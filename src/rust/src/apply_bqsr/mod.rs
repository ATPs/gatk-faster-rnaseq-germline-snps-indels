use anyhow::{bail, Context, Result};
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use rust_htslib::bam;
use rust_htslib::bam::record::{Aux, Cigar, Record};
use rust_htslib::bam::Read;
use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::OnceLock;

const BASE_SUBSTITUTION_EVENT: &str = "M";
const DEFAULT_PRESERVE_QSCORES_LESS_THAN: u8 = 6;
const DEFAULT_LOW_QUALITY_TAIL: u8 = 2;
const DEFAULT_MISMATCHES_CONTEXT_SIZE: usize = 2;
const DEFAULT_MAXIMUM_CYCLE_VALUE: i32 = 500;
const MAX_REASONABLE_Q_SCORE: i32 = 60;
const MAX_RECALIBRATED_Q_SCORE: i32 = 93;
const MAX_NUMBER_OF_OBSERVATIONS: u64 = 2_147_483_646;
const APPLY_BQSR_BATCH_RECORDS: usize = 16_384;

#[derive(Debug, Clone)]
pub struct ApplyBqsrConfig {
    pub input_bam: PathBuf,
    pub recal_table: PathBuf,
    pub output_bam: PathBuf,
    pub output_index: Option<PathBuf>,
    pub threads: usize,
    pub use_original_qualities: bool,
    pub allow_missing_read_groups: bool,
    pub use_report_quantization: bool,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct ApplyBqsrStats {
    pub records: u64,
    pub bases: u64,
}

pub fn apply_bqsr(config: &ApplyBqsrConfig) -> Result<ApplyBqsrStats> {
    if !config.input_bam.exists() {
        bail!("input BAM not found: {}", config.input_bam.display());
    }
    if !config.recal_table.exists() {
        bail!(
            "recalibration table not found: {}",
            config.recal_table.display()
        );
    }
    if let Some(parent) = config.output_bam.parent() {
        if !parent.as_os_str().is_empty() {
            fs::create_dir_all(parent).with_context(|| format!("creating {}", parent.display()))?;
        }
    }

    let model = RecalibrationModel::from_path(&config.recal_table, config.use_report_quantization)?;
    let mut reader = bam::Reader::from_path(&config.input_bam)
        .with_context(|| format!("opening {}", config.input_bam.display()))?;
    let total_threads = config.threads.max(1);
    let htslib_threads = if total_threads >= 8 {
        4
    } else {
        total_threads.min(2)
    };
    if htslib_threads > 1 {
        reader
            .set_threads(htslib_threads)
            .context("enabling BAM reader threads")?;
    }

    let read_group_ids = read_group_identifiers_from_header(reader.header().as_bytes())?;
    let header = bam::Header::from_template(reader.header());
    let mut writer = bam::Writer::from_path(&config.output_bam, &header, bam::Format::Bam)
        .with_context(|| format!("creating {}", config.output_bam.display()))?;
    if htslib_threads > 1 {
        writer
            .set_threads(htslib_threads)
            .context("enabling BAM writer threads")?;
    }

    let worker_threads = total_threads.saturating_sub(htslib_threads).max(1);
    let thread_pool = ThreadPoolBuilder::new()
        .num_threads(worker_threads)
        .build()
        .context("creating ApplyBQSR worker thread pool")?;
    let mut stats = ApplyBqsrStats {
        records: 0,
        bases: 0,
    };

    let mut batch = Vec::with_capacity(APPLY_BQSR_BATCH_RECORDS);
    for record_result in reader.records() {
        batch.push(record_result.context("reading BAM record")?);
        if batch.len() == APPLY_BQSR_BATCH_RECORDS {
            process_apply_bqsr_batch(
                &mut batch,
                &thread_pool,
                &model,
                &read_group_ids,
                config.use_original_qualities,
                config.allow_missing_read_groups,
                &mut writer,
                &mut stats,
            )?;
        }
    }
    process_apply_bqsr_batch(
        &mut batch,
        &thread_pool,
        &model,
        &read_group_ids,
        config.use_original_qualities,
        config.allow_missing_read_groups,
        &mut writer,
        &mut stats,
    )?;

    drop(writer);
    let output_index = config
        .output_index
        .clone()
        .unwrap_or_else(|| bam_index_path(&config.output_bam));
    bam::index::build(
        &config.output_bam,
        Some(&output_index),
        bam::index::Type::Bai,
        total_threads as u32,
    )
    .with_context(|| format!("creating BAM index {}", output_index.display()))?;

    Ok(stats)
}

fn process_apply_bqsr_batch(
    batch: &mut Vec<Record>,
    thread_pool: &rayon::ThreadPool,
    model: &RecalibrationModel,
    read_group_ids: &HashMap<String, String>,
    use_original_qualities: bool,
    allow_missing_read_groups: bool,
    writer: &mut bam::Writer,
    stats: &mut ApplyBqsrStats,
) -> Result<()> {
    if batch.is_empty() {
        return Ok(());
    }

    let records = std::mem::replace(batch, Vec::with_capacity(APPLY_BQSR_BATCH_RECORDS));
    let transformed = thread_pool.install(|| {
        records
            .into_par_iter()
            .map(|mut record| {
                if !passes_apply_bqsr_read_filters(&record) {
                    return Ok(None);
                }
                transform_record(
                    &mut record,
                    model,
                    read_group_ids,
                    use_original_qualities,
                    allow_missing_read_groups,
                )
                .map(|()| Some(record))
            })
            .collect::<Vec<_>>()
    });

    for record_result in transformed {
        let Some(record) = record_result? else {
            continue;
        };
        stats.records += 1;
        stats.bases += record.seq_len() as u64;
        writer.write(&record).context("writing BAM record")?;
    }

    Ok(())
}

fn passes_apply_bqsr_read_filters(record: &Record) -> bool {
    !record
        .cigar()
        .iter()
        .any(|op| matches!(op, Cigar::RefSkip(_)))
}

pub fn bam_index_path(bam: &Path) -> PathBuf {
    if bam.extension().and_then(|ext| ext.to_str()) == Some("bam") {
        bam.with_extension("bai")
    } else {
        PathBuf::from(format!("{}.bai", bam.display()))
    }
}

fn transform_record(
    record: &mut Record,
    model: &RecalibrationModel,
    read_group_ids: &HashMap<String, String>,
    use_original_qualities: bool,
    allow_missing_read_groups: bool,
) -> Result<()> {
    let rg_id = aux_string(record, b"RG")?
        .with_context(|| format!("read {} has no RG tag", qname_for_error(record)))?;
    let read_group_identifier = read_group_ids
        .get(&rg_id)
        .with_context(|| format!("read group {rg_id} is not present in the BAM header"))?;

    let bases = record.seq().as_bytes();
    let mut qualities = if use_original_qualities {
        match aux_string(record, b"OQ")? {
            Some(oq) => fastq_to_phred(oq.as_bytes())
                .with_context(|| format!("invalid OQ tag on read {}", qname_for_error(record)))?,
            None => record.qual().to_vec(),
        }
    } else {
        record.qual().to_vec()
    };

    if qualities.len() != bases.len() {
        bail!(
            "read {} has {} bases but {} qualities",
            qname_for_error(record),
            bases.len(),
            qualities.len()
        );
    }

    qualities = model.recalibrate_qualities(
        read_group_identifier,
        &bases,
        &qualities,
        record.flags(),
        allow_missing_read_groups,
    )?;

    remove_aux_if_present(record, b"BI")?;
    remove_aux_if_present(record, b"BD")?;

    let qname = record.qname().to_vec();
    let cigar = record.cigar().take();
    record.set(&qname, Some(&cigar), &bases, &qualities);
    Ok(())
}

fn qname_for_error(record: &Record) -> String {
    String::from_utf8_lossy(record.qname()).into_owned()
}

fn aux_string(record: &Record, tag: &[u8]) -> Result<Option<String>> {
    match record.aux(tag) {
        Ok(Aux::String(value)) => Ok(Some(value.to_owned())),
        Ok(_) => bail!(
            "aux tag {} is present but is not a string",
            String::from_utf8_lossy(tag)
        ),
        Err(_) => Ok(None),
    }
}

fn remove_aux_if_present(record: &mut Record, tag: &[u8]) -> Result<()> {
    if record.aux(tag).is_ok() {
        record
            .remove_aux(tag)
            .with_context(|| format!("removing aux tag {}", String::from_utf8_lossy(tag)))?;
    }
    Ok(())
}

fn fastq_to_phred(value: &[u8]) -> Result<Vec<u8>> {
    value
        .iter()
        .map(|&b| {
            if b < 33 {
                bail!("FASTQ quality byte below '!': {b}");
            }
            Ok(b - 33)
        })
        .collect()
}

fn read_group_identifiers_from_header(header: &[u8]) -> Result<HashMap<String, String>> {
    let header = std::str::from_utf8(header).context("BAM header is not UTF-8")?;
    let mut read_groups = HashMap::new();
    for line in header.lines() {
        if !line.starts_with("@RG\t") {
            continue;
        }
        let mut id = None;
        let mut platform_unit = None;
        for field in line.split('\t').skip(1) {
            if let Some(value) = field.strip_prefix("ID:") {
                id = Some(value.to_owned());
            } else if let Some(value) = field.strip_prefix("PU:") {
                platform_unit = Some(value.to_owned());
            }
        }
        let id = id.with_context(|| format!("malformed @RG line without ID: {line}"))?;
        let identifier = platform_unit.unwrap_or_else(|| id.clone());
        read_groups.insert(id, identifier);
    }
    Ok(read_groups)
}

include!("model.rs");
include!("tests_block.rs");
