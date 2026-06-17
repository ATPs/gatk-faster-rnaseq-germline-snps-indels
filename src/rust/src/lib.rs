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

pub mod assembly;
pub mod haplotype_caller;
pub mod hc_tools;
pub mod interval_tools;
pub mod pair_hmm;
pub mod smith_waterman;
pub mod split_n_cigar;

pub fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}

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

#[derive(Debug)]
struct RecalibrationModel {
    args: RecalibrationArgs,
    read_groups: HashMap<String, usize>,
    read_group_table: Vec<Option<RecalDatum>>,
    quality_score_table: HashMap<(usize, u8), RecalDatum>,
    cycle_table: HashMap<(usize, u8, i32), RecalDatum>,
    context_table: HashMap<(usize, u8, u32), RecalDatum>,
    quantized_quals: Vec<u8>,
}

#[derive(Debug, Clone, Copy)]
struct RecalibrationArgs {
    preserve_qscores_less_than: u8,
    mismatches_context_size: usize,
    low_quality_tail: u8,
    maximum_cycle_value: i32,
}

impl Default for RecalibrationArgs {
    fn default() -> Self {
        Self {
            preserve_qscores_less_than: DEFAULT_PRESERVE_QSCORES_LESS_THAN,
            mismatches_context_size: DEFAULT_MISMATCHES_CONTEXT_SIZE,
            low_quality_tail: DEFAULT_LOW_QUALITY_TAIL,
            maximum_cycle_value: DEFAULT_MAXIMUM_CYCLE_VALUE,
        }
    }
}

impl RecalibrationModel {
    fn from_path(path: &Path, use_report_quantization: bool) -> Result<Self> {
        let text = fs::read_to_string(path)
            .with_context(|| format!("reading recalibration table {}", path.display()))?;
        Self::from_str(&text, use_report_quantization)
    }

    fn from_str(text: &str, use_report_quantization: bool) -> Result<Self> {
        let mut parser = RecalibrationTableParser::new(use_report_quantization);
        parser.parse(text)?;
        parser.finish()
    }

    fn recalibrate_qualities(
        &self,
        read_group_identifier: &str,
        bases: &[u8],
        qualities: &[u8],
        flags: u16,
        allow_missing_read_groups: bool,
    ) -> Result<Vec<u8>> {
        let Some(&rg_key) = self.read_groups.get(read_group_identifier) else {
            if allow_missing_read_groups {
                return Ok(qualities
                    .iter()
                    .map(|&q| self.quantized_quals[q as usize])
                    .collect());
            }
            bail!("read group {read_group_identifier} not found in recalibration table");
        };

        let read_group_datum = self
            .read_group_table
            .get(rg_key)
            .and_then(Option::as_ref)
            .with_context(|| {
                format!("read group {read_group_identifier} has no substitution row")
            })?;
        let prior_quality = read_group_datum.reported_quality;
        let empirical_quality_for_read_group = read_group_datum.empirical_quality(prior_quality);
        let contexts = context_keys_for_read(
            bases,
            qualities,
            flags,
            self.args.mismatches_context_size,
            self.args.low_quality_tail,
        );

        let mut recalibrated = qualities.to_vec();
        for offset in 0..qualities.len() {
            let reported_quality = qualities[offset];
            if reported_quality < self.args.preserve_qscores_less_than {
                continue;
            }

            let posterior_quality = self
                .quality_score_table
                .get(&(rg_key, reported_quality))
                .map(|datum| datum.empirical_quality(empirical_quality_for_read_group))
                .unwrap_or(empirical_quality_for_read_group);

            let mut delta_special_covariates = 0.0;
            let cycle = cycle_for_offset(
                offset,
                qualities.len(),
                flags,
                self.args.maximum_cycle_value,
            )?;
            if let Some(datum) = self.cycle_table.get(&(rg_key, reported_quality, cycle)) {
                delta_special_covariates +=
                    datum.empirical_quality(posterior_quality) - posterior_quality;
            }
            if let Some(context_key) = contexts[offset] {
                if let Some(datum) =
                    self.context_table
                        .get(&(rg_key, reported_quality, context_key))
                {
                    delta_special_covariates +=
                        datum.empirical_quality(posterior_quality) - posterior_quality;
                }
            }

            let raw_quality = posterior_quality + delta_special_covariates;
            let bounded_quality = bound_qual(fast_round(raw_quality), MAX_RECALIBRATED_Q_SCORE);
            recalibrated[offset] = self.quantized_quals[bounded_quality as usize];
        }
        Ok(recalibrated)
    }
}

struct RecalibrationTableParser {
    args: RecalibrationArgs,
    read_groups: HashMap<String, usize>,
    read_group_table: Vec<Option<RecalDatum>>,
    quality_score_table: HashMap<(usize, u8), RecalDatum>,
    cycle_table: HashMap<(usize, u8, i32), RecalDatum>,
    context_table: HashMap<(usize, u8, u32), RecalDatum>,
    report_quantized_quals: Vec<Option<u8>>,
    use_report_quantization: bool,
}

impl RecalibrationTableParser {
    fn new(use_report_quantization: bool) -> Self {
        Self {
            args: RecalibrationArgs::default(),
            read_groups: HashMap::new(),
            read_group_table: Vec::new(),
            quality_score_table: HashMap::new(),
            cycle_table: HashMap::new(),
            context_table: HashMap::new(),
            report_quantized_quals: vec![None; (MAX_RECALIBRATED_Q_SCORE + 1) as usize],
            use_report_quantization,
        }
    }

    fn parse(&mut self, text: &str) -> Result<()> {
        let mut current_table: Option<&str> = None;
        let mut expect_header = false;

        for raw_line in text.lines() {
            let line = raw_line.trim_end();
            if line.is_empty() {
                continue;
            }
            if let Some(rest) = line.strip_prefix("#:GATKTable:") {
                let table_name = rest.split(':').next().unwrap_or_default();
                current_table = match table_name {
                    "Arguments" | "Quantized" | "RecalTable0" | "RecalTable1" | "RecalTable2" => {
                        Some(table_name)
                    }
                    _ => None,
                };
                expect_header = current_table.is_some();
                continue;
            }
            if line.starts_with('#') {
                continue;
            }
            let Some(table) = current_table else {
                continue;
            };
            if expect_header {
                expect_header = false;
                continue;
            }

            match table {
                "Arguments" => self.parse_argument_row(line)?,
                "Quantized" => self.parse_quantized_row(line)?,
                "RecalTable0" => self.parse_read_group_row(line)?,
                "RecalTable1" => self.parse_quality_score_row(line)?,
                "RecalTable2" => self.parse_covariate_row(line)?,
                _ => {}
            }
        }

        Ok(())
    }

    fn finish(self) -> Result<RecalibrationModel> {
        if self.read_groups.is_empty() {
            bail!("recalibration table has no read groups");
        }

        let quantized_quals = if self.use_report_quantization {
            self.report_quantized_quals
                .into_iter()
                .enumerate()
                .map(|(q, value)| value.with_context(|| format!("missing quantized row for Q{q}")))
                .collect::<Result<Vec<_>>>()?
        } else {
            (0..=MAX_RECALIBRATED_Q_SCORE).map(|q| q as u8).collect()
        };

        Ok(RecalibrationModel {
            args: self.args,
            read_groups: self.read_groups,
            read_group_table: self.read_group_table,
            quality_score_table: self.quality_score_table,
            cycle_table: self.cycle_table,
            context_table: self.context_table,
            quantized_quals,
        })
    }

    fn parse_argument_row(&mut self, line: &str) -> Result<()> {
        let mut fields = line.split_whitespace();
        let Some(argument) = fields.next() else {
            return Ok(());
        };
        let value = fields.next().unwrap_or("null");
        match argument {
            "covariate" => {
                let expected =
                    "ReadGroupCovariate,QualityScoreCovariate,ContextCovariate,CycleCovariate";
                if value != expected {
                    bail!(
                        "unsupported covariate list {value}; only standard BQSR covariates are supported"
                    );
                }
            }
            "no_standard_covs" if value == "true" => {
                bail!("recalibration tables with no_standard_covs=true are not supported");
            }
            "mismatches_context_size" => {
                self.args.mismatches_context_size = parse_field(value, argument)?;
            }
            "low_quality_tail" => {
                self.args.low_quality_tail = parse_field(value, argument)?;
            }
            "maximum_cycle_value" => {
                self.args.maximum_cycle_value = parse_field(value, argument)?;
            }
            _ => {}
        }
        Ok(())
    }

    fn parse_quantized_row(&mut self, line: &str) -> Result<()> {
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 3 {
            bail!("malformed Quantized row: {line}");
        }
        let quality: usize = parse_field(fields[0], "QualityScore")?;
        let quantized: u8 = parse_field(fields[2], "QuantizedScore")?;
        if quality < self.report_quantized_quals.len() {
            self.report_quantized_quals[quality] = Some(quantized);
        }
        Ok(())
    }

    fn parse_read_group_row(&mut self, line: &str) -> Result<()> {
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 6 {
            bail!("malformed RecalTable0 row: {line}");
        }
        if fields[1] != BASE_SUBSTITUTION_EVENT {
            return Ok(());
        }
        let rg = self.read_group_key(fields[0]);
        let reported_quality = parse_field(fields[3], "EstimatedQReported")?;
        let observations = parse_field(fields[4], "Observations")?;
        let errors = parse_field(fields[5], "Errors")?;
        self.read_group_table[rg] = Some(RecalDatum::new(observations, errors, reported_quality));
        Ok(())
    }

    fn parse_quality_score_row(&mut self, line: &str) -> Result<()> {
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 6 {
            bail!("malformed RecalTable1 row: {line}");
        }
        if fields[2] != BASE_SUBSTITUTION_EVENT {
            return Ok(());
        }
        let rg = self.read_group_key(fields[0]);
        let quality = parse_field(fields[1], "QualityScore")?;
        let observations = parse_field(fields[4], "Observations")?;
        let errors = parse_field(fields[5], "Errors")?;
        self.quality_score_table.insert(
            (rg, quality),
            RecalDatum::new(observations, errors, f64::from(quality)),
        );
        Ok(())
    }

    fn parse_covariate_row(&mut self, line: &str) -> Result<()> {
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 8 {
            bail!("malformed RecalTable2 row: {line}");
        }
        if fields[4] != BASE_SUBSTITUTION_EVENT {
            return Ok(());
        }
        let rg = self.read_group_key(fields[0]);
        let quality = parse_field(fields[1], "QualityScore")?;
        let observations = parse_field(fields[6], "Observations")?;
        let errors = parse_field(fields[7], "Errors")?;
        let datum = RecalDatum::new(observations, errors, f64::from(quality));
        match fields[3] {
            "Cycle" => {
                let cycle = parse_field(fields[2], "Cycle")?;
                self.cycle_table.insert((rg, quality, cycle), datum);
            }
            "Context" => {
                let context_key = key_from_context(fields[2].as_bytes())
                    .with_context(|| format!("invalid Context covariate {}", fields[2]))?;
                self.context_table.insert((rg, quality, context_key), datum);
            }
            covariate => bail!("unsupported covariate {covariate} in RecalTable2"),
        }
        Ok(())
    }

    fn read_group_key(&mut self, read_group: &str) -> usize {
        if let Some(&key) = self.read_groups.get(read_group) {
            return key;
        }
        let key = self.read_groups.len();
        self.read_groups.insert(read_group.to_owned(), key);
        self.read_group_table.push(None);
        key
    }
}

fn parse_field<T>(value: &str, name: &str) -> Result<T>
where
    T: std::str::FromStr,
    T::Err: std::error::Error + Send + Sync + 'static,
{
    value
        .parse()
        .with_context(|| format!("parsing {name} value {value}"))
}

#[derive(Debug, Clone)]
struct RecalDatum {
    observations: u64,
    errors: f64,
    reported_quality: f64,
    empirical_quality: OnceLock<i32>,
}

impl RecalDatum {
    fn new(observations: u64, errors: f64, reported_quality: f64) -> Self {
        Self {
            observations,
            errors,
            reported_quality,
            empirical_quality: OnceLock::new(),
        }
    }

    fn empirical_quality(&self, prior_quality_score: f64) -> f64 {
        f64::from(*self.empirical_quality.get_or_init(|| {
            let mut observations = self.observations + 2;
            let mut errors = (self.errors + 0.5) as u64 + 1;
            if observations > MAX_NUMBER_OF_OBSERVATIONS {
                let fraction = MAX_NUMBER_OF_OBSERVATIONS as f64 / observations as f64;
                errors = (errors as f64 * fraction).round() as u64;
                observations = MAX_NUMBER_OF_OBSERVATIONS;
            }

            bayesian_estimate_of_empirical_quality(observations, errors, prior_quality_score)
                .min(MAX_RECALIBRATED_Q_SCORE)
        }))
    }
}

fn bayesian_estimate_of_empirical_quality(
    observations: u64,
    errors: u64,
    prior_quality_score: f64,
) -> i32 {
    let mut best_quality = 0;
    let mut best_log_posterior = f64::NEG_INFINITY;
    for quality in 0..=MAX_REASONABLE_Q_SCORE {
        let log_posterior = log_prior(quality, prior_quality_score)
            + log_binomial_likelihood(quality, observations, errors);
        if log_posterior > best_log_posterior {
            best_log_posterior = log_posterior;
            best_quality = quality;
        }
    }
    best_quality
}

fn log_prior(quality: i32, prior_quality_score: f64) -> f64 {
    let difference = ((f64::from(quality) - prior_quality_score) as i32)
        .abs()
        .min(40);
    let sigma = 0.5;
    -0.5 * (f64::from(difference) / sigma).powi(2)
}

fn log_binomial_likelihood(quality: i32, observations: u64, errors: u64) -> f64 {
    if observations == 0 {
        return 0.0;
    }
    let p_error = 10_f64.powf(f64::from(quality) / -10.0);
    if p_error == 1.0 {
        return if errors == observations {
            0.0
        } else {
            f64::NEG_INFINITY
        };
    }
    if p_error == 0.0 {
        return if errors == 0 { 0.0 } else { f64::NEG_INFINITY };
    }

    errors as f64 * p_error.ln() + (observations - errors) as f64 * (-p_error).ln_1p()
}

fn fast_round(value: f64) -> i32 {
    if value > 0.0 {
        (value + 0.5) as i32
    } else {
        (value - 0.5) as i32
    }
}

fn bound_qual(quality: i32, max_quality: i32) -> i32 {
    quality.clamp(1, max_quality)
}

fn context_keys_for_read(
    bases: &[u8],
    qualities: &[u8],
    flags: u16,
    context_size: usize,
    low_quality_tail: u8,
) -> Vec<Option<u32>> {
    let read_len = bases.len();
    let mut contexts = vec![None; read_len];
    if read_len == 0 || context_size == 0 {
        return contexts;
    }

    let Some(mut stranded_bases) = stranded_low_quality_tail_clipped_bases(
        bases,
        qualities,
        is_reverse(flags),
        low_quality_tail,
    ) else {
        return contexts;
    };

    if stranded_bases.len() != read_len {
        return contexts;
    }

    for stranded_offset in 0..read_len {
        let read_offset = if is_reverse(flags) {
            read_len - stranded_offset - 1
        } else {
            stranded_offset
        };
        if stranded_offset + 1 < context_size {
            continue;
        }
        let start = stranded_offset + 1 - context_size;
        contexts[read_offset] = key_from_context(&stranded_bases[start..=stranded_offset]);
    }

    stranded_bases.clear();
    contexts
}

fn stranded_low_quality_tail_clipped_bases(
    bases: &[u8],
    qualities: &[u8],
    reverse: bool,
    low_quality_tail: u8,
) -> Option<Vec<u8>> {
    let read_len = bases.len();
    let mut right = read_len;
    while right > 0 && qualities[right - 1] <= low_quality_tail {
        right -= 1;
    }
    let mut left = 0;
    while left < read_len && qualities[left] <= low_quality_tail {
        left += 1;
    }
    if left >= right {
        return None;
    }

    let mut clipped = bases.to_vec();
    for base in &mut clipped[..left] {
        *base = b'N';
    }
    for base in &mut clipped[right..] {
        *base = b'N';
    }

    if reverse {
        clipped = clipped.into_iter().rev().map(simple_complement).collect();
    }
    Some(clipped)
}

fn cycle_for_offset(
    offset: usize,
    read_len: usize,
    flags: u16,
    maximum_cycle_value: i32,
) -> Result<i32> {
    let read_order_factor = if is_paired(flags) && is_second_in_pair(flags) {
        -1
    } else {
        1
    };
    let (cycle_start, increment) = if is_reverse(flags) {
        (read_len as i32 * read_order_factor, -read_order_factor)
    } else {
        (read_order_factor, read_order_factor)
    };
    let cycle = cycle_start + offset as i32 * increment;
    if cycle.abs() > maximum_cycle_value {
        bail!(
            "cycle {} exceeds maximum_cycle_value {}",
            cycle,
            maximum_cycle_value
        );
    }
    Ok(cycle)
}

fn is_paired(flags: u16) -> bool {
    flags & 0x1 != 0
}

fn is_reverse(flags: u16) -> bool {
    flags & 0x10 != 0
}

fn is_second_in_pair(flags: u16) -> bool {
    flags & 0x80 != 0
}

fn key_from_context(context: &[u8]) -> Option<u32> {
    let mut key = context.len() as u32;
    let mut bit_offset = 4;
    for &base in context {
        let base_index = simple_base_to_index(base)?;
        key |= base_index << bit_offset;
        bit_offset += 2;
    }
    Some(key)
}

fn simple_base_to_index(base: u8) -> Option<u32> {
    match base {
        b'A' | b'a' | b'*' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

fn simple_complement(base: u8) -> u8 {
    match base {
        b'A' | b'a' => b'T',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        b'T' | b't' => b'A',
        _ => base,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::{Cigar, CigarString};

    #[test]
    fn synthetic_quality_recalibration_preserves_low_qualities() {
        let model = RecalibrationModel::from_str(
            &synthetic_table(
                &["rg1  M  30.0000  30.0000  1000000  1000.00"],
                &["rg1  30  M  35.0000  1000000  316.00"],
                &[],
            ),
            false,
        )
        .unwrap();

        let recalibrated = model
            .recalibrate_qualities("rg1", b"AAC", &[30, 30, 5], 0, false)
            .unwrap();
        assert_eq!(recalibrated, vec![34, 34, 5]);
    }

    #[test]
    fn synthetic_cycle_and_context_covariates_are_applied() {
        let model = RecalibrationModel::from_str(
            &synthetic_table(
                &["rg1  M  30.0000  30.0000  1000000  1000.00"],
                &["rg1  30  M  30.0000  1000000  1000.00"],
                &[
                    "rg1  30  1   Cycle    M  40.0000  1000000  100.00",
                    "rg1  30  AA  Context  M  40.0000  1000000  100.00",
                ],
            ),
            false,
        )
        .unwrap();

        let recalibrated = model
            .recalibrate_qualities("rg1", b"AA", &[30, 30], 0, false)
            .unwrap();
        assert_eq!(recalibrated, vec![37, 37]);
    }

    #[test]
    fn header_read_group_identifier_prefers_platform_unit() {
        let header = b"@HD\tVN:1.6\n@RG\tID:rg1\tSM:s1\tPU:unit1\n@RG\tID:rg2\tSM:s1\n";
        let read_groups = read_group_identifiers_from_header(header).unwrap();
        assert_eq!(read_groups.get("rg1").unwrap(), "unit1");
        assert_eq!(read_groups.get("rg2").unwrap(), "rg2");
    }

    #[test]
    fn transform_record_rewrites_qualities_and_preserves_alignment_fields() {
        let model = RecalibrationModel::from_str(
            &synthetic_table(
                &["rg1  M  30.0000  30.0000  1000000  1000.00"],
                &["rg1  30  M  35.0000  1000000  316.00"],
                &[],
            ),
            false,
        )
        .unwrap();
        let read_groups = HashMap::from([("rg1".to_owned(), "rg1".to_owned())]);
        let cigar = CigarString(vec![Cigar::Match(3)]);
        let mut record = Record::new();
        record.set(b"read1", Some(&cigar), b"AAC", &[10, 10, 10]);
        record.set_tid(2);
        record.set_pos(123);
        record.set_mtid(2);
        record.set_mpos(456);
        record.set_insert_size(333);
        record.set_flags(0x41);
        record.set_mapq(60);
        record.push_aux(b"RG", Aux::String("rg1")).unwrap();
        record.push_aux(b"OQ", Aux::String("???")).unwrap();
        record.push_aux(b"NM", Aux::I32(1)).unwrap();
        record.push_aux(b"BI", Aux::String("old")).unwrap();
        record.push_aux(b"BD", Aux::String("old")).unwrap();

        transform_record(&mut record, &model, &read_groups, true, false).unwrap();

        assert_eq!(record.qname(), b"read1");
        assert_eq!(record.seq().as_bytes(), b"AAC");
        assert_eq!(record.qual(), &[34, 34, 34]);
        assert_eq!(record.tid(), 2);
        assert_eq!(record.pos(), 123);
        assert_eq!(record.mtid(), 2);
        assert_eq!(record.mpos(), 456);
        assert_eq!(record.insert_size(), 333);
        assert_eq!(record.flags(), 0x41);
        assert_eq!(record.mapq(), 60);
        assert!(matches!(record.aux(b"RG").unwrap(), Aux::String("rg1")));
        assert!(matches!(record.aux(b"OQ").unwrap(), Aux::String("???")));
        assert!(matches!(record.aux(b"NM").unwrap(), Aux::I32(1)));
        assert!(record.aux(b"BI").is_err());
        assert!(record.aux(b"BD").is_err());
    }

    #[test]
    fn apply_bqsr_read_filter_rejects_cigar_n_records() {
        let mut record = Record::new();
        record.set(
            b"read1",
            Some(&CigarString(vec![
                Cigar::Match(2),
                Cigar::RefSkip(10),
                Cigar::Match(2),
            ])),
            b"AACC",
            &[30, 30, 30, 30],
        );

        assert!(!passes_apply_bqsr_read_filters(&record));
    }

    fn synthetic_table(
        read_group_rows: &[&str],
        quality_rows: &[&str],
        covariate_rows: &[&str],
    ) -> String {
        let mut text = String::new();
        text.push_str("#:GATKReport.v1.1:5\n");
        text.push_str("#:GATKTable:2:17:%s:%s:;\n");
        text.push_str(
            "#:GATKTable:Arguments:Recalibration argument collection values used in this run\n",
        );
        text.push_str("Argument                    Value\n");
        text.push_str("covariate                   ReadGroupCovariate,QualityScoreCovariate,ContextCovariate,CycleCovariate\n");
        text.push_str("low_quality_tail            2\n");
        text.push_str("maximum_cycle_value         500\n");
        text.push_str("mismatches_context_size     2\n");
        text.push_str("no_standard_covs            false\n\n");

        text.push_str("#:GATKTable:3:94:%d:%d:%d:;\n");
        text.push_str("#:GATKTable:Quantized:Quality quantization map\n");
        text.push_str("QualityScore  Count      QuantizedScore\n");
        for quality in 0..=MAX_RECALIBRATED_Q_SCORE {
            text.push_str(&format!("{quality}  0  {quality}\n"));
        }
        text.push('\n');

        text.push_str("#:GATKTable:6:1:%s:%s:%.4f:%.4f:%d:%.2f:;\n");
        text.push_str("#:GATKTable:RecalTable0:\n");
        text.push_str("ReadGroup        EventType  EmpiricalQuality  EstimatedQReported  Observations  Errors\n");
        for row in read_group_rows {
            text.push_str(row);
            text.push('\n');
        }
        text.push('\n');

        text.push_str("#:GATKTable:6:1:%s:%d:%s:%.4f:%d:%.2f:;\n");
        text.push_str("#:GATKTable:RecalTable1:\n");
        text.push_str(
            "ReadGroup        QualityScore  EventType  EmpiricalQuality  Observations  Errors\n",
        );
        for row in quality_rows {
            text.push_str(row);
            text.push('\n');
        }
        text.push('\n');

        text.push_str("#:GATKTable:8:1:%s:%d:%s:%s:%s:%.4f:%d:%.2f:;\n");
        text.push_str("#:GATKTable:RecalTable2:\n");
        text.push_str("ReadGroup        QualityScore  CovariateValue  CovariateName  EventType  EmpiricalQuality  Observations  Errors\n");
        for row in covariate_rows {
            text.push_str(row);
            text.push('\n');
        }
        text
    }
}
