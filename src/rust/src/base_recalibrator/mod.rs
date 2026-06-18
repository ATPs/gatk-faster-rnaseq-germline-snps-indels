use anyhow::{anyhow, bail, Context, Result};
use clap::Parser;
use rust_htslib::bam;
use rust_htslib::bam::record::{Aux, Cigar};
use rust_htslib::bam::Read;
use rust_htslib::bcf;
use rust_htslib::bcf::Read as BcfRead;
use rust_htslib::faidx;
use std::cmp::Ordering;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering as AtomicOrdering};
use std::sync::Arc;
use std::thread;

const MAX_SAM_QUAL_SCORE: usize = 93;
const MAX_REASONABLE_Q_SCORE: i32 = 60;
const MIN_USABLE_Q_SCORE: u8 = 6;
const MAX_GATK_USABLE_Q_SCORE: i32 = 40;

#[derive(Clone, Debug, Parser)]
#[command(
    name = "rust_base_recalibrator",
    about = "Mismatch-only GATK BaseRecalibrator replacement for RNA-seq preprocessing"
)]
struct Args {
    #[arg(short = 'R', long = "ref")]
    reference: PathBuf,

    #[arg(short = 'I', long = "input-bam")]
    input_bam: PathBuf,

    #[arg(long = "known-sites", required = true)]
    known_sites: Vec<PathBuf>,

    #[arg(short = 'O', long = "output-table")]
    output_table: PathBuf,

    #[arg(long = "use-original-qualities")]
    use_original_qualities: bool,

    #[arg(long = "mismatches-context-size", default_value_t = 2)]
    context_size: usize,

    #[arg(long = "low-quality-tail", default_value_t = 2)]
    low_quality_tail: u8,

    #[arg(long = "maximum-cycle-value", default_value_t = 500)]
    max_cycle: i32,

    #[arg(long = "quantizing-levels", default_value_t = 16)]
    quantizing_levels: usize,

    #[arg(long = "known-sites-chunk-size", default_value_t = 1_000_000)]
    chunk_size: u64,

    #[arg(long = "threads", default_value_t = 1)]
    threads: usize,

    #[arg(long = "region-bases", default_value_t = 25_000_000)]
    region_bases: u64,
}

#[derive(Clone, Debug)]
struct ReadGroupInfo {
    identifiers: Vec<String>,
    id_to_identifier: HashMap<String, String>,
}

#[derive(Clone, Debug)]
struct PreparedRead {
    contig: String,
    start0: u64,
    end0: u64,
    bases: Vec<u8>,
    quals: Vec<u8>,
    cigar: Vec<SimpleCigar>,
    rg_index: usize,
    is_reverse: bool,
    is_second_in_pair: bool,
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum SimpleCigar {
    Match(u32),
    Ins(u32),
    Del(u32),
    RefSkip(u32),
    Equal(u32),
    Diff(u32),
}

#[derive(Clone, Debug)]
struct KnownInterval {
    start0: u64,
    end0: u64,
}

#[derive(Clone, Debug)]
struct WorkUnit {
    tid: u32,
    contig: String,
    start0: u64,
    end0_exclusive: u64,
}

struct ReferenceCache {
    reader: faidx::Reader,
    contig: Option<String>,
    start0: u64,
    end0: u64,
    bases: Vec<u8>,
    chunk_size: u64,
}

struct KnownSiteReader {
    path: PathBuf,
    reader: bcf::IndexedReader,
    unfetchable_contigs: HashSet<String>,
}

struct KnownSiteCache {
    readers: Vec<KnownSiteReader>,
    contig: Option<String>,
    start0: u64,
    end0: u64,
    intervals: Vec<KnownInterval>,
    chunk_size: u64,
}

#[derive(Clone, Debug)]
struct Datum {
    observations: u64,
    errors: f64,
    reported_quality: f64,
}

impl Datum {
    fn new(reported_quality: u8) -> Self {
        Self {
            observations: 0,
            errors: 0.0,
            reported_quality: f64::from(reported_quality),
        }
    }

    fn increment(&mut self, is_error: f64) {
        self.observations += 1;
        self.errors += is_error;
    }

    fn combine(&mut self, other: &Datum) {
        let expected_errors = self.expected_errors() + other.expected_errors();
        self.observations += other.observations;
        self.errors += other.errors;
        self.reported_quality = -10.0 * (expected_errors / self.observations as f64).log10();
    }

    fn add_same_key(&mut self, other: &Datum) {
        self.observations += other.observations;
        self.errors += other.errors;
    }

    fn expected_errors(&self) -> f64 {
        self.observations as f64 * qual_to_error_prob(self.reported_quality)
    }
}

#[derive(Default)]
struct RecalTables {
    quality: BTreeMap<(usize, u8), Datum>,
    cycle: BTreeMap<(usize, u8, i32), Datum>,
    context: BTreeMap<(usize, u8, String), Datum>,
}

impl RecalTables {
    fn merge_from(&mut self, other: RecalTables) {
        merge_table(&mut self.quality, other.quality);
        merge_table(&mut self.cycle, other.cycle);
        merge_table(&mut self.context, other.context);
    }
}

fn merge_table<K: Ord>(target: &mut BTreeMap<K, Datum>, source: BTreeMap<K, Datum>) {
    for (key, datum) in source {
        target
            .entry(key)
            .and_modify(|existing| existing.add_same_key(&datum))
            .or_insert(datum);
    }
}

pub fn run_cli() -> Result<()> {
    let args = Args::parse();
    run(args)
}

fn run(args: Args) -> Result<()> {
    if args.context_size == 0 || args.context_size > 13 {
        bail!("--mismatches-context-size must be between 1 and 13");
    }
    if args.max_cycle < 1 {
        bail!("--maximum-cycle-value must be at least 1");
    }
    if args.quantizing_levels == 0 {
        bail!("--quantizing-levels must be at least 1");
    }
    if args.chunk_size == 0 {
        bail!("--known-sites-chunk-size must be at least 1");
    }
    if args.threads == 0 {
        bail!("--threads must be at least 1");
    }
    if args.region_bases == 0 {
        bail!("--region-bases must be at least 1");
    }

    let mut bam_reader = bam::Reader::from_path(&args.input_bam)
        .with_context(|| format!("failed to open BAM {}", args.input_bam.display()))?;
    let header = bam_reader.header().clone();
    let rg_info = read_groups_from_header(&header)?;
    let tables = if args.threads == 1 {
        run_sequential(&args, &mut bam_reader, &header, &rg_info)?
    } else {
        drop(bam_reader);
        run_parallel(&args, &header, &rg_info)?
    };

    let read_group_table = collapse_read_group_table(&tables.quality);
    let quantized = build_quantization_table(&tables.quality, args.quantizing_levels);
    write_report(&args, &rg_info, &read_group_table, &tables, &quantized)?;
    Ok(())
}

fn run_sequential(
    args: &Args,
    bam_reader: &mut bam::Reader,
    header: &bam::HeaderView,
    rg_info: &ReadGroupInfo,
) -> Result<RecalTables> {
    let mut reference = ReferenceCache::new(&args.reference, args.chunk_size)?;
    let mut known_sites = KnownSiteCache::new(&args.known_sites, args.chunk_size)?;
    let mut tables = RecalTables::default();
    for result in bam_reader.records() {
        let record =
            result.with_context(|| format!("failed to read BAM {}", args.input_bam.display()))?;
        process_record(
            record,
            header,
            rg_info,
            &mut reference,
            &mut known_sites,
            args,
            &mut tables,
        )?;
    }
    Ok(tables)
}

fn run_parallel(
    args: &Args,
    header: &bam::HeaderView,
    rg_info: &ReadGroupInfo,
) -> Result<RecalTables> {
    let work_units = Arc::new(build_work_units(header, args.region_bases)?);
    if work_units.is_empty() {
        return Ok(RecalTables::default());
    }
    let next = Arc::new(AtomicUsize::new(0));
    let worker_count = args.threads.min(work_units.len());
    let mut handles = Vec::with_capacity(worker_count);

    for worker_index in 0..worker_count {
        let args = args.clone();
        let rg_info = rg_info.clone();
        let work_units = Arc::clone(&work_units);
        let next = Arc::clone(&next);
        handles.push(thread::spawn(move || -> Result<RecalTables> {
            let mut bam_reader =
                bam::IndexedReader::from_path(&args.input_bam).with_context(|| {
                    format!("failed to open indexed BAM {}", args.input_bam.display())
                })?;
            let header = bam_reader.header().clone();
            let mut reference = ReferenceCache::new(&args.reference, args.chunk_size)?;
            let mut known_sites = KnownSiteCache::new(&args.known_sites, args.chunk_size)?;
            let mut tables = RecalTables::default();
            loop {
                let work_index = next.fetch_add(1, AtomicOrdering::Relaxed);
                if work_index >= work_units.len() {
                    break;
                }
                process_work_unit(
                    &args,
                    worker_index,
                    &work_units[work_index],
                    &mut bam_reader,
                    &header,
                    &rg_info,
                    &mut reference,
                    &mut known_sites,
                    &mut tables,
                )?;
            }
            Ok(tables)
        }));
    }

    let mut merged = RecalTables::default();
    for handle in handles {
        let worker_tables = handle
            .join()
            .map_err(|_| anyhow!("rust_base_recalibrator worker thread panicked"))??;
        merged.merge_from(worker_tables);
    }
    Ok(merged)
}

include!("work_units.rs");
include!("reporting.rs");
include!("tests_block.rs");
