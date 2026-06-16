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

fn main() -> Result<()> {
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

fn build_work_units(header: &bam::HeaderView, region_bases: u64) -> Result<Vec<WorkUnit>> {
    let mut work_units = Vec::new();
    for tid in 0..header.target_count() {
        let contig = String::from_utf8(header.tid2name(tid).to_vec())
            .with_context(|| format!("BAM header target id {tid} is not valid UTF-8"))?;
        let Some(contig_len) = header.target_len(tid) else {
            continue;
        };
        let mut start0 = 0_u64;
        while start0 < contig_len {
            let end0_exclusive = start0.saturating_add(region_bases).min(contig_len);
            work_units.push(WorkUnit {
                tid,
                contig: contig.clone(),
                start0,
                end0_exclusive,
            });
            start0 = end0_exclusive;
        }
    }
    Ok(work_units)
}

#[allow(clippy::too_many_arguments)]
fn process_work_unit(
    args: &Args,
    worker_index: usize,
    work: &WorkUnit,
    bam_reader: &mut bam::IndexedReader,
    header: &bam::HeaderView,
    rg_info: &ReadGroupInfo,
    reference: &mut ReferenceCache,
    known_sites: &mut KnownSiteCache,
    tables: &mut RecalTables,
) -> Result<()> {
    let start = i64::try_from(work.start0)
        .with_context(|| format!("{} start coordinate exceeds i64", work.contig))?;
    let end = i64::try_from(work.end0_exclusive)
        .with_context(|| format!("{} end coordinate exceeds i64", work.contig))?;
    bam_reader
        .fetch((work.tid as i32, start, end))
        .with_context(|| {
            format!(
                "worker {worker_index} failed to fetch BAM region {}:{}-{}",
                work.contig,
                work.start0 + 1,
                work.end0_exclusive
            )
        })?;

    for result in bam_reader.records() {
        let record = result.with_context(|| {
            format!(
                "worker {worker_index} failed to read BAM region {}:{}-{}",
                work.contig,
                work.start0 + 1,
                work.end0_exclusive
            )
        })?;
        if record.tid() != work.tid as i32 || record.pos() < start || record.pos() >= end {
            continue;
        }
        process_record(
            record,
            header,
            rg_info,
            reference,
            known_sites,
            args,
            tables,
        )?;
    }
    Ok(())
}

fn process_record(
    record: bam::Record,
    header: &bam::HeaderView,
    rg_info: &ReadGroupInfo,
    reference: &mut ReferenceCache,
    known_sites: &mut KnownSiteCache,
    args: &Args,
    tables: &mut RecalTables,
) -> Result<()> {
    if !passes_bqsr_filters(&record) {
        return Ok(());
    }
    let Some(read) = prepare_read(&record, header, rg_info, args)? else {
        return Ok(());
    };
    let reference_bases = reference.bases_for(&read.contig, read.start0, read.end0)?;
    let known = known_sites.intervals_for(&read.contig, read.start0, read.end0)?;
    process_read(&read, reference_bases, known, args, tables)
}

fn passes_bqsr_filters(record: &bam::Record) -> bool {
    !record.is_unmapped()
        && !record.is_secondary()
        && !record.is_duplicate()
        && !record.is_quality_check_failed()
        && record.mapq() != 0
        && record.mapq() != 255
        && record.tid() >= 0
        && record.pos() >= 0
}

fn read_groups_from_header(header: &bam::HeaderView) -> Result<ReadGroupInfo> {
    let text = String::from_utf8_lossy(header.as_bytes());
    let mut identifiers = Vec::new();
    let mut id_to_identifier = HashMap::new();
    let mut identifier_to_index = HashMap::new();
    for line in text.lines() {
        if !line.starts_with("@RG\t") {
            continue;
        }
        let mut id = None;
        let mut platform_unit = None;
        for field in line.split('\t').skip(1) {
            if let Some(value) = field.strip_prefix("ID:") {
                id = Some(value.to_string());
            } else if let Some(value) = field.strip_prefix("PU:") {
                platform_unit = Some(value.to_string());
            }
        }
        let id = id.context("encountered @RG header line without ID")?;
        let identifier = platform_unit.unwrap_or_else(|| id.clone());
        if !identifier_to_index.contains_key(&identifier) {
            identifier_to_index.insert(identifier.clone(), identifiers.len());
            identifiers.push(identifier.clone());
        }
        id_to_identifier.insert(id, identifier);
    }

    if identifiers.is_empty() {
        bail!("input BAM header has no read groups");
    }
    Ok(ReadGroupInfo {
        identifiers,
        id_to_identifier,
    })
}

fn prepare_read(
    record: &bam::Record,
    header: &bam::HeaderView,
    rg_info: &ReadGroupInfo,
    args: &Args,
) -> Result<Option<PreparedRead>> {
    let contig = String::from_utf8_lossy(header.tid2name(record.tid() as u32)).to_string();
    let rg_id = aux_string(record, b"RG")?;
    let rg_identifier = rg_info
        .id_to_identifier
        .get(&rg_id)
        .with_context(|| format!("read references unknown read group '{rg_id}'"))?;
    let rg_index = rg_info
        .identifiers
        .iter()
        .position(|identifier| identifier == rg_identifier)
        .with_context(|| format!("read group '{rg_identifier}' missing from lookup"))?;

    let bases = record.seq().as_bytes();
    let quals = if args.use_original_qualities {
        original_qualities(record)?.unwrap_or_else(|| record.qual().to_vec())
    } else {
        record.qual().to_vec()
    };
    if quals.len() != bases.len() {
        bail!(
            "read {} has {} bases but {} qualities",
            String::from_utf8_lossy(record.qname()),
            bases.len(),
            quals.len()
        );
    }

    let mut clipped_bases = Vec::with_capacity(bases.len());
    let mut clipped_quals = Vec::with_capacity(quals.len());
    let mut cigar = Vec::new();
    let mut read_pos = 0_usize;
    let mut ref_len = 0_u64;

    for op in &record.cigar() {
        match *op {
            Cigar::Match(len) => {
                copy_read_segment(
                    &bases,
                    &quals,
                    read_pos,
                    len,
                    &mut clipped_bases,
                    &mut clipped_quals,
                );
                read_pos += len as usize;
                ref_len += u64::from(len);
                cigar.push(SimpleCigar::Match(len));
            }
            Cigar::Equal(len) => {
                copy_read_segment(
                    &bases,
                    &quals,
                    read_pos,
                    len,
                    &mut clipped_bases,
                    &mut clipped_quals,
                );
                read_pos += len as usize;
                ref_len += u64::from(len);
                cigar.push(SimpleCigar::Equal(len));
            }
            Cigar::Diff(len) => {
                copy_read_segment(
                    &bases,
                    &quals,
                    read_pos,
                    len,
                    &mut clipped_bases,
                    &mut clipped_quals,
                );
                read_pos += len as usize;
                ref_len += u64::from(len);
                cigar.push(SimpleCigar::Diff(len));
            }
            Cigar::Ins(len) => {
                copy_read_segment(
                    &bases,
                    &quals,
                    read_pos,
                    len,
                    &mut clipped_bases,
                    &mut clipped_quals,
                );
                read_pos += len as usize;
                cigar.push(SimpleCigar::Ins(len));
            }
            Cigar::Del(len) => {
                ref_len += u64::from(len);
                cigar.push(SimpleCigar::Del(len));
            }
            Cigar::RefSkip(len) => {
                ref_len += u64::from(len);
                cigar.push(SimpleCigar::RefSkip(len));
            }
            Cigar::SoftClip(len) => {
                read_pos += len as usize;
            }
            Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }

    if clipped_bases.is_empty() || ref_len == 0 {
        return Ok(None);
    }
    let start0 = record.pos() as u64;
    let end0 = start0 + ref_len - 1;

    Ok(Some(PreparedRead {
        contig,
        start0,
        end0,
        bases: clipped_bases,
        quals: clipped_quals,
        cigar,
        rg_index,
        is_reverse: record.is_reverse(),
        is_second_in_pair: record.is_paired() && record.is_last_in_template(),
    }))
}

fn copy_read_segment(
    bases: &[u8],
    quals: &[u8],
    read_pos: usize,
    len: u32,
    clipped_bases: &mut Vec<u8>,
    clipped_quals: &mut Vec<u8>,
) {
    let end = read_pos + len as usize;
    clipped_bases.extend_from_slice(&bases[read_pos..end]);
    clipped_quals.extend_from_slice(&quals[read_pos..end]);
}

fn aux_string(record: &bam::Record, tag: &[u8]) -> Result<String> {
    match record.aux(tag) {
        Ok(Aux::String(value)) => Ok(value.to_string()),
        Ok(Aux::Char(value)) => Ok((value as char).to_string()),
        Ok(other) => bail!(
            "aux tag {} has unsupported value type {:?}",
            String::from_utf8_lossy(tag),
            other
        ),
        Err(_) => bail!(
            "read {} is missing required aux tag {}",
            String::from_utf8_lossy(record.qname()),
            String::from_utf8_lossy(tag)
        ),
    }
}

fn original_qualities(record: &bam::Record) -> Result<Option<Vec<u8>>> {
    match record.aux(b"OQ") {
        Ok(Aux::String(value)) => Ok(Some(value.bytes().map(|b| b.saturating_sub(33)).collect())),
        Ok(_) => bail!(
            "read {} has non-string OQ tag",
            String::from_utf8_lossy(record.qname())
        ),
        Err(_) => Ok(None),
    }
}

impl ReferenceCache {
    fn new(reference: &PathBuf, chunk_size: u64) -> Result<Self> {
        let reader = faidx::Reader::from_path(reference)
            .with_context(|| format!("failed to open FASTA index for {}", reference.display()))?;
        Ok(Self {
            reader,
            contig: None,
            start0: 0,
            end0: 0,
            bases: Vec::new(),
            chunk_size,
        })
    }

    fn bases_for(&mut self, contig: &str, start0: u64, end0: u64) -> Result<&[u8]> {
        let cache_hit =
            self.contig.as_deref() == Some(contig) && start0 >= self.start0 && end0 <= self.end0;
        if !cache_hit {
            self.load_chunk(contig, start0, end0)?;
        }
        let offset = usize::try_from(start0 - self.start0)
            .with_context(|| format!("reference offset overflow for {contig}:{start0}-{end0}"))?;
        let length = usize::try_from(end0 - start0 + 1)
            .with_context(|| format!("reference length overflow for {contig}:{start0}-{end0}"))?;
        Ok(&self.bases[offset..offset + length])
    }

    fn load_chunk(&mut self, contig: &str, start0: u64, end0: u64) -> Result<()> {
        self.contig = Some(contig.to_string());
        self.start0 = start0;
        let contig_len = self.reader.fetch_seq_len(contig);
        if contig_len == 0 {
            bail!("contig '{contig}' is not present in the reference FASTA");
        }
        self.end0 = end0
            .max(start0.saturating_add(self.chunk_size).saturating_sub(1))
            .min(contig_len - 1);
        self.bases = self
            .reader
            .fetch_seq_string(contig, self.start0 as usize, self.end0 as usize)
            .with_context(|| {
                format!(
                    "failed to fetch reference {}:{}-{}",
                    contig,
                    self.start0 + 1,
                    self.end0 + 1
                )
            })?
            .into_bytes();
        Ok(())
    }
}

impl KnownSiteCache {
    fn new(paths: &[PathBuf], chunk_size: u64) -> Result<Self> {
        let mut readers = Vec::new();
        for path in paths {
            let reader = bcf::IndexedReader::from_path(path).with_context(|| {
                format!("failed to open indexed known-sites VCF {}", path.display())
            })?;
            readers.push(KnownSiteReader {
                path: path.clone(),
                reader,
                unfetchable_contigs: HashSet::new(),
            });
        }
        Ok(Self {
            readers,
            contig: None,
            start0: 0,
            end0: 0,
            intervals: Vec::new(),
            chunk_size,
        })
    }

    fn intervals_for(&mut self, contig: &str, start0: u64, end0: u64) -> Result<&[KnownInterval]> {
        let cache_hit =
            self.contig.as_deref() == Some(contig) && start0 >= self.start0 && end0 <= self.end0;
        if !cache_hit {
            self.load_chunk(contig, start0, end0)?;
        }
        Ok(&self.intervals)
    }

    fn load_chunk(&mut self, contig: &str, start0: u64, end0: u64) -> Result<()> {
        self.contig = Some(contig.to_string());
        self.start0 = start0;
        self.end0 = end0.max(start0.saturating_add(self.chunk_size).saturating_sub(1));
        self.intervals.clear();

        for known_reader in &mut self.readers {
            if known_reader.unfetchable_contigs.contains(contig) {
                continue;
            }
            let rid = match known_reader.reader.header().name2rid(contig.as_bytes()) {
                Ok(rid) => rid,
                Err(_) => continue,
            };
            if let Err(err) = known_reader
                .reader
                .fetch(rid, self.start0, Some(self.end0 + 1))
            {
                eprintln!(
                    "warning: skipping known-sites VCF {} on unfetchable contig {}: {}",
                    known_reader.path.display(),
                    contig,
                    err
                );
                known_reader.unfetchable_contigs.insert(contig.to_string());
                continue;
            }
            for record_result in known_reader.reader.records() {
                let record = record_result.with_context(|| {
                    format!(
                        "failed to read known-sites VCF {}",
                        known_reader.path.display()
                    )
                })?;
                let rec_start = record.pos().max(0) as u64;
                let rec_len = record.rlen().max(1) as u64;
                let rec_end = rec_start + rec_len - 1;
                if rec_end >= self.start0 && rec_start <= self.end0 {
                    self.intervals.push(KnownInterval {
                        start0: rec_start,
                        end0: rec_end,
                    });
                }
            }
        }
        self.intervals
            .sort_by(|a, b| a.start0.cmp(&b.start0).then(a.end0.cmp(&b.end0)));
        Ok(())
    }
}

fn process_read(
    read: &PreparedRead,
    reference_bases: &[u8],
    known_intervals: &[KnownInterval],
    args: &Args,
    tables: &mut RecalTables,
) -> Result<()> {
    let contexts = context_values(read, args.context_size, args.low_quality_tail);
    let mut read_pos = 0_usize;
    let mut ref_pos = 0_usize;
    let mut known_index = 0_usize;

    for op in &read.cigar {
        match *op {
            SimpleCigar::Match(len) | SimpleCigar::Equal(len) | SimpleCigar::Diff(len) => {
                for _ in 0..len {
                    let base = read.bases[read_pos];
                    let qual = read.quals[read_pos].min(MAX_SAM_QUAL_SCORE as u8);
                    let ref_base = reference_bases.get(ref_pos).copied().with_context(|| {
                        format!(
                            "reference slice shorter than CIGAR span for {}:{}-{}",
                            read.contig,
                            read.start0 + 1,
                            read.end0 + 1
                        )
                    })?;
                    let absolute_ref = read.start0 + ref_pos as u64;
                    if is_regular_base(base)
                        && qual >= MIN_USABLE_Q_SCORE
                        && !is_known_site(absolute_ref, known_intervals, &mut known_index)
                    {
                        let is_error = if bases_equal(base, ref_base) {
                            0.0
                        } else {
                            1.0
                        };
                        update_tables(read, read_pos, qual, is_error, &contexts, args, tables)?;
                    }
                    read_pos += 1;
                    ref_pos += 1;
                }
            }
            SimpleCigar::Ins(len) => {
                read_pos += len as usize;
            }
            SimpleCigar::Del(len) | SimpleCigar::RefSkip(len) => {
                ref_pos += len as usize;
            }
        }
    }
    Ok(())
}

fn update_tables(
    read: &PreparedRead,
    read_pos: usize,
    qual: u8,
    is_error: f64,
    contexts: &[Option<String>],
    args: &Args,
    tables: &mut RecalTables,
) -> Result<()> {
    tables
        .quality
        .entry((read.rg_index, qual))
        .or_insert_with(|| Datum::new(qual))
        .increment(is_error);
    let cycle = cycle_value(
        read_pos,
        read.bases.len(),
        read.is_reverse,
        read.is_second_in_pair,
        args.max_cycle,
    )?;
    tables
        .cycle
        .entry((read.rg_index, qual, cycle))
        .or_insert_with(|| Datum::new(qual))
        .increment(is_error);
    if let Some(context) = &contexts[read_pos] {
        tables
            .context
            .entry((read.rg_index, qual, context.clone()))
            .or_insert_with(|| Datum::new(qual))
            .increment(is_error);
    }
    Ok(())
}

fn collapse_read_group_table(
    quality_table: &BTreeMap<(usize, u8), Datum>,
) -> BTreeMap<usize, Datum> {
    let mut read_groups: BTreeMap<usize, Datum> = BTreeMap::new();
    for ((rg, _qual), datum) in quality_table {
        read_groups
            .entry(*rg)
            .and_modify(|existing| existing.combine(datum))
            .or_insert_with(|| datum.clone());
    }
    read_groups
}

fn is_known_site(position0: u64, intervals: &[KnownInterval], index: &mut usize) -> bool {
    while *index < intervals.len() && intervals[*index].end0 < position0 {
        *index += 1;
    }
    *index < intervals.len()
        && intervals[*index].start0 <= position0
        && intervals[*index].end0 >= position0
}

fn context_values(
    read: &PreparedRead,
    context_size: usize,
    low_quality_tail: u8,
) -> Vec<Option<String>> {
    let mut bases = read.bases.clone();
    let mut left = 0_usize;
    while left < bases.len() && read.quals[left] <= low_quality_tail {
        bases[left] = b'N';
        left += 1;
    }
    let mut right = bases.len();
    while right > left && read.quals[right - 1] <= low_quality_tail {
        bases[right - 1] = b'N';
        right -= 1;
    }
    if read.is_reverse {
        bases = reverse_complement(&bases);
    }

    let mut stranded_contexts = vec![None; bases.len()];
    if bases.len() >= context_size {
        for i in (context_size - 1)..bases.len() {
            let start = i + 1 - context_size;
            if bases[start..=i].iter().all(|base| is_regular_base(*base)) {
                stranded_contexts[i] = Some(
                    bases[start..=i]
                        .iter()
                        .map(|base| normalize_base(*base) as char)
                        .collect(),
                );
            }
        }
    }

    let mut contexts = vec![None; bases.len()];
    for (stranded_offset, context) in stranded_contexts.into_iter().enumerate() {
        let read_offset = if read.is_reverse {
            bases.len() - stranded_offset - 1
        } else {
            stranded_offset
        };
        contexts[read_offset] = context;
    }
    contexts
}

fn cycle_value(
    base_number: usize,
    read_length: usize,
    is_reverse: bool,
    is_second_in_pair: bool,
    max_cycle: i32,
) -> Result<i32> {
    let read_order_factor = if is_second_in_pair { -1 } else { 1 };
    let (mut cycle, increment) = if is_reverse {
        (read_length as i32 * read_order_factor, -read_order_factor)
    } else {
        (read_order_factor, read_order_factor)
    };
    cycle += base_number as i32 * increment;
    if cycle.abs() > max_cycle {
        bail!(
            "cycle {} exceeds --maximum-cycle-value {}",
            cycle,
            max_cycle
        );
    }
    Ok(cycle)
}

fn is_regular_base(base: u8) -> bool {
    matches!(normalize_base(base), b'A' | b'C' | b'G' | b'T')
}

fn bases_equal(left: u8, right: u8) -> bool {
    normalize_base(left) == normalize_base(right)
}

fn normalize_base(base: u8) -> u8 {
    base.to_ascii_uppercase()
}

fn reverse_complement(bases: &[u8]) -> Vec<u8> {
    bases
        .iter()
        .rev()
        .map(|base| match normalize_base(*base) {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => b'N',
        })
        .collect()
}

fn build_quantization_table(
    quality_table: &BTreeMap<(usize, u8), Datum>,
    quantizing_levels: usize,
) -> Vec<(usize, u64, u8)> {
    let mut counts = vec![0_u64; MAX_SAM_QUAL_SCORE + 1];
    for datum in quality_table.values() {
        let empirical = empirical_quality(datum).round() as usize;
        counts[empirical.min(MAX_SAM_QUAL_SCORE)] += datum.observations;
    }
    let quantized =
        quantize_quality_scores(&counts, quantizing_levels, usize::from(MIN_USABLE_Q_SCORE));
    (0..=MAX_SAM_QUAL_SCORE)
        .map(|qual| (qual, counts[qual], quantized[qual]))
        .collect()
}

#[derive(Clone, Debug)]
struct QualInterval {
    q_start: usize,
    q_end: usize,
    fixed_qual: Option<u8>,
    observations: u64,
    errors: u64,
    sub_intervals: Vec<QualInterval>,
}

impl QualInterval {
    fn leaf(q: usize, observations: u64) -> Self {
        let errors = (observations as f64 * qual_to_error_prob(q as f64)).floor() as u64;
        Self {
            q_start: q,
            q_end: q,
            fixed_qual: Some(q as u8),
            observations,
            errors,
            sub_intervals: Vec::new(),
        }
    }

    fn merge(left: &Self, right: &Self) -> Self {
        Self {
            q_start: left.q_start,
            q_end: right.q_end,
            fixed_qual: None,
            observations: left.observations + right.observations,
            errors: left.errors + right.errors,
            sub_intervals: vec![left.clone(), right.clone()],
        }
    }

    fn error_rate(&self) -> f64 {
        if let Some(qual) = self.fixed_qual {
            qual_to_error_prob(f64::from(qual))
        } else if self.observations == 0 {
            0.0
        } else {
            (self.errors + 1) as f64 / (self.observations + 1) as f64
        }
    }

    fn qual(&self) -> u8 {
        if let Some(qual) = self.fixed_qual {
            qual
        } else {
            error_prob_to_qual(self.error_rate())
        }
    }

    fn penalty(&self, min_interesting_qual: usize) -> f64 {
        self.calc_penalty(self.error_rate(), min_interesting_qual)
    }

    fn calc_penalty(&self, global_error_rate: f64, min_interesting_qual: usize) -> f64 {
        if global_error_rate == 0.0 {
            return 0.0;
        }
        if self.sub_intervals.is_empty() {
            if self.q_end <= min_interesting_qual {
                0.0
            } else {
                (self.error_rate().log10() - global_error_rate.log10()).abs()
                    * self.observations as f64
            }
        } else {
            self.sub_intervals
                .iter()
                .map(|interval| interval.calc_penalty(global_error_rate, min_interesting_qual))
                .sum()
        }
    }
}

fn quantize_quality_scores(
    counts: &[u64],
    n_levels: usize,
    min_interesting_qual: usize,
) -> Vec<u8> {
    let mut intervals: Vec<QualInterval> = counts
        .iter()
        .enumerate()
        .map(|(qual, count)| QualInterval::leaf(qual, *count))
        .collect();

    while intervals.len() > n_levels {
        let mut best_index = 0_usize;
        let mut best_penalty = f64::INFINITY;
        for index in 0..(intervals.len() - 1) {
            let merged = QualInterval::merge(&intervals[index], &intervals[index + 1]);
            let penalty = merged.penalty(min_interesting_qual);
            if penalty < best_penalty {
                best_index = index;
                best_penalty = penalty;
            }
        }
        let merged = QualInterval::merge(&intervals[best_index], &intervals[best_index + 1]);
        intervals.splice(best_index..=best_index + 1, [merged]);
    }

    let mut map = vec![0_u8; counts.len()];
    for interval in intervals {
        for qual in interval.q_start..=interval.q_end {
            map[qual] = interval.qual();
        }
    }
    map
}

fn empirical_quality(datum: &Datum) -> f64 {
    let mismatches = (datum.errors + 0.5) as u64 + 1;
    let observations = datum.observations + 2;
    bayesian_estimate_of_empirical_quality(observations, mismatches, datum.reported_quality)
        .min(MAX_SAM_QUAL_SCORE as i32) as f64
}

fn bayesian_estimate_of_empirical_quality(
    observations: u64,
    errors: u64,
    prior_mean_quality: f64,
) -> i32 {
    let mut best_q = 0_i32;
    let mut best_value = f64::NEG_INFINITY;
    for q in 0..=MAX_REASONABLE_Q_SCORE {
        let value = log_prior(q as f64, prior_mean_quality)
            + log_binomial_likelihood(q as f64, observations, errors);
        if value > best_value {
            best_value = value;
            best_q = q;
        }
    }
    best_q
}

fn log_prior(quality_score: f64, prior_quality_score: f64) -> f64 {
    let difference = ((quality_score - prior_quality_score) as i32)
        .abs()
        .min(MAX_GATK_USABLE_Q_SCORE);
    let sigma = 0.5_f64;
    -((difference as f64).powi(2)) / (2.0 * sigma.powi(2))
}

fn log_binomial_likelihood(quality_score: f64, observations: u64, mut errors: u64) -> f64 {
    if observations == 0 {
        return 0.0;
    }
    let max_observations = i32::MAX as u64 - 1;
    let mut scaled_observations = observations;
    if scaled_observations > max_observations {
        let fraction = max_observations as f64 / scaled_observations as f64;
        errors = (errors as f64 * fraction).round() as u64;
        scaled_observations = max_observations;
    }

    let p = qual_to_error_prob(quality_score);
    if p == 1.0 {
        return if errors == scaled_observations {
            0.0
        } else {
            f64::NEG_INFINITY
        };
    }
    errors as f64 * p.ln() + (scaled_observations - errors) as f64 * (1.0 - p).ln()
}

fn qual_to_error_prob(qual: f64) -> f64 {
    10.0_f64.powf(qual / -10.0)
}

fn error_prob_to_qual(error_rate: f64) -> u8 {
    let qual = (-10.0 * error_rate.log10()).round() as i32;
    qual.clamp(1, MAX_SAM_QUAL_SCORE as i32) as u8
}

fn write_report(
    args: &Args,
    rg_info: &ReadGroupInfo,
    read_group_table: &BTreeMap<usize, Datum>,
    tables: &RecalTables,
    quantized: &[(usize, u64, u8)],
) -> Result<()> {
    if let Some(parent) = args.output_table.parent() {
        if !parent.as_os_str().is_empty() {
            fs::create_dir_all(parent)
                .with_context(|| format!("failed to create {}", parent.display()))?;
        }
    }
    let file = File::create(&args.output_table)
        .with_context(|| format!("failed to create {}", args.output_table.display()))?;
    let mut writer = BufWriter::new(file);
    writeln!(writer, "#:GATKReport.v1.1:5")?;

    write_table(
        &mut writer,
        "Arguments",
        "Recalibration argument collection values used in this run",
        &[("Argument", "%s"), ("Value", "%s")],
        argument_rows(args),
    )?;

    let quant_rows: Vec<Vec<String>> = quantized
        .iter()
        .map(|(qual, count, quantized_score)| {
            vec![
                qual.to_string(),
                count.to_string(),
                quantized_score.to_string(),
            ]
        })
        .collect();
    write_table(
        &mut writer,
        "Quantized",
        "Quality quantization map",
        &[
            ("QualityScore", "%d"),
            ("Count", "%d"),
            ("QuantizedScore", "%d"),
        ],
        quant_rows,
    )?;

    let rg_rows: Vec<Vec<String>> = read_group_table
        .iter()
        .map(|(rg, datum)| {
            vec![
                rg_info.identifiers[*rg].clone(),
                "M".to_string(),
                format!("{:.4}", empirical_quality(datum)),
                format!("{:.4}", datum.reported_quality),
                datum.observations.to_string(),
                format!("{:.2}", datum.errors),
            ]
        })
        .collect();
    write_table(
        &mut writer,
        "RecalTable0",
        "",
        &[
            ("ReadGroup", "%s"),
            ("EventType", "%s"),
            ("EmpiricalQuality", "%.4f"),
            ("EstimatedQReported", "%.4f"),
            ("Observations", "%d"),
            ("Errors", "%.2f"),
        ],
        rg_rows,
    )?;

    let qual_rows: Vec<Vec<String>> = tables
        .quality
        .iter()
        .map(|((rg, qual), datum)| {
            vec![
                rg_info.identifiers[*rg].clone(),
                qual.to_string(),
                "M".to_string(),
                format!("{:.4}", empirical_quality(datum)),
                datum.observations.to_string(),
                format!("{:.2}", datum.errors),
            ]
        })
        .collect();
    write_table(
        &mut writer,
        "RecalTable1",
        "",
        &[
            ("ReadGroup", "%s"),
            ("QualityScore", "%d"),
            ("EventType", "%s"),
            ("EmpiricalQuality", "%.4f"),
            ("Observations", "%d"),
            ("Errors", "%.2f"),
        ],
        qual_rows,
    )?;

    let mut cov_rows = Vec::with_capacity(tables.cycle.len() + tables.context.len());
    for ((rg, qual, cycle), datum) in &tables.cycle {
        cov_rows.push(vec![
            rg_info.identifiers[*rg].clone(),
            qual.to_string(),
            cycle.to_string(),
            "Cycle".to_string(),
            "M".to_string(),
            format!("{:.4}", empirical_quality(datum)),
            datum.observations.to_string(),
            format!("{:.2}", datum.errors),
        ]);
    }
    for ((rg, qual, context), datum) in &tables.context {
        cov_rows.push(vec![
            rg_info.identifiers[*rg].clone(),
            qual.to_string(),
            context.clone(),
            "Context".to_string(),
            "M".to_string(),
            format!("{:.4}", empirical_quality(datum)),
            datum.observations.to_string(),
            format!("{:.2}", datum.errors),
        ]);
    }
    cov_rows.sort_by(compare_recal_table2_rows);
    write_table(
        &mut writer,
        "RecalTable2",
        "",
        &[
            ("ReadGroup", "%s"),
            ("QualityScore", "%d"),
            ("CovariateValue", "%s"),
            ("CovariateName", "%s"),
            ("EventType", "%s"),
            ("EmpiricalQuality", "%.4f"),
            ("Observations", "%d"),
            ("Errors", "%.2f"),
        ],
        cov_rows,
    )?;
    writer.flush()?;
    Ok(())
}

fn argument_rows(args: &Args) -> Vec<Vec<String>> {
    vec![
        vec!["binary_tag_name".to_string(), "null".to_string()],
        vec![
            "covariate".to_string(),
            "ReadGroupCovariate,QualityScoreCovariate,ContextCovariate,CycleCovariate".to_string(),
        ],
        vec!["default_platform".to_string(), "null".to_string()],
        vec!["deletions_default_quality".to_string(), "45".to_string()],
        vec!["force_platform".to_string(), "null".to_string()],
        vec!["indels_context_size".to_string(), "3".to_string()],
        vec!["insertions_default_quality".to_string(), "45".to_string()],
        vec![
            "low_quality_tail".to_string(),
            args.low_quality_tail.to_string(),
        ],
        vec![
            "maximum_cycle_value".to_string(),
            args.max_cycle.to_string(),
        ],
        vec![
            "mismatches_context_size".to_string(),
            args.context_size.to_string(),
        ],
        vec!["mismatches_default_quality".to_string(), "-1".to_string()],
        vec!["no_standard_covs".to_string(), "false".to_string()],
        vec![
            "quantizing_levels".to_string(),
            args.quantizing_levels.to_string(),
        ],
        vec!["recalibration_report".to_string(), "null".to_string()],
        vec!["run_without_dbsnp".to_string(), "false".to_string()],
        vec![
            "solid_nocall_strategy".to_string(),
            "THROW_EXCEPTION".to_string(),
        ],
        vec!["solid_recal_mode".to_string(), "SET_Q_ZERO".to_string()],
    ]
}

fn compare_recal_table2_rows(left: &Vec<String>, right: &Vec<String>) -> Ordering {
    left[0]
        .cmp(&right[0])
        .then_with(|| numeric_string_cmp(&left[1], &right[1]))
        .then_with(|| left[2].cmp(&right[2]))
        .then_with(|| left[3].cmp(&right[3]))
}

fn numeric_string_cmp(left: &str, right: &str) -> Ordering {
    match (left.parse::<i32>(), right.parse::<i32>()) {
        (Ok(l), Ok(r)) => l.cmp(&r),
        _ => left.cmp(right),
    }
}

fn write_table<W: Write>(
    writer: &mut W,
    name: &str,
    description: &str,
    columns: &[(&str, &str)],
    rows: Vec<Vec<String>>,
) -> std::io::Result<()> {
    writeln!(
        writer,
        "#:GATKTable:{}:{}:{}:;",
        columns.len(),
        rows.len(),
        columns
            .iter()
            .map(|(_, format)| *format)
            .collect::<Vec<&str>>()
            .join(":")
    )?;
    writeln!(writer, "#:GATKTable:{name}:{description}")?;

    let mut widths: Vec<usize> = columns.iter().map(|(column, _)| column.len()).collect();
    for row in &rows {
        for (idx, value) in row.iter().enumerate() {
            widths[idx] = widths[idx].max(value.len());
        }
    }
    write_fixed_row(
        writer,
        &columns
            .iter()
            .map(|(column, _)| column.to_string())
            .collect::<Vec<String>>(),
        &widths,
        true,
    )?;
    for row in &rows {
        write_fixed_row(writer, row, &widths, false)?;
    }
    writeln!(writer)?;
    Ok(())
}

fn write_fixed_row<W: Write>(
    writer: &mut W,
    row: &[String],
    widths: &[usize],
    header: bool,
) -> std::io::Result<()> {
    for (idx, value) in row.iter().enumerate() {
        if idx > 0 {
            write!(writer, "  ")?;
        }
        if header || !right_align(value) {
            write!(writer, "{:<width$}", value, width = widths[idx])?;
        } else {
            write!(writer, "{:>width$}", value, width = widths[idx])?;
        }
    }
    writeln!(writer)
}

fn right_align(value: &str) -> bool {
    value == "null" || value == "NA" || value.parse::<i64>().is_ok() || value.parse::<f64>().is_ok()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::BTreeSet;

    fn test_read(reverse: bool) -> PreparedRead {
        PreparedRead {
            contig: "chr1".to_string(),
            start0: 0,
            end0: 3,
            bases: b"ACGT".to_vec(),
            quals: vec![30, 30, 30, 30],
            cigar: vec![SimpleCigar::Match(4)],
            rg_index: 0,
            is_reverse: reverse,
            is_second_in_pair: false,
        }
    }

    #[test]
    fn cycle_matches_gatk_orientation_rules() {
        assert_eq!(cycle_value(0, 5, false, false, 500).unwrap(), 1);
        assert_eq!(cycle_value(4, 5, false, false, 500).unwrap(), 5);
        assert_eq!(cycle_value(0, 5, true, false, 500).unwrap(), 5);
        assert_eq!(cycle_value(4, 5, true, false, 500).unwrap(), 1);
        assert_eq!(cycle_value(0, 5, false, true, 500).unwrap(), -1);
        assert_eq!(cycle_value(4, 5, true, true, 500).unwrap(), -1);
    }

    #[test]
    fn context_uses_preceding_base_and_reverse_complements_negative_reads() {
        assert_eq!(
            context_values(&test_read(false), 2, 2),
            vec![
                None,
                Some("AC".to_string()),
                Some("CG".to_string()),
                Some("GT".to_string())
            ]
        );
        assert_eq!(
            context_values(&test_read(true), 2, 2),
            vec![
                Some("GT".to_string()),
                Some("CG".to_string()),
                Some("AC".to_string()),
                None
            ]
        );
    }

    #[test]
    fn known_site_pointer_masks_overlapping_positions() {
        let intervals = vec![
            KnownInterval {
                start0: 10,
                end0: 12,
            },
            KnownInterval {
                start0: 20,
                end0: 20,
            },
        ];
        let mut index = 0;
        assert!(!is_known_site(9, &intervals, &mut index));
        assert!(is_known_site(10, &intervals, &mut index));
        assert!(is_known_site(12, &intervals, &mut index));
        assert!(!is_known_site(13, &intervals, &mut index));
        assert!(is_known_site(20, &intervals, &mut index));
    }

    #[test]
    fn empirical_quality_moves_down_for_high_error_bins() {
        let mut clean = Datum::new(30);
        for _ in 0..100 {
            clean.increment(0.0);
        }
        let mut noisy = Datum::new(30);
        for _ in 0..100 {
            noisy.increment(1.0);
        }
        assert!(empirical_quality(&clean) > empirical_quality(&noisy));
    }

    #[test]
    fn quantizer_keeps_requested_number_of_levels() {
        let mut counts = vec![0_u64; MAX_SAM_QUAL_SCORE + 1];
        counts[10] = 100;
        counts[20] = 100;
        counts[30] = 100;
        let map = quantize_quality_scores(&counts, 4, usize::from(MIN_USABLE_Q_SCORE));
        let levels: BTreeSet<u8> = map.into_iter().collect();
        assert!(levels.len() <= 4);
    }
}
