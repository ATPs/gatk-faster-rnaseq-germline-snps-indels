use anyhow::{bail, Context, Result};
use clap::Parser;
use rust_htslib::bam::header::HeaderRecord;
use rust_htslib::bam::record::Aux;
use rust_htslib::{bam, bam::Read};
use std::cmp::Ordering;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

const DUPLICATE_FLAG: u16 = 0x400;
const MIN_BASE_QUAL: u8 = 15;
const SCORE_CAP_PER_END: i32 = i16::MAX as i32 / 2;
const QC_FAIL_PENALTY: i32 = i16::MIN as i32 / 2;
const UNKNOWN_LIBRARY: &str = "Unknown Library";

#[derive(Debug, Parser)]
#[command(
    author,
    version,
    about = "Mark PCR duplicates in a coordinate-sorted BAM"
)]
struct Args {
    #[arg(long)]
    input_bam: PathBuf,

    #[arg(long)]
    output_bam: PathBuf,

    #[arg(long)]
    output_metrics: PathBuf,

    #[arg(
        long,
        default_value = "/data/p/samtools/samtools-1.22.1_installed/bin/samtools"
    )]
    samtools: PathBuf,

    #[arg(long, default_value_t = 1)]
    threads: usize,

    #[arg(long, default_value_t = true, action = clap::ArgAction::Set)]
    create_index: bool,
}

#[derive(Clone, Debug)]
struct HeaderLibraries {
    by_read_group: HashMap<Vec<u8>, String>,
    library_ids: HashMap<String, u16>,
}

#[derive(Clone, Debug)]
struct Template {
    reads: Vec<ReadSummary>,
}

#[derive(Clone, Debug)]
struct ReadSummary {
    tid: i32,
    pos: i64,
    reverse: bool,
    unmapped: bool,
    has_mapped_mate: bool,
    secondary_or_supplementary: bool,
    library_id: u16,
    score: i32,
    leading_clip: i64,
    trailing_clip: i64,
    reference_len: i64,
}

#[derive(Clone, Debug, Eq, Hash, PartialEq)]
enum DuplicateKey {
    Fragment(EndKeyWithLibrary),
    Pair(EndKeyWithLibrary, EndKey),
}

#[derive(Clone, Debug, Eq, Hash, PartialEq)]
struct EndKeyWithLibrary {
    tid: i32,
    stranded_unclipped_start: i64,
    reverse: bool,
    library_id: u16,
}

#[derive(Clone, Debug, Eq, Hash, PartialEq)]
struct EndKey {
    tid: i32,
    stranded_unclipped_start: i64,
    reverse: bool,
}

#[derive(Clone, Debug)]
struct Candidate {
    name: Vec<u8>,
    score: i32,
}

#[derive(Clone, Debug)]
struct GroupState {
    best: Candidate,
}

#[derive(Clone, Debug, Default)]
struct Metrics {
    unpaired_reads_examined: u64,
    read_pairs_examined_raw: u64,
    secondary_or_supplementary_reads: u64,
    unmapped_reads: u64,
    unpaired_read_duplicates: u64,
    read_pair_duplicates_raw: u64,
}

fn main() -> Result<()> {
    let args = Args::parse();
    validate_args(&args)?;
    mark_duplicates(&args)
}

fn validate_args(args: &Args) -> Result<()> {
    if !args.input_bam.exists() {
        bail!("input BAM not found: {}", args.input_bam.display());
    }
    if !args.samtools.exists() {
        bail!("samtools not found: {}", args.samtools.display());
    }
    if args.input_bam == args.output_bam {
        bail!("--input-bam and --output-bam must be different paths");
    }
    if args.threads == 0 {
        bail!("--threads must be at least 1");
    }
    Ok(())
}

fn mark_duplicates(args: &Args) -> Result<()> {
    create_parent_dir(&args.output_bam)?;
    create_parent_dir(&args.output_metrics)?;

    let mut reader = bam::Reader::from_path(&args.input_bam)
        .with_context(|| format!("failed to open {}", args.input_bam.display()))?;
    if args.threads > 1 {
        reader.set_threads(args.threads).with_context(|| {
            format!(
                "failed to set reader threads for {}",
                args.input_bam.display()
            )
        })?;
    }
    let header_view = reader.header().clone();
    let header_libraries = HeaderLibraries::from_header(&header_view);
    let output_header = output_header(&header_view);
    let duplicate_names = collect_duplicate_names(reader, &header_libraries)?;

    let mut reader = bam::Reader::from_path(&args.input_bam)
        .with_context(|| format!("failed to reopen {}", args.input_bam.display()))?;
    if args.threads > 1 {
        reader.set_threads(args.threads).with_context(|| {
            format!(
                "failed to set reader threads for {}",
                args.input_bam.display()
            )
        })?;
    }
    let mut writer = bam::Writer::from_path(&args.output_bam, &output_header, bam::Format::Bam)
        .with_context(|| format!("failed to create {}", args.output_bam.display()))?;
    if args.threads > 1 {
        writer.set_threads(args.threads).with_context(|| {
            format!(
                "failed to set writer threads for {}",
                args.output_bam.display()
            )
        })?;
    }

    let mut metrics_by_library = header_libraries.empty_metrics_by_library();
    for record_result in reader.records() {
        let mut record = record_result
            .with_context(|| format!("failed to read {}", args.input_bam.display()))?;
        let should_mark_duplicate = !record.is_unmapped()
            && !record.is_secondary()
            && !record.is_supplementary()
            && duplicate_names.contains(record.qname());
        if should_mark_duplicate {
            record.set_flags(record.flags() | DUPLICATE_FLAG);
        } else {
            record.set_flags(record.flags() & !DUPLICATE_FLAG);
        }
        update_metrics(&mut metrics_by_library, &header_libraries, &record);
        writer
            .write(&record)
            .with_context(|| format!("failed to write {}", args.output_bam.display()))?;
    }
    drop(writer);

    write_metrics(&args.output_metrics, &metrics_by_library)?;
    if args.create_index {
        index_bam(
            &args.samtools,
            &args.output_bam,
            &bam_index_path(&args.output_bam),
            args.threads,
        )?;
    }
    Ok(())
}

fn collect_duplicate_names(
    mut reader: bam::Reader,
    header_libraries: &HeaderLibraries,
) -> Result<HashSet<Vec<u8>>> {
    let mut templates: HashMap<Vec<u8>, Template> = HashMap::new();
    for record_result in reader.records() {
        let record = record_result.context("failed to read BAM record")?;
        if record.is_secondary() || record.is_supplementary() || record.is_unmapped() {
            continue;
        }
        let summary = summarize_record(&record, header_libraries);
        templates
            .entry(record.qname().to_vec())
            .or_insert_with(|| Template {
                reads: Vec::with_capacity(2),
            })
            .reads
            .push(summary);
    }

    let mut group_states: HashMap<DuplicateKey, GroupState> = HashMap::new();
    let mut duplicate_names = HashSet::new();
    for (name, template) in templates {
        for (key, candidate) in candidates_for_template(&name, &template)? {
            update_group_state(&mut group_states, &mut duplicate_names, key, candidate);
        }
    }
    Ok(duplicate_names)
}

fn candidates_for_template(
    name: &[u8],
    template: &Template,
) -> Result<Vec<(DuplicateKey, Candidate)>> {
    let primary_reads: Vec<&ReadSummary> = template
        .reads
        .iter()
        .filter(|read| !read.unmapped && !read.secondary_or_supplementary)
        .collect();
    if primary_reads.is_empty() {
        return Ok(Vec::new());
    }
    if primary_reads.len() > 2 {
        bail!(
            "read '{}' has more than two primary mapped records; v1 only supports singleton fragments and pairs",
            String::from_utf8_lossy(name)
        );
    }

    let mut out = Vec::new();
    for read in primary_reads.iter().filter(|read| !read.has_mapped_mate) {
        out.push((
            DuplicateKey::Fragment(read.end_key_with_library()),
            Candidate {
                name: name.to_vec(),
                score: read.score,
            },
        ));
    }

    let paired: Vec<&ReadSummary> = primary_reads
        .iter()
        .copied()
        .filter(|read| read.has_mapped_mate)
        .collect();
    if paired.len() == 2 {
        let (first, second) = order_pair(paired[0], paired[1]);
        out.push((
            DuplicateKey::Pair(first.end_key_with_library(), second.end_key()),
            Candidate {
                name: name.to_vec(),
                score: first.score + second.score,
            },
        ));
    }
    Ok(out)
}

fn update_group_state(
    group_states: &mut HashMap<DuplicateKey, GroupState>,
    duplicate_names: &mut HashSet<Vec<u8>>,
    key: DuplicateKey,
    candidate: Candidate,
) {
    match group_states.get_mut(&key) {
        None => {
            group_states.insert(key, GroupState { best: candidate });
        }
        Some(state) => {
            if better_candidate(&candidate, &state.best) {
                duplicate_names.insert(state.best.name.clone());
                state.best = candidate;
            } else {
                duplicate_names.insert(candidate.name);
            }
        }
    }
}

fn better_candidate(left: &Candidate, right: &Candidate) -> bool {
    match left.score.cmp(&right.score) {
        Ordering::Greater => true,
        Ordering::Less => false,
        Ordering::Equal => left.name < right.name,
    }
}

fn order_pair<'a>(
    left: &'a ReadSummary,
    right: &'a ReadSummary,
) -> (&'a ReadSummary, &'a ReadSummary) {
    let left_key = (left.tid, left.stranded_unclipped_start());
    let right_key = (right.tid, right.stranded_unclipped_start());
    let (mut first, mut second) = if left_key <= right_key {
        (left, right)
    } else {
        (right, left)
    };
    if first.tid == second.tid
        && first.stranded_unclipped_start() == second.stranded_unclipped_start()
        && first.reverse
        && !second.reverse
    {
        std::mem::swap(&mut first, &mut second);
    }
    (first, second)
}

fn summarize_record(record: &bam::Record, header_libraries: &HeaderLibraries) -> ReadSummary {
    let (_, library_id) = header_libraries.library_for_record(record);
    let cigar = record.cigar();
    ReadSummary {
        tid: record.tid(),
        pos: record.pos(),
        reverse: record.is_reverse(),
        unmapped: record.is_unmapped(),
        has_mapped_mate: record.is_paired() && !record.is_mate_unmapped(),
        secondary_or_supplementary: record.is_secondary() || record.is_supplementary(),
        library_id,
        score: duplicate_score(record),
        leading_clip: cigar.leading_softclips() + cigar.leading_hardclips(),
        trailing_clip: cigar.trailing_softclips() + cigar.trailing_hardclips(),
        reference_len: cigar.end_pos() - record.pos(),
    }
}

impl ReadSummary {
    fn unclipped_start(&self) -> i64 {
        self.pos - self.leading_clip + 1
    }

    fn unclipped_end(&self) -> i64 {
        self.pos + self.reference_len + self.trailing_clip
    }

    fn stranded_unclipped_start(&self) -> i64 {
        if self.reverse {
            self.unclipped_end()
        } else {
            self.unclipped_start()
        }
    }

    fn end_key_with_library(&self) -> EndKeyWithLibrary {
        EndKeyWithLibrary {
            tid: self.tid,
            stranded_unclipped_start: self.stranded_unclipped_start(),
            reverse: self.reverse,
            library_id: self.library_id,
        }
    }

    fn end_key(&self) -> EndKey {
        EndKey {
            tid: self.tid,
            stranded_unclipped_start: self.stranded_unclipped_start(),
            reverse: self.reverse,
        }
    }
}

fn duplicate_score(record: &bam::Record) -> i32 {
    let quality_sum: i32 = record
        .qual()
        .iter()
        .filter(|quality| **quality >= MIN_BASE_QUAL)
        .map(|quality| i32::from(*quality))
        .sum::<i32>()
        .min(SCORE_CAP_PER_END);
    quality_sum
        + if record.is_quality_check_failed() {
            QC_FAIL_PENALTY
        } else {
            0
        }
}

impl HeaderLibraries {
    fn from_header(header: &bam::HeaderView) -> Self {
        let header_text = String::from_utf8_lossy(header.as_bytes());
        let mut by_read_group = HashMap::new();
        let mut library_ids = HashMap::new();
        for line in header_text.lines().filter(|line| line.starts_with("@RG\t")) {
            let mut id: Option<Vec<u8>> = None;
            let mut library: Option<String> = None;
            for field in line.split('\t').skip(1) {
                if let Some(value) = field.strip_prefix("ID:") {
                    id = Some(value.as_bytes().to_vec());
                } else if let Some(value) = field.strip_prefix("LB:") {
                    library = Some(value.to_string());
                }
            }
            if let Some(id) = id {
                let library = library.unwrap_or_else(|| UNKNOWN_LIBRARY.to_string());
                by_read_group.insert(id, library.clone());
                next_library_id(&mut library_ids, &library);
            }
        }
        next_library_id(&mut library_ids, UNKNOWN_LIBRARY);
        Self {
            by_read_group,
            library_ids,
        }
    }

    fn library_for_record(&self, record: &bam::Record) -> (String, u16) {
        let library = match record.aux(b"RG") {
            Ok(Aux::String(read_group)) => self
                .by_read_group
                .get(read_group.as_bytes())
                .cloned()
                .unwrap_or_else(|| UNKNOWN_LIBRARY.to_string()),
            _ => UNKNOWN_LIBRARY.to_string(),
        };
        let library_id = *self
            .library_ids
            .get(&library)
            .or_else(|| self.library_ids.get(UNKNOWN_LIBRARY))
            .unwrap_or(&0);
        (library, library_id)
    }

    fn empty_metrics_by_library(&self) -> BTreeMap<String, Metrics> {
        self.by_read_group
            .values()
            .map(|library| (library.clone(), Metrics::default()))
            .collect()
    }
}

fn next_library_id(library_ids: &mut HashMap<String, u16>, library: &str) -> u16 {
    if let Some(id) = library_ids.get(library) {
        return *id;
    }
    let id = u16::try_from(library_ids.len()).unwrap_or(u16::MAX);
    library_ids.insert(library.to_string(), id);
    id
}

fn output_header(input_header: &bam::HeaderView) -> bam::Header {
    let mut header = bam::Header::from_template(input_header);
    let mut pg = HeaderRecord::new(b"PG");
    pg.push_tag(b"ID", "rust_mark_duplicates")
        .push_tag(b"PN", "rust_mark_duplicates")
        .push_tag(b"VN", env!("CARGO_PKG_VERSION"));
    header.push_record(&pg);
    header
}

fn update_metrics(
    metrics_by_library: &mut BTreeMap<String, Metrics>,
    header_libraries: &HeaderLibraries,
    record: &bam::Record,
) {
    let (library, _) = header_libraries.library_for_record(record);
    let metrics = metrics_by_library.entry(library).or_default();
    let secondary_or_supplementary = record.is_secondary() || record.is_supplementary();
    let has_mapped_mate = record.is_paired() && !record.is_mate_unmapped();

    if record.is_unmapped() {
        metrics.unmapped_reads += 1;
    } else if secondary_or_supplementary {
        metrics.secondary_or_supplementary_reads += 1;
    } else if !has_mapped_mate {
        metrics.unpaired_reads_examined += 1;
    } else {
        metrics.read_pairs_examined_raw += 1;
    }

    if record.is_duplicate() && !secondary_or_supplementary && !record.is_unmapped() {
        if !has_mapped_mate {
            metrics.unpaired_read_duplicates += 1;
        } else {
            metrics.read_pair_duplicates_raw += 1;
        }
    }
}

fn write_metrics(path: &Path, metrics_by_library: &BTreeMap<String, Metrics>) -> Result<()> {
    let file =
        File::create(path).with_context(|| format!("failed to create {}", path.display()))?;
    let mut writer = BufWriter::new(file);
    writeln!(writer, "## METRICS CLASS\tpicard.sam.DuplicationMetrics")?;
    writeln!(
        writer,
        "LIBRARY\tUNPAIRED_READS_EXAMINED\tREAD_PAIRS_EXAMINED\tSECONDARY_OR_SUPPLEMENTARY_RDS\tUNMAPPED_READS\tUNPAIRED_READ_DUPLICATES\tREAD_PAIR_DUPLICATES\tREAD_PAIR_OPTICAL_DUPLICATES\tPERCENT_DUPLICATION\tESTIMATED_LIBRARY_SIZE"
    )?;
    for (library, metrics) in metrics_by_library {
        let read_pairs_examined = metrics.read_pairs_examined_raw / 2;
        let read_pair_duplicates = metrics.read_pair_duplicates_raw / 2;
        let denominator = metrics.unpaired_reads_examined + 2 * read_pairs_examined;
        let numerator = metrics.unpaired_read_duplicates + 2 * read_pair_duplicates;
        let percent_duplication = if denominator == 0 {
            0.0
        } else {
            numerator as f64 / denominator as f64
        };
        let estimated_library_size = estimate_library_size(
            read_pairs_examined,
            read_pairs_examined - read_pair_duplicates,
        )
        .map(|value| value.to_string())
        .unwrap_or_default();
        writeln!(
            writer,
            "{library}\t{}\t{}\t{}\t{}\t{}\t{}\t0\t{:.6}\t{}",
            metrics.unpaired_reads_examined,
            read_pairs_examined,
            metrics.secondary_or_supplementary_reads,
            metrics.unmapped_reads,
            metrics.unpaired_read_duplicates,
            read_pair_duplicates,
            percent_duplication,
            estimated_library_size
        )?;
    }
    Ok(())
}

fn estimate_library_size(read_pairs: u64, unique_read_pairs: u64) -> Option<u64> {
    if read_pairs == 0 || unique_read_pairs == 0 || unique_read_pairs >= read_pairs {
        return None;
    }

    let c = unique_read_pairs as f64;
    let n = read_pairs as f64;
    let mut lower = 1.0;
    let mut upper = 100.0;
    if library_size_function(lower * c, c, n) < 0.0 {
        return None;
    }
    while library_size_function(upper * c, c, n) > 0.0 {
        upper *= 10.0;
    }
    for _ in 0..40 {
        let mid = (lower + upper) / 2.0;
        let value = library_size_function(mid * c, c, n);
        if value == 0.0 {
            lower = mid;
            upper = mid;
            break;
        }
        if value > 0.0 {
            lower = mid;
        } else {
            upper = mid;
        }
    }
    Some((unique_read_pairs as f64 * (lower + upper) / 2.0) as u64)
}

fn library_size_function(x: f64, c: f64, n: f64) -> f64 {
    c / x - 1.0 + (-n / x).exp()
}

fn index_bam(samtools: &Path, bam: &Path, bai: &Path, threads: usize) -> Result<()> {
    if bai.exists() {
        fs::remove_file(bai).with_context(|| format!("failed to remove {}", bai.display()))?;
    }
    let status = Command::new(samtools)
        .arg("index")
        .arg("-@")
        .arg(threads.saturating_sub(1).to_string())
        .arg("-o")
        .arg(bai)
        .arg(bam)
        .status()
        .with_context(|| format!("failed to run {} index", samtools.display()))?;
    if !status.success() {
        bail!("{} index failed with status {status}", samtools.display());
    }
    Ok(())
}

fn bam_index_path(bam: &Path) -> PathBuf {
    if bam.extension().is_some_and(|ext| ext == "bam") {
        bam.with_extension("bai")
    } else {
        PathBuf::from(format!("{}.bai", bam.display()))
    }
}

fn create_parent_dir(path: &Path) -> Result<()> {
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            fs::create_dir_all(parent)
                .with_context(|| format!("failed to create directory {}", parent.display()))?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn read(tid: i32, pos: i64, reverse: bool, score: i32) -> ReadSummary {
        ReadSummary {
            tid,
            pos,
            reverse,
            unmapped: false,
            has_mapped_mate: true,
            secondary_or_supplementary: false,
            library_id: 0,
            score,
            leading_clip: 0,
            trailing_clip: 0,
            reference_len: 100,
        }
    }

    #[test]
    fn paired_duplicate_group_marks_lower_scoring_template() {
        let first = Template {
            reads: vec![read(0, 99, false, 3000), read(0, 499, true, 3000)],
        };
        let second = Template {
            reads: vec![read(0, 99, false, 1000), read(0, 499, true, 1000)],
        };
        let mut states = HashMap::new();
        let mut duplicates = HashSet::new();
        for (name, template) in [
            (b"winner".as_slice(), &first),
            (b"loser".as_slice(), &second),
        ] {
            for (key, candidate) in candidates_for_template(name, template).unwrap() {
                update_group_state(&mut states, &mut duplicates, key, candidate);
            }
        }
        assert!(duplicates.contains(b"loser".as_slice()));
        assert!(!duplicates.contains(b"winner".as_slice()));
    }

    #[test]
    fn better_late_candidate_replaces_previous_best() {
        let first = Template {
            reads: vec![read(0, 99, false, 1000), read(0, 499, true, 1000)],
        };
        let second = Template {
            reads: vec![read(0, 99, false, 3000), read(0, 499, true, 3000)],
        };
        let mut states = HashMap::new();
        let mut duplicates = HashSet::new();
        for (name, template) in [
            (b"old_best".as_slice(), &first),
            (b"new_best".as_slice(), &second),
        ] {
            for (key, candidate) in candidates_for_template(name, template).unwrap() {
                update_group_state(&mut states, &mut duplicates, key, candidate);
            }
        }
        assert!(duplicates.contains(b"old_best".as_slice()));
        assert!(!duplicates.contains(b"new_best".as_slice()));
    }

    #[test]
    fn reverse_strand_uses_unclipped_end_as_stranded_start() {
        let mut reverse = read(0, 99, true, 1000);
        reverse.leading_clip = 5;
        reverse.trailing_clip = 7;
        reverse.reference_len = 50;
        assert_eq!(reverse.unclipped_start(), 95);
        assert_eq!(reverse.unclipped_end(), 156);
        assert_eq!(reverse.stranded_unclipped_start(), 156);
    }

    #[test]
    fn estimated_library_size_matches_picard() {
        assert_eq!(
            estimate_library_size(7_360_339, 7_360_339 - 602_223),
            Some(42_490_397)
        );
        assert_eq!(estimate_library_size(10, 10), None);
    }
}
