use anyhow::{bail, Context, Result};
use clap::{Parser, ValueEnum};
use rust_htslib::bam::pileup::Indel;
use rust_htslib::{bam, bam::Read, faidx};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

#[derive(Debug, Parser)]
#[command(
    author,
    version,
    about = "Build HaplotypeCaller candidate intervals from BAM pileup evidence"
)]
struct Args {
    #[command(subcommand)]
    command: Command,
}

#[derive(Debug, clap::Subcommand)]
enum Command {
    CandidateIntervals(CandidateIntervalsArgs),
}

#[derive(Copy, Clone, Debug, Eq, PartialEq, ValueEnum)]
enum EmptyBehavior {
    Error,
    Input,
}

#[derive(Debug, Parser)]
struct CandidateIntervalsArgs {
    #[arg(long)]
    input_bam: PathBuf,

    #[arg(long = "ref")]
    reference: PathBuf,

    #[arg(long)]
    input_interval_list: PathBuf,

    #[arg(long)]
    output_interval_list: PathBuf,

    #[arg(long)]
    output_summary: Option<PathBuf>,

    #[arg(long)]
    output_bed: Option<PathBuf>,

    #[arg(long, default_value_t = 20)]
    min_mapq: u8,

    #[arg(long, default_value_t = 10)]
    min_baseq: u8,

    #[arg(long, default_value_t = 1)]
    min_alt_count: u32,

    #[arg(long, default_value_t = 1)]
    min_indel_count: u32,

    #[arg(long, default_value_t = 0.0)]
    min_alt_fraction: f64,

    #[arg(long, default_value_t = 150)]
    padding: u64,

    #[arg(long, default_value_t = 100_000)]
    max_depth: u32,

    #[arg(long, default_value_t = 1)]
    threads: usize,

    #[arg(long, default_value_t = false)]
    exclude_supplementary: bool,

    #[arg(long, value_enum, default_value_t = EmptyBehavior::Input)]
    empty_behavior: EmptyBehavior,
}

#[derive(Clone, Debug)]
struct DictRecord {
    length: u64,
}

#[derive(Clone, Debug)]
struct SequenceDict {
    header_lines: Vec<String>,
    records: Vec<DictRecord>,
    index_by_name: HashMap<String, usize>,
}

impl SequenceDict {
    fn order(&self, contig: &str) -> Option<usize> {
        self.index_by_name.get(contig).copied()
    }

    fn contig_length(&self, contig: &str) -> Option<u64> {
        self.order(contig).map(|idx| self.records[idx].length)
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct Interval {
    contig: String,
    start: u64,
    end: u64,
}

impl Interval {
    fn len(&self) -> u64 {
        self.end - self.start + 1
    }
}

#[derive(Default)]
struct Summary {
    input_intervals: u64,
    input_bases: u64,
    pileup_sites: u64,
    covered_sites: u64,
    candidate_sites: u64,
    total_depth: u128,
    max_depth: u32,
    output_intervals: u64,
    output_bases: u64,
    used_input_fallback: bool,
}

#[derive(Copy, Clone)]
struct Evidence {
    depth: u32,
    alt_count: u32,
    indel_count: u32,
}

pub fn run_cli() -> Result<()> {
    match Args::parse().command {
        Command::CandidateIntervals(args) => run_candidate_intervals(&args),
    }
}

fn run_candidate_intervals(args: &CandidateIntervalsArgs) -> Result<()> {
    validate_args(args)?;
    let (dict, input_intervals) = read_interval_list(&args.input_interval_list)?;
    let mut bam = bam::IndexedReader::from_path(&args.input_bam)
        .with_context(|| format!("failed to open indexed BAM {}", args.input_bam.display()))?;
    if args.threads > 1 {
        bam.set_threads(args.threads).with_context(|| {
            format!(
                "failed to set BAM reader threads for {}",
                args.input_bam.display()
            )
        })?;
    }
    let bam_tid_by_name = bam_tid_by_name(bam.header())?;
    let reference = faidx::Reader::from_path(&args.reference).with_context(|| {
        format!(
            "failed to open FASTA index for {}",
            args.reference.display()
        )
    })?;

    let mut summary = Summary {
        input_intervals: input_intervals.len() as u64,
        input_bases: input_intervals.iter().map(Interval::len).sum(),
        ..Summary::default()
    };
    let mut candidates = Vec::new();

    for interval in &input_intervals {
        let tid = *bam_tid_by_name.get(&interval.contig).with_context(|| {
            format!(
                "contig '{}' from {} is not present in BAM header",
                interval.contig,
                args.input_interval_list.display()
            )
        })?;
        let contig_len = dict.contig_length(&interval.contig).with_context(|| {
            format!(
                "contig '{}' from interval body is missing from interval_list header",
                interval.contig
            )
        })?;
        let ref_len = reference.fetch_seq_len(&interval.contig);
        if ref_len == 0 {
            bail!(
                "contig '{}' is not present in reference FASTA {}",
                interval.contig,
                args.reference.display()
            );
        }
        if interval.end > ref_len {
            bail!(
                "interval {}:{}-{} extends past FASTA contig length {}",
                interval.contig,
                interval.start,
                interval.end,
                ref_len
            );
        }

        let ref_bases = reference
            .fetch_seq(
                &interval.contig,
                (interval.start - 1) as usize,
                (interval.end - 1) as usize,
            )
            .with_context(|| {
                format!(
                    "failed to fetch reference sequence {}:{}-{}",
                    interval.contig, interval.start, interval.end
                )
            })?;

        bam.fetch((tid as i32, (interval.start - 1) as i64, interval.end as i64))
            .with_context(|| {
                format!(
                    "failed to fetch BAM region {}:{}-{}",
                    interval.contig, interval.start, interval.end
                )
            })?;

        let mut pileups = bam.pileup();
        pileups.set_max_depth(args.max_depth);
        for pileup in pileups {
            let pileup = pileup.with_context(|| {
                format!(
                    "failed while reading pileup for {}:{}-{}",
                    interval.contig, interval.start, interval.end
                )
            })?;
            if pileup.tid() != tid {
                continue;
            }
            let pos0 = u64::from(pileup.pos());
            if pos0 < interval.start - 1 || pos0 >= interval.end {
                continue;
            }
            summary.pileup_sites += 1;
            let ref_base = normalize_base(ref_bases[(pos0 - (interval.start - 1)) as usize]);
            let evidence = pileup_evidence(&pileup, ref_base, args);
            if evidence.depth > 0 {
                summary.covered_sites += 1;
                summary.total_depth += u128::from(evidence.depth);
                summary.max_depth = summary.max_depth.max(evidence.depth);
            }
            if is_candidate(evidence, args) {
                summary.candidate_sites += 1;
                candidates.push(padded_interval(
                    &interval.contig,
                    pos0 + 1,
                    args.padding,
                    contig_len,
                ));
            }
        }
    }

    if candidates.is_empty() && args.empty_behavior == EmptyBehavior::Input {
        summary.used_input_fallback = true;
        candidates = input_intervals.clone();
    }
    if candidates.is_empty() {
        bail!("no candidate intervals were found; rerun with --empty-behavior input or looser thresholds");
    }

    sort_and_merge(&mut candidates, &dict)?;
    if !summary.used_input_fallback {
        let mut input_bounds = input_intervals.clone();
        sort_and_merge(&mut input_bounds, &dict)?;
        candidates = intersect_intervals(&candidates, &input_bounds, &dict)?;
        if candidates.is_empty() {
            bail!("candidate intervals were removed after clipping to the input intervals");
        }
        sort_and_merge(&mut candidates, &dict)?;
    }
    summary.output_intervals = candidates.len() as u64;
    summary.output_bases = candidates.iter().map(Interval::len).sum();

    write_interval_list(&args.output_interval_list, &dict, &candidates)?;
    if let Some(path) = &args.output_bed {
        write_bed(path, &candidates)?;
    }
    if let Some(path) = &args.output_summary {
        write_summary(path, args, &summary)?;
    }
    Ok(())
}

fn validate_args(args: &CandidateIntervalsArgs) -> Result<()> {
    if args.min_alt_count == 0 {
        bail!("--min-alt-count must be at least 1");
    }
    if args.min_indel_count == 0 {
        bail!("--min-indel-count must be at least 1");
    }
    if !(0.0..=1.0).contains(&args.min_alt_fraction) {
        bail!("--min-alt-fraction must be between 0 and 1");
    }
    if args.max_depth == 0 {
        bail!("--max-depth must be at least 1");
    }
    if args.threads == 0 {
        bail!("--threads must be at least 1");
    }
    Ok(())
}

fn pileup_evidence(
    pileup: &bam::pileup::Pileup,
    ref_base: u8,
    args: &CandidateIntervalsArgs,
) -> Evidence {
    let mut evidence = Evidence {
        depth: 0,
        alt_count: 0,
        indel_count: 0,
    };
    for alignment in pileup.alignments() {
        let record = alignment.record();
        if !read_is_usable(&record, args) || alignment.is_refskip() {
            continue;
        }

        let base_quality_ok = alignment
            .qpos()
            .and_then(|qpos| record.qual().get(qpos).copied())
            .is_some_and(|quality| quality >= args.min_baseq);

        if let Some(qpos) = alignment.qpos() {
            if base_quality_ok {
                evidence.depth += 1;
                let read_base = normalize_base(record.seq()[qpos]);
                if is_acgt(read_base) && is_acgt(ref_base) && read_base != ref_base {
                    evidence.alt_count += 1;
                }
            }
        } else if alignment.is_del() {
            evidence.depth += 1;
        }

        if alignment.indel() != Indel::None && base_quality_ok {
            evidence.indel_count += 1;
        }
    }
    evidence
}

fn read_is_usable(record: &bam::Record, args: &CandidateIntervalsArgs) -> bool {
    const UNMAPPED: u16 = 0x4;
    const SECONDARY: u16 = 0x100;
    const QCFAIL: u16 = 0x200;
    const DUPLICATE: u16 = 0x400;
    const SUPPLEMENTARY: u16 = 0x800;

    let mut excluded = UNMAPPED | SECONDARY | QCFAIL | DUPLICATE;
    if args.exclude_supplementary {
        excluded |= SUPPLEMENTARY;
    }
    record.flags() & excluded == 0 && record.mapq() >= args.min_mapq
}

fn is_candidate(evidence: Evidence, args: &CandidateIntervalsArgs) -> bool {
    if evidence.indel_count >= args.min_indel_count {
        return true;
    }
    if evidence.depth == 0 || evidence.alt_count < args.min_alt_count {
        return false;
    }
    f64::from(evidence.alt_count) / f64::from(evidence.depth) >= args.min_alt_fraction
}

fn padded_interval(contig: &str, pos1: u64, padding: u64, contig_len: u64) -> Interval {
    Interval {
        contig: contig.to_string(),
        start: pos1.saturating_sub(padding).max(1),
        end: pos1.saturating_add(padding).min(contig_len),
    }
}

fn read_interval_list(path: &Path) -> Result<(SequenceDict, Vec<Interval>)> {
    let file = File::open(path).with_context(|| format!("failed to open {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut header_lines = Vec::new();
    let mut body_lines = Vec::new();
    for line in reader.lines() {
        let line = line.with_context(|| format!("failed to read {}", path.display()))?;
        if line.starts_with('@') {
            header_lines.push(line);
        } else {
            body_lines.push(line);
        }
    }
    let dict = parse_dict_lines(&header_lines, path)?;
    let mut intervals = Vec::new();
    for (line_idx, line) in body_lines.iter().enumerate() {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        let fields: Vec<&str> = trimmed.split_whitespace().collect();
        if fields.len() < 3 {
            bail!(
                "{}: interval row has fewer than 3 columns near body line {}",
                path.display(),
                line_idx + 1
            );
        }
        let contig = fields[0].to_string();
        let start = parse_u64(fields[1], path, line_idx + 1, "interval start")?;
        let end = parse_u64(fields[2], path, line_idx + 1, "interval end")?;
        validate_interval(path, line_idx + 1, &contig, start, end, &dict)?;
        intervals.push(Interval { contig, start, end });
    }
    if intervals.is_empty() {
        bail!("no intervals found in {}", path.display());
    }
    Ok((dict, intervals))
}

fn parse_dict_lines(lines: &[String], path: &Path) -> Result<SequenceDict> {
    let mut records = Vec::new();
    let mut index_by_name = HashMap::new();
    for line in lines {
        if !line.starts_with("@SQ\t") {
            continue;
        }
        let mut name = None;
        let mut length = None;
        for field in line.split('\t').skip(1) {
            if let Some(value) = field.strip_prefix("SN:") {
                name = Some(value.to_string());
            } else if let Some(value) = field.strip_prefix("LN:") {
                length =
                    Some(value.parse::<u64>().with_context(|| {
                        format!("invalid LN value in {}: {value}", path.display())
                    })?);
            }
        }
        let name =
            name.with_context(|| format!("missing SN field in @SQ line in {}", path.display()))?;
        let length = length
            .with_context(|| format!("missing LN field in @SQ line in {}", path.display()))?;
        if length == 0 {
            bail!("zero-length contig '{name}' in {}", path.display());
        }
        if index_by_name.contains_key(&name) {
            bail!("duplicate contig '{name}' in {}", path.display());
        }
        index_by_name.insert(name, records.len());
        records.push(DictRecord { length });
    }
    if records.is_empty() {
        bail!("no @SQ records found in {}", path.display());
    }
    Ok(SequenceDict {
        header_lines: lines.to_vec(),
        records,
        index_by_name,
    })
}

fn parse_u64(value: &str, path: &Path, line_no: usize, label: &str) -> Result<u64> {
    value.parse::<u64>().with_context(|| {
        format!(
            "{}:{line_no}: invalid {label} value '{value}'",
            path.display()
        )
    })
}

fn validate_interval(
    path: &Path,
    line_no: usize,
    contig: &str,
    start: u64,
    end: u64,
    dict: &SequenceDict,
) -> Result<()> {
    if start == 0 {
        bail!(
            "{}:{line_no}: interval start must be at least 1",
            path.display()
        );
    }
    if start > end {
        bail!(
            "{}:{line_no}: interval start is greater than end",
            path.display()
        );
    }
    let length = dict.contig_length(contig).with_context(|| {
        format!(
            "{}:{line_no}: contig '{contig}' is not present in the sequence dictionary",
            path.display()
        )
    })?;
    if end > length {
        bail!(
            "{}:{line_no}: interval {contig}:{start}-{end} extends past contig length {length}",
            path.display()
        );
    }
    Ok(())
}

fn bam_tid_by_name(header: &bam::HeaderView) -> Result<HashMap<String, u32>> {
    let mut tids = HashMap::new();
    for tid in 0..header.target_count() {
        let name = String::from_utf8(header.tid2name(tid).to_vec())
            .with_context(|| format!("BAM header target id {tid} is not valid UTF-8"))?;
        tids.insert(name, tid);
    }
    Ok(tids)
}

fn sort_and_merge(intervals: &mut Vec<Interval>, dict: &SequenceDict) -> Result<()> {
    for interval in intervals.iter() {
        if dict.order(&interval.contig).is_none() {
            bail!(
                "contig '{}' is not present in the sequence dictionary",
                interval.contig
            );
        }
    }
    intervals.sort_by(|a, b| {
        dict.order(&a.contig)
            .cmp(&dict.order(&b.contig))
            .then(a.start.cmp(&b.start))
            .then(a.end.cmp(&b.end))
    });

    let mut merged: Vec<Interval> = Vec::with_capacity(intervals.len());
    for interval in intervals.drain(..) {
        if let Some(current) = merged.last_mut() {
            if current.contig == interval.contig && interval.start <= current.end.saturating_add(1)
            {
                current.end = current.end.max(interval.end);
                continue;
            }
        }
        merged.push(interval);
    }
    *intervals = merged;
    Ok(())
}

fn intersect_intervals(
    candidates: &[Interval],
    allowed: &[Interval],
    dict: &SequenceDict,
) -> Result<Vec<Interval>> {
    let mut intersections = Vec::new();
    let mut candidate_idx = 0;
    let mut allowed_idx = 0;

    while candidate_idx < candidates.len() && allowed_idx < allowed.len() {
        let candidate = &candidates[candidate_idx];
        let allowed_interval = &allowed[allowed_idx];
        let candidate_order = dict.order(&candidate.contig).with_context(|| {
            format!(
                "contig '{}' is not present in the sequence dictionary",
                candidate.contig
            )
        })?;
        let allowed_order = dict.order(&allowed_interval.contig).with_context(|| {
            format!(
                "contig '{}' is not present in the sequence dictionary",
                allowed_interval.contig
            )
        })?;

        if candidate_order < allowed_order
            || (candidate_order == allowed_order && candidate.end < allowed_interval.start)
        {
            candidate_idx += 1;
            continue;
        }
        if allowed_order < candidate_order
            || (candidate_order == allowed_order && allowed_interval.end < candidate.start)
        {
            allowed_idx += 1;
            continue;
        }

        let start = candidate.start.max(allowed_interval.start);
        let end = candidate.end.min(allowed_interval.end);
        if start <= end {
            intersections.push(Interval {
                contig: candidate.contig.clone(),
                start,
                end,
            });
        }

        if candidate.end < allowed_interval.end {
            candidate_idx += 1;
        } else {
            allowed_idx += 1;
        }
    }

    Ok(intersections)
}

fn write_interval_list(path: &Path, dict: &SequenceDict, intervals: &[Interval]) -> Result<()> {
    create_parent_dir(path)?;
    let file =
        File::create(path).with_context(|| format!("failed to create {}", path.display()))?;
    let mut writer = BufWriter::new(file);
    for line in &dict.header_lines {
        writeln!(writer, "{line}")
            .with_context(|| format!("failed to write {}", path.display()))?;
    }
    for interval in intervals {
        writeln!(
            writer,
            "{}\t{}\t{}\t+\t.",
            interval.contig, interval.start, interval.end
        )
        .with_context(|| format!("failed to write {}", path.display()))?;
    }
    Ok(())
}

fn write_bed(path: &Path, intervals: &[Interval]) -> Result<()> {
    create_parent_dir(path)?;
    let file =
        File::create(path).with_context(|| format!("failed to create {}", path.display()))?;
    let mut writer = BufWriter::new(file);
    for interval in intervals {
        writeln!(
            writer,
            "{}\t{}\t{}",
            interval.contig,
            interval.start - 1,
            interval.end
        )
        .with_context(|| format!("failed to write {}", path.display()))?;
    }
    Ok(())
}

fn write_summary(path: &Path, args: &CandidateIntervalsArgs, summary: &Summary) -> Result<()> {
    create_parent_dir(path)?;
    let file =
        File::create(path).with_context(|| format!("failed to create {}", path.display()))?;
    let mut writer = BufWriter::new(file);
    let mean_depth = if summary.covered_sites == 0 {
        0.0
    } else {
        summary.total_depth as f64 / summary.covered_sites as f64
    };
    let reduction_fraction = if summary.input_bases == 0 {
        0.0
    } else {
        1.0 - (summary.output_bases as f64 / summary.input_bases as f64)
    };
    writeln!(writer, "metric\tvalue")?;
    writeln!(writer, "input_bam\t{}", args.input_bam.display())?;
    writeln!(writer, "reference\t{}", args.reference.display())?;
    writeln!(
        writer,
        "input_interval_list\t{}",
        args.input_interval_list.display()
    )?;
    writeln!(writer, "min_mapq\t{}", args.min_mapq)?;
    writeln!(writer, "min_baseq\t{}", args.min_baseq)?;
    writeln!(writer, "min_alt_count\t{}", args.min_alt_count)?;
    writeln!(writer, "min_indel_count\t{}", args.min_indel_count)?;
    writeln!(writer, "min_alt_fraction\t{:.6}", args.min_alt_fraction)?;
    writeln!(writer, "padding\t{}", args.padding)?;
    writeln!(writer, "max_depth_setting\t{}", args.max_depth)?;
    writeln!(
        writer,
        "exclude_supplementary\t{}",
        args.exclude_supplementary
    )?;
    writeln!(writer, "clip_to_input_intervals\ttrue")?;
    writeln!(writer, "input_intervals\t{}", summary.input_intervals)?;
    writeln!(writer, "input_bases\t{}", summary.input_bases)?;
    writeln!(writer, "pileup_sites\t{}", summary.pileup_sites)?;
    writeln!(writer, "covered_sites\t{}", summary.covered_sites)?;
    writeln!(writer, "candidate_sites\t{}", summary.candidate_sites)?;
    writeln!(writer, "mean_filtered_depth\t{mean_depth:.3}")?;
    writeln!(writer, "max_filtered_depth\t{}", summary.max_depth)?;
    writeln!(writer, "output_intervals\t{}", summary.output_intervals)?;
    writeln!(writer, "output_bases\t{}", summary.output_bases)?;
    writeln!(writer, "reduction_fraction\t{reduction_fraction:.6}")?;
    writeln!(
        writer,
        "used_input_fallback\t{}",
        summary.used_input_fallback
    )?;
    Ok(())
}

fn create_parent_dir(path: &Path) -> Result<()> {
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent)
                .with_context(|| format!("failed to create directory {}", parent.display()))?;
        }
    }
    Ok(())
}

fn normalize_base(base: u8) -> u8 {
    base.to_ascii_uppercase()
}

fn is_acgt(base: u8) -> bool {
    matches!(base, b'A' | b'C' | b'G' | b'T')
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_dict() -> SequenceDict {
        let lines = vec![
            "@HD\tVN:1.6\tSO:coordinate".to_string(),
            "@SQ\tSN:chr2\tLN:100".to_string(),
            "@SQ\tSN:chr1\tLN:200".to_string(),
        ];
        parse_dict_lines(&lines, Path::new("test.interval_list")).unwrap()
    }

    fn test_args() -> CandidateIntervalsArgs {
        CandidateIntervalsArgs {
            input_bam: PathBuf::from("in.bam"),
            reference: PathBuf::from("ref.fa"),
            input_interval_list: PathBuf::from("in.interval_list"),
            output_interval_list: PathBuf::from("out.interval_list"),
            output_summary: None,
            output_bed: None,
            min_mapq: 20,
            min_baseq: 10,
            min_alt_count: 2,
            min_indel_count: 1,
            min_alt_fraction: 0.2,
            padding: 150,
            max_depth: 100_000,
            threads: 1,
            exclude_supplementary: false,
            empty_behavior: EmptyBehavior::Input,
        }
    }

    #[test]
    fn padded_interval_clamps_to_contig_bounds() {
        assert_eq!(
            padded_interval("chr1", 3, 10, 20),
            Interval {
                contig: "chr1".to_string(),
                start: 1,
                end: 13,
            }
        );
        assert_eq!(
            padded_interval("chr1", 18, 10, 20),
            Interval {
                contig: "chr1".to_string(),
                start: 8,
                end: 20,
            }
        );
    }

    #[test]
    fn intervals_sort_and_merge_by_dictionary_order() {
        let dict = test_dict();
        let mut intervals = vec![
            Interval {
                contig: "chr1".to_string(),
                start: 50,
                end: 80,
            },
            Interval {
                contig: "chr2".to_string(),
                start: 5,
                end: 10,
            },
            Interval {
                contig: "chr1".to_string(),
                start: 81,
                end: 90,
            },
        ];
        sort_and_merge(&mut intervals, &dict).unwrap();
        assert_eq!(
            intervals,
            vec![
                Interval {
                    contig: "chr2".to_string(),
                    start: 5,
                    end: 10,
                },
                Interval {
                    contig: "chr1".to_string(),
                    start: 50,
                    end: 90,
                },
            ]
        );
    }

    #[test]
    fn candidate_intervals_are_clipped_to_input_intervals() {
        let dict = test_dict();
        let candidates = vec![
            Interval {
                contig: "chr2".to_string(),
                start: 1,
                end: 20,
            },
            Interval {
                contig: "chr1".to_string(),
                start: 40,
                end: 100,
            },
        ];
        let allowed = vec![
            Interval {
                contig: "chr2".to_string(),
                start: 5,
                end: 10,
            },
            Interval {
                contig: "chr1".to_string(),
                start: 50,
                end: 60,
            },
            Interval {
                contig: "chr1".to_string(),
                start: 80,
                end: 90,
            },
        ];

        let clipped = intersect_intervals(&candidates, &allowed, &dict).unwrap();

        assert_eq!(
            clipped,
            vec![
                Interval {
                    contig: "chr2".to_string(),
                    start: 5,
                    end: 10,
                },
                Interval {
                    contig: "chr1".to_string(),
                    start: 50,
                    end: 60,
                },
                Interval {
                    contig: "chr1".to_string(),
                    start: 80,
                    end: 90,
                },
            ]
        );
    }

    #[test]
    fn candidate_requires_alt_threshold_or_indel_threshold() {
        let args = test_args();
        assert!(!is_candidate(
            Evidence {
                depth: 10,
                alt_count: 1,
                indel_count: 0,
            },
            &args
        ));
        assert!(is_candidate(
            Evidence {
                depth: 10,
                alt_count: 2,
                indel_count: 0,
            },
            &args
        ));
        assert!(is_candidate(
            Evidence {
                depth: 10,
                alt_count: 0,
                indel_count: 1,
            },
            &args
        ));
    }
}
