use anyhow::{bail, Context, Result};
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use rust_htslib::bam::pileup::Indel;
use rust_htslib::bam::record::Cigar;
use rust_htslib::{bam, bam::Read, faidx};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

#[derive(Debug, Clone)]
pub struct ActiveRegionDiscoveryConfig {
    pub input_bam: PathBuf,
    pub reference: PathBuf,
    pub input_interval_list: PathBuf,
    pub output_active_bed: PathBuf,
    pub output_summary: PathBuf,
    pub min_mapq: u8,
    pub min_baseq: u8,
    pub min_alt_count: u32,
    pub min_indel_count: u32,
    pub min_alt_fraction: f64,
    pub active_region_padding: u64,
    pub max_depth: u32,
    pub threads: usize,
    pub exclude_supplementary: bool,
}

#[derive(Debug, Clone)]
pub struct HaplotypeCallerConfig {
    pub input_bam: PathBuf,
    pub reference: PathBuf,
    pub input_interval_list: PathBuf,
    pub output_vcf: PathBuf,
    pub dbsnp: Option<PathBuf>,
    pub dont_use_soft_clipped_bases: bool,
    pub standard_min_confidence_threshold_for_calling: f64,
    pub threads: usize,
    pub memory_gb: usize,
    pub native_pair_hmm_threads: usize,
    pub pair_hmm_implementation: String,
}

#[derive(Debug, Clone, Copy, Default, Eq, PartialEq)]
pub struct ActiveRegionSummary {
    pub input_intervals: u64,
    pub input_bases: u64,
    pub pileup_sites: u64,
    pub covered_sites: u64,
    pub active_sites: u64,
    pub active_regions: u64,
    pub active_bases: u64,
    pub max_filtered_depth: u32,
}

pub fn call_variants(config: &HaplotypeCallerConfig) -> Result<()> {
    validate_haplotype_caller_config(config)?;
    bail!(
        "rust_haplotype_caller call has a validated CLI surface, but assembly, PairHMM, genotyping, annotations, and VCF writing are not implemented yet. Use discover-active-regions for the current implemented subsystem."
    );
}

#[derive(Clone, Debug)]
struct DictRecord {
    length: u64,
}

#[derive(Clone, Debug)]
struct SequenceDict {
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

#[derive(Clone, Debug, Eq, PartialEq)]
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

#[derive(Copy, Clone, Debug, Default, Eq, PartialEq)]
struct SiteEvidence {
    depth: u32,
    alt_count: u32,
    indel_count: u32,
}

#[derive(Default)]
struct WorkerOutput {
    summary: ActiveRegionSummary,
    active_regions: Vec<Interval>,
}

pub fn discover_active_regions(
    config: &ActiveRegionDiscoveryConfig,
) -> Result<ActiveRegionSummary> {
    validate_active_region_config(config)?;
    let (dict, intervals) = read_interval_list(&config.input_interval_list)?;
    let partitions = partition_intervals_by_bases(&intervals, config.threads.max(1));

    let thread_pool = ThreadPoolBuilder::new()
        .num_threads(config.threads.max(1))
        .build()
        .context("creating HaplotypeCaller active-region thread pool")?;
    let worker_outputs: Vec<Result<WorkerOutput>> = thread_pool.install(|| {
        partitions
            .into_par_iter()
            .map(|partition| scan_partition(config, &dict, &partition))
            .collect()
    });

    let mut summary = ActiveRegionSummary {
        input_intervals: intervals.len() as u64,
        input_bases: intervals.iter().map(Interval::len).sum(),
        ..ActiveRegionSummary::default()
    };
    let mut active_regions = Vec::new();
    for worker_output in worker_outputs {
        let worker_output = worker_output?;
        summary.pileup_sites += worker_output.summary.pileup_sites;
        summary.covered_sites += worker_output.summary.covered_sites;
        summary.active_sites += worker_output.summary.active_sites;
        summary.max_filtered_depth = summary
            .max_filtered_depth
            .max(worker_output.summary.max_filtered_depth);
        active_regions.extend(worker_output.active_regions);
    }

    sort_and_merge(&mut active_regions, &dict)?;
    summary.active_regions = active_regions.len() as u64;
    summary.active_bases = active_regions.iter().map(Interval::len).sum();

    write_bed(&config.output_active_bed, &active_regions)?;
    write_active_region_summary(&config.output_summary, config, &summary)?;
    Ok(summary)
}

fn validate_active_region_config(config: &ActiveRegionDiscoveryConfig) -> Result<()> {
    if !config.input_bam.exists() {
        bail!("input BAM not found: {}", config.input_bam.display());
    }
    if !config.reference.exists() {
        bail!("reference FASTA not found: {}", config.reference.display());
    }
    if !config.input_interval_list.exists() {
        bail!(
            "input interval_list not found: {}",
            config.input_interval_list.display()
        );
    }
    if config.min_alt_count == 0 {
        bail!("--min-alt-count must be at least 1");
    }
    if config.min_indel_count == 0 {
        bail!("--min-indel-count must be at least 1");
    }
    if !(0.0..=1.0).contains(&config.min_alt_fraction) {
        bail!("--min-alt-fraction must be between 0 and 1");
    }
    if config.max_depth == 0 {
        bail!("--max-depth must be at least 1");
    }
    if config.threads == 0 {
        bail!("--threads must be at least 1");
    }
    Ok(())
}

fn validate_haplotype_caller_config(config: &HaplotypeCallerConfig) -> Result<()> {
    if !config.input_bam.exists() {
        bail!("input BAM not found: {}", config.input_bam.display());
    }
    if !config.reference.exists() {
        bail!("reference FASTA not found: {}", config.reference.display());
    }
    if !config.input_interval_list.exists() {
        bail!(
            "input interval_list not found: {}",
            config.input_interval_list.display()
        );
    }
    if let Some(dbsnp) = &config.dbsnp {
        if !dbsnp.exists() {
            bail!("dbsnp VCF not found: {}", dbsnp.display());
        }
    }
    if config.standard_min_confidence_threshold_for_calling < 0.0 {
        bail!("--standard-min-confidence-threshold-for-calling must be non-negative");
    }
    if config.threads == 0 {
        bail!("--threads must be at least 1");
    }
    if config.memory_gb == 0 {
        bail!("--memory-gb must be at least 1");
    }
    if config.native_pair_hmm_threads == 0 {
        bail!("--native-pair-hmm-threads must be at least 1");
    }
    if config.pair_hmm_implementation != "rust" && config.pair_hmm_implementation != "native" {
        bail!("--pair-hmm-implementation must be rust or native");
    }
    Ok(())
}

fn scan_partition(
    config: &ActiveRegionDiscoveryConfig,
    dict: &SequenceDict,
    intervals: &[Interval],
) -> Result<WorkerOutput> {
    let mut output = WorkerOutput::default();
    let mut bam = bam::IndexedReader::from_path(&config.input_bam)
        .with_context(|| format!("opening indexed BAM {}", config.input_bam.display()))?;
    let bam_tid_by_name = bam_tid_by_name(bam.header())?;
    let reference = faidx::Reader::from_path(&config.reference)
        .with_context(|| format!("opening reference FASTA {}", config.reference.display()))?;

    for interval in intervals {
        scan_interval(
            config,
            dict,
            &bam_tid_by_name,
            &reference,
            &mut bam,
            interval,
            &mut output,
        )?;
    }
    Ok(output)
}

fn scan_interval(
    config: &ActiveRegionDiscoveryConfig,
    dict: &SequenceDict,
    bam_tid_by_name: &HashMap<String, u32>,
    reference: &faidx::Reader,
    bam: &mut bam::IndexedReader,
    interval: &Interval,
    output: &mut WorkerOutput,
) -> Result<()> {
    let tid = *bam_tid_by_name.get(&interval.contig).with_context(|| {
        format!(
            "contig '{}' from {} is not present in BAM header",
            interval.contig,
            config.input_interval_list.display()
        )
    })?;
    let contig_len = dict.contig_length(&interval.contig).with_context(|| {
        format!(
            "contig '{}' missing from interval_list header",
            interval.contig
        )
    })?;
    let ref_len = reference.fetch_seq_len(&interval.contig);
    if ref_len == 0 {
        bail!(
            "contig '{}' is not present in reference FASTA {}",
            interval.contig,
            config.reference.display()
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
                "fetching reference sequence {}:{}-{}",
                interval.contig, interval.start, interval.end
            )
        })?;

    bam.fetch((tid as i32, (interval.start - 1) as i64, interval.end as i64))
        .with_context(|| {
            format!(
                "fetching BAM region {}:{}-{}",
                interval.contig, interval.start, interval.end
            )
        })?;

    let mut pileups = bam.pileup();
    pileups.set_max_depth(config.max_depth);
    for pileup in pileups {
        let pileup = pileup.with_context(|| {
            format!(
                "reading pileup for {}:{}-{}",
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

        output.summary.pileup_sites += 1;
        let ref_base = normalize_base(ref_bases[(pos0 - (interval.start - 1)) as usize]);
        let evidence = pileup_evidence(&pileup, ref_base, config);
        if evidence.depth > 0 {
            output.summary.covered_sites += 1;
            output.summary.max_filtered_depth =
                output.summary.max_filtered_depth.max(evidence.depth);
        }
        if is_active_site(evidence, config) {
            output.summary.active_sites += 1;
            output.active_regions.push(padded_interval(
                &interval.contig,
                pos0 + 1,
                config.active_region_padding,
                contig_len,
            ));
        }
    }
    Ok(())
}

fn pileup_evidence(
    pileup: &bam::pileup::Pileup,
    ref_base: u8,
    config: &ActiveRegionDiscoveryConfig,
) -> SiteEvidence {
    let mut evidence = SiteEvidence::default();
    for alignment in pileup.alignments() {
        let record = alignment.record();
        if !standard_hc_read_filter(&record, config) || alignment.is_refskip() {
            continue;
        }

        let base_quality_ok = alignment
            .qpos()
            .and_then(|qpos| record.qual().get(qpos).copied())
            .is_some_and(|quality| quality >= config.min_baseq);

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

pub fn standard_hc_read_filter(record: &bam::Record, config: &ActiveRegionDiscoveryConfig) -> bool {
    const UNMAPPED: u16 = 0x4;
    const SECONDARY: u16 = 0x100;
    const QCFAIL: u16 = 0x200;
    const DUPLICATE: u16 = 0x400;
    const SUPPLEMENTARY: u16 = 0x800;

    let mut excluded = UNMAPPED | SECONDARY | QCFAIL | DUPLICATE;
    if config.exclude_supplementary {
        excluded |= SUPPLEMENTARY;
    }

    record.flags() & excluded == 0
        && record.mapq() >= config.min_mapq
        && cigar_has_reference_bases(record)
}

fn cigar_has_reference_bases(record: &bam::Record) -> bool {
    record.cigar().iter().any(|op| {
        matches!(
            op,
            Cigar::Match(_) | Cigar::Equal(_) | Cigar::Diff(_) | Cigar::Del(_) | Cigar::RefSkip(_)
        )
    })
}

fn is_active_site(evidence: SiteEvidence, config: &ActiveRegionDiscoveryConfig) -> bool {
    if evidence.indel_count >= config.min_indel_count {
        return true;
    }
    if evidence.depth == 0 || evidence.alt_count < config.min_alt_count {
        return false;
    }
    f64::from(evidence.alt_count) / f64::from(evidence.depth) >= config.min_alt_fraction
}

fn padded_interval(contig: &str, pos1: u64, padding: u64, contig_len: u64) -> Interval {
    Interval {
        contig: contig.to_string(),
        start: pos1.saturating_sub(padding).max(1),
        end: pos1.saturating_add(padding).min(contig_len),
    }
}

fn partition_intervals_by_bases(intervals: &[Interval], threads: usize) -> Vec<Vec<Interval>> {
    if intervals.is_empty() {
        return Vec::new();
    }
    let workers = threads.min(intervals.len()).max(1);
    let total_bases: u64 = intervals.iter().map(Interval::len).sum();
    let target_bases = total_bases.div_ceil(workers as u64).max(1);

    let mut partitions = Vec::with_capacity(workers);
    let mut current = Vec::new();
    let mut current_bases = 0_u64;
    for interval in intervals {
        if !current.is_empty() && partitions.len() + 1 < workers && current_bases >= target_bases {
            partitions.push(current);
            current = Vec::new();
            current_bases = 0;
        }
        current_bases += interval.len();
        current.push(interval.clone());
    }
    if !current.is_empty() {
        partitions.push(current);
    }
    partitions
}

fn read_interval_list(path: &Path) -> Result<(SequenceDict, Vec<Interval>)> {
    let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut header_lines = Vec::new();
    let mut body_lines = Vec::new();
    for line in reader.lines() {
        let line = line.with_context(|| format!("reading {}", path.display()))?;
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

fn write_bed(path: &Path, intervals: &[Interval]) -> Result<()> {
    create_parent_dir(path)?;
    let file = File::create(path).with_context(|| format!("creating {}", path.display()))?;
    let mut writer = BufWriter::new(file);
    for interval in intervals {
        writeln!(
            writer,
            "{}\t{}\t{}",
            interval.contig,
            interval.start - 1,
            interval.end
        )
        .with_context(|| format!("writing {}", path.display()))?;
    }
    Ok(())
}

fn write_active_region_summary(
    path: &Path,
    config: &ActiveRegionDiscoveryConfig,
    summary: &ActiveRegionSummary,
) -> Result<()> {
    create_parent_dir(path)?;
    let file = File::create(path).with_context(|| format!("creating {}", path.display()))?;
    let mut writer = BufWriter::new(file);
    writeln!(writer, "metric\tvalue")?;
    writeln!(writer, "mode\tactive_region_discovery")?;
    writeln!(writer, "input_bam\t{}", config.input_bam.display())?;
    writeln!(writer, "reference\t{}", config.reference.display())?;
    writeln!(
        writer,
        "input_interval_list\t{}",
        config.input_interval_list.display()
    )?;
    writeln!(writer, "threads\t{}", config.threads)?;
    writeln!(writer, "min_mapq\t{}", config.min_mapq)?;
    writeln!(writer, "min_baseq\t{}", config.min_baseq)?;
    writeln!(writer, "min_alt_count\t{}", config.min_alt_count)?;
    writeln!(writer, "min_indel_count\t{}", config.min_indel_count)?;
    writeln!(writer, "min_alt_fraction\t{:.6}", config.min_alt_fraction)?;
    writeln!(
        writer,
        "active_region_padding\t{}",
        config.active_region_padding
    )?;
    writeln!(writer, "max_depth_setting\t{}", config.max_depth)?;
    writeln!(
        writer,
        "exclude_supplementary\t{}",
        config.exclude_supplementary
    )?;
    writeln!(writer, "input_intervals\t{}", summary.input_intervals)?;
    writeln!(writer, "input_bases\t{}", summary.input_bases)?;
    writeln!(writer, "pileup_sites\t{}", summary.pileup_sites)?;
    writeln!(writer, "covered_sites\t{}", summary.covered_sites)?;
    writeln!(writer, "active_sites\t{}", summary.active_sites)?;
    writeln!(writer, "active_regions\t{}", summary.active_regions)?;
    writeln!(writer, "active_bases\t{}", summary.active_bases)?;
    writeln!(writer, "max_filtered_depth\t{}", summary.max_filtered_depth)?;
    Ok(())
}

fn create_parent_dir(path: &Path) -> Result<()> {
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent)
                .with_context(|| format!("creating directory {}", parent.display()))?;
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

    fn test_config() -> ActiveRegionDiscoveryConfig {
        ActiveRegionDiscoveryConfig {
            input_bam: PathBuf::from("in.bam"),
            reference: PathBuf::from("ref.fa"),
            input_interval_list: PathBuf::from("in.interval_list"),
            output_active_bed: PathBuf::from("active.bed"),
            output_summary: PathBuf::from("summary.tsv"),
            min_mapq: 20,
            min_baseq: 10,
            min_alt_count: 2,
            min_indel_count: 1,
            min_alt_fraction: 0.2,
            active_region_padding: 150,
            max_depth: 100_000,
            threads: 4,
            exclude_supplementary: false,
        }
    }

    #[test]
    fn active_site_requires_alt_or_indel_threshold() {
        let config = test_config();
        assert!(!is_active_site(
            SiteEvidence {
                depth: 10,
                alt_count: 1,
                indel_count: 0,
            },
            &config
        ));
        assert!(is_active_site(
            SiteEvidence {
                depth: 10,
                alt_count: 2,
                indel_count: 0,
            },
            &config
        ));
        assert!(is_active_site(
            SiteEvidence {
                depth: 0,
                alt_count: 0,
                indel_count: 1,
            },
            &config
        ));
    }

    #[test]
    fn active_region_padding_clamps_to_contig_bounds() {
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
    fn active_regions_merge_by_dictionary_order() {
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
    fn interval_partitioning_preserves_all_bases() {
        let intervals = vec![
            Interval {
                contig: "chr1".to_string(),
                start: 1,
                end: 10,
            },
            Interval {
                contig: "chr1".to_string(),
                start: 11,
                end: 30,
            },
            Interval {
                contig: "chr1".to_string(),
                start: 31,
                end: 60,
            },
        ];
        let total_bases: u64 = intervals.iter().map(Interval::len).sum();
        let partitions = partition_intervals_by_bases(&intervals, 2);
        let partition_bases: u64 = partitions.iter().flatten().map(Interval::len).sum();
        assert_eq!(partition_bases, total_bases);
        assert_eq!(
            partitions.iter().map(Vec::len).sum::<usize>(),
            intervals.len()
        );
    }
}
