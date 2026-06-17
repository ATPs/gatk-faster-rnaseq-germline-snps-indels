use anyhow::{bail, Context, Result};
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use rust_htslib::bam::pileup::Indel;
use rust_htslib::bam::record::Cigar;
use rust_htslib::tbx::Read as TbxRead;
use rust_htslib::{bam, bam::Read, bgzf, faidx, htslib, tbx};
use std::collections::HashMap;
use std::ffi::CString;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::os::unix::ffi::OsStrExt;
use std::path::{Path, PathBuf};
use std::ptr;
use crate::pair_hmm;

const PIPELINE_MIN_MAPQ: u8 = 20;
const PIPELINE_MIN_BASEQ: u8 = 10;
const PIPELINE_MAX_DEPTH: u32 = 100_000;
// GATK: minTailQuality = minBaseQualityScore - 1 for the non-error-correction path.
// See AssemblyBasedCallerUtils.finalizeRegion and HaplotypeCallerEngine.
const PIPELINE_MIN_TAIL_QUALITY: u8 = PIPELINE_MIN_BASEQ - 1;
const FETCH_WINDOW_GAP: u64 = 1_000;
const FETCH_WINDOW_MAX_BASES: u64 = 2_000_000;
const MAX_BOOTSTRAP_INDEL_LEN: u32 = 200;
const MAX_SAM_QUAL: u8 = 60;
const DEFAULT_INDEL_QUAL: u8 = 30;
const HALF_DEFAULT_PCR_SNV_QUAL: u8 = 20;
const FISHER_STRAND_TARGET_TABLE_SIZE: f64 = 200.0;
const FISHER_STRAND_MIN_PVALUE: f64 = 1e-320;
const CALL_PARTITIONS_PER_THREAD: usize = 8;
const ACTIVE_REGION_MAX_GAP: u64 = 50;
const ACTIVE_REGION_PADDING: u64 = 50;

// GATK's isActive() uses a minimal confidence threshold of 4.0 for
// active-region discovery, much lower than the standard calling threshold.
const ACTIVE_REGION_DISCOVERY_CONFIDENCE: f64 = 4.0;

// GATK: reads shorter than this after trimming are removed before genotyping.
// See AssemblyBasedCallerUtils.MINIMUM_READ_LENGTH_AFTER_TRIMMING.
const MIN_READ_LENGTH_AFTER_TRIMMING: usize = 10;

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

#[derive(Debug, Clone)]
pub struct HaplotypeReplayConfig {
    pub caller: HaplotypeCallerConfig,
    pub output_prefix: PathBuf,
}

#[derive(Debug, Clone, Default, Eq, PartialEq)]
pub struct HaplotypeReplayStats {
    pub active_regions: usize,
    pub active_loci: usize,
    pub read_observations: usize,
    pub candidate_events: usize,
    pub genotype_rows: usize,
}

/// Simplified ref-vs-any active probability, mirroring GATK's
/// ReferenceConfidenceModel.calcGenotypeLikelihoodsOfRefVsAny +
/// MinimalGenotypingEngine.calculateSingleSampleRefVsAnyActiveStateProfileValue.
///
/// Returns true if the locus is "active" — i.e., likely to harbor variation.

fn is_active_indel(evidence: &IndelEvidence) -> (bool, u32) {
    if evidence.counts.depth < 2 {
        return (false, 0);
    }
    let mut log10_likelihoods = [0.0_f64; 3];
    for observation in &evidence.observations {
        let error = phred_error_probability(observation.quality);
        let is_ref = matches!(observation.allele, IndelObservationAllele::Ref);
        let ref_prob = snp_observation_probability(is_ref, error);
        let alt_prob = if !is_ref {
            1.0 - error
        } else {
            error / 3.0
        }
        .max(f64::MIN_POSITIVE);
        add_diploid_observation(&mut log10_likelihoods, ref_prob, alt_prob);
    }
    let denominator = log10_sum_exp(&log10_likelihoods);
    let ref_posterior_log10 = log10_likelihoods[0] - denominator;
    let qual = phred_from_log10_probability(ref_posterior_log10);
    (f64::from(qual) >= ACTIVE_REGION_DISCOVERY_CONFIDENCE, qual)
}

fn is_active_locus(ref_index: Option<usize>, evidence: &SnpEvidence, depth: u32) -> (bool, u32) {
    let ref_index = match ref_index {
        Some(idx) => idx,
        None => return (false, 0),
    };
    if depth < 2 {
        return (false, 0);
    }

    // Compute ref-vs-any log10 likelihoods for diploid genotype states.
    let mut log10_likelihoods = [0.0_f64; 3]; // hom-ref, het, hom-alt
    for observation in &evidence.observations {
        let error = phred_error_probability(observation.quality);
        let ref_prob = snp_observation_probability(observation.base_index == ref_index, error);
        // For ref-vs-any, any non-reference base is considered alt evidence.
        let alt_prob = if observation.base_index != ref_index {
            1.0 - error
        } else {
            error / 3.0
        }
        .max(f64::MIN_POSITIVE);
        add_diploid_observation(&mut log10_likelihoods, ref_prob, alt_prob);
    }

    // Use a flat prior (no population prior for active detection, like GATK
    // uses a very low effective confidence threshold).
    let denominator = log10_sum_exp(&log10_likelihoods);
    let ref_posterior_log10 = log10_likelihoods[0] - denominator;
    let qual = phred_from_log10_probability(ref_posterior_log10);

    // Active if QUAL exceeds the low discovery threshold (GATK uses 4.0).
    (f64::from(qual) >= ACTIVE_REGION_DISCOVERY_CONFIDENCE, qual)
}

/// Compute the effective clipped span of a read for evidence collection.
/// Returns (seq_start_qpos, seq_end_qpos_exclusive) where the span excludes:
///  - soft-clipped bases (hard-clipped when dont_use_soft_clipped_bases)
///  - low-quality tail bases (quality < min_tail_quality)
///
/// Mirrors GATK's AssemblyBasedCallerUtils.finalizeRegion() clipping logic:
/// ReadClipper.hardClipSoftClippedBases + ReadClipper.hardClipLowQualEnds.
fn clip_read_for_evidence(
    record: &bam::Record,
    min_tail_quality: u8,
    dont_use_soft_clipped_bases: bool,
) -> Option<(usize, usize)> {
    let seq = record.seq();
    let qual = record.qual();
    if seq.is_empty() {
        return None;
    }

    let seq_len = seq.len();
    let mut qpos = 0_usize;
    let mut align_start: Option<usize> = None;
    let mut align_end: Option<usize> = None;

    for op in record.cigar().iter() {
        let len = match op {
            Cigar::Match(l)
            | Cigar::Equal(l)
            | Cigar::Diff(l)
            | Cigar::Ins(l)
            | Cigar::SoftClip(l) => *l as usize,
            Cigar::Del(_) | Cigar::RefSkip(_) => {
                continue;
            }
            Cigar::HardClip(_) | Cigar::Pad(_) => continue,
        };

        match op {
            Cigar::SoftClip(_) => {
                if dont_use_soft_clipped_bases {
                    qpos += len;
                } else {
                    if align_start.is_none() {
                        align_start = Some(qpos);
                    }
                    qpos += len;
                    align_end = Some(qpos);
                }
            }
            Cigar::Match(_) | Cigar::Equal(_) | Cigar::Diff(_) | Cigar::Ins(_) => {
                if align_start.is_none() {
                    align_start = Some(qpos);
                }
                qpos += len;
                align_end = Some(qpos);
            }
            _ => {}
        }
    }

    let mut start = align_start?;
    let mut end = align_end?;
    if end <= start || start >= seq_len {
        return None;
    }
    end = end.min(seq_len);

    while start < end && qual.get(start).copied().unwrap_or(0) < min_tail_quality {
        start += 1;
    }
    while end > start && qual.get(end - 1).copied().unwrap_or(0) < min_tail_quality {
        end -= 1;
    }

    if end - start < MIN_READ_LENGTH_AFTER_TRIMMING {
        return None;
    }
    Some((start, end))
}

pub fn call_variants(config: &HaplotypeCallerConfig) -> Result<()> {
    validate_haplotype_caller_config(config)?;
    let (dict, mut intervals) = read_interval_list(&config.input_interval_list)?;
    sort_intervals(&mut intervals, &dict)?;
    let fetch_windows = coalesce_fetch_windows(&intervals);
    let partitions = partition_fetch_windows_by_bases(
        &fetch_windows,
        config
            .threads
            .max(1)
            .saturating_mul(CALL_PARTITIONS_PER_THREAD),
    );
    let sample_name = sample_name_from_bam(&config.input_bam)?;

    let thread_pool = ThreadPoolBuilder::new()
        .num_threads(config.threads.max(1))
        .build()
        .context("creating HaplotypeCaller call thread pool")?;
    let worker_outputs: Vec<Result<CallWorkerOutput>> = thread_pool.install(|| {
        partitions
            .into_par_iter()
            .map(|partition| scan_call_partition(config, &partition))
            .collect()
    });

    let mut variants = Vec::new();
    for worker_output in worker_outputs {
        variants.extend(worker_output?.variants);
    }
    sort_variant_calls(&mut variants, &dict)?;
    dedup_variant_calls(&mut variants);
    if let Some(dbsnp) = &config.dbsnp {
        annotate_dbsnp(dbsnp, &mut variants)?;
    }

    write_bootstrap_vcf(&config.output_vcf, config, &dict, &sample_name, &variants)?;
    if is_gzip_path(&config.output_vcf) {
        write_tabix_index(&config.output_vcf, config.threads)?;
    }
    Ok(())
}

pub fn replay_regions(config: &HaplotypeReplayConfig) -> Result<HaplotypeReplayStats> {
    validate_haplotype_caller_config(&config.caller)?;
    let (dict, mut intervals) = read_interval_list(&config.caller.input_interval_list)?;
    sort_intervals(&mut intervals, &dict)?;
    let fetch_windows = coalesce_fetch_windows(&intervals);
    let partitions = partition_fetch_windows_by_bases(
        &fetch_windows,
        config
            .caller
            .threads
            .max(1)
            .saturating_mul(CALL_PARTITIONS_PER_THREAD),
    );

    let thread_pool = ThreadPoolBuilder::new()
        .num_threads(config.caller.threads.max(1))
        .build()
        .context("creating HaplotypeCaller replay thread pool")?;
    let worker_outputs: Vec<Result<ReplayWorkerOutput>> = thread_pool.install(|| {
        partitions
            .into_par_iter()
            .map(|partition| scan_replay_partition(&config.caller, &partition))
            .collect()
    });

    let mut output = ReplayWorkerOutput::default();
    for worker_output in worker_outputs {
        output.extend(worker_output?);
    }
    sort_replay_rows(&mut output, &dict)?;
    sort_variant_calls(&mut output.variants, &dict)?;
    dedup_variant_calls(&mut output.variants);
    if let Some(dbsnp) = &config.caller.dbsnp {
        annotate_dbsnp(dbsnp, &mut output.variants)?;
    }

    let genotype_rows: Vec<ReplayGenotypeRow> = output
        .variants
        .iter()
        .map(ReplayGenotypeRow::from)
        .collect();
    let stats = HaplotypeReplayStats {
        active_regions: output.active_regions.len(),
        active_loci: output.active_loci.len(),
        read_observations: output.read_observations.len(),
        candidate_events: output.events.len(),
        genotype_rows: genotype_rows.len(),
    };
    write_replay_tables(&config.output_prefix, &output, &genotype_rows)?;
    Ok(stats)
}

#[derive(Clone, Debug)]
struct DictRecord {
    name: String,
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

#[derive(Clone, Debug, Eq, PartialEq)]
struct FetchWindow {
    contig: String,
    start: u64,
    end: u64,
    intervals: Vec<Interval>,
}

impl FetchWindow {
    fn len(&self) -> u64 {
        self.end - self.start + 1
    }
}

#[derive(Clone, Debug, PartialEq)]
struct VariantCall {
    contig: String,
    pos: u64,
    id: Option<String>,
    db: bool,
    ref_allele: Vec<u8>,
    alt_allele: Vec<u8>,
    depth: u32,
    ref_count: u32,
    alt_count: u32,
    qual: u32,
    fs: f64,
    pl: [u32; 3],
    genotype_index: usize,
}

impl VariantCall {
    fn genotype(&self) -> &'static str {
        match self.genotype_index {
            0 => "0/0",
            1 => "0/1",
            2 => "1/1",
            _ => "0/1",
        }
    }

    fn alt_allele_count(&self) -> u32 {
        match self.genotype_index {
            0 => 0,
            1 => 1,
            2 => 2,
            _ => 1,
        }
    }

    fn gq(&self) -> u32 {
        self.pl
            .iter()
            .enumerate()
            .filter(|(idx, _)| *idx != self.genotype_index)
            .map(|(_, pl)| *pl)
            .min()
            .unwrap_or(0)
            .min(99)
    }
}

#[derive(Copy, Clone, Debug, Default, Eq, PartialEq)]
struct BaseCounts {
    counts: [u32; 4],
    depth: u32,
}

#[derive(Copy, Clone, Debug, Default, Eq, PartialEq)]
struct StrandCounts {
    forward: u32,
    reverse: u32,
}

impl StrandCounts {
    fn increment(&mut self, is_reverse: bool) {
        if is_reverse {
            self.reverse += 1;
        } else {
            self.forward += 1;
        }
    }
}

#[derive(Copy, Clone, Debug, Eq, PartialEq)]
struct BaseObservation {
    base_index: usize,
    quality: u8,
    is_reverse: bool,
}

#[derive(Clone, Debug, Default, Eq, PartialEq)]
struct SnpEvidence {
    counts: BaseCounts,
    strands: [StrandCounts; 4],
    observations: Vec<BaseObservation>,
}

#[derive(Clone, Debug, Eq, Hash, PartialEq)]
enum IndelAllele {
    Insertion(Vec<u8>),
    Deletion(u32),
}

#[derive(Clone, Debug, Default, Eq, PartialEq)]
struct IndelCounts {
    depth: u32,
    ref_count: u32,
    counts: HashMap<IndelAllele, u32>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
enum IndelObservationAllele {
    Ref,
    Alt(IndelAllele),
}

#[derive(Clone, Debug, Eq, PartialEq)]
struct IndelObservation {
    allele: IndelObservationAllele,
    quality: u8,
    is_reverse: bool,
}

#[derive(Clone, Debug, Default, PartialEq)]
struct IndelEvidence {
    counts: IndelCounts,
    ref_strand: StrandCounts,
    alt_strands: HashMap<IndelAllele, StrandCounts>,
    observations: Vec<IndelObservation>,
}

#[derive(Clone, Debug)]
struct VariantModel {
    qual: u32,
    pl: [u32; 3],
    genotype_index: usize,
}

#[derive(Default)]
struct CallWorkerOutput {
    variants: Vec<VariantCall>,
}

#[derive(Clone, Debug, Default)]
struct ReplayWorkerOutput {
    variants: Vec<VariantCall>,
    active_regions: Vec<ReplayActiveRegionRow>,
    active_loci: Vec<ReplayActiveLocusRow>,
    read_observations: Vec<ReplayReadObservationRow>,
    events: Vec<ReplayEventRow>,
    haplotypes: Vec<ReplayHaplotypeRow>,
    pairhmms: Vec<ReplayPairHmmRow>,
}

impl ReplayWorkerOutput {
    fn extend(&mut self, other: ReplayWorkerOutput) {
        self.variants.extend(other.variants);
        self.active_regions.extend(other.active_regions);
        self.active_loci.extend(other.active_loci);
        self.read_observations.extend(other.read_observations);
        self.events.extend(other.events);
        self.haplotypes.extend(other.haplotypes);
        self.pairhmms.extend(other.pairhmms);
    }
}

#[derive(Clone, Debug, Default)]
struct ReplayActiveRegionRow {
    contig: String,
    start: u64,
    end: u64,
    region: String,
    observed_loci: u64,
    active_loci: u64,
    candidate_events: u64,
    max_alt_fraction: f64,
    mean_alt_fraction: f64,
}

#[derive(Clone, Debug)]
struct ReplayActiveLocusRow {
    contig: String,
    pos: u64,
    region: String,
    ref_base: u8,
    depth: u32,
    snp_alt_count: u32,
    snp_best_alt: String,
    indel_alt_count: u32,
    indel_best_alt: String,
    alt_fraction: f64,
    active_probability_proxy: f64,
}

#[derive(Clone, Debug)]
struct ReplayReadObservationRow {
    region: String,
    read: String,
    kind: &'static str,
    pos: u64,
    qpos: usize,
    allele: String,
    adjusted_quality: u8,
    mapq: u8,
    strand: &'static str,
}

#[derive(Clone, Debug)]
struct ReplayEventRow {
    region: String,
    event: String,
    chrom: String,
    pos: u64,
    event_type: &'static str,
    alleles: String,
    raw: String,
    depth: u32,
    ref_count: u32,
    alt_count: u32,
    qual: u32,
    gt: &'static str,
}

#[derive(Clone, Debug)]
struct ReplayGenotypeRow {
    chrom: String,
    pos: u64,
    ref_allele: String,
    alt: String,
    qual: u32,
    filter: &'static str,
    gt: &'static str,
    gq: u32,
    dp: u32,
    ad_ref: u32,
    ad_alt: u32,
    fs: f64,
    qd: f64,
    pl: String,
    db: bool,
}

#[derive(Clone, Debug)]
struct ReplayHaplotypeRow {
    region: String,
    stage: &'static str,
    haplotype: usize,
    span_start: u64,
    span_end: u64,
    kmer: u32,
    length: u32,
    cigar: String,
    is_ref: bool,
    bases: String,
}

#[derive(Clone, Debug)]
struct ReplayPairHmmRow {
    region: String,
    read: String,
    haplotype: usize,
    read_index: usize,
    cigar: String,
    mapq: u8,
    loc: u64,
    unclipped_loc: u64,
    length: u32,
    score: f64,
}

impl From<&VariantCall> for ReplayGenotypeRow {
    fn from(variant: &VariantCall) -> Self {
        let qd = if variant.depth == 0 {
            0.0
        } else {
            fix_too_high_qd(f64::from(variant.qual) / f64::from(variant.depth))
        };
        ReplayGenotypeRow {
            chrom: variant.contig.clone(),
            pos: variant.pos,
            ref_allele: String::from_utf8_lossy(&variant.ref_allele).into_owned(),
            alt: String::from_utf8_lossy(&variant.alt_allele).into_owned(),
            qual: variant.qual,
            filter: "PASS",
            gt: variant.genotype(),
            gq: variant.gq(),
            dp: variant.depth,
            ad_ref: variant.ref_count,
            ad_alt: variant.alt_count,
            fs: variant.fs,
            qd,
            pl: genotype_likelihoods(variant),
            db: variant.db,
        }
    }
}

#[derive(Clone, Debug)]
struct ReplayRowContext<'a> {
    region: &'a str,
    pos: u64,
}

#[derive(Clone, Debug)]
struct NamedBaseObservation {
    read_name: Vec<u8>,
    qpos: usize,
    mapq: u8,
    observation: BaseObservation,
}

#[derive(Clone, Debug)]
struct NamedIndelObservation {
    read_name: Vec<u8>,
    qpos: usize,
    mapq: u8,
    observation: IndelObservation,
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

fn scan_call_partition(
    config: &HaplotypeCallerConfig,
    windows: &[FetchWindow],
) -> Result<CallWorkerOutput> {
    let mut output = CallWorkerOutput::default();
    let mut bam = bam::IndexedReader::from_path(&config.input_bam)
        .with_context(|| format!("opening indexed BAM {}", config.input_bam.display()))?;
    let bam_tid_by_name = bam_tid_by_name(bam.header())?;
    let reference = faidx::Reader::from_path(&config.reference)
        .with_context(|| format!("opening reference FASTA {}", config.reference.display()))?;

    for window in windows {
        scan_call_window(
            config,
            &bam_tid_by_name,
            &reference,
            &mut bam,
            window,
            &mut output,
        )?;
    }
    Ok(output)
}

fn scan_replay_partition(
    config: &HaplotypeCallerConfig,
    windows: &[FetchWindow],
) -> Result<ReplayWorkerOutput> {
    let mut output = ReplayWorkerOutput::default();
    let mut bam = bam::IndexedReader::from_path(&config.input_bam)
        .with_context(|| format!("opening indexed BAM {}", config.input_bam.display()))?;
    let bam_tid_by_name = bam_tid_by_name(bam.header())?;
    let reference = faidx::Reader::from_path(&config.reference)
        .with_context(|| format!("opening reference FASTA {}", config.reference.display()))?;

    for window in windows {
        scan_replay_window(
            config,
            &bam_tid_by_name,
            &reference,
            &mut bam,
            window,
            &mut output,
        )?;
    }
    Ok(output)
}


fn discover_active_regions(
    config: &HaplotypeCallerConfig,
    tid: u32,
    window: &FetchWindow,
    ref_bases: &[u8],
    bam: &mut bam::IndexedReader,
) -> Result<Vec<Interval>> {
    let mut active_regions: Vec<Interval> = Vec::new();
    let mut activity_scores: std::collections::HashMap<u64, u32> = std::collections::HashMap::new();

    bam.fetch((tid as i32, (window.start - 1) as i64, window.end as i64))
        .with_context(|| {
            format!(
                "fetching BAM region {}:{}-{} for active region discovery",
                window.contig, window.start, window.end
            )
        })?;
    let mut pileups = bam.pileup();
    pileups.set_max_depth(PIPELINE_MAX_DEPTH);
    let mut interval_cursor = 0_usize;
    let mut current_region: Option<Interval> = None;

    for pileup in pileups {
        let pileup = pileup.with_context(|| {
            format!(
                "reading pileup for {}:{}-{} active discovery",
                window.contig, window.start, window.end
            )
        })?;
        if pileup.tid() != tid {
            continue;
        }
        let pos0 = u64::from(pileup.pos());
        if pos0 < window.start - 1 || pos0 >= window.end {
            continue;
        }
        let pos1 = pos0 + 1;
        if !position_is_requested(&window.intervals, pos1, &mut interval_cursor) {
            continue;
        }

        let ref_base = normalize_base(ref_bases[(pos0 - (window.start - 1)) as usize]);
        let snp_evidence = pileup_snp_evidence(
            &pileup,
            PIPELINE_MIN_BASEQ,
            PIPELINE_MIN_MAPQ,
            PIPELINE_MIN_TAIL_QUALITY,
            config.dont_use_soft_clipped_bases,
        );
        let indel_evidence = pileup_indel_evidence(
            &pileup,
            PIPELINE_MIN_BASEQ,
            PIPELINE_MIN_MAPQ,
            PIPELINE_MIN_TAIL_QUALITY,
            config.dont_use_soft_clipped_bases,
        );
        let depth = snp_evidence.counts.depth.max(indel_evidence.counts.depth);

        let (snp_active, snp_qual) = is_active_locus(base_index(ref_base), &snp_evidence, depth);
        let (indel_active, indel_qual) = is_active_indel(&indel_evidence);
        let is_active = snp_active || indel_active;
        let qual = snp_qual.max(indel_qual);

        if is_active {
            activity_scores.insert(pos1, qual);
            if let Some(ref mut reg) = current_region {
                if pos1 <= reg.end + ACTIVE_REGION_MAX_GAP {
                    reg.end = pos1;
                } else {
                    active_regions.push(reg.clone());
                    *reg = Interval {
                        contig: window.contig.clone(),
                        start: pos1,
                        end: pos1,
                    };
                }
            } else {
                current_region = Some(Interval {
                    contig: window.contig.clone(),
                    start: pos1,
                    end: pos1,
                });
            }
        }
    }
    if let Some(reg) = current_region {
        active_regions.push(reg);
    }

    let mut split_active_regions = Vec::new();
    for reg in active_regions {
        let mut start = reg.start;
        let end = reg.end;
        while end - start + 1 > 300 {
            let search_start = start + 50;
            let search_end = (start + 300).min(end - 50);

            let mut best_cut = search_start;
            let mut min_score = u32::MAX;

            for p in search_start..=search_end {
                let score = *activity_scores.get(&p).unwrap_or(&0);
                if score < min_score {
                    min_score = score;
                    best_cut = p;
                } else if score == min_score {
                    let mid = (search_start + search_end) / 2;
                    if (p as i64 - mid as i64).abs() < (best_cut as i64 - mid as i64).abs() {
                        best_cut = p;
                    }
                }
            }

            split_active_regions.push(Interval {
                contig: reg.contig.clone(),
                start,
                end: best_cut - 1,
            });
            start = best_cut;
        }
        split_active_regions.push(Interval {
            contig: reg.contig.clone(),
            start,
            end,
        });
    }

    let mut coalesced: Vec<Interval> = Vec::new();
    for mut reg in split_active_regions {
        reg.start = reg.start.saturating_sub(ACTIVE_REGION_PADDING).max(window.start);
        reg.end = reg.end.saturating_add(ACTIVE_REGION_PADDING).min(window.end);
        coalesced.push(reg);
    }

    Ok(coalesced)
}

fn scan_call_window(
    config: &HaplotypeCallerConfig,
    bam_tid_by_name: &HashMap<String, u32>,
    reference: &faidx::Reader,
    bam: &mut bam::IndexedReader,
    window: &FetchWindow,
    output: &mut CallWorkerOutput,
) -> Result<()> {
    let tid = *bam_tid_by_name.get(&window.contig).with_context(|| {
        format!(
            "contig '{}' from {} is not present in BAM header",
            window.contig,
            config.input_interval_list.display()
        )
    })?;
    let ref_len = reference.fetch_seq_len(&window.contig);
    if ref_len == 0 {
        bail!("contig '{}' is not present in reference FASTA {}", window.contig, config.reference.display());
    }
    if window.end > ref_len {
        bail!("interval {}:{}-{} extends past FASTA contig length {}", window.contig, window.start, window.end, ref_len);
    }

    let ref_end = window.end.saturating_add(u64::from(MAX_BOOTSTRAP_INDEL_LEN)).min(ref_len);
    let ref_bases = reference.fetch_seq(&window.contig, (window.start - 1) as usize, (ref_end - 1) as usize)
        .with_context(|| format!("fetching reference sequence {}:{}-{}", window.contig, window.start, ref_end))?;

    let active_regions = discover_active_regions(config, tid, window, &ref_bases, bam)?;

    for region in active_regions {
        bam.fetch((tid as i32, (region.start - 1) as i64, region.end as i64))
            .with_context(|| format!("fetching BAM region {}:{}-{}", region.contig, region.start, region.end))?;

        let mut read_bases_list = Vec::new();
        let mut read_quals_list = Vec::new();
        let mut read_is_reverse_list = Vec::new();
        
        for r in bam.records() {
            let record = r?;
            if record.is_unmapped() || record.is_secondary() || record.is_supplementary() || record.is_duplicate() || record.is_quality_check_failed() {
                continue;
            }
            if record.mapq() < PIPELINE_MIN_MAPQ {
                continue;
            }

            let mapq = record.mapq();
            let mut r_bases = Vec::new();
            let mut r_quals = Vec::new();
            let seq = record.seq();
            let quals = record.qual();
            let mut seq_idx = 0;

            for view in record.cigar().iter() {
                use rust_htslib::bam::record::Cigar::*;
                let len = view.len() as usize;
                match view {
                    Match(_) | Equal(_) | Diff(_) | Ins(_) => {
                        for i in seq_idx..seq_idx + len {
                            if i < seq.len() {
                                r_bases.push(seq[i]);
                                r_quals.push(quals[i].min(mapq).max(18));
                            }
                        }
                        seq_idx += len;
                    }
                    SoftClip(_) => {
                        seq_idx += len;
                    }
                    Del(_) | RefSkip(_) => {}
                    HardClip(_) | Pad(_) => {}
                }
            }

            if !r_bases.is_empty() {
                read_bases_list.push(r_bases);
                read_quals_list.push(r_quals);
                read_is_reverse_list.push(record.is_reverse());
            }
        }

        if read_bases_list.is_empty() {
            continue;
        }

        let ref_region_start_offset = (region.start - window.start) as usize;
        let ref_region_end_offset = (region.end - window.start) as usize;
        let local_ref_bases = &ref_bases[ref_region_start_offset..=ref_region_end_offset];

        // 1. Assemble haplotypes instead of pileup!
        let max_mnp_distance = 0; // Default GATK
        // Cap reads for assembly to avoid blowup on deep regions
        let assembly_reads: Vec<&Vec<u8>> = read_bases_list.iter().take(500).collect();
        let assembly_reads_owned: Vec<Vec<u8>> = assembly_reads.iter().map(|r| (*r).clone()).collect();
        let (local_haplotypes, valid_events) = assemble_haplotypes(
            &region.contig,
            region.start,
            local_ref_bases,
            &assembly_reads_owned,
            &[10, 25],
            max_mnp_distance,
        );

        if valid_events.is_empty() {
            continue;
        }

        let n_reads = read_bases_list.len();
        let mut read_haplotype_likelihoods: Vec<Vec<f64>> = Vec::with_capacity(n_reads);

        for i in 0..n_reads {
            let r_bases = &read_bases_list[i];
            let r_quals = &read_quals_list[i];
            let read_ins_quals = vec![45; r_bases.len()];
            let read_del_quals = vec![45; r_bases.len()];
            let gcp = 10;

            let mut hap_likelihoods = Vec::with_capacity(local_haplotypes.len());
            for hap in &local_haplotypes {
                let score = pair_hmm::compute_read_likelihood_given_haplotype(
                    &hap.bases,
                    r_bases,
                    r_quals,
                    &read_ins_quals,
                    &read_del_quals,
                    gcp,
                );
                hap_likelihoods.push(score);
            }
            read_haplotype_likelihoods.push(hap_likelihoods);
        }

        if read_haplotype_likelihoods.is_empty() {
            continue;
        }

        for (event_idx, event) in valid_events.iter().enumerate() {
            let mut hap_contains_allele = vec![false; local_haplotypes.len()];
            for (hap_idx, hap) in local_haplotypes.iter().enumerate() {
                if hap.event_indices.contains(&event_idx) {
                    hap_contains_allele[hap_idx] = true;
                }
            }

            let mut read_allele_likelihoods_ref = vec![0.0; n_reads];
            let mut read_allele_likelihoods_alt = vec![0.0; n_reads];

            for read_idx in 0..n_reads {
                let mut best_ref_lhood = f64::NEG_INFINITY;
                let mut best_alt_lhood = f64::NEG_INFINITY;
                for hap_idx in 0..local_haplotypes.len() {
                    let lhood = read_haplotype_likelihoods[read_idx][hap_idx];
                    if hap_contains_allele[hap_idx] {
                        if lhood > best_alt_lhood { best_alt_lhood = lhood; }
                    } else {
                        if lhood > best_ref_lhood { best_ref_lhood = lhood; }
                    }
                }
                read_allele_likelihoods_ref[read_idx] = best_ref_lhood;
                read_allele_likelihoods_alt[read_idx] = best_alt_lhood;
            }

            let mut ref_count = 0;
            let mut alt_count = 0;
            let mut ref_strand = StrandCounts::default();
            let mut alt_strand = StrandCounts::default();

            for read_idx in 0..n_reads {
                let rL = read_allele_likelihoods_ref[read_idx];
                let aL = read_allele_likelihoods_alt[read_idx];
                let is_rev = read_is_reverse_list[read_idx];
                if rL - aL > 0.2 {
                    ref_count += 1;
                    ref_strand.increment(is_rev);
                } else if aL - rL > 0.2 {
                    alt_count += 1;
                    alt_strand.increment(is_rev);
                }
            }
            let depth = ref_count + alt_count;
            let fs = fisher_strand_score(ref_strand, alt_strand);

            let mut log10_likelihoods = [0.0; 3];
            for read_idx in 0..n_reads {
                let rL = read_allele_likelihoods_ref[read_idx];
                let aL = read_allele_likelihoods_alt[read_idx];
                log10_likelihoods[0] += rL;
                log10_likelihoods[1] += log10_sum_exp(&[rL, aL]) - 0.3010299956639812;
                log10_likelihoods[2] += aL;
            }

            let heterozygosity: f64 = 1e-3;
            let log10_priors = [
                (1.0 - 1.5 * heterozygosity).log10(),
                heterozygosity.log10(),
                (0.5 * heterozygosity).log10(),
            ];

            let log10_posteriors = [
                log10_likelihoods[0] + log10_priors[0],
                log10_likelihoods[1] + log10_priors[1],
                log10_likelihoods[2] + log10_priors[2],
            ];
            
            let genotype_index = if log10_posteriors[0] >= log10_posteriors[1] && log10_posteriors[0] >= log10_posteriors[2] { 0 }
            else if log10_posteriors[1] >= log10_posteriors[0] && log10_posteriors[1] >= log10_posteriors[2] { 1 }
            else { 2 };
            
            let denominator = log10_sum_exp(&log10_posteriors);
            let ref_posterior_log10 = log10_posteriors[0] - denominator;
            let qual = phred_from_log10_probability(ref_posterior_log10);
            let best_posterior = log10_posteriors[genotype_index];
            
            let mut pl = [
                phred_likelihood_delta(log10_posteriors[0], best_posterior),
                phred_likelihood_delta(log10_posteriors[1], best_posterior),
                phred_likelihood_delta(log10_posteriors[2], best_posterior),
            ];
            
            if pl[0] == 0 && pl[1] == 0 && pl[2] == 0 {
                for i in 0..3 {
                    if i != genotype_index { pl[i] = 9999; }
                }
            }
            
            let min_confidence = config.standard_min_confidence_threshold_for_calling;
            if (qual as f64) >= min_confidence {
                let mut final_call = event.clone();
                final_call.pl = pl;
                final_call.genotype_index = genotype_index;
                final_call.qual = qual;
                final_call.depth = depth;
                final_call.ref_count = ref_count;
                final_call.alt_count = alt_count;
                final_call.fs = fs;
                
                output.variants.push(final_call);
            }
        }
    }
    Ok(())
}

fn scan_replay_window(
    config: &HaplotypeCallerConfig,
    bam_tid_by_name: &HashMap<String, u32>,
    reference: &faidx::Reader,
    bam: &mut bam::IndexedReader,
    window: &FetchWindow,
    output: &mut ReplayWorkerOutput,
) -> Result<()> {
    let tid = *bam_tid_by_name.get(&window.contig).with_context(|| {
        format!(
            "contig '{}' from {} is not present in BAM header",
            window.contig,
            config.input_interval_list.display()
        )
    })?;
    let ref_len = reference.fetch_seq_len(&window.contig);
    if ref_len == 0 {
        bail!(
            "contig '{}' is not present in reference FASTA {}",
            window.contig,
            config.reference.display()
        );
    }
    if window.end > ref_len {
        bail!(
            "interval {}:{}-{} extends past FASTA contig length {}",
            window.contig,
            window.start,
            window.end,
            ref_len
        );
    }

    let ref_end = window
        .end
        .saturating_add(u64::from(MAX_BOOTSTRAP_INDEL_LEN))
        .min(ref_len);
    let ref_bases = reference
        .fetch_seq(
            &window.contig,
            (window.start - 1) as usize,
            (ref_end - 1) as usize,
        )
        .with_context(|| {
            format!(
                "fetching reference sequence {}:{}-{}",
                window.contig, window.start, ref_end
            )
        })?;

    let active_regions = discover_active_regions(config, tid, window, &ref_bases, bam)?;

    let mut interval_cursor = 0_usize;
    for region_interval in active_regions {
        bam.fetch((tid as i32, (region_interval.start - 1) as i64, region_interval.end as i64))
            .with_context(|| {
                format!(
                    "fetching BAM region {}:{}-{}",
                    region_interval.contig, region_interval.start, region_interval.end
                )
            })?;

        let region = region_name(&region_interval.contig, region_interval.start, region_interval.end);
        let mut active_region = ReplayActiveRegionRow {
            contig: region_interval.contig.clone(),
            start: region_interval.start,
            end: region_interval.end,
            region: region.clone(),
            ..ReplayActiveRegionRow::default()
        };        let mut region_events = Vec::new();

        let mut pileups = bam.pileup();
        pileups.set_max_depth(PIPELINE_MAX_DEPTH);
        for pileup in pileups {
            let pileup = pileup.with_context(|| {
                format!(
                    "reading pileup for {}:{}-{}",
                    region_interval.contig, region_interval.start, region_interval.end
                )
            })?;
            if pileup.tid() != tid {
                continue;
            }
            let pos0 = u64::from(pileup.pos());
            if pos0 < region_interval.start - 1 || pos0 >= region_interval.end {
                continue;
            }
            let pos1 = pos0 + 1;
            if !position_is_requested(&window.intervals, pos1, &mut interval_cursor) {
                continue;
            }

            let ref_base = normalize_base(ref_bases[(pos0 - (window.start - 1)) as usize]);
            let row_context = ReplayRowContext {
                region: &region,
                pos: pos1,
            };
            let (snp_evidence, snp_rows) = pileup_snp_evidence_with_rows(
                &pileup,
                PIPELINE_MIN_BASEQ,
                PIPELINE_MIN_MAPQ,
                PIPELINE_MIN_TAIL_QUALITY,
                config.dont_use_soft_clipped_bases,
                Some(&row_context),
            );
            let (indel_evidence, indel_rows) = pileup_indel_evidence_with_rows(
                &pileup,
                PIPELINE_MIN_BASEQ,
                PIPELINE_MIN_MAPQ,
                PIPELINE_MIN_TAIL_QUALITY,
                config.dont_use_soft_clipped_bases,
                Some(&row_context),
            );
            output.read_observations.extend(snp_rows);
            output.read_observations.extend(indel_rows);

            let ref_index = base_index(ref_base);
            let snp_alt = best_snp_alt(ref_index, &snp_evidence);
            let indel_alt = best_indel_alt(&indel_evidence);
            let snp_alt_count = snp_alt.map(|(_, count)| count).unwrap_or(0);
            let indel_alt_count = indel_alt.map(|(_, count)| *count).unwrap_or(0);
            let depth = snp_evidence.counts.depth.max(indel_evidence.counts.depth);
            let best_alt_count = snp_alt_count.max(indel_alt_count);
            let alt_fraction = if depth == 0 {
                0.0
            } else {
                f64::from(best_alt_count) / f64::from(depth)
            };
            let is_active = is_active_locus(ref_index, &snp_evidence, depth).0 || is_active_indel(&indel_evidence).0;
            active_region.observed_loci += 1;
            active_region.max_alt_fraction = active_region.max_alt_fraction.max(alt_fraction);
            active_region.mean_alt_fraction += alt_fraction;

            if is_active {
                active_region.active_loci += 1;
                output.active_loci.push(ReplayActiveLocusRow {
                    contig: region_interval.contig.clone(),
                    pos: pos1,
                    region: region.clone(),
                    ref_base,
                    depth,
                    snp_alt_count,
                    snp_best_alt: snp_alt
                        .map(|(base_index, _)| base_from_index(base_index) as char)
                        .map(|base| base.to_string())
                        .unwrap_or_default(),
                    indel_alt_count,
                    indel_best_alt: indel_alt
                        .as_ref()
                        .map(|(allele, _)| indel_allele_label(allele))
                        .unwrap_or_default(),
                    alt_fraction,
                    active_probability_proxy: is_active as u8 as f64,
                });
            }

            if let Some(call) = best_snp_call(
                &region_interval.contig,
                pos1,
                ref_base,
                snp_evidence,
                config.standard_min_confidence_threshold_for_calling,
            ) {
                active_region.candidate_events += 1;
                output.events.push(replay_event_row(&region, &call)?);
                region_events.push(call);
            }
            if let Some(call) = best_indel_call(
                &region_interval.contig,
                pos1,
                window.start,
                &ref_bases,
                indel_evidence,
                config.standard_min_confidence_threshold_for_calling,
            ) {
                active_region.candidate_events += 1;
                output.events.push(replay_event_row(&region, &call)?);
                region_events.push(call);
            }
        }

        let ref_region_start_offset = (region_interval.start - window.start) as usize;
        let ref_region_end_offset = (region_interval.end - window.start) as usize;
        let local_ref_bases = &ref_bases[ref_region_start_offset..=ref_region_end_offset];
        
        let local_haplotypes = build_local_haplotypes(
            &region_interval.contig,
            region_interval.start,
            region_interval.end,
            local_ref_bases,
            &region_events,
            128,
        );

        // Fetch reads for this active region for PairHMM
        bam.fetch((tid as i32, (region_interval.start - 1) as i64, region_interval.end as i64))
            .with_context(|| format!("fetching BAM region for PairHMM: {}:{}-{}", region_interval.contig, region_interval.start, region_interval.end))?;
            
        let mut read_index = 0;
        for r in bam.records() {
            let record = r?;
            if record.is_unmapped() || record.is_secondary() || record.is_supplementary() || record.is_duplicate() || record.is_quality_check_failed() {
                continue;
            }
            if record.mapq() < PIPELINE_MIN_MAPQ {
                continue;
            }

            let mapq = record.mapq();
            let mut read_bases = Vec::new();
            let mut read_quals = Vec::new();
            let seq = record.seq();
            let quals = record.qual();
            let mut seq_idx = 0;

            for view in record.cigar().iter() {
                use rust_htslib::bam::record::Cigar::*;
                let len = view.len() as usize;
                match view {
                    Match(_) | Equal(_) | Diff(_) | Ins(_) => {
                        for i in seq_idx..seq_idx + len {
                            read_bases.push(seq[i]);
                            read_quals.push(quals[i].min(mapq).max(18));
                        }
                        seq_idx += len;
                    }
                    SoftClip(_) => {
                        seq_idx += len;
                    }
                    Del(_) | RefSkip(_) => {}
                    HardClip(_) | Pad(_) => {}
                }
            }

            if read_bases.is_empty() {
                continue;
            }

            let read_ins_quals = vec![45; read_bases.len()];
            let read_del_quals = vec![45; read_bases.len()];
            let gcp = 10;

            let read_name = String::from_utf8_lossy(record.qname()).into_owned();
            let unclipped_loc = (record.pos() + 1) as u64; // Approximation for debug output

            for (i, hap) in local_haplotypes.iter().enumerate() {
                let score = pair_hmm::compute_read_likelihood_given_haplotype(
                    &hap.bases,
                    &read_bases,
                    &read_quals,
                    &read_ins_quals,
                    &read_del_quals,
                    gcp,
                );

                output.pairhmms.push(ReplayPairHmmRow {
                    region: region.clone(),
                    read: read_name.clone(),
                    haplotype: i,
                    read_index,
                    cigar: record.cigar().to_string(),
                    mapq,
                    loc: unclipped_loc, // Approximation
                    unclipped_loc,
                    length: read_bases.len() as u32,
                    score,
                });
            }
            read_index += 1;
        }

        for (i, hap) in local_haplotypes.into_iter().enumerate() {
            output.haplotypes.push(ReplayHaplotypeRow {
                region: region.clone(),
                stage: "unclipped",
                haplotype: i,
                span_start: region_interval.start,
                span_end: region_interval.end,
                kmer: 0,
                length: hap.bases.len() as u32,
                cigar: hap.cigar,
                is_ref: hap.is_ref,
                bases: String::from_utf8_lossy(&hap.bases).into_owned(),
            });
        }

        output.variants.extend(region_events);
        if active_region.observed_loci > 0 {
            active_region.mean_alt_fraction /= active_region.observed_loci as f64;
        }
        output.active_regions.push(active_region);
    }
    Ok(())
}

fn pileup_snp_evidence(
    pileup: &bam::pileup::Pileup,
    min_baseq: u8,
    min_mapq: u8,
    min_tail_quality: u8,
    dont_use_soft_clipped_bases: bool,
) -> SnpEvidence {
    pileup_snp_evidence_with_rows(
        pileup,
        min_baseq,
        min_mapq,
        min_tail_quality,
        dont_use_soft_clipped_bases,
        None,
    )
    .0
}

fn pileup_snp_evidence_with_rows(
    pileup: &bam::pileup::Pileup,
    min_baseq: u8,
    min_mapq: u8,
    min_tail_quality: u8,
    dont_use_soft_clipped_bases: bool,
    row_context: Option<&ReplayRowContext<'_>>,
) -> (SnpEvidence, Vec<ReplayReadObservationRow>) {
    let mut observations_by_fragment: HashMap<Vec<u8>, Vec<NamedBaseObservation>> = HashMap::new();
    for alignment in pileup.alignments() {
        let record = alignment.record();
        if !read_passes_hc_filter(&record, min_mapq, false) || alignment.is_refskip() {
            continue;
        }
        let Some(qpos) = alignment.qpos() else {
            continue;
        };
        if record
            .qual()
            .get(qpos)
            .is_none_or(|quality| *quality < min_baseq)
        {
            continue;
        }

        // Apply GATK-like read clipping (hard-clip soft-clips, low-quality tails).
        if let Some((clip_start, clip_end)) =
            clip_read_for_evidence(&record, min_tail_quality, dont_use_soft_clipped_bases)
        {
            if qpos < clip_start || qpos >= clip_end {
                continue;
            }
        } else {
            continue;
        }

        let base = normalize_base(record.seq()[qpos]);
        let Some(index) = base_index(base) else {
            continue;
        };
        observations_by_fragment
            .entry(record.qname().to_vec())
            .or_default()
            .push(NamedBaseObservation {
                read_name: record.qname().to_vec(),
                qpos,
                mapq: record.mapq(),
                observation: BaseObservation {
                    base_index: index,
                    quality: record.qual()[qpos],
                    is_reverse: record.is_reverse(),
                },
            });
    }

    let mut evidence = SnpEvidence::default();
    let mut rows = Vec::new();
    for observations in observations_by_fragment.into_values() {
        for named in adjust_named_base_observations(&observations) {
            let observation = named.observation;
            if observation.quality < min_baseq {
                continue;
            }
            evidence.counts.counts[observation.base_index] += 1;
            evidence.counts.depth += 1;
            evidence.strands[observation.base_index].increment(observation.is_reverse);
            evidence.observations.push(observation);
            if let Some(context) = row_context {
                rows.push(ReplayReadObservationRow {
                    region: context.region.to_string(),
                    read: String::from_utf8_lossy(&named.read_name).into_owned(),
                    kind: "snp",
                    pos: context.pos,
                    qpos: named.qpos,
                    allele: (base_from_index(observation.base_index) as char).to_string(),
                    adjusted_quality: observation.quality,
                    mapq: named.mapq,
                    strand: strand_label(observation.is_reverse),
                });
            }
        }
    }
    (evidence, rows)
}

#[cfg(test)]
fn adjust_fragment_base_observations(observations: &[BaseObservation]) -> Vec<BaseObservation> {
    if observations.len() <= 1 {
        return observations.to_vec();
    }

    let first_base_index = observations[0].base_index;
    if observations
        .iter()
        .all(|observation| observation.base_index == first_base_index)
    {
        observations
            .iter()
            .map(|observation| BaseObservation {
                quality: observation.quality.min(HALF_DEFAULT_PCR_SNV_QUAL),
                ..*observation
            })
            .collect()
    } else {
        observations
            .iter()
            .map(|observation| BaseObservation {
                quality: 0,
                ..*observation
            })
            .collect()
    }
}

fn adjust_named_base_observations(
    observations: &[NamedBaseObservation],
) -> Vec<NamedBaseObservation> {
    if observations.len() <= 1 {
        return observations.to_vec();
    }

    let first_base_index = observations[0].observation.base_index;
    let all_same_base = observations
        .iter()
        .all(|observation| observation.observation.base_index == first_base_index);
    observations
        .iter()
        .map(|observation| {
            let mut adjusted = observation.clone();
            adjusted.observation.quality = if all_same_base {
                adjusted.observation.quality.min(HALF_DEFAULT_PCR_SNV_QUAL)
            } else {
                0
            };
            adjusted
        })
        .collect()
}

fn pileup_indel_evidence(
    pileup: &bam::pileup::Pileup,
    min_baseq: u8,
    min_mapq: u8,
    min_tail_quality: u8,
    dont_use_soft_clipped_bases: bool,
) -> IndelEvidence {
    pileup_indel_evidence_with_rows(
        pileup,
        min_baseq,
        min_mapq,
        min_tail_quality,
        dont_use_soft_clipped_bases,
        None,
    )
    .0
}

fn pileup_indel_evidence_with_rows(
    pileup: &bam::pileup::Pileup,
    min_baseq: u8,
    min_mapq: u8,
    min_tail_quality: u8,
    dont_use_soft_clipped_bases: bool,
    row_context: Option<&ReplayRowContext<'_>>,
) -> (IndelEvidence, Vec<ReplayReadObservationRow>) {
    let mut observations_by_fragment: HashMap<Vec<u8>, Vec<NamedIndelObservation>> = HashMap::new();
    for alignment in pileup.alignments() {
        let record = alignment.record();
        if !read_passes_hc_filter(&record, min_mapq, false) || alignment.is_refskip() {
            continue;
        }
        let Some(qpos) = alignment.qpos() else {
            continue;
        };
        if record
            .qual()
            .get(qpos)
            .is_none_or(|quality| *quality < min_baseq)
        {
            continue;
        }

        // Apply GATK-like read clipping.
        if let Some((clip_start, clip_end)) =
            clip_read_for_evidence(&record, min_tail_quality, dont_use_soft_clipped_bases)
        {
            if qpos < clip_start || qpos >= clip_end {
                continue;
            }
        } else {
            continue;
        }

        let allele = match alignment.indel() {
            Indel::None => Some(IndelObservationAllele::Ref),
            Indel::Ins(len) if len <= MAX_BOOTSTRAP_INDEL_LEN => {
                let inserted = inserted_bases(&record, qpos, len);
                if inserted.is_empty() {
                    None
                } else {
                    Some(IndelObservationAllele::Alt(IndelAllele::Insertion(
                        inserted,
                    )))
                }
            }
            Indel::Del(len) if len <= MAX_BOOTSTRAP_INDEL_LEN => {
                Some(IndelObservationAllele::Alt(IndelAllele::Deletion(len)))
            }
            _ => None,
        };
        let Some(allele) = allele else {
            continue;
        };
        observations_by_fragment
            .entry(record.qname().to_vec())
            .or_default()
            .push(NamedIndelObservation {
                read_name: record.qname().to_vec(),
                qpos,
                mapq: record.mapq(),
                observation: IndelObservation {
                    allele,
                    quality: indel_observation_quality(record.qual()[qpos]),
                    is_reverse: record.is_reverse(),
                },
            });
    }

    let mut evidence = IndelEvidence::default();
    let mut rows = Vec::new();
    for observations in observations_by_fragment.into_values() {
        for named in adjust_named_indel_observations(&observations) {
            let observation = named.observation;
            if observation.quality < min_baseq {
                continue;
            }
            evidence.counts.depth += 1;
            match &observation.allele {
                IndelObservationAllele::Ref => {
                    evidence.counts.ref_count += 1;
                    evidence.ref_strand.increment(observation.is_reverse);
                }
                IndelObservationAllele::Alt(allele) => {
                    *evidence.counts.counts.entry(allele.clone()).or_insert(0) += 1;
                    evidence
                        .alt_strands
                        .entry(allele.clone())
                        .or_default()
                        .increment(observation.is_reverse);
                }
            }
            if let Some(context) = row_context {
                rows.push(ReplayReadObservationRow {
                    region: context.region.to_string(),
                    read: String::from_utf8_lossy(&named.read_name).into_owned(),
                    kind: "indel",
                    pos: context.pos,
                    qpos: named.qpos,
                    allele: indel_observation_allele_label(&observation.allele),
                    adjusted_quality: observation.quality,
                    mapq: named.mapq,
                    strand: strand_label(observation.is_reverse),
                });
            }
            evidence.observations.push(observation);
        }
    }
    (evidence, rows)
}

fn adjust_named_indel_observations(
    observations: &[NamedIndelObservation],
) -> Vec<NamedIndelObservation> {
    observations.to_vec()
}

fn indel_observation_quality(base_quality: u8) -> u8 {
    base_quality.min(DEFAULT_INDEL_QUAL)
}

fn best_snp_call(
    contig: &str,
    pos: u64,
    ref_base: u8,
    evidence: SnpEvidence,
    min_qual: f64,
) -> Option<VariantCall> {
    let ref_index = base_index(ref_base)?;
    let mut best_alt_index = None;
    let mut best_alt_count = 0_u32;
    for (index, count) in evidence.counts.counts.iter().copied().enumerate() {
        if index == ref_index || count <= best_alt_count {
            continue;
        }
        best_alt_index = Some(index);
        best_alt_count = count;
    }
    let alt_index = best_alt_index?;
    if !alt_support_passes(evidence.counts.depth, best_alt_count) {
        return None;
    }
    let model = snp_variant_model(&evidence.observations, ref_index, alt_index);
    if f64::from(model.qual) < min_qual {
        return None;
    }

    Some(VariantCall {
        contig: contig.to_string(),
        pos,
        id: None,
        db: false,
        ref_allele: vec![ref_base],
        alt_allele: vec![base_from_index(alt_index)],
        depth: evidence.counts.depth,
        ref_count: evidence.counts.counts[ref_index],
        alt_count: best_alt_count,
        qual: model.qual,
        fs: fisher_strand_score(evidence.strands[ref_index], evidence.strands[alt_index]),
        pl: model.pl,
        genotype_index: model.genotype_index,
    })
}

fn best_snp_alt(ref_index: Option<usize>, evidence: &SnpEvidence) -> Option<(usize, u32)> {
    let ref_index = ref_index?;
    let mut best_alt_index = None;
    let mut best_alt_count = 0_u32;
    for (index, count) in evidence.counts.counts.iter().copied().enumerate() {
        if index == ref_index || count <= best_alt_count {
            continue;
        }
        best_alt_index = Some(index);
        best_alt_count = count;
    }
    best_alt_index.map(|index| (index, best_alt_count))
}

fn best_indel_call(
    contig: &str,
    pos: u64,
    ref_start: u64,
    ref_bases: &[u8],
    evidence: IndelEvidence,
    min_qual: f64,
) -> Option<VariantCall> {
    let (best_allele, best_alt_count) =
        evidence
            .counts
            .counts
            .iter()
            .max_by(|(allele_a, count_a), (allele_b, count_b)| {
                count_a.cmp(count_b).then_with(|| {
                    indel_allele_sort_key(allele_a).cmp(&indel_allele_sort_key(allele_b))
                })
            })?;
    if !alt_support_passes(evidence.counts.depth, *best_alt_count) {
        return None;
    }
    let model = indel_variant_model(&evidence.observations, best_allele);
    if f64::from(model.qual) < min_qual {
        return None;
    }

    let offset = usize::try_from(pos.checked_sub(ref_start)?).ok()?;
    let anchor = normalize_base(*ref_bases.get(offset)?);
    if !is_acgt(anchor) {
        return None;
    }
    let (ref_allele, alt_allele) = match best_allele {
        IndelAllele::Insertion(inserted) => {
            if inserted.iter().any(|base| !is_acgt(*base)) {
                return None;
            }
            let mut alt_allele = Vec::with_capacity(inserted.len() + 1);
            alt_allele.push(anchor);
            alt_allele.extend_from_slice(inserted);
            (vec![anchor], alt_allele)
        }
        IndelAllele::Deletion(len) => {
            let delete_len = usize::try_from(*len).ok()?;
            let end_offset = offset.checked_add(delete_len)?;
            let deleted = ref_bases.get(offset..=end_offset)?;
            let ref_allele: Vec<u8> = deleted.iter().map(|base| normalize_base(*base)).collect();
            if ref_allele.iter().any(|base| !is_acgt(*base)) {
                return None;
            }
            (ref_allele, vec![anchor])
        }
    };
    let (pos, ref_allele, alt_allele) =
        left_normalize_indel(pos, ref_start, ref_bases, ref_allele, alt_allele);

    Some(VariantCall {
        contig: contig.to_string(),
        pos,
        id: None,
        db: false,
        ref_allele,
        alt_allele,
        depth: evidence.counts.depth,
        ref_count: evidence.counts.ref_count,
        alt_count: *best_alt_count,
        qual: model.qual,
        fs: fisher_strand_score(
            evidence.ref_strand,
            evidence
                .alt_strands
                .get(best_allele)
                .copied()
                .unwrap_or_default(),
        ),
        pl: model.pl,
        genotype_index: model.genotype_index,
    })
}

fn best_indel_alt(evidence: &IndelEvidence) -> Option<(&IndelAllele, &u32)> {
    evidence
        .counts
        .counts
        .iter()
        .max_by(|(allele_a, count_a), (allele_b, count_b)| {
            count_a
                .cmp(count_b)
                .then_with(|| indel_allele_sort_key(allele_a).cmp(&indel_allele_sort_key(allele_b)))
        })
}

fn replay_event_row(region: &str, call: &VariantCall) -> Result<ReplayEventRow> {
    let event_type = if call.ref_allele.len() == 1 && call.alt_allele.len() == 1 {
        "SNP"
    } else {
        "INDEL"
    };
    let ref_allele = allele_string(&call.ref_allele)?.to_string();
    let alt_allele = allele_string(&call.alt_allele)?.to_string();
    let event = format!(
        "{}:{}:{}:{}*,{}",
        call.contig, call.pos, event_type, ref_allele, alt_allele
    );
    let alleles = format!("{}*,{}", ref_allele, alt_allele);
    let raw = format!(
        "{} depth={} ref_count={} alt_count={} qual={} gt={}",
        event,
        call.depth,
        call.ref_count,
        call.alt_count,
        call.qual,
        call.genotype()
    );
    Ok(ReplayEventRow {
        region: region.to_string(),
        event,
        chrom: call.contig.clone(),
        pos: call.pos,
        event_type,
        alleles,
        raw,
        depth: call.depth,
        ref_count: call.ref_count,
        alt_count: call.alt_count,
        qual: call.qual,
        gt: call.genotype(),
    })
}

fn inserted_bases(record: &bam::Record, qpos: usize, len: u32) -> Vec<u8> {
    let len = len as usize;
    let start = qpos.saturating_add(1);
    let end = start.saturating_add(len);
    if end > record.seq_len() {
        return Vec::new();
    }
    (start..end)
        .map(|idx| normalize_base(record.seq()[idx]))
        .collect()
}

fn indel_allele_sort_key(allele: &IndelAllele) -> (u8, u32, Vec<u8>) {
    match allele {
        IndelAllele::Insertion(bases) => (0, bases.len() as u32, bases.clone()),
        IndelAllele::Deletion(len) => (1, *len, Vec::new()),
    }
}

fn left_normalize_indel(
    mut pos: u64,
    ref_start: u64,
    ref_bases: &[u8],
    mut ref_allele: Vec<u8>,
    mut alt_allele: Vec<u8>,
) -> (u64, Vec<u8>, Vec<u8>) {
    if ref_allele.len() == alt_allele.len() {
        return (pos, ref_allele, alt_allele);
    }

    loop {
        let mut changed = false;
        while ref_allele.len() > 1 && alt_allele.len() > 1 && ref_allele.last() == alt_allele.last()
        {
            ref_allele.pop();
            alt_allele.pop();
            changed = true;
        }
        while ref_allele.len() > 1
            && alt_allele.len() > 1
            && ref_allele.first() == alt_allele.first()
        {
            ref_allele.remove(0);
            alt_allele.remove(0);
            pos += 1;
            changed = true;
        }
        if ref_allele.last() == alt_allele.last() && pos > ref_start {
            let Some(prev_base) = reference_base_at(ref_bases, ref_start, pos - 1) else {
                break;
            };
            if !is_acgt(prev_base) {
                break;
            }
            ref_allele.insert(0, prev_base);
            alt_allele.insert(0, prev_base);
            pos -= 1;
            changed = true;
        }
        if !changed {
            break;
        }
    }

    (pos, ref_allele, alt_allele)
}

fn reference_base_at(ref_bases: &[u8], ref_start: u64, pos: u64) -> Option<u8> {
    let offset = usize::try_from(pos.checked_sub(ref_start)?).ok()?;
    ref_bases.get(offset).map(|base| normalize_base(*base))
}

fn indel_observation_allele_label(allele: &IndelObservationAllele) -> String {
    match allele {
        IndelObservationAllele::Ref => "REF".to_string(),
        IndelObservationAllele::Alt(allele) => indel_allele_label(allele),
    }
}

fn indel_allele_label(allele: &IndelAllele) -> String {
    match allele {
        IndelAllele::Insertion(bases) => {
            format!("INS:{}", String::from_utf8_lossy(bases))
        }
        IndelAllele::Deletion(len) => format!("DEL:{len}"),
    }
}

fn strand_label(is_reverse: bool) -> &'static str {
    if is_reverse {
        "-"
    } else {
        "+"
    }
}

fn region_name(contig: &str, start: u64, end: u64) -> String {
    format!("{contig}:{start}-{end}")
}

fn alt_support_passes(depth: u32, alt_count: u32) -> bool {
    depth > 0 && alt_count > 0
}

fn snp_variant_model(
    observations: &[BaseObservation],
    ref_index: usize,
    alt_index: usize,
) -> VariantModel {
    let mut log10_likelihoods = [0.0_f64; 3];
    for observation in observations {
        let error = phred_error_probability(observation.quality);
        let ref_prob = snp_observation_probability(observation.base_index == ref_index, error);
        let alt_prob = snp_observation_probability(observation.base_index == alt_index, error);
        add_diploid_observation(&mut log10_likelihoods, ref_prob, alt_prob);
    }

    let snp_het = 1e-3_f64.log10() - 3.0_f64.log10();
    variant_model_from_log10(log10_likelihoods, [0.0, snp_het, snp_het * 2.0])
}

fn indel_variant_model(
    observations: &[IndelObservation],
    alt_allele: &IndelAllele,
) -> VariantModel {
    let mut log10_likelihoods = [0.0_f64; 3];
    for observation in observations {
        let error = phred_error_probability(observation.quality);
        let ref_prob = match &observation.allele {
            IndelObservationAllele::Ref => 1.0 - error,
            IndelObservationAllele::Alt(_) => error,
        }
        .max(f64::MIN_POSITIVE);
        let alt_prob = match &observation.allele {
            IndelObservationAllele::Alt(allele) if allele == alt_allele => 1.0 - error,
            _ => error,
        }
        .max(f64::MIN_POSITIVE);
        add_diploid_observation(&mut log10_likelihoods, ref_prob, alt_prob);
    }

    let indel_het = (1.0_f64 / 8_000.0).log10();
    variant_model_from_log10(log10_likelihoods, [0.0, indel_het, indel_het * 2.0])
}

fn snp_observation_probability(matches_allele: bool, error: f64) -> f64 {
    if matches_allele {
        1.0 - error
    } else {
        error / 3.0
    }
    .max(f64::MIN_POSITIVE)
}

fn add_diploid_observation(log10_likelihoods: &mut [f64; 3], ref_prob: f64, alt_prob: f64) {
    log10_likelihoods[0] += ref_prob.log10();
    log10_likelihoods[1] += (0.5 * ref_prob + 0.5 * alt_prob)
        .max(f64::MIN_POSITIVE)
        .log10();
    log10_likelihoods[2] += alt_prob.log10();
}

fn phred_error_probability(quality: u8) -> f64 {
    10_f64.powf(-f64::from(quality.min(MAX_SAM_QUAL)) / 10.0)
}

fn variant_model_from_log10(log10_likelihoods: [f64; 3], log10_priors: [f64; 3]) -> VariantModel {
    let log10_posteriors = [
        log10_likelihoods[0] + log10_priors[0],
        log10_likelihoods[1] + log10_priors[1],
        log10_likelihoods[2] + log10_priors[2],
    ];
    let genotype_index = max_index(&log10_posteriors);
    let denominator = log10_sum_exp(&log10_posteriors);
    let ref_posterior_log10 = log10_posteriors[0] - denominator;
    let qual = phred_from_log10_probability(ref_posterior_log10);
    let best_posterior = log10_posteriors[genotype_index];
    let pl = [
        phred_likelihood_delta(log10_posteriors[0], best_posterior),
        phred_likelihood_delta(log10_posteriors[1], best_posterior),
        phred_likelihood_delta(log10_posteriors[2], best_posterior),
    ];

    VariantModel {
        qual,
        pl,
        genotype_index,
    }
}

fn max_index(values: &[f64; 3]) -> usize {
    if values[2] > values[1] && values[2] > values[0] {
        2
    } else if values[1] > values[0] {
        1
    } else {
        0
    }
}

fn log10_sum_exp(values: &[f64]) -> f64 {
    let max = values
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, |acc, value| acc.max(value));
    if !max.is_finite() {
        return max;
    }
    let sum: f64 = values.iter().map(|value| 10_f64.powf(value - max)).sum();
    max + sum.log10()
}

fn phred_from_log10_probability(log10_probability: f64) -> u32 {
    if log10_probability <= -999.9 {
        9_999
    } else {
        (-10.0 * log10_probability).round().clamp(0.0, 9_999.0) as u32
    }
}

fn phred_likelihood_delta(log10_value: f64, best_log10_value: f64) -> u32 {
    (-10.0 * (log10_value - best_log10_value))
        .round()
        .clamp(0.0, 999.0) as u32
}

fn fisher_strand_score(ref_strand: StrandCounts, alt_strand: StrandCounts) -> f64 {
    let pvalue = fisher_exact_pvalue(ref_strand, alt_strand);
    -10.0 * pvalue.max(FISHER_STRAND_MIN_PVALUE).log10()
}

fn fisher_exact_pvalue(ref_strand: StrandCounts, alt_strand: StrandCounts) -> f64 {
    let table = normalize_fisher_table([
        [ref_strand.forward, ref_strand.reverse],
        [alt_strand.forward, alt_strand.reverse],
    ]);
    let row_ref = table[0][0] + table[0][1];
    let row_alt = table[1][0] + table[1][1];
    let col_forward = table[0][0] + table[1][0];
    let total = row_ref + row_alt;
    let lo = col_forward.saturating_sub(row_alt);
    let hi = col_forward.min(row_ref);
    if hi <= lo {
        return 1.0;
    }

    let observed = table[0][0];
    let observed_log = hypergeometric_log_probability(observed, total, row_ref, col_forward);
    let mut pvalue = 0.0;
    for value in lo..=hi {
        let logp = hypergeometric_log_probability(value, total, row_ref, col_forward);
        if logp <= observed_log + 1e-12 {
            pvalue += logp.exp();
        }
    }
    pvalue.min(1.0)
}

fn normalize_fisher_table(mut table: [[u32; 2]; 2]) -> [[u32; 2]; 2] {
    let sum = table[0][0] + table[0][1] + table[1][0] + table[1][1];
    if f64::from(sum) <= FISHER_STRAND_TARGET_TABLE_SIZE * 2.0 {
        return table;
    }
    let factor = f64::from(sum) / FISHER_STRAND_TARGET_TABLE_SIZE;
    for row in &mut table {
        for value in row {
            *value = (f64::from(*value) / factor) as u32;
        }
    }
    table
}

fn hypergeometric_log_probability(k: u32, population: u32, success_states: u32, draws: u32) -> f64 {
    log_choose(success_states, k) + log_choose(population - success_states, draws - k)
        - log_choose(population, draws)
}

fn log_choose(n: u32, k: u32) -> f64 {
    if k > n {
        return f64::NEG_INFINITY;
    }
    log_factorial(n) - log_factorial(k) - log_factorial(n - k)
}

fn log_factorial(n: u32) -> f64 {
    (2..=n).map(|value| f64::from(value).ln()).sum()
}

fn read_passes_hc_filter(record: &bam::Record, min_mapq: u8, exclude_supplementary: bool) -> bool {
    const UNMAPPED: u16 = 0x4;
    const SECONDARY: u16 = 0x100;
    const QCFAIL: u16 = 0x200;
    const DUPLICATE: u16 = 0x400;
    const SUPPLEMENTARY: u16 = 0x800;

    let mut excluded = UNMAPPED | SECONDARY | QCFAIL | DUPLICATE;
    if exclude_supplementary {
        excluded |= SUPPLEMENTARY;
    }

    record.flags() & excluded == 0 && record.mapq() >= min_mapq && cigar_has_reference_bases(record)
}

fn cigar_has_reference_bases(record: &bam::Record) -> bool {
    record.cigar().iter().any(|op| {
        matches!(
            op,
            Cigar::Match(_) | Cigar::Equal(_) | Cigar::Diff(_) | Cigar::Del(_) | Cigar::RefSkip(_)
        )
    })
}

fn position_is_requested(intervals: &[Interval], pos1: u64, cursor: &mut usize) -> bool {
    while *cursor < intervals.len() && intervals[*cursor].end < pos1 {
        *cursor += 1;
    }
    intervals
        .get(*cursor)
        .is_some_and(|interval| interval.start <= pos1 && pos1 <= interval.end)
}

fn sort_intervals(intervals: &mut [Interval], dict: &SequenceDict) -> Result<()> {
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
    Ok(())
}

fn coalesce_fetch_windows(intervals: &[Interval]) -> Vec<FetchWindow> {
    let mut windows: Vec<FetchWindow> = Vec::new();
    for interval in intervals {
        if let Some(current) = windows.last_mut() {
            let same_contig = current.contig == interval.contig;
            let close_enough = interval.start <= current.end.saturating_add(FETCH_WINDOW_GAP + 1);
            let merged_len = interval.end.saturating_sub(current.start).saturating_add(1);
            if same_contig && close_enough && merged_len <= FETCH_WINDOW_MAX_BASES {
                current.end = current.end.max(interval.end);
                current.intervals.push(interval.clone());
                continue;
            }
        }
        windows.push(FetchWindow {
            contig: interval.contig.clone(),
            start: interval.start,
            end: interval.end,
            intervals: vec![interval.clone()],
        });
    }
    windows
}

fn partition_fetch_windows_by_bases(
    windows: &[FetchWindow],
    threads: usize,
) -> Vec<Vec<FetchWindow>> {
    if windows.is_empty() {
        return Vec::new();
    }
    let workers = threads.min(windows.len()).max(1);
    let total_bases: u64 = windows.iter().map(FetchWindow::len).sum();
    let target_bases = total_bases.div_ceil(workers as u64).max(1);

    let mut partitions = Vec::with_capacity(workers);
    let mut current = Vec::new();
    let mut current_bases = 0_u64;
    for window in windows {
        if !current.is_empty() && partitions.len() + 1 < workers && current_bases >= target_bases {
            partitions.push(current);
            current = Vec::new();
            current_bases = 0;
        }
        current_bases += window.len();
        current.push(window.clone());
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
        index_by_name.insert(name.clone(), records.len());
        records.push(DictRecord { name, length });
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

fn annotate_dbsnp(path: &Path, variants: &mut [VariantCall]) -> Result<()> {
    if variants.is_empty() {
        return Ok(());
    }
    let mut reader = tbx::Reader::from_path(path)
        .with_context(|| format!("opening dbSNP {}", path.display()))?;
    let mut tid_cache: HashMap<String, Option<u64>> = HashMap::new();
    let mut record = Vec::new();
    for variant in variants {
        let tid = match tid_cache.get(&variant.contig) {
            Some(tid) => *tid,
            None => {
                let tid = reader.tid(&variant.contig).ok();
                tid_cache.insert(variant.contig.clone(), tid);
                tid
            }
        };
        let Some(tid) = tid else {
            continue;
        };
        let start0 = variant.pos.saturating_sub(1);
        let end0 = variant
            .pos
            .saturating_add(variant.ref_allele.len() as u64)
            .max(start0 + 1);
        if reader.fetch(tid, start0, end0).is_err() {
            continue;
        }
        while TbxRead::read(&mut reader, &mut record)
            .with_context(|| format!("reading dbSNP {}", path.display()))?
        {
            if dbsnp_record_matches(&record, variant)? {
                let id = dbsnp_record_id(&record)?;
                if !id.is_empty() {
                    variant.id = Some(id);
                }
                variant.db = true;
                break;
            }
        }
    }
    Ok(())
}

fn dbsnp_record_matches(record: &[u8], variant: &VariantCall) -> Result<bool> {
    if record.starts_with(b"#") {
        return Ok(false);
    }
    let line = std::str::from_utf8(record).context("dbSNP record is not UTF-8")?;
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 5 {
        return Ok(false);
    }
    let pos = fields[1]
        .parse::<u64>()
        .with_context(|| format!("invalid dbSNP POS value '{}'", fields[1]))?;
    if fields[0] != variant.contig || pos != variant.pos {
        return Ok(false);
    }
    if fields[3].as_bytes() != variant.ref_allele {
        return Ok(false);
    }
    Ok(fields[4]
        .split(',')
        .any(|alt| alt.as_bytes() == variant.alt_allele))
}

fn dbsnp_record_id(record: &[u8]) -> Result<String> {
    let line = std::str::from_utf8(record).context("dbSNP record is not UTF-8")?;
    let Some(id) = line.split('\t').nth(2) else {
        return Ok(String::new());
    };
    if id == "." {
        Ok(String::new())
    } else {
        Ok(id.to_string())
    }
}

fn sort_variant_calls(variants: &mut [VariantCall], dict: &SequenceDict) -> Result<()> {
    for variant in variants.iter() {
        if dict.order(&variant.contig).is_none() {
            bail!(
                "contig '{}' is not present in the sequence dictionary",
                variant.contig
            );
        }
    }
    variants.sort_by(|a, b| {
        dict.order(&a.contig)
            .cmp(&dict.order(&b.contig))
            .then(a.pos.cmp(&b.pos))
            .then(a.ref_allele.cmp(&b.ref_allele))
            .then(a.alt_allele.cmp(&b.alt_allele))
    });
    Ok(())
}

fn dedup_variant_calls(variants: &mut Vec<VariantCall>) {
    variants.dedup_by(|a, b| {
        a.contig == b.contig
            && a.pos == b.pos
            && a.ref_allele == b.ref_allele
            && a.alt_allele == b.alt_allele
    });
}

fn write_bootstrap_vcf(
    path: &Path,
    config: &HaplotypeCallerConfig,
    dict: &SequenceDict,
    sample_name: &str,
    variants: &[VariantCall],
) -> Result<()> {
    create_parent_dir(path)?;
    let mut writer = open_vcf_writer(path)?;
    writeln!(writer, "##fileformat=VCFv4.2")?;
    writeln!(writer, "##source=rust_haplotype_caller_pipeline_bootstrap")?;
    writeln!(
        writer,
        "##rust_hc_pipeline_bootstrap=Pipeline-only pileup SNP/short-indel caller; local assembly, PairHMM, phasing, and GVCF are not implemented"
    )?;
    writeln!(writer, "##reference={}", config.reference.display())?;
    for record in &dict.records {
        writeln!(
            writer,
            "##contig=<ID={},length={}>",
            record.name, record.length
        )?;
    }
    writeln!(
        writer,
        "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">"
    )?;
    writeln!(
        writer,
        "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency in genotypes\">"
    )?;
    writeln!(writer, "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">")?;
    writeln!(writer, "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Filtered fragment depth used by the pipeline-only Rust caller\">")?;
    writeln!(
        writer,
        "##INFO=<ID=DB,Number=0,Type=Flag,Description=\"dbSNP exact REF/ALT match\">"
    )?;
    writeln!(writer, "##INFO=<ID=FS,Number=1,Type=Float,Description=\"FisherStrand-style strand-bias phred score from fragment evidence\">")?;
    writeln!(writer, "##INFO=<ID=QD,Number=1,Type=Float,Description=\"QUAL divided by AD depth for the variant sample\">")?;
    writeln!(
        writer,
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    )?;
    writeln!(writer, "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for ref and alt alleles\">")?;
    writeln!(
        writer,
        "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Filtered read depth\">"
    )?;
    writeln!(
        writer,
        "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality from diploid fragment likelihoods\">"
    )?;
    writeln!(writer, "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Normalized phred-scaled genotype likelihoods from fragment evidence\">")?;
    writeln!(
        writer,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}",
        sample_name
    )?;
    for variant in variants {
        write_variant_record(&mut writer, variant)?;
    }
    writer.flush()?;
    Ok(())
}

fn sort_replay_rows(output: &mut ReplayWorkerOutput, dict: &SequenceDict) -> Result<()> {
    for row in &output.active_regions {
        if dict.order(&row.contig).is_none() {
            bail!(
                "contig '{}' is not present in the sequence dictionary",
                row.contig
            );
        }
    }
    for row in &output.active_loci {
        if dict.order(&row.contig).is_none() {
            bail!(
                "contig '{}' is not present in the sequence dictionary",
                row.contig
            );
        }
    }
    for row in &output.events {
        if dict.order(&row.chrom).is_none() {
            bail!(
                "contig '{}' is not present in the sequence dictionary",
                row.chrom
            );
        }
    }
    output.active_regions.sort_by(|a, b| {
        dict.order(&a.contig)
            .cmp(&dict.order(&b.contig))
            .then(a.start.cmp(&b.start))
            .then(a.end.cmp(&b.end))
    });
    output.active_loci.sort_by(|a, b| {
        dict.order(&a.contig)
            .cmp(&dict.order(&b.contig))
            .then(a.pos.cmp(&b.pos))
            .then(a.region.cmp(&b.region))
    });
    output.read_observations.sort_by(|a, b| {
        a.region
            .cmp(&b.region)
            .then(a.pos.cmp(&b.pos))
            .then(a.read.cmp(&b.read))
            .then(a.kind.cmp(&b.kind))
            .then(a.qpos.cmp(&b.qpos))
    });
    output.events.sort_by(|a, b| {
        dict.order(&a.chrom)
            .cmp(&dict.order(&b.chrom))
            .then(a.pos.cmp(&b.pos))
            .then(a.event_type.cmp(&b.event_type))
            .then(a.alleles.cmp(&b.alleles))
    });
    output.haplotypes.sort_by(|a, b| {
        a.region
            .cmp(&b.region)
            .then(a.stage.cmp(&b.stage))
            .then(a.haplotype.cmp(&b.haplotype))
    });
    output.pairhmms.sort_by(|a, b| {
        a.region
            .cmp(&b.region)
            .then(a.read_index.cmp(&b.read_index))
            .then(a.haplotype.cmp(&b.haplotype))
    });
    Ok(())
}

fn write_replay_tables(
    prefix: &Path,
    output: &ReplayWorkerOutput,
    genotype_rows: &[ReplayGenotypeRow],
) -> Result<()> {
    write_replay_active_regions(&replay_prefixed_path(prefix, "active_regions.tsv"), output)?;
    write_replay_active_loci(&replay_prefixed_path(prefix, "active_loci.tsv"), output)?;
    write_replay_read_observations(
        &replay_prefixed_path(prefix, "read_observations.tsv"),
        output,
    )?;
    write_replay_events(&replay_prefixed_path(prefix, "events.tsv"), output)?;
    write_replay_genotypes(
        &replay_prefixed_path(prefix, "genotypes.tsv"),
        genotype_rows,
    )?;
    write_replay_haplotypes(&replay_prefixed_path(prefix, "haplotypes.tsv"), output)?;
    write_replay_pairhmms(&replay_prefixed_path(prefix, "pairhmm.tsv"), output)?;
    write_empty_allele_likelihoods(prefix)?;
    Ok(())
}

fn write_replay_active_regions(path: &Path, output: &ReplayWorkerOutput) -> Result<()> {
    create_parent_dir(path)?;
    let mut writer =
        BufWriter::new(File::create(path).with_context(|| format!("creating {}", path.display()))?);
    writeln!(
        writer,
        "contig\tstart\tend\tregion\tobserved_loci\tactive_loci\tcandidate_events\tmax_alt_fraction\tmean_alt_fraction"
    )?;
    for row in &output.active_regions {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.6}\t{:.6}",
            row.contig,
            row.start,
            row.end,
            row.region,
            row.observed_loci,
            row.active_loci,
            row.candidate_events,
            row.max_alt_fraction,
            row.mean_alt_fraction
        )?;
    }
    Ok(())
}

fn write_replay_active_loci(path: &Path, output: &ReplayWorkerOutput) -> Result<()> {
    create_parent_dir(path)?;
    let mut writer =
        BufWriter::new(File::create(path).with_context(|| format!("creating {}", path.display()))?);
    writeln!(
        writer,
        "contig\tpos\tregion\tref_base\tdepth\tsnp_alt_count\tsnp_best_alt\tindel_alt_count\tindel_best_alt\talt_fraction\tactive_probability_proxy"
    )?;
    for row in &output.active_loci {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.6}\t{:.6}",
            row.contig,
            row.pos,
            row.region,
            row.ref_base as char,
            row.depth,
            row.snp_alt_count,
            row.snp_best_alt,
            row.indel_alt_count,
            row.indel_best_alt,
            row.alt_fraction,
            row.active_probability_proxy
        )?;
    }
    Ok(())
}

fn write_replay_read_observations(path: &Path, output: &ReplayWorkerOutput) -> Result<()> {
    create_parent_dir(path)?;
    let mut writer =
        BufWriter::new(File::create(path).with_context(|| format!("creating {}", path.display()))?);
    writeln!(
        writer,
        "region\tread\tkind\tpos\tqpos\tallele\tadjusted_quality\tmapq\tstrand"
    )?;
    for row in &output.read_observations {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            row.region,
            row.read,
            row.kind,
            row.pos,
            row.qpos,
            row.allele,
            row.adjusted_quality,
            row.mapq,
            row.strand
        )?;
    }
    Ok(())
}

fn write_replay_events(path: &Path, output: &ReplayWorkerOutput) -> Result<()> {
    create_parent_dir(path)?;
    let mut writer =
        BufWriter::new(File::create(path).with_context(|| format!("creating {}", path.display()))?);
    writeln!(
        writer,
        "region\tevent\tchrom\tpos\ttype\talleles\traw\tdepth\tref_count\talt_count\tqual\tgt"
    )?;
    for row in &output.events {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            row.region,
            row.event,
            row.chrom,
            row.pos,
            row.event_type,
            row.alleles,
            row.raw,
            row.depth,
            row.ref_count,
            row.alt_count,
            row.qual,
            row.gt
        )?;
    }
    Ok(())
}

fn write_replay_genotypes(path: &Path, rows: &[ReplayGenotypeRow]) -> Result<()> {
    create_parent_dir(path)?;
    let mut writer =
        BufWriter::new(File::create(path).with_context(|| format!("creating {}", path.display()))?);
    writeln!(
        writer,
        "chrom\tpos\tref\talt\tqual\tfilter\tgt\tgq\tdp\tad_ref\tad_alt\tfs\tqd\tpl\tdb"
    )?;
    for row in rows {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.3}\t{:.2}\t{}\t{}",
            row.chrom,
            row.pos,
            row.ref_allele,
            row.alt,
            row.qual,
            row.filter,
            row.gt,
            row.gq,
            row.dp,
            row.ad_ref,
            row.ad_alt,
            row.fs,
            row.qd,
            row.pl,
            row.db
        )?;
    }
    Ok(())
}
fn write_replay_haplotypes(path: &Path, output: &ReplayWorkerOutput) -> Result<()> {
    create_parent_dir(path)?;
    let mut writer =
        BufWriter::new(File::create(path).with_context(|| format!("creating {}", path.display()))?);
    writeln!(
        writer,
        "region\tstage\thaplotype\tspan_start\tspan_end\tkmer\tlength\tcigar\tis_ref\tbases"
    )?;
    for row in &output.haplotypes {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            row.region,
            row.stage,
            row.haplotype,
            row.span_start,
            row.span_end,
            row.kmer,
            row.length,
            row.cigar,
            row.is_ref,
            row.bases
        )?;
    }
    Ok(())
}

#[derive(Clone, Debug, Eq, PartialEq, Hash)]
struct LocalHaplotype {
    bases: Vec<u8>,
    is_ref: bool,
    cigar: String,
    event_indices: Vec<usize>,
}

fn push_cigar(cigar: &mut Vec<(u32, char)>, len: u32, op: char) {
    if len == 0 {
        return;
    }
    if let Some(last) = cigar.last_mut() {
        if last.1 == op {
            last.0 += len;
            return;
        }
    }
    cigar.push((len, op));
}

fn format_cigar(cigar: &[(u32, char)]) -> String {
    let mut s = String::with_capacity(cigar.len() * 4);
    for &(len, op) in cigar {
        s.push_str(&len.to_string());
        s.push(op);
    }
    s
}

fn build_local_haplotypes(
    _contig: &str,
    region_start: u64,
    region_end: u64,
    ref_bases: &[u8],
    candidate_events: &[VariantCall],
    max_haplotypes: usize,
) -> Vec<LocalHaplotype> {
    let ref_hap = LocalHaplotype {
        bases: ref_bases.to_vec(),
        is_ref: true,
        cigar: format!("{}M", ref_bases.len()),
        event_indices: Vec::new(),
    };

    let mut haplotypes = vec![ref_hap];

    let mut valid_events = Vec::new();
    for event in candidate_events {
        let event_end = event.pos + event.ref_allele.len() as u64 - 1;
        if event.pos >= region_start && event_end <= region_end {
            valid_events.push(event);
        }
    }

    if valid_events.len() > 7 {
        valid_events.truncate(7);
    }

    let n_events = valid_events.len();
    for mask in 1..(1 << n_events) {
        if haplotypes.len() >= max_haplotypes {
            break;
        }

        let mut overlap = false;
        let mut last_end = 0;
        let mut selected_events = Vec::new();
        let mut event_indices = Vec::new();
        for i in 0..n_events {
            if (mask & (1 << i)) != 0 {
                let ev = valid_events[i];
                if ev.pos <= last_end {
                    overlap = true;
                    break;
                }
                last_end = ev.pos + ev.ref_allele.len() as u64 - 1;
                selected_events.push(ev);
                event_indices.push(i);
            }
        }

        if overlap {
            continue;
        }

        let mut bases = Vec::with_capacity(ref_bases.len());
        let mut cigar_ops = Vec::new();
        let mut ref_offset = 0;
        let mut current_pos = region_start;

        for ev in &selected_events {
            let dist = ev.pos.saturating_sub(current_pos) as usize;
            if dist > 0 {
                bases.extend_from_slice(&ref_bases[ref_offset..ref_offset + dist]);
                push_cigar(&mut cigar_ops, dist as u32, 'M');
                ref_offset += dist;
                current_pos += dist as u64;
            }

            bases.extend_from_slice(&ev.alt_allele);
            let match_len = (ev.alt_allele.len() as u32).min(ev.ref_allele.len() as u32);
            push_cigar(&mut cigar_ops, match_len, 'M');
            
            if ev.alt_allele.len() > ev.ref_allele.len() {
                push_cigar(&mut cigar_ops, (ev.alt_allele.len() - ev.ref_allele.len()) as u32, 'I');
            } else if ev.ref_allele.len() > ev.alt_allele.len() {
                push_cigar(&mut cigar_ops, (ev.ref_allele.len() - ev.alt_allele.len()) as u32, 'D');
            }

            ref_offset += ev.ref_allele.len();
            current_pos += ev.ref_allele.len() as u64;
        }

        let rem = ref_bases.len().saturating_sub(ref_offset);
        if rem > 0 {
            bases.extend_from_slice(&ref_bases[ref_offset..]);
            push_cigar(&mut cigar_ops, rem as u32, 'M');
        }

        haplotypes.push(LocalHaplotype {
            bases,
            is_ref: false,
            cigar: format_cigar(&cigar_ops),
            event_indices,
        });
    }

    haplotypes
}

fn write_replay_pairhmms(path: &Path, output: &ReplayWorkerOutput) -> Result<()> {
    create_parent_dir(path)?;
    let mut writer =
        BufWriter::new(File::create(path).with_context(|| format!("creating {}", path.display()))?);
    writeln!(
        writer,
        "region\tread\thaplotype\tread_index\tcigar\tmapq\tloc\tunclipped_loc\tlength\tscore"
    )?;
    for row in &output.pairhmms {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            row.region,
            row.read,
            row.haplotype,
            row.read_index,
            row.cigar,
            row.mapq,
            row.loc,
            row.unclipped_loc,
            row.length,
            row.score
        )?;
    }
    Ok(())
}

fn write_empty_allele_likelihoods(prefix: &Path) -> Result<()> {
    let mut allele_likelihoods = BufWriter::new(File::create(replay_prefixed_path(
        prefix,
        "allele_likelihoods.tsv",
    ))?);
    writeln!(
        allele_likelihoods,
        "region\tevent\tmatrix\tread\tread_index\tallele\tscore"
    )?;
    Ok(())
}

fn replay_prefixed_path(prefix: &Path, suffix: &str) -> PathBuf {
    PathBuf::from(format!("{}.{}", prefix.display(), suffix))
}

fn open_vcf_writer(path: &Path) -> Result<Box<dyn Write>> {
    if is_gzip_path(path) {
        let writer = bgzf::Writer::from_path(path)
            .with_context(|| format!("creating bgzipped VCF {}", path.display()))?;
        Ok(Box::new(writer))
    } else {
        let file = File::create(path).with_context(|| format!("creating {}", path.display()))?;
        Ok(Box::new(BufWriter::new(file)))
    }
}

fn write_variant_record(writer: &mut dyn Write, variant: &VariantCall) -> Result<()> {
    let ac = variant.alt_allele_count();
    let af = f64::from(ac) / 2.0;
    let qd = if variant.depth == 0 {
        0.0
    } else {
        fix_too_high_qd(f64::from(variant.qual) / f64::from(variant.depth))
    };
    let gt = variant.genotype();
    let pl = genotype_likelihoods(variant);
    let ref_allele = allele_string(&variant.ref_allele)?;
    let alt_allele = allele_string(&variant.alt_allele)?;
    let id = variant.id.as_deref().unwrap_or(".");
    let db = if variant.db { ";DB" } else { "" };
    writeln!(
        writer,
        "{}\t{}\t{}\t{}\t{}\t{}\tPASS\tAC={};AF={:.3};AN=2;DP={}{};FS={:.3};QD={:.2}\tGT:AD:DP:GQ:PL\t{}:{},{}:{}:{}:{}",
        variant.contig,
        variant.pos,
        id,
        ref_allele,
        alt_allele,
        variant.qual,
        ac,
        af,
        variant.depth,
        db,
        variant.fs,
        qd,
        gt,
        variant.ref_count,
        variant.alt_count,
        variant.depth,
        variant.gq(),
        pl
    )?;
    Ok(())
}

fn genotype_likelihoods(variant: &VariantCall) -> String {
    format!("{},{},{}", variant.pl[0], variant.pl[1], variant.pl[2])
}

fn fix_too_high_qd(qd: f64) -> f64 {
    if qd < 35.0 {
        qd
    } else {
        30.0
    }
}

fn allele_string(allele: &[u8]) -> Result<&str> {
    std::str::from_utf8(allele).context("allele contains non-UTF-8 bases")
}

fn sample_name_from_bam(path: &Path) -> Result<String> {
    let reader =
        bam::Reader::from_path(path).with_context(|| format!("opening {}", path.display()))?;
    let header = String::from_utf8_lossy(reader.header().as_bytes());
    for line in header.lines() {
        if !line.starts_with("@RG\t") {
            continue;
        }
        for field in line.split('\t').skip(1) {
            if let Some(sample) = field.strip_prefix("SM:") {
                if !sample.is_empty() {
                    return Ok(sample.to_string());
                }
            }
        }
    }
    Ok("SAMPLE".to_string())
}

fn is_gzip_path(path: &Path) -> bool {
    path.extension().is_some_and(|extension| extension == "gz")
}

fn write_tabix_index(path: &Path, threads: usize) -> Result<()> {
    let c_path = CString::new(path.as_os_str().as_bytes())?;
    let thread_count = threads.clamp(1, i32::MAX as usize) as i32;
    let result = unsafe {
        htslib::tbx_index_build3(
            c_path.as_ptr(),
            ptr::null(),
            0,
            thread_count,
            &htslib::tbx_conf_vcf,
        )
    };
    if result != 0 {
        bail!("failed to create tabix index for {}", path.display());
    }
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

fn base_index(base: u8) -> Option<usize> {
    match base {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

fn base_from_index(index: usize) -> u8 {
    match index {
        0 => b'A',
        1 => b'C',
        2 => b'G',
        3 => b'T',
        _ => unreachable!("invalid base index"),
    }
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

    fn snp_evidence(
        ref_index: usize,
        ref_count: u32,
        alt_index: usize,
        alt_count: u32,
        quality: u8,
    ) -> SnpEvidence {
        let mut evidence = SnpEvidence::default();
        for _ in 0..ref_count {
            evidence.counts.counts[ref_index] += 1;
            evidence.counts.depth += 1;
            evidence.strands[ref_index].increment(false);
            evidence.observations.push(BaseObservation {
                base_index: ref_index,
                quality,
                is_reverse: false,
            });
        }
        for idx in 0..alt_count {
            let is_reverse = idx % 2 == 1;
            evidence.counts.counts[alt_index] += 1;
            evidence.counts.depth += 1;
            evidence.strands[alt_index].increment(is_reverse);
            evidence.observations.push(BaseObservation {
                base_index: alt_index,
                quality,
                is_reverse,
            });
        }
        evidence
    }

    fn indel_evidence(
        ref_count: u32,
        alt_allele: IndelAllele,
        alt_count: u32,
        quality: u8,
    ) -> IndelEvidence {
        let mut evidence = IndelEvidence::default();
        evidence.counts.ref_count = ref_count;
        evidence.counts.depth = ref_count + alt_count;
        evidence.counts.counts.insert(alt_allele.clone(), alt_count);
        for _ in 0..ref_count {
            evidence.ref_strand.increment(false);
            evidence.observations.push(IndelObservation {
                allele: IndelObservationAllele::Ref,
                quality,
                is_reverse: false,
            });
        }
        for idx in 0..alt_count {
            let is_reverse = idx % 2 == 1;
            evidence
                .alt_strands
                .entry(alt_allele.clone())
                .or_default()
                .increment(is_reverse);
            evidence.observations.push(IndelObservation {
                allele: IndelObservationAllele::Alt(alt_allele.clone()),
                quality,
                is_reverse,
            });
        }
        evidence
    }

    fn test_variant(contig: &str, pos: u64, ref_allele: &[u8], alt_allele: &[u8]) -> VariantCall {
        VariantCall {
            contig: contig.to_string(),
            pos,
            id: None,
            db: false,
            ref_allele: ref_allele.to_vec(),
            alt_allele: alt_allele.to_vec(),
            depth: 10,
            ref_count: 8,
            alt_count: 2,
            qual: 20,
            fs: 0.0,
            pl: [20, 0, 20],
            genotype_index: 1,
        }
    }

    #[test]
    fn intervals_sort_by_dictionary_order() {
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
        sort_intervals(&mut intervals, &dict).unwrap();
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
                    end: 80,
                },
                Interval {
                    contig: "chr1".to_string(),
                    start: 81,
                    end: 90,
                },
            ]
        );
    }

    #[test]
    fn fetch_windows_coalesce_nearby_intervals() {
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
                start: 2_000,
                end: 2_010,
            },
        ];
        let windows = coalesce_fetch_windows(&intervals);
        assert_eq!(windows.len(), 2);
        assert_eq!(windows[0].start, 1);
        assert_eq!(windows[0].end, 30);
        assert_eq!(windows[0].intervals.len(), 2);
        assert_eq!(windows[1].start, 2_000);
    }

    #[test]
    fn fetch_window_partitioning_preserves_all_bases() {
        let windows = vec![
            FetchWindow {
                contig: "chr1".to_string(),
                start: 1,
                end: 10,
                intervals: Vec::new(),
            },
            FetchWindow {
                contig: "chr1".to_string(),
                start: 11,
                end: 30,
                intervals: Vec::new(),
            },
            FetchWindow {
                contig: "chr1".to_string(),
                start: 31,
                end: 60,
                intervals: Vec::new(),
            },
        ];
        let total_bases: u64 = windows.iter().map(FetchWindow::len).sum();
        let partitions = partition_fetch_windows_by_bases(&windows, 2);
        let partition_bases: u64 = partitions.iter().flatten().map(FetchWindow::len).sum();
        assert_eq!(partition_bases, total_bases);
        assert_eq!(
            partitions.iter().map(Vec::len).sum::<usize>(),
            windows.len()
        );
    }

    #[test]
    fn requested_position_uses_sorted_interval_cursor() {
        let intervals = vec![
            Interval {
                contig: "chr1".to_string(),
                start: 10,
                end: 20,
            },
            Interval {
                contig: "chr1".to_string(),
                start: 30,
                end: 40,
            },
        ];
        let mut cursor = 0;
        assert!(!position_is_requested(&intervals, 9, &mut cursor));
        assert!(position_is_requested(&intervals, 10, &mut cursor));
        assert!(position_is_requested(&intervals, 20, &mut cursor));
        assert!(!position_is_requested(&intervals, 25, &mut cursor));
        assert!(position_is_requested(&intervals, 35, &mut cursor));
    }

    #[test]
    fn snp_call_uses_likelihood_quality_threshold() {
        let evidence = snp_evidence(0, 5, 1, 5, 30);
        assert!(best_snp_call("chr1", 10, b'A', evidence.clone(), 200.0).is_none());
        let call = best_snp_call("chr1", 10, b'A', evidence, 20.0).unwrap();
        assert_eq!(call.ref_allele, b"A");
        assert_eq!(call.alt_allele, b"C");
        assert!(call.qual >= 90);
        assert_eq!(call.genotype(), "0/1");
        assert_eq!(call.pl[1], 0);
    }

    #[test]
    fn snp_call_uses_hom_alt_when_likelihoods_overcome_prior() {
        let evidence = snp_evidence(0, 0, 1, 20, 30);
        let call = best_snp_call("chr1", 10, b'A', evidence, 20.0).unwrap();
        assert_eq!(call.genotype(), "1/1");
        assert_eq!(call.alt_allele_count(), 2);
        assert_eq!(call.pl[2], 0);
    }

    #[test]
    fn snp_call_rejects_low_alt_fraction_noise_after_likelihoods() {
        let evidence = snp_evidence(0, 29, 1, 2, 30);
        assert!(best_snp_call("chr1", 10, b'A', evidence, 20.0).is_none());
    }

    #[test]
    fn overlapping_same_base_fragment_caps_both_observations() {
        let observations = vec![
            BaseObservation {
                base_index: 1,
                quality: 35,
                is_reverse: false,
            },
            BaseObservation {
                base_index: 1,
                quality: 33,
                is_reverse: true,
            },
        ];
        let kept = adjust_fragment_base_observations(&observations);
        assert_eq!(kept.len(), 2);
        assert!(kept.iter().all(|observation| observation.base_index == 1));
        assert!(kept
            .iter()
            .all(|observation| observation.quality == HALF_DEFAULT_PCR_SNV_QUAL));
    }

    #[test]
    fn overlapping_discordant_fragment_zeroes_both_observations() {
        let observations = vec![
            BaseObservation {
                base_index: 1,
                quality: 35,
                is_reverse: false,
            },
            BaseObservation {
                base_index: 2,
                quality: 33,
                is_reverse: true,
            },
        ];
        let kept = adjust_fragment_base_observations(&observations);
        assert_eq!(kept.len(), 2);
        assert!(kept.iter().all(|observation| observation.quality == 0));
    }

    #[test]
    fn fisher_strand_score_detects_one_sided_alt_support() {
        let fs = fisher_strand_score(
            StrandCounts {
                forward: 10,
                reverse: 0,
            },
            StrandCounts {
                forward: 0,
                reverse: 10,
            },
        );
        assert!(fs > 30.0);
    }

    #[test]
    fn insertion_call_uses_left_anchor() {
        let evidence = indel_evidence(6, IndelAllele::Insertion(b"TG".to_vec()), 6, 30);
        let call = best_indel_call("chr1", 100, 100, b"ACGT", evidence, 20.0).unwrap();
        assert_eq!(call.pos, 100);
        assert_eq!(call.ref_allele, b"A");
        assert_eq!(call.alt_allele, b"ATG");
        assert_eq!(call.ref_count, 6);
        assert_eq!(call.alt_count, 6);
    }

    #[test]
    fn deletion_call_uses_left_anchor_and_deleted_reference() {
        let evidence = indel_evidence(6, IndelAllele::Deletion(2), 6, 30);
        let call = best_indel_call("chr1", 100, 100, b"ACGT", evidence, 20.0).unwrap();
        assert_eq!(call.pos, 100);
        assert_eq!(call.ref_allele, b"ACG");
        assert_eq!(call.alt_allele, b"A");
    }

    #[test]
    fn indel_normalization_left_aligns_homopolymer_insertion() {
        let (pos, ref_allele, alt_allele) =
            left_normalize_indel(102, 100, b"ATTTG", b"T".to_vec(), b"TT".to_vec());
        assert_eq!(pos, 100);
        assert_eq!(ref_allele, b"A");
        assert_eq!(alt_allele, b"AT");
    }

    #[test]
    fn indel_normalization_left_aligns_homopolymer_deletion() {
        let (pos, ref_allele, alt_allele) =
            left_normalize_indel(102, 100, b"ATTTG", b"TT".to_vec(), b"T".to_vec());
        assert_eq!(pos, 100);
        assert_eq!(ref_allele, b"AT");
        assert_eq!(alt_allele, b"A");
    }

    #[test]
    fn replay_event_row_uses_gatk_like_event_key() {
        let call = test_variant("chr1", 100, b"A", b"C");
        let row = replay_event_row("chr1:90-110", &call).unwrap();
        assert_eq!(row.region, "chr1:90-110");
        assert_eq!(row.event, "chr1:100:SNP:A*,C");
        assert_eq!(row.event_type, "SNP");
        assert_eq!(row.alleles, "A*,C");
        assert_eq!(row.gt, "0/1");
    }

    #[test]
    fn dbsnp_record_match_requires_exact_ref_and_alt() {
        let variant = test_variant("chr1", 100, b"A", b"ATG");
        let record = b"chr1\t100\trs1\tA\tC,ATG\t.\t.\t.";
        assert!(dbsnp_record_matches(record, &variant).unwrap());
        assert_eq!(dbsnp_record_id(record).unwrap(), "rs1");
        let non_match = b"chr1\t100\trs2\tA\tC,G\t.\t.\t.";
        assert!(!dbsnp_record_matches(non_match, &variant).unwrap());
    }

    #[test]
    fn variant_calls_sort_by_dictionary_order() {
        let dict = test_dict();
        let mut variants = vec![
            test_variant("chr1", 20, b"A", b"C"),
            test_variant("chr2", 30, b"G", b"T"),
            test_variant("chr2", 10, b"A", b"G"),
        ];
        sort_variant_calls(&mut variants, &dict).unwrap();
        assert_eq!(variants[0].contig, "chr2");
        assert_eq!(variants[0].pos, 10);
        assert_eq!(variants[1].contig, "chr2");
        assert_eq!(variants[1].pos, 30);
        assert_eq!(variants[2].contig, "chr1");
    }
}

fn is_regular_allele(allele: &[u8]) -> bool {
    allele.iter().all(|&b| matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't'))
}

pub fn extract_variants_from_cigar(
    contig: &str,
    ref_bases: &[u8],
    alt_bases: &[u8],
    cigar: &rust_htslib::bam::record::CigarString,
    alignment_offset: i32,
    region_start: u64,
    max_mnp_distance: usize,
) -> Vec<VariantCall> {
    use rust_htslib::bam::record::Cigar::*;
    
    let mut ref_pos = alignment_offset;
    if ref_pos < 0 {
        return Vec::new();
    }
    
    let mut alignment_pos = 0;
    let mut proposed_events = Vec::new();
    let num_cigar_elements = cigar.len();

    for (cigar_index, ce) in cigar.iter().enumerate() {
        match ce {
            Ins(len) => {
                let element_length = *len as usize;
                if element_length <= 10 && ref_pos > 0 && cigar_index > 0 && cigar_index < num_cigar_elements - 1 {
                    let insertion_start = region_start + ref_pos as u64 - 1;
                    let ref_byte = ref_bases[ref_pos as usize - 1];
                    let mut insertion_bases = vec![ref_byte];
                    insertion_bases.extend_from_slice(&alt_bases[alignment_pos..alignment_pos + element_length]);
                    
                    if is_regular_allele(&[ref_byte]) && is_regular_allele(&insertion_bases) {
                        proposed_events.push(VariantCall {
                            contig: contig.to_string(),
                            pos: insertion_start,
                            id: None,
                            db: false,
                            ref_allele: vec![ref_byte],
                            alt_allele: insertion_bases,
                            depth: 0,
                            ref_count: 0,
                            alt_count: 0,
                            qual: 0,
                            fs: 0.0,
                            pl: [0, 0, 0],
                            genotype_index: 0,
                        });
                    }
                }
                alignment_pos += element_length;
            }
            SoftClip(len) => {
                alignment_pos += *len as usize;
            }
            Del(len) => {
                let element_length = *len as usize;
                if element_length <= 10 && ref_pos > 0 {
                    let deletion_start = region_start + ref_pos as u64 - 1;
                    let ref_byte = ref_bases[ref_pos as usize - 1];
                    let mut deletion_bases = vec![ref_byte];
                    deletion_bases.extend_from_slice(&ref_bases[ref_pos as usize..ref_pos as usize + element_length]);
                    
                    if is_regular_allele(&deletion_bases) && is_regular_allele(&[ref_byte]) {
                        proposed_events.push(VariantCall {
                            contig: contig.to_string(),
                            pos: deletion_start,
                            id: None,
                            db: false,
                            ref_allele: deletion_bases,
                            alt_allele: vec![ref_byte],
                            depth: 0,
                            ref_count: 0,
                            alt_count: 0,
                            qual: 0,
                            fs: 0.0,
                            pl: [0, 0, 0],
                            genotype_index: 0,
                        });
                    }
                }
                ref_pos += element_length as i32;
            }
            Match(len) | Equal(len) | Diff(len) => {
                let element_length = *len as usize;
                let mut mismatch_offsets = std::collections::VecDeque::new();
                
                for offset in 0..element_length {
                    let r_idx = ref_pos as usize + offset;
                    let a_idx = alignment_pos + offset;
                    if r_idx < ref_bases.len() && a_idx < alt_bases.len() {
                        let ref_byte = ref_bases[r_idx];
                        let alt_byte = alt_bases[a_idx];
                        // we ignore N vs N mismatches in practice, but keeping simple
                        if ref_byte != alt_byte {
                            mismatch_offsets.push_back(offset);
                        }
                    }
                }

                while let Some(start) = mismatch_offsets.pop_front() {
                    let mut end = start;
                    while let Some(&next) = mismatch_offsets.front() {
                        if next - end <= max_mnp_distance {
                            end = mismatch_offsets.pop_front().unwrap();
                        } else {
                            break;
                        }
                    }
                    
                    let ref_allele = ref_bases[ref_pos as usize + start..=ref_pos as usize + end].to_vec();
                    let alt_allele = alt_bases[alignment_pos + start..=alignment_pos + end].to_vec();
                    
                    if is_regular_allele(&ref_allele) && is_regular_allele(&alt_allele) {
                        proposed_events.push(VariantCall {
                            contig: contig.to_string(),
                            pos: region_start + ref_pos as u64 + start as u64,
                            id: None,
                            db: false,
                            ref_allele,
                            alt_allele,
                            depth: 0,
                            ref_count: 0,
                            alt_count: 0,
                            qual: 0,
                            fs: 0.0,
                            pl: [0, 0, 0],
                            genotype_index: 0,
                        });
                    }
                }

                ref_pos += element_length as i32;
                alignment_pos += element_length;
            }
            _ => {
                // skip others for now
            }
        }
    }
    
    proposed_events
}

pub fn assemble_haplotypes(
    contig: &str,
    region_start: u64,
    ref_bases: &[u8],
    reads_bases: &[Vec<u8>],
    kmer_sizes: &[usize],
    max_mnp_distance: usize,
) -> (Vec<LocalHaplotype>, Vec<VariantCall>) {
    use std::collections::HashSet;
    let mut assembled_haplotypes_set = HashSet::new();
    let mut local_haps = Vec::new();
    
    let ref_hap = LocalHaplotype {
        bases: ref_bases.to_vec(),
        is_ref: true,
        cigar: format!("{}M", ref_bases.len()),
        event_indices: Vec::new(),
    };
    local_haps.push(ref_hap);
    assembled_haplotypes_set.insert(ref_bases.to_vec());
    
    // We'll store events per haplotype as Vec<Vec<VariantCall>>
    let mut hap_events = vec![Vec::new()]; // first is ref, no events
    
    for &kmer_size in kmer_sizes {
        if ref_bases.len() < kmer_size {
            continue;
        }
        let mut graph = crate::assembly::ReadThreadingGraph::new(kmer_size);
        graph.add_sequence(ref_bases, true);
        
        for read_bases in reads_bases {
            if read_bases.len() >= kmer_size {
                graph.add_sequence(read_bases, false);
            }
        }
        
        let source_kmer = &ref_bases[0..kmer_size];
        let sink_kmer = &ref_bases[ref_bases.len() - kmer_size..];
        
        let source_idx = graph.get_or_create_vertex(source_kmer);
        let sink_idx = graph.get_or_create_vertex(sink_kmer);
        
        let best_paths = graph.find_best_haplotypes(source_idx, sink_idx, 10);
        
        let mut found_nonref = false;
        for path in best_paths {
            let seq = graph.reconstruct_sequence(&path);
            if !assembled_haplotypes_set.contains(&seq) {
                assembled_haplotypes_set.insert(seq.clone());
                
                let sw_params = crate::smith_waterman::SWParameters::default();
                let align_result = crate::smith_waterman::align(
                    ref_bases,
                    &seq,
                    &sw_params,
                    crate::smith_waterman::SWOverhangStrategy::SoftClip,
                );
                
                let events = extract_variants_from_cigar(
                    contig,
                    ref_bases,
                    &seq,
                    &align_result.cigar,
                    align_result.alignment_offset,
                    region_start,
                    max_mnp_distance,
                );
                
                hap_events.push(events);
                
                local_haps.push(LocalHaplotype {
                    bases: seq,
                    is_ref: false,
                    cigar: align_result.cigar.to_string(),
                    event_indices: Vec::new(), // To be filled later
                });
                found_nonref = true;
            }
        }
        // Like GATK: if this kmer size produced non-ref haplotypes, stop
        if found_nonref {
            break;
        }
    }
    
    // Deduplicate events across all haplotypes
    let mut unique_events: Vec<VariantCall> = Vec::new();
    for events in &hap_events {
        for event in events {
            if !unique_events.iter().any(|e| e.pos == event.pos && e.ref_allele == event.ref_allele && e.alt_allele == event.alt_allele) {
                unique_events.push(event.clone());
            }
        }
    }
    
    // Now map event indices back to each haplotype
    for (hap_idx, events) in hap_events.iter().enumerate() {
        for event in events {
            if let Some(idx) = unique_events.iter().position(|e| e.pos == event.pos && e.ref_allele == event.ref_allele && e.alt_allele == event.alt_allele) {
                if !local_haps[hap_idx].event_indices.contains(&idx) {
                    local_haps[hap_idx].event_indices.push(idx);
                }
            }
        }
    }
    
    (local_haps, unique_events)
}
