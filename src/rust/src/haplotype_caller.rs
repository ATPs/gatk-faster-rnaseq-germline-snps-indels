use crate::pair_hmm;
use anyhow::{bail, Context, Result};
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use rust_htslib::bam::ext::BamRecordExtensions;
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
// GATK PairHMM squashes base qualities below this threshold down to the minimum
// usable quality instead of raising them to the threshold itself.
const PAIR_HMM_BASE_QUALITY_SCORE_THRESHOLD: u8 = 18;
const PAIR_HMM_MIN_USABLE_Q_SCORE: u8 = 6;
const FISHER_STRAND_TARGET_TABLE_SIZE: f64 = 200.0;
const FISHER_STRAND_MIN_PVALUE: f64 = 1e-320;
const CALL_PARTITIONS_PER_THREAD: usize = 8;
const ACTIVE_REGION_MAX_GAP: u64 = 50;
const ACTIVE_REGION_PADDING: u64 = 100;
const SNP_CLUSTER_WINDOW: u64 = 35;
const WEAK_SUPPLEMENTAL_CLUSTER_SNP_MAX_ALT_COUNT: u32 = 4;
const STRONG_SINGLE_SNP_RESCUE_MIN_ALT_COUNT: u32 = 8;
const STRONG_SINGLE_SNP_RESCUE_MIN_ALT_FRACTION: f64 = 0.35;

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
    pub exclude_supplementary: bool,
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
        let alt_prob = if !is_ref { 1.0 - error } else { error / 3.0 }.max(f64::MIN_POSITIVE);
        add_diploid_observation(&mut log10_likelihoods, ref_prob, alt_prob);
    }
    let denominator = log10_sum_exp(&log10_likelihoods);
    let ref_posterior_log10 = log10_likelihoods[0] - denominator;
    let qual = phred_from_log10_probability(ref_posterior_log10);
    (f64::from(qual) >= ACTIVE_REGION_DISCOVERY_CONFIDENCE, qual)
}

fn is_active_locus(ref_index: Option<usize>, evidence: &SnpEvidence, _depth: u32) -> (bool, u32) {
    let ref_index = match ref_index {
        Some(idx) => idx,
        None => return (false, 0),
    };
    if evidence.active_observations.len() < 2 {
        return (false, 0);
    }

    // Compute ref-vs-any log10 likelihoods for diploid genotype states.
    let mut log10_likelihoods = [0.0_f64; 3]; // hom-ref, het, hom-alt
    for observation in &evidence.active_observations {
        let error = phred_error_probability(observation.quality);
        let is_ref = observation.base_index == Some(ref_index);
        let ref_prob = snp_observation_probability(is_ref, error);
        // For ref-vs-any, Java treats any non-reference read base, including N,
        // as alt evidence during active-region discovery.
        let alt_prob = if !is_ref { 1.0 - error } else { error / 3.0 }.max(f64::MIN_POSITIVE);
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

struct PreparedHmmRead {
    bases: Vec<u8>,
    quals: Vec<u8>,
    ins_quals: Vec<u8>,
    del_quals: Vec<u8>,
    assembly_segments: Vec<Vec<u8>>,
    ref_span: (u64, u64),
}

fn append_prepared_read_slice(
    bases: &mut Vec<u8>,
    quals_out: &mut Vec<u8>,
    ins_quals_out: &mut Vec<u8>,
    del_quals_out: &mut Vec<u8>,
    assembly_quals_out: &mut Vec<u8>,
    seq: &rust_htslib::bam::record::Seq<'_>,
    quals: &[u8],
    bi_bytes: Option<&[u8]>,
    bd_bytes: Option<&[u8]>,
    mapq: u8,
    start: usize,
    end: usize,
) {
    for i in start..end {
        if i >= seq.len() {
            break;
        }
        bases.push(seq[i]);
        quals_out.push(pair_hmm_base_quality(quals[i], mapq));
        assembly_quals_out.push(quals[i]);
        ins_quals_out.push(
            bi_bytes
                .map(|bytes| pair_hmm_indel_open_quality(bytes[i].saturating_sub(33)))
                .unwrap_or(45),
        );
        del_quals_out.push(
            bd_bytes
                .map(|bytes| pair_hmm_indel_open_quality(bytes[i].saturating_sub(33)))
                .unwrap_or(45),
        );
    }
}

fn prepared_read_reference_start(
    record: &bam::Record,
    dont_use_soft_clipped_bases: bool,
) -> Option<u64> {
    let aligned_start = u64::try_from(record.pos()).ok()?.saturating_add(1);
    if dont_use_soft_clipped_bases {
        return Some(aligned_start);
    }

    let leading_soft_clip = record
        .cigar()
        .iter()
        .find_map(|op| match op {
            Cigar::SoftClip(len) => Some(u64::from(*len)),
            Cigar::HardClip(_) => None,
            _ => Some(0),
        })
        .unwrap_or(0);
    Some(aligned_start.saturating_sub(leading_soft_clip))
}

#[allow(clippy::too_many_arguments)]
fn append_prepared_ref_consuming_segment(
    bases: &mut Vec<u8>,
    quals_out: &mut Vec<u8>,
    ins_quals_out: &mut Vec<u8>,
    del_quals_out: &mut Vec<u8>,
    assembly_quals_out: &mut Vec<u8>,
    ref_span_start: &mut Option<u64>,
    ref_span_end: &mut Option<u64>,
    seq: &rust_htslib::bam::record::Seq<'_>,
    quals: &[u8],
    bi_bytes: Option<&[u8]>,
    bd_bytes: Option<&[u8]>,
    mapq: u8,
    clip_start: usize,
    clip_end: usize,
    qpos: usize,
    len: usize,
    ref_pos: u64,
    region_start: u64,
    region_end_exclusive: u64,
) {
    let op_start = qpos;
    let op_end = qpos + len;
    let ref_end = ref_pos + len as u64;

    let region_clip_start = if ref_pos < region_start {
        op_start + (region_start - ref_pos) as usize
    } else {
        op_start
    };
    let region_clip_end = if ref_end > region_end_exclusive {
        op_start + (region_end_exclusive.saturating_sub(ref_pos)) as usize
    } else {
        op_end
    };

    let kept_start = clip_start.max(region_clip_start);
    let kept_end = clip_end.min(region_clip_end);
    if kept_start >= kept_end {
        return;
    }

    append_prepared_read_slice(
        bases,
        quals_out,
        ins_quals_out,
        del_quals_out,
        assembly_quals_out,
        seq,
        quals,
        bi_bytes,
        bd_bytes,
        mapq,
        kept_start,
        kept_end,
    );
    let segment_start = ref_pos + (kept_start - op_start) as u64;
    let segment_end = ref_pos + (kept_end - op_start) as u64 - 1;
    ref_span_start.get_or_insert(segment_start);
    *ref_span_end = Some(segment_end);
}

fn base_is_usable_for_assembly(base: u8, qual: u8) -> bool {
    normalize_base(base) != b'N' && qual >= PIPELINE_MIN_BASEQ
}

fn assembly_read_segments(bases: &[u8], quals: &[u8]) -> Vec<Vec<u8>> {
    let mut segments = Vec::new();
    let mut start = None;

    for end in 0..=bases.len() {
        let usable = end < bases.len() && base_is_usable_for_assembly(bases[end], quals[end]);
        if usable {
            if start.is_none() {
                start = Some(end);
            }
        } else if let Some(segment_start) = start.take() {
            if segment_start < end {
                segments.push(bases[segment_start..end].to_vec());
            }
        }
    }

    segments
}

fn prepare_hmm_read(
    record: &bam::Record,
    min_tail_quality: u8,
    dont_use_soft_clipped_bases: bool,
    region_start: u64,
    region_end: u64,
) -> Option<PreparedHmmRead> {
    let (clip_start, clip_end) =
        clip_read_for_evidence(record, min_tail_quality, dont_use_soft_clipped_bases)?;
    let mapq = record.mapq();
    let seq = record.seq();
    let quals = record.qual();

    let bi_bytes: Option<&[u8]> = match record.aux(b"BI") {
        Ok(rust_htslib::bam::record::Aux::String(s)) if s.len() == seq.len() => Some(s.as_bytes()),
        _ => None,
    };
    let bd_bytes: Option<&[u8]> = match record.aux(b"BD") {
        Ok(rust_htslib::bam::record::Aux::String(s)) if s.len() == seq.len() => Some(s.as_bytes()),
        _ => None,
    };

    let mut bases = Vec::new();
    let mut read_quals = Vec::new();
    let mut ins_quals = Vec::new();
    let mut del_quals = Vec::new();
    let mut assembly_quals = Vec::new();
    let mut qpos = 0_usize;
    let mut ref_pos = prepared_read_reference_start(record, dont_use_soft_clipped_bases)?;
    let mut ref_span_start = None;
    let mut ref_span_end = None;
    let region_end_exclusive = region_end.saturating_add(1);

    for view in record.cigar().iter() {
        use rust_htslib::bam::record::Cigar::*;
        let len = view.len() as usize;
        match view {
            Match(_) | Equal(_) | Diff(_) => {
                append_prepared_ref_consuming_segment(
                    &mut bases,
                    &mut read_quals,
                    &mut ins_quals,
                    &mut del_quals,
                    &mut assembly_quals,
                    &mut ref_span_start,
                    &mut ref_span_end,
                    &seq,
                    quals,
                    bi_bytes,
                    bd_bytes,
                    mapq,
                    clip_start,
                    clip_end,
                    qpos,
                    len,
                    ref_pos,
                    region_start,
                    region_end_exclusive,
                );
                qpos += len;
                ref_pos += len as u64;
            }
            Ins(_) => {
                let kept_start = clip_start.max(qpos);
                let kept_end = clip_end.min(qpos + len);
                if kept_start < kept_end
                    && ref_pos >= region_start
                    && ref_pos <= region_end_exclusive
                {
                    append_prepared_read_slice(
                        &mut bases,
                        &mut read_quals,
                        &mut ins_quals,
                        &mut del_quals,
                        &mut assembly_quals,
                        &seq,
                        quals,
                        bi_bytes,
                        bd_bytes,
                        mapq,
                        kept_start,
                        kept_end,
                    );
                }
                qpos += len;
            }
            SoftClip(_) => {
                if !dont_use_soft_clipped_bases {
                    append_prepared_ref_consuming_segment(
                        &mut bases,
                        &mut read_quals,
                        &mut ins_quals,
                        &mut del_quals,
                        &mut assembly_quals,
                        &mut ref_span_start,
                        &mut ref_span_end,
                        &seq,
                        quals,
                        bi_bytes,
                        bd_bytes,
                        mapq,
                        clip_start,
                        clip_end,
                        qpos,
                        len,
                        ref_pos,
                        region_start,
                        region_end_exclusive,
                    );
                    ref_pos += len as u64;
                }
                qpos += len;
            }
            Del(_) | RefSkip(_) => {
                ref_pos += len as u64;
            }
            HardClip(_) | Pad(_) => {}
        }
    }

    if bases.len() < MIN_READ_LENGTH_AFTER_TRIMMING {
        return None;
    }
    let (ref_span_start, ref_span_end) = match (ref_span_start, ref_span_end) {
        (Some(start), Some(end)) => (start, end),
        _ => return None,
    };

    Some(PreparedHmmRead {
        assembly_segments: assembly_read_segments(&bases, &assembly_quals),
        bases,
        quals: read_quals,
        ins_quals,
        del_quals,
        ref_span: (ref_span_start, ref_span_end),
    })
}

pub fn call_variants(config: &HaplotypeCallerConfig) -> Result<()> {
    validate_haplotype_caller_config(config)?;
    let (dict, mut intervals) = read_interval_list(&config.input_interval_list)?;
    sort_intervals(&mut intervals, &dict)?;
    let fetch_windows = coalesce_fetch_windows(&intervals, &dict);
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
    let fetch_windows = coalesce_fetch_windows(&intervals, &dict);
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

#[derive(Copy, Clone, Debug, Eq, PartialEq)]
struct ActiveBaseObservation {
    base_index: Option<usize>,
    quality: u8,
}

#[derive(Clone, Debug, Default, Eq, PartialEq)]
struct SnpEvidence {
    counts: BaseCounts,
    strands: [StrandCounts; 4],
    observations: Vec<BaseObservation>,
    active_observations: Vec<ActiveBaseObservation>,
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
    active_start: u64,
    active_end: u64,
    padded_start: u64,
    padded_end: u64,
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

#[derive(Clone, Debug, Eq, PartialEq)]
struct ActiveRegionSpan {
    active: Interval,
    padded: Interval,
}

#[derive(Clone, Debug)]
struct NamedSnpObservation {
    read_name: Vec<u8>,
    qpos: usize,
    mapq: u8,
    base_index: Option<usize>,
    quality: u8,
    is_reverse: bool,
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

fn select_assembly_reads<'a>(
    read_assembly_segments_list: &'a [Vec<Vec<u8>>],
    read_use_for_assembly_list: &[bool],
) -> Vec<&'a Vec<u8>> {
    read_assembly_segments_list
        .iter()
        .zip(read_use_for_assembly_list.iter().copied())
        .flat_map(|(segments, use_for_assembly)| {
            use_for_assembly
                .then_some(segments.as_slice())
                .unwrap_or(&[])
                .iter()
        })
        .collect()
}

fn discover_active_regions(
    config: &HaplotypeCallerConfig,
    tid: u32,
    window: &FetchWindow,
    ref_bases: &[u8],
    bam: &mut bam::IndexedReader,
) -> Result<Vec<ActiveRegionSpan>> {
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
            config.exclude_supplementary,
            config.dont_use_soft_clipped_bases,
        );
        let indel_evidence = pileup_indel_evidence(
            &pileup,
            PIPELINE_MIN_BASEQ,
            PIPELINE_MIN_MAPQ,
            PIPELINE_MIN_TAIL_QUALITY,
            config.exclude_supplementary,
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

    let mut coalesced = Vec::new();
    for reg in split_active_regions {
        let padded = pad_active_region(&reg, window);
        coalesced.push(ActiveRegionSpan {
            active: reg,
            padded,
        });
    }

    Ok(coalesced)
}

fn pad_active_region(interval: &Interval, window: &FetchWindow) -> Interval {
    Interval {
        contig: interval.contig.clone(),
        start: interval
            .start
            .saturating_sub(ACTIVE_REGION_PADDING)
            .max(window.start),
        end: interval
            .end
            .saturating_add(ACTIVE_REGION_PADDING)
            .min(window.end),
    }
}

fn collect_call_active_loci_rows(
    config: &HaplotypeCallerConfig,
    tid: u32,
    window: &FetchWindow,
    active_span: &Interval,
    region: &Interval,
    ref_bases: &[u8],
    bam: &mut bam::IndexedReader,
) -> Result<Vec<ReplayActiveLocusRow>> {
    let mut active_loci_rows = Vec::new();
    let region_label = region_name(&region.contig, region.start, region.end);
    bam.fetch((
        tid as i32,
        (active_span.start - 1) as i64,
        active_span.end as i64,
    ))
    .with_context(|| {
        format!(
            "fetching BAM region for pileup fallback {}:{}-{}",
            active_span.contig, active_span.start, active_span.end
        )
    })?;
    for pileup_result in bam.pileup() {
        let pileup = pileup_result?;
        let pos0 = pileup.pos() as u64;
        let pos1 = pos0 + 1;
        if pos1 < active_span.start || pos1 > active_span.end {
            continue;
        }
        let ref_base = normalize_base(ref_bases[(pos0 - (window.start - 1)) as usize]);
        let snp_evidence = pileup_snp_evidence(
            &pileup,
            PIPELINE_MIN_BASEQ,
            PIPELINE_MIN_MAPQ,
            PIPELINE_MIN_TAIL_QUALITY,
            config.exclude_supplementary,
            config.dont_use_soft_clipped_bases,
        );
        let indel_evidence = pileup_indel_evidence(
            &pileup,
            PIPELINE_MIN_BASEQ,
            PIPELINE_MIN_MAPQ,
            PIPELINE_MIN_TAIL_QUALITY,
            config.exclude_supplementary,
            config.dont_use_soft_clipped_bases,
        );
        let ref_index = base_index(ref_base);
        let snp_alt = best_snp_alt(ref_index, &snp_evidence);
        let indel_alt = best_indel_alt(&indel_evidence);
        let snp_alt_count = snp_alt.map(|(_, count)| count).unwrap_or(0);
        let indel_alt_count = indel_alt.map(|(_, count)| *count).unwrap_or(0);
        let depth = snp_evidence.counts.depth.max(indel_evidence.counts.depth);
        let is_active = is_active_locus(ref_index, &snp_evidence, depth).0
            || is_active_indel(&indel_evidence).0;
        if is_active {
            let best_alt_count = snp_alt_count.max(indel_alt_count);
            let alt_fraction = if depth == 0 {
                0.0
            } else {
                f64::from(best_alt_count) / f64::from(depth)
            };
            active_loci_rows.push(ReplayActiveLocusRow {
                contig: region.contig.clone(),
                pos: pos1,
                region: region_label.clone(),
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
                active_probability_proxy: 1.0,
            });
        }
    }

    Ok(active_loci_rows)
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

    for region in active_regions {
        let active_span = region.active.clone();
        let region = region.padded;
        let active_loci_rows = collect_call_active_loci_rows(
            config,
            tid,
            window,
            &active_span,
            &region,
            &ref_bases,
            bam,
        )?;
        bam.fetch((tid as i32, (region.start - 1) as i64, region.end as i64))
            .with_context(|| {
                format!(
                    "fetching BAM region {}:{}-{}",
                    region.contig, region.start, region.end
                )
            })?;

        let mut read_bases_list = Vec::new();
        let mut read_quals_list = Vec::new();
        let mut read_ins_quals_list = Vec::new();
        let mut read_del_quals_list = Vec::new();
        let mut read_is_reverse_list = Vec::new();
        let mut read_ref_spans = Vec::new();
        let mut read_use_for_assembly_list = Vec::new();
        let mut read_assembly_segments_list = Vec::new();
        let mut reads_by_start = std::collections::HashMap::new();

        for r in bam.records() {
            let record = r?;
            if !read_passes_hc_filter(&record, PIPELINE_MIN_MAPQ, config.exclude_supplementary) {
                continue;
            }
            let Some(prepared_read) = prepare_hmm_read(
                &record,
                PIPELINE_MIN_TAIL_QUALITY,
                config.dont_use_soft_clipped_bases,
                region.start,
                region.end,
            ) else {
                continue;
            };

            let start_pos = record.pos();
            let count = reads_by_start.entry(start_pos).or_insert(0);
            let use_for_assembly = *count < 50;
            if use_for_assembly {
                *count += 1;
            }

            read_assembly_segments_list.push(prepared_read.assembly_segments);
            read_bases_list.push(prepared_read.bases);
            read_quals_list.push(prepared_read.quals);
            read_ins_quals_list.push(prepared_read.ins_quals);
            read_del_quals_list.push(prepared_read.del_quals);
            read_is_reverse_list.push(record.is_reverse());
            read_ref_spans.push(prepared_read.ref_span);
            read_use_for_assembly_list.push(use_for_assembly);
        }

        if read_bases_list.is_empty() {
            continue;
        }

        let ref_region_start_offset = (region.start - window.start) as usize;
        let ref_region_end_offset = (region.end - window.start) as usize;
        let local_ref_bases = &ref_bases[ref_region_start_offset..=ref_region_end_offset];

        // 1. Assemble haplotypes instead of pileup!
        let max_mnp_distance = 0; // Default GATK
        let assembly_reads =
            select_assembly_reads(&read_assembly_segments_list, &read_use_for_assembly_list);
        let assembly_reads_owned: Vec<Vec<u8>> =
            assembly_reads.iter().map(|r| (*r).clone()).collect();
        let (mut local_haplotypes, mut valid_events) = assemble_haplotypes(
            &region.contig,
            region.start,
            local_ref_bases,
            &assembly_reads_owned,
            &[10, 25],
            max_mnp_distance,
        );
        supplement_missing_pileup_events(
            &region.contig,
            region.start,
            local_ref_bases,
            &active_loci_rows,
            config.standard_min_confidence_threshold_for_calling,
            &mut local_haplotypes,
            &mut valid_events,
        );
        filter_non_acgt_haplotypes_for_single_snp_region(&mut local_haplotypes, &valid_events);

        if valid_events.is_empty() {
            continue;
        }

        let n_reads = read_bases_list.len();
        let mut read_haplotype_likelihoods: Vec<Vec<f64>> = Vec::with_capacity(n_reads);

        for i in 0..n_reads {
            let r_bases = &read_bases_list[i];
            let r_quals = &read_quals_list[i];
            let read_ins_quals = &read_ins_quals_list[i];
            let read_del_quals = &read_del_quals_list[i];
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

        let mut final_calls = genotype_assembled_events(
            &local_haplotypes,
            &valid_events,
            &read_haplotype_likelihoods,
            &read_is_reverse_list,
            &read_ref_spans,
            config.standard_min_confidence_threshold_for_calling,
        );
        prune_unsupported_simple_snp_calls_in_dense_clusters(&mut final_calls, &active_loci_rows);
        merge_missing_strong_snp_cluster_rescues_from_pileup(
            &mut final_calls,
            &region.contig,
            region.start,
            local_ref_bases,
            &active_loci_rows,
            &valid_events,
            config.standard_min_confidence_threshold_for_calling,
        );
        if final_calls.is_empty() {
            final_calls = rescue_collapsed_strong_snp_cluster_from_pileup(
                &region.contig,
                region.start,
                local_ref_bases,
                &active_loci_rows,
                &valid_events,
                config.standard_min_confidence_threshold_for_calling,
            );
        }
        output.variants.extend(final_calls);
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
    for discovered_region in active_regions {
        let active_span = &discovered_region.active;
        let region_interval = &discovered_region.padded;
        bam.fetch((
            tid as i32,
            (region_interval.start - 1) as i64,
            region_interval.end as i64,
        ))
        .with_context(|| {
            format!(
                "fetching BAM region {}:{}-{}",
                region_interval.contig, region_interval.start, region_interval.end
            )
        })?;

        let region = region_name(
            &region_interval.contig,
            region_interval.start,
            region_interval.end,
        );
        let mut active_region = ReplayActiveRegionRow {
            contig: region_interval.contig.clone(),
            start: region_interval.start,
            end: region_interval.end,
            region: region.clone(),
            active_start: active_span.start,
            active_end: active_span.end,
            padded_start: region_interval.start,
            padded_end: region_interval.end,
            ..ReplayActiveRegionRow::default()
        };

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
                config.exclude_supplementary,
                config.dont_use_soft_clipped_bases,
                Some(&row_context),
            );
            let (indel_evidence, indel_rows) = pileup_indel_evidence_with_rows(
                &pileup,
                PIPELINE_MIN_BASEQ,
                PIPELINE_MIN_MAPQ,
                PIPELINE_MIN_TAIL_QUALITY,
                config.exclude_supplementary,
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
            let is_active = is_active_locus(ref_index, &snp_evidence, depth).0
                || is_active_indel(&indel_evidence).0;
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

            // We removed best_snp_call and best_indel_call here
            // since we will use assemble_haplotypes below.
        }

        let call_active_loci_rows = collect_call_active_loci_rows(
            config,
            tid,
            window,
            active_span,
            region_interval,
            &ref_bases,
            bam,
        )?;

        let ref_region_start_offset = (region_interval.start - window.start) as usize;
        let ref_region_end_offset = (region_interval.end - window.start) as usize;
        let local_ref_bases = &ref_bases[ref_region_start_offset..=ref_region_end_offset];

        bam.fetch((
            tid as i32,
            (region_interval.start - 1) as i64,
            region_interval.end as i64,
        ))
        .with_context(|| {
            format!(
                "fetching BAM region for PairHMM: {}:{}-{}",
                region_interval.contig, region_interval.start, region_interval.end
            )
        })?;

        let mut read_bases_list = Vec::new();
        let mut read_quals_list = Vec::new();
        let mut read_ins_quals_list = Vec::new();
        let mut read_del_quals_list = Vec::new();
        let mut read_use_for_assembly_list = Vec::new();
        let mut read_assembly_segments_list = Vec::new();
        let mut reads_by_start = std::collections::HashMap::new();
        let mut read_names_list = Vec::new();
        let mut read_is_reverse_list = Vec::new();
        let mut read_ref_spans = Vec::new();
        let mut mapq_list = Vec::new();
        let mut unclipped_loc_list = Vec::new();
        let mut cigar_string_list = Vec::new();

        for r in bam.records() {
            let record = r?;
            if !read_passes_hc_filter(&record, PIPELINE_MIN_MAPQ, config.exclude_supplementary) {
                continue;
            }
            let Some(prepared_read) = prepare_hmm_read(
                &record,
                PIPELINE_MIN_TAIL_QUALITY,
                config.dont_use_soft_clipped_bases,
                region_interval.start,
                region_interval.end,
            ) else {
                continue;
            };

            let start_pos = record.pos();
            let count = reads_by_start.entry(start_pos).or_insert(0);
            let use_for_assembly = *count < 50;
            if use_for_assembly {
                *count += 1;
            }

            read_assembly_segments_list.push(prepared_read.assembly_segments);
            read_bases_list.push(prepared_read.bases);
            read_quals_list.push(prepared_read.quals);
            read_ins_quals_list.push(prepared_read.ins_quals);
            read_del_quals_list.push(prepared_read.del_quals);
            read_use_for_assembly_list.push(use_for_assembly);
            read_names_list.push(String::from_utf8_lossy(record.qname()).into_owned());
            read_is_reverse_list.push(record.is_reverse());
            mapq_list.push(record.mapq());
            unclipped_loc_list.push((record.pos() + 1) as u64);
            cigar_string_list.push(record.cigar().to_string());
            read_ref_spans.push(prepared_read.ref_span);
        }

        let assembly_reads =
            select_assembly_reads(&read_assembly_segments_list, &read_use_for_assembly_list);
        let assembly_reads_owned: Vec<Vec<u8>> =
            assembly_reads.iter().map(|r| (*r).clone()).collect();
        let max_mnp_distance = 0; // Default GATK
        let (mut local_haplotypes, mut valid_events) = assemble_haplotypes(
            &region_interval.contig,
            region_interval.start,
            local_ref_bases,
            &assembly_reads_owned,
            &[10, 25],
            max_mnp_distance,
        );
        supplement_missing_pileup_events(
            &region_interval.contig,
            region_interval.start,
            local_ref_bases,
            &call_active_loci_rows,
            config.standard_min_confidence_threshold_for_calling,
            &mut local_haplotypes,
            &mut valid_events,
        );
        filter_non_acgt_haplotypes_for_single_snp_region(&mut local_haplotypes, &valid_events);

        for call in &valid_events {
            active_region.candidate_events += 1;
            output.events.push(replay_event_row(&region, call)?);
        }

        for (hap_idx, hap) in local_haplotypes.iter().enumerate() {
            output.haplotypes.push(ReplayHaplotypeRow {
                region: region.clone(),
                stage: "assembled",
                haplotype: hap_idx,
                span_start: region_interval.start,
                span_end: region_interval.end,
                kmer: 0,
                length: hap.bases.len() as u32,
                cigar: hap.cigar.clone(),
                is_ref: hap.is_ref,
                bases: String::from_utf8_lossy(&hap.bases).into_owned(),
            });
        }

        let n_reads = read_bases_list.len();
        let mut read_haplotype_likelihoods: Vec<Vec<f64>> = Vec::with_capacity(n_reads);
        for i in 0..n_reads {
            let read_bases = &read_bases_list[i];
            let read_quals = &read_quals_list[i];
            let read_ins_quals = &read_ins_quals_list[i];
            let read_del_quals = &read_del_quals_list[i];
            let gcp = 10;
            let mut hap_likelihoods = Vec::with_capacity(local_haplotypes.len());

            for (hap_idx, hap) in local_haplotypes.iter().enumerate() {
                let score = pair_hmm::compute_read_likelihood_given_haplotype(
                    &hap.bases,
                    read_bases,
                    read_quals,
                    read_ins_quals,
                    read_del_quals,
                    gcp,
                );

                output.pairhmms.push(ReplayPairHmmRow {
                    region: region.clone(),
                    read: read_names_list[i].clone(),
                    haplotype: hap_idx,
                    read_index: i,
                    cigar: cigar_string_list[i].clone(),
                    mapq: mapq_list[i],
                    loc: unclipped_loc_list[i],
                    unclipped_loc: unclipped_loc_list[i],
                    length: read_bases.len() as u32,
                    score,
                });
                hap_likelihoods.push(score);
            }
            read_haplotype_likelihoods.push(hap_likelihoods);
        }

        let mut final_calls = genotype_assembled_events(
            &local_haplotypes,
            &valid_events,
            &read_haplotype_likelihoods,
            &read_is_reverse_list,
            &read_ref_spans,
            config.standard_min_confidence_threshold_for_calling,
        );
        prune_unsupported_simple_snp_calls_in_dense_clusters(
            &mut final_calls,
            &call_active_loci_rows,
        );
        merge_missing_strong_snp_cluster_rescues_from_pileup(
            &mut final_calls,
            &region_interval.contig,
            region_interval.start,
            local_ref_bases,
            &call_active_loci_rows,
            &valid_events,
            config.standard_min_confidence_threshold_for_calling,
        );
        if final_calls.is_empty() {
            final_calls = rescue_collapsed_strong_snp_cluster_from_pileup(
                &region_interval.contig,
                region_interval.start,
                local_ref_bases,
                &call_active_loci_rows,
                &valid_events,
                config.standard_min_confidence_threshold_for_calling,
            );
        }
        output.variants.extend(final_calls);

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
    exclude_supplementary: bool,
    dont_use_soft_clipped_bases: bool,
) -> SnpEvidence {
    pileup_snp_evidence_with_rows(
        pileup,
        min_baseq,
        min_mapq,
        min_tail_quality,
        exclude_supplementary,
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
    exclude_supplementary: bool,
    dont_use_soft_clipped_bases: bool,
    row_context: Option<&ReplayRowContext<'_>>,
) -> (SnpEvidence, Vec<ReplayReadObservationRow>) {
    let mut observations_by_fragment: HashMap<Vec<u8>, Vec<NamedSnpObservation>> = HashMap::new();
    for alignment in pileup.alignments() {
        let record = alignment.record();
        if !read_passes_hc_filter(&record, min_mapq, exclude_supplementary)
            || alignment.is_refskip()
        {
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
        observations_by_fragment
            .entry(record.qname().to_vec())
            .or_default()
            .push(NamedSnpObservation {
                read_name: record.qname().to_vec(),
                qpos,
                mapq: record.mapq(),
                base_index: base_index(base),
                quality: record.qual()[qpos],
                is_reverse: record.is_reverse(),
            });
    }

    let mut evidence = SnpEvidence::default();
    let mut rows = Vec::new();
    for observations in observations_by_fragment.into_values() {
        for named in adjust_named_snp_observations(&observations) {
            if named.quality < min_baseq {
                continue;
            }
            evidence.active_observations.push(ActiveBaseObservation {
                base_index: named.base_index,
                quality: named.quality,
            });
            if let Some(index) = named.base_index {
                let observation = BaseObservation {
                    base_index: index,
                    quality: named.quality,
                    is_reverse: named.is_reverse,
                };
                evidence.counts.counts[index] += 1;
                evidence.counts.depth += 1;
                evidence.strands[index].increment(named.is_reverse);
                evidence.observations.push(observation);
                if let Some(context) = row_context {
                    rows.push(ReplayReadObservationRow {
                        region: context.region.to_string(),
                        read: String::from_utf8_lossy(&named.read_name).into_owned(),
                        kind: "snp",
                        pos: context.pos,
                        qpos: named.qpos,
                        allele: (base_from_index(index) as char).to_string(),
                        adjusted_quality: named.quality,
                        mapq: named.mapq,
                        strand: strand_label(named.is_reverse),
                    });
                }
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

fn adjust_named_snp_observations(observations: &[NamedSnpObservation]) -> Vec<NamedSnpObservation> {
    if observations.len() <= 1 {
        return observations.to_vec();
    }

    let first_base_index = observations[0].base_index;
    let all_same_base = observations
        .iter()
        .all(|observation| observation.base_index == first_base_index);
    observations
        .iter()
        .map(|observation| {
            let mut adjusted = observation.clone();
            adjusted.quality = if all_same_base {
                adjusted.quality.min(HALF_DEFAULT_PCR_SNV_QUAL)
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
    exclude_supplementary: bool,
    dont_use_soft_clipped_bases: bool,
) -> IndelEvidence {
    pileup_indel_evidence_with_rows(
        pileup,
        min_baseq,
        min_mapq,
        min_tail_quality,
        exclude_supplementary,
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
    exclude_supplementary: bool,
    dont_use_soft_clipped_bases: bool,
    row_context: Option<&ReplayRowContext<'_>>,
) -> (IndelEvidence, Vec<ReplayReadObservationRow>) {
    let mut observations_by_fragment: HashMap<Vec<u8>, Vec<NamedIndelObservation>> = HashMap::new();
    for alignment in pileup.alignments() {
        let record = alignment.record();
        if !read_passes_hc_filter(&record, min_mapq, exclude_supplementary)
            || alignment.is_refskip()
        {
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

fn pair_hmm_base_quality(base_quality: u8, mapq: u8) -> u8 {
    let capped = base_quality.min(mapq);
    if capped < PAIR_HMM_BASE_QUALITY_SCORE_THRESHOLD {
        PAIR_HMM_MIN_USABLE_Q_SCORE
    } else {
        capped
    }
}

fn pair_hmm_indel_open_quality(quality: u8) -> u8 {
    quality.max(PAIR_HMM_MIN_USABLE_Q_SCORE)
}

fn filter_non_acgt_haplotypes_for_single_snp_region(
    local_haplotypes: &mut Vec<LocalHaplotype>,
    valid_events: &[VariantCall],
) {
    if valid_events.len() != 1 {
        return;
    }
    let event = &valid_events[0];
    if event.ref_allele.len() != 1 || event.alt_allele.len() != 1 {
        return;
    }

    let has_regular_alt = local_haplotypes
        .iter()
        .any(|hap| !hap.is_ref && hap.event_indices.contains(&0) && is_regular_bases(&hap.bases));
    if !has_regular_alt {
        return;
    }

    local_haplotypes.retain(|hap| hap.is_ref || is_regular_bases(&hap.bases));
}

fn genotype_assembled_events(
    local_haplotypes: &[LocalHaplotype],
    valid_events: &[VariantCall],
    read_haplotype_likelihoods: &[Vec<f64>],
    read_is_reverse_list: &[bool],
    read_ref_spans: &[(u64, u64)],
    min_confidence: f64,
) -> Vec<VariantCall> {
    if read_haplotype_likelihoods.is_empty() {
        return Vec::new();
    }

    let use_pair_genotyping = overlapping_event_mask(valid_events);
    let overlapping_event_indices = overlapping_event_indices(valid_events);
    let pair_context = use_pair_genotyping
        .iter()
        .copied()
        .any(|use_pair| use_pair)
        .then(|| build_pair_genotyping_context(local_haplotypes.len(), read_haplotype_likelihoods));
    let mut final_calls = Vec::new();

    for (event_idx, event) in valid_events.iter().enumerate() {
        let maybe_call = if use_pair_genotyping[event_idx] {
            genotype_overlapping_assembled_event(
                local_haplotypes,
                event_idx,
                event,
                read_haplotype_likelihoods,
                read_is_reverse_list,
                read_ref_spans,
                min_confidence,
                &overlapping_event_indices[event_idx],
                pair_context
                    .as_ref()
                    .expect("pair genotyping context must exist for overlapping events"),
            )
        } else {
            genotype_isolated_assembled_event(
                local_haplotypes,
                event_idx,
                event,
                read_haplotype_likelihoods,
                read_is_reverse_list,
                read_ref_spans,
                min_confidence,
            )
        };
        if let Some(final_call) = maybe_call {
            final_calls.push(final_call);
        }
    }

    final_calls
}

struct PairGenotypingContext {
    haplotype_pairs: Vec<(usize, usize)>,
    pair_log10_likelihoods: Vec<f64>,
    best_pair: (usize, usize),
}

fn build_pair_genotyping_context(
    n_haplotypes: usize,
    read_haplotype_likelihoods: &[Vec<f64>],
) -> PairGenotypingContext {
    let haplotype_pairs = enumerate_haplotype_pairs(n_haplotypes);
    let pair_log10_likelihoods =
        compute_pair_log10_likelihoods(&haplotype_pairs, read_haplotype_likelihoods);
    let best_pair_idx = pair_log10_likelihoods
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
        .map(|(idx, _)| idx)
        .unwrap_or(0);
    let best_pair = haplotype_pairs[best_pair_idx];

    PairGenotypingContext {
        haplotype_pairs,
        pair_log10_likelihoods,
        best_pair,
    }
}

fn genotype_isolated_assembled_event(
    local_haplotypes: &[LocalHaplotype],
    event_idx: usize,
    event: &VariantCall,
    read_haplotype_likelihoods: &[Vec<f64>],
    read_is_reverse_list: &[bool],
    read_ref_spans: &[(u64, u64)],
    min_confidence: f64,
) -> Option<VariantCall> {
    let hap_contains_allele = event_haplotype_mask(local_haplotypes, event_idx);
    let (ref_haps, alt_haps) = split_event_haplotype_indices(&hap_contains_allele);
    if alt_haps.is_empty() {
        return None;
    }

    let (depth, ref_count, alt_count, fs) = count_event_support(
        read_haplotype_likelihoods,
        read_is_reverse_list,
        read_ref_spans,
        event,
        &ref_haps,
        &alt_haps,
    );

    let mut log10_likelihoods = [0.0; 3];
    for read_idx in 0..read_haplotype_likelihoods.len() {
        if !read_span_overlaps_event(read_ref_spans[read_idx], event) {
            continue;
        }
        let ref_log10 =
            marginalize_haplotype_indices(read_haplotype_likelihoods, read_idx, &ref_haps);
        let alt_log10 =
            marginalize_haplotype_indices(read_haplotype_likelihoods, read_idx, &alt_haps);
        log10_likelihoods[0] += ref_log10;
        log10_likelihoods[1] += log10_sum_exp(&[ref_log10, alt_log10]) - 0.3010299956639812;
        log10_likelihoods[2] += alt_log10;
    }

    finalize_assembled_event_call(
        event,
        variant_model_from_log10(log10_likelihoods, assembled_event_log10_priors()),
        depth,
        ref_count,
        alt_count,
        fs,
        min_confidence,
    )
}

fn genotype_overlapping_assembled_event(
    local_haplotypes: &[LocalHaplotype],
    event_idx: usize,
    event: &VariantCall,
    read_haplotype_likelihoods: &[Vec<f64>],
    read_is_reverse_list: &[bool],
    read_ref_spans: &[(u64, u64)],
    min_confidence: f64,
    overlapping_event_indices: &[usize],
    pair_context: &PairGenotypingContext,
) -> Option<VariantCall> {
    let hap_contains_allele = event_haplotype_mask(local_haplotypes, event_idx);
    if hap_contains_allele.iter().all(|present| !present) {
        return None;
    }

    let mut log10_likelihoods = [f64::NEG_INFINITY; 3];
    for (pair_idx, pair) in pair_context.haplotype_pairs.iter().copied().enumerate() {
        let genotype_index = pair_event_genotype_index(pair, &hap_contains_allele);
        if genotype_index == 0
            && pair_contains_competing_overlapping_event(
                pair,
                local_haplotypes,
                overlapping_event_indices,
            )
        {
            continue;
        }
        log10_likelihoods[genotype_index] = log10_sum_exp(&[
            log10_likelihoods[genotype_index],
            pair_context.pair_log10_likelihoods[pair_idx],
        ]);
    }

    let genotype_index = pair_event_genotype_index(pair_context.best_pair, &hap_contains_allele);
    if genotype_index == 0 {
        return None;
    }

    let (best_ref_haps, best_alt_haps) =
        best_pair_hap_sets(pair_context.best_pair, &hap_contains_allele);
    let (depth, ref_count, alt_count, fs) = count_event_support(
        read_haplotype_likelihoods,
        read_is_reverse_list,
        read_ref_spans,
        event,
        &best_ref_haps,
        &best_alt_haps,
    );

    let mut model = variant_model_from_log10(log10_likelihoods, assembled_event_log10_priors());
    model.genotype_index = genotype_index;
    normalize_variant_model_pl(&mut model);

    finalize_assembled_event_call(
        event,
        model,
        depth,
        ref_count,
        alt_count,
        fs,
        min_confidence,
    )
}

fn event_haplotype_mask(local_haplotypes: &[LocalHaplotype], event_idx: usize) -> Vec<bool> {
    local_haplotypes
        .iter()
        .map(|hap| hap.event_indices.contains(&event_idx))
        .collect()
}

fn split_event_haplotype_indices(hap_contains_allele: &[bool]) -> (Vec<usize>, Vec<usize>) {
    let mut ref_haps = Vec::new();
    let mut alt_haps = Vec::new();
    for (hap_idx, contains_allele) in hap_contains_allele.iter().copied().enumerate() {
        if contains_allele {
            alt_haps.push(hap_idx);
        } else {
            ref_haps.push(hap_idx);
        }
    }
    (ref_haps, alt_haps)
}

fn count_event_support(
    read_haplotype_likelihoods: &[Vec<f64>],
    read_is_reverse_list: &[bool],
    read_ref_spans: &[(u64, u64)],
    event: &VariantCall,
    ref_haps: &[usize],
    alt_haps: &[usize],
) -> (u32, u32, u32, f64) {
    let mut ref_count = 0_u32;
    let mut alt_count = 0_u32;
    let mut ref_strand = StrandCounts::default();
    let mut alt_strand = StrandCounts::default();

    for read_idx in 0..read_haplotype_likelihoods.len() {
        if !read_span_overlaps_event(read_ref_spans[read_idx], event) {
            continue;
        }
        let ref_log10 =
            marginalize_haplotype_indices(read_haplotype_likelihoods, read_idx, ref_haps);
        let alt_log10 =
            marginalize_haplotype_indices(read_haplotype_likelihoods, read_idx, alt_haps);
        let is_reverse = read_is_reverse_list[read_idx];
        if ref_log10 - alt_log10 > 0.2 {
            ref_count += 1;
            ref_strand.increment(is_reverse);
        } else if alt_log10 - ref_log10 > 0.2 {
            alt_count += 1;
            alt_strand.increment(is_reverse);
        }
    }

    let depth = ref_count + alt_count;
    let fs = fisher_strand_score(ref_strand, alt_strand);
    (depth, ref_count, alt_count, fs)
}

fn finalize_assembled_event_call(
    event: &VariantCall,
    mut model: VariantModel,
    depth: u32,
    ref_count: u32,
    alt_count: u32,
    fs: f64,
    min_confidence: f64,
) -> Option<VariantCall> {
    normalize_variant_model_pl(&mut model);
    if model.genotype_index == 0 || f64::from(model.qual) < min_confidence {
        return None;
    }

    let mut final_call = event.clone();
    final_call.pl = model.pl;
    final_call.genotype_index = model.genotype_index;
    final_call.qual = model.qual;
    final_call.depth = depth;
    final_call.ref_count = ref_count;
    final_call.alt_count = alt_count;
    final_call.fs = fs;
    Some(final_call)
}

fn normalize_variant_model_pl(model: &mut VariantModel) {
    if model.pl[0] == 0 && model.pl[1] == 0 && model.pl[2] == 0 {
        for idx in 0..3 {
            if idx != model.genotype_index {
                model.pl[idx] = 9999;
            }
        }
    }
}

fn assembled_event_log10_priors() -> [f64; 3] {
    let heterozygosity: f64 = 1e-3;
    [
        (1.0 - 1.5 * heterozygosity).log10(),
        heterozygosity.log10(),
        (0.5 * heterozygosity).log10(),
    ]
}

fn overlapping_event_mask(valid_events: &[VariantCall]) -> Vec<bool> {
    let mut mask = vec![false; valid_events.len()];
    for left_idx in 0..valid_events.len() {
        for right_idx in left_idx + 1..valid_events.len() {
            if events_overlap(&valid_events[left_idx], &valid_events[right_idx]) {
                mask[left_idx] = true;
                mask[right_idx] = true;
            }
        }
    }
    mask
}

fn overlapping_event_indices(valid_events: &[VariantCall]) -> Vec<Vec<usize>> {
    let mut overlaps = vec![Vec::new(); valid_events.len()];
    for left_idx in 0..valid_events.len() {
        for right_idx in left_idx + 1..valid_events.len() {
            if events_overlap(&valid_events[left_idx], &valid_events[right_idx]) {
                overlaps[left_idx].push(right_idx);
                overlaps[right_idx].push(left_idx);
            }
        }
    }
    overlaps
}

fn events_overlap(left: &VariantCall, right: &VariantCall) -> bool {
    let left_end = left.pos + left.ref_allele.len() as u64 - 1;
    let right_end = right.pos + right.ref_allele.len() as u64 - 1;
    left.pos <= right_end && right.pos <= left_end
}

fn pair_contains_competing_overlapping_event(
    pair: (usize, usize),
    local_haplotypes: &[LocalHaplotype],
    overlapping_event_indices: &[usize],
) -> bool {
    if overlapping_event_indices.is_empty() {
        return false;
    }
    [pair.0, pair.1].into_iter().any(|hap_idx| {
        local_haplotypes[hap_idx]
            .event_indices
            .iter()
            .any(|event_idx| overlapping_event_indices.contains(event_idx))
    })
}

fn read_reference_span(record: &bam::Record) -> (u64, u64) {
    let start = (record.pos() + 1) as u64;
    let end = record.reference_end() as u64;
    (start, end)
}

fn read_reference_span_from_start_and_cigar(start: u64, cigar: &str) -> (u64, u64) {
    let mut reference_len = 0_u64;
    let mut current_len = 0_u64;
    for byte in cigar.bytes() {
        if byte.is_ascii_digit() {
            current_len = current_len
                .saturating_mul(10)
                .saturating_add(u64::from(byte - b'0'));
            continue;
        }
        match byte as char {
            'M' | 'D' | 'N' | '=' | 'X' => {
                reference_len = reference_len.saturating_add(current_len);
            }
            'I' | 'S' | 'H' | 'P' => {}
            _ => {}
        }
        current_len = 0;
    }
    let end = start.saturating_add(reference_len.saturating_sub(1));
    (start, end)
}

fn read_span_overlaps_event(read_span: (u64, u64), event: &VariantCall) -> bool {
    let event_end = event.pos + event.ref_allele.len() as u64 - 1;
    read_span.0 <= event_end && event.pos <= read_span.1
}

fn enumerate_haplotype_pairs(n_haplotypes: usize) -> Vec<(usize, usize)> {
    let mut pairs = Vec::new();
    for left in 0..n_haplotypes {
        for right in left..n_haplotypes {
            pairs.push((left, right));
        }
    }
    pairs
}

fn compute_pair_log10_likelihoods(
    haplotype_pairs: &[(usize, usize)],
    read_haplotype_likelihoods: &[Vec<f64>],
) -> Vec<f64> {
    let mut pair_log10_likelihoods = vec![0.0; haplotype_pairs.len()];
    for (pair_idx, (left, right)) in haplotype_pairs.iter().copied().enumerate() {
        let mut total = 0.0;
        for read_likelihoods in read_haplotype_likelihoods {
            total += if left == right {
                read_likelihoods[left]
            } else {
                log10_sum_exp(&[read_likelihoods[left], read_likelihoods[right]])
                    - 0.3010299956639812
            };
        }
        pair_log10_likelihoods[pair_idx] = total;
    }
    pair_log10_likelihoods
}

fn pair_event_genotype_index(pair: (usize, usize), hap_contains_allele: &[bool]) -> usize {
    let left = hap_contains_allele[pair.0];
    let right = hap_contains_allele[pair.1];
    match (left, right) {
        (false, false) => 0,
        (true, true) => 2,
        _ => 1,
    }
}

fn best_pair_hap_sets(
    pair: (usize, usize),
    hap_contains_allele: &[bool],
) -> (Vec<usize>, Vec<usize>) {
    let mut ref_haps = Vec::new();
    let mut alt_haps = Vec::new();
    for hap_idx in [pair.0, pair.1] {
        if hap_contains_allele[hap_idx] {
            if !alt_haps.contains(&hap_idx) {
                alt_haps.push(hap_idx);
            }
        } else if !ref_haps.contains(&hap_idx) {
            ref_haps.push(hap_idx);
        }
    }
    (ref_haps, alt_haps)
}

fn marginalize_haplotype_indices(
    read_haplotype_likelihoods: &[Vec<f64>],
    read_idx: usize,
    hap_indices: &[usize],
) -> f64 {
    if hap_indices.is_empty() {
        return f64::NEG_INFINITY;
    }
    let mut values = Vec::with_capacity(hap_indices.len());
    for hap_idx in hap_indices {
        values.push(read_haplotype_likelihoods[read_idx][*hap_idx]);
    }
    marginalize_allele_likelihoods(&values)
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

fn collect_pileup_fallback_events(
    contig: &str,
    region_start: u64,
    local_ref_bases: &[u8],
    active_loci: &[ReplayActiveLocusRow],
    min_qual: f64,
) -> Vec<VariantCall> {
    let region_end = region_start + local_ref_bases.len().saturating_sub(1) as u64;
    let mut events = Vec::new();
    for row in active_loci {
        if row.contig != contig || row.pos < region_start || row.pos > region_end {
            continue;
        }
        if row.depth < 10 {
            continue;
        }
        if row.snp_alt_count >= 3 {
            let alt_base = row.snp_best_alt.as_bytes().first().copied();
            let ref_base = base_index(row.ref_base);
            let alt_index = alt_base.and_then(base_index);
            let offset = (row.pos - region_start) as usize;
            if let (Some(ref_index), Some(alt_index)) = (ref_base, alt_index) {
                let mut observations = Vec::with_capacity(row.depth as usize);
                for _ in 0..row.depth.saturating_sub(row.snp_alt_count) {
                    observations.push(BaseObservation {
                        base_index: ref_index,
                        quality: 30,
                        is_reverse: false,
                    });
                }
                for _ in 0..row.snp_alt_count {
                    observations.push(BaseObservation {
                        base_index: alt_index,
                        quality: 30,
                        is_reverse: false,
                    });
                }
                let model = snp_variant_model(&observations, ref_index, alt_index);
                if f64::from(model.qual) >= min_qual {
                    let alt_base = base_from_index(alt_index);
                    if let Some(ref_byte) = local_ref_bases.get(offset).copied() {
                        events.push(VariantCall {
                            contig: contig.to_string(),
                            pos: row.pos,
                            id: None,
                            db: false,
                            ref_allele: vec![normalize_base(ref_byte)],
                            alt_allele: vec![alt_base],
                            depth: row.depth,
                            ref_count: row.depth.saturating_sub(row.snp_alt_count),
                            alt_count: row.snp_alt_count,
                            qual: model.qual,
                            fs: 0.0,
                            pl: model.pl,
                            genotype_index: model.genotype_index,
                        });
                    }
                }
            }
        }
        if row.indel_alt_count >= 3 {
            if let Some(label) = row.indel_best_alt.strip_prefix("INS:") {
                let inserted = label.as_bytes().to_vec();
                let offset = (row.pos - region_start) as usize;
                if let Some(ref_byte) = local_ref_bases.get(offset).copied() {
                    let mut observations = Vec::with_capacity(row.depth as usize);
                    for _ in 0..row.depth.saturating_sub(row.indel_alt_count) {
                        observations.push(IndelObservation {
                            allele: IndelObservationAllele::Ref,
                            quality: 30,
                            is_reverse: false,
                        });
                    }
                    let alt_allele = IndelAllele::Insertion(inserted.clone());
                    for _ in 0..row.indel_alt_count {
                        observations.push(IndelObservation {
                            allele: IndelObservationAllele::Alt(alt_allele.clone()),
                            quality: 30,
                            is_reverse: false,
                        });
                    }
                    let model = indel_variant_model(&observations, &alt_allele);
                    if f64::from(model.qual) >= min_qual {
                        let anchor = normalize_base(ref_byte);
                        let mut alt_bases = vec![anchor];
                        alt_bases.extend_from_slice(&inserted);
                        events.push(VariantCall {
                            contig: contig.to_string(),
                            pos: row.pos,
                            id: None,
                            db: false,
                            ref_allele: vec![anchor],
                            alt_allele: alt_bases,
                            depth: row.depth,
                            ref_count: row.depth.saturating_sub(row.indel_alt_count),
                            alt_count: row.indel_alt_count,
                            qual: model.qual,
                            fs: 0.0,
                            pl: model.pl,
                            genotype_index: model.genotype_index,
                        });
                    }
                }
            }
        }
    }
    events.sort_by(|a, b| {
        a.pos
            .cmp(&b.pos)
            .then_with(|| a.ref_allele.cmp(&b.ref_allele))
            .then_with(|| a.alt_allele.cmp(&b.alt_allele))
    });
    events.dedup_by(|a, b| {
        a.pos == b.pos && a.ref_allele == b.ref_allele && a.alt_allele == b.alt_allele
    });
    events
}

fn collect_zero_candidate_simple_snp_seed_events(
    contig: &str,
    region_start: u64,
    local_ref_bases: &[u8],
    active_loci: &[ReplayActiveLocusRow],
) -> Vec<VariantCall> {
    collect_pileup_fallback_events(contig, region_start, local_ref_bases, active_loci, 0.0)
        .into_iter()
        .filter(|event| event.ref_allele.len() == 1 && event.alt_allele.len() == 1)
        .filter(|event| event.genotype_index == 0)
        .filter(|event| {
            active_loci
                .iter()
                .any(|row| active_locus_exact_simple_snp_support_without_indel(row, event))
        })
        .collect()
}

fn same_event_key(left: &VariantCall, right: &VariantCall) -> bool {
    left.pos == right.pos
        && left.ref_allele == right.ref_allele
        && left.alt_allele == right.alt_allele
}

fn merge_supplemental_haplotype(
    local_haplotypes: &mut Vec<LocalHaplotype>,
    mut haplotype: LocalHaplotype,
) {
    if haplotype.is_ref {
        return;
    }
    haplotype.event_indices.sort_unstable();
    haplotype.event_indices.dedup();
    if let Some(existing) = local_haplotypes
        .iter_mut()
        .find(|existing| existing.bases == haplotype.bases)
    {
        existing
            .event_indices
            .extend(haplotype.event_indices.iter().copied());
        existing.event_indices.sort_unstable();
        existing.event_indices.dedup();
        return;
    }
    local_haplotypes.push(haplotype);
}

fn haplotype_base_index_for_reference_pos(
    region_start: u64,
    haplotype: &LocalHaplotype,
    ref_pos: u64,
) -> Option<usize> {
    let mut current_ref = region_start;
    let mut current_hap = 0_usize;
    let mut op_len = 0_usize;

    for byte in haplotype.cigar.bytes() {
        if byte.is_ascii_digit() {
            op_len = op_len
                .saturating_mul(10)
                .saturating_add(usize::from(byte - b'0'));
            continue;
        }

        match byte as char {
            'M' | '=' | 'X' => {
                let op_end = current_ref + op_len as u64;
                if current_ref <= ref_pos && ref_pos < op_end {
                    return Some(current_hap + (ref_pos - current_ref) as usize);
                }
                current_ref = op_end;
                current_hap += op_len;
            }
            'I' | 'S' => {
                current_hap += op_len;
            }
            'D' | 'N' => {
                let op_end = current_ref + op_len as u64;
                if current_ref <= ref_pos && ref_pos < op_end {
                    return None;
                }
                current_ref = op_end;
            }
            'H' | 'P' => {}
            _ => {}
        }
        op_len = 0;
    }

    None
}

fn overlay_supplemental_snp_on_haplotype(
    region_start: u64,
    haplotype: &LocalHaplotype,
    existing_events: &[VariantCall],
    event: &VariantCall,
    event_idx: usize,
) -> Option<LocalHaplotype> {
    if haplotype.is_ref || event.ref_allele.len() != 1 || event.alt_allele.len() != 1 {
        return None;
    }
    if !haplotype.event_indices.iter().all(|existing_idx| {
        existing_events
            .get(*existing_idx)
            .is_some_and(|existing_event| {
                existing_event.ref_allele.len() == 1
                    && existing_event.alt_allele.len() == 1
                    && !events_overlap(existing_event, event)
            })
    }) {
        return None;
    }

    let base_idx = haplotype_base_index_for_reference_pos(region_start, haplotype, event.pos)?;
    if haplotype.bases.get(base_idx).copied()? != event.ref_allele[0] {
        return None;
    }

    let mut overlaid = haplotype.clone();
    overlaid.bases[base_idx] = event.alt_allele[0];
    overlaid.event_indices.push(event_idx);
    Some(overlaid)
}

fn supplement_missing_pileup_events(
    contig: &str,
    region_start: u64,
    local_ref_bases: &[u8],
    active_loci: &[ReplayActiveLocusRow],
    min_qual: f64,
    local_haplotypes: &mut Vec<LocalHaplotype>,
    valid_events: &mut Vec<VariantCall>,
) {
    let mut fallback_events = collect_pileup_fallback_events(
        contig,
        region_start,
        local_ref_bases,
        active_loci,
        min_qual,
    );
    if fallback_events.is_empty() && valid_events.is_empty() {
        // In some RNA-seq repeat contexts the pileup model under-scores a real
        // simple SNP, but PairHMM can recover it once the candidate exists.
        // Only relax seeding for zero-candidate regions, simple SNPs, and
        // loci without competing indel support.
        fallback_events = collect_zero_candidate_simple_snp_seed_events(
            contig,
            region_start,
            local_ref_bases,
            active_loci,
        );
    }
    if fallback_events.is_empty() {
        return;
    }

    if valid_events.is_empty() {
        *local_haplotypes = haplotypes_from_candidate_events(
            contig,
            region_start,
            local_ref_bases,
            &fallback_events,
        );
        *valid_events = fallback_events;
        return;
    }

    // Keep the original zero-candidate fallback behavior for both SNPs and
    // indels, but narrow the non-empty-region supplement path to simple SNPs.
    // Round8 full-call evidence showed that the Java-only gain here was SNP
    // dominated, while the new Rust-only regression included a large added
    // indel class.
    fallback_events.retain(|event| event.ref_allele.len() == 1 && event.alt_allele.len() == 1);
    if fallback_events.is_empty() {
        return;
    }

    let missing_events: Vec<VariantCall> = fallback_events
        .into_iter()
        .filter(|event| {
            !valid_events
                .iter()
                .any(|existing| same_event_key(existing, event))
                && !should_skip_weak_supplemental_snp_in_dense_snp_cluster(event, valid_events)
        })
        .collect();
    if missing_events.is_empty() {
        return;
    }

    let base_event_idx = valid_events.len();
    let existing_alt_haplotypes: Vec<LocalHaplotype> = local_haplotypes
        .iter()
        .filter(|hap| !hap.is_ref)
        .cloned()
        .collect();

    // Single-event synthetic haplotypes can underfit nearby multi-SNP reads.
    // When a pileup-strong missing SNP sits on the same reads as an already
    // assembled nearby ALT haplotype, also overlay that SNP onto the existing
    // ALT haplotype so PairHMM can score the combined sequence.
    for (event_offset, event) in missing_events.iter().enumerate() {
        let event_idx = base_event_idx + event_offset;
        for haplotype in &existing_alt_haplotypes {
            if let Some(overlaid) = overlay_supplemental_snp_on_haplotype(
                region_start,
                haplotype,
                valid_events,
                event,
                event_idx,
            ) {
                merge_supplemental_haplotype(local_haplotypes, overlaid);
            }
        }
    }

    let supplemental_haplotypes =
        haplotypes_from_candidate_events(contig, region_start, local_ref_bases, &missing_events);
    valid_events.extend(missing_events.iter().cloned());

    for mut haplotype in supplemental_haplotypes
        .into_iter()
        .filter(|hap| !hap.is_ref)
    {
        for event_idx in &mut haplotype.event_indices {
            *event_idx += base_event_idx;
        }
        merge_supplemental_haplotype(local_haplotypes, haplotype);
    }
}

fn should_skip_weak_supplemental_snp_in_dense_snp_cluster(
    event: &VariantCall,
    valid_events: &[VariantCall],
) -> bool {
    if event.ref_allele.len() != 1
        || event.alt_allele.len() != 1
        || event.alt_count > WEAK_SUPPLEMENTAL_CLUSTER_SNP_MAX_ALT_COUNT
    {
        return false;
    }

    let mut nearby_positions = Vec::with_capacity(2);
    for existing in valid_events {
        if existing.ref_allele.len() != 1 || existing.alt_allele.len() != 1 {
            continue;
        }
        if existing.pos.abs_diff(event.pos) > SNP_CLUSTER_WINDOW {
            continue;
        }
        if !nearby_positions.contains(&existing.pos) {
            nearby_positions.push(existing.pos);
            if nearby_positions.len() >= 2 {
                return true;
            }
        }
    }

    false
}

fn rescue_collapsed_strong_snp_cluster_from_pileup(
    contig: &str,
    region_start: u64,
    local_ref_bases: &[u8],
    active_loci: &[ReplayActiveLocusRow],
    valid_events: &[VariantCall],
    min_qual: f64,
) -> Vec<VariantCall> {
    let fallback_events = collect_pileup_fallback_events(
        contig,
        region_start,
        local_ref_bases,
        active_loci,
        min_qual,
    );
    let rescued = exact_strong_simple_snp_pileup_matches_from_fallback_events(
        &fallback_events,
        active_loci,
        valid_events,
    );

    if rescued.len() >= 2 {
        return rescued;
    }

    if rescued.len() == 1 && fallback_events.len() == 1 {
        return rescued;
    }

    if rescued.len() == 1 && has_high_confidence_single_snp_rescue_support(active_loci, &rescued[0])
    {
        return rescued;
    }

    Vec::new()
}

fn exact_strong_simple_snp_pileup_matches_from_fallback_events(
    fallback_events: &[VariantCall],
    active_loci: &[ReplayActiveLocusRow],
    valid_events: &[VariantCall],
) -> Vec<VariantCall> {
    fallback_events
        .iter()
        .filter(|event| event.ref_allele.len() == 1 && event.alt_allele.len() == 1)
        .filter(|event| {
            valid_events
                .iter()
                .any(|valid| same_event_key(valid, event))
                && active_loci
                    .iter()
                    .any(|row| active_locus_exact_simple_snp_support_without_indel(row, event))
        })
        .cloned()
        .collect()
}

fn exact_strong_simple_snp_pileup_matches(
    contig: &str,
    region_start: u64,
    local_ref_bases: &[u8],
    active_loci: &[ReplayActiveLocusRow],
    valid_events: &[VariantCall],
    min_qual: f64,
) -> Vec<VariantCall> {
    let fallback_events = collect_pileup_fallback_events(
        contig,
        region_start,
        local_ref_bases,
        active_loci,
        min_qual,
    );
    exact_strong_simple_snp_pileup_matches_from_fallback_events(
        &fallback_events,
        active_loci,
        valid_events,
    )
}

fn active_locus_exact_simple_snp_support(row: &ReplayActiveLocusRow, event: &VariantCall) -> bool {
    event.ref_allele.len() == 1
        && event.alt_allele.len() == 1
        && event.pos == row.pos
        && event.ref_allele[0] == row.ref_base
        && row
            .snp_best_alt
            .as_bytes()
            .first()
            .is_some_and(|alt| event.alt_allele[0] == *alt)
}

fn active_locus_exact_simple_snp_support_without_indel(
    row: &ReplayActiveLocusRow,
    event: &VariantCall,
) -> bool {
    row.indel_alt_count == 0 && active_locus_exact_simple_snp_support(row, event)
}

fn active_locus_high_confidence_single_snp_rescue_support(
    row: &ReplayActiveLocusRow,
    event: &VariantCall,
) -> bool {
    active_locus_exact_simple_snp_support_without_indel(row, event)
        && row.snp_alt_count >= STRONG_SINGLE_SNP_RESCUE_MIN_ALT_COUNT
        && row.alt_fraction >= STRONG_SINGLE_SNP_RESCUE_MIN_ALT_FRACTION
}

fn has_high_confidence_single_snp_rescue_support(
    active_loci: &[ReplayActiveLocusRow],
    event: &VariantCall,
) -> bool {
    active_loci
        .iter()
        .any(|row| active_locus_high_confidence_single_snp_rescue_support(row, event))
}

fn prune_unsupported_simple_snp_calls_in_dense_clusters(
    final_calls: &mut Vec<VariantCall>,
    active_loci: &[ReplayActiveLocusRow],
) {
    if final_calls.len() < 3 {
        return;
    }

    let dense_unsupported_snp_keys: Vec<(u64, Vec<u8>, Vec<u8>)> = final_calls
        .iter()
        .filter(|call| call.ref_allele.len() == 1 && call.alt_allele.len() == 1)
        .filter(|call| {
            !active_loci
                .iter()
                .any(|row| active_locus_exact_simple_snp_support(row, call))
        })
        .filter(|call| {
            let mut nearby_positions = Vec::with_capacity(2);
            for other in final_calls.iter() {
                if same_event_key(call, other)
                    || other.ref_allele.len() != 1
                    || other.alt_allele.len() != 1
                    || other.pos.abs_diff(call.pos) > SNP_CLUSTER_WINDOW
                {
                    continue;
                }
                if !nearby_positions.contains(&other.pos) {
                    nearby_positions.push(other.pos);
                    if nearby_positions.len() >= 2 {
                        return true;
                    }
                }
            }
            false
        })
        .map(|call| (call.pos, call.ref_allele.clone(), call.alt_allele.clone()))
        .collect();

    if dense_unsupported_snp_keys.is_empty() {
        return;
    }

    final_calls.retain(|call| {
        !dense_unsupported_snp_keys
            .iter()
            .any(|(pos, ref_allele, alt_allele)| {
                call.pos == *pos && call.ref_allele == *ref_allele && call.alt_allele == *alt_allele
            })
    });
}

fn merge_missing_strong_snp_cluster_rescues_from_pileup(
    final_calls: &mut Vec<VariantCall>,
    contig: &str,
    region_start: u64,
    local_ref_bases: &[u8],
    active_loci: &[ReplayActiveLocusRow],
    valid_events: &[VariantCall],
    min_qual: f64,
) {
    let rescued = exact_strong_simple_snp_pileup_matches(
        contig,
        region_start,
        local_ref_bases,
        active_loci,
        valid_events,
        min_qual,
    );
    if rescued.len() < 2
        && !rescued
            .first()
            .is_some_and(|event| has_high_confidence_single_snp_rescue_support(active_loci, event))
    {
        return;
    }

    for rescued_call in rescued {
        if !final_calls
            .iter()
            .any(|final_call| same_event_key(final_call, &rescued_call))
        {
            final_calls.push(rescued_call);
        }
    }
}

fn haplotypes_from_candidate_events(
    contig: &str,
    region_start: u64,
    local_ref_bases: &[u8],
    candidate_events: &[VariantCall],
) -> Vec<LocalHaplotype> {
    let region_end = region_start + local_ref_bases.len().saturating_sub(1) as u64;
    build_local_haplotypes(
        contig,
        region_start,
        region_end,
        local_ref_bases,
        candidate_events,
        128,
    )
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

fn marginalize_allele_likelihoods(values: &[f64]) -> f64 {
    if values.is_empty() {
        f64::NEG_INFINITY
    } else {
        values.iter().copied().fold(f64::NEG_INFINITY, f64::max)
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

fn padded_fetch_bounds(interval: &Interval, dict: &SequenceDict) -> (u64, u64) {
    let contig_len = dict
        .contig_length(&interval.contig)
        .expect("interval contig validated by sort_intervals");
    let start = interval.start.saturating_sub(ACTIVE_REGION_PADDING).max(1);
    let end = interval
        .end
        .saturating_add(ACTIVE_REGION_PADDING)
        .min(contig_len);
    (start, end)
}

fn coalesce_fetch_windows(intervals: &[Interval], dict: &SequenceDict) -> Vec<FetchWindow> {
    let mut windows: Vec<FetchWindow> = Vec::new();
    for interval in intervals {
        let (window_start, window_end) = padded_fetch_bounds(interval, dict);
        if let Some(current) = windows.last_mut() {
            let same_contig = current.contig == interval.contig;
            let close_enough = window_start <= current.end.saturating_add(FETCH_WINDOW_GAP + 1);
            let merged_len = window_end.saturating_sub(current.start).saturating_add(1);
            if same_contig && close_enough && merged_len <= FETCH_WINDOW_MAX_BASES {
                current.end = current.end.max(window_end);
                current.intervals.push(interval.clone());
                continue;
            }
        }
        windows.push(FetchWindow {
            contig: interval.contig.clone(),
            start: window_start,
            end: window_end,
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
        "contig\tstart\tend\tregion\tactive_start\tactive_end\tpadded_start\tpadded_end\tobserved_loci\tactive_loci\tcandidate_events\tmax_alt_fraction\tmean_alt_fraction"
    )?;
    for row in &output.active_regions {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.6}\t{:.6}",
            row.contig,
            row.start,
            row.end,
            row.region,
            row.active_start,
            row.active_end,
            row.padded_start,
            row.padded_end,
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
                push_cigar(
                    &mut cigar_ops,
                    (ev.alt_allele.len() - ev.ref_allele.len()) as u32,
                    'I',
                );
            } else if ev.ref_allele.len() > ev.alt_allele.len() {
                push_cigar(
                    &mut cigar_ops,
                    (ev.ref_allele.len() - ev.alt_allele.len()) as u32,
                    'D',
                );
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

    fn test_bam_record(flags: u16, cigar: Vec<rust_htslib::bam::record::Cigar>) -> bam::Record {
        let cigar = rust_htslib::bam::record::CigarString(cigar);
        let read_len = cigar
            .iter()
            .map(|op| match op {
                rust_htslib::bam::record::Cigar::Match(len)
                | rust_htslib::bam::record::Cigar::Equal(len)
                | rust_htslib::bam::record::Cigar::Diff(len)
                | rust_htslib::bam::record::Cigar::Ins(len)
                | rust_htslib::bam::record::Cigar::SoftClip(len) => *len as usize,
                _ => 0,
            })
            .sum::<usize>();
        let bases = vec![b'A'; read_len];
        let quals = vec![30_u8; read_len];
        let mut record = bam::Record::new();
        record.set(b"read1", Some(&cigar), &bases, &quals);
        record.set_flags(flags);
        record.set_mapq(60);
        record
    }

    fn test_bam_record_with_bases_quals(
        cigar: Vec<rust_htslib::bam::record::Cigar>,
        bases: &[u8],
        quals: &[u8],
        pos0: i64,
    ) -> bam::Record {
        let cigar = rust_htslib::bam::record::CigarString(cigar);
        let mut record = bam::Record::new();
        record.set(b"read1", Some(&cigar), bases, quals);
        record.set_flags(0);
        record.set_mapq(60);
        record.set_pos(pos0);
        record
    }

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
            evidence.active_observations.push(ActiveBaseObservation {
                base_index: Some(ref_index),
                quality,
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
            evidence.active_observations.push(ActiveBaseObservation {
                base_index: Some(alt_index),
                quality,
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
        let dict = parse_dict_lines(
            &[
                "@HD\tVN:1.6\tSO:coordinate".to_string(),
                "@SQ\tSN:chr1\tLN:5000".to_string(),
            ],
            Path::new("test.interval_list"),
        )
        .unwrap();
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
        let windows = coalesce_fetch_windows(&intervals, &dict);
        assert_eq!(windows.len(), 2);
        assert_eq!(windows[0].start, 1);
        assert_eq!(windows[0].end, 130);
        assert_eq!(windows[0].intervals.len(), 2);
        assert_eq!(windows[1].start, 1_900);
        assert_eq!(windows[1].end, 2_110);
    }

    #[test]
    fn fetch_windows_pad_requested_intervals_and_clip_to_contig_edges() {
        let dict = test_dict();
        let left_windows = coalesce_fetch_windows(
            &[Interval {
                contig: "chr1".to_string(),
                start: 10,
                end: 20,
            }],
            &dict,
        );
        assert_eq!(left_windows.len(), 1);
        assert_eq!(left_windows[0].start, 1);
        assert_eq!(left_windows[0].end, 120);

        let right_windows = coalesce_fetch_windows(
            &[Interval {
                contig: "chr1".to_string(),
                start: 181,
                end: 200,
            }],
            &dict,
        );
        assert_eq!(right_windows.len(), 1);
        assert_eq!(right_windows[0].start, 81);
        assert_eq!(right_windows[0].end, 200);
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
    fn select_assembly_reads_keeps_only_flagged_segments_in_input_order() {
        let reads = vec![
            vec![b"AAA".to_vec(), b"TTT".to_vec()],
            vec![b"CCC".to_vec()],
            vec![b"GGG".to_vec()],
        ];
        let selected = select_assembly_reads(&reads, &[true, false, true]);
        let selected_owned: Vec<Vec<u8>> = selected.into_iter().cloned().collect();
        assert_eq!(
            selected_owned,
            vec![b"AAA".to_vec(), b"TTT".to_vec(), b"GGG".to_vec()]
        );
    }

    #[test]
    fn assembly_read_segments_split_on_n_and_low_base_quality() {
        let segments = assembly_read_segments(
            b"AAAANCCCTGGG",
            &[30, 30, 30, 30, 30, 30, 30, 30, 9, 30, 30, 30],
        );

        assert_eq!(
            segments,
            vec![b"AAAA".to_vec(), b"CCC".to_vec(), b"GGG".to_vec()]
        );
    }

    #[test]
    fn allele_marginalization_uses_max_across_supporting_haplotypes() {
        let marginalized = marginalize_allele_likelihoods(&[-2.0, -2.0]);
        assert_eq!(marginalized, -2.0);
    }

    #[test]
    fn allele_marginalization_returns_negative_infinity_for_empty_support() {
        assert_eq!(marginalize_allele_likelihoods(&[]), f64::NEG_INFINITY);
    }

    #[test]
    fn pair_hmm_base_quality_caps_by_mapq_then_squashes_low_values() {
        assert_eq!(pair_hmm_base_quality(30, 40), 30);
        assert_eq!(pair_hmm_base_quality(30, 20), 20);
        assert_eq!(pair_hmm_base_quality(17, 40), PAIR_HMM_MIN_USABLE_Q_SCORE);
        assert_eq!(pair_hmm_base_quality(30, 10), PAIR_HMM_MIN_USABLE_Q_SCORE);
        assert_eq!(pair_hmm_base_quality(18, 40), 18);
    }

    #[test]
    fn pair_hmm_indel_open_quality_uses_min_usable_floor() {
        assert_eq!(pair_hmm_indel_open_quality(4), PAIR_HMM_MIN_USABLE_Q_SCORE);
        assert_eq!(pair_hmm_indel_open_quality(6), 6);
        assert_eq!(pair_hmm_indel_open_quality(25), 25);
    }

    #[test]
    fn pad_active_region_clips_to_fetch_window_and_preserves_active_span() {
        let window = FetchWindow {
            contig: "chr1".to_string(),
            start: 100,
            end: 200,
            intervals: vec![Interval {
                contig: "chr1".to_string(),
                start: 100,
                end: 200,
            }],
        };
        let interval = Interval {
            contig: "chr1".to_string(),
            start: 110,
            end: 150,
        };
        let padded = pad_active_region(&interval, &window);
        assert_eq!(padded.start, 100);
        assert_eq!(padded.end, 200);
        assert_eq!(interval.start, 110);
        assert_eq!(interval.end, 150);
    }

    #[test]
    fn genotype_assembled_events_emits_final_calls_from_pairhmm_matrix() {
        let local_haplotypes = vec![
            LocalHaplotype {
                bases: b"AAA".to_vec(),
                is_ref: true,
                cigar: "3M".to_string(),
                event_indices: vec![],
            },
            LocalHaplotype {
                bases: b"ACA".to_vec(),
                is_ref: false,
                cigar: "3M".to_string(),
                event_indices: vec![0],
            },
        ];
        let valid_events = vec![test_variant("chr1", 10, b"A", b"C")];
        let read_haplotype_likelihoods = vec![vec![-10.0, 0.0], vec![-10.0, 0.0], vec![-10.0, 0.0]];
        let read_is_reverse_list = vec![false, true, false];
        let read_ref_spans = vec![(10, 10), (10, 10), (10, 10)];

        let final_calls = genotype_assembled_events(
            &local_haplotypes,
            &valid_events,
            &read_haplotype_likelihoods,
            &read_is_reverse_list,
            &read_ref_spans,
            20.0,
        );

        assert_eq!(final_calls.len(), 1);
        let call = &final_calls[0];
        assert_eq!(call.genotype(), "1/1");
        assert_eq!(call.depth, 3);
        assert_eq!(call.ref_count, 0);
        assert_eq!(call.alt_count, 3);
        assert!(call.qual >= 20);
    }

    #[test]
    fn genotype_assembled_events_ignores_overlapping_alt_haplotypes_as_ref_evidence() {
        let local_haplotypes = vec![
            LocalHaplotype {
                bases: b"AAA".to_vec(),
                is_ref: true,
                cigar: "3M".to_string(),
                event_indices: vec![],
            },
            LocalHaplotype {
                bases: b"ACA".to_vec(),
                is_ref: false,
                cigar: "3M".to_string(),
                event_indices: vec![0],
            },
            LocalHaplotype {
                bases: b"AGA".to_vec(),
                is_ref: false,
                cigar: "3M".to_string(),
                event_indices: vec![1],
            },
        ];
        let valid_events = vec![
            test_variant("chr1", 10, b"A", b"C"),
            test_variant("chr1", 10, b"A", b"G"),
        ];
        let mut read_haplotype_likelihoods = Vec::new();
        read_haplotype_likelihoods.extend(vec![vec![-10.0, 0.0, -10.0]; 2]);
        read_haplotype_likelihoods.extend(vec![vec![-10.0, -10.0, 0.0]; 8]);
        let read_is_reverse_list = vec![false; read_haplotype_likelihoods.len()];
        let read_ref_spans = vec![(10, 10); read_haplotype_likelihoods.len()];

        let final_calls = genotype_assembled_events(
            &local_haplotypes,
            &valid_events,
            &read_haplotype_likelihoods,
            &read_is_reverse_list,
            &read_ref_spans,
            0.0,
        );

        let c_call = final_calls
            .iter()
            .find(|call| call.ref_allele == b"A" && call.alt_allele == b"C")
            .unwrap();
        let g_call = final_calls
            .iter()
            .find(|call| call.ref_allele == b"A" && call.alt_allele == b"G")
            .unwrap();
        assert_ne!(c_call.genotype(), "0/0");
        assert_ne!(g_call.genotype(), "0/0");
    }

    #[test]
    fn genotype_assembled_events_pair_genotyping_ignores_competing_overlap_pairs_as_ref() {
        let local_haplotypes = vec![
            LocalHaplotype {
                bases: b"AAAA".to_vec(),
                is_ref: true,
                cigar: "4M".to_string(),
                event_indices: vec![],
            },
            LocalHaplotype {
                bases: b"ACAA".to_vec(),
                is_ref: false,
                cigar: "4M".to_string(),
                event_indices: vec![0],
            },
            LocalHaplotype {
                bases: b"AGAA".to_vec(),
                is_ref: false,
                cigar: "4M".to_string(),
                event_indices: vec![1],
            },
            LocalHaplotype {
                bases: b"ACCA".to_vec(),
                is_ref: false,
                cigar: "4M".to_string(),
                event_indices: vec![0],
            },
            LocalHaplotype {
                bases: b"ACGA".to_vec(),
                is_ref: false,
                cigar: "4M".to_string(),
                event_indices: vec![0],
            },
        ];
        let valid_events = vec![
            test_variant("chr1", 10, b"A", b"C"),
            test_variant("chr1", 10, b"A", b"G"),
        ];
        let read_haplotype_likelihoods = vec![
            vec![-10.0, 0.0, -12.0, -1.0, -1.0],
            vec![-10.0, -12.0, 0.0, -12.0, -12.0],
            vec![0.0, -10.0, -10.0, -10.0, -10.0],
        ];
        let read_is_reverse_list = vec![false; read_haplotype_likelihoods.len()];
        let read_ref_spans = vec![(10, 10); read_haplotype_likelihoods.len()];

        let final_calls = genotype_assembled_events(
            &local_haplotypes,
            &valid_events,
            &read_haplotype_likelihoods,
            &read_is_reverse_list,
            &read_ref_spans,
            0.0,
        );

        let c_call = final_calls
            .iter()
            .find(|call| call.ref_allele == b"A" && call.alt_allele == b"C")
            .unwrap();
        let g_call = final_calls
            .iter()
            .find(|call| call.ref_allele == b"A" && call.alt_allele == b"G")
            .unwrap();
        assert_eq!(c_call.genotype(), "0/1");
        assert_eq!(g_call.genotype(), "0/1");
    }

    #[test]
    fn genotype_assembled_events_uses_isolated_event_counts_for_simple_snp() {
        let local_haplotypes = vec![
            LocalHaplotype {
                bases: b"AAAA".to_vec(),
                is_ref: true,
                cigar: "4M".to_string(),
                event_indices: vec![],
            },
            LocalHaplotype {
                bases: b"ACAA".to_vec(),
                is_ref: false,
                cigar: "4M".to_string(),
                event_indices: vec![0],
            },
            LocalHaplotype {
                bases: b"AAGA".to_vec(),
                is_ref: false,
                cigar: "4M".to_string(),
                event_indices: vec![1],
            },
        ];
        let valid_events = vec![
            test_variant("chr1", 10, b"A", b"C"),
            test_variant("chr1", 12, b"A", b"G"),
        ];
        let read_haplotype_likelihoods = vec![
            vec![-10.0, 0.0, -12.0],
            vec![-10.0, 0.0, -12.0],
            vec![-10.0, -12.0, 0.0],
            vec![-10.0, -12.0, 0.0],
        ];
        let read_is_reverse_list = vec![false, true, false, true];
        let read_ref_spans = vec![(10, 10), (10, 10), (12, 12), (12, 12)];

        let final_calls = genotype_assembled_events(
            &local_haplotypes,
            &valid_events,
            &read_haplotype_likelihoods,
            &read_is_reverse_list,
            &read_ref_spans,
            0.0,
        );

        let c_call = final_calls
            .iter()
            .find(|call| call.ref_allele == b"A" && call.alt_allele == b"C")
            .unwrap();
        let g_call = final_calls
            .iter()
            .find(|call| call.ref_allele == b"A" && call.alt_allele == b"G")
            .unwrap();
        assert_eq!(c_call.depth, 2);
        assert_eq!(c_call.alt_count, 2);
        assert_eq!(g_call.depth, 2);
        assert_eq!(g_call.alt_count, 2);
    }

    #[test]
    fn genotype_assembled_events_ignores_reads_that_do_not_span_the_event() {
        let local_haplotypes = vec![
            LocalHaplotype {
                bases: b"AAAA".to_vec(),
                is_ref: true,
                cigar: "4M".to_string(),
                event_indices: vec![],
            },
            LocalHaplotype {
                bases: b"ACGA".to_vec(),
                is_ref: false,
                cigar: "4M".to_string(),
                event_indices: vec![0, 1],
            },
        ];
        let valid_events = vec![
            test_variant("chr1", 10, b"A", b"C"),
            test_variant("chr1", 12, b"A", b"G"),
        ];
        let read_haplotype_likelihoods = vec![vec![-10.0, 0.0], vec![-10.0, 0.0], vec![-10.0, 0.0]];
        let read_is_reverse_list = vec![false, false, true];
        let read_ref_spans = vec![(10, 10), (12, 12), (12, 12)];

        let final_calls = genotype_assembled_events(
            &local_haplotypes,
            &valid_events,
            &read_haplotype_likelihoods,
            &read_is_reverse_list,
            &read_ref_spans,
            0.0,
        );

        let c_call = final_calls
            .iter()
            .find(|call| call.ref_allele == b"A" && call.alt_allele == b"C")
            .unwrap();
        let g_call = final_calls
            .iter()
            .find(|call| call.ref_allele == b"A" && call.alt_allele == b"G")
            .unwrap();
        assert_eq!(c_call.depth, 1);
        assert_eq!(c_call.alt_count, 1);
        assert_eq!(g_call.depth, 2);
        assert_eq!(g_call.alt_count, 2);
    }

    #[test]
    fn genotype_assembled_events_marks_same_position_events_for_pair_genotyping() {
        let valid_events = vec![
            test_variant("chr1", 10, b"A", b"C"),
            test_variant("chr1", 10, b"A", b"G"),
            test_variant("chr1", 12, b"T", b"C"),
        ];

        assert_eq!(
            overlapping_event_mask(&valid_events),
            vec![true, true, false]
        );
    }

    #[test]
    fn read_reference_span_from_start_and_cigar_uses_reference_consuming_ops_only() {
        assert_eq!(
            read_reference_span_from_start_and_cigar(101, "10M1I10M"),
            (101, 120)
        );
        assert_eq!(
            read_reference_span_from_start_and_cigar(101, "10M100N10M1D5M"),
            (101, 226)
        );
    }

    #[test]
    fn prepare_hmm_read_trims_low_quality_tails_and_updates_ref_span() {
        let record = test_bam_record_with_bases_quals(
            vec![rust_htslib::bam::record::Cigar::Match(12)],
            b"AACCGGTTAACC",
            &[5, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 5],
            99,
        );

        let prepared = prepare_hmm_read(&record, 10, true, 100, 111).unwrap();

        assert_eq!(prepared.bases, b"ACCGGTTAAC");
        assert_eq!(prepared.ref_span, (101, 110));
        assert_eq!(prepared.quals, vec![30; 10]);
        assert_eq!(prepared.assembly_segments, vec![b"ACCGGTTAAC".to_vec()]);
    }

    #[test]
    fn prepare_hmm_read_excludes_soft_clips_from_bases_and_ref_span() {
        let record = test_bam_record_with_bases_quals(
            vec![
                rust_htslib::bam::record::Cigar::SoftClip(2),
                rust_htslib::bam::record::Cigar::Match(10),
                rust_htslib::bam::record::Cigar::SoftClip(2),
            ],
            b"TTACCGGTTAACAA",
            &[30; 14],
            199,
        );

        let prepared = prepare_hmm_read(&record, 10, true, 200, 209).unwrap();

        assert_eq!(prepared.bases, b"ACCGGTTAAC");
        assert_eq!(prepared.ref_span, (200, 209));
        assert_eq!(prepared.quals, vec![30; 10]);
    }

    #[test]
    fn prepare_hmm_read_clips_match_bases_to_region_span() {
        let record = test_bam_record_with_bases_quals(
            vec![rust_htslib::bam::record::Cigar::Match(14)],
            b"AACCGGTTAACCGG",
            &[30; 14],
            99,
        );

        let prepared = prepare_hmm_read(&record, 10, true, 103, 112).unwrap();

        assert_eq!(prepared.bases, b"CGGTTAACCG");
        assert_eq!(prepared.ref_span, (103, 112));
        assert_eq!(prepared.quals, vec![30; 10]);
    }

    #[test]
    fn prepare_hmm_read_keeps_boundary_insertion_when_clipping_left_tail() {
        let record = test_bam_record_with_bases_quals(
            vec![
                rust_htslib::bam::record::Cigar::Match(3),
                rust_htslib::bam::record::Cigar::Ins(2),
                rust_htslib::bam::record::Cigar::Match(8),
            ],
            b"AAAGGTTTTTTTT",
            &[30; 13],
            99,
        );

        let prepared = prepare_hmm_read(&record, 10, true, 103, 110).unwrap();

        assert_eq!(prepared.bases, b"GGTTTTTTTT");
        assert_eq!(prepared.ref_span, (103, 110));
        assert_eq!(prepared.quals, vec![30; 10]);
    }

    #[test]
    fn hc_filter_includes_supplementary_reads_by_default() {
        let record = test_bam_record(0x800, vec![rust_htslib::bam::record::Cigar::Match(5)]);
        assert!(read_passes_hc_filter(&record, 20, false));
    }

    #[test]
    fn hc_filter_can_exclude_supplementary_reads_for_debugging() {
        let record = test_bam_record(0x800, vec![rust_htslib::bam::record::Cigar::Match(5)]);
        assert!(!read_passes_hc_filter(&record, 20, true));
    }

    #[test]
    fn hc_filter_rejects_reads_without_reference_consuming_cigar_ops() {
        let record = test_bam_record(0, vec![rust_htslib::bam::record::Cigar::SoftClip(5)]);
        assert!(!read_passes_hc_filter(&record, 20, false));
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
    fn active_locus_treats_non_acgt_bases_as_non_ref_only_for_discovery() {
        let mut evidence = SnpEvidence::default();
        evidence.active_observations.push(ActiveBaseObservation {
            base_index: None,
            quality: 30,
        });
        evidence.active_observations.push(ActiveBaseObservation {
            base_index: None,
            quality: 30,
        });

        let (active, qual) = is_active_locus(Some(0), &evidence, evidence.counts.depth);
        assert!(active);
        assert!(qual >= ACTIVE_REGION_DISCOVERY_CONFIDENCE as u32);
        assert_eq!(evidence.counts.depth, 0);
        assert!(best_snp_alt(Some(0), &evidence).is_none());
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

    #[test]
    fn pileup_fallback_collects_strong_snp_candidate() {
        let active_loci = vec![ReplayActiveLocusRow {
            contig: "chr1".to_string(),
            pos: 101,
            region: "chr1:100-102".to_string(),
            ref_base: b'A',
            depth: 16,
            snp_alt_count: 8,
            snp_best_alt: "G".to_string(),
            indel_alt_count: 0,
            indel_best_alt: String::new(),
            alt_fraction: 0.5,
            active_probability_proxy: 1.0,
        }];
        let events = collect_pileup_fallback_events("chr1", 100, b"AAT", &active_loci, 20.0);
        assert_eq!(events.len(), 1);
        let event = &events[0];
        assert_eq!(event.pos, 101);
        assert_eq!(event.ref_allele, b"A");
        assert_eq!(event.alt_allele, b"G");
        assert!(event.qual >= 20);
    }

    #[test]
    fn pileup_fallback_requires_strong_support() {
        let active_loci = vec![
            ReplayActiveLocusRow {
                contig: "chr1".to_string(),
                pos: 101,
                region: "chr1:100-102".to_string(),
                ref_base: b'A',
                depth: 9,
                snp_alt_count: 8,
                snp_best_alt: "G".to_string(),
                indel_alt_count: 0,
                indel_best_alt: String::new(),
                alt_fraction: 0.5,
                active_probability_proxy: 1.0,
            },
            ReplayActiveLocusRow {
                contig: "chr1".to_string(),
                pos: 102,
                region: "chr1:100-102".to_string(),
                ref_base: b'A',
                depth: 16,
                snp_alt_count: 2,
                snp_best_alt: "G".to_string(),
                indel_alt_count: 0,
                indel_best_alt: String::new(),
                alt_fraction: 0.125,
                active_probability_proxy: 1.0,
            },
        ];
        let events = collect_pileup_fallback_events("chr1", 100, b"AAT", &active_loci, 20.0);
        assert!(events.is_empty());
    }

    #[test]
    fn rescue_collapsed_strong_snp_cluster_from_pileup_recovers_matching_candidates() {
        let valid_events = vec![
            test_variant("chr1", 100, b"A", b"G"),
            test_variant("chr1", 101, b"A", b"C"),
        ];
        let active_loci = vec![
            ReplayActiveLocusRow {
                contig: "chr1".to_string(),
                pos: 100,
                region: "chr1:100-102".to_string(),
                ref_base: b'A',
                depth: 16,
                snp_alt_count: 8,
                snp_best_alt: "G".to_string(),
                indel_alt_count: 0,
                indel_best_alt: String::new(),
                alt_fraction: 0.5,
                active_probability_proxy: 1.0,
            },
            ReplayActiveLocusRow {
                contig: "chr1".to_string(),
                pos: 101,
                region: "chr1:100-102".to_string(),
                ref_base: b'A',
                depth: 16,
                snp_alt_count: 8,
                snp_best_alt: "C".to_string(),
                indel_alt_count: 0,
                indel_best_alt: String::new(),
                alt_fraction: 0.5,
                active_probability_proxy: 1.0,
            },
        ];

        let rescued = rescue_collapsed_strong_snp_cluster_from_pileup(
            "chr1",
            100,
            b"AAT",
            &active_loci,
            &valid_events,
            20.0,
        );

        assert_eq!(rescued.len(), 2);
        assert!(rescued
            .iter()
            .any(|event| event.pos == 100 && event.ref_allele == b"A" && event.alt_allele == b"G"));
        assert!(rescued
            .iter()
            .any(|event| event.pos == 101 && event.ref_allele == b"A" && event.alt_allele == b"C"));
    }

    #[test]
    fn rescue_collapsed_strong_snp_cluster_from_pileup_does_not_recover_low_confidence_single_match(
    ) {
        let valid_events = vec![test_variant("chr1", 100, b"A", b"G")];
        let active_loci = vec![
            ReplayActiveLocusRow {
                contig: "chr1".to_string(),
                pos: 100,
                region: "chr1:100-102".to_string(),
                ref_base: b'A',
                depth: 16,
                snp_alt_count: 4,
                snp_best_alt: "G".to_string(),
                indel_alt_count: 0,
                indel_best_alt: String::new(),
                alt_fraction: 0.25,
                active_probability_proxy: 1.0,
            },
            ReplayActiveLocusRow {
                contig: "chr1".to_string(),
                pos: 101,
                region: "chr1:100-102".to_string(),
                ref_base: b'A',
                depth: 16,
                snp_alt_count: 8,
                snp_best_alt: "C".to_string(),
                indel_alt_count: 0,
                indel_best_alt: String::new(),
                alt_fraction: 0.5,
                active_probability_proxy: 1.0,
            },
        ];

        let rescued = rescue_collapsed_strong_snp_cluster_from_pileup(
            "chr1",
            100,
            b"AAT",
            &active_loci,
            &valid_events,
            20.0,
        );

        assert!(rescued.is_empty());
    }

    #[test]
    fn rescue_collapsed_strong_snp_cluster_from_pileup_recovers_single_exact_match_for_isolated_snp_locus(
    ) {
        let valid_events = vec![
            test_variant("chr1", 100, b"A", b"T"),
            test_variant("chr1", 99, b"G", b"GA"),
        ];
        let active_loci = vec![ReplayActiveLocusRow {
            contig: "chr1".to_string(),
            pos: 100,
            region: "chr1:99-101".to_string(),
            ref_base: b'A',
            depth: 16,
            snp_alt_count: 8,
            snp_best_alt: "T".to_string(),
            indel_alt_count: 0,
            indel_best_alt: String::new(),
            alt_fraction: 0.5,
            active_probability_proxy: 1.0,
        }];

        let rescued = rescue_collapsed_strong_snp_cluster_from_pileup(
            "chr1",
            99,
            b"GAA",
            &active_loci,
            &valid_events,
            20.0,
        );

        assert_eq!(rescued.len(), 1);
        assert_eq!(rescued[0].pos, 100);
        assert_eq!(rescued[0].ref_allele, b"A");
        assert_eq!(rescued[0].alt_allele, b"T");
    }

    #[test]
    fn rescue_collapsed_strong_snp_cluster_from_pileup_recovers_single_exact_match_with_only_weak_other_active_loci(
    ) {
        let valid_events = vec![
            test_variant("chr1", 100, b"A", b"T"),
            test_variant("chr1", 99, b"G", b"GA"),
        ];
        let active_loci = vec![
            ReplayActiveLocusRow {
                contig: "chr1".to_string(),
                pos: 100,
                region: "chr1:99-106".to_string(),
                ref_base: b'A',
                depth: 16,
                snp_alt_count: 8,
                snp_best_alt: "T".to_string(),
                indel_alt_count: 0,
                indel_best_alt: String::new(),
                alt_fraction: 0.5,
                active_probability_proxy: 1.0,
            },
            ReplayActiveLocusRow {
                contig: "chr1".to_string(),
                pos: 103,
                region: "chr1:99-106".to_string(),
                ref_base: b'C',
                depth: 6,
                snp_alt_count: 1,
                snp_best_alt: "A".to_string(),
                indel_alt_count: 0,
                indel_best_alt: String::new(),
                alt_fraction: 0.166667,
                active_probability_proxy: 1.0,
            },
            ReplayActiveLocusRow {
                contig: "chr1".to_string(),
                pos: 105,
                region: "chr1:99-106".to_string(),
                ref_base: b'G',
                depth: 6,
                snp_alt_count: 1,
                snp_best_alt: "C".to_string(),
                indel_alt_count: 0,
                indel_best_alt: String::new(),
                alt_fraction: 0.166667,
                active_probability_proxy: 1.0,
            },
        ];

        let rescued = rescue_collapsed_strong_snp_cluster_from_pileup(
            "chr1",
            99,
            b"GAAACCG",
            &active_loci,
            &valid_events,
            20.0,
        );

        assert_eq!(rescued.len(), 1);
        assert_eq!(rescued[0].pos, 100);
        assert_eq!(rescued[0].ref_allele, b"A");
        assert_eq!(rescued[0].alt_allele, b"T");
    }

    #[test]
    fn rescue_collapsed_strong_snp_cluster_from_pileup_does_not_recover_single_match_with_indel_evidence(
    ) {
        let valid_events = vec![
            test_variant("chr1", 100, b"A", b"T"),
            test_variant("chr1", 99, b"G", b"GA"),
        ];
        let active_loci = vec![ReplayActiveLocusRow {
            contig: "chr1".to_string(),
            pos: 100,
            region: "chr1:99-101".to_string(),
            ref_base: b'A',
            depth: 16,
            snp_alt_count: 8,
            snp_best_alt: "T".to_string(),
            indel_alt_count: 4,
            indel_best_alt: "INS:A".to_string(),
            alt_fraction: 0.5,
            active_probability_proxy: 1.0,
        }];

        let rescued = rescue_collapsed_strong_snp_cluster_from_pileup(
            "chr1",
            99,
            b"GAA",
            &active_loci,
            &valid_events,
            20.0,
        );

        assert!(rescued.is_empty());
    }

    #[test]
    fn merge_missing_strong_snp_cluster_rescues_from_pileup_adds_missing_exact_cluster_call() {
        let valid_events = vec![
            test_variant("chr1", 100, b"A", b"G"),
            test_variant("chr1", 101, b"A", b"C"),
            test_variant("chr1", 102, b"A", b"T"),
        ];
        let active_loci = vec![
            ReplayActiveLocusRow {
                contig: "chr1".to_string(),
                pos: 100,
                region: "chr1:100-103".to_string(),
                ref_base: b'A',
                depth: 16,
                snp_alt_count: 8,
                snp_best_alt: "G".to_string(),
                indel_alt_count: 0,
                indel_best_alt: String::new(),
                alt_fraction: 0.5,
                active_probability_proxy: 1.0,
            },
            ReplayActiveLocusRow {
                contig: "chr1".to_string(),
                pos: 101,
                region: "chr1:100-103".to_string(),
                ref_base: b'A',
                depth: 16,
                snp_alt_count: 7,
                snp_best_alt: "C".to_string(),
                indel_alt_count: 0,
                indel_best_alt: String::new(),
                alt_fraction: 0.4375,
                active_probability_proxy: 1.0,
            },
            ReplayActiveLocusRow {
                contig: "chr1".to_string(),
                pos: 102,
                region: "chr1:100-103".to_string(),
                ref_base: b'A',
                depth: 16,
                snp_alt_count: 6,
                snp_best_alt: "T".to_string(),
                indel_alt_count: 0,
                indel_best_alt: String::new(),
                alt_fraction: 0.375,
                active_probability_proxy: 1.0,
            },
        ];
        let mut final_calls = vec![
            test_variant("chr1", 100, b"A", b"G"),
            test_variant("chr1", 102, b"A", b"T"),
        ];

        merge_missing_strong_snp_cluster_rescues_from_pileup(
            &mut final_calls,
            "chr1",
            100,
            b"AAAT",
            &active_loci,
            &valid_events,
            20.0,
        );

        assert_eq!(final_calls.len(), 3);
        assert!(final_calls
            .iter()
            .any(|event| event.pos == 101 && event.ref_allele == b"A" && event.alt_allele == b"C"));
    }

    #[test]
    fn merge_missing_strong_snp_cluster_rescues_from_pileup_does_not_add_single_isolated_match() {
        let valid_events = vec![
            test_variant("chr1", 100, b"A", b"T"),
            test_variant("chr1", 100, b"A", b"AG"),
        ];
        let active_loci = vec![ReplayActiveLocusRow {
            contig: "chr1".to_string(),
            pos: 100,
            region: "chr1:100-101".to_string(),
            ref_base: b'A',
            depth: 16,
            snp_alt_count: 4,
            snp_best_alt: "T".to_string(),
            indel_alt_count: 0,
            indel_best_alt: String::new(),
            alt_fraction: 0.25,
            active_probability_proxy: 1.0,
        }];
        let mut final_calls = vec![test_variant("chr1", 100, b"A", b"AG")];

        merge_missing_strong_snp_cluster_rescues_from_pileup(
            &mut final_calls,
            "chr1",
            100,
            b"AA",
            &active_loci,
            &valid_events,
            20.0,
        );

        assert_eq!(final_calls.len(), 1);
        assert_eq!(final_calls[0].alt_allele, b"AG");
    }

    #[test]
    fn merge_missing_strong_snp_cluster_rescues_from_pileup_adds_single_high_confidence_match() {
        let valid_events = vec![
            test_variant("chr1", 100, b"A", b"T"),
            test_variant("chr1", 100, b"A", b"AG"),
        ];
        let active_loci = vec![ReplayActiveLocusRow {
            contig: "chr1".to_string(),
            pos: 100,
            region: "chr1:100-101".to_string(),
            ref_base: b'A',
            depth: 16,
            snp_alt_count: 12,
            snp_best_alt: "T".to_string(),
            indel_alt_count: 0,
            indel_best_alt: String::new(),
            alt_fraction: 0.75,
            active_probability_proxy: 1.0,
        }];
        let mut final_calls = vec![test_variant("chr1", 100, b"A", b"AG")];

        merge_missing_strong_snp_cluster_rescues_from_pileup(
            &mut final_calls,
            "chr1",
            100,
            b"AA",
            &active_loci,
            &valid_events,
            20.0,
        );

        assert_eq!(final_calls.len(), 2);
        assert!(final_calls
            .iter()
            .any(|event| event.pos == 100 && event.ref_allele == b"A" && event.alt_allele == b"T"));
    }

    #[test]
    fn prune_unsupported_simple_snp_calls_in_dense_clusters_drops_dense_non_active_snps() {
        let active_loci = vec![
            ReplayActiveLocusRow {
                contig: "chr1".to_string(),
                pos: 110,
                region: "chr1:100-140".to_string(),
                ref_base: b'A',
                depth: 20,
                snp_alt_count: 10,
                snp_best_alt: "C".to_string(),
                indel_alt_count: 0,
                indel_best_alt: String::new(),
                alt_fraction: 0.5,
                active_probability_proxy: 1.0,
            },
            ReplayActiveLocusRow {
                contig: "chr1".to_string(),
                pos: 136,
                region: "chr1:100-140".to_string(),
                ref_base: b'G',
                depth: 18,
                snp_alt_count: 9,
                snp_best_alt: "T".to_string(),
                indel_alt_count: 0,
                indel_best_alt: String::new(),
                alt_fraction: 0.5,
                active_probability_proxy: 1.0,
            },
        ];
        let mut final_calls = vec![
            test_variant("chr1", 100, b"A", b"G"),
            test_variant("chr1", 110, b"A", b"C"),
            test_variant("chr1", 118, b"T", b"TA"),
            test_variant("chr1", 120, b"T", b"G"),
            test_variant("chr1", 136, b"G", b"T"),
        ];

        prune_unsupported_simple_snp_calls_in_dense_clusters(&mut final_calls, &active_loci);

        assert_eq!(final_calls.len(), 3);
        assert!(final_calls
            .iter()
            .any(|call| call.pos == 110 && call.alt_allele == b"C"));
        assert!(final_calls
            .iter()
            .any(|call| call.pos == 118 && call.alt_allele == b"TA"));
        assert!(final_calls
            .iter()
            .any(|call| call.pos == 136 && call.alt_allele == b"T"));
        assert!(!final_calls
            .iter()
            .any(|call| call.pos == 100 && call.alt_allele == b"G"));
        assert!(!final_calls
            .iter()
            .any(|call| call.pos == 120 && call.alt_allele == b"G"));
    }

    #[test]
    fn prune_unsupported_simple_snp_calls_in_dense_clusters_keeps_isolated_unsupported_snp() {
        let active_loci = vec![ReplayActiveLocusRow {
            contig: "chr1".to_string(),
            pos: 110,
            region: "chr1:100-160".to_string(),
            ref_base: b'A',
            depth: 20,
            snp_alt_count: 10,
            snp_best_alt: "C".to_string(),
            indel_alt_count: 0,
            indel_best_alt: String::new(),
            alt_fraction: 0.5,
            active_probability_proxy: 1.0,
        }];
        let mut final_calls = vec![
            test_variant("chr1", 110, b"A", b"C"),
            test_variant("chr1", 140, b"T", b"G"),
        ];

        prune_unsupported_simple_snp_calls_in_dense_clusters(&mut final_calls, &active_loci);

        assert_eq!(final_calls.len(), 2);
        assert!(final_calls
            .iter()
            .any(|call| call.pos == 110 && call.alt_allele == b"C"));
        assert!(final_calls
            .iter()
            .any(|call| call.pos == 140 && call.alt_allele == b"G"));
    }

    #[test]
    fn supplement_missing_pileup_events_adds_missing_event_and_haplotype() {
        let mut valid_events = vec![test_variant("chr1", 100, b"A", b"G")];
        let mut local_haplotypes =
            haplotypes_from_candidate_events("chr1", 100, b"AAT", &valid_events);
        let active_loci = vec![
            ReplayActiveLocusRow {
                contig: "chr1".to_string(),
                pos: 100,
                region: "chr1:100-102".to_string(),
                ref_base: b'A',
                depth: 16,
                snp_alt_count: 8,
                snp_best_alt: "G".to_string(),
                indel_alt_count: 0,
                indel_best_alt: String::new(),
                alt_fraction: 0.5,
                active_probability_proxy: 1.0,
            },
            ReplayActiveLocusRow {
                contig: "chr1".to_string(),
                pos: 101,
                region: "chr1:100-102".to_string(),
                ref_base: b'A',
                depth: 16,
                snp_alt_count: 8,
                snp_best_alt: "C".to_string(),
                indel_alt_count: 0,
                indel_best_alt: String::new(),
                alt_fraction: 0.5,
                active_probability_proxy: 1.0,
            },
        ];

        supplement_missing_pileup_events(
            "chr1",
            100,
            b"AAT",
            &active_loci,
            20.0,
            &mut local_haplotypes,
            &mut valid_events,
        );

        assert_eq!(valid_events.len(), 2);
        assert_eq!(valid_events[1].pos, 101);
        assert_eq!(valid_events[1].ref_allele, b"A");
        assert_eq!(valid_events[1].alt_allele, b"C");
        assert!(local_haplotypes
            .iter()
            .any(|hap| hap.event_indices == vec![1] && hap.bases == b"ACT"));
    }

    #[test]
    fn supplement_missing_pileup_events_seeds_low_qual_zero_candidate_simple_snp() {
        let mut valid_events = Vec::new();
        let mut local_haplotypes = haplotypes_from_candidate_events("chr1", 100, b"A", &[]);
        let active_loci = vec![ReplayActiveLocusRow {
            contig: "chr1".to_string(),
            pos: 100,
            region: "chr1:100-100".to_string(),
            ref_base: b'A',
            depth: 24,
            snp_alt_count: 3,
            snp_best_alt: "G".to_string(),
            indel_alt_count: 0,
            indel_best_alt: String::new(),
            alt_fraction: 0.125,
            active_probability_proxy: 1.0,
        }];

        supplement_missing_pileup_events(
            "chr1",
            100,
            b"A",
            &active_loci,
            20.0,
            &mut local_haplotypes,
            &mut valid_events,
        );

        assert_eq!(valid_events.len(), 1);
        assert_eq!(valid_events[0].pos, 100);
        assert_eq!(valid_events[0].ref_allele, b"A");
        assert_eq!(valid_events[0].alt_allele, b"G");
        assert!(valid_events[0].qual < 20);
        assert!(local_haplotypes
            .iter()
            .any(|hap| !hap.is_ref && hap.event_indices == vec![0] && hap.bases == b"G"));
    }

    #[test]
    fn supplement_missing_pileup_events_skips_low_qual_pileup_het_seed() {
        let mut valid_events = Vec::new();
        let mut local_haplotypes = haplotypes_from_candidate_events("chr1", 100, b"G", &[]);
        let active_loci = vec![ReplayActiveLocusRow {
            contig: "chr1".to_string(),
            pos: 100,
            region: "chr1:100-100".to_string(),
            ref_base: b'G',
            depth: 17,
            snp_alt_count: 3,
            snp_best_alt: "C".to_string(),
            indel_alt_count: 0,
            indel_best_alt: String::new(),
            alt_fraction: 0.176471,
            active_probability_proxy: 1.0,
        }];

        supplement_missing_pileup_events(
            "chr1",
            100,
            b"G",
            &active_loci,
            20.0,
            &mut local_haplotypes,
            &mut valid_events,
        );

        assert!(valid_events.is_empty());
        assert!(local_haplotypes.iter().all(|hap| hap.is_ref));
    }

    #[test]
    fn supplement_missing_pileup_events_overlays_missing_snp_on_existing_alt_haplotype() {
        let mut valid_events = vec![test_variant("chr1", 100, b"A", b"G")];
        let mut local_haplotypes =
            haplotypes_from_candidate_events("chr1", 100, b"AAT", &valid_events);
        let active_loci = vec![
            ReplayActiveLocusRow {
                contig: "chr1".to_string(),
                pos: 100,
                region: "chr1:100-102".to_string(),
                ref_base: b'A',
                depth: 16,
                snp_alt_count: 8,
                snp_best_alt: "G".to_string(),
                indel_alt_count: 0,
                indel_best_alt: String::new(),
                alt_fraction: 0.5,
                active_probability_proxy: 1.0,
            },
            ReplayActiveLocusRow {
                contig: "chr1".to_string(),
                pos: 101,
                region: "chr1:100-102".to_string(),
                ref_base: b'A',
                depth: 16,
                snp_alt_count: 8,
                snp_best_alt: "C".to_string(),
                indel_alt_count: 0,
                indel_best_alt: String::new(),
                alt_fraction: 0.5,
                active_probability_proxy: 1.0,
            },
        ];

        supplement_missing_pileup_events(
            "chr1",
            100,
            b"AAT",
            &active_loci,
            20.0,
            &mut local_haplotypes,
            &mut valid_events,
        );

        assert_eq!(valid_events.len(), 2);
        assert!(local_haplotypes
            .iter()
            .any(|hap| hap.event_indices == vec![0, 1] && hap.bases == b"GCT"));
    }

    #[test]
    fn supplement_missing_pileup_events_does_not_overlay_missing_snp_on_indel_haplotype() {
        let mut valid_events = vec![test_variant("chr1", 100, b"T", b"TA")];
        let mut local_haplotypes =
            haplotypes_from_candidate_events("chr1", 100, b"TT", &valid_events);
        let active_loci = vec![
            ReplayActiveLocusRow {
                contig: "chr1".to_string(),
                pos: 100,
                region: "chr1:100-101".to_string(),
                ref_base: b'T',
                depth: 16,
                snp_alt_count: 0,
                snp_best_alt: String::new(),
                indel_alt_count: 8,
                indel_best_alt: "INS:A".to_string(),
                alt_fraction: 0.5,
                active_probability_proxy: 1.0,
            },
            ReplayActiveLocusRow {
                contig: "chr1".to_string(),
                pos: 100,
                region: "chr1:100-101".to_string(),
                ref_base: b'T',
                depth: 16,
                snp_alt_count: 8,
                snp_best_alt: "C".to_string(),
                indel_alt_count: 8,
                indel_best_alt: "INS:A".to_string(),
                alt_fraction: 0.5,
                active_probability_proxy: 1.0,
            },
        ];

        supplement_missing_pileup_events(
            "chr1",
            100,
            b"TT",
            &active_loci,
            20.0,
            &mut local_haplotypes,
            &mut valid_events,
        );

        assert_eq!(valid_events.len(), 2);
        assert!(local_haplotypes
            .iter()
            .any(|hap| hap.event_indices == vec![0]));
        assert!(local_haplotypes
            .iter()
            .any(|hap| hap.event_indices == vec![1]));
        assert!(!local_haplotypes
            .iter()
            .any(|hap| hap.event_indices == vec![0, 1]));
    }

    #[test]
    fn supplement_missing_pileup_events_skips_indels_when_region_already_has_events() {
        let mut valid_events = vec![test_variant("chr1", 100, b"A", b"G")];
        let mut local_haplotypes =
            haplotypes_from_candidate_events("chr1", 100, b"AAAT", &valid_events);
        let active_loci = vec![ReplayActiveLocusRow {
            contig: "chr1".to_string(),
            pos: 101,
            region: "chr1:100-103".to_string(),
            ref_base: b'A',
            depth: 16,
            snp_alt_count: 0,
            snp_best_alt: String::new(),
            indel_alt_count: 8,
            indel_best_alt: "INS:A".to_string(),
            alt_fraction: 0.5,
            active_probability_proxy: 1.0,
        }];

        supplement_missing_pileup_events(
            "chr1",
            100,
            b"AAAT",
            &active_loci,
            20.0,
            &mut local_haplotypes,
            &mut valid_events,
        );

        assert_eq!(valid_events.len(), 1);
        assert_eq!(valid_events[0].ref_allele, b"A");
        assert_eq!(valid_events[0].alt_allele, b"G");
    }

    #[test]
    fn supplement_missing_pileup_events_skips_weak_snp_in_dense_existing_snp_cluster() {
        let local_ref_bases = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        let mut valid_events = vec![
            test_variant("chr1", 100, b"A", b"G"),
            test_variant("chr1", 121, b"A", b"C"),
        ];
        let mut local_haplotypes =
            haplotypes_from_candidate_events("chr1", 100, local_ref_bases, &valid_events);
        let active_loci = vec![
            ReplayActiveLocusRow {
                contig: "chr1".to_string(),
                pos: 100,
                region: "chr1:100-129".to_string(),
                ref_base: b'A',
                depth: 16,
                snp_alt_count: 8,
                snp_best_alt: "G".to_string(),
                indel_alt_count: 0,
                indel_best_alt: String::new(),
                alt_fraction: 0.5,
                active_probability_proxy: 1.0,
            },
            ReplayActiveLocusRow {
                contig: "chr1".to_string(),
                pos: 121,
                region: "chr1:100-129".to_string(),
                ref_base: b'A',
                depth: 16,
                snp_alt_count: 8,
                snp_best_alt: "C".to_string(),
                indel_alt_count: 0,
                indel_best_alt: String::new(),
                alt_fraction: 0.5,
                active_probability_proxy: 1.0,
            },
            ReplayActiveLocusRow {
                contig: "chr1".to_string(),
                pos: 127,
                region: "chr1:100-129".to_string(),
                ref_base: b'A',
                depth: 14,
                snp_alt_count: 4,
                snp_best_alt: "T".to_string(),
                indel_alt_count: 0,
                indel_best_alt: String::new(),
                alt_fraction: 0.285714,
                active_probability_proxy: 1.0,
            },
        ];

        supplement_missing_pileup_events(
            "chr1",
            100,
            local_ref_bases,
            &active_loci,
            20.0,
            &mut local_haplotypes,
            &mut valid_events,
        );

        assert_eq!(valid_events.len(), 2);
        assert!(!valid_events.iter().any(|event| event.pos == 127));
    }

    #[test]
    fn supplement_missing_pileup_events_keeps_stronger_snp_in_dense_existing_snp_cluster() {
        let local_ref_bases = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        let mut valid_events = vec![
            test_variant("chr1", 100, b"A", b"G"),
            test_variant("chr1", 121, b"A", b"C"),
        ];
        let mut local_haplotypes =
            haplotypes_from_candidate_events("chr1", 100, local_ref_bases, &valid_events);
        let active_loci = vec![
            ReplayActiveLocusRow {
                contig: "chr1".to_string(),
                pos: 100,
                region: "chr1:100-129".to_string(),
                ref_base: b'A',
                depth: 16,
                snp_alt_count: 8,
                snp_best_alt: "G".to_string(),
                indel_alt_count: 0,
                indel_best_alt: String::new(),
                alt_fraction: 0.5,
                active_probability_proxy: 1.0,
            },
            ReplayActiveLocusRow {
                contig: "chr1".to_string(),
                pos: 121,
                region: "chr1:100-129".to_string(),
                ref_base: b'A',
                depth: 16,
                snp_alt_count: 8,
                snp_best_alt: "C".to_string(),
                indel_alt_count: 0,
                indel_best_alt: String::new(),
                alt_fraction: 0.5,
                active_probability_proxy: 1.0,
            },
            ReplayActiveLocusRow {
                contig: "chr1".to_string(),
                pos: 127,
                region: "chr1:100-129".to_string(),
                ref_base: b'A',
                depth: 16,
                snp_alt_count: 5,
                snp_best_alt: "T".to_string(),
                indel_alt_count: 0,
                indel_best_alt: String::new(),
                alt_fraction: 0.3125,
                active_probability_proxy: 1.0,
            },
        ];

        supplement_missing_pileup_events(
            "chr1",
            100,
            local_ref_bases,
            &active_loci,
            20.0,
            &mut local_haplotypes,
            &mut valid_events,
        );

        assert_eq!(valid_events.len(), 3);
        assert!(valid_events.iter().any(|event| event.pos == 127));
    }

    #[test]
    fn assemble_haplotypes_skips_non_acgt_haplotype_sequences() {
        let (haplotypes, events) = assemble_haplotypes(
            "chr1",
            100,
            b"AAACC",
            &[
                b"AAGCC".to_vec(),
                b"AAGCC".to_vec(),
                b"AANCC".to_vec(),
                b"AANCC".to_vec(),
            ],
            &[2],
            0,
        );

        assert!(events
            .iter()
            .any(|event| event.pos == 102 && event.ref_allele == b"A" && event.alt_allele == b"G"));
        assert!(haplotypes
            .iter()
            .any(|hap| !hap.is_ref && hap.bases == b"AAGCC"));
        assert!(haplotypes.iter().all(|hap| is_regular_bases(&hap.bases)));
    }

    #[test]
    fn align_haplotype_to_reference_rejects_softclipped_path_alignments() {
        let align_result = align_haplotype_to_reference(b"ACGTACGT", b"ACGTACGTGGGG");
        assert!(align_result.is_none());
    }

    #[test]
    fn align_haplotype_to_reference_keeps_simple_indel_alignments() {
        let align_result = align_haplotype_to_reference(b"ACGTACGT", b"ACGTTACGT").unwrap();
        assert_eq!(align_result.alignment_offset, 0);
        assert!(align_result
            .cigar
            .iter()
            .all(|ce| !matches!(ce, rust_htslib::bam::record::Cigar::SoftClip(_))));
        assert_eq!(
            align_result
                .cigar
                .iter()
                .filter(|ce| matches!(ce, rust_htslib::bam::record::Cigar::Ins(_)))
                .count(),
            1
        );
    }
}

fn is_regular_allele(allele: &[u8]) -> bool {
    allele
        .iter()
        .all(|&b| matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't'))
}

fn is_regular_bases(bases: &[u8]) -> bool {
    bases.iter().all(|&b| is_acgt(normalize_base(b)))
}

fn haplotype_to_reference_sw_parameters() -> crate::smith_waterman::SWParameters {
    crate::smith_waterman::SWParameters {
        match_value: 200,
        mismatch_penalty: -150,
        gap_open_penalty: -260,
        gap_extend_penalty: -11,
    }
}

fn align_haplotype_to_reference(
    ref_bases: &[u8],
    haplotype_bases: &[u8],
) -> Option<crate::smith_waterman::SWAlignmentResult> {
    let align_result = crate::smith_waterman::align(
        ref_bases,
        haplotype_bases,
        &haplotype_to_reference_sw_parameters(),
        crate::smith_waterman::SWOverhangStrategy::SoftClip,
    );

    if align_result.alignment_offset > 0
        || align_result
            .cigar
            .iter()
            .any(|ce| matches!(ce, rust_htslib::bam::record::Cigar::SoftClip(_)))
    {
        return None;
    }

    Some(align_result)
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
                if element_length <= 10
                    && ref_pos > 0
                    && cigar_index > 0
                    && cigar_index < num_cigar_elements - 1
                {
                    let insertion_start = region_start + ref_pos as u64 - 1;
                    let ref_byte = ref_bases[ref_pos as usize - 1];
                    let mut insertion_bases = vec![ref_byte];
                    insertion_bases.extend_from_slice(
                        &alt_bases[alignment_pos..alignment_pos + element_length],
                    );

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
                    deletion_bases.extend_from_slice(
                        &ref_bases[ref_pos as usize..ref_pos as usize + element_length],
                    );

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

                    let ref_allele =
                        ref_bases[ref_pos as usize + start..=ref_pos as usize + end].to_vec();
                    let alt_allele =
                        alt_bases[alignment_pos + start..=alignment_pos + end].to_vec();

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

        graph.prune(2); // min_prune_factor = 2

        let best_paths = graph.find_best_haplotypes(source_idx, sink_idx, 128);

        let mut found_nonref = false;
        for path in best_paths {
            let seq = graph.reconstruct_sequence(&path);
            if !is_regular_bases(&seq) {
                continue;
            }
            if !assembled_haplotypes_set.contains(&seq) {
                assembled_haplotypes_set.insert(seq.clone());

                let Some(align_result) = align_haplotype_to_reference(ref_bases, &seq) else {
                    continue;
                };

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
            if !unique_events.iter().any(|e| {
                e.pos == event.pos
                    && e.ref_allele == event.ref_allele
                    && e.alt_allele == event.alt_allele
            }) {
                unique_events.push(event.clone());
            }
        }
    }

    // Now map event indices back to each haplotype
    for (hap_idx, events) in hap_events.iter().enumerate() {
        for event in events {
            if let Some(idx) = unique_events.iter().position(|e| {
                e.pos == event.pos
                    && e.ref_allele == event.ref_allele
                    && e.alt_allele == event.alt_allele
            }) {
                if !local_haps[hap_idx].event_indices.contains(&idx) {
                    local_haps[hap_idx].event_indices.push(idx);
                }
            }
        }
    }

    (local_haps, unique_events)
}
