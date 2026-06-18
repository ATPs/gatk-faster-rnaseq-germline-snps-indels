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

include!("activity.rs");
include!("evidence.rs");
include!("genotyping.rs");
include!("candidate_events.rs");
include!("replay_io.rs");
include!("tests_block.rs");
include!("alignment.rs");
