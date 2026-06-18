use anyhow::{bail, Context, Result};
use flate2::read::MultiGzDecoder;
use serde_json::Value;
use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Read, Write};
use std::path::{Path, PathBuf};

#[derive(Debug, Clone, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct VcfKey {
    pub chrom: String,
    pub pos: u64,
    pub ref_allele: String,
    pub alt: String,
}

impl VcfKey {
    pub fn display(&self) -> String {
        format!(
            "{}:{}:{}>{}",
            self.chrom, self.pos, self.ref_allele, self.alt
        )
    }
}

#[derive(Debug, Clone)]
pub struct VcfRecord {
    pub key: VcfKey,
    pub filter: String,
    pub gt: String,
    pub qual: String,
    pub info: BTreeMap<String, String>,
}

impl VcfRecord {
    pub fn is_pass(&self) -> bool {
        self.filter == "PASS"
    }

    pub fn variant_type(&self) -> &'static str {
        variant_type(&self.key.ref_allele, &self.key.alt)
    }
}

#[derive(Debug, Clone, Default)]
pub struct VcfDatasetSummary {
    pub total: usize,
    pub pass: usize,
    pub nonpass: usize,
    pub filters: BTreeMap<String, usize>,
    pub all_types: BTreeMap<String, usize>,
    pub pass_types: BTreeMap<String, usize>,
}

#[derive(Debug, Clone, Default)]
pub struct VcfSetComparison {
    pub a_count: usize,
    pub b_count: usize,
    pub shared: usize,
    pub a_private: usize,
    pub b_private: usize,
    pub a_sensitivity: f64,
    pub b_precision_vs_a: f64,
    pub shared_types: BTreeMap<String, usize>,
    pub a_private_types: BTreeMap<String, usize>,
    pub b_private_types: BTreeMap<String, usize>,
    pub gt_same: usize,
    pub gt_diff: usize,
}

#[derive(Debug, Clone)]
pub struct VcfComparison {
    pub a_label: String,
    pub b_label: String,
    pub a_summary: VcfDatasetSummary,
    pub b_summary: VcfDatasetSummary,
    pub all_records: VcfSetComparison,
    pub pass_records: VcfSetComparison,
    pub detail_rows: Vec<VcfDetailRow>,
}

#[derive(Debug, Clone)]
pub struct VcfDetailRow {
    pub category: String,
    pub key: VcfKey,
    pub variant_type: &'static str,
    pub a_filter: String,
    pub b_filter: String,
    pub a_gt: String,
    pub b_gt: String,
}

#[derive(Debug, Clone)]
pub struct SelectedRegion {
    pub category: String,
    pub key: VcfKey,
    pub variant_type: &'static str,
    pub start: u64,
    pub end: u64,
    pub a_gt: String,
    pub b_gt: String,
}

#[derive(Debug, Clone)]
pub struct StageDiffConfig {
    pub java_path: PathBuf,
    pub rust_path: PathBuf,
    pub key_columns: Vec<String>,
    pub numeric_tolerance: f64,
    pub output_prefix: PathBuf,
    pub stage_name: String,
}

#[derive(Debug, Clone, Default)]
pub struct StageDiffSummary {
    pub java_rows: usize,
    pub rust_rows: usize,
    pub shared_rows: usize,
    pub java_private_rows: usize,
    pub rust_private_rows: usize,
    pub shared_rows_with_diffs: usize,
    pub field_diffs: usize,
}

#[derive(Debug, Clone)]
pub struct GatkDebugExtractConfig {
    pub genotyper_debug: Option<PathBuf>,
    pub assembly_state: Option<PathBuf>,
    pub output_prefix: PathBuf,
}

#[derive(Debug, Clone, Default)]
pub struct GatkDebugExtractStats {
    pub genotyper_haplotypes: usize,
    pub pairhmm_scores: usize,
    pub event_allele_links: usize,
    pub allele_likelihoods: usize,
    pub read_quality_rows: usize,
    pub assembly_reads: usize,
    pub assembly_haplotypes: usize,
}

#[derive(Debug, Clone)]
struct StageRow {
    fields: BTreeMap<String, String>,
}

#[derive(Debug, Clone)]
struct DebugReadMeta {
    index: usize,
    name: String,
    cigar: String,
    mapq: String,
    loc: String,
    unclipped_loc: String,
    length: String,
}

#[derive(Debug, Clone)]
struct DebugHaplotypeMeta {
    start: String,
    end: String,
    kmer: String,
    length: String,
    cigar: String,
    is_ref: bool,
}

include!("vcf_compare.rs");
include!("debug_tools.rs");
include!("vcf_helpers.rs");
include!("stage_diff.rs");
include!("tests_block.rs");
