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

pub fn read_vcf_records(path: &Path) -> Result<Vec<VcfRecord>> {
    let mut records = Vec::new();
    let mut reader = text_reader(path)?;
    let mut line = String::new();
    while reader.read_line(&mut line)? != 0 {
        if line.starts_with('#') || line.trim().is_empty() {
            line.clear();
            continue;
        }
        records.extend(
            parse_vcf_records(line.trim_end())
                .with_context(|| format!("parsing VCF record from {}", path.display()))?,
        );
        line.clear();
    }
    Ok(records)
}

pub fn compare_vcfs(
    a_path: &Path,
    b_path: &Path,
    a_label: &str,
    b_label: &str,
) -> Result<VcfComparison> {
    let a_records =
        read_vcf_records(a_path).with_context(|| format!("reading {}", a_path.display()))?;
    let b_records =
        read_vcf_records(b_path).with_context(|| format!("reading {}", b_path.display()))?;
    let a_map = record_map(&a_records);
    let b_map = record_map(&b_records);
    let detail_rows = build_detail_rows(&a_map, &b_map);
    Ok(VcfComparison {
        a_label: a_label.to_string(),
        b_label: b_label.to_string(),
        a_summary: summarize_vcf_records(&a_records),
        b_summary: summarize_vcf_records(&b_records),
        all_records: compare_record_maps(&a_map, &b_map, false),
        pass_records: compare_record_maps(&a_map, &b_map, true),
        detail_rows,
    })
}

pub fn write_vcf_comparison(prefix: &Path, comparison: &VcfComparison) -> Result<()> {
    create_output_parent(prefix)?;
    write_vcf_summary_markdown(&prefixed_path(prefix, "summary.md"), comparison)?;
    write_vcf_summary_tsv(&prefixed_path(prefix, "summary.tsv"), comparison)?;
    write_vcf_summary_json(&prefixed_path(prefix, "summary.json"), comparison)?;
    write_vcf_detail_tsv(&prefixed_path(prefix, "details.tsv"), comparison)?;
    Ok(())
}

pub fn select_regions(
    a_path: &Path,
    b_path: &Path,
    padding: u64,
    max_per_category: usize,
    pass_only: bool,
) -> Result<Vec<SelectedRegion>> {
    let a_records =
        read_vcf_records(a_path).with_context(|| format!("reading {}", a_path.display()))?;
    let b_records =
        read_vcf_records(b_path).with_context(|| format!("reading {}", b_path.display()))?;
    let a_map = record_map(&a_records);
    let b_map = record_map(&b_records);
    let a_keys = filtered_keys(&a_map, pass_only);
    let b_keys = filtered_keys(&b_map, pass_only);
    let mut rows = Vec::new();

    push_key_category(
        &mut rows,
        "a_private_snp",
        a_keys
            .difference(&b_keys)
            .filter(|key| variant_type(&key.ref_allele, &key.alt) == "SNP"),
        &a_map,
        &b_map,
        padding,
        max_per_category,
    );
    push_key_category(
        &mut rows,
        "a_private_indel",
        a_keys
            .difference(&b_keys)
            .filter(|key| variant_type(&key.ref_allele, &key.alt) != "SNP"),
        &a_map,
        &b_map,
        padding,
        max_per_category,
    );
    push_key_category(
        &mut rows,
        "b_private_snp",
        b_keys
            .difference(&a_keys)
            .filter(|key| variant_type(&key.ref_allele, &key.alt) == "SNP"),
        &a_map,
        &b_map,
        padding,
        max_per_category,
    );
    push_key_category(
        &mut rows,
        "b_private_indel",
        b_keys
            .difference(&a_keys)
            .filter(|key| variant_type(&key.ref_allele, &key.alt) != "SNP"),
        &a_map,
        &b_map,
        padding,
        max_per_category,
    );
    push_key_category(
        &mut rows,
        "shared_gt_diff",
        a_keys.intersection(&b_keys).filter(|key| {
            let a_gt = a_map.get(*key).map(|r| r.gt.as_str()).unwrap_or("");
            let b_gt = b_map.get(*key).map(|r| r.gt.as_str()).unwrap_or("");
            !a_gt.is_empty() && !b_gt.is_empty() && a_gt != b_gt
        }),
        &a_map,
        &b_map,
        padding,
        max_per_category,
    );
    push_key_category(
        &mut rows,
        "shared_indel",
        a_keys
            .intersection(&b_keys)
            .filter(|key| variant_type(&key.ref_allele, &key.alt) != "SNP"),
        &a_map,
        &b_map,
        padding,
        max_per_category,
    );
    push_key_category(
        &mut rows,
        "shared_match",
        a_keys.intersection(&b_keys).filter(|key| {
            let a_gt = a_map.get(*key).map(|r| r.gt.as_str()).unwrap_or("");
            let b_gt = b_map.get(*key).map(|r| r.gt.as_str()).unwrap_or("");
            a_gt == b_gt
        }),
        &a_map,
        &b_map,
        padding,
        max_per_category,
    );

    rows.sort_by(|a, b| {
        a.key
            .chrom
            .cmp(&b.key.chrom)
            .then_with(|| a.start.cmp(&b.start))
            .then_with(|| a.end.cmp(&b.end))
            .then_with(|| a.category.cmp(&b.category))
    });
    Ok(rows)
}

pub fn write_selected_regions(
    prefix: &Path,
    rows: &[SelectedRegion],
    interval_list_template: Option<&Path>,
) -> Result<()> {
    create_output_parent(prefix)?;
    write_region_manifest(&prefixed_path(prefix, "manifest.tsv"), rows)?;
    write_region_bed(&prefixed_path(prefix, "regions.bed"), rows)?;
    if let Some(template) = interval_list_template {
        write_region_interval_list(
            &prefixed_path(prefix, "regions.interval_list"),
            rows,
            template,
        )?;
    }
    Ok(())
}

pub fn run_stage_diff(config: &StageDiffConfig) -> Result<StageDiffSummary> {
    let java_rows = read_stage_rows(&config.java_path)
        .with_context(|| format!("reading Java stage file {}", config.java_path.display()))?;
    let rust_rows = read_stage_rows(&config.rust_path)
        .with_context(|| format!("reading Rust stage file {}", config.rust_path.display()))?;
    let java_map = stage_row_map(&java_rows, &config.key_columns, "java")?;
    let rust_map = stage_row_map(&rust_rows, &config.key_columns, "rust")?;
    let mut summary = StageDiffSummary {
        java_rows: java_map.len(),
        rust_rows: rust_map.len(),
        ..StageDiffSummary::default()
    };
    let mut detail = Vec::new();
    let keys: BTreeSet<_> = java_map.keys().chain(rust_map.keys()).cloned().collect();
    for key in keys {
        match (java_map.get(&key), rust_map.get(&key)) {
            (Some(java), Some(rust)) => {
                summary.shared_rows += 1;
                let diffs =
                    diff_stage_fields(java, rust, &config.key_columns, config.numeric_tolerance);
                if !diffs.is_empty() {
                    summary.shared_rows_with_diffs += 1;
                    summary.field_diffs += diffs.len();
                    for (field, java_value, rust_value) in diffs {
                        detail.push(vec![
                            "field_diff".to_string(),
                            key.clone(),
                            field,
                            java_value,
                            rust_value,
                        ]);
                    }
                }
            }
            (Some(_), None) => {
                summary.java_private_rows += 1;
                detail.push(vec![
                    "java_private".to_string(),
                    key.clone(),
                    String::new(),
                    String::new(),
                    String::new(),
                ]);
            }
            (None, Some(_)) => {
                summary.rust_private_rows += 1;
                detail.push(vec![
                    "rust_private".to_string(),
                    key.clone(),
                    String::new(),
                    String::new(),
                    String::new(),
                ]);
            }
            (None, None) => {}
        }
    }

    create_output_parent(&config.output_prefix)?;
    write_stage_summary_markdown(
        &prefixed_path(&config.output_prefix, "summary.md"),
        config,
        &summary,
    )?;
    write_stage_summary_tsv(
        &prefixed_path(&config.output_prefix, "summary.tsv"),
        config,
        &summary,
    )?;
    write_stage_detail_tsv(
        &prefixed_path(&config.output_prefix, "details.tsv"),
        &detail,
    )?;
    Ok(summary)
}

pub fn write_acceptance_report(output: &Path, inputs: &[PathBuf], title: &str) -> Result<()> {
    if let Some(parent) = output.parent() {
        if !parent.as_os_str().is_empty() {
            fs::create_dir_all(parent).with_context(|| format!("creating {}", parent.display()))?;
        }
    }
    let mut out = File::create(output).with_context(|| format!("creating {}", output.display()))?;
    writeln!(out, "# {title}")?;
    writeln!(out)?;
    writeln!(out, "Generated from {} input report(s).", inputs.len())?;
    for path in inputs {
        writeln!(out)?;
        writeln!(out, "## {}", path.display())?;
        writeln!(out)?;
        let text =
            fs::read_to_string(path).with_context(|| format!("reading {}", path.display()))?;
        out.write_all(text.as_bytes())?;
        if !text.ends_with('\n') {
            writeln!(out)?;
        }
    }
    Ok(())
}

pub fn write_vcf_genotype_table(vcf: &Path, output: &Path) -> Result<usize> {
    if let Some(parent) = output.parent() {
        if !parent.as_os_str().is_empty() {
            fs::create_dir_all(parent).with_context(|| format!("creating {}", parent.display()))?;
        }
    }
    let mut out = File::create(output).with_context(|| format!("creating {}", output.display()))?;
    writeln!(
        out,
        "chrom\tpos\tref\talt\tqual\tfilter\tgt\tgq\tdp\tad_ref\tad_alt\tfs\tqd\tpl\tdb"
    )?;
    let mut rows = 0usize;
    let mut reader = text_reader(vcf)?;
    let mut line = String::new();
    while reader.read_line(&mut line)? != 0 {
        if line.starts_with('#') || line.trim().is_empty() {
            line.clear();
            continue;
        }
        let fields: Vec<&str> = line.trim_end().split('\t').collect();
        if fields.len() < 8 {
            bail!("VCF record has fewer than 8 columns in {}", vcf.display());
        }
        let info = parse_info(fields[7]);
        let sample = if fields.len() >= 10 {
            parse_sample_value_map(fields[8], fields[9])
        } else {
            BTreeMap::new()
        };
        let ad = sample.get("AD").map(|value| value.as_str()).unwrap_or("");
        let (ad_ref, ad_alt) = split_ad(ad);
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            fields[0],
            fields[1],
            fields[3],
            fields[4],
            fields[5],
            fields[6],
            sample.get("GT").map(|value| value.as_str()).unwrap_or(""),
            sample.get("GQ").map(|value| value.as_str()).unwrap_or(""),
            sample
                .get("DP")
                .map(|value| value.as_str())
                .or_else(|| info.get("DP").map(|value| value.as_str()))
                .unwrap_or(""),
            ad_ref,
            ad_alt,
            info.get("FS").map(|value| value.as_str()).unwrap_or(""),
            info.get("QD").map(|value| value.as_str()).unwrap_or(""),
            sample.get("PL").map(|value| value.as_str()).unwrap_or(""),
            info.contains_key("DB"),
        )?;
        rows += 1;
        line.clear();
    }
    Ok(rows)
}

pub fn extract_gatk_debug_tables(config: &GatkDebugExtractConfig) -> Result<GatkDebugExtractStats> {
    create_output_parent(&config.output_prefix)?;
    let mut stats = GatkDebugExtractStats::default();
    if let Some(path) = &config.genotyper_debug {
        extract_genotyper_debug(path, &config.output_prefix, &mut stats)
            .with_context(|| format!("extracting {}", path.display()))?;
    }
    if let Some(path) = &config.assembly_state {
        extract_assembly_state(path, &config.output_prefix, &mut stats)
            .with_context(|| format!("extracting {}", path.display()))?;
    }
    Ok(stats)
}

pub fn prefixed_path(prefix: &Path, suffix: &str) -> PathBuf {
    PathBuf::from(format!("{}.{}", prefix.display(), suffix))
}

fn extract_genotyper_debug(
    path: &Path,
    prefix: &Path,
    stats: &mut GatkDebugExtractStats,
) -> Result<()> {
    let mut haplotypes = File::create(prefixed_path(prefix, "genotyper_haplotypes.tsv"))?;
    let mut pairhmm = File::create(prefixed_path(prefix, "pairhmm.tsv"))?;
    let mut read_qualities = File::create(prefixed_path(prefix, "read_qualities.tsv"))?;
    let mut events = File::create(prefixed_path(prefix, "events.tsv"))?;
    let mut event_haps = File::create(prefixed_path(prefix, "event_allele_haps.tsv"))?;
    let mut allele_likelihoods = File::create(prefixed_path(prefix, "allele_likelihoods.tsv"))?;

    writeln!(
        haplotypes,
        "region\tstage\thaplotype\tspan_start\tspan_end\tkmer\tlength\tcigar\tis_ref\tbases"
    )?;
    writeln!(
        pairhmm,
        "region\tread\thaplotype\tread_index\tcigar\tmapq\tloc\tunclipped_loc\tlength\tscore"
    )?;
    writeln!(
        read_qualities,
        "region\tread\tread_index\tcigar\tmapq\tloc\tunclipped_loc\tlength\tqualities"
    )?;
    writeln!(events, "region\tevent\tchrom\tpos\ttype\talleles\traw")?;
    writeln!(event_haps, "region\tevent\tallele\thaplotypes")?;
    writeln!(
        allele_likelihoods,
        "region\tevent\tmatrix\tread\tread_index\tallele\tscore"
    )?;

    let reader = text_reader(path)?;
    let mut current_region = String::new();
    let mut hap_stage = String::new();
    let mut hap_index = 0usize;
    let mut pending_haplotype: Option<DebugHaplotypeMeta> = None;
    let mut pending_read: Option<DebugReadMeta> = None;
    let mut current_event = String::new();
    let mut matrix_mode = String::new();
    let mut matrix_alleles: Vec<String> = Vec::new();
    let mut awaiting_matrix_header = false;

    for line_result in reader.lines() {
        let line = line_result?;
        let trimmed = line.trim();
        if trimmed.is_empty() {
            pending_haplotype = None;
            pending_read = None;
            awaiting_matrix_header = false;
            matrix_mode.clear();
            continue;
        }

        if let Some(region) = trimmed.strip_prefix("calling for region: ") {
            current_region = region.to_string();
            continue;
        }
        if let Some(region) = trimmed.strip_prefix("assemblyRegion: ") {
            current_region = region.to_string();
            continue;
        }
        if trimmed.starts_with("Unclipped Haplotypes(") {
            hap_stage = "unclipped".to_string();
            hap_index = 0;
            pending_haplotype = None;
            continue;
        }
        if trimmed.starts_with("Clipped Haplot") {
            hap_stage = "clipped".to_string();
            hap_index = 0;
            pending_haplotype = None;
            continue;
        }
        if !hap_stage.is_empty() && trimmed.starts_with('[') && trimmed.contains(" k=") {
            pending_haplotype = parse_haplotype_meta(trimmed);
            continue;
        }
        if let Some(meta) = pending_haplotype.take() {
            writeln!(
                haplotypes,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                current_region,
                hap_stage,
                hap_index,
                meta.start,
                meta.end,
                meta.kmer,
                meta.length,
                meta.cigar,
                meta.is_ref,
                trimmed
            )?;
            stats.genotyper_haplotypes += 1;
            hap_index += 1;
            continue;
        }

        if let Some(event_raw) = trimmed.strip_prefix("Event at: ") {
            let (event_id, chrom, pos, event_type, alleles) = parse_event_line(event_raw);
            current_event = event_id;
            writeln!(
                events,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                current_region, current_event, chrom, pos, event_type, alleles, event_raw
            )?;
            matrix_mode.clear();
            awaiting_matrix_header = false;
            continue;
        }
        if let Some(rest) = trimmed.strip_prefix("Allele: ") {
            if !current_event.is_empty() {
                let (allele, haps) = parse_allele_hap_line(rest);
                writeln!(
                    event_haps,
                    "{}\t{}\t{}\t{}",
                    current_region, current_event, allele, haps
                )?;
                stats.event_allele_links += 1;
            }
            continue;
        }
        if trimmed == "Read-allele matrix:" {
            matrix_mode = "raw".to_string();
            matrix_alleles.clear();
            awaiting_matrix_header = true;
            continue;
        }
        if trimmed == "Normalized Read-Allele matrix:" {
            matrix_mode = "normalized".to_string();
            matrix_alleles.clear();
            awaiting_matrix_header = true;
            continue;
        }
        if awaiting_matrix_header {
            matrix_alleles = trimmed.split_whitespace().map(|s| s.to_string()).collect();
            awaiting_matrix_header = false;
            continue;
        }
        if !matrix_mode.is_empty() {
            if let Some((read_index, read_name, scores)) = parse_read_allele_matrix_row(trimmed) {
                for (allele_index, score) in scores.iter().enumerate() {
                    let allele = matrix_alleles
                        .get(allele_index)
                        .cloned()
                        .unwrap_or_else(|| allele_index.to_string());
                    writeln!(
                        allele_likelihoods,
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        current_region,
                        current_event,
                        matrix_mode,
                        read_name,
                        read_index,
                        allele,
                        score
                    )?;
                    stats.allele_likelihoods += 1;
                }
                continue;
            }
            matrix_mode.clear();
        }

        if trimmed.starts_with("read ") && trimmed.contains(" cigar: ") {
            pending_read = parse_read_meta(trimmed);
            continue;
        }
        if let Some(read_meta) = pending_read.take() {
            if trimmed.starts_with('[') {
                writeln!(
                    read_qualities,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    current_region,
                    read_meta.name,
                    read_meta.index,
                    read_meta.cigar,
                    read_meta.mapq,
                    read_meta.loc,
                    read_meta.unclipped_loc,
                    read_meta.length,
                    trimmed
                )?;
                stats.read_quality_rows += 1;
                continue;
            }
            if trimmed.starts_with(',') {
                for (haplotype_index, score) in trimmed
                    .split(',')
                    .filter(|value| !value.is_empty())
                    .enumerate()
                {
                    writeln!(
                        pairhmm,
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        current_region,
                        read_meta.name,
                        haplotype_index,
                        read_meta.index,
                        read_meta.cigar,
                        read_meta.mapq,
                        read_meta.loc,
                        read_meta.unclipped_loc,
                        read_meta.length,
                        score
                    )?;
                    stats.pairhmm_scores += 1;
                }
            }
        }
    }
    Ok(())
}

fn extract_assembly_state(
    path: &Path,
    prefix: &Path,
    stats: &mut GatkDebugExtractStats,
) -> Result<()> {
    let mut regions = File::create(prefixed_path(prefix, "assembly_regions.tsv"))?;
    let mut reads = File::create(prefixed_path(prefix, "assembly_reads.tsv"))?;
    let mut haplotypes = File::create(prefixed_path(prefix, "assembly_haplotypes.tsv"))?;
    writeln!(regions, "region\tmetric\tvalue")?;
    writeln!(reads, "region\tread\tflags")?;
    writeln!(haplotypes, "region\thaplotype\tbases")?;

    let reader = text_reader(path)?;
    let mut current_region = String::new();
    let mut in_reads = false;
    let mut in_haplotypes = false;
    let mut hap_index = 0usize;
    for line_result in reader.lines() {
        let line = line_result?;
        let trimmed = line.trim();
        if trimmed.is_empty() {
            in_reads = false;
            in_haplotypes = false;
            continue;
        }
        if looks_like_interval(trimmed) {
            current_region = trimmed.to_string();
            in_reads = false;
            in_haplotypes = false;
            hap_index = 0;
            continue;
        }
        if let Some(count) = parse_number_of_reads(trimmed) {
            writeln!(regions, "{}\treads\t{}", current_region, count)?;
            in_reads = true;
            in_haplotypes = false;
            continue;
        }
        if let Some(count) = parse_number_of_haplotypes(trimmed) {
            writeln!(regions, "{}\thaplotypes\t{}", current_region, count)?;
            in_reads = false;
            in_haplotypes = true;
            hap_index = 0;
            continue;
        }
        if in_reads {
            let mut parts = trimmed.split_whitespace();
            if let (Some(read), Some(flags)) = (parts.next(), parts.next()) {
                writeln!(reads, "{}\t{}\t{}", current_region, read, flags)?;
                stats.assembly_reads += 1;
            }
            continue;
        }
        if in_haplotypes {
            writeln!(haplotypes, "{}\t{}\t{}", current_region, hap_index, trimmed)?;
            stats.assembly_haplotypes += 1;
            hap_index += 1;
        }
    }
    Ok(())
}

fn parse_vcf_records(line: &str) -> Result<Vec<VcfRecord>> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 8 {
        bail!("VCF record has fewer than 8 columns: {line}");
    }
    let pos = fields[1]
        .parse::<u64>()
        .with_context(|| format!("invalid VCF position {}", fields[1]))?;
    let gt = if fields.len() >= 10 {
        parse_gt(fields[8], fields[9])
    } else {
        String::new()
    };
    let info = parse_info(fields[7]);
    Ok(fields[4]
        .split(',')
        .map(|alt| VcfRecord {
            key: VcfKey {
                chrom: fields[0].to_string(),
                pos,
                ref_allele: fields[3].to_string(),
                alt: alt.to_string(),
            },
            qual: fields[5].to_string(),
            filter: fields[6].to_string(),
            info: info.clone(),
            gt: gt.clone(),
        })
        .collect())
}

fn parse_gt(format: &str, sample: &str) -> String {
    let keys: Vec<&str> = format.split(':').collect();
    let Some(gt_index) = keys.iter().position(|key| *key == "GT") else {
        return String::new();
    };
    sample.split(':').nth(gt_index).unwrap_or("").to_string()
}

fn parse_sample_value_map(format: &str, sample: &str) -> BTreeMap<String, String> {
    let keys: Vec<&str> = format.split(':').collect();
    let values: Vec<&str> = sample.split(':').collect();
    let mut parsed = BTreeMap::new();
    for (index, key) in keys.iter().enumerate() {
        parsed.insert(
            (*key).to_string(),
            values.get(index).copied().unwrap_or("").to_string(),
        );
    }
    parsed
}

fn split_ad(ad: &str) -> (String, String) {
    if ad.is_empty() {
        return (String::new(), String::new());
    }
    let values: Vec<&str> = ad.split(',').collect();
    let ad_ref = values.first().copied().unwrap_or("").to_string();
    let ad_alt = if values.len() <= 1 {
        String::new()
    } else {
        values[1..].join(",")
    };
    (ad_ref, ad_alt)
}

fn parse_info(info: &str) -> BTreeMap<String, String> {
    let mut parsed = BTreeMap::new();
    if info == "." {
        return parsed;
    }
    for item in info.split(';') {
        if item.is_empty() {
            continue;
        }
        if let Some((key, value)) = item.split_once('=') {
            parsed.insert(key.to_string(), value.to_string());
        } else {
            parsed.insert(item.to_string(), "true".to_string());
        }
    }
    parsed
}

fn text_reader(path: &Path) -> Result<Box<dyn BufRead>> {
    let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let reader: Box<dyn Read> = if is_gzip_path(path) {
        Box::new(MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    };
    Ok(Box::new(BufReader::new(reader)))
}

fn is_gzip_path(path: &Path) -> bool {
    let name = path.file_name().and_then(|n| n.to_str()).unwrap_or("");
    name.ends_with(".gz") || name.ends_with(".bgz")
}

fn record_map(records: &[VcfRecord]) -> BTreeMap<VcfKey, VcfRecord> {
    let mut map = BTreeMap::new();
    for record in records {
        map.insert(record.key.clone(), record.clone());
    }
    map
}

fn summarize_vcf_records(records: &[VcfRecord]) -> VcfDatasetSummary {
    let mut summary = VcfDatasetSummary {
        total: records.len(),
        ..VcfDatasetSummary::default()
    };
    for record in records {
        *summary.filters.entry(record.filter.clone()).or_default() += 1;
        *summary
            .all_types
            .entry(record.variant_type().to_string())
            .or_default() += 1;
        if record.is_pass() {
            summary.pass += 1;
            *summary
                .pass_types
                .entry(record.variant_type().to_string())
                .or_default() += 1;
        }
    }
    summary.nonpass = summary.total - summary.pass;
    summary
}

fn compare_record_maps(
    a_map: &BTreeMap<VcfKey, VcfRecord>,
    b_map: &BTreeMap<VcfKey, VcfRecord>,
    pass_only: bool,
) -> VcfSetComparison {
    let a_keys = filtered_keys(a_map, pass_only);
    let b_keys = filtered_keys(b_map, pass_only);
    let shared_keys: BTreeSet<_> = a_keys.intersection(&b_keys).cloned().collect();
    let a_private_keys: BTreeSet<_> = a_keys.difference(&b_keys).cloned().collect();
    let b_private_keys: BTreeSet<_> = b_keys.difference(&a_keys).cloned().collect();
    let mut comparison = VcfSetComparison {
        a_count: a_keys.len(),
        b_count: b_keys.len(),
        shared: shared_keys.len(),
        a_private: a_private_keys.len(),
        b_private: b_private_keys.len(),
        a_sensitivity: percent(shared_keys.len(), a_keys.len()),
        b_precision_vs_a: percent(shared_keys.len(), b_keys.len()),
        ..VcfSetComparison::default()
    };
    for key in &shared_keys {
        *comparison
            .shared_types
            .entry(variant_type(&key.ref_allele, &key.alt).to_string())
            .or_default() += 1;
        let a_gt = a_map
            .get(key)
            .map(|record| record.gt.as_str())
            .unwrap_or("");
        let b_gt = b_map
            .get(key)
            .map(|record| record.gt.as_str())
            .unwrap_or("");
        if a_gt == b_gt {
            comparison.gt_same += 1;
        } else {
            comparison.gt_diff += 1;
        }
    }
    for key in &a_private_keys {
        *comparison
            .a_private_types
            .entry(variant_type(&key.ref_allele, &key.alt).to_string())
            .or_default() += 1;
    }
    for key in &b_private_keys {
        *comparison
            .b_private_types
            .entry(variant_type(&key.ref_allele, &key.alt).to_string())
            .or_default() += 1;
    }
    comparison
}

fn build_detail_rows(
    a_map: &BTreeMap<VcfKey, VcfRecord>,
    b_map: &BTreeMap<VcfKey, VcfRecord>,
) -> Vec<VcfDetailRow> {
    let keys: BTreeSet<_> = a_map.keys().chain(b_map.keys()).cloned().collect();
    let mut rows = Vec::with_capacity(keys.len());
    for key in keys {
        let a = a_map.get(&key);
        let b = b_map.get(&key);
        let category = match (a, b) {
            (Some(a_record), Some(b_record)) if a_record.gt == b_record.gt => "shared_gt_same",
            (Some(_), Some(_)) => "shared_gt_diff",
            (Some(_), None) => "a_private",
            (None, Some(_)) => "b_private",
            (None, None) => unreachable!(),
        };
        rows.push(VcfDetailRow {
            category: category.to_string(),
            variant_type: variant_type(&key.ref_allele, &key.alt),
            a_filter: a.map(|record| record.filter.clone()).unwrap_or_default(),
            b_filter: b.map(|record| record.filter.clone()).unwrap_or_default(),
            a_gt: a.map(|record| record.gt.clone()).unwrap_or_default(),
            b_gt: b.map(|record| record.gt.clone()).unwrap_or_default(),
            key,
        });
    }
    rows
}

fn filtered_keys(map: &BTreeMap<VcfKey, VcfRecord>, pass_only: bool) -> BTreeSet<VcfKey> {
    map.iter()
        .filter_map(|(key, record)| {
            if !pass_only || record.is_pass() {
                Some(key.clone())
            } else {
                None
            }
        })
        .collect()
}

pub fn variant_type(ref_allele: &str, alt: &str) -> &'static str {
    let alts: Vec<&str> = alt.split(',').collect();
    if ref_allele.len() == 1 && alts.iter().all(|allele| allele.len() == 1) {
        "SNP"
    } else {
        "INDEL_OR_COMPLEX"
    }
}

fn percent(numerator: usize, denominator: usize) -> f64 {
    if denominator == 0 {
        0.0
    } else {
        numerator as f64 / denominator as f64 * 100.0
    }
}

fn write_vcf_summary_markdown(path: &Path, comparison: &VcfComparison) -> Result<()> {
    let mut out = File::create(path).with_context(|| format!("creating {}", path.display()))?;
    writeln!(
        out,
        "# VCF comparison: {} vs {}",
        comparison.a_label, comparison.b_label
    )?;
    writeln!(out)?;
    writeln!(out, "## Counts")?;
    writeln!(out)?;
    writeln!(
        out,
        "| dataset | total allele keys | PASS allele keys | non-PASS allele keys | PASS SNP | PASS indel/complex | all SNP | all indel/complex | filters |"
    )?;
    writeln!(out, "|---|---:|---:|---:|---:|---:|---:|---:|---|")?;
    write_summary_count_row(&mut out, &comparison.a_label, &comparison.a_summary)?;
    write_summary_count_row(&mut out, &comparison.b_label, &comparison.b_summary)?;
    writeln!(out)?;
    writeln!(out, "## Exact Allele-Key Comparison")?;
    writeln!(out)?;
    writeln!(
        out,
        "| set | A count | B count | shared | A private | B private | A sensitivity | B precision vs A | shared SNP | shared indel/complex | A-private SNP | A-private indel/complex | B-private SNP | B-private indel/complex | GT same | GT diff |"
    )?;
    writeln!(
        out,
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|"
    )?;
    write_comparison_row(&mut out, "all allele keys", &comparison.all_records)?;
    write_comparison_row(&mut out, "PASS only", &comparison.pass_records)?;
    Ok(())
}

fn write_summary_count_row(
    out: &mut dyn Write,
    label: &str,
    summary: &VcfDatasetSummary,
) -> Result<()> {
    writeln!(
        out,
        "| {label} | {} | {} | {} | {} | {} | {} | {} | {} |",
        summary.total,
        summary.pass,
        summary.nonpass,
        map_get(&summary.pass_types, "SNP"),
        map_get(&summary.pass_types, "INDEL_OR_COMPLEX"),
        map_get(&summary.all_types, "SNP"),
        map_get(&summary.all_types, "INDEL_OR_COMPLEX"),
        format_counter(&summary.filters)
    )?;
    Ok(())
}

fn write_comparison_row(out: &mut dyn Write, label: &str, comp: &VcfSetComparison) -> Result<()> {
    writeln!(
        out,
        "| {label} | {} | {} | {} | {} | {} | {:.3}% | {:.3}% | {} | {} | {} | {} | {} | {} | {} | {} |",
        comp.a_count,
        comp.b_count,
        comp.shared,
        comp.a_private,
        comp.b_private,
        comp.a_sensitivity,
        comp.b_precision_vs_a,
        map_get(&comp.shared_types, "SNP"),
        map_get(&comp.shared_types, "INDEL_OR_COMPLEX"),
        map_get(&comp.a_private_types, "SNP"),
        map_get(&comp.a_private_types, "INDEL_OR_COMPLEX"),
        map_get(&comp.b_private_types, "SNP"),
        map_get(&comp.b_private_types, "INDEL_OR_COMPLEX"),
        comp.gt_same,
        comp.gt_diff,
    )?;
    Ok(())
}

fn write_vcf_summary_tsv(path: &Path, comparison: &VcfComparison) -> Result<()> {
    let mut out = File::create(path).with_context(|| format!("creating {}", path.display()))?;
    writeln!(
        out,
        "set\ta_count\tb_count\tshared\ta_private\tb_private\ta_sensitivity\tb_precision_vs_a\tshared_snp\tshared_indel_or_complex\ta_private_snp\ta_private_indel_or_complex\tb_private_snp\tb_private_indel_or_complex\tgt_same\tgt_diff"
    )?;
    write_comparison_tsv(&mut out, "all", &comparison.all_records)?;
    write_comparison_tsv(&mut out, "pass", &comparison.pass_records)?;
    Ok(())
}

fn write_comparison_tsv(out: &mut dyn Write, label: &str, comp: &VcfSetComparison) -> Result<()> {
    writeln!(
        out,
        "{label}\t{}\t{}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        comp.a_count,
        comp.b_count,
        comp.shared,
        comp.a_private,
        comp.b_private,
        comp.a_sensitivity,
        comp.b_precision_vs_a,
        map_get(&comp.shared_types, "SNP"),
        map_get(&comp.shared_types, "INDEL_OR_COMPLEX"),
        map_get(&comp.a_private_types, "SNP"),
        map_get(&comp.a_private_types, "INDEL_OR_COMPLEX"),
        map_get(&comp.b_private_types, "SNP"),
        map_get(&comp.b_private_types, "INDEL_OR_COMPLEX"),
        comp.gt_same,
        comp.gt_diff,
    )?;
    Ok(())
}

fn write_vcf_summary_json(path: &Path, comparison: &VcfComparison) -> Result<()> {
    let mut out = File::create(path).with_context(|| format!("creating {}", path.display()))?;
    writeln!(out, "{{")?;
    writeln!(
        out,
        "  \"a_label\": \"{}\",",
        json_escape(&comparison.a_label)
    )?;
    writeln!(
        out,
        "  \"b_label\": \"{}\",",
        json_escape(&comparison.b_label)
    )?;
    write_json_comparison(&mut out, "all_records", &comparison.all_records, true)?;
    write_json_comparison(&mut out, "pass_records", &comparison.pass_records, false)?;
    writeln!(out, "}}")?;
    Ok(())
}

fn write_json_comparison(
    out: &mut dyn Write,
    name: &str,
    comp: &VcfSetComparison,
    trailing_comma: bool,
) -> Result<()> {
    writeln!(out, "  \"{name}\": {{")?;
    writeln!(out, "    \"a_count\": {},", comp.a_count)?;
    writeln!(out, "    \"b_count\": {},", comp.b_count)?;
    writeln!(out, "    \"shared\": {},", comp.shared)?;
    writeln!(out, "    \"a_private\": {},", comp.a_private)?;
    writeln!(out, "    \"b_private\": {},", comp.b_private)?;
    writeln!(out, "    \"a_sensitivity\": {:.6},", comp.a_sensitivity)?;
    writeln!(
        out,
        "    \"b_precision_vs_a\": {:.6},",
        comp.b_precision_vs_a
    )?;
    writeln!(out, "    \"gt_same\": {},", comp.gt_same)?;
    writeln!(out, "    \"gt_diff\": {}", comp.gt_diff)?;
    writeln!(out, "  }}{}", if trailing_comma { "," } else { "" })?;
    Ok(())
}

fn write_vcf_detail_tsv(path: &Path, comparison: &VcfComparison) -> Result<()> {
    let mut out = File::create(path).with_context(|| format!("creating {}", path.display()))?;
    writeln!(
        out,
        "category\tchrom\tpos\tref\talt\ttype\ta_filter\tb_filter\ta_gt\tb_gt"
    )?;
    for row in &comparison.detail_rows {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            row.category,
            row.key.chrom,
            row.key.pos,
            row.key.ref_allele,
            row.key.alt,
            row.variant_type,
            row.a_filter,
            row.b_filter,
            row.a_gt,
            row.b_gt
        )?;
    }
    Ok(())
}

fn map_get(map: &BTreeMap<String, usize>, key: &str) -> usize {
    map.get(key).copied().unwrap_or(0)
}

fn format_counter(map: &BTreeMap<String, usize>) -> String {
    map.iter()
        .map(|(key, value)| format!("{key}:{value}"))
        .collect::<Vec<_>>()
        .join(", ")
}

fn json_escape(value: &str) -> String {
    value
        .chars()
        .flat_map(|c| match c {
            '"' => "\\\"".chars().collect::<Vec<_>>(),
            '\\' => "\\\\".chars().collect(),
            '\n' => "\\n".chars().collect(),
            '\r' => "\\r".chars().collect(),
            '\t' => "\\t".chars().collect(),
            _ => vec![c],
        })
        .collect()
}

fn push_key_category<'a, I>(
    rows: &mut Vec<SelectedRegion>,
    category: &str,
    keys: I,
    a_map: &BTreeMap<VcfKey, VcfRecord>,
    b_map: &BTreeMap<VcfKey, VcfRecord>,
    padding: u64,
    max_per_category: usize,
) where
    I: Iterator<Item = &'a VcfKey>,
{
    for key in keys.take(max_per_category) {
        let allele_len = key
            .alt
            .split(',')
            .map(str::len)
            .chain(std::iter::once(key.ref_allele.len()))
            .max()
            .unwrap_or(1) as u64;
        rows.push(SelectedRegion {
            category: category.to_string(),
            variant_type: variant_type(&key.ref_allele, &key.alt),
            start: key.pos.saturating_sub(padding).max(1),
            end: key.pos + allele_len.saturating_sub(1) + padding,
            a_gt: a_map
                .get(key)
                .map(|record| record.gt.clone())
                .unwrap_or_default(),
            b_gt: b_map
                .get(key)
                .map(|record| record.gt.clone())
                .unwrap_or_default(),
            key: key.clone(),
        });
    }
}

fn write_region_manifest(path: &Path, rows: &[SelectedRegion]) -> Result<()> {
    let mut out = File::create(path).with_context(|| format!("creating {}", path.display()))?;
    writeln!(
        out,
        "category\tchrom\tstart\tend\tpos\tref\talt\ttype\ta_gt\tb_gt"
    )?;
    for row in rows {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            row.category,
            row.key.chrom,
            row.start,
            row.end,
            row.key.pos,
            row.key.ref_allele,
            row.key.alt,
            row.variant_type,
            row.a_gt,
            row.b_gt
        )?;
    }
    Ok(())
}

fn write_region_bed(path: &Path, rows: &[SelectedRegion]) -> Result<()> {
    let mut out = File::create(path).with_context(|| format!("creating {}", path.display()))?;
    for row in rows {
        writeln!(
            out,
            "{}\t{}\t{}\t{}_{}",
            row.key.chrom,
            row.start.saturating_sub(1),
            row.end,
            row.category,
            row.key.display()
        )?;
    }
    Ok(())
}

fn write_region_interval_list(path: &Path, rows: &[SelectedRegion], template: &Path) -> Result<()> {
    let mut out = File::create(path).with_context(|| format!("creating {}", path.display()))?;
    let template_reader = text_reader(template)?;
    for line in template_reader.lines() {
        let line = line?;
        if line.starts_with('@') {
            writeln!(out, "{line}")?;
        }
    }
    for row in rows {
        writeln!(
            out,
            "{}\t{}\t{}\t+\t{}",
            row.key.chrom, row.start, row.end, row.category
        )?;
    }
    Ok(())
}

fn read_stage_rows(path: &Path) -> Result<Vec<StageRow>> {
    let mut reader = text_reader(path)?;
    let mut first = String::new();
    loop {
        first.clear();
        if reader.read_line(&mut first)? == 0 {
            return Ok(Vec::new());
        }
        if !first.trim().is_empty() {
            break;
        }
    }

    if first.trim_start().starts_with('{') {
        read_jsonl_rows(first, reader)
    } else {
        read_tsv_rows(first, reader)
    }
}

fn read_jsonl_rows(first: String, mut reader: Box<dyn BufRead>) -> Result<Vec<StageRow>> {
    let mut rows = Vec::new();
    parse_jsonl_row(&first, &mut rows)?;
    let mut line = String::new();
    while reader.read_line(&mut line)? != 0 {
        if !line.trim().is_empty() {
            parse_jsonl_row(&line, &mut rows)?;
        }
        line.clear();
    }
    Ok(rows)
}

fn parse_jsonl_row(line: &str, rows: &mut Vec<StageRow>) -> Result<()> {
    let value: Value = serde_json::from_str(line.trim()).context("parsing JSONL stage row")?;
    let Some(object) = value.as_object() else {
        bail!("JSONL stage row is not an object");
    };
    let mut fields = BTreeMap::new();
    for (key, value) in object {
        let text = match value {
            Value::Null => String::new(),
            Value::Bool(v) => v.to_string(),
            Value::Number(v) => v.to_string(),
            Value::String(v) => v.clone(),
            Value::Array(_) | Value::Object(_) => value.to_string(),
        };
        fields.insert(key.clone(), text);
    }
    rows.push(StageRow { fields });
    Ok(())
}

fn read_tsv_rows(first: String, mut reader: Box<dyn BufRead>) -> Result<Vec<StageRow>> {
    let header: Vec<String> = first
        .trim_end()
        .split('\t')
        .map(|field| field.to_string())
        .collect();
    let mut rows = Vec::new();
    let mut line = String::new();
    while reader.read_line(&mut line)? != 0 {
        let trimmed = line.trim_end();
        if trimmed.is_empty() {
            line.clear();
            continue;
        }
        let values: Vec<&str> = trimmed.split('\t').collect();
        let mut fields = BTreeMap::new();
        for (index, key) in header.iter().enumerate() {
            fields.insert(
                key.clone(),
                values.get(index).copied().unwrap_or("").to_string(),
            );
        }
        rows.push(StageRow { fields });
        line.clear();
    }
    Ok(rows)
}

fn stage_row_map(
    rows: &[StageRow],
    key_columns: &[String],
    label: &str,
) -> Result<BTreeMap<String, StageRow>> {
    if key_columns.is_empty() {
        bail!("at least one key column is required");
    }
    let mut map = BTreeMap::new();
    for row in rows {
        let mut key_parts = Vec::new();
        for column in key_columns {
            let Some(value) = row.fields.get(column) else {
                bail!("{label} row is missing key column {column}");
            };
            key_parts.push(value.clone());
        }
        map.insert(key_parts.join("\x1f"), row.clone());
    }
    Ok(map)
}

fn diff_stage_fields(
    java: &StageRow,
    rust: &StageRow,
    key_columns: &[String],
    tolerance: f64,
) -> Vec<(String, String, String)> {
    let key_set: BTreeSet<_> = key_columns.iter().cloned().collect();
    let fields: BTreeSet<_> = java
        .fields
        .keys()
        .chain(rust.fields.keys())
        .filter(|field| !key_set.contains(*field))
        .cloned()
        .collect();
    let mut diffs = Vec::new();
    for field in fields {
        let java_value = java.fields.get(&field).cloned().unwrap_or_default();
        let rust_value = rust.fields.get(&field).cloned().unwrap_or_default();
        if values_match(&java_value, &rust_value, tolerance) {
            continue;
        }
        diffs.push((field, java_value, rust_value));
    }
    diffs
}

fn values_match(a: &str, b: &str, tolerance: f64) -> bool {
    if a == b {
        return true;
    }
    match (a.parse::<f64>(), b.parse::<f64>()) {
        (Ok(a_value), Ok(b_value)) => {
            if a_value.is_nan() && b_value.is_nan() {
                true
            } else {
                (a_value - b_value).abs() <= tolerance
            }
        }
        _ => false,
    }
}

fn write_stage_summary_markdown(
    path: &Path,
    config: &StageDiffConfig,
    summary: &StageDiffSummary,
) -> Result<()> {
    let mut out = File::create(path).with_context(|| format!("creating {}", path.display()))?;
    writeln!(out, "# {} stage diff", config.stage_name)?;
    writeln!(out)?;
    writeln!(out, "- Java file: `{}`", config.java_path.display())?;
    writeln!(out, "- Rust file: `{}`", config.rust_path.display())?;
    writeln!(out, "- key columns: `{}`", config.key_columns.join(","))?;
    writeln!(out, "- numeric tolerance: `{}`", config.numeric_tolerance)?;
    writeln!(out)?;
    writeln!(out, "| metric | value |")?;
    writeln!(out, "|---|---:|")?;
    writeln!(out, "| Java rows | {} |", summary.java_rows)?;
    writeln!(out, "| Rust rows | {} |", summary.rust_rows)?;
    writeln!(out, "| shared rows | {} |", summary.shared_rows)?;
    writeln!(out, "| Java-private rows | {} |", summary.java_private_rows)?;
    writeln!(out, "| Rust-private rows | {} |", summary.rust_private_rows)?;
    writeln!(
        out,
        "| shared rows with field diffs | {} |",
        summary.shared_rows_with_diffs
    )?;
    writeln!(out, "| field diffs | {} |", summary.field_diffs)?;
    Ok(())
}

fn write_stage_summary_tsv(
    path: &Path,
    config: &StageDiffConfig,
    summary: &StageDiffSummary,
) -> Result<()> {
    let mut out = File::create(path).with_context(|| format!("creating {}", path.display()))?;
    writeln!(
        out,
        "stage\tjava_rows\trust_rows\tshared_rows\tjava_private_rows\trust_private_rows\tshared_rows_with_diffs\tfield_diffs"
    )?;
    writeln!(
        out,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        config.stage_name,
        summary.java_rows,
        summary.rust_rows,
        summary.shared_rows,
        summary.java_private_rows,
        summary.rust_private_rows,
        summary.shared_rows_with_diffs,
        summary.field_diffs
    )?;
    Ok(())
}

fn write_stage_detail_tsv(path: &Path, rows: &[Vec<String>]) -> Result<()> {
    let mut out = File::create(path).with_context(|| format!("creating {}", path.display()))?;
    writeln!(out, "category\trow_key\tfield\tjava_value\trust_value")?;
    for row in rows {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}",
            row.first().cloned().unwrap_or_default(),
            row.get(1).cloned().unwrap_or_default(),
            row.get(2).cloned().unwrap_or_default(),
            row.get(3).cloned().unwrap_or_default(),
            row.get(4).cloned().unwrap_or_default()
        )?;
    }
    Ok(())
}

fn create_output_parent(prefix: &Path) -> Result<()> {
    if let Some(parent) = prefix.parent() {
        if !parent.as_os_str().is_empty() {
            fs::create_dir_all(parent).with_context(|| format!("creating {}", parent.display()))?;
        }
    }
    Ok(())
}

fn parse_haplotype_meta(line: &str) -> Option<DebugHaplotypeMeta> {
    let tokens: Vec<&str> = line.split_whitespace().collect();
    if tokens.len() < 5 {
        return None;
    }
    let span = tokens[0].trim_start_matches('[').trim_end_matches(']');
    let (start, end) = span.split_once('-')?;
    let kmer = tokens[1].strip_prefix("k=")?.to_string();
    let length = tokens[3].to_string();
    let cigar_token = tokens[4];
    let is_ref = cigar_token.ends_with("ref");
    let cigar = cigar_token.trim_end_matches("ref").to_string();
    Some(DebugHaplotypeMeta {
        start: start.to_string(),
        end: end.to_string(),
        kmer,
        length,
        cigar,
        is_ref,
    })
}

fn parse_read_meta(line: &str) -> Option<DebugReadMeta> {
    let rest = line.strip_prefix("read ")?;
    let (index_text, rest) = rest.split_once(": ")?;
    let index = index_text.parse::<usize>().ok()?;
    let (name, rest) = rest.split_once(" cigar: ")?;
    let (cigar, rest) = rest.split_once(" mapQ: ")?;
    let (mapq, rest) = rest.split_once(" loc: ")?;
    let (loc, rest) = rest.split_once(" unclippedloc: ")?;
    let (unclipped_loc, length) = if let Some((unclipped, length)) = rest.split_once(" length:") {
        (unclipped, length)
    } else {
        (rest, "")
    };
    Some(DebugReadMeta {
        index,
        name: name.to_string(),
        cigar: cigar.to_string(),
        mapq: mapq.to_string(),
        loc: loc.to_string(),
        unclipped_loc: unclipped_loc.to_string(),
        length: length.to_string(),
    })
}

fn parse_event_line(line: &str) -> (String, String, String, String, String) {
    let locus = extract_between(line, "@ ", " Q.").unwrap_or_default();
    let (chrom, pos) = locus
        .split_once(':')
        .map(|(chrom, pos)| (chrom.to_string(), pos.to_string()))
        .unwrap_or_else(|| (String::new(), String::new()));
    let event_type = extract_between(line, "type=", " alleles=").unwrap_or_default();
    let alleles = extract_between(line, "alleles=[", "]").unwrap_or_default();
    let event_id = format!(
        "{}:{}:{}:{}",
        chrom,
        pos,
        event_type,
        alleles.replace(' ', "")
    );
    (event_id, chrom, pos, event_type, alleles)
}

fn parse_allele_hap_line(line: &str) -> (String, String) {
    if let Some((allele, haps)) = line.split_once(" Haps: ") {
        (allele.to_string(), haps.to_string())
    } else {
        (line.to_string(), String::new())
    }
}

fn parse_read_allele_matrix_row(line: &str) -> Option<(usize, String, Vec<String>)> {
    let rest = line.strip_prefix("read: ")?;
    let mut parts = rest.split_whitespace();
    let read_index = parts.next()?.parse::<usize>().ok()?;
    let read_name = parts.next()?.to_string();
    let scores = parts.map(|part| part.to_string()).collect::<Vec<_>>();
    Some((read_index, read_name, scores))
}

fn looks_like_interval(line: &str) -> bool {
    line.contains(':')
        && line.contains('-')
        && !line.contains(' ')
        && line
            .split_once(':')
            .and_then(|(_, rest)| rest.split_once('-'))
            .is_some()
}

fn parse_number_of_reads(line: &str) -> Option<String> {
    let rest = line.strip_prefix("Number of reads in region: ")?;
    rest.split_whitespace()
        .next()
        .map(|value| value.to_string())
}

fn parse_number_of_haplotypes(line: &str) -> Option<String> {
    let rest = line.strip_prefix("There were ")?;
    rest.split_whitespace()
        .next()
        .map(|value| value.to_string())
}

fn extract_between(line: &str, start: &str, end: &str) -> Option<String> {
    let after_start = line.split_once(start)?.1;
    let value = after_start.split_once(end)?.0;
    Some(value.to_string())
}

#[allow(dead_code)]
fn compare_contig_pos(a: &SelectedRegion, b: &SelectedRegion) -> Ordering {
    a.key
        .chrom
        .cmp(&b.key.chrom)
        .then_with(|| a.start.cmp(&b.start))
        .then_with(|| a.end.cmp(&b.end))
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

    #[test]
    fn vcf_compare_counts_private_shared_and_gt_diff() {
        let dir = tempdir().unwrap();
        let a = dir.path().join("a.vcf");
        let b = dir.path().join("b.vcf");
        fs::write(
            &a,
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\nchr1\t10\t.\tA\tG\t50\tPASS\tDP=10\tGT:DP\t0/1:10\nchr1\t20\t.\tAT\tA\t50\tPASS\tDP=5\tGT:DP\t1/1:5\n",
        )
        .unwrap();
        fs::write(
            &b,
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\nchr1\t10\t.\tA\tG\t50\tPASS\tDP=10\tGT:DP\t1/1:10\nchr1\t30\t.\tC\tT\t50\tPASS\tDP=6\tGT:DP\t0/1:6\n",
        )
        .unwrap();
        let comparison = compare_vcfs(&a, &b, "a", "b").unwrap();
        assert_eq!(comparison.pass_records.a_count, 2);
        assert_eq!(comparison.pass_records.b_count, 2);
        assert_eq!(comparison.pass_records.shared, 1);
        assert_eq!(comparison.pass_records.a_private, 1);
        assert_eq!(comparison.pass_records.b_private, 1);
        assert_eq!(comparison.pass_records.gt_diff, 1);
        assert_eq!(
            map_get(&comparison.pass_records.a_private_types, "INDEL_OR_COMPLEX"),
            1
        );
    }

    #[test]
    fn vcf_compare_splits_multiallelic_alt_keys() {
        let dir = tempdir().unwrap();
        let a = dir.path().join("a.vcf");
        let b = dir.path().join("b.vcf");
        fs::write(
            &a,
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\nchr1\t10\t.\tA\tC,G\t50\tPASS\t.\tGT\t1/2\n",
        )
        .unwrap();
        fs::write(
            &b,
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\nchr1\t10\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n",
        )
        .unwrap();
        let comparison = compare_vcfs(&a, &b, "a", "b").unwrap();
        assert_eq!(comparison.pass_records.a_count, 2);
        assert_eq!(comparison.pass_records.b_count, 1);
        assert_eq!(comparison.pass_records.shared, 1);
        assert_eq!(comparison.pass_records.a_private, 1);
        assert_eq!(comparison.pass_records.b_private, 0);
        assert_eq!(map_get(&comparison.pass_records.shared_types, "SNP"), 1);
        assert_eq!(map_get(&comparison.pass_records.a_private_types, "SNP"), 1);
    }

    #[test]
    fn region_selection_splits_categories() {
        let dir = tempdir().unwrap();
        let a = dir.path().join("a.vcf");
        let b = dir.path().join("b.vcf");
        fs::write(
            &a,
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\nchr1\t10\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\nchr1\t20\t.\tAT\tA\t50\tPASS\t.\tGT\t1/1\n",
        )
        .unwrap();
        fs::write(
            &b,
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\nchr1\t10\t.\tA\tG\t50\tPASS\t.\tGT\t1/1\nchr1\t30\t.\tC\tT\t50\tPASS\t.\tGT\t0/1\n",
        )
        .unwrap();
        let rows = select_regions(&a, &b, 5, 10, true).unwrap();
        let categories: BTreeSet<_> = rows.iter().map(|row| row.category.as_str()).collect();
        assert!(categories.contains("a_private_indel"));
        assert!(categories.contains("b_private_snp"));
        assert!(categories.contains("shared_gt_diff"));
    }

    #[test]
    fn stage_diff_accepts_numeric_tolerance() {
        let dir = tempdir().unwrap();
        let java = dir.path().join("java.tsv");
        let rust = dir.path().join("rust.tsv");
        fs::write(&java, "region\tread\tscore\nr1\tread1\t-10.0001\n").unwrap();
        fs::write(&rust, "region\tread\tscore\nr1\tread1\t-10.0002\n").unwrap();
        let config = StageDiffConfig {
            java_path: java,
            rust_path: rust,
            key_columns: vec!["region".to_string(), "read".to_string()],
            numeric_tolerance: 0.001,
            output_prefix: dir.path().join("diff"),
            stage_name: "pairhmm".to_string(),
        };
        let summary = run_stage_diff(&config).unwrap();
        assert_eq!(summary.shared_rows, 1);
        assert_eq!(summary.field_diffs, 0);
    }

    #[test]
    fn vcf_genotype_table_extracts_shared_fields() {
        let dir = tempdir().unwrap();
        let vcf = dir.path().join("calls.vcf");
        let output = dir.path().join("calls.tsv");
        fs::write(
            &vcf,
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\nchr1\t10\t.\tA\tG\t42.5\tPASS\tDP=12;FS=3.0;QD=4.2;DB\tGT:AD:DP:GQ:PL\t0/1:7,5:12:99:100,0,90\n",
        )
        .unwrap();
        let rows = write_vcf_genotype_table(&vcf, &output).unwrap();
        assert_eq!(rows, 1);
        let table = fs::read_to_string(output).unwrap();
        assert!(table
            .contains("chr1\t10\tA\tG\t42.5\tPASS\t0/1\t99\t12\t7\t5\t3.0\t4.2\t100,0,90\ttrue"));
    }
}
