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

