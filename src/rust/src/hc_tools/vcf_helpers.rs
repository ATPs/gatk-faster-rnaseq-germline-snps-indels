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

