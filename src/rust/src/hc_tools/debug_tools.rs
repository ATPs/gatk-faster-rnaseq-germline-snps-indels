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

