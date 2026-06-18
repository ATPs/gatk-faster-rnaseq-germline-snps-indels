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

