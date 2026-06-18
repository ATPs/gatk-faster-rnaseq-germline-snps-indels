fn process_read(
    read: &PreparedRead,
    reference_bases: &[u8],
    known_intervals: &[KnownInterval],
    args: &Args,
    tables: &mut RecalTables,
) -> Result<()> {
    let contexts = context_values(read, args.context_size, args.low_quality_tail);
    let mut read_pos = 0_usize;
    let mut ref_pos = 0_usize;
    let mut known_index = 0_usize;

    for op in &read.cigar {
        match *op {
            SimpleCigar::Match(len) | SimpleCigar::Equal(len) | SimpleCigar::Diff(len) => {
                for _ in 0..len {
                    let base = read.bases[read_pos];
                    let qual = read.quals[read_pos].min(MAX_SAM_QUAL_SCORE as u8);
                    let ref_base = reference_bases.get(ref_pos).copied().with_context(|| {
                        format!(
                            "reference slice shorter than CIGAR span for {}:{}-{}",
                            read.contig,
                            read.start0 + 1,
                            read.end0 + 1
                        )
                    })?;
                    let absolute_ref = read.start0 + ref_pos as u64;
                    if is_regular_base(base)
                        && qual >= MIN_USABLE_Q_SCORE
                        && !is_known_site(absolute_ref, known_intervals, &mut known_index)
                    {
                        let is_error = if bases_equal(base, ref_base) {
                            0.0
                        } else {
                            1.0
                        };
                        update_tables(read, read_pos, qual, is_error, &contexts, args, tables)?;
                    }
                    read_pos += 1;
                    ref_pos += 1;
                }
            }
            SimpleCigar::Ins(len) => {
                read_pos += len as usize;
            }
            SimpleCigar::Del(len) | SimpleCigar::RefSkip(len) => {
                ref_pos += len as usize;
            }
        }
    }
    Ok(())
}

fn update_tables(
    read: &PreparedRead,
    read_pos: usize,
    qual: u8,
    is_error: f64,
    contexts: &[Option<String>],
    args: &Args,
    tables: &mut RecalTables,
) -> Result<()> {
    tables
        .quality
        .entry((read.rg_index, qual))
        .or_insert_with(|| Datum::new(qual))
        .increment(is_error);
    let cycle = cycle_value(
        read_pos,
        read.bases.len(),
        read.is_reverse,
        read.is_second_in_pair,
        args.max_cycle,
    )?;
    tables
        .cycle
        .entry((read.rg_index, qual, cycle))
        .or_insert_with(|| Datum::new(qual))
        .increment(is_error);
    if let Some(context) = &contexts[read_pos] {
        tables
            .context
            .entry((read.rg_index, qual, context.clone()))
            .or_insert_with(|| Datum::new(qual))
            .increment(is_error);
    }
    Ok(())
}

fn collapse_read_group_table(
    quality_table: &BTreeMap<(usize, u8), Datum>,
) -> BTreeMap<usize, Datum> {
    let mut read_groups: BTreeMap<usize, Datum> = BTreeMap::new();
    for ((rg, _qual), datum) in quality_table {
        read_groups
            .entry(*rg)
            .and_modify(|existing| existing.combine(datum))
            .or_insert_with(|| datum.clone());
    }
    read_groups
}

fn is_known_site(position0: u64, intervals: &[KnownInterval], index: &mut usize) -> bool {
    while *index < intervals.len() && intervals[*index].end0 < position0 {
        *index += 1;
    }
    *index < intervals.len()
        && intervals[*index].start0 <= position0
        && intervals[*index].end0 >= position0
}

fn context_values(
    read: &PreparedRead,
    context_size: usize,
    low_quality_tail: u8,
) -> Vec<Option<String>> {
    let mut bases = read.bases.clone();
    let mut left = 0_usize;
    while left < bases.len() && read.quals[left] <= low_quality_tail {
        bases[left] = b'N';
        left += 1;
    }
    let mut right = bases.len();
    while right > left && read.quals[right - 1] <= low_quality_tail {
        bases[right - 1] = b'N';
        right -= 1;
    }
    if read.is_reverse {
        bases = reverse_complement(&bases);
    }

    let mut stranded_contexts = vec![None; bases.len()];
    if bases.len() >= context_size {
        for i in (context_size - 1)..bases.len() {
            let start = i + 1 - context_size;
            if bases[start..=i].iter().all(|base| is_regular_base(*base)) {
                stranded_contexts[i] = Some(
                    bases[start..=i]
                        .iter()
                        .map(|base| normalize_base(*base) as char)
                        .collect(),
                );
            }
        }
    }

    let mut contexts = vec![None; bases.len()];
    for (stranded_offset, context) in stranded_contexts.into_iter().enumerate() {
        let read_offset = if read.is_reverse {
            bases.len() - stranded_offset - 1
        } else {
            stranded_offset
        };
        contexts[read_offset] = context;
    }
    contexts
}

fn cycle_value(
    base_number: usize,
    read_length: usize,
    is_reverse: bool,
    is_second_in_pair: bool,
    max_cycle: i32,
) -> Result<i32> {
    let read_order_factor = if is_second_in_pair { -1 } else { 1 };
    let (mut cycle, increment) = if is_reverse {
        (read_length as i32 * read_order_factor, -read_order_factor)
    } else {
        (read_order_factor, read_order_factor)
    };
    cycle += base_number as i32 * increment;
    if cycle.abs() > max_cycle {
        bail!(
            "cycle {} exceeds --maximum-cycle-value {}",
            cycle,
            max_cycle
        );
    }
    Ok(cycle)
}

fn is_regular_base(base: u8) -> bool {
    matches!(normalize_base(base), b'A' | b'C' | b'G' | b'T')
}

fn bases_equal(left: u8, right: u8) -> bool {
    normalize_base(left) == normalize_base(right)
}

fn normalize_base(base: u8) -> u8 {
    base.to_ascii_uppercase()
}

fn reverse_complement(bases: &[u8]) -> Vec<u8> {
    bases
        .iter()
        .rev()
        .map(|base| match normalize_base(*base) {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => b'N',
        })
        .collect()
}

fn build_quantization_table(
    quality_table: &BTreeMap<(usize, u8), Datum>,
    quantizing_levels: usize,
) -> Vec<(usize, u64, u8)> {
    let mut counts = vec![0_u64; MAX_SAM_QUAL_SCORE + 1];
    for datum in quality_table.values() {
        let empirical = empirical_quality(datum).round() as usize;
        counts[empirical.min(MAX_SAM_QUAL_SCORE)] += datum.observations;
    }
    let quantized =
        quantize_quality_scores(&counts, quantizing_levels, usize::from(MIN_USABLE_Q_SCORE));
    (0..=MAX_SAM_QUAL_SCORE)
        .map(|qual| (qual, counts[qual], quantized[qual]))
        .collect()
}

#[derive(Clone, Debug)]
struct QualInterval {
    q_start: usize,
    q_end: usize,
    fixed_qual: Option<u8>,
    observations: u64,
    errors: u64,
    sub_intervals: Vec<QualInterval>,
}

impl QualInterval {
    fn leaf(q: usize, observations: u64) -> Self {
        let errors = (observations as f64 * qual_to_error_prob(q as f64)).floor() as u64;
        Self {
            q_start: q,
            q_end: q,
            fixed_qual: Some(q as u8),
            observations,
            errors,
            sub_intervals: Vec::new(),
        }
    }

    fn merge(left: &Self, right: &Self) -> Self {
        Self {
            q_start: left.q_start,
            q_end: right.q_end,
            fixed_qual: None,
            observations: left.observations + right.observations,
            errors: left.errors + right.errors,
            sub_intervals: vec![left.clone(), right.clone()],
        }
    }

    fn error_rate(&self) -> f64 {
        if let Some(qual) = self.fixed_qual {
            qual_to_error_prob(f64::from(qual))
        } else if self.observations == 0 {
            0.0
        } else {
            (self.errors + 1) as f64 / (self.observations + 1) as f64
        }
    }

    fn qual(&self) -> u8 {
        if let Some(qual) = self.fixed_qual {
            qual
        } else {
            error_prob_to_qual(self.error_rate())
        }
    }

    fn penalty(&self, min_interesting_qual: usize) -> f64 {
        self.calc_penalty(self.error_rate(), min_interesting_qual)
    }

    fn calc_penalty(&self, global_error_rate: f64, min_interesting_qual: usize) -> f64 {
        if global_error_rate == 0.0 {
            return 0.0;
        }
        if self.sub_intervals.is_empty() {
            if self.q_end <= min_interesting_qual {
                0.0
            } else {
                (self.error_rate().log10() - global_error_rate.log10()).abs()
                    * self.observations as f64
            }
        } else {
            self.sub_intervals
                .iter()
                .map(|interval| interval.calc_penalty(global_error_rate, min_interesting_qual))
                .sum()
        }
    }
}

fn quantize_quality_scores(
    counts: &[u64],
    n_levels: usize,
    min_interesting_qual: usize,
) -> Vec<u8> {
    let mut intervals: Vec<QualInterval> = counts
        .iter()
        .enumerate()
        .map(|(qual, count)| QualInterval::leaf(qual, *count))
        .collect();

    while intervals.len() > n_levels {
        let mut best_index = 0_usize;
        let mut best_penalty = f64::INFINITY;
        for index in 0..(intervals.len() - 1) {
            let merged = QualInterval::merge(&intervals[index], &intervals[index + 1]);
            let penalty = merged.penalty(min_interesting_qual);
            if penalty < best_penalty {
                best_index = index;
                best_penalty = penalty;
            }
        }
        let merged = QualInterval::merge(&intervals[best_index], &intervals[best_index + 1]);
        intervals.splice(best_index..=best_index + 1, [merged]);
    }

    let mut map = vec![0_u8; counts.len()];
    for interval in intervals {
        for qual in interval.q_start..=interval.q_end {
            map[qual] = interval.qual();
        }
    }
    map
}

fn empirical_quality(datum: &Datum) -> f64 {
    let mismatches = (datum.errors + 0.5) as u64 + 1;
    let observations = datum.observations + 2;
    bayesian_estimate_of_empirical_quality(observations, mismatches, datum.reported_quality)
        .min(MAX_SAM_QUAL_SCORE as i32) as f64
}

fn bayesian_estimate_of_empirical_quality(
    observations: u64,
    errors: u64,
    prior_mean_quality: f64,
) -> i32 {
    let mut best_q = 0_i32;
    let mut best_value = f64::NEG_INFINITY;
    for q in 0..=MAX_REASONABLE_Q_SCORE {
        let value = log_prior(q as f64, prior_mean_quality)
            + log_binomial_likelihood(q as f64, observations, errors);
        if value > best_value {
            best_value = value;
            best_q = q;
        }
    }
    best_q
}

fn log_prior(quality_score: f64, prior_quality_score: f64) -> f64 {
    let difference = ((quality_score - prior_quality_score) as i32)
        .abs()
        .min(MAX_GATK_USABLE_Q_SCORE);
    let sigma = 0.5_f64;
    -((difference as f64).powi(2)) / (2.0 * sigma.powi(2))
}

fn log_binomial_likelihood(quality_score: f64, observations: u64, mut errors: u64) -> f64 {
    if observations == 0 {
        return 0.0;
    }
    let max_observations = i32::MAX as u64 - 1;
    let mut scaled_observations = observations;
    if scaled_observations > max_observations {
        let fraction = max_observations as f64 / scaled_observations as f64;
        errors = (errors as f64 * fraction).round() as u64;
        scaled_observations = max_observations;
    }

    let p = qual_to_error_prob(quality_score);
    if p == 1.0 {
        return if errors == scaled_observations {
            0.0
        } else {
            f64::NEG_INFINITY
        };
    }
    errors as f64 * p.ln() + (scaled_observations - errors) as f64 * (1.0 - p).ln()
}

fn qual_to_error_prob(qual: f64) -> f64 {
    10.0_f64.powf(qual / -10.0)
}

fn error_prob_to_qual(error_rate: f64) -> u8 {
    let qual = (-10.0 * error_rate.log10()).round() as i32;
    qual.clamp(1, MAX_SAM_QUAL_SCORE as i32) as u8
}

fn write_report(
    args: &Args,
    rg_info: &ReadGroupInfo,
    read_group_table: &BTreeMap<usize, Datum>,
    tables: &RecalTables,
    quantized: &[(usize, u64, u8)],
) -> Result<()> {
    if let Some(parent) = args.output_table.parent() {
        if !parent.as_os_str().is_empty() {
            fs::create_dir_all(parent)
                .with_context(|| format!("failed to create {}", parent.display()))?;
        }
    }
    let file = File::create(&args.output_table)
        .with_context(|| format!("failed to create {}", args.output_table.display()))?;
    let mut writer = BufWriter::new(file);
    writeln!(writer, "#:GATKReport.v1.1:5")?;

    write_table(
        &mut writer,
        "Arguments",
        "Recalibration argument collection values used in this run",
        &[("Argument", "%s"), ("Value", "%s")],
        argument_rows(args),
    )?;

    let quant_rows: Vec<Vec<String>> = quantized
        .iter()
        .map(|(qual, count, quantized_score)| {
            vec![
                qual.to_string(),
                count.to_string(),
                quantized_score.to_string(),
            ]
        })
        .collect();
    write_table(
        &mut writer,
        "Quantized",
        "Quality quantization map",
        &[
            ("QualityScore", "%d"),
            ("Count", "%d"),
            ("QuantizedScore", "%d"),
        ],
        quant_rows,
    )?;

    let rg_rows: Vec<Vec<String>> = read_group_table
        .iter()
        .map(|(rg, datum)| {
            vec![
                rg_info.identifiers[*rg].clone(),
                "M".to_string(),
                format!("{:.4}", empirical_quality(datum)),
                format!("{:.4}", datum.reported_quality),
                datum.observations.to_string(),
                format!("{:.2}", datum.errors),
            ]
        })
        .collect();
    write_table(
        &mut writer,
        "RecalTable0",
        "",
        &[
            ("ReadGroup", "%s"),
            ("EventType", "%s"),
            ("EmpiricalQuality", "%.4f"),
            ("EstimatedQReported", "%.4f"),
            ("Observations", "%d"),
            ("Errors", "%.2f"),
        ],
        rg_rows,
    )?;

    let qual_rows: Vec<Vec<String>> = tables
        .quality
        .iter()
        .map(|((rg, qual), datum)| {
            vec![
                rg_info.identifiers[*rg].clone(),
                qual.to_string(),
                "M".to_string(),
                format!("{:.4}", empirical_quality(datum)),
                datum.observations.to_string(),
                format!("{:.2}", datum.errors),
            ]
        })
        .collect();
    write_table(
        &mut writer,
        "RecalTable1",
        "",
        &[
            ("ReadGroup", "%s"),
            ("QualityScore", "%d"),
            ("EventType", "%s"),
            ("EmpiricalQuality", "%.4f"),
            ("Observations", "%d"),
            ("Errors", "%.2f"),
        ],
        qual_rows,
    )?;

    let mut cov_rows = Vec::with_capacity(tables.cycle.len() + tables.context.len());
    for ((rg, qual, cycle), datum) in &tables.cycle {
        cov_rows.push(vec![
            rg_info.identifiers[*rg].clone(),
            qual.to_string(),
            cycle.to_string(),
            "Cycle".to_string(),
            "M".to_string(),
            format!("{:.4}", empirical_quality(datum)),
            datum.observations.to_string(),
            format!("{:.2}", datum.errors),
        ]);
    }
    for ((rg, qual, context), datum) in &tables.context {
        cov_rows.push(vec![
            rg_info.identifiers[*rg].clone(),
            qual.to_string(),
            context.clone(),
            "Context".to_string(),
            "M".to_string(),
            format!("{:.4}", empirical_quality(datum)),
            datum.observations.to_string(),
            format!("{:.2}", datum.errors),
        ]);
    }
    cov_rows.sort_by(compare_recal_table2_rows);
    write_table(
        &mut writer,
        "RecalTable2",
        "",
        &[
            ("ReadGroup", "%s"),
            ("QualityScore", "%d"),
            ("CovariateValue", "%s"),
            ("CovariateName", "%s"),
            ("EventType", "%s"),
            ("EmpiricalQuality", "%.4f"),
            ("Observations", "%d"),
            ("Errors", "%.2f"),
        ],
        cov_rows,
    )?;
    writer.flush()?;
    Ok(())
}

fn argument_rows(args: &Args) -> Vec<Vec<String>> {
    vec![
        vec!["binary_tag_name".to_string(), "null".to_string()],
        vec![
            "covariate".to_string(),
            "ReadGroupCovariate,QualityScoreCovariate,ContextCovariate,CycleCovariate".to_string(),
        ],
        vec!["default_platform".to_string(), "null".to_string()],
        vec!["deletions_default_quality".to_string(), "45".to_string()],
        vec!["force_platform".to_string(), "null".to_string()],
        vec!["indels_context_size".to_string(), "3".to_string()],
        vec!["insertions_default_quality".to_string(), "45".to_string()],
        vec![
            "low_quality_tail".to_string(),
            args.low_quality_tail.to_string(),
        ],
        vec![
            "maximum_cycle_value".to_string(),
            args.max_cycle.to_string(),
        ],
        vec![
            "mismatches_context_size".to_string(),
            args.context_size.to_string(),
        ],
        vec!["mismatches_default_quality".to_string(), "-1".to_string()],
        vec!["no_standard_covs".to_string(), "false".to_string()],
        vec![
            "quantizing_levels".to_string(),
            args.quantizing_levels.to_string(),
        ],
        vec!["recalibration_report".to_string(), "null".to_string()],
        vec!["run_without_dbsnp".to_string(), "false".to_string()],
        vec![
            "solid_nocall_strategy".to_string(),
            "THROW_EXCEPTION".to_string(),
        ],
        vec!["solid_recal_mode".to_string(), "SET_Q_ZERO".to_string()],
    ]
}

fn compare_recal_table2_rows(left: &Vec<String>, right: &Vec<String>) -> Ordering {
    left[0]
        .cmp(&right[0])
        .then_with(|| numeric_string_cmp(&left[1], &right[1]))
        .then_with(|| left[2].cmp(&right[2]))
        .then_with(|| left[3].cmp(&right[3]))
}

fn numeric_string_cmp(left: &str, right: &str) -> Ordering {
    match (left.parse::<i32>(), right.parse::<i32>()) {
        (Ok(l), Ok(r)) => l.cmp(&r),
        _ => left.cmp(right),
    }
}

fn write_table<W: Write>(
    writer: &mut W,
    name: &str,
    description: &str,
    columns: &[(&str, &str)],
    rows: Vec<Vec<String>>,
) -> std::io::Result<()> {
    writeln!(
        writer,
        "#:GATKTable:{}:{}:{}:;",
        columns.len(),
        rows.len(),
        columns
            .iter()
            .map(|(_, format)| *format)
            .collect::<Vec<&str>>()
            .join(":")
    )?;
    writeln!(writer, "#:GATKTable:{name}:{description}")?;

    let mut widths: Vec<usize> = columns.iter().map(|(column, _)| column.len()).collect();
    for row in &rows {
        for (idx, value) in row.iter().enumerate() {
            widths[idx] = widths[idx].max(value.len());
        }
    }
    write_fixed_row(
        writer,
        &columns
            .iter()
            .map(|(column, _)| column.to_string())
            .collect::<Vec<String>>(),
        &widths,
        true,
    )?;
    for row in &rows {
        write_fixed_row(writer, row, &widths, false)?;
    }
    writeln!(writer)?;
    Ok(())
}

fn write_fixed_row<W: Write>(
    writer: &mut W,
    row: &[String],
    widths: &[usize],
    header: bool,
) -> std::io::Result<()> {
    for (idx, value) in row.iter().enumerate() {
        if idx > 0 {
            write!(writer, "  ")?;
        }
        if header || !right_align(value) {
            write!(writer, "{:<width$}", value, width = widths[idx])?;
        } else {
            write!(writer, "{:>width$}", value, width = widths[idx])?;
        }
    }
    writeln!(writer)
}

fn right_align(value: &str) -> bool {
    value == "null" || value == "NA" || value.parse::<i64>().is_ok() || value.parse::<f64>().is_ok()
}

