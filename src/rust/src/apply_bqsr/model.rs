#[derive(Debug)]
struct RecalibrationModel {
    args: RecalibrationArgs,
    read_groups: HashMap<String, usize>,
    read_group_table: Vec<Option<RecalDatum>>,
    quality_score_table: HashMap<(usize, u8), RecalDatum>,
    cycle_table: HashMap<(usize, u8, i32), RecalDatum>,
    context_table: HashMap<(usize, u8, u32), RecalDatum>,
    quantized_quals: Vec<u8>,
}

#[derive(Debug, Clone, Copy)]
struct RecalibrationArgs {
    preserve_qscores_less_than: u8,
    mismatches_context_size: usize,
    low_quality_tail: u8,
    maximum_cycle_value: i32,
}

impl Default for RecalibrationArgs {
    fn default() -> Self {
        Self {
            preserve_qscores_less_than: DEFAULT_PRESERVE_QSCORES_LESS_THAN,
            mismatches_context_size: DEFAULT_MISMATCHES_CONTEXT_SIZE,
            low_quality_tail: DEFAULT_LOW_QUALITY_TAIL,
            maximum_cycle_value: DEFAULT_MAXIMUM_CYCLE_VALUE,
        }
    }
}

impl RecalibrationModel {
    fn from_path(path: &Path, use_report_quantization: bool) -> Result<Self> {
        let text = fs::read_to_string(path)
            .with_context(|| format!("reading recalibration table {}", path.display()))?;
        Self::from_str(&text, use_report_quantization)
    }

    fn from_str(text: &str, use_report_quantization: bool) -> Result<Self> {
        let mut parser = RecalibrationTableParser::new(use_report_quantization);
        parser.parse(text)?;
        parser.finish()
    }

    fn recalibrate_qualities(
        &self,
        read_group_identifier: &str,
        bases: &[u8],
        qualities: &[u8],
        flags: u16,
        allow_missing_read_groups: bool,
    ) -> Result<Vec<u8>> {
        let Some(&rg_key) = self.read_groups.get(read_group_identifier) else {
            if allow_missing_read_groups {
                return Ok(qualities
                    .iter()
                    .map(|&q| self.quantized_quals[q as usize])
                    .collect());
            }
            bail!("read group {read_group_identifier} not found in recalibration table");
        };

        let read_group_datum = self
            .read_group_table
            .get(rg_key)
            .and_then(Option::as_ref)
            .with_context(|| {
                format!("read group {read_group_identifier} has no substitution row")
            })?;
        let prior_quality = read_group_datum.reported_quality;
        let empirical_quality_for_read_group = read_group_datum.empirical_quality(prior_quality);
        let contexts = context_keys_for_read(
            bases,
            qualities,
            flags,
            self.args.mismatches_context_size,
            self.args.low_quality_tail,
        );

        let mut recalibrated = qualities.to_vec();
        for offset in 0..qualities.len() {
            let reported_quality = qualities[offset];
            if reported_quality < self.args.preserve_qscores_less_than {
                continue;
            }

            let posterior_quality = self
                .quality_score_table
                .get(&(rg_key, reported_quality))
                .map(|datum| datum.empirical_quality(empirical_quality_for_read_group))
                .unwrap_or(empirical_quality_for_read_group);

            let mut delta_special_covariates = 0.0;
            let cycle = cycle_for_offset(
                offset,
                qualities.len(),
                flags,
                self.args.maximum_cycle_value,
            )?;
            if let Some(datum) = self.cycle_table.get(&(rg_key, reported_quality, cycle)) {
                delta_special_covariates +=
                    datum.empirical_quality(posterior_quality) - posterior_quality;
            }
            if let Some(context_key) = contexts[offset] {
                if let Some(datum) =
                    self.context_table
                        .get(&(rg_key, reported_quality, context_key))
                {
                    delta_special_covariates +=
                        datum.empirical_quality(posterior_quality) - posterior_quality;
                }
            }

            let raw_quality = posterior_quality + delta_special_covariates;
            let bounded_quality = bound_qual(fast_round(raw_quality), MAX_RECALIBRATED_Q_SCORE);
            recalibrated[offset] = self.quantized_quals[bounded_quality as usize];
        }
        Ok(recalibrated)
    }
}

struct RecalibrationTableParser {
    args: RecalibrationArgs,
    read_groups: HashMap<String, usize>,
    read_group_table: Vec<Option<RecalDatum>>,
    quality_score_table: HashMap<(usize, u8), RecalDatum>,
    cycle_table: HashMap<(usize, u8, i32), RecalDatum>,
    context_table: HashMap<(usize, u8, u32), RecalDatum>,
    report_quantized_quals: Vec<Option<u8>>,
    use_report_quantization: bool,
}

impl RecalibrationTableParser {
    fn new(use_report_quantization: bool) -> Self {
        Self {
            args: RecalibrationArgs::default(),
            read_groups: HashMap::new(),
            read_group_table: Vec::new(),
            quality_score_table: HashMap::new(),
            cycle_table: HashMap::new(),
            context_table: HashMap::new(),
            report_quantized_quals: vec![None; (MAX_RECALIBRATED_Q_SCORE + 1) as usize],
            use_report_quantization,
        }
    }

    fn parse(&mut self, text: &str) -> Result<()> {
        let mut current_table: Option<&str> = None;
        let mut expect_header = false;

        for raw_line in text.lines() {
            let line = raw_line.trim_end();
            if line.is_empty() {
                continue;
            }
            if let Some(rest) = line.strip_prefix("#:GATKTable:") {
                let table_name = rest.split(':').next().unwrap_or_default();
                current_table = match table_name {
                    "Arguments" | "Quantized" | "RecalTable0" | "RecalTable1" | "RecalTable2" => {
                        Some(table_name)
                    }
                    _ => None,
                };
                expect_header = current_table.is_some();
                continue;
            }
            if line.starts_with('#') {
                continue;
            }
            let Some(table) = current_table else {
                continue;
            };
            if expect_header {
                expect_header = false;
                continue;
            }

            match table {
                "Arguments" => self.parse_argument_row(line)?,
                "Quantized" => self.parse_quantized_row(line)?,
                "RecalTable0" => self.parse_read_group_row(line)?,
                "RecalTable1" => self.parse_quality_score_row(line)?,
                "RecalTable2" => self.parse_covariate_row(line)?,
                _ => {}
            }
        }

        Ok(())
    }

    fn finish(self) -> Result<RecalibrationModel> {
        if self.read_groups.is_empty() {
            bail!("recalibration table has no read groups");
        }

        let quantized_quals = if self.use_report_quantization {
            self.report_quantized_quals
                .into_iter()
                .enumerate()
                .map(|(q, value)| value.with_context(|| format!("missing quantized row for Q{q}")))
                .collect::<Result<Vec<_>>>()?
        } else {
            (0..=MAX_RECALIBRATED_Q_SCORE).map(|q| q as u8).collect()
        };

        Ok(RecalibrationModel {
            args: self.args,
            read_groups: self.read_groups,
            read_group_table: self.read_group_table,
            quality_score_table: self.quality_score_table,
            cycle_table: self.cycle_table,
            context_table: self.context_table,
            quantized_quals,
        })
    }

    fn parse_argument_row(&mut self, line: &str) -> Result<()> {
        let mut fields = line.split_whitespace();
        let Some(argument) = fields.next() else {
            return Ok(());
        };
        let value = fields.next().unwrap_or("null");
        match argument {
            "covariate" => {
                let expected =
                    "ReadGroupCovariate,QualityScoreCovariate,ContextCovariate,CycleCovariate";
                if value != expected {
                    bail!(
                        "unsupported covariate list {value}; only standard BQSR covariates are supported"
                    );
                }
            }
            "no_standard_covs" if value == "true" => {
                bail!("recalibration tables with no_standard_covs=true are not supported");
            }
            "mismatches_context_size" => {
                self.args.mismatches_context_size = parse_field(value, argument)?;
            }
            "low_quality_tail" => {
                self.args.low_quality_tail = parse_field(value, argument)?;
            }
            "maximum_cycle_value" => {
                self.args.maximum_cycle_value = parse_field(value, argument)?;
            }
            _ => {}
        }
        Ok(())
    }

    fn parse_quantized_row(&mut self, line: &str) -> Result<()> {
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 3 {
            bail!("malformed Quantized row: {line}");
        }
        let quality: usize = parse_field(fields[0], "QualityScore")?;
        let quantized: u8 = parse_field(fields[2], "QuantizedScore")?;
        if quality < self.report_quantized_quals.len() {
            self.report_quantized_quals[quality] = Some(quantized);
        }
        Ok(())
    }

    fn parse_read_group_row(&mut self, line: &str) -> Result<()> {
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 6 {
            bail!("malformed RecalTable0 row: {line}");
        }
        if fields[1] != BASE_SUBSTITUTION_EVENT {
            return Ok(());
        }
        let rg = self.read_group_key(fields[0]);
        let reported_quality = parse_field(fields[3], "EstimatedQReported")?;
        let observations = parse_field(fields[4], "Observations")?;
        let errors = parse_field(fields[5], "Errors")?;
        self.read_group_table[rg] = Some(RecalDatum::new(observations, errors, reported_quality));
        Ok(())
    }

    fn parse_quality_score_row(&mut self, line: &str) -> Result<()> {
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 6 {
            bail!("malformed RecalTable1 row: {line}");
        }
        if fields[2] != BASE_SUBSTITUTION_EVENT {
            return Ok(());
        }
        let rg = self.read_group_key(fields[0]);
        let quality = parse_field(fields[1], "QualityScore")?;
        let observations = parse_field(fields[4], "Observations")?;
        let errors = parse_field(fields[5], "Errors")?;
        self.quality_score_table.insert(
            (rg, quality),
            RecalDatum::new(observations, errors, f64::from(quality)),
        );
        Ok(())
    }

    fn parse_covariate_row(&mut self, line: &str) -> Result<()> {
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 8 {
            bail!("malformed RecalTable2 row: {line}");
        }
        if fields[4] != BASE_SUBSTITUTION_EVENT {
            return Ok(());
        }
        let rg = self.read_group_key(fields[0]);
        let quality = parse_field(fields[1], "QualityScore")?;
        let observations = parse_field(fields[6], "Observations")?;
        let errors = parse_field(fields[7], "Errors")?;
        let datum = RecalDatum::new(observations, errors, f64::from(quality));
        match fields[3] {
            "Cycle" => {
                let cycle = parse_field(fields[2], "Cycle")?;
                self.cycle_table.insert((rg, quality, cycle), datum);
            }
            "Context" => {
                let context_key = key_from_context(fields[2].as_bytes())
                    .with_context(|| format!("invalid Context covariate {}", fields[2]))?;
                self.context_table.insert((rg, quality, context_key), datum);
            }
            covariate => bail!("unsupported covariate {covariate} in RecalTable2"),
        }
        Ok(())
    }

    fn read_group_key(&mut self, read_group: &str) -> usize {
        if let Some(&key) = self.read_groups.get(read_group) {
            return key;
        }
        let key = self.read_groups.len();
        self.read_groups.insert(read_group.to_owned(), key);
        self.read_group_table.push(None);
        key
    }
}

fn parse_field<T>(value: &str, name: &str) -> Result<T>
where
    T: std::str::FromStr,
    T::Err: std::error::Error + Send + Sync + 'static,
{
    value
        .parse()
        .with_context(|| format!("parsing {name} value {value}"))
}

#[derive(Debug, Clone)]
struct RecalDatum {
    observations: u64,
    errors: f64,
    reported_quality: f64,
    empirical_quality: OnceLock<i32>,
}

impl RecalDatum {
    fn new(observations: u64, errors: f64, reported_quality: f64) -> Self {
        Self {
            observations,
            errors,
            reported_quality,
            empirical_quality: OnceLock::new(),
        }
    }

    fn empirical_quality(&self, prior_quality_score: f64) -> f64 {
        f64::from(*self.empirical_quality.get_or_init(|| {
            let mut observations = self.observations + 2;
            let mut errors = (self.errors + 0.5) as u64 + 1;
            if observations > MAX_NUMBER_OF_OBSERVATIONS {
                let fraction = MAX_NUMBER_OF_OBSERVATIONS as f64 / observations as f64;
                errors = (errors as f64 * fraction).round() as u64;
                observations = MAX_NUMBER_OF_OBSERVATIONS;
            }

            bayesian_estimate_of_empirical_quality(observations, errors, prior_quality_score)
                .min(MAX_RECALIBRATED_Q_SCORE)
        }))
    }
}

fn bayesian_estimate_of_empirical_quality(
    observations: u64,
    errors: u64,
    prior_quality_score: f64,
) -> i32 {
    let mut best_quality = 0;
    let mut best_log_posterior = f64::NEG_INFINITY;
    for quality in 0..=MAX_REASONABLE_Q_SCORE {
        let log_posterior = log_prior(quality, prior_quality_score)
            + log_binomial_likelihood(quality, observations, errors);
        if log_posterior > best_log_posterior {
            best_log_posterior = log_posterior;
            best_quality = quality;
        }
    }
    best_quality
}

fn log_prior(quality: i32, prior_quality_score: f64) -> f64 {
    let difference = ((f64::from(quality) - prior_quality_score) as i32)
        .abs()
        .min(40);
    let sigma = 0.5;
    -0.5 * (f64::from(difference) / sigma).powi(2)
}

fn log_binomial_likelihood(quality: i32, observations: u64, errors: u64) -> f64 {
    if observations == 0 {
        return 0.0;
    }
    let p_error = 10_f64.powf(f64::from(quality) / -10.0);
    if p_error == 1.0 {
        return if errors == observations {
            0.0
        } else {
            f64::NEG_INFINITY
        };
    }
    if p_error == 0.0 {
        return if errors == 0 { 0.0 } else { f64::NEG_INFINITY };
    }

    errors as f64 * p_error.ln() + (observations - errors) as f64 * (-p_error).ln_1p()
}

fn fast_round(value: f64) -> i32 {
    if value > 0.0 {
        (value + 0.5) as i32
    } else {
        (value - 0.5) as i32
    }
}

fn bound_qual(quality: i32, max_quality: i32) -> i32 {
    quality.clamp(1, max_quality)
}

fn context_keys_for_read(
    bases: &[u8],
    qualities: &[u8],
    flags: u16,
    context_size: usize,
    low_quality_tail: u8,
) -> Vec<Option<u32>> {
    let read_len = bases.len();
    let mut contexts = vec![None; read_len];
    if read_len == 0 || context_size == 0 {
        return contexts;
    }

    let Some(mut stranded_bases) = stranded_low_quality_tail_clipped_bases(
        bases,
        qualities,
        is_reverse(flags),
        low_quality_tail,
    ) else {
        return contexts;
    };

    if stranded_bases.len() != read_len {
        return contexts;
    }

    for stranded_offset in 0..read_len {
        let read_offset = if is_reverse(flags) {
            read_len - stranded_offset - 1
        } else {
            stranded_offset
        };
        if stranded_offset + 1 < context_size {
            continue;
        }
        let start = stranded_offset + 1 - context_size;
        contexts[read_offset] = key_from_context(&stranded_bases[start..=stranded_offset]);
    }

    stranded_bases.clear();
    contexts
}

fn stranded_low_quality_tail_clipped_bases(
    bases: &[u8],
    qualities: &[u8],
    reverse: bool,
    low_quality_tail: u8,
) -> Option<Vec<u8>> {
    let read_len = bases.len();
    let mut right = read_len;
    while right > 0 && qualities[right - 1] <= low_quality_tail {
        right -= 1;
    }
    let mut left = 0;
    while left < read_len && qualities[left] <= low_quality_tail {
        left += 1;
    }
    if left >= right {
        return None;
    }

    let mut clipped = bases.to_vec();
    for base in &mut clipped[..left] {
        *base = b'N';
    }
    for base in &mut clipped[right..] {
        *base = b'N';
    }

    if reverse {
        clipped = clipped.into_iter().rev().map(simple_complement).collect();
    }
    Some(clipped)
}

fn cycle_for_offset(
    offset: usize,
    read_len: usize,
    flags: u16,
    maximum_cycle_value: i32,
) -> Result<i32> {
    let read_order_factor = if is_paired(flags) && is_second_in_pair(flags) {
        -1
    } else {
        1
    };
    let (cycle_start, increment) = if is_reverse(flags) {
        (read_len as i32 * read_order_factor, -read_order_factor)
    } else {
        (read_order_factor, read_order_factor)
    };
    let cycle = cycle_start + offset as i32 * increment;
    if cycle.abs() > maximum_cycle_value {
        bail!(
            "cycle {} exceeds maximum_cycle_value {}",
            cycle,
            maximum_cycle_value
        );
    }
    Ok(cycle)
}

fn is_paired(flags: u16) -> bool {
    flags & 0x1 != 0
}

fn is_reverse(flags: u16) -> bool {
    flags & 0x10 != 0
}

fn is_second_in_pair(flags: u16) -> bool {
    flags & 0x80 != 0
}

fn key_from_context(context: &[u8]) -> Option<u32> {
    let mut key = context.len() as u32;
    let mut bit_offset = 4;
    for &base in context {
        let base_index = simple_base_to_index(base)?;
        key |= base_index << bit_offset;
        bit_offset += 2;
    }
    Some(key)
}

fn simple_base_to_index(base: u8) -> Option<u32> {
    match base {
        b'A' | b'a' | b'*' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

fn simple_complement(base: u8) -> u8 {
    match base {
        b'A' | b'a' => b'T',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        b'T' | b't' => b'A',
        _ => base,
    }
}
