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

