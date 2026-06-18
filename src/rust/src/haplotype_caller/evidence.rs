fn pileup_snp_evidence(
    pileup: &bam::pileup::Pileup,
    min_baseq: u8,
    min_mapq: u8,
    min_tail_quality: u8,
    exclude_supplementary: bool,
    dont_use_soft_clipped_bases: bool,
) -> SnpEvidence {
    pileup_snp_evidence_with_rows(
        pileup,
        min_baseq,
        min_mapq,
        min_tail_quality,
        exclude_supplementary,
        dont_use_soft_clipped_bases,
        None,
    )
    .0
}

fn pileup_snp_evidence_with_rows(
    pileup: &bam::pileup::Pileup,
    min_baseq: u8,
    min_mapq: u8,
    min_tail_quality: u8,
    exclude_supplementary: bool,
    dont_use_soft_clipped_bases: bool,
    row_context: Option<&ReplayRowContext<'_>>,
) -> (SnpEvidence, Vec<ReplayReadObservationRow>) {
    let mut observations_by_fragment: HashMap<Vec<u8>, Vec<NamedSnpObservation>> = HashMap::new();
    for alignment in pileup.alignments() {
        let record = alignment.record();
        if !read_passes_hc_filter(&record, min_mapq, exclude_supplementary)
            || alignment.is_refskip()
        {
            continue;
        }
        let Some(qpos) = alignment.qpos() else {
            continue;
        };
        if record
            .qual()
            .get(qpos)
            .is_none_or(|quality| *quality < min_baseq)
        {
            continue;
        }

        // Apply GATK-like read clipping (hard-clip soft-clips, low-quality tails).
        if let Some((clip_start, clip_end)) =
            clip_read_for_evidence(&record, min_tail_quality, dont_use_soft_clipped_bases)
        {
            if qpos < clip_start || qpos >= clip_end {
                continue;
            }
        } else {
            continue;
        }

        let base = normalize_base(record.seq()[qpos]);
        observations_by_fragment
            .entry(record.qname().to_vec())
            .or_default()
            .push(NamedSnpObservation {
                read_name: record.qname().to_vec(),
                qpos,
                mapq: record.mapq(),
                base_index: base_index(base),
                quality: record.qual()[qpos],
                is_reverse: record.is_reverse(),
            });
    }

    let mut evidence = SnpEvidence::default();
    let mut rows = Vec::new();
    for observations in observations_by_fragment.into_values() {
        for named in adjust_named_snp_observations(&observations) {
            if named.quality < min_baseq {
                continue;
            }
            evidence.active_observations.push(ActiveBaseObservation {
                base_index: named.base_index,
                quality: named.quality,
            });
            if let Some(index) = named.base_index {
                let observation = BaseObservation {
                    base_index: index,
                    quality: named.quality,
                    is_reverse: named.is_reverse,
                };
                evidence.counts.counts[index] += 1;
                evidence.counts.depth += 1;
                evidence.strands[index].increment(named.is_reverse);
                evidence.observations.push(observation);
                if let Some(context) = row_context {
                    rows.push(ReplayReadObservationRow {
                        region: context.region.to_string(),
                        read: String::from_utf8_lossy(&named.read_name).into_owned(),
                        kind: "snp",
                        pos: context.pos,
                        qpos: named.qpos,
                        allele: (base_from_index(index) as char).to_string(),
                        adjusted_quality: named.quality,
                        mapq: named.mapq,
                        strand: strand_label(named.is_reverse),
                    });
                }
            }
        }
    }
    (evidence, rows)
}

#[cfg(test)]
fn adjust_fragment_base_observations(observations: &[BaseObservation]) -> Vec<BaseObservation> {
    if observations.len() <= 1 {
        return observations.to_vec();
    }

    let first_base_index = observations[0].base_index;
    if observations
        .iter()
        .all(|observation| observation.base_index == first_base_index)
    {
        observations
            .iter()
            .map(|observation| BaseObservation {
                quality: observation.quality.min(HALF_DEFAULT_PCR_SNV_QUAL),
                ..*observation
            })
            .collect()
    } else {
        observations
            .iter()
            .map(|observation| BaseObservation {
                quality: 0,
                ..*observation
            })
            .collect()
    }
}

fn adjust_named_snp_observations(observations: &[NamedSnpObservation]) -> Vec<NamedSnpObservation> {
    if observations.len() <= 1 {
        return observations.to_vec();
    }

    let first_base_index = observations[0].base_index;
    let all_same_base = observations
        .iter()
        .all(|observation| observation.base_index == first_base_index);
    observations
        .iter()
        .map(|observation| {
            let mut adjusted = observation.clone();
            adjusted.quality = if all_same_base {
                adjusted.quality.min(HALF_DEFAULT_PCR_SNV_QUAL)
            } else {
                0
            };
            adjusted
        })
        .collect()
}

fn pileup_indel_evidence(
    pileup: &bam::pileup::Pileup,
    min_baseq: u8,
    min_mapq: u8,
    min_tail_quality: u8,
    exclude_supplementary: bool,
    dont_use_soft_clipped_bases: bool,
) -> IndelEvidence {
    pileup_indel_evidence_with_rows(
        pileup,
        min_baseq,
        min_mapq,
        min_tail_quality,
        exclude_supplementary,
        dont_use_soft_clipped_bases,
        None,
    )
    .0
}

fn pileup_indel_evidence_with_rows(
    pileup: &bam::pileup::Pileup,
    min_baseq: u8,
    min_mapq: u8,
    min_tail_quality: u8,
    exclude_supplementary: bool,
    dont_use_soft_clipped_bases: bool,
    row_context: Option<&ReplayRowContext<'_>>,
) -> (IndelEvidence, Vec<ReplayReadObservationRow>) {
    let mut observations_by_fragment: HashMap<Vec<u8>, Vec<NamedIndelObservation>> = HashMap::new();
    for alignment in pileup.alignments() {
        let record = alignment.record();
        if !read_passes_hc_filter(&record, min_mapq, exclude_supplementary)
            || alignment.is_refskip()
        {
            continue;
        }
        let Some(qpos) = alignment.qpos() else {
            continue;
        };
        if record
            .qual()
            .get(qpos)
            .is_none_or(|quality| *quality < min_baseq)
        {
            continue;
        }

        // Apply GATK-like read clipping.
        if let Some((clip_start, clip_end)) =
            clip_read_for_evidence(&record, min_tail_quality, dont_use_soft_clipped_bases)
        {
            if qpos < clip_start || qpos >= clip_end {
                continue;
            }
        } else {
            continue;
        }

        let allele = match alignment.indel() {
            Indel::None => Some(IndelObservationAllele::Ref),
            Indel::Ins(len) if len <= MAX_BOOTSTRAP_INDEL_LEN => {
                let inserted = inserted_bases(&record, qpos, len);
                if inserted.is_empty() {
                    None
                } else {
                    Some(IndelObservationAllele::Alt(IndelAllele::Insertion(
                        inserted,
                    )))
                }
            }
            Indel::Del(len) if len <= MAX_BOOTSTRAP_INDEL_LEN => {
                Some(IndelObservationAllele::Alt(IndelAllele::Deletion(len)))
            }
            _ => None,
        };
        let Some(allele) = allele else {
            continue;
        };
        observations_by_fragment
            .entry(record.qname().to_vec())
            .or_default()
            .push(NamedIndelObservation {
                read_name: record.qname().to_vec(),
                qpos,
                mapq: record.mapq(),
                observation: IndelObservation {
                    allele,
                    quality: indel_observation_quality(record.qual()[qpos]),
                    is_reverse: record.is_reverse(),
                },
            });
    }

    let mut evidence = IndelEvidence::default();
    let mut rows = Vec::new();
    for observations in observations_by_fragment.into_values() {
        for named in adjust_named_indel_observations(&observations) {
            let observation = named.observation;
            if observation.quality < min_baseq {
                continue;
            }
            evidence.counts.depth += 1;
            match &observation.allele {
                IndelObservationAllele::Ref => {
                    evidence.counts.ref_count += 1;
                    evidence.ref_strand.increment(observation.is_reverse);
                }
                IndelObservationAllele::Alt(allele) => {
                    *evidence.counts.counts.entry(allele.clone()).or_insert(0) += 1;
                    evidence
                        .alt_strands
                        .entry(allele.clone())
                        .or_default()
                        .increment(observation.is_reverse);
                }
            }
            if let Some(context) = row_context {
                rows.push(ReplayReadObservationRow {
                    region: context.region.to_string(),
                    read: String::from_utf8_lossy(&named.read_name).into_owned(),
                    kind: "indel",
                    pos: context.pos,
                    qpos: named.qpos,
                    allele: indel_observation_allele_label(&observation.allele),
                    adjusted_quality: observation.quality,
                    mapq: named.mapq,
                    strand: strand_label(observation.is_reverse),
                });
            }
            evidence.observations.push(observation);
        }
    }
    (evidence, rows)
}

fn adjust_named_indel_observations(
    observations: &[NamedIndelObservation],
) -> Vec<NamedIndelObservation> {
    observations.to_vec()
}

fn indel_observation_quality(base_quality: u8) -> u8 {
    base_quality.min(DEFAULT_INDEL_QUAL)
}

fn pair_hmm_base_quality(base_quality: u8, mapq: u8) -> u8 {
    let capped = base_quality.min(mapq);
    if capped < PAIR_HMM_BASE_QUALITY_SCORE_THRESHOLD {
        PAIR_HMM_MIN_USABLE_Q_SCORE
    } else {
        capped
    }
}

fn pair_hmm_indel_open_quality(quality: u8) -> u8 {
    quality.max(PAIR_HMM_MIN_USABLE_Q_SCORE)
}

fn filter_non_acgt_haplotypes_for_single_snp_region(
    local_haplotypes: &mut Vec<LocalHaplotype>,
    valid_events: &[VariantCall],
) {
    if valid_events.len() != 1 {
        return;
    }
    let event = &valid_events[0];
    if event.ref_allele.len() != 1 || event.alt_allele.len() != 1 {
        return;
    }

    let has_regular_alt = local_haplotypes
        .iter()
        .any(|hap| !hap.is_ref && hap.event_indices.contains(&0) && is_regular_bases(&hap.bases));
    if !has_regular_alt {
        return;
    }

    local_haplotypes.retain(|hap| hap.is_ref || is_regular_bases(&hap.bases));
}

fn genotype_assembled_events(
    local_haplotypes: &[LocalHaplotype],
    valid_events: &[VariantCall],
    read_haplotype_likelihoods: &[Vec<f64>],
    read_is_reverse_list: &[bool],
    read_ref_spans: &[(u64, u64)],
    min_confidence: f64,
) -> Vec<VariantCall> {
    if read_haplotype_likelihoods.is_empty() {
        return Vec::new();
    }

    let use_pair_genotyping = overlapping_event_mask(valid_events);
    let overlapping_event_indices = overlapping_event_indices(valid_events);
    let pair_context = use_pair_genotyping
        .iter()
        .copied()
        .any(|use_pair| use_pair)
        .then(|| build_pair_genotyping_context(local_haplotypes.len(), read_haplotype_likelihoods));
    let mut final_calls = Vec::new();

    for (event_idx, event) in valid_events.iter().enumerate() {
        let maybe_call = if use_pair_genotyping[event_idx] {
            genotype_overlapping_assembled_event(
                local_haplotypes,
                event_idx,
                event,
                read_haplotype_likelihoods,
                read_is_reverse_list,
                read_ref_spans,
                min_confidence,
                &overlapping_event_indices[event_idx],
                pair_context
                    .as_ref()
                    .expect("pair genotyping context must exist for overlapping events"),
            )
        } else {
            genotype_isolated_assembled_event(
                local_haplotypes,
                event_idx,
                event,
                read_haplotype_likelihoods,
                read_is_reverse_list,
                read_ref_spans,
                min_confidence,
            )
        };
        if let Some(final_call) = maybe_call {
            final_calls.push(final_call);
        }
    }

    final_calls
}

struct PairGenotypingContext {
    haplotype_pairs: Vec<(usize, usize)>,
    pair_log10_likelihoods: Vec<f64>,
    best_pair: (usize, usize),
}

fn build_pair_genotyping_context(
    n_haplotypes: usize,
    read_haplotype_likelihoods: &[Vec<f64>],
) -> PairGenotypingContext {
    let haplotype_pairs = enumerate_haplotype_pairs(n_haplotypes);
    let pair_log10_likelihoods =
        compute_pair_log10_likelihoods(&haplotype_pairs, read_haplotype_likelihoods);
    let best_pair_idx = pair_log10_likelihoods
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
        .map(|(idx, _)| idx)
        .unwrap_or(0);
    let best_pair = haplotype_pairs[best_pair_idx];

    PairGenotypingContext {
        haplotype_pairs,
        pair_log10_likelihoods,
        best_pair,
    }
}

fn genotype_isolated_assembled_event(
    local_haplotypes: &[LocalHaplotype],
    event_idx: usize,
    event: &VariantCall,
    read_haplotype_likelihoods: &[Vec<f64>],
    read_is_reverse_list: &[bool],
    read_ref_spans: &[(u64, u64)],
    min_confidence: f64,
) -> Option<VariantCall> {
    let hap_contains_allele = event_haplotype_mask(local_haplotypes, event_idx);
    let (ref_haps, alt_haps) = split_event_haplotype_indices(&hap_contains_allele);
    if alt_haps.is_empty() {
        return None;
    }

    let (depth, ref_count, alt_count, fs) = count_event_support(
        read_haplotype_likelihoods,
        read_is_reverse_list,
        read_ref_spans,
        event,
        &ref_haps,
        &alt_haps,
    );

    let mut log10_likelihoods = [0.0; 3];
    for read_idx in 0..read_haplotype_likelihoods.len() {
        if !read_span_overlaps_event(read_ref_spans[read_idx], event) {
            continue;
        }
        let ref_log10 =
            marginalize_haplotype_indices(read_haplotype_likelihoods, read_idx, &ref_haps);
        let alt_log10 =
            marginalize_haplotype_indices(read_haplotype_likelihoods, read_idx, &alt_haps);
        log10_likelihoods[0] += ref_log10;
        log10_likelihoods[1] += log10_sum_exp(&[ref_log10, alt_log10]) - 0.3010299956639812;
        log10_likelihoods[2] += alt_log10;
    }

    finalize_assembled_event_call(
        event,
        variant_model_from_log10(log10_likelihoods, assembled_event_log10_priors()),
        depth,
        ref_count,
        alt_count,
        fs,
        min_confidence,
    )
}

fn genotype_overlapping_assembled_event(
    local_haplotypes: &[LocalHaplotype],
    event_idx: usize,
    event: &VariantCall,
    read_haplotype_likelihoods: &[Vec<f64>],
    read_is_reverse_list: &[bool],
    read_ref_spans: &[(u64, u64)],
    min_confidence: f64,
    overlapping_event_indices: &[usize],
    pair_context: &PairGenotypingContext,
) -> Option<VariantCall> {
    let hap_contains_allele = event_haplotype_mask(local_haplotypes, event_idx);
    if hap_contains_allele.iter().all(|present| !present) {
        return None;
    }

    let mut log10_likelihoods = [f64::NEG_INFINITY; 3];
    for (pair_idx, pair) in pair_context.haplotype_pairs.iter().copied().enumerate() {
        let genotype_index = pair_event_genotype_index(pair, &hap_contains_allele);
        if genotype_index == 0
            && pair_contains_competing_overlapping_event(
                pair,
                local_haplotypes,
                overlapping_event_indices,
            )
        {
            continue;
        }
        log10_likelihoods[genotype_index] = log10_sum_exp(&[
            log10_likelihoods[genotype_index],
            pair_context.pair_log10_likelihoods[pair_idx],
        ]);
    }

    let genotype_index = pair_event_genotype_index(pair_context.best_pair, &hap_contains_allele);
    if genotype_index == 0 {
        return None;
    }

    let (best_ref_haps, best_alt_haps) =
        best_pair_hap_sets(pair_context.best_pair, &hap_contains_allele);
    let (depth, ref_count, alt_count, fs) = count_event_support(
        read_haplotype_likelihoods,
        read_is_reverse_list,
        read_ref_spans,
        event,
        &best_ref_haps,
        &best_alt_haps,
    );

    let mut model = variant_model_from_log10(log10_likelihoods, assembled_event_log10_priors());
    model.genotype_index = genotype_index;
    normalize_variant_model_pl(&mut model);

    finalize_assembled_event_call(
        event,
        model,
        depth,
        ref_count,
        alt_count,
        fs,
        min_confidence,
    )
}

fn event_haplotype_mask(local_haplotypes: &[LocalHaplotype], event_idx: usize) -> Vec<bool> {
    local_haplotypes
        .iter()
        .map(|hap| hap.event_indices.contains(&event_idx))
        .collect()
}

fn split_event_haplotype_indices(hap_contains_allele: &[bool]) -> (Vec<usize>, Vec<usize>) {
    let mut ref_haps = Vec::new();
    let mut alt_haps = Vec::new();
    for (hap_idx, contains_allele) in hap_contains_allele.iter().copied().enumerate() {
        if contains_allele {
            alt_haps.push(hap_idx);
        } else {
            ref_haps.push(hap_idx);
        }
    }
    (ref_haps, alt_haps)
}

fn count_event_support(
    read_haplotype_likelihoods: &[Vec<f64>],
    read_is_reverse_list: &[bool],
    read_ref_spans: &[(u64, u64)],
    event: &VariantCall,
    ref_haps: &[usize],
    alt_haps: &[usize],
) -> (u32, u32, u32, f64) {
    let mut ref_count = 0_u32;
    let mut alt_count = 0_u32;
    let mut ref_strand = StrandCounts::default();
    let mut alt_strand = StrandCounts::default();

    for read_idx in 0..read_haplotype_likelihoods.len() {
        if !read_span_overlaps_event(read_ref_spans[read_idx], event) {
            continue;
        }
        let ref_log10 =
            marginalize_haplotype_indices(read_haplotype_likelihoods, read_idx, ref_haps);
        let alt_log10 =
            marginalize_haplotype_indices(read_haplotype_likelihoods, read_idx, alt_haps);
        let is_reverse = read_is_reverse_list[read_idx];
        if ref_log10 - alt_log10 > 0.2 {
            ref_count += 1;
            ref_strand.increment(is_reverse);
        } else if alt_log10 - ref_log10 > 0.2 {
            alt_count += 1;
            alt_strand.increment(is_reverse);
        }
    }

    let depth = ref_count + alt_count;
    let fs = fisher_strand_score(ref_strand, alt_strand);
    (depth, ref_count, alt_count, fs)
}

fn finalize_assembled_event_call(
    event: &VariantCall,
    mut model: VariantModel,
    depth: u32,
    ref_count: u32,
    alt_count: u32,
    fs: f64,
    min_confidence: f64,
) -> Option<VariantCall> {
    normalize_variant_model_pl(&mut model);
    if model.genotype_index == 0 || f64::from(model.qual) < min_confidence {
        return None;
    }

    let mut final_call = event.clone();
    final_call.pl = model.pl;
    final_call.genotype_index = model.genotype_index;
    final_call.qual = model.qual;
    final_call.depth = depth;
    final_call.ref_count = ref_count;
    final_call.alt_count = alt_count;
    final_call.fs = fs;
    Some(final_call)
}

fn normalize_variant_model_pl(model: &mut VariantModel) {
    if model.pl[0] == 0 && model.pl[1] == 0 && model.pl[2] == 0 {
        for idx in 0..3 {
            if idx != model.genotype_index {
                model.pl[idx] = 9999;
            }
        }
    }
}

fn assembled_event_log10_priors() -> [f64; 3] {
    let heterozygosity: f64 = 1e-3;
    [
        (1.0 - 1.5 * heterozygosity).log10(),
        heterozygosity.log10(),
        (0.5 * heterozygosity).log10(),
    ]
}

fn overlapping_event_mask(valid_events: &[VariantCall]) -> Vec<bool> {
    let mut mask = vec![false; valid_events.len()];
    for left_idx in 0..valid_events.len() {
        for right_idx in left_idx + 1..valid_events.len() {
            if events_overlap(&valid_events[left_idx], &valid_events[right_idx]) {
                mask[left_idx] = true;
                mask[right_idx] = true;
            }
        }
    }
    mask
}

fn overlapping_event_indices(valid_events: &[VariantCall]) -> Vec<Vec<usize>> {
    let mut overlaps = vec![Vec::new(); valid_events.len()];
    for left_idx in 0..valid_events.len() {
        for right_idx in left_idx + 1..valid_events.len() {
            if events_overlap(&valid_events[left_idx], &valid_events[right_idx]) {
                overlaps[left_idx].push(right_idx);
                overlaps[right_idx].push(left_idx);
            }
        }
    }
    overlaps
}

fn events_overlap(left: &VariantCall, right: &VariantCall) -> bool {
    let left_end = left.pos + left.ref_allele.len() as u64 - 1;
    let right_end = right.pos + right.ref_allele.len() as u64 - 1;
    left.pos <= right_end && right.pos <= left_end
}

fn pair_contains_competing_overlapping_event(
    pair: (usize, usize),
    local_haplotypes: &[LocalHaplotype],
    overlapping_event_indices: &[usize],
) -> bool {
    if overlapping_event_indices.is_empty() {
        return false;
    }
    [pair.0, pair.1].into_iter().any(|hap_idx| {
        local_haplotypes[hap_idx]
            .event_indices
            .iter()
            .any(|event_idx| overlapping_event_indices.contains(event_idx))
    })
}

fn read_reference_span(record: &bam::Record) -> (u64, u64) {
    let start = (record.pos() + 1) as u64;
    let end = record.reference_end() as u64;
    (start, end)
}

fn read_reference_span_from_start_and_cigar(start: u64, cigar: &str) -> (u64, u64) {
    let mut reference_len = 0_u64;
    let mut current_len = 0_u64;
    for byte in cigar.bytes() {
        if byte.is_ascii_digit() {
            current_len = current_len
                .saturating_mul(10)
                .saturating_add(u64::from(byte - b'0'));
            continue;
        }
        match byte as char {
            'M' | 'D' | 'N' | '=' | 'X' => {
                reference_len = reference_len.saturating_add(current_len);
            }
            'I' | 'S' | 'H' | 'P' => {}
            _ => {}
        }
        current_len = 0;
    }
    let end = start.saturating_add(reference_len.saturating_sub(1));
    (start, end)
}

fn read_span_overlaps_event(read_span: (u64, u64), event: &VariantCall) -> bool {
    let event_end = event.pos + event.ref_allele.len() as u64 - 1;
    read_span.0 <= event_end && event.pos <= read_span.1
}

fn enumerate_haplotype_pairs(n_haplotypes: usize) -> Vec<(usize, usize)> {
    let mut pairs = Vec::new();
    for left in 0..n_haplotypes {
        for right in left..n_haplotypes {
            pairs.push((left, right));
        }
    }
    pairs
}

fn compute_pair_log10_likelihoods(
    haplotype_pairs: &[(usize, usize)],
    read_haplotype_likelihoods: &[Vec<f64>],
) -> Vec<f64> {
    let mut pair_log10_likelihoods = vec![0.0; haplotype_pairs.len()];
    for (pair_idx, (left, right)) in haplotype_pairs.iter().copied().enumerate() {
        let mut total = 0.0;
        for read_likelihoods in read_haplotype_likelihoods {
            total += if left == right {
                read_likelihoods[left]
            } else {
                log10_sum_exp(&[read_likelihoods[left], read_likelihoods[right]])
                    - 0.3010299956639812
            };
        }
        pair_log10_likelihoods[pair_idx] = total;
    }
    pair_log10_likelihoods
}

fn pair_event_genotype_index(pair: (usize, usize), hap_contains_allele: &[bool]) -> usize {
    let left = hap_contains_allele[pair.0];
    let right = hap_contains_allele[pair.1];
    match (left, right) {
        (false, false) => 0,
        (true, true) => 2,
        _ => 1,
    }
}

fn best_pair_hap_sets(
    pair: (usize, usize),
    hap_contains_allele: &[bool],
) -> (Vec<usize>, Vec<usize>) {
    let mut ref_haps = Vec::new();
    let mut alt_haps = Vec::new();
    for hap_idx in [pair.0, pair.1] {
        if hap_contains_allele[hap_idx] {
            if !alt_haps.contains(&hap_idx) {
                alt_haps.push(hap_idx);
            }
        } else if !ref_haps.contains(&hap_idx) {
            ref_haps.push(hap_idx);
        }
    }
    (ref_haps, alt_haps)
}

fn marginalize_haplotype_indices(
    read_haplotype_likelihoods: &[Vec<f64>],
    read_idx: usize,
    hap_indices: &[usize],
) -> f64 {
    if hap_indices.is_empty() {
        return f64::NEG_INFINITY;
    }
    let mut values = Vec::with_capacity(hap_indices.len());
    for hap_idx in hap_indices {
        values.push(read_haplotype_likelihoods[read_idx][*hap_idx]);
    }
    marginalize_allele_likelihoods(&values)
}

fn best_snp_call(
    contig: &str,
    pos: u64,
    ref_base: u8,
    evidence: SnpEvidence,
    min_qual: f64,
) -> Option<VariantCall> {
    let ref_index = base_index(ref_base)?;
    let mut best_alt_index = None;
    let mut best_alt_count = 0_u32;
    for (index, count) in evidence.counts.counts.iter().copied().enumerate() {
        if index == ref_index || count <= best_alt_count {
            continue;
        }
        best_alt_index = Some(index);
        best_alt_count = count;
    }
    let alt_index = best_alt_index?;
    if !alt_support_passes(evidence.counts.depth, best_alt_count) {
        return None;
    }
    let model = snp_variant_model(&evidence.observations, ref_index, alt_index);
    if f64::from(model.qual) < min_qual {
        return None;
    }

    Some(VariantCall {
        contig: contig.to_string(),
        pos,
        id: None,
        db: false,
        ref_allele: vec![ref_base],
        alt_allele: vec![base_from_index(alt_index)],
        depth: evidence.counts.depth,
        ref_count: evidence.counts.counts[ref_index],
        alt_count: best_alt_count,
        qual: model.qual,
        fs: fisher_strand_score(evidence.strands[ref_index], evidence.strands[alt_index]),
        pl: model.pl,
        genotype_index: model.genotype_index,
    })
}

