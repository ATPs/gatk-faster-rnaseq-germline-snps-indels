fn best_snp_alt(ref_index: Option<usize>, evidence: &SnpEvidence) -> Option<(usize, u32)> {
    let ref_index = ref_index?;
    let mut best_alt_index = None;
    let mut best_alt_count = 0_u32;
    for (index, count) in evidence.counts.counts.iter().copied().enumerate() {
        if index == ref_index || count <= best_alt_count {
            continue;
        }
        best_alt_index = Some(index);
        best_alt_count = count;
    }
    best_alt_index.map(|index| (index, best_alt_count))
}

fn best_indel_call(
    contig: &str,
    pos: u64,
    ref_start: u64,
    ref_bases: &[u8],
    evidence: IndelEvidence,
    min_qual: f64,
) -> Option<VariantCall> {
    let (best_allele, best_alt_count) =
        evidence
            .counts
            .counts
            .iter()
            .max_by(|(allele_a, count_a), (allele_b, count_b)| {
                count_a.cmp(count_b).then_with(|| {
                    indel_allele_sort_key(allele_a).cmp(&indel_allele_sort_key(allele_b))
                })
            })?;
    if !alt_support_passes(evidence.counts.depth, *best_alt_count) {
        return None;
    }
    let model = indel_variant_model(&evidence.observations, best_allele);
    if f64::from(model.qual) < min_qual {
        return None;
    }

    let offset = usize::try_from(pos.checked_sub(ref_start)?).ok()?;
    let anchor = normalize_base(*ref_bases.get(offset)?);
    if !is_acgt(anchor) {
        return None;
    }
    let (ref_allele, alt_allele) = match best_allele {
        IndelAllele::Insertion(inserted) => {
            if inserted.iter().any(|base| !is_acgt(*base)) {
                return None;
            }
            let mut alt_allele = Vec::with_capacity(inserted.len() + 1);
            alt_allele.push(anchor);
            alt_allele.extend_from_slice(inserted);
            (vec![anchor], alt_allele)
        }
        IndelAllele::Deletion(len) => {
            let delete_len = usize::try_from(*len).ok()?;
            let end_offset = offset.checked_add(delete_len)?;
            let deleted = ref_bases.get(offset..=end_offset)?;
            let ref_allele: Vec<u8> = deleted.iter().map(|base| normalize_base(*base)).collect();
            if ref_allele.iter().any(|base| !is_acgt(*base)) {
                return None;
            }
            (ref_allele, vec![anchor])
        }
    };
    let (pos, ref_allele, alt_allele) =
        left_normalize_indel(pos, ref_start, ref_bases, ref_allele, alt_allele);

    Some(VariantCall {
        contig: contig.to_string(),
        pos,
        id: None,
        db: false,
        ref_allele,
        alt_allele,
        depth: evidence.counts.depth,
        ref_count: evidence.counts.ref_count,
        alt_count: *best_alt_count,
        qual: model.qual,
        fs: fisher_strand_score(
            evidence.ref_strand,
            evidence
                .alt_strands
                .get(best_allele)
                .copied()
                .unwrap_or_default(),
        ),
        pl: model.pl,
        genotype_index: model.genotype_index,
    })
}

fn best_indel_alt(evidence: &IndelEvidence) -> Option<(&IndelAllele, &u32)> {
    evidence
        .counts
        .counts
        .iter()
        .max_by(|(allele_a, count_a), (allele_b, count_b)| {
            count_a
                .cmp(count_b)
                .then_with(|| indel_allele_sort_key(allele_a).cmp(&indel_allele_sort_key(allele_b)))
        })
}

fn collect_pileup_fallback_events(
    contig: &str,
    region_start: u64,
    local_ref_bases: &[u8],
    active_loci: &[ReplayActiveLocusRow],
    min_qual: f64,
) -> Vec<VariantCall> {
    let region_end = region_start + local_ref_bases.len().saturating_sub(1) as u64;
    let mut events = Vec::new();
    for row in active_loci {
        if row.contig != contig || row.pos < region_start || row.pos > region_end {
            continue;
        }
        if row.depth < 10 {
            continue;
        }
        if row.snp_alt_count >= 3 {
            let alt_base = row.snp_best_alt.as_bytes().first().copied();
            let ref_base = base_index(row.ref_base);
            let alt_index = alt_base.and_then(base_index);
            let offset = (row.pos - region_start) as usize;
            if let (Some(ref_index), Some(alt_index)) = (ref_base, alt_index) {
                let mut observations = Vec::with_capacity(row.depth as usize);
                for _ in 0..row.depth.saturating_sub(row.snp_alt_count) {
                    observations.push(BaseObservation {
                        base_index: ref_index,
                        quality: 30,
                        is_reverse: false,
                    });
                }
                for _ in 0..row.snp_alt_count {
                    observations.push(BaseObservation {
                        base_index: alt_index,
                        quality: 30,
                        is_reverse: false,
                    });
                }
                let model = snp_variant_model(&observations, ref_index, alt_index);
                if f64::from(model.qual) >= min_qual {
                    let alt_base = base_from_index(alt_index);
                    if let Some(ref_byte) = local_ref_bases.get(offset).copied() {
                        events.push(VariantCall {
                            contig: contig.to_string(),
                            pos: row.pos,
                            id: None,
                            db: false,
                            ref_allele: vec![normalize_base(ref_byte)],
                            alt_allele: vec![alt_base],
                            depth: row.depth,
                            ref_count: row.depth.saturating_sub(row.snp_alt_count),
                            alt_count: row.snp_alt_count,
                            qual: model.qual,
                            fs: 0.0,
                            pl: model.pl,
                            genotype_index: model.genotype_index,
                        });
                    }
                }
            }
        }
        if row.indel_alt_count >= 3 {
            if let Some(label) = row.indel_best_alt.strip_prefix("INS:") {
                let inserted = label.as_bytes().to_vec();
                let offset = (row.pos - region_start) as usize;
                if let Some(ref_byte) = local_ref_bases.get(offset).copied() {
                    let mut observations = Vec::with_capacity(row.depth as usize);
                    for _ in 0..row.depth.saturating_sub(row.indel_alt_count) {
                        observations.push(IndelObservation {
                            allele: IndelObservationAllele::Ref,
                            quality: 30,
                            is_reverse: false,
                        });
                    }
                    let alt_allele = IndelAllele::Insertion(inserted.clone());
                    for _ in 0..row.indel_alt_count {
                        observations.push(IndelObservation {
                            allele: IndelObservationAllele::Alt(alt_allele.clone()),
                            quality: 30,
                            is_reverse: false,
                        });
                    }
                    let model = indel_variant_model(&observations, &alt_allele);
                    if f64::from(model.qual) >= min_qual {
                        let anchor = normalize_base(ref_byte);
                        let mut alt_bases = vec![anchor];
                        alt_bases.extend_from_slice(&inserted);
                        events.push(VariantCall {
                            contig: contig.to_string(),
                            pos: row.pos,
                            id: None,
                            db: false,
                            ref_allele: vec![anchor],
                            alt_allele: alt_bases,
                            depth: row.depth,
                            ref_count: row.depth.saturating_sub(row.indel_alt_count),
                            alt_count: row.indel_alt_count,
                            qual: model.qual,
                            fs: 0.0,
                            pl: model.pl,
                            genotype_index: model.genotype_index,
                        });
                    }
                }
            }
        }
    }
    events.sort_by(|a, b| {
        a.pos
            .cmp(&b.pos)
            .then_with(|| a.ref_allele.cmp(&b.ref_allele))
            .then_with(|| a.alt_allele.cmp(&b.alt_allele))
    });
    events.dedup_by(|a, b| {
        a.pos == b.pos && a.ref_allele == b.ref_allele && a.alt_allele == b.alt_allele
    });
    events
}

fn collect_zero_candidate_simple_snp_seed_events(
    contig: &str,
    region_start: u64,
    local_ref_bases: &[u8],
    active_loci: &[ReplayActiveLocusRow],
) -> Vec<VariantCall> {
    collect_pileup_fallback_events(contig, region_start, local_ref_bases, active_loci, 0.0)
        .into_iter()
        .filter(|event| event.ref_allele.len() == 1 && event.alt_allele.len() == 1)
        .filter(|event| event.genotype_index == 0)
        .filter(|event| {
            active_loci
                .iter()
                .any(|row| active_locus_exact_simple_snp_support_without_indel(row, event))
        })
        .collect()
}

fn same_event_key(left: &VariantCall, right: &VariantCall) -> bool {
    left.pos == right.pos
        && left.ref_allele == right.ref_allele
        && left.alt_allele == right.alt_allele
}

fn merge_supplemental_haplotype(
    local_haplotypes: &mut Vec<LocalHaplotype>,
    mut haplotype: LocalHaplotype,
) {
    if haplotype.is_ref {
        return;
    }
    haplotype.event_indices.sort_unstable();
    haplotype.event_indices.dedup();
    if let Some(existing) = local_haplotypes
        .iter_mut()
        .find(|existing| existing.bases == haplotype.bases)
    {
        existing
            .event_indices
            .extend(haplotype.event_indices.iter().copied());
        existing.event_indices.sort_unstable();
        existing.event_indices.dedup();
        return;
    }
    local_haplotypes.push(haplotype);
}

fn haplotype_base_index_for_reference_pos(
    region_start: u64,
    haplotype: &LocalHaplotype,
    ref_pos: u64,
) -> Option<usize> {
    let mut current_ref = region_start;
    let mut current_hap = 0_usize;
    let mut op_len = 0_usize;

    for byte in haplotype.cigar.bytes() {
        if byte.is_ascii_digit() {
            op_len = op_len
                .saturating_mul(10)
                .saturating_add(usize::from(byte - b'0'));
            continue;
        }

        match byte as char {
            'M' | '=' | 'X' => {
                let op_end = current_ref + op_len as u64;
                if current_ref <= ref_pos && ref_pos < op_end {
                    return Some(current_hap + (ref_pos - current_ref) as usize);
                }
                current_ref = op_end;
                current_hap += op_len;
            }
            'I' | 'S' => {
                current_hap += op_len;
            }
            'D' | 'N' => {
                let op_end = current_ref + op_len as u64;
                if current_ref <= ref_pos && ref_pos < op_end {
                    return None;
                }
                current_ref = op_end;
            }
            'H' | 'P' => {}
            _ => {}
        }
        op_len = 0;
    }

    None
}

fn overlay_supplemental_snp_on_haplotype(
    region_start: u64,
    haplotype: &LocalHaplotype,
    existing_events: &[VariantCall],
    event: &VariantCall,
    event_idx: usize,
) -> Option<LocalHaplotype> {
    if haplotype.is_ref || event.ref_allele.len() != 1 || event.alt_allele.len() != 1 {
        return None;
    }
    if !haplotype.event_indices.iter().all(|existing_idx| {
        existing_events
            .get(*existing_idx)
            .is_some_and(|existing_event| {
                existing_event.ref_allele.len() == 1
                    && existing_event.alt_allele.len() == 1
                    && !events_overlap(existing_event, event)
            })
    }) {
        return None;
    }

    let base_idx = haplotype_base_index_for_reference_pos(region_start, haplotype, event.pos)?;
    if haplotype.bases.get(base_idx).copied()? != event.ref_allele[0] {
        return None;
    }

    let mut overlaid = haplotype.clone();
    overlaid.bases[base_idx] = event.alt_allele[0];
    overlaid.event_indices.push(event_idx);
    Some(overlaid)
}

fn supplement_missing_pileup_events(
    contig: &str,
    region_start: u64,
    local_ref_bases: &[u8],
    active_loci: &[ReplayActiveLocusRow],
    min_qual: f64,
    local_haplotypes: &mut Vec<LocalHaplotype>,
    valid_events: &mut Vec<VariantCall>,
) {
    let mut fallback_events = collect_pileup_fallback_events(
        contig,
        region_start,
        local_ref_bases,
        active_loci,
        min_qual,
    );
    if fallback_events.is_empty() && valid_events.is_empty() {
        // In some RNA-seq repeat contexts the pileup model under-scores a real
        // simple SNP, but PairHMM can recover it once the candidate exists.
        // Only relax seeding for zero-candidate regions, simple SNPs, and
        // loci without competing indel support.
        fallback_events = collect_zero_candidate_simple_snp_seed_events(
            contig,
            region_start,
            local_ref_bases,
            active_loci,
        );
    }
    if fallback_events.is_empty() {
        return;
    }

    if valid_events.is_empty() {
        *local_haplotypes = haplotypes_from_candidate_events(
            contig,
            region_start,
            local_ref_bases,
            &fallback_events,
        );
        *valid_events = fallback_events;
        return;
    }

    // Keep the original zero-candidate fallback behavior for both SNPs and
    // indels, but narrow the non-empty-region supplement path to simple SNPs.
    // Round8 full-call evidence showed that the Java-only gain here was SNP
    // dominated, while the new Rust-only regression included a large added
    // indel class.
    fallback_events.retain(|event| event.ref_allele.len() == 1 && event.alt_allele.len() == 1);
    if fallback_events.is_empty() {
        return;
    }

    let missing_events: Vec<VariantCall> = fallback_events
        .into_iter()
        .filter(|event| {
            !valid_events
                .iter()
                .any(|existing| same_event_key(existing, event))
                && !should_skip_weak_supplemental_snp_in_dense_snp_cluster(event, valid_events)
        })
        .collect();
    if missing_events.is_empty() {
        return;
    }

    let base_event_idx = valid_events.len();
    let existing_alt_haplotypes: Vec<LocalHaplotype> = local_haplotypes
        .iter()
        .filter(|hap| !hap.is_ref)
        .cloned()
        .collect();

    // Single-event synthetic haplotypes can underfit nearby multi-SNP reads.
    // When a pileup-strong missing SNP sits on the same reads as an already
    // assembled nearby ALT haplotype, also overlay that SNP onto the existing
    // ALT haplotype so PairHMM can score the combined sequence.
    for (event_offset, event) in missing_events.iter().enumerate() {
        let event_idx = base_event_idx + event_offset;
        for haplotype in &existing_alt_haplotypes {
            if let Some(overlaid) = overlay_supplemental_snp_on_haplotype(
                region_start,
                haplotype,
                valid_events,
                event,
                event_idx,
            ) {
                merge_supplemental_haplotype(local_haplotypes, overlaid);
            }
        }
    }

    let supplemental_haplotypes =
        haplotypes_from_candidate_events(contig, region_start, local_ref_bases, &missing_events);
    valid_events.extend(missing_events.iter().cloned());

    for mut haplotype in supplemental_haplotypes
        .into_iter()
        .filter(|hap| !hap.is_ref)
    {
        for event_idx in &mut haplotype.event_indices {
            *event_idx += base_event_idx;
        }
        merge_supplemental_haplotype(local_haplotypes, haplotype);
    }
}

fn should_skip_weak_supplemental_snp_in_dense_snp_cluster(
    event: &VariantCall,
    valid_events: &[VariantCall],
) -> bool {
    if event.ref_allele.len() != 1
        || event.alt_allele.len() != 1
        || event.alt_count > WEAK_SUPPLEMENTAL_CLUSTER_SNP_MAX_ALT_COUNT
    {
        return false;
    }

    let mut nearby_positions = Vec::with_capacity(2);
    for existing in valid_events {
        if existing.ref_allele.len() != 1 || existing.alt_allele.len() != 1 {
            continue;
        }
        if existing.pos.abs_diff(event.pos) > SNP_CLUSTER_WINDOW {
            continue;
        }
        if !nearby_positions.contains(&existing.pos) {
            nearby_positions.push(existing.pos);
            if nearby_positions.len() >= 2 {
                return true;
            }
        }
    }

    false
}

fn rescue_collapsed_strong_snp_cluster_from_pileup(
    contig: &str,
    region_start: u64,
    local_ref_bases: &[u8],
    active_loci: &[ReplayActiveLocusRow],
    valid_events: &[VariantCall],
    min_qual: f64,
) -> Vec<VariantCall> {
    let fallback_events = collect_pileup_fallback_events(
        contig,
        region_start,
        local_ref_bases,
        active_loci,
        min_qual,
    );
    let rescued = exact_strong_simple_snp_pileup_matches_from_fallback_events(
        &fallback_events,
        active_loci,
        valid_events,
    );

    if rescued.len() >= 2 {
        return rescued;
    }

    if rescued.len() == 1 && fallback_events.len() == 1 {
        return rescued;
    }

    if rescued.len() == 1 && has_high_confidence_single_snp_rescue_support(active_loci, &rescued[0])
    {
        return rescued;
    }

    Vec::new()
}

fn exact_strong_simple_snp_pileup_matches_from_fallback_events(
    fallback_events: &[VariantCall],
    active_loci: &[ReplayActiveLocusRow],
    valid_events: &[VariantCall],
) -> Vec<VariantCall> {
    fallback_events
        .iter()
        .filter(|event| event.ref_allele.len() == 1 && event.alt_allele.len() == 1)
        .filter(|event| {
            valid_events
                .iter()
                .any(|valid| same_event_key(valid, event))
                && active_loci
                    .iter()
                    .any(|row| active_locus_exact_simple_snp_support_without_indel(row, event))
        })
        .cloned()
        .collect()
}

fn exact_strong_simple_snp_pileup_matches(
    contig: &str,
    region_start: u64,
    local_ref_bases: &[u8],
    active_loci: &[ReplayActiveLocusRow],
    valid_events: &[VariantCall],
    min_qual: f64,
) -> Vec<VariantCall> {
    let fallback_events = collect_pileup_fallback_events(
        contig,
        region_start,
        local_ref_bases,
        active_loci,
        min_qual,
    );
    exact_strong_simple_snp_pileup_matches_from_fallback_events(
        &fallback_events,
        active_loci,
        valid_events,
    )
}

fn active_locus_exact_simple_snp_support(row: &ReplayActiveLocusRow, event: &VariantCall) -> bool {
    event.ref_allele.len() == 1
        && event.alt_allele.len() == 1
        && event.pos == row.pos
        && event.ref_allele[0] == row.ref_base
        && row
            .snp_best_alt
            .as_bytes()
            .first()
            .is_some_and(|alt| event.alt_allele[0] == *alt)
}

fn active_locus_exact_simple_snp_support_without_indel(
    row: &ReplayActiveLocusRow,
    event: &VariantCall,
) -> bool {
    row.indel_alt_count == 0 && active_locus_exact_simple_snp_support(row, event)
}

fn active_locus_high_confidence_single_snp_rescue_support(
    row: &ReplayActiveLocusRow,
    event: &VariantCall,
) -> bool {
    active_locus_exact_simple_snp_support_without_indel(row, event)
        && row.snp_alt_count >= STRONG_SINGLE_SNP_RESCUE_MIN_ALT_COUNT
        && row.alt_fraction >= STRONG_SINGLE_SNP_RESCUE_MIN_ALT_FRACTION
}

fn has_high_confidence_single_snp_rescue_support(
    active_loci: &[ReplayActiveLocusRow],
    event: &VariantCall,
) -> bool {
    active_loci
        .iter()
        .any(|row| active_locus_high_confidence_single_snp_rescue_support(row, event))
}

fn prune_unsupported_simple_snp_calls_in_dense_clusters(
    final_calls: &mut Vec<VariantCall>,
    active_loci: &[ReplayActiveLocusRow],
) {
    if final_calls.len() < 3 {
        return;
    }

    let dense_unsupported_snp_keys: Vec<(u64, Vec<u8>, Vec<u8>)> = final_calls
        .iter()
        .filter(|call| call.ref_allele.len() == 1 && call.alt_allele.len() == 1)
        .filter(|call| {
            !active_loci
                .iter()
                .any(|row| active_locus_exact_simple_snp_support(row, call))
        })
        .filter(|call| {
            let mut nearby_positions = Vec::with_capacity(2);
            for other in final_calls.iter() {
                if same_event_key(call, other)
                    || other.ref_allele.len() != 1
                    || other.alt_allele.len() != 1
                    || other.pos.abs_diff(call.pos) > SNP_CLUSTER_WINDOW
                {
                    continue;
                }
                if !nearby_positions.contains(&other.pos) {
                    nearby_positions.push(other.pos);
                    if nearby_positions.len() >= 2 {
                        return true;
                    }
                }
            }
            false
        })
        .map(|call| (call.pos, call.ref_allele.clone(), call.alt_allele.clone()))
        .collect();

    if dense_unsupported_snp_keys.is_empty() {
        return;
    }

    final_calls.retain(|call| {
        !dense_unsupported_snp_keys
            .iter()
            .any(|(pos, ref_allele, alt_allele)| {
                call.pos == *pos && call.ref_allele == *ref_allele && call.alt_allele == *alt_allele
            })
    });
}

fn merge_missing_strong_snp_cluster_rescues_from_pileup(
    final_calls: &mut Vec<VariantCall>,
    contig: &str,
    region_start: u64,
    local_ref_bases: &[u8],
    active_loci: &[ReplayActiveLocusRow],
    valid_events: &[VariantCall],
    min_qual: f64,
) {
    let rescued = exact_strong_simple_snp_pileup_matches(
        contig,
        region_start,
        local_ref_bases,
        active_loci,
        valid_events,
        min_qual,
    );
    if rescued.len() < 2
        && !rescued
            .first()
            .is_some_and(|event| has_high_confidence_single_snp_rescue_support(active_loci, event))
    {
        return;
    }

    for rescued_call in rescued {
        if !final_calls
            .iter()
            .any(|final_call| same_event_key(final_call, &rescued_call))
        {
            final_calls.push(rescued_call);
        }
    }
}

