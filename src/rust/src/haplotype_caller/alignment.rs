fn is_regular_allele(allele: &[u8]) -> bool {
    allele
        .iter()
        .all(|&b| matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't'))
}

fn is_regular_bases(bases: &[u8]) -> bool {
    bases.iter().all(|&b| is_acgt(normalize_base(b)))
}

fn haplotype_to_reference_sw_parameters() -> crate::smith_waterman::SWParameters {
    crate::smith_waterman::SWParameters {
        match_value: 200,
        mismatch_penalty: -150,
        gap_open_penalty: -260,
        gap_extend_penalty: -11,
    }
}

fn align_haplotype_to_reference(
    ref_bases: &[u8],
    haplotype_bases: &[u8],
) -> Option<crate::smith_waterman::SWAlignmentResult> {
    let align_result = crate::smith_waterman::align(
        ref_bases,
        haplotype_bases,
        &haplotype_to_reference_sw_parameters(),
        crate::smith_waterman::SWOverhangStrategy::SoftClip,
    );

    if align_result.alignment_offset > 0
        || align_result
            .cigar
            .iter()
            .any(|ce| matches!(ce, rust_htslib::bam::record::Cigar::SoftClip(_)))
    {
        return None;
    }

    Some(align_result)
}

pub fn extract_variants_from_cigar(
    contig: &str,
    ref_bases: &[u8],
    alt_bases: &[u8],
    cigar: &rust_htslib::bam::record::CigarString,
    alignment_offset: i32,
    region_start: u64,
    max_mnp_distance: usize,
) -> Vec<VariantCall> {
    use rust_htslib::bam::record::Cigar::*;

    let mut ref_pos = alignment_offset;
    if ref_pos < 0 {
        return Vec::new();
    }

    let mut alignment_pos = 0;
    let mut proposed_events = Vec::new();
    let num_cigar_elements = cigar.len();

    for (cigar_index, ce) in cigar.iter().enumerate() {
        match ce {
            Ins(len) => {
                let element_length = *len as usize;
                if element_length <= 10
                    && ref_pos > 0
                    && cigar_index > 0
                    && cigar_index < num_cigar_elements - 1
                {
                    let insertion_start = region_start + ref_pos as u64 - 1;
                    let ref_byte = ref_bases[ref_pos as usize - 1];
                    let mut insertion_bases = vec![ref_byte];
                    insertion_bases.extend_from_slice(
                        &alt_bases[alignment_pos..alignment_pos + element_length],
                    );

                    if is_regular_allele(&[ref_byte]) && is_regular_allele(&insertion_bases) {
                        proposed_events.push(VariantCall {
                            contig: contig.to_string(),
                            pos: insertion_start,
                            id: None,
                            db: false,
                            ref_allele: vec![ref_byte],
                            alt_allele: insertion_bases,
                            depth: 0,
                            ref_count: 0,
                            alt_count: 0,
                            qual: 0,
                            fs: 0.0,
                            pl: [0, 0, 0],
                            genotype_index: 0,
                        });
                    }
                }
                alignment_pos += element_length;
            }
            SoftClip(len) => {
                alignment_pos += *len as usize;
            }
            Del(len) => {
                let element_length = *len as usize;
                if element_length <= 10 && ref_pos > 0 {
                    let deletion_start = region_start + ref_pos as u64 - 1;
                    let ref_byte = ref_bases[ref_pos as usize - 1];
                    let mut deletion_bases = vec![ref_byte];
                    deletion_bases.extend_from_slice(
                        &ref_bases[ref_pos as usize..ref_pos as usize + element_length],
                    );

                    if is_regular_allele(&deletion_bases) && is_regular_allele(&[ref_byte]) {
                        proposed_events.push(VariantCall {
                            contig: contig.to_string(),
                            pos: deletion_start,
                            id: None,
                            db: false,
                            ref_allele: deletion_bases,
                            alt_allele: vec![ref_byte],
                            depth: 0,
                            ref_count: 0,
                            alt_count: 0,
                            qual: 0,
                            fs: 0.0,
                            pl: [0, 0, 0],
                            genotype_index: 0,
                        });
                    }
                }
                ref_pos += element_length as i32;
            }
            Match(len) | Equal(len) | Diff(len) => {
                let element_length = *len as usize;
                let mut mismatch_offsets = std::collections::VecDeque::new();

                for offset in 0..element_length {
                    let r_idx = ref_pos as usize + offset;
                    let a_idx = alignment_pos + offset;
                    if r_idx < ref_bases.len() && a_idx < alt_bases.len() {
                        let ref_byte = ref_bases[r_idx];
                        let alt_byte = alt_bases[a_idx];
                        // we ignore N vs N mismatches in practice, but keeping simple
                        if ref_byte != alt_byte {
                            mismatch_offsets.push_back(offset);
                        }
                    }
                }

                while let Some(start) = mismatch_offsets.pop_front() {
                    let mut end = start;
                    while let Some(&next) = mismatch_offsets.front() {
                        if next - end <= max_mnp_distance {
                            end = mismatch_offsets.pop_front().unwrap();
                        } else {
                            break;
                        }
                    }

                    let ref_allele =
                        ref_bases[ref_pos as usize + start..=ref_pos as usize + end].to_vec();
                    let alt_allele =
                        alt_bases[alignment_pos + start..=alignment_pos + end].to_vec();

                    if is_regular_allele(&ref_allele) && is_regular_allele(&alt_allele) {
                        proposed_events.push(VariantCall {
                            contig: contig.to_string(),
                            pos: region_start + ref_pos as u64 + start as u64,
                            id: None,
                            db: false,
                            ref_allele,
                            alt_allele,
                            depth: 0,
                            ref_count: 0,
                            alt_count: 0,
                            qual: 0,
                            fs: 0.0,
                            pl: [0, 0, 0],
                            genotype_index: 0,
                        });
                    }
                }

                ref_pos += element_length as i32;
                alignment_pos += element_length;
            }
            _ => {
                // skip others for now
            }
        }
    }

    proposed_events
}

pub fn assemble_haplotypes(
    contig: &str,
    region_start: u64,
    ref_bases: &[u8],
    reads_bases: &[Vec<u8>],
    kmer_sizes: &[usize],
    max_mnp_distance: usize,
) -> (Vec<LocalHaplotype>, Vec<VariantCall>) {
    use std::collections::HashSet;
    let mut assembled_haplotypes_set = HashSet::new();
    let mut local_haps = Vec::new();

    let ref_hap = LocalHaplotype {
        bases: ref_bases.to_vec(),
        is_ref: true,
        cigar: format!("{}M", ref_bases.len()),
        event_indices: Vec::new(),
    };
    local_haps.push(ref_hap);
    assembled_haplotypes_set.insert(ref_bases.to_vec());

    // We'll store events per haplotype as Vec<Vec<VariantCall>>
    let mut hap_events = vec![Vec::new()]; // first is ref, no events

    for &kmer_size in kmer_sizes {
        if ref_bases.len() < kmer_size {
            continue;
        }
        let mut graph = crate::assembly::ReadThreadingGraph::new(kmer_size);
        graph.add_sequence(ref_bases, true);

        for read_bases in reads_bases {
            if read_bases.len() >= kmer_size {
                graph.add_sequence(read_bases, false);
            }
        }

        let source_kmer = &ref_bases[0..kmer_size];
        let sink_kmer = &ref_bases[ref_bases.len() - kmer_size..];

        let source_idx = graph.get_or_create_vertex(source_kmer);
        let sink_idx = graph.get_or_create_vertex(sink_kmer);

        graph.prune(2); // min_prune_factor = 2

        let best_paths = graph.find_best_haplotypes(source_idx, sink_idx, 128);

        let mut found_nonref = false;
        for path in best_paths {
            let seq = graph.reconstruct_sequence(&path);
            if !is_regular_bases(&seq) {
                continue;
            }
            if !assembled_haplotypes_set.contains(&seq) {
                assembled_haplotypes_set.insert(seq.clone());

                let Some(align_result) = align_haplotype_to_reference(ref_bases, &seq) else {
                    continue;
                };

                let events = extract_variants_from_cigar(
                    contig,
                    ref_bases,
                    &seq,
                    &align_result.cigar,
                    align_result.alignment_offset,
                    region_start,
                    max_mnp_distance,
                );

                hap_events.push(events);

                local_haps.push(LocalHaplotype {
                    bases: seq,
                    is_ref: false,
                    cigar: align_result.cigar.to_string(),
                    event_indices: Vec::new(), // To be filled later
                });
                found_nonref = true;
            }
        }
        // Like GATK: if this kmer size produced non-ref haplotypes, stop
        if found_nonref {
            break;
        }
    }

    // Deduplicate events across all haplotypes
    let mut unique_events: Vec<VariantCall> = Vec::new();
    for events in &hap_events {
        for event in events {
            if !unique_events.iter().any(|e| {
                e.pos == event.pos
                    && e.ref_allele == event.ref_allele
                    && e.alt_allele == event.alt_allele
            }) {
                unique_events.push(event.clone());
            }
        }
    }

    // Now map event indices back to each haplotype
    for (hap_idx, events) in hap_events.iter().enumerate() {
        for event in events {
            if let Some(idx) = unique_events.iter().position(|e| {
                e.pos == event.pos
                    && e.ref_allele == event.ref_allele
                    && e.alt_allele == event.alt_allele
            }) {
                if !local_haps[hap_idx].event_indices.contains(&idx) {
                    local_haps[hap_idx].event_indices.push(idx);
                }
            }
        }
    }

    (local_haps, unique_events)
}
