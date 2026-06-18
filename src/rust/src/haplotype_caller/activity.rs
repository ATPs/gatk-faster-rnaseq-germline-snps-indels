fn validate_haplotype_caller_config(config: &HaplotypeCallerConfig) -> Result<()> {
    if !config.input_bam.exists() {
        bail!("input BAM not found: {}", config.input_bam.display());
    }
    if !config.reference.exists() {
        bail!("reference FASTA not found: {}", config.reference.display());
    }
    if !config.input_interval_list.exists() {
        bail!(
            "input interval_list not found: {}",
            config.input_interval_list.display()
        );
    }
    if let Some(dbsnp) = &config.dbsnp {
        if !dbsnp.exists() {
            bail!("dbsnp VCF not found: {}", dbsnp.display());
        }
    }
    if config.standard_min_confidence_threshold_for_calling < 0.0 {
        bail!("--standard-min-confidence-threshold-for-calling must be non-negative");
    }
    if config.threads == 0 {
        bail!("--threads must be at least 1");
    }
    if config.memory_gb == 0 {
        bail!("--memory-gb must be at least 1");
    }
    if config.native_pair_hmm_threads == 0 {
        bail!("--native-pair-hmm-threads must be at least 1");
    }
    if config.pair_hmm_implementation != "rust" && config.pair_hmm_implementation != "native" {
        bail!("--pair-hmm-implementation must be rust or native");
    }
    Ok(())
}

fn scan_call_partition(
    config: &HaplotypeCallerConfig,
    windows: &[FetchWindow],
) -> Result<CallWorkerOutput> {
    let mut output = CallWorkerOutput::default();
    let mut bam = bam::IndexedReader::from_path(&config.input_bam)
        .with_context(|| format!("opening indexed BAM {}", config.input_bam.display()))?;
    let bam_tid_by_name = bam_tid_by_name(bam.header())?;
    let reference = faidx::Reader::from_path(&config.reference)
        .with_context(|| format!("opening reference FASTA {}", config.reference.display()))?;

    for window in windows {
        scan_call_window(
            config,
            &bam_tid_by_name,
            &reference,
            &mut bam,
            window,
            &mut output,
        )?;
    }
    Ok(output)
}

fn scan_replay_partition(
    config: &HaplotypeCallerConfig,
    windows: &[FetchWindow],
) -> Result<ReplayWorkerOutput> {
    let mut output = ReplayWorkerOutput::default();
    let mut bam = bam::IndexedReader::from_path(&config.input_bam)
        .with_context(|| format!("opening indexed BAM {}", config.input_bam.display()))?;
    let bam_tid_by_name = bam_tid_by_name(bam.header())?;
    let reference = faidx::Reader::from_path(&config.reference)
        .with_context(|| format!("opening reference FASTA {}", config.reference.display()))?;

    for window in windows {
        scan_replay_window(
            config,
            &bam_tid_by_name,
            &reference,
            &mut bam,
            window,
            &mut output,
        )?;
    }
    Ok(output)
}

fn select_assembly_reads<'a>(
    read_assembly_segments_list: &'a [Vec<Vec<u8>>],
    read_use_for_assembly_list: &[bool],
) -> Vec<&'a Vec<u8>> {
    read_assembly_segments_list
        .iter()
        .zip(read_use_for_assembly_list.iter().copied())
        .flat_map(|(segments, use_for_assembly)| {
            use_for_assembly
                .then_some(segments.as_slice())
                .unwrap_or(&[])
                .iter()
        })
        .collect()
}

fn discover_active_regions(
    config: &HaplotypeCallerConfig,
    tid: u32,
    window: &FetchWindow,
    ref_bases: &[u8],
    bam: &mut bam::IndexedReader,
) -> Result<Vec<ActiveRegionSpan>> {
    let mut active_regions: Vec<Interval> = Vec::new();
    let mut activity_scores: std::collections::HashMap<u64, u32> = std::collections::HashMap::new();

    bam.fetch((tid as i32, (window.start - 1) as i64, window.end as i64))
        .with_context(|| {
            format!(
                "fetching BAM region {}:{}-{} for active region discovery",
                window.contig, window.start, window.end
            )
        })?;
    let mut pileups = bam.pileup();
    pileups.set_max_depth(PIPELINE_MAX_DEPTH);
    let mut interval_cursor = 0_usize;
    let mut current_region: Option<Interval> = None;

    for pileup in pileups {
        let pileup = pileup.with_context(|| {
            format!(
                "reading pileup for {}:{}-{} active discovery",
                window.contig, window.start, window.end
            )
        })?;
        if pileup.tid() != tid {
            continue;
        }
        let pos0 = u64::from(pileup.pos());
        if pos0 < window.start - 1 || pos0 >= window.end {
            continue;
        }
        let pos1 = pos0 + 1;
        if !position_is_requested(&window.intervals, pos1, &mut interval_cursor) {
            continue;
        }

        let ref_base = normalize_base(ref_bases[(pos0 - (window.start - 1)) as usize]);
        let snp_evidence = pileup_snp_evidence(
            &pileup,
            PIPELINE_MIN_BASEQ,
            PIPELINE_MIN_MAPQ,
            PIPELINE_MIN_TAIL_QUALITY,
            config.exclude_supplementary,
            config.dont_use_soft_clipped_bases,
        );
        let indel_evidence = pileup_indel_evidence(
            &pileup,
            PIPELINE_MIN_BASEQ,
            PIPELINE_MIN_MAPQ,
            PIPELINE_MIN_TAIL_QUALITY,
            config.exclude_supplementary,
            config.dont_use_soft_clipped_bases,
        );
        let depth = snp_evidence.counts.depth.max(indel_evidence.counts.depth);

        let (snp_active, snp_qual) = is_active_locus(base_index(ref_base), &snp_evidence, depth);
        let (indel_active, indel_qual) = is_active_indel(&indel_evidence);
        let is_active = snp_active || indel_active;
        let qual = snp_qual.max(indel_qual);

        if is_active {
            activity_scores.insert(pos1, qual);
            if let Some(ref mut reg) = current_region {
                if pos1 <= reg.end + ACTIVE_REGION_MAX_GAP {
                    reg.end = pos1;
                } else {
                    active_regions.push(reg.clone());
                    *reg = Interval {
                        contig: window.contig.clone(),
                        start: pos1,
                        end: pos1,
                    };
                }
            } else {
                current_region = Some(Interval {
                    contig: window.contig.clone(),
                    start: pos1,
                    end: pos1,
                });
            }
        }
    }
    if let Some(reg) = current_region {
        active_regions.push(reg);
    }

    let mut split_active_regions = Vec::new();
    for reg in active_regions {
        let mut start = reg.start;
        let end = reg.end;
        while end - start + 1 > 300 {
            let search_start = start + 50;
            let search_end = (start + 300).min(end - 50);

            let mut best_cut = search_start;
            let mut min_score = u32::MAX;

            for p in search_start..=search_end {
                let score = *activity_scores.get(&p).unwrap_or(&0);
                if score < min_score {
                    min_score = score;
                    best_cut = p;
                } else if score == min_score {
                    let mid = (search_start + search_end) / 2;
                    if (p as i64 - mid as i64).abs() < (best_cut as i64 - mid as i64).abs() {
                        best_cut = p;
                    }
                }
            }

            split_active_regions.push(Interval {
                contig: reg.contig.clone(),
                start,
                end: best_cut - 1,
            });
            start = best_cut;
        }
        split_active_regions.push(Interval {
            contig: reg.contig.clone(),
            start,
            end,
        });
    }

    let mut coalesced = Vec::new();
    for reg in split_active_regions {
        let padded = pad_active_region(&reg, window);
        coalesced.push(ActiveRegionSpan {
            active: reg,
            padded,
        });
    }

    Ok(coalesced)
}

fn pad_active_region(interval: &Interval, window: &FetchWindow) -> Interval {
    Interval {
        contig: interval.contig.clone(),
        start: interval
            .start
            .saturating_sub(ACTIVE_REGION_PADDING)
            .max(window.start),
        end: interval
            .end
            .saturating_add(ACTIVE_REGION_PADDING)
            .min(window.end),
    }
}

fn collect_call_active_loci_rows(
    config: &HaplotypeCallerConfig,
    tid: u32,
    window: &FetchWindow,
    active_span: &Interval,
    region: &Interval,
    ref_bases: &[u8],
    bam: &mut bam::IndexedReader,
) -> Result<Vec<ReplayActiveLocusRow>> {
    let mut active_loci_rows = Vec::new();
    let region_label = region_name(&region.contig, region.start, region.end);
    bam.fetch((
        tid as i32,
        (active_span.start - 1) as i64,
        active_span.end as i64,
    ))
    .with_context(|| {
        format!(
            "fetching BAM region for pileup fallback {}:{}-{}",
            active_span.contig, active_span.start, active_span.end
        )
    })?;
    for pileup_result in bam.pileup() {
        let pileup = pileup_result?;
        let pos0 = pileup.pos() as u64;
        let pos1 = pos0 + 1;
        if pos1 < active_span.start || pos1 > active_span.end {
            continue;
        }
        let ref_base = normalize_base(ref_bases[(pos0 - (window.start - 1)) as usize]);
        let snp_evidence = pileup_snp_evidence(
            &pileup,
            PIPELINE_MIN_BASEQ,
            PIPELINE_MIN_MAPQ,
            PIPELINE_MIN_TAIL_QUALITY,
            config.exclude_supplementary,
            config.dont_use_soft_clipped_bases,
        );
        let indel_evidence = pileup_indel_evidence(
            &pileup,
            PIPELINE_MIN_BASEQ,
            PIPELINE_MIN_MAPQ,
            PIPELINE_MIN_TAIL_QUALITY,
            config.exclude_supplementary,
            config.dont_use_soft_clipped_bases,
        );
        let ref_index = base_index(ref_base);
        let snp_alt = best_snp_alt(ref_index, &snp_evidence);
        let indel_alt = best_indel_alt(&indel_evidence);
        let snp_alt_count = snp_alt.map(|(_, count)| count).unwrap_or(0);
        let indel_alt_count = indel_alt.map(|(_, count)| *count).unwrap_or(0);
        let depth = snp_evidence.counts.depth.max(indel_evidence.counts.depth);
        let is_active = is_active_locus(ref_index, &snp_evidence, depth).0
            || is_active_indel(&indel_evidence).0;
        if is_active {
            let best_alt_count = snp_alt_count.max(indel_alt_count);
            let alt_fraction = if depth == 0 {
                0.0
            } else {
                f64::from(best_alt_count) / f64::from(depth)
            };
            active_loci_rows.push(ReplayActiveLocusRow {
                contig: region.contig.clone(),
                pos: pos1,
                region: region_label.clone(),
                ref_base,
                depth,
                snp_alt_count,
                snp_best_alt: snp_alt
                    .map(|(base_index, _)| base_from_index(base_index) as char)
                    .map(|base| base.to_string())
                    .unwrap_or_default(),
                indel_alt_count,
                indel_best_alt: indel_alt
                    .as_ref()
                    .map(|(allele, _)| indel_allele_label(allele))
                    .unwrap_or_default(),
                alt_fraction,
                active_probability_proxy: 1.0,
            });
        }
    }

    Ok(active_loci_rows)
}

fn scan_call_window(
    config: &HaplotypeCallerConfig,
    bam_tid_by_name: &HashMap<String, u32>,
    reference: &faidx::Reader,
    bam: &mut bam::IndexedReader,
    window: &FetchWindow,
    output: &mut CallWorkerOutput,
) -> Result<()> {
    let tid = *bam_tid_by_name.get(&window.contig).with_context(|| {
        format!(
            "contig '{}' from {} is not present in BAM header",
            window.contig,
            config.input_interval_list.display()
        )
    })?;
    let ref_len = reference.fetch_seq_len(&window.contig);
    if ref_len == 0 {
        bail!(
            "contig '{}' is not present in reference FASTA {}",
            window.contig,
            config.reference.display()
        );
    }
    if window.end > ref_len {
        bail!(
            "interval {}:{}-{} extends past FASTA contig length {}",
            window.contig,
            window.start,
            window.end,
            ref_len
        );
    }

    let ref_end = window
        .end
        .saturating_add(u64::from(MAX_BOOTSTRAP_INDEL_LEN))
        .min(ref_len);
    let ref_bases = reference
        .fetch_seq(
            &window.contig,
            (window.start - 1) as usize,
            (ref_end - 1) as usize,
        )
        .with_context(|| {
            format!(
                "fetching reference sequence {}:{}-{}",
                window.contig, window.start, ref_end
            )
        })?;

    let active_regions = discover_active_regions(config, tid, window, &ref_bases, bam)?;

    for region in active_regions {
        let active_span = region.active.clone();
        let region = region.padded;
        let active_loci_rows = collect_call_active_loci_rows(
            config,
            tid,
            window,
            &active_span,
            &region,
            &ref_bases,
            bam,
        )?;
        bam.fetch((tid as i32, (region.start - 1) as i64, region.end as i64))
            .with_context(|| {
                format!(
                    "fetching BAM region {}:{}-{}",
                    region.contig, region.start, region.end
                )
            })?;

        let mut read_bases_list = Vec::new();
        let mut read_quals_list = Vec::new();
        let mut read_ins_quals_list = Vec::new();
        let mut read_del_quals_list = Vec::new();
        let mut read_is_reverse_list = Vec::new();
        let mut read_ref_spans = Vec::new();
        let mut read_use_for_assembly_list = Vec::new();
        let mut read_assembly_segments_list = Vec::new();
        let mut reads_by_start = std::collections::HashMap::new();

        for r in bam.records() {
            let record = r?;
            if !read_passes_hc_filter(&record, PIPELINE_MIN_MAPQ, config.exclude_supplementary) {
                continue;
            }
            let Some(prepared_read) = prepare_hmm_read(
                &record,
                PIPELINE_MIN_TAIL_QUALITY,
                config.dont_use_soft_clipped_bases,
                region.start,
                region.end,
            ) else {
                continue;
            };

            let start_pos = record.pos();
            let count = reads_by_start.entry(start_pos).or_insert(0);
            let use_for_assembly = *count < 50;
            if use_for_assembly {
                *count += 1;
            }

            read_assembly_segments_list.push(prepared_read.assembly_segments);
            read_bases_list.push(prepared_read.bases);
            read_quals_list.push(prepared_read.quals);
            read_ins_quals_list.push(prepared_read.ins_quals);
            read_del_quals_list.push(prepared_read.del_quals);
            read_is_reverse_list.push(record.is_reverse());
            read_ref_spans.push(prepared_read.ref_span);
            read_use_for_assembly_list.push(use_for_assembly);
        }

        if read_bases_list.is_empty() {
            continue;
        }

        let ref_region_start_offset = (region.start - window.start) as usize;
        let ref_region_end_offset = (region.end - window.start) as usize;
        let local_ref_bases = &ref_bases[ref_region_start_offset..=ref_region_end_offset];

        // 1. Assemble haplotypes instead of pileup!
        let max_mnp_distance = 0; // Default GATK
        let assembly_reads =
            select_assembly_reads(&read_assembly_segments_list, &read_use_for_assembly_list);
        let assembly_reads_owned: Vec<Vec<u8>> =
            assembly_reads.iter().map(|r| (*r).clone()).collect();
        let (mut local_haplotypes, mut valid_events) = assemble_haplotypes(
            &region.contig,
            region.start,
            local_ref_bases,
            &assembly_reads_owned,
            &[10, 25],
            max_mnp_distance,
        );
        supplement_missing_pileup_events(
            &region.contig,
            region.start,
            local_ref_bases,
            &active_loci_rows,
            config.standard_min_confidence_threshold_for_calling,
            &mut local_haplotypes,
            &mut valid_events,
        );
        filter_non_acgt_haplotypes_for_single_snp_region(&mut local_haplotypes, &valid_events);

        if valid_events.is_empty() {
            continue;
        }

        let n_reads = read_bases_list.len();
        let mut read_haplotype_likelihoods: Vec<Vec<f64>> = Vec::with_capacity(n_reads);

        for i in 0..n_reads {
            let r_bases = &read_bases_list[i];
            let r_quals = &read_quals_list[i];
            let read_ins_quals = &read_ins_quals_list[i];
            let read_del_quals = &read_del_quals_list[i];
            let gcp = 10;

            let mut hap_likelihoods = Vec::with_capacity(local_haplotypes.len());
            for hap in &local_haplotypes {
                let score = pair_hmm::compute_read_likelihood_given_haplotype(
                    &hap.bases,
                    r_bases,
                    r_quals,
                    &read_ins_quals,
                    &read_del_quals,
                    gcp,
                );
                hap_likelihoods.push(score);
            }
            read_haplotype_likelihoods.push(hap_likelihoods);
        }

        if read_haplotype_likelihoods.is_empty() {
            continue;
        }

        let mut final_calls = genotype_assembled_events(
            &local_haplotypes,
            &valid_events,
            &read_haplotype_likelihoods,
            &read_is_reverse_list,
            &read_ref_spans,
            config.standard_min_confidence_threshold_for_calling,
        );
        prune_unsupported_simple_snp_calls_in_dense_clusters(&mut final_calls, &active_loci_rows);
        merge_missing_strong_snp_cluster_rescues_from_pileup(
            &mut final_calls,
            &region.contig,
            region.start,
            local_ref_bases,
            &active_loci_rows,
            &valid_events,
            config.standard_min_confidence_threshold_for_calling,
        );
        if final_calls.is_empty() {
            final_calls = rescue_collapsed_strong_snp_cluster_from_pileup(
                &region.contig,
                region.start,
                local_ref_bases,
                &active_loci_rows,
                &valid_events,
                config.standard_min_confidence_threshold_for_calling,
            );
        }
        output.variants.extend(final_calls);
    }
    Ok(())
}

fn scan_replay_window(
    config: &HaplotypeCallerConfig,
    bam_tid_by_name: &HashMap<String, u32>,
    reference: &faidx::Reader,
    bam: &mut bam::IndexedReader,
    window: &FetchWindow,
    output: &mut ReplayWorkerOutput,
) -> Result<()> {
    let tid = *bam_tid_by_name.get(&window.contig).with_context(|| {
        format!(
            "contig '{}' from {} is not present in BAM header",
            window.contig,
            config.input_interval_list.display()
        )
    })?;
    let ref_len = reference.fetch_seq_len(&window.contig);
    if ref_len == 0 {
        bail!(
            "contig '{}' is not present in reference FASTA {}",
            window.contig,
            config.reference.display()
        );
    }
    if window.end > ref_len {
        bail!(
            "interval {}:{}-{} extends past FASTA contig length {}",
            window.contig,
            window.start,
            window.end,
            ref_len
        );
    }

    let ref_end = window
        .end
        .saturating_add(u64::from(MAX_BOOTSTRAP_INDEL_LEN))
        .min(ref_len);
    let ref_bases = reference
        .fetch_seq(
            &window.contig,
            (window.start - 1) as usize,
            (ref_end - 1) as usize,
        )
        .with_context(|| {
            format!(
                "fetching reference sequence {}:{}-{}",
                window.contig, window.start, ref_end
            )
        })?;

    let active_regions = discover_active_regions(config, tid, window, &ref_bases, bam)?;

    let mut interval_cursor = 0_usize;
    for discovered_region in active_regions {
        let active_span = &discovered_region.active;
        let region_interval = &discovered_region.padded;
        bam.fetch((
            tid as i32,
            (region_interval.start - 1) as i64,
            region_interval.end as i64,
        ))
        .with_context(|| {
            format!(
                "fetching BAM region {}:{}-{}",
                region_interval.contig, region_interval.start, region_interval.end
            )
        })?;

        let region = region_name(
            &region_interval.contig,
            region_interval.start,
            region_interval.end,
        );
        let mut active_region = ReplayActiveRegionRow {
            contig: region_interval.contig.clone(),
            start: region_interval.start,
            end: region_interval.end,
            region: region.clone(),
            active_start: active_span.start,
            active_end: active_span.end,
            padded_start: region_interval.start,
            padded_end: region_interval.end,
            ..ReplayActiveRegionRow::default()
        };

        let mut pileups = bam.pileup();
        pileups.set_max_depth(PIPELINE_MAX_DEPTH);
        for pileup in pileups {
            let pileup = pileup.with_context(|| {
                format!(
                    "reading pileup for {}:{}-{}",
                    region_interval.contig, region_interval.start, region_interval.end
                )
            })?;
            if pileup.tid() != tid {
                continue;
            }
            let pos0 = u64::from(pileup.pos());
            if pos0 < region_interval.start - 1 || pos0 >= region_interval.end {
                continue;
            }
            let pos1 = pos0 + 1;
            if !position_is_requested(&window.intervals, pos1, &mut interval_cursor) {
                continue;
            }

            let ref_base = normalize_base(ref_bases[(pos0 - (window.start - 1)) as usize]);
            let row_context = ReplayRowContext {
                region: &region,
                pos: pos1,
            };
            let (snp_evidence, snp_rows) = pileup_snp_evidence_with_rows(
                &pileup,
                PIPELINE_MIN_BASEQ,
                PIPELINE_MIN_MAPQ,
                PIPELINE_MIN_TAIL_QUALITY,
                config.exclude_supplementary,
                config.dont_use_soft_clipped_bases,
                Some(&row_context),
            );
            let (indel_evidence, indel_rows) = pileup_indel_evidence_with_rows(
                &pileup,
                PIPELINE_MIN_BASEQ,
                PIPELINE_MIN_MAPQ,
                PIPELINE_MIN_TAIL_QUALITY,
                config.exclude_supplementary,
                config.dont_use_soft_clipped_bases,
                Some(&row_context),
            );
            output.read_observations.extend(snp_rows);
            output.read_observations.extend(indel_rows);

            let ref_index = base_index(ref_base);
            let snp_alt = best_snp_alt(ref_index, &snp_evidence);
            let indel_alt = best_indel_alt(&indel_evidence);
            let snp_alt_count = snp_alt.map(|(_, count)| count).unwrap_or(0);
            let indel_alt_count = indel_alt.map(|(_, count)| *count).unwrap_or(0);
            let depth = snp_evidence.counts.depth.max(indel_evidence.counts.depth);
            let best_alt_count = snp_alt_count.max(indel_alt_count);
            let alt_fraction = if depth == 0 {
                0.0
            } else {
                f64::from(best_alt_count) / f64::from(depth)
            };
            let is_active = is_active_locus(ref_index, &snp_evidence, depth).0
                || is_active_indel(&indel_evidence).0;
            active_region.observed_loci += 1;
            active_region.max_alt_fraction = active_region.max_alt_fraction.max(alt_fraction);
            active_region.mean_alt_fraction += alt_fraction;

            if is_active {
                active_region.active_loci += 1;
                output.active_loci.push(ReplayActiveLocusRow {
                    contig: region_interval.contig.clone(),
                    pos: pos1,
                    region: region.clone(),
                    ref_base,
                    depth,
                    snp_alt_count,
                    snp_best_alt: snp_alt
                        .map(|(base_index, _)| base_from_index(base_index) as char)
                        .map(|base| base.to_string())
                        .unwrap_or_default(),
                    indel_alt_count,
                    indel_best_alt: indel_alt
                        .as_ref()
                        .map(|(allele, _)| indel_allele_label(allele))
                        .unwrap_or_default(),
                    alt_fraction,
                    active_probability_proxy: is_active as u8 as f64,
                });
            }

            // We removed best_snp_call and best_indel_call here
            // since we will use assemble_haplotypes below.
        }

        let call_active_loci_rows = collect_call_active_loci_rows(
            config,
            tid,
            window,
            active_span,
            region_interval,
            &ref_bases,
            bam,
        )?;

        let ref_region_start_offset = (region_interval.start - window.start) as usize;
        let ref_region_end_offset = (region_interval.end - window.start) as usize;
        let local_ref_bases = &ref_bases[ref_region_start_offset..=ref_region_end_offset];

        bam.fetch((
            tid as i32,
            (region_interval.start - 1) as i64,
            region_interval.end as i64,
        ))
        .with_context(|| {
            format!(
                "fetching BAM region for PairHMM: {}:{}-{}",
                region_interval.contig, region_interval.start, region_interval.end
            )
        })?;

        let mut read_bases_list = Vec::new();
        let mut read_quals_list = Vec::new();
        let mut read_ins_quals_list = Vec::new();
        let mut read_del_quals_list = Vec::new();
        let mut read_use_for_assembly_list = Vec::new();
        let mut read_assembly_segments_list = Vec::new();
        let mut reads_by_start = std::collections::HashMap::new();
        let mut read_names_list = Vec::new();
        let mut read_is_reverse_list = Vec::new();
        let mut read_ref_spans = Vec::new();
        let mut mapq_list = Vec::new();
        let mut unclipped_loc_list = Vec::new();
        let mut cigar_string_list = Vec::new();

        for r in bam.records() {
            let record = r?;
            if !read_passes_hc_filter(&record, PIPELINE_MIN_MAPQ, config.exclude_supplementary) {
                continue;
            }
            let Some(prepared_read) = prepare_hmm_read(
                &record,
                PIPELINE_MIN_TAIL_QUALITY,
                config.dont_use_soft_clipped_bases,
                region_interval.start,
                region_interval.end,
            ) else {
                continue;
            };

            let start_pos = record.pos();
            let count = reads_by_start.entry(start_pos).or_insert(0);
            let use_for_assembly = *count < 50;
            if use_for_assembly {
                *count += 1;
            }

            read_assembly_segments_list.push(prepared_read.assembly_segments);
            read_bases_list.push(prepared_read.bases);
            read_quals_list.push(prepared_read.quals);
            read_ins_quals_list.push(prepared_read.ins_quals);
            read_del_quals_list.push(prepared_read.del_quals);
            read_use_for_assembly_list.push(use_for_assembly);
            read_names_list.push(String::from_utf8_lossy(record.qname()).into_owned());
            read_is_reverse_list.push(record.is_reverse());
            mapq_list.push(record.mapq());
            unclipped_loc_list.push((record.pos() + 1) as u64);
            cigar_string_list.push(record.cigar().to_string());
            read_ref_spans.push(prepared_read.ref_span);
        }

        let assembly_reads =
            select_assembly_reads(&read_assembly_segments_list, &read_use_for_assembly_list);
        let assembly_reads_owned: Vec<Vec<u8>> =
            assembly_reads.iter().map(|r| (*r).clone()).collect();
        let max_mnp_distance = 0; // Default GATK
        let (mut local_haplotypes, mut valid_events) = assemble_haplotypes(
            &region_interval.contig,
            region_interval.start,
            local_ref_bases,
            &assembly_reads_owned,
            &[10, 25],
            max_mnp_distance,
        );
        supplement_missing_pileup_events(
            &region_interval.contig,
            region_interval.start,
            local_ref_bases,
            &call_active_loci_rows,
            config.standard_min_confidence_threshold_for_calling,
            &mut local_haplotypes,
            &mut valid_events,
        );
        filter_non_acgt_haplotypes_for_single_snp_region(&mut local_haplotypes, &valid_events);

        for call in &valid_events {
            active_region.candidate_events += 1;
            output.events.push(replay_event_row(&region, call)?);
        }

        for (hap_idx, hap) in local_haplotypes.iter().enumerate() {
            output.haplotypes.push(ReplayHaplotypeRow {
                region: region.clone(),
                stage: "assembled",
                haplotype: hap_idx,
                span_start: region_interval.start,
                span_end: region_interval.end,
                kmer: 0,
                length: hap.bases.len() as u32,
                cigar: hap.cigar.clone(),
                is_ref: hap.is_ref,
                bases: String::from_utf8_lossy(&hap.bases).into_owned(),
            });
        }

        let n_reads = read_bases_list.len();
        let mut read_haplotype_likelihoods: Vec<Vec<f64>> = Vec::with_capacity(n_reads);
        for i in 0..n_reads {
            let read_bases = &read_bases_list[i];
            let read_quals = &read_quals_list[i];
            let read_ins_quals = &read_ins_quals_list[i];
            let read_del_quals = &read_del_quals_list[i];
            let gcp = 10;
            let mut hap_likelihoods = Vec::with_capacity(local_haplotypes.len());

            for (hap_idx, hap) in local_haplotypes.iter().enumerate() {
                let score = pair_hmm::compute_read_likelihood_given_haplotype(
                    &hap.bases,
                    read_bases,
                    read_quals,
                    read_ins_quals,
                    read_del_quals,
                    gcp,
                );

                output.pairhmms.push(ReplayPairHmmRow {
                    region: region.clone(),
                    read: read_names_list[i].clone(),
                    haplotype: hap_idx,
                    read_index: i,
                    cigar: cigar_string_list[i].clone(),
                    mapq: mapq_list[i],
                    loc: unclipped_loc_list[i],
                    unclipped_loc: unclipped_loc_list[i],
                    length: read_bases.len() as u32,
                    score,
                });
                hap_likelihoods.push(score);
            }
            read_haplotype_likelihoods.push(hap_likelihoods);
        }

        let mut final_calls = genotype_assembled_events(
            &local_haplotypes,
            &valid_events,
            &read_haplotype_likelihoods,
            &read_is_reverse_list,
            &read_ref_spans,
            config.standard_min_confidence_threshold_for_calling,
        );
        prune_unsupported_simple_snp_calls_in_dense_clusters(
            &mut final_calls,
            &call_active_loci_rows,
        );
        merge_missing_strong_snp_cluster_rescues_from_pileup(
            &mut final_calls,
            &region_interval.contig,
            region_interval.start,
            local_ref_bases,
            &call_active_loci_rows,
            &valid_events,
            config.standard_min_confidence_threshold_for_calling,
        );
        if final_calls.is_empty() {
            final_calls = rescue_collapsed_strong_snp_cluster_from_pileup(
                &region_interval.contig,
                region_interval.start,
                local_ref_bases,
                &call_active_loci_rows,
                &valid_events,
                config.standard_min_confidence_threshold_for_calling,
            );
        }
        output.variants.extend(final_calls);

        for (i, hap) in local_haplotypes.into_iter().enumerate() {
            output.haplotypes.push(ReplayHaplotypeRow {
                region: region.clone(),
                stage: "unclipped",
                haplotype: i,
                span_start: region_interval.start,
                span_end: region_interval.end,
                kmer: 0,
                length: hap.bases.len() as u32,
                cigar: hap.cigar,
                is_ref: hap.is_ref,
                bases: String::from_utf8_lossy(&hap.bases).into_owned(),
            });
        }
        if active_region.observed_loci > 0 {
            active_region.mean_alt_fraction /= active_region.observed_loci as f64;
        }
        output.active_regions.push(active_region);
    }
    Ok(())
}

