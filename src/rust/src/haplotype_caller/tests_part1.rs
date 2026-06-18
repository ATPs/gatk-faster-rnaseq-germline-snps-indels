    use super::*;

    fn test_bam_record(flags: u16, cigar: Vec<rust_htslib::bam::record::Cigar>) -> bam::Record {
        let cigar = rust_htslib::bam::record::CigarString(cigar);
        let read_len = cigar
            .iter()
            .map(|op| match op {
                rust_htslib::bam::record::Cigar::Match(len)
                | rust_htslib::bam::record::Cigar::Equal(len)
                | rust_htslib::bam::record::Cigar::Diff(len)
                | rust_htslib::bam::record::Cigar::Ins(len)
                | rust_htslib::bam::record::Cigar::SoftClip(len) => *len as usize,
                _ => 0,
            })
            .sum::<usize>();
        let bases = vec![b'A'; read_len];
        let quals = vec![30_u8; read_len];
        let mut record = bam::Record::new();
        record.set(b"read1", Some(&cigar), &bases, &quals);
        record.set_flags(flags);
        record.set_mapq(60);
        record
    }

    fn test_bam_record_with_bases_quals(
        cigar: Vec<rust_htslib::bam::record::Cigar>,
        bases: &[u8],
        quals: &[u8],
        pos0: i64,
    ) -> bam::Record {
        let cigar = rust_htslib::bam::record::CigarString(cigar);
        let mut record = bam::Record::new();
        record.set(b"read1", Some(&cigar), bases, quals);
        record.set_flags(0);
        record.set_mapq(60);
        record.set_pos(pos0);
        record
    }

    fn test_dict() -> SequenceDict {
        let lines = vec![
            "@HD\tVN:1.6\tSO:coordinate".to_string(),
            "@SQ\tSN:chr2\tLN:100".to_string(),
            "@SQ\tSN:chr1\tLN:200".to_string(),
        ];
        parse_dict_lines(&lines, Path::new("test.interval_list")).unwrap()
    }

    fn snp_evidence(
        ref_index: usize,
        ref_count: u32,
        alt_index: usize,
        alt_count: u32,
        quality: u8,
    ) -> SnpEvidence {
        let mut evidence = SnpEvidence::default();
        for _ in 0..ref_count {
            evidence.counts.counts[ref_index] += 1;
            evidence.counts.depth += 1;
            evidence.strands[ref_index].increment(false);
            evidence.observations.push(BaseObservation {
                base_index: ref_index,
                quality,
                is_reverse: false,
            });
            evidence.active_observations.push(ActiveBaseObservation {
                base_index: Some(ref_index),
                quality,
            });
        }
        for idx in 0..alt_count {
            let is_reverse = idx % 2 == 1;
            evidence.counts.counts[alt_index] += 1;
            evidence.counts.depth += 1;
            evidence.strands[alt_index].increment(is_reverse);
            evidence.observations.push(BaseObservation {
                base_index: alt_index,
                quality,
                is_reverse,
            });
            evidence.active_observations.push(ActiveBaseObservation {
                base_index: Some(alt_index),
                quality,
            });
        }
        evidence
    }

    fn indel_evidence(
        ref_count: u32,
        alt_allele: IndelAllele,
        alt_count: u32,
        quality: u8,
    ) -> IndelEvidence {
        let mut evidence = IndelEvidence::default();
        evidence.counts.ref_count = ref_count;
        evidence.counts.depth = ref_count + alt_count;
        evidence.counts.counts.insert(alt_allele.clone(), alt_count);
        for _ in 0..ref_count {
            evidence.ref_strand.increment(false);
            evidence.observations.push(IndelObservation {
                allele: IndelObservationAllele::Ref,
                quality,
                is_reverse: false,
            });
        }
        for idx in 0..alt_count {
            let is_reverse = idx % 2 == 1;
            evidence
                .alt_strands
                .entry(alt_allele.clone())
                .or_default()
                .increment(is_reverse);
            evidence.observations.push(IndelObservation {
                allele: IndelObservationAllele::Alt(alt_allele.clone()),
                quality,
                is_reverse,
            });
        }
        evidence
    }

    fn test_variant(contig: &str, pos: u64, ref_allele: &[u8], alt_allele: &[u8]) -> VariantCall {
        VariantCall {
            contig: contig.to_string(),
            pos,
            id: None,
            db: false,
            ref_allele: ref_allele.to_vec(),
            alt_allele: alt_allele.to_vec(),
            depth: 10,
            ref_count: 8,
            alt_count: 2,
            qual: 20,
            fs: 0.0,
            pl: [20, 0, 20],
            genotype_index: 1,
        }
    }

    #[test]
    fn intervals_sort_by_dictionary_order() {
        let dict = test_dict();
        let mut intervals = vec![
            Interval {
                contig: "chr1".to_string(),
                start: 50,
                end: 80,
            },
            Interval {
                contig: "chr2".to_string(),
                start: 5,
                end: 10,
            },
            Interval {
                contig: "chr1".to_string(),
                start: 81,
                end: 90,
            },
        ];
        sort_intervals(&mut intervals, &dict).unwrap();
        assert_eq!(
            intervals,
            vec![
                Interval {
                    contig: "chr2".to_string(),
                    start: 5,
                    end: 10,
                },
                Interval {
                    contig: "chr1".to_string(),
                    start: 50,
                    end: 80,
                },
                Interval {
                    contig: "chr1".to_string(),
                    start: 81,
                    end: 90,
                },
            ]
        );
    }

    #[test]
    fn fetch_windows_coalesce_nearby_intervals() {
        let dict = parse_dict_lines(
            &[
                "@HD\tVN:1.6\tSO:coordinate".to_string(),
                "@SQ\tSN:chr1\tLN:5000".to_string(),
            ],
            Path::new("test.interval_list"),
        )
        .unwrap();
        let intervals = vec![
            Interval {
                contig: "chr1".to_string(),
                start: 1,
                end: 10,
            },
            Interval {
                contig: "chr1".to_string(),
                start: 11,
                end: 30,
            },
            Interval {
                contig: "chr1".to_string(),
                start: 2_000,
                end: 2_010,
            },
        ];
        let windows = coalesce_fetch_windows(&intervals, &dict);
        assert_eq!(windows.len(), 2);
        assert_eq!(windows[0].start, 1);
        assert_eq!(windows[0].end, 130);
        assert_eq!(windows[0].intervals.len(), 2);
        assert_eq!(windows[1].start, 1_900);
        assert_eq!(windows[1].end, 2_110);
    }
    #[test]
    fn fetch_windows_pad_requested_intervals_and_clip_to_contig_edges() {
        let dict = test_dict();
        let left_windows = coalesce_fetch_windows(
            &[Interval {
                contig: "chr1".to_string(),
                start: 10,
                end: 20,
            }],
            &dict,
        );
        assert_eq!(left_windows.len(), 1);
        assert_eq!(left_windows[0].start, 1);
        assert_eq!(left_windows[0].end, 120);

        let right_windows = coalesce_fetch_windows(
            &[Interval {
                contig: "chr1".to_string(),
                start: 181,
                end: 200,
            }],
            &dict,
        );
        assert_eq!(right_windows.len(), 1);
        assert_eq!(right_windows[0].start, 81);
        assert_eq!(right_windows[0].end, 200);
    }

    #[test]
    fn fetch_window_partitioning_preserves_all_bases() {
        let windows = vec![
            FetchWindow {
                contig: "chr1".to_string(),
                start: 1,
                end: 10,
                intervals: Vec::new(),
            },
            FetchWindow {
                contig: "chr1".to_string(),
                start: 11,
                end: 30,
                intervals: Vec::new(),
            },
            FetchWindow {
                contig: "chr1".to_string(),
                start: 31,
                end: 60,
                intervals: Vec::new(),
            },
        ];
        let total_bases: u64 = windows.iter().map(FetchWindow::len).sum();
        let partitions = partition_fetch_windows_by_bases(&windows, 2);
        let partition_bases: u64 = partitions.iter().flatten().map(FetchWindow::len).sum();
        assert_eq!(partition_bases, total_bases);
        assert_eq!(
            partitions.iter().map(Vec::len).sum::<usize>(),
            windows.len()
        );
    }

    #[test]
    fn requested_position_uses_sorted_interval_cursor() {
        let intervals = vec![
            Interval {
                contig: "chr1".to_string(),
                start: 10,
                end: 20,
            },
            Interval {
                contig: "chr1".to_string(),
                start: 30,
                end: 40,
            },
        ];
        let mut cursor = 0;
        assert!(!position_is_requested(&intervals, 9, &mut cursor));
        assert!(position_is_requested(&intervals, 10, &mut cursor));
        assert!(position_is_requested(&intervals, 20, &mut cursor));
        assert!(!position_is_requested(&intervals, 25, &mut cursor));
        assert!(position_is_requested(&intervals, 35, &mut cursor));
    }

    #[test]
    fn select_assembly_reads_keeps_only_flagged_segments_in_input_order() {
        let reads = vec![
            vec![b"AAA".to_vec(), b"TTT".to_vec()],
            vec![b"CCC".to_vec()],
            vec![b"GGG".to_vec()],
        ];
        let selected = select_assembly_reads(&reads, &[true, false, true]);
        let selected_owned: Vec<Vec<u8>> = selected.into_iter().cloned().collect();
        assert_eq!(
            selected_owned,
            vec![b"AAA".to_vec(), b"TTT".to_vec(), b"GGG".to_vec()]
        );
    }

    #[test]
    fn assembly_read_segments_split_on_n_and_low_base_quality() {
        let segments = assembly_read_segments(
            b"AAAANCCCTGGG",
            &[30, 30, 30, 30, 30, 30, 30, 30, 9, 30, 30, 30],
        );

        assert_eq!(
            segments,
            vec![b"AAAA".to_vec(), b"CCC".to_vec(), b"GGG".to_vec()]
        );
    }

    #[test]
    fn allele_marginalization_uses_max_across_supporting_haplotypes() {
        let marginalized = marginalize_allele_likelihoods(&[-2.0, -2.0]);
        assert_eq!(marginalized, -2.0);
    }

    #[test]
    fn allele_marginalization_returns_negative_infinity_for_empty_support() {
        assert_eq!(marginalize_allele_likelihoods(&[]), f64::NEG_INFINITY);
    }

    #[test]
    fn pair_hmm_base_quality_caps_by_mapq_then_squashes_low_values() {
        assert_eq!(pair_hmm_base_quality(30, 40), 30);
        assert_eq!(pair_hmm_base_quality(30, 20), 20);
        assert_eq!(pair_hmm_base_quality(17, 40), PAIR_HMM_MIN_USABLE_Q_SCORE);
        assert_eq!(pair_hmm_base_quality(30, 10), PAIR_HMM_MIN_USABLE_Q_SCORE);
        assert_eq!(pair_hmm_base_quality(18, 40), 18);
    }

    #[test]
    fn pair_hmm_indel_open_quality_uses_min_usable_floor() {
        assert_eq!(pair_hmm_indel_open_quality(4), PAIR_HMM_MIN_USABLE_Q_SCORE);
        assert_eq!(pair_hmm_indel_open_quality(6), 6);
        assert_eq!(pair_hmm_indel_open_quality(25), 25);
    }

    #[test]
    fn pad_active_region_clips_to_fetch_window_and_preserves_active_span() {
        let window = FetchWindow {
            contig: "chr1".to_string(),
            start: 100,
            end: 200,
            intervals: vec![Interval {
                contig: "chr1".to_string(),
                start: 100,
                end: 200,
            }],
        };
        let interval = Interval {
            contig: "chr1".to_string(),
            start: 110,
            end: 150,
        };
        let padded = pad_active_region(&interval, &window);
        assert_eq!(padded.start, 100);
        assert_eq!(padded.end, 200);
        assert_eq!(interval.start, 110);
        assert_eq!(interval.end, 150);
    }

    #[test]
    fn genotype_assembled_events_emits_final_calls_from_pairhmm_matrix() {
        let local_haplotypes = vec![
            LocalHaplotype {
                bases: b"AAA".to_vec(),
                is_ref: true,
                cigar: "3M".to_string(),
                event_indices: vec![],
            },
            LocalHaplotype {
                bases: b"ACA".to_vec(),
                is_ref: false,
                cigar: "3M".to_string(),
                event_indices: vec![0],
            },
        ];
        let valid_events = vec![test_variant("chr1", 10, b"A", b"C")];
        let read_haplotype_likelihoods = vec![vec![-10.0, 0.0], vec![-10.0, 0.0], vec![-10.0, 0.0]];
        let read_is_reverse_list = vec![false, true, false];
        let read_ref_spans = vec![(10, 10), (10, 10), (10, 10)];

        let final_calls = genotype_assembled_events(
            &local_haplotypes,
            &valid_events,
            &read_haplotype_likelihoods,
            &read_is_reverse_list,
            &read_ref_spans,
            20.0,
        );

        assert_eq!(final_calls.len(), 1);
        let call = &final_calls[0];
        assert_eq!(call.genotype(), "1/1");
        assert_eq!(call.depth, 3);
        assert_eq!(call.ref_count, 0);
        assert_eq!(call.alt_count, 3);
        assert!(call.qual >= 20);
    }

    #[test]
    fn genotype_assembled_events_ignores_overlapping_alt_haplotypes_as_ref_evidence() {
        let local_haplotypes = vec![
            LocalHaplotype {
                bases: b"AAA".to_vec(),
                is_ref: true,
                cigar: "3M".to_string(),
                event_indices: vec![],
            },
            LocalHaplotype {
                bases: b"ACA".to_vec(),
                is_ref: false,
                cigar: "3M".to_string(),
                event_indices: vec![0],
            },
            LocalHaplotype {
                bases: b"AGA".to_vec(),
                is_ref: false,
                cigar: "3M".to_string(),
                event_indices: vec![1],
            },
        ];
        let valid_events = vec![
            test_variant("chr1", 10, b"A", b"C"),
            test_variant("chr1", 10, b"A", b"G"),
        ];
        let mut read_haplotype_likelihoods = Vec::new();
        read_haplotype_likelihoods.extend(vec![vec![-10.0, 0.0, -10.0]; 2]);
        read_haplotype_likelihoods.extend(vec![vec![-10.0, -10.0, 0.0]; 8]);
        let read_is_reverse_list = vec![false; read_haplotype_likelihoods.len()];
        let read_ref_spans = vec![(10, 10); read_haplotype_likelihoods.len()];

        let final_calls = genotype_assembled_events(
            &local_haplotypes,
            &valid_events,
            &read_haplotype_likelihoods,
            &read_is_reverse_list,
            &read_ref_spans,
            0.0,
        );

        let c_call = final_calls
            .iter()
            .find(|call| call.ref_allele == b"A" && call.alt_allele == b"C")
            .unwrap();
        let g_call = final_calls
            .iter()
            .find(|call| call.ref_allele == b"A" && call.alt_allele == b"G")
            .unwrap();
        assert_ne!(c_call.genotype(), "0/0");
        assert_ne!(g_call.genotype(), "0/0");
    }

    #[test]
    fn genotype_assembled_events_pair_genotyping_ignores_competing_overlap_pairs_as_ref() {
        let local_haplotypes = vec![
            LocalHaplotype {
                bases: b"AAAA".to_vec(),
                is_ref: true,
                cigar: "4M".to_string(),
                event_indices: vec![],
            },
            LocalHaplotype {
                bases: b"ACAA".to_vec(),
                is_ref: false,
                cigar: "4M".to_string(),
                event_indices: vec![0],
            },
            LocalHaplotype {
                bases: b"AGAA".to_vec(),
                is_ref: false,
                cigar: "4M".to_string(),
                event_indices: vec![1],
            },
            LocalHaplotype {
                bases: b"ACCA".to_vec(),
                is_ref: false,
                cigar: "4M".to_string(),
                event_indices: vec![0],
            },
            LocalHaplotype {
                bases: b"ACGA".to_vec(),
                is_ref: false,
                cigar: "4M".to_string(),
                event_indices: vec![0],
            },
        ];
        let valid_events = vec![
            test_variant("chr1", 10, b"A", b"C"),
            test_variant("chr1", 10, b"A", b"G"),
        ];
        let read_haplotype_likelihoods = vec![
            vec![-10.0, 0.0, -12.0, -1.0, -1.0],
            vec![-10.0, -12.0, 0.0, -12.0, -12.0],
            vec![0.0, -10.0, -10.0, -10.0, -10.0],
        ];
        let read_is_reverse_list = vec![false; read_haplotype_likelihoods.len()];
        let read_ref_spans = vec![(10, 10); read_haplotype_likelihoods.len()];

        let final_calls = genotype_assembled_events(
            &local_haplotypes,
            &valid_events,
            &read_haplotype_likelihoods,
            &read_is_reverse_list,
            &read_ref_spans,
            0.0,
        );

        let c_call = final_calls
            .iter()
            .find(|call| call.ref_allele == b"A" && call.alt_allele == b"C")
            .unwrap();
        let g_call = final_calls
            .iter()
            .find(|call| call.ref_allele == b"A" && call.alt_allele == b"G")
            .unwrap();
        assert_eq!(c_call.genotype(), "0/1");
        assert_eq!(g_call.genotype(), "0/1");
    }

    #[test]
    fn genotype_assembled_events_uses_isolated_event_counts_for_simple_snp() {
        let local_haplotypes = vec![
            LocalHaplotype {
                bases: b"AAAA".to_vec(),
                is_ref: true,
                cigar: "4M".to_string(),
                event_indices: vec![],
            },
            LocalHaplotype {
                bases: b"ACAA".to_vec(),
                is_ref: false,
                cigar: "4M".to_string(),
                event_indices: vec![0],
            },
            LocalHaplotype {
                bases: b"AAGA".to_vec(),
                is_ref: false,
                cigar: "4M".to_string(),
                event_indices: vec![1],
            },
        ];
        let valid_events = vec![
            test_variant("chr1", 10, b"A", b"C"),
            test_variant("chr1", 12, b"A", b"G"),
        ];
        let read_haplotype_likelihoods = vec![
            vec![-10.0, 0.0, -12.0],
            vec![-10.0, 0.0, -12.0],
            vec![-10.0, -12.0, 0.0],
            vec![-10.0, -12.0, 0.0],
        ];
        let read_is_reverse_list = vec![false, true, false, true];
        let read_ref_spans = vec![(10, 10), (10, 10), (12, 12), (12, 12)];

        let final_calls = genotype_assembled_events(
            &local_haplotypes,
            &valid_events,
            &read_haplotype_likelihoods,
            &read_is_reverse_list,
            &read_ref_spans,
            0.0,
        );

        let c_call = final_calls
            .iter()
            .find(|call| call.ref_allele == b"A" && call.alt_allele == b"C")
            .unwrap();
        let g_call = final_calls
            .iter()
            .find(|call| call.ref_allele == b"A" && call.alt_allele == b"G")
            .unwrap();
        assert_eq!(c_call.depth, 2);
        assert_eq!(c_call.alt_count, 2);
        assert_eq!(g_call.depth, 2);
        assert_eq!(g_call.alt_count, 2);
    }

    #[test]
    fn genotype_assembled_events_ignores_reads_that_do_not_span_the_event() {
        let local_haplotypes = vec![
            LocalHaplotype {
                bases: b"AAAA".to_vec(),
                is_ref: true,
                cigar: "4M".to_string(),
                event_indices: vec![],
            },
            LocalHaplotype {
                bases: b"ACGA".to_vec(),
                is_ref: false,
                cigar: "4M".to_string(),
                event_indices: vec![0, 1],
            },
        ];
        let valid_events = vec![
            test_variant("chr1", 10, b"A", b"C"),
            test_variant("chr1", 12, b"A", b"G"),
        ];
        let read_haplotype_likelihoods = vec![vec![-10.0, 0.0], vec![-10.0, 0.0], vec![-10.0, 0.0]];
        let read_is_reverse_list = vec![false, false, true];
        let read_ref_spans = vec![(10, 10), (12, 12), (12, 12)];

        let final_calls = genotype_assembled_events(
            &local_haplotypes,
            &valid_events,
            &read_haplotype_likelihoods,
            &read_is_reverse_list,
            &read_ref_spans,
            0.0,
        );

        let c_call = final_calls
            .iter()
            .find(|call| call.ref_allele == b"A" && call.alt_allele == b"C")
            .unwrap();
        let g_call = final_calls
            .iter()
            .find(|call| call.ref_allele == b"A" && call.alt_allele == b"G")
            .unwrap();
        assert_eq!(c_call.depth, 1);
        assert_eq!(c_call.alt_count, 1);
        assert_eq!(g_call.depth, 2);
        assert_eq!(g_call.alt_count, 2);
    }

    #[test]
    fn genotype_assembled_events_marks_same_position_events_for_pair_genotyping() {
        let valid_events = vec![
            test_variant("chr1", 10, b"A", b"C"),
            test_variant("chr1", 10, b"A", b"G"),
            test_variant("chr1", 12, b"T", b"C"),
        ];

        assert_eq!(
            overlapping_event_mask(&valid_events),
            vec![true, true, false]
        );
    }

    #[test]
    fn read_reference_span_from_start_and_cigar_uses_reference_consuming_ops_only() {
        assert_eq!(
            read_reference_span_from_start_and_cigar(101, "10M1I10M"),
            (101, 120)
        );
        assert_eq!(
            read_reference_span_from_start_and_cigar(101, "10M100N10M1D5M"),
            (101, 226)
        );
    }

    #[test]
    fn prepare_hmm_read_trims_low_quality_tails_and_updates_ref_span() {
        let record = test_bam_record_with_bases_quals(
            vec![rust_htslib::bam::record::Cigar::Match(12)],
            b"AACCGGTTAACC",
            &[5, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 5],
            99,
        );

        let prepared = prepare_hmm_read(&record, 10, true, 100, 111).unwrap();

        assert_eq!(prepared.bases, b"ACCGGTTAAC");
        assert_eq!(prepared.ref_span, (101, 110));
        assert_eq!(prepared.quals, vec![30; 10]);
        assert_eq!(prepared.assembly_segments, vec![b"ACCGGTTAAC".to_vec()]);
    }

    #[test]
    fn prepare_hmm_read_excludes_soft_clips_from_bases_and_ref_span() {
        let record = test_bam_record_with_bases_quals(
            vec![
                rust_htslib::bam::record::Cigar::SoftClip(2),
                rust_htslib::bam::record::Cigar::Match(10),
                rust_htslib::bam::record::Cigar::SoftClip(2),
            ],
            b"TTACCGGTTAACAA",
            &[30; 14],
            199,
        );

        let prepared = prepare_hmm_read(&record, 10, true, 200, 209).unwrap();

        assert_eq!(prepared.bases, b"ACCGGTTAAC");
        assert_eq!(prepared.ref_span, (200, 209));
        assert_eq!(prepared.quals, vec![30; 10]);
    }

    #[test]
    fn prepare_hmm_read_clips_match_bases_to_region_span() {
        let record = test_bam_record_with_bases_quals(
            vec![rust_htslib::bam::record::Cigar::Match(14)],
            b"AACCGGTTAACCGG",
            &[30; 14],
            99,
        );

        let prepared = prepare_hmm_read(&record, 10, true, 103, 112).unwrap();

        assert_eq!(prepared.bases, b"CGGTTAACCG");
        assert_eq!(prepared.ref_span, (103, 112));
        assert_eq!(prepared.quals, vec![30; 10]);
    }

    #[test]
    fn prepare_hmm_read_keeps_boundary_insertion_when_clipping_left_tail() {
        let record = test_bam_record_with_bases_quals(
            vec![
                rust_htslib::bam::record::Cigar::Match(3),
                rust_htslib::bam::record::Cigar::Ins(2),
                rust_htslib::bam::record::Cigar::Match(8),
            ],
            b"AAAGGTTTTTTTT",
            &[30; 13],
            99,
        );

        let prepared = prepare_hmm_read(&record, 10, true, 103, 110).unwrap();

        assert_eq!(prepared.bases, b"GGTTTTTTTT");
        assert_eq!(prepared.ref_span, (103, 110));
        assert_eq!(prepared.quals, vec![30; 10]);
    }

    #[test]
    fn hc_filter_includes_supplementary_reads_by_default() {
        let record = test_bam_record(0x800, vec![rust_htslib::bam::record::Cigar::Match(5)]);
        assert!(read_passes_hc_filter(&record, 20, false));
    }

    #[test]
    fn hc_filter_can_exclude_supplementary_reads_for_debugging() {
        let record = test_bam_record(0x800, vec![rust_htslib::bam::record::Cigar::Match(5)]);
        assert!(!read_passes_hc_filter(&record, 20, true));
    }

    #[test]
    fn hc_filter_rejects_reads_without_reference_consuming_cigar_ops() {
        let record = test_bam_record(0, vec![rust_htslib::bam::record::Cigar::SoftClip(5)]);
        assert!(!read_passes_hc_filter(&record, 20, false));
    }

    #[test]
    fn snp_call_uses_likelihood_quality_threshold() {
        let evidence = snp_evidence(0, 5, 1, 5, 30);
        assert!(best_snp_call("chr1", 10, b'A', evidence.clone(), 200.0).is_none());
        let call = best_snp_call("chr1", 10, b'A', evidence, 20.0).unwrap();
        assert_eq!(call.ref_allele, b"A");
        assert_eq!(call.alt_allele, b"C");
        assert!(call.qual >= 90);
        assert_eq!(call.genotype(), "0/1");
        assert_eq!(call.pl[1], 0);
    }

    #[test]
    fn snp_call_uses_hom_alt_when_likelihoods_overcome_prior() {
        let evidence = snp_evidence(0, 0, 1, 20, 30);
        let call = best_snp_call("chr1", 10, b'A', evidence, 20.0).unwrap();
        assert_eq!(call.genotype(), "1/1");
        assert_eq!(call.alt_allele_count(), 2);
        assert_eq!(call.pl[2], 0);
    }

    #[test]
    fn snp_call_rejects_low_alt_fraction_noise_after_likelihoods() {
        let evidence = snp_evidence(0, 29, 1, 2, 30);
        assert!(best_snp_call("chr1", 10, b'A', evidence, 20.0).is_none());
    }

    #[test]
    fn active_locus_treats_non_acgt_bases_as_non_ref_only_for_discovery() {
        let mut evidence = SnpEvidence::default();
        evidence.active_observations.push(ActiveBaseObservation {
            base_index: None,
            quality: 30,
        });
        evidence.active_observations.push(ActiveBaseObservation {
            base_index: None,
            quality: 30,
        });

        let (active, qual) = is_active_locus(Some(0), &evidence, evidence.counts.depth);
        assert!(active);
        assert!(qual >= ACTIVE_REGION_DISCOVERY_CONFIDENCE as u32);
        assert_eq!(evidence.counts.depth, 0);
        assert!(best_snp_alt(Some(0), &evidence).is_none());
        assert!(best_snp_call("chr1", 10, b'A', evidence, 20.0).is_none());
    }

    #[test]
    fn overlapping_same_base_fragment_caps_both_observations() {
        let observations = vec![
            BaseObservation {
                base_index: 1,
                quality: 35,
                is_reverse: false,
            },
            BaseObservation {
                base_index: 1,
                quality: 33,
                is_reverse: true,
            },
        ];
        let kept = adjust_fragment_base_observations(&observations);
        assert_eq!(kept.len(), 2);
        assert!(kept.iter().all(|observation| observation.base_index == 1));
        assert!(kept
            .iter()
            .all(|observation| observation.quality == HALF_DEFAULT_PCR_SNV_QUAL));
    }

    #[test]
    fn overlapping_discordant_fragment_zeroes_both_observations() {
        let observations = vec![
            BaseObservation {
                base_index: 1,
                quality: 35,
                is_reverse: false,
            },
            BaseObservation {
                base_index: 2,
                quality: 33,
                is_reverse: true,
            },
        ];
        let kept = adjust_fragment_base_observations(&observations);
        assert_eq!(kept.len(), 2);
        assert!(kept.iter().all(|observation| observation.quality == 0));
    }

    #[test]
    fn fisher_strand_score_detects_one_sided_alt_support() {
        let fs = fisher_strand_score(
            StrandCounts {
                forward: 10,
                reverse: 0,
            },
            StrandCounts {
                forward: 0,
                reverse: 10,
            },
        );
        assert!(fs > 30.0);
    }

    #[test]
    fn insertion_call_uses_left_anchor() {
        let evidence = indel_evidence(6, IndelAllele::Insertion(b"TG".to_vec()), 6, 30);
        let call = best_indel_call("chr1", 100, 100, b"ACGT", evidence, 20.0).unwrap();
        assert_eq!(call.pos, 100);
        assert_eq!(call.ref_allele, b"A");
        assert_eq!(call.alt_allele, b"ATG");
        assert_eq!(call.ref_count, 6);
        assert_eq!(call.alt_count, 6);
    }

    #[test]
    fn deletion_call_uses_left_anchor_and_deleted_reference() {
        let evidence = indel_evidence(6, IndelAllele::Deletion(2), 6, 30);
        let call = best_indel_call("chr1", 100, 100, b"ACGT", evidence, 20.0).unwrap();
        assert_eq!(call.pos, 100);
        assert_eq!(call.ref_allele, b"ACG");
        assert_eq!(call.alt_allele, b"A");
    }

    #[test]
    fn indel_normalization_left_aligns_homopolymer_insertion() {
        let (pos, ref_allele, alt_allele) =
            left_normalize_indel(102, 100, b"ATTTG", b"T".to_vec(), b"TT".to_vec());
        assert_eq!(pos, 100);
        assert_eq!(ref_allele, b"A");
        assert_eq!(alt_allele, b"AT");
    }

    #[test]
    fn indel_normalization_left_aligns_homopolymer_deletion() {
        let (pos, ref_allele, alt_allele) =
            left_normalize_indel(102, 100, b"ATTTG", b"TT".to_vec(), b"T".to_vec());
        assert_eq!(pos, 100);
        assert_eq!(ref_allele, b"AT");
        assert_eq!(alt_allele, b"A");
    }

    #[test]
    fn replay_event_row_uses_gatk_like_event_key() {
        let call = test_variant("chr1", 100, b"A", b"C");
        let row = replay_event_row("chr1:90-110", &call).unwrap();
        assert_eq!(row.region, "chr1:90-110");
        assert_eq!(row.event, "chr1:100:SNP:A*,C");
        assert_eq!(row.event_type, "SNP");
        assert_eq!(row.alleles, "A*,C");
        assert_eq!(row.gt, "0/1");
    }

    #[test]
    fn dbsnp_record_match_requires_exact_ref_and_alt() {
        let variant = test_variant("chr1", 100, b"A", b"ATG");
        let record = b"chr1\t100\trs1\tA\tC,ATG\t.\t.\t.";
        assert!(dbsnp_record_matches(record, &variant).unwrap());
        assert_eq!(dbsnp_record_id(record).unwrap(), "rs1");
        let non_match = b"chr1\t100\trs2\tA\tC,G\t.\t.\t.";
        assert!(!dbsnp_record_matches(non_match, &variant).unwrap());
    }

    #[test]
    fn variant_calls_sort_by_dictionary_order() {
        let dict = test_dict();
        let mut variants = vec![
            test_variant("chr1", 20, b"A", b"C"),
            test_variant("chr2", 30, b"G", b"T"),
            test_variant("chr2", 10, b"A", b"G"),
        ];
        sort_variant_calls(&mut variants, &dict).unwrap();
        assert_eq!(variants[0].contig, "chr2");
        assert_eq!(variants[0].pos, 10);
        assert_eq!(variants[1].contig, "chr2");
        assert_eq!(variants[1].pos, 30);
        assert_eq!(variants[2].contig, "chr1");
    }

    #[test]
    fn pileup_fallback_collects_strong_snp_candidate() {
        let active_loci = vec![ReplayActiveLocusRow {
            contig: "chr1".to_string(),
            pos: 101,
            region: "chr1:100-102".to_string(),
            ref_base: b'A',
            depth: 16,
            snp_alt_count: 8,
            snp_best_alt: "G".to_string(),
            indel_alt_count: 0,
            indel_best_alt: String::new(),
            alt_fraction: 0.5,
            active_probability_proxy: 1.0,
        }];
        let events = collect_pileup_fallback_events("chr1", 100, b"AAT", &active_loci, 20.0);
        assert_eq!(events.len(), 1);
        let event = &events[0];
        assert_eq!(event.pos, 101);
        assert_eq!(event.ref_allele, b"A");
        assert_eq!(event.alt_allele, b"G");
        assert!(event.qual >= 20);
    }

    #[test]
    fn pileup_fallback_requires_strong_support() {
        let active_loci = vec![
            ReplayActiveLocusRow {
                contig: "chr1".to_string(),
                pos: 101,
                region: "chr1:100-102".to_string(),
                ref_base: b'A',
                depth: 9,
                snp_alt_count: 8,
                snp_best_alt: "G".to_string(),
                indel_alt_count: 0,
                indel_best_alt: String::new(),
                alt_fraction: 0.5,
                active_probability_proxy: 1.0,
            },
            ReplayActiveLocusRow {
                contig: "chr1".to_string(),
                pos: 102,
                region: "chr1:100-102".to_string(),
                ref_base: b'A',
                depth: 16,
                snp_alt_count: 2,
                snp_best_alt: "G".to_string(),
                indel_alt_count: 0,
                indel_best_alt: String::new(),
                alt_fraction: 0.125,
                active_probability_proxy: 1.0,
            },
        ];
        let events = collect_pileup_fallback_events("chr1", 100, b"AAT", &active_loci, 20.0);
        assert!(events.is_empty());
    }
