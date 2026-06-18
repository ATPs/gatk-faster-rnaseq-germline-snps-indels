    #[test]
    fn rescue_collapsed_strong_snp_cluster_from_pileup_recovers_matching_candidates() {
        let valid_events = vec![
            test_variant("chr1", 100, b"A", b"G"),
            test_variant("chr1", 101, b"A", b"C"),
        ];
        let active_loci = vec![
            ReplayActiveLocusRow {
                contig: "chr1".to_string(),
                pos: 100,
                region: "chr1:100-102".to_string(),
                ref_base: b'A',
                depth: 16,
                snp_alt_count: 8,
                snp_best_alt: "G".to_string(),
                indel_alt_count: 0,
                indel_best_alt: String::new(),
                alt_fraction: 0.5,
                active_probability_proxy: 1.0,
            },
            ReplayActiveLocusRow {
                contig: "chr1".to_string(),
                pos: 101,
                region: "chr1:100-102".to_string(),
                ref_base: b'A',
                depth: 16,
                snp_alt_count: 8,
                snp_best_alt: "C".to_string(),
                indel_alt_count: 0,
                indel_best_alt: String::new(),
                alt_fraction: 0.5,
                active_probability_proxy: 1.0,
            },
        ];

        let rescued = rescue_collapsed_strong_snp_cluster_from_pileup(
            "chr1",
            100,
            b"AAT",
            &active_loci,
            &valid_events,
            20.0,
        );

        assert_eq!(rescued.len(), 2);
        assert!(rescued
            .iter()
            .any(|event| event.pos == 100 && event.ref_allele == b"A" && event.alt_allele == b"G"));
        assert!(rescued
            .iter()
            .any(|event| event.pos == 101 && event.ref_allele == b"A" && event.alt_allele == b"C"));
    }

    #[test]
    fn rescue_collapsed_strong_snp_cluster_from_pileup_does_not_recover_low_confidence_single_match(
    ) {
        let valid_events = vec![test_variant("chr1", 100, b"A", b"G")];
        let active_loci = vec![
            ReplayActiveLocusRow {
                contig: "chr1".to_string(),
                pos: 100,
                region: "chr1:100-102".to_string(),
                ref_base: b'A',
                depth: 16,
                snp_alt_count: 4,
                snp_best_alt: "G".to_string(),
                indel_alt_count: 0,
                indel_best_alt: String::new(),
                alt_fraction: 0.25,
                active_probability_proxy: 1.0,
            },
            ReplayActiveLocusRow {
                contig: "chr1".to_string(),
                pos: 101,
                region: "chr1:100-102".to_string(),
                ref_base: b'A',
                depth: 16,
                snp_alt_count: 8,
                snp_best_alt: "C".to_string(),
                indel_alt_count: 0,
                indel_best_alt: String::new(),
                alt_fraction: 0.5,
                active_probability_proxy: 1.0,
            },
        ];

        let rescued = rescue_collapsed_strong_snp_cluster_from_pileup(
            "chr1",
            100,
            b"AAT",
            &active_loci,
            &valid_events,
            20.0,
        );

        assert!(rescued.is_empty());
    }

    #[test]
    fn rescue_collapsed_strong_snp_cluster_from_pileup_recovers_single_exact_match_for_isolated_snp_locus(
    ) {
        let valid_events = vec![
            test_variant("chr1", 100, b"A", b"T"),
            test_variant("chr1", 99, b"G", b"GA"),
        ];
        let active_loci = vec![ReplayActiveLocusRow {
            contig: "chr1".to_string(),
            pos: 100,
            region: "chr1:99-101".to_string(),
            ref_base: b'A',
            depth: 16,
            snp_alt_count: 8,
            snp_best_alt: "T".to_string(),
            indel_alt_count: 0,
            indel_best_alt: String::new(),
            alt_fraction: 0.5,
            active_probability_proxy: 1.0,
        }];

        let rescued = rescue_collapsed_strong_snp_cluster_from_pileup(
            "chr1",
            99,
            b"GAA",
            &active_loci,
            &valid_events,
            20.0,
        );

        assert_eq!(rescued.len(), 1);
        assert_eq!(rescued[0].pos, 100);
        assert_eq!(rescued[0].ref_allele, b"A");
        assert_eq!(rescued[0].alt_allele, b"T");
    }

    #[test]
    fn rescue_collapsed_strong_snp_cluster_from_pileup_recovers_single_exact_match_with_only_weak_other_active_loci(
    ) {
        let valid_events = vec![
            test_variant("chr1", 100, b"A", b"T"),
            test_variant("chr1", 99, b"G", b"GA"),
        ];
        let active_loci = vec![
            ReplayActiveLocusRow {
                contig: "chr1".to_string(),
                pos: 100,
                region: "chr1:99-106".to_string(),
                ref_base: b'A',
                depth: 16,
                snp_alt_count: 8,
                snp_best_alt: "T".to_string(),
                indel_alt_count: 0,
                indel_best_alt: String::new(),
                alt_fraction: 0.5,
                active_probability_proxy: 1.0,
            },
            ReplayActiveLocusRow {
                contig: "chr1".to_string(),
                pos: 103,
                region: "chr1:99-106".to_string(),
                ref_base: b'C',
                depth: 6,
                snp_alt_count: 1,
                snp_best_alt: "A".to_string(),
                indel_alt_count: 0,
                indel_best_alt: String::new(),
                alt_fraction: 0.166667,
                active_probability_proxy: 1.0,
            },
            ReplayActiveLocusRow {
                contig: "chr1".to_string(),
                pos: 105,
                region: "chr1:99-106".to_string(),
                ref_base: b'G',
                depth: 6,
                snp_alt_count: 1,
                snp_best_alt: "C".to_string(),
                indel_alt_count: 0,
                indel_best_alt: String::new(),
                alt_fraction: 0.166667,
                active_probability_proxy: 1.0,
            },
        ];

        let rescued = rescue_collapsed_strong_snp_cluster_from_pileup(
            "chr1",
            99,
            b"GAAACCG",
            &active_loci,
            &valid_events,
            20.0,
        );

        assert_eq!(rescued.len(), 1);
        assert_eq!(rescued[0].pos, 100);
        assert_eq!(rescued[0].ref_allele, b"A");
        assert_eq!(rescued[0].alt_allele, b"T");
    }

    #[test]
    fn rescue_collapsed_strong_snp_cluster_from_pileup_does_not_recover_single_match_with_indel_evidence(
    ) {
        let valid_events = vec![
            test_variant("chr1", 100, b"A", b"T"),
            test_variant("chr1", 99, b"G", b"GA"),
        ];
        let active_loci = vec![ReplayActiveLocusRow {
            contig: "chr1".to_string(),
            pos: 100,
            region: "chr1:99-101".to_string(),
            ref_base: b'A',
            depth: 16,
            snp_alt_count: 8,
            snp_best_alt: "T".to_string(),
            indel_alt_count: 4,
            indel_best_alt: "INS:A".to_string(),
            alt_fraction: 0.5,
            active_probability_proxy: 1.0,
        }];

        let rescued = rescue_collapsed_strong_snp_cluster_from_pileup(
            "chr1",
            99,
            b"GAA",
            &active_loci,
            &valid_events,
            20.0,
        );

        assert!(rescued.is_empty());
    }

    #[test]
    fn merge_missing_strong_snp_cluster_rescues_from_pileup_adds_missing_exact_cluster_call() {
        let valid_events = vec![
            test_variant("chr1", 100, b"A", b"G"),
            test_variant("chr1", 101, b"A", b"C"),
            test_variant("chr1", 102, b"A", b"T"),
        ];
        let active_loci = vec![
            ReplayActiveLocusRow {
                contig: "chr1".to_string(),
                pos: 100,
                region: "chr1:100-103".to_string(),
                ref_base: b'A',
                depth: 16,
                snp_alt_count: 8,
                snp_best_alt: "G".to_string(),
                indel_alt_count: 0,
                indel_best_alt: String::new(),
                alt_fraction: 0.5,
                active_probability_proxy: 1.0,
            },
            ReplayActiveLocusRow {
                contig: "chr1".to_string(),
                pos: 101,
                region: "chr1:100-103".to_string(),
                ref_base: b'A',
                depth: 16,
                snp_alt_count: 7,
                snp_best_alt: "C".to_string(),
                indel_alt_count: 0,
                indel_best_alt: String::new(),
                alt_fraction: 0.4375,
                active_probability_proxy: 1.0,
            },
            ReplayActiveLocusRow {
                contig: "chr1".to_string(),
                pos: 102,
                region: "chr1:100-103".to_string(),
                ref_base: b'A',
                depth: 16,
                snp_alt_count: 6,
                snp_best_alt: "T".to_string(),
                indel_alt_count: 0,
                indel_best_alt: String::new(),
                alt_fraction: 0.375,
                active_probability_proxy: 1.0,
            },
        ];
        let mut final_calls = vec![
            test_variant("chr1", 100, b"A", b"G"),
            test_variant("chr1", 102, b"A", b"T"),
        ];

        merge_missing_strong_snp_cluster_rescues_from_pileup(
            &mut final_calls,
            "chr1",
            100,
            b"AAAT",
            &active_loci,
            &valid_events,
            20.0,
        );

        assert_eq!(final_calls.len(), 3);
        assert!(final_calls
            .iter()
            .any(|event| event.pos == 101 && event.ref_allele == b"A" && event.alt_allele == b"C"));
    }

    #[test]
    fn merge_missing_strong_snp_cluster_rescues_from_pileup_does_not_add_single_isolated_match() {
        let valid_events = vec![
            test_variant("chr1", 100, b"A", b"T"),
            test_variant("chr1", 100, b"A", b"AG"),
        ];
        let active_loci = vec![ReplayActiveLocusRow {
            contig: "chr1".to_string(),
            pos: 100,
            region: "chr1:100-101".to_string(),
            ref_base: b'A',
            depth: 16,
            snp_alt_count: 4,
            snp_best_alt: "T".to_string(),
            indel_alt_count: 0,
            indel_best_alt: String::new(),
            alt_fraction: 0.25,
            active_probability_proxy: 1.0,
        }];
        let mut final_calls = vec![test_variant("chr1", 100, b"A", b"AG")];

        merge_missing_strong_snp_cluster_rescues_from_pileup(
            &mut final_calls,
            "chr1",
            100,
            b"AA",
            &active_loci,
            &valid_events,
            20.0,
        );

        assert_eq!(final_calls.len(), 1);
        assert_eq!(final_calls[0].alt_allele, b"AG");
    }

    #[test]
    fn merge_missing_strong_snp_cluster_rescues_from_pileup_adds_single_high_confidence_match() {
        let valid_events = vec![
            test_variant("chr1", 100, b"A", b"T"),
            test_variant("chr1", 100, b"A", b"AG"),
        ];
        let active_loci = vec![ReplayActiveLocusRow {
            contig: "chr1".to_string(),
            pos: 100,
            region: "chr1:100-101".to_string(),
            ref_base: b'A',
            depth: 16,
            snp_alt_count: 12,
            snp_best_alt: "T".to_string(),
            indel_alt_count: 0,
            indel_best_alt: String::new(),
            alt_fraction: 0.75,
            active_probability_proxy: 1.0,
        }];
        let mut final_calls = vec![test_variant("chr1", 100, b"A", b"AG")];

        merge_missing_strong_snp_cluster_rescues_from_pileup(
            &mut final_calls,
            "chr1",
            100,
            b"AA",
            &active_loci,
            &valid_events,
            20.0,
        );

        assert_eq!(final_calls.len(), 2);
        assert!(final_calls
            .iter()
            .any(|event| event.pos == 100 && event.ref_allele == b"A" && event.alt_allele == b"T"));
    }

    #[test]
    fn prune_unsupported_simple_snp_calls_in_dense_clusters_drops_dense_non_active_snps() {
        let active_loci = vec![
            ReplayActiveLocusRow {
                contig: "chr1".to_string(),
                pos: 110,
                region: "chr1:100-140".to_string(),
                ref_base: b'A',
                depth: 20,
                snp_alt_count: 10,
                snp_best_alt: "C".to_string(),
                indel_alt_count: 0,
                indel_best_alt: String::new(),
                alt_fraction: 0.5,
                active_probability_proxy: 1.0,
            },
            ReplayActiveLocusRow {
                contig: "chr1".to_string(),
                pos: 136,
                region: "chr1:100-140".to_string(),
                ref_base: b'G',
                depth: 18,
                snp_alt_count: 9,
                snp_best_alt: "T".to_string(),
                indel_alt_count: 0,
                indel_best_alt: String::new(),
                alt_fraction: 0.5,
                active_probability_proxy: 1.0,
            },
        ];
        let mut final_calls = vec![
            test_variant("chr1", 100, b"A", b"G"),
            test_variant("chr1", 110, b"A", b"C"),
            test_variant("chr1", 118, b"T", b"TA"),
            test_variant("chr1", 120, b"T", b"G"),
            test_variant("chr1", 136, b"G", b"T"),
        ];

        prune_unsupported_simple_snp_calls_in_dense_clusters(&mut final_calls, &active_loci);

        assert_eq!(final_calls.len(), 3);
        assert!(final_calls
            .iter()
            .any(|call| call.pos == 110 && call.alt_allele == b"C"));
        assert!(final_calls
            .iter()
            .any(|call| call.pos == 118 && call.alt_allele == b"TA"));
        assert!(final_calls
            .iter()
            .any(|call| call.pos == 136 && call.alt_allele == b"T"));
        assert!(!final_calls
            .iter()
            .any(|call| call.pos == 100 && call.alt_allele == b"G"));
        assert!(!final_calls
            .iter()
            .any(|call| call.pos == 120 && call.alt_allele == b"G"));
    }

    #[test]
    fn prune_unsupported_simple_snp_calls_in_dense_clusters_keeps_isolated_unsupported_snp() {
        let active_loci = vec![ReplayActiveLocusRow {
            contig: "chr1".to_string(),
            pos: 110,
            region: "chr1:100-160".to_string(),
            ref_base: b'A',
            depth: 20,
            snp_alt_count: 10,
            snp_best_alt: "C".to_string(),
            indel_alt_count: 0,
            indel_best_alt: String::new(),
            alt_fraction: 0.5,
            active_probability_proxy: 1.0,
        }];
        let mut final_calls = vec![
            test_variant("chr1", 110, b"A", b"C"),
            test_variant("chr1", 140, b"T", b"G"),
        ];

        prune_unsupported_simple_snp_calls_in_dense_clusters(&mut final_calls, &active_loci);

        assert_eq!(final_calls.len(), 2);
        assert!(final_calls
            .iter()
            .any(|call| call.pos == 110 && call.alt_allele == b"C"));
        assert!(final_calls
            .iter()
            .any(|call| call.pos == 140 && call.alt_allele == b"G"));
    }

    #[test]
    fn supplement_missing_pileup_events_adds_missing_event_and_haplotype() {
        let mut valid_events = vec![test_variant("chr1", 100, b"A", b"G")];
        let mut local_haplotypes =
            haplotypes_from_candidate_events("chr1", 100, b"AAT", &valid_events);
        let active_loci = vec![
            ReplayActiveLocusRow {
                contig: "chr1".to_string(),
                pos: 100,
                region: "chr1:100-102".to_string(),
                ref_base: b'A',
                depth: 16,
                snp_alt_count: 8,
                snp_best_alt: "G".to_string(),
                indel_alt_count: 0,
                indel_best_alt: String::new(),
                alt_fraction: 0.5,
                active_probability_proxy: 1.0,
            },
            ReplayActiveLocusRow {
                contig: "chr1".to_string(),
                pos: 101,
                region: "chr1:100-102".to_string(),
                ref_base: b'A',
                depth: 16,
                snp_alt_count: 8,
                snp_best_alt: "C".to_string(),
                indel_alt_count: 0,
                indel_best_alt: String::new(),
                alt_fraction: 0.5,
                active_probability_proxy: 1.0,
            },
        ];

        supplement_missing_pileup_events(
            "chr1",
            100,
            b"AAT",
            &active_loci,
            20.0,
            &mut local_haplotypes,
            &mut valid_events,
        );

        assert_eq!(valid_events.len(), 2);
        assert_eq!(valid_events[1].pos, 101);
        assert_eq!(valid_events[1].ref_allele, b"A");
        assert_eq!(valid_events[1].alt_allele, b"C");
        assert!(local_haplotypes
            .iter()
            .any(|hap| hap.event_indices == vec![1] && hap.bases == b"ACT"));
    }

    #[test]
    fn supplement_missing_pileup_events_seeds_low_qual_zero_candidate_simple_snp() {
        let mut valid_events = Vec::new();
        let mut local_haplotypes = haplotypes_from_candidate_events("chr1", 100, b"A", &[]);
        let active_loci = vec![ReplayActiveLocusRow {
            contig: "chr1".to_string(),
            pos: 100,
            region: "chr1:100-100".to_string(),
            ref_base: b'A',
            depth: 24,
            snp_alt_count: 3,
            snp_best_alt: "G".to_string(),
            indel_alt_count: 0,
            indel_best_alt: String::new(),
            alt_fraction: 0.125,
            active_probability_proxy: 1.0,
        }];

        supplement_missing_pileup_events(
            "chr1",
            100,
            b"A",
            &active_loci,
            20.0,
            &mut local_haplotypes,
            &mut valid_events,
        );

        assert_eq!(valid_events.len(), 1);
        assert_eq!(valid_events[0].pos, 100);
        assert_eq!(valid_events[0].ref_allele, b"A");
        assert_eq!(valid_events[0].alt_allele, b"G");
        assert!(valid_events[0].qual < 20);
        assert!(local_haplotypes
            .iter()
            .any(|hap| !hap.is_ref && hap.event_indices == vec![0] && hap.bases == b"G"));
    }

    #[test]
    fn supplement_missing_pileup_events_skips_low_qual_pileup_het_seed() {
        let mut valid_events = Vec::new();
        let mut local_haplotypes = haplotypes_from_candidate_events("chr1", 100, b"G", &[]);
        let active_loci = vec![ReplayActiveLocusRow {
            contig: "chr1".to_string(),
            pos: 100,
            region: "chr1:100-100".to_string(),
            ref_base: b'G',
            depth: 17,
            snp_alt_count: 3,
            snp_best_alt: "C".to_string(),
            indel_alt_count: 0,
            indel_best_alt: String::new(),
            alt_fraction: 0.176471,
            active_probability_proxy: 1.0,
        }];

        supplement_missing_pileup_events(
            "chr1",
            100,
            b"G",
            &active_loci,
            20.0,
            &mut local_haplotypes,
            &mut valid_events,
        );

        assert!(valid_events.is_empty());
        assert!(local_haplotypes.iter().all(|hap| hap.is_ref));
    }

    #[test]
    fn supplement_missing_pileup_events_overlays_missing_snp_on_existing_alt_haplotype() {
        let mut valid_events = vec![test_variant("chr1", 100, b"A", b"G")];
        let mut local_haplotypes =
            haplotypes_from_candidate_events("chr1", 100, b"AAT", &valid_events);
        let active_loci = vec![
            ReplayActiveLocusRow {
                contig: "chr1".to_string(),
                pos: 100,
                region: "chr1:100-102".to_string(),
                ref_base: b'A',
                depth: 16,
                snp_alt_count: 8,
                snp_best_alt: "G".to_string(),
                indel_alt_count: 0,
                indel_best_alt: String::new(),
                alt_fraction: 0.5,
                active_probability_proxy: 1.0,
            },
            ReplayActiveLocusRow {
                contig: "chr1".to_string(),
                pos: 101,
                region: "chr1:100-102".to_string(),
                ref_base: b'A',
                depth: 16,
                snp_alt_count: 8,
                snp_best_alt: "C".to_string(),
                indel_alt_count: 0,
                indel_best_alt: String::new(),
                alt_fraction: 0.5,
                active_probability_proxy: 1.0,
            },
        ];

        supplement_missing_pileup_events(
            "chr1",
            100,
            b"AAT",
            &active_loci,
            20.0,
            &mut local_haplotypes,
            &mut valid_events,
        );

        assert_eq!(valid_events.len(), 2);
        assert!(local_haplotypes
            .iter()
            .any(|hap| hap.event_indices == vec![0, 1] && hap.bases == b"GCT"));
    }

    #[test]
    fn supplement_missing_pileup_events_does_not_overlay_missing_snp_on_indel_haplotype() {
        let mut valid_events = vec![test_variant("chr1", 100, b"T", b"TA")];
        let mut local_haplotypes =
            haplotypes_from_candidate_events("chr1", 100, b"TT", &valid_events);
        let active_loci = vec![
            ReplayActiveLocusRow {
                contig: "chr1".to_string(),
                pos: 100,
                region: "chr1:100-101".to_string(),
                ref_base: b'T',
                depth: 16,
                snp_alt_count: 0,
                snp_best_alt: String::new(),
                indel_alt_count: 8,
                indel_best_alt: "INS:A".to_string(),
                alt_fraction: 0.5,
                active_probability_proxy: 1.0,
            },
            ReplayActiveLocusRow {
                contig: "chr1".to_string(),
                pos: 100,
                region: "chr1:100-101".to_string(),
                ref_base: b'T',
                depth: 16,
                snp_alt_count: 8,
                snp_best_alt: "C".to_string(),
                indel_alt_count: 8,
                indel_best_alt: "INS:A".to_string(),
                alt_fraction: 0.5,
                active_probability_proxy: 1.0,
            },
        ];

        supplement_missing_pileup_events(
            "chr1",
            100,
            b"TT",
            &active_loci,
            20.0,
            &mut local_haplotypes,
            &mut valid_events,
        );

        assert_eq!(valid_events.len(), 2);
        assert!(local_haplotypes
            .iter()
            .any(|hap| hap.event_indices == vec![0]));
        assert!(local_haplotypes
            .iter()
            .any(|hap| hap.event_indices == vec![1]));
        assert!(!local_haplotypes
            .iter()
            .any(|hap| hap.event_indices == vec![0, 1]));
    }

    #[test]
    fn supplement_missing_pileup_events_skips_indels_when_region_already_has_events() {
        let mut valid_events = vec![test_variant("chr1", 100, b"A", b"G")];
        let mut local_haplotypes =
            haplotypes_from_candidate_events("chr1", 100, b"AAAT", &valid_events);
        let active_loci = vec![ReplayActiveLocusRow {
            contig: "chr1".to_string(),
            pos: 101,
            region: "chr1:100-103".to_string(),
            ref_base: b'A',
            depth: 16,
            snp_alt_count: 0,
            snp_best_alt: String::new(),
            indel_alt_count: 8,
            indel_best_alt: "INS:A".to_string(),
            alt_fraction: 0.5,
            active_probability_proxy: 1.0,
        }];

        supplement_missing_pileup_events(
            "chr1",
            100,
            b"AAAT",
            &active_loci,
            20.0,
            &mut local_haplotypes,
            &mut valid_events,
        );

        assert_eq!(valid_events.len(), 1);
        assert_eq!(valid_events[0].ref_allele, b"A");
        assert_eq!(valid_events[0].alt_allele, b"G");
    }

    #[test]
    fn supplement_missing_pileup_events_skips_weak_snp_in_dense_existing_snp_cluster() {
        let local_ref_bases = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        let mut valid_events = vec![
            test_variant("chr1", 100, b"A", b"G"),
            test_variant("chr1", 121, b"A", b"C"),
        ];
        let mut local_haplotypes =
            haplotypes_from_candidate_events("chr1", 100, local_ref_bases, &valid_events);
        let active_loci = vec![
            ReplayActiveLocusRow {
                contig: "chr1".to_string(),
                pos: 100,
                region: "chr1:100-129".to_string(),
                ref_base: b'A',
                depth: 16,
                snp_alt_count: 8,
                snp_best_alt: "G".to_string(),
                indel_alt_count: 0,
                indel_best_alt: String::new(),
                alt_fraction: 0.5,
                active_probability_proxy: 1.0,
            },
            ReplayActiveLocusRow {
                contig: "chr1".to_string(),
                pos: 121,
                region: "chr1:100-129".to_string(),
                ref_base: b'A',
                depth: 16,
                snp_alt_count: 8,
                snp_best_alt: "C".to_string(),
                indel_alt_count: 0,
                indel_best_alt: String::new(),
                alt_fraction: 0.5,
                active_probability_proxy: 1.0,
            },
            ReplayActiveLocusRow {
                contig: "chr1".to_string(),
                pos: 127,
                region: "chr1:100-129".to_string(),
                ref_base: b'A',
                depth: 14,
                snp_alt_count: 4,
                snp_best_alt: "T".to_string(),
                indel_alt_count: 0,
                indel_best_alt: String::new(),
                alt_fraction: 0.285714,
                active_probability_proxy: 1.0,
            },
        ];

        supplement_missing_pileup_events(
            "chr1",
            100,
            local_ref_bases,
            &active_loci,
            20.0,
            &mut local_haplotypes,
            &mut valid_events,
        );

        assert_eq!(valid_events.len(), 2);
        assert!(!valid_events.iter().any(|event| event.pos == 127));
    }

    #[test]
    fn supplement_missing_pileup_events_keeps_stronger_snp_in_dense_existing_snp_cluster() {
        let local_ref_bases = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        let mut valid_events = vec![
            test_variant("chr1", 100, b"A", b"G"),
            test_variant("chr1", 121, b"A", b"C"),
        ];
        let mut local_haplotypes =
            haplotypes_from_candidate_events("chr1", 100, local_ref_bases, &valid_events);
        let active_loci = vec![
            ReplayActiveLocusRow {
                contig: "chr1".to_string(),
                pos: 100,
                region: "chr1:100-129".to_string(),
                ref_base: b'A',
                depth: 16,
                snp_alt_count: 8,
                snp_best_alt: "G".to_string(),
                indel_alt_count: 0,
                indel_best_alt: String::new(),
                alt_fraction: 0.5,
                active_probability_proxy: 1.0,
            },
            ReplayActiveLocusRow {
                contig: "chr1".to_string(),
                pos: 121,
                region: "chr1:100-129".to_string(),
                ref_base: b'A',
                depth: 16,
                snp_alt_count: 8,
                snp_best_alt: "C".to_string(),
                indel_alt_count: 0,
                indel_best_alt: String::new(),
                alt_fraction: 0.5,
                active_probability_proxy: 1.0,
            },
            ReplayActiveLocusRow {
                contig: "chr1".to_string(),
                pos: 127,
                region: "chr1:100-129".to_string(),
                ref_base: b'A',
                depth: 16,
                snp_alt_count: 5,
                snp_best_alt: "T".to_string(),
                indel_alt_count: 0,
                indel_best_alt: String::new(),
                alt_fraction: 0.3125,
                active_probability_proxy: 1.0,
            },
        ];

        supplement_missing_pileup_events(
            "chr1",
            100,
            local_ref_bases,
            &active_loci,
            20.0,
            &mut local_haplotypes,
            &mut valid_events,
        );

        assert_eq!(valid_events.len(), 3);
        assert!(valid_events.iter().any(|event| event.pos == 127));
    }

    #[test]
    fn assemble_haplotypes_skips_non_acgt_haplotype_sequences() {
        let (haplotypes, events) = assemble_haplotypes(
            "chr1",
            100,
            b"AAACC",
            &[
                b"AAGCC".to_vec(),
                b"AAGCC".to_vec(),
                b"AANCC".to_vec(),
                b"AANCC".to_vec(),
            ],
            &[2],
            0,
        );

        assert!(events
            .iter()
            .any(|event| event.pos == 102 && event.ref_allele == b"A" && event.alt_allele == b"G"));
        assert!(haplotypes
            .iter()
            .any(|hap| !hap.is_ref && hap.bases == b"AAGCC"));
        assert!(haplotypes.iter().all(|hap| is_regular_bases(&hap.bases)));
    }

    #[test]
    fn align_haplotype_to_reference_rejects_softclipped_path_alignments() {
        let align_result = align_haplotype_to_reference(b"ACGTACGT", b"ACGTACGTGGGG");
        assert!(align_result.is_none());
    }

    #[test]
    fn align_haplotype_to_reference_keeps_simple_indel_alignments() {
        let align_result = align_haplotype_to_reference(b"ACGTACGT", b"ACGTTACGT").unwrap();
        assert_eq!(align_result.alignment_offset, 0);
        assert!(align_result
            .cigar
            .iter()
            .all(|ce| !matches!(ce, rust_htslib::bam::record::Cigar::SoftClip(_))));
        assert_eq!(
            align_result
                .cigar
                .iter()
                .filter(|ce| matches!(ce, rust_htslib::bam::record::Cigar::Ins(_)))
                .count(),
            1
        );
    }
