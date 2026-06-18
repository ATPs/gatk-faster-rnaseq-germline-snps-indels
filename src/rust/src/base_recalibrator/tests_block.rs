#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::BTreeSet;

    fn test_read(reverse: bool) -> PreparedRead {
        PreparedRead {
            contig: "chr1".to_string(),
            start0: 0,
            end0: 3,
            bases: b"ACGT".to_vec(),
            quals: vec![30, 30, 30, 30],
            cigar: vec![SimpleCigar::Match(4)],
            rg_index: 0,
            is_reverse: reverse,
            is_second_in_pair: false,
        }
    }

    #[test]
    fn cycle_matches_gatk_orientation_rules() {
        assert_eq!(cycle_value(0, 5, false, false, 500).unwrap(), 1);
        assert_eq!(cycle_value(4, 5, false, false, 500).unwrap(), 5);
        assert_eq!(cycle_value(0, 5, true, false, 500).unwrap(), 5);
        assert_eq!(cycle_value(4, 5, true, false, 500).unwrap(), 1);
        assert_eq!(cycle_value(0, 5, false, true, 500).unwrap(), -1);
        assert_eq!(cycle_value(4, 5, true, true, 500).unwrap(), -1);
    }

    #[test]
    fn context_uses_preceding_base_and_reverse_complements_negative_reads() {
        assert_eq!(
            context_values(&test_read(false), 2, 2),
            vec![
                None,
                Some("AC".to_string()),
                Some("CG".to_string()),
                Some("GT".to_string())
            ]
        );
        assert_eq!(
            context_values(&test_read(true), 2, 2),
            vec![
                Some("GT".to_string()),
                Some("CG".to_string()),
                Some("AC".to_string()),
                None
            ]
        );
    }

    #[test]
    fn known_site_pointer_masks_overlapping_positions() {
        let intervals = vec![
            KnownInterval {
                start0: 10,
                end0: 12,
            },
            KnownInterval {
                start0: 20,
                end0: 20,
            },
        ];
        let mut index = 0;
        assert!(!is_known_site(9, &intervals, &mut index));
        assert!(is_known_site(10, &intervals, &mut index));
        assert!(is_known_site(12, &intervals, &mut index));
        assert!(!is_known_site(13, &intervals, &mut index));
        assert!(is_known_site(20, &intervals, &mut index));
    }

    #[test]
    fn empirical_quality_moves_down_for_high_error_bins() {
        let mut clean = Datum::new(30);
        for _ in 0..100 {
            clean.increment(0.0);
        }
        let mut noisy = Datum::new(30);
        for _ in 0..100 {
            noisy.increment(1.0);
        }
        assert!(empirical_quality(&clean) > empirical_quality(&noisy));
    }

    #[test]
    fn quantizer_keeps_requested_number_of_levels() {
        let mut counts = vec![0_u64; MAX_SAM_QUAL_SCORE + 1];
        counts[10] = 100;
        counts[20] = 100;
        counts[30] = 100;
        let map = quantize_quality_scores(&counts, 4, usize::from(MIN_USABLE_Q_SCORE));
        let levels: BTreeSet<u8> = map.into_iter().collect();
        assert!(levels.len() <= 4);
    }
}
