#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::Aux;

    fn cigar_text(cigar: &CigarString) -> String {
        cigar.to_string()
    }

    fn record_with_cigar(cigar: Vec<Cigar>, flags: u16) -> bam::Record {
        let cigar = CigarString(cigar);
        let read_len = cigar.iter().copied().map(read_len).sum::<u32>() as usize;
        let bases = vec![b'A'; read_len];
        let quals = vec![30u8; read_len];
        let mut record = bam::Record::new();
        record.set(b"read1", Some(&cigar), &bases, &quals);
        record.set_tid(0);
        record.set_pos(100);
        record.set_mtid(0);
        record.set_mpos(200);
        record.set_insert_size(150);
        record.set_mapq(255);
        record.set_flags(flags);
        record
    }

    #[test]
    fn splits_simple_n_cigar() {
        let plans =
            split_cigar_at_ref_skips(&[Cigar::Match(5), Cigar::RefSkip(10), Cigar::Match(5)])
                .unwrap();

        assert_eq!(plans.len(), 2);
        assert_eq!(plans[0].ref_offset, 0);
        assert_eq!(cigar_text(&plans[0].cigar), "5M5S");
        assert_eq!(plans[1].ref_offset, 15);
        assert_eq!(cigar_text(&plans[1].cigar), "5S5M");
    }

    #[test]
    fn preserves_soft_and_hard_clips() {
        let plans = split_cigar_at_ref_skips(&[
            Cigar::HardClip(1),
            Cigar::SoftClip(2),
            Cigar::Match(3),
            Cigar::RefSkip(8),
            Cigar::Match(4),
            Cigar::SoftClip(5),
            Cigar::HardClip(6),
        ])
        .unwrap();

        assert_eq!(plans.len(), 2);
        assert_eq!(cigar_text(&plans[0].cigar), "1H2S3M9S6H");
        assert_eq!(plans[1].ref_offset, 11);
        assert_eq!(cigar_text(&plans[1].cigar), "1H5S4M5S6H");
    }

    #[test]
    fn handles_insertions_and_deletions() {
        let plans = split_cigar_at_ref_skips(&[
            Cigar::HardClip(1),
            Cigar::Match(2),
            Cigar::Del(2),
            Cigar::Match(1),
            Cigar::RefSkip(2),
            Cigar::Match(1),
            Cigar::Ins(2),
            Cigar::RefSkip(1),
            Cigar::Match(1),
            Cigar::SoftClip(2),
        ])
        .unwrap();

        assert_eq!(plans.len(), 3);
        assert_eq!(cigar_text(&plans[0].cigar), "1H2M2D1M6S");
        assert_eq!(plans[1].ref_offset, 7);
        assert_eq!(cigar_text(&plans[1].cigar), "1H3S1M2I3S");
        assert_eq!(plans[2].ref_offset, 9);
        assert_eq!(cigar_text(&plans[2].cigar), "1H6S1M2S");
    }

    #[test]
    fn trims_deletions_at_split_edges() {
        let plans = split_cigar_at_ref_skips(&[
            Cigar::Match(4),
            Cigar::Del(3),
            Cigar::RefSkip(5),
            Cigar::Del(2),
            Cigar::Match(4),
        ])
        .unwrap();

        assert_eq!(plans.len(), 2);
        assert_eq!(cigar_text(&plans[0].cigar), "4M4S");
        assert_eq!(plans[1].ref_offset, 14);
        assert_eq!(cigar_text(&plans[1].cigar), "4S4M");
    }

    #[test]
    fn leaves_bogus_leading_n_only_split_unchanged() {
        let plans = split_cigar_at_ref_skips(&[
            Cigar::SoftClip(1),
            Cigar::RefSkip(3),
            Cigar::Match(2),
            Cigar::HardClip(4),
        ])
        .unwrap();

        assert!(plans.is_empty());
    }

    #[test]
    fn transform_splits_paired_read_and_removes_stale_tags() {
        let mut record = record_with_cigar(
            vec![Cigar::Match(5), Cigar::RefSkip(10), Cigar::Match(5)],
            0x1 | 0x2 | 0x40,
        );
        record.push_aux(b"NM", Aux::I32(1)).unwrap();
        record.push_aux(b"MD", Aux::String("5")).unwrap();
        record.push_aux(b"NH", Aux::I32(1)).unwrap();
        record.push_aux(b"MC", Aux::String("10M")).unwrap();

        let records = transform_record(&record, SplitOptions::default()).unwrap();

        assert_eq!(records.len(), 2);
        assert_eq!(records[0].mapq(), 60);
        assert_eq!(records[0].pos(), 100);
        assert_eq!(records[1].pos(), 115);
        assert_eq!(records[0].flags() & 0x1, 0x1);
        assert_eq!(records[1].flags() & SUPPLEMENTARY_FLAG, SUPPLEMENTARY_FLAG);
        assert_eq!(records[1].mtid(), 0);
        assert_eq!(records[1].mpos(), 200);
        assert!(records[0].aux(b"NM").is_err());
        assert!(records[0].aux(b"MD").is_err());
        assert!(records[0].aux(b"NH").is_err());
        assert!(records[0].aux(b"MC").is_err());
    }

    #[test]
    fn compatibility_mode_repairs_sa_tags_for_split_family() {
        let mut record = record_with_cigar(
            vec![Cigar::Match(5), Cigar::RefSkip(10), Cigar::Match(5)],
            0x1 | 0x2 | 0x40,
        );
        record
            .push_aux(b"SA", Aux::String("chr9,7,-,3M,40,2;"))
            .unwrap();
        let records = transform_record(
            &record,
            SplitOptions {
                mode: SplitMode::Compatibility,
                ..SplitOptions::default()
            },
        )
        .unwrap();

        assert_eq!(records.len(), 2);
        assert_eq!(records[0].flags() & SUPPLEMENTARY_FLAG, 0);
        assert_eq!(records[1].flags() & SUPPLEMENTARY_FLAG, SUPPLEMENTARY_FLAG);
        assert_eq!(
            records[0].aux(b"SA").unwrap(),
            Aux::String("chr9,7,-,3M,40,2;1,116,+,5S5M,60,*;")
        );
        assert_eq!(
            records[1].aux(b"SA").unwrap(),
            Aux::String("1,101,+,5M5S,60,*;chr9,7,-,3M,40,2;")
        );
    }

    #[test]
    fn compatibility_mode_uses_header_contig_names_in_sa_tags() {
        let record = record_with_cigar(
            vec![Cigar::Match(5), Cigar::RefSkip(10), Cigar::Match(5)],
            0,
        );
        let contigs = vec!["chr1".to_string()];
        let records = transform_record_with_contig_names(
            &record,
            SplitOptions {
                mode: SplitMode::Compatibility,
                ..SplitOptions::default()
            },
            Some(&contigs),
        )
        .unwrap();

        assert_eq!(
            records[0].aux(b"SA").unwrap(),
            Aux::String("chr1,116,+,5S5M,60,*;")
        );
    }

    #[test]
    fn compatibility_mode_repairs_mate_cigar_to_first_split_segment() {
        let mut record = record_with_cigar(
            vec![Cigar::Match(5), Cigar::RefSkip(10), Cigar::Match(5)],
            0x1 | 0x2 | 0x40,
        );
        record.push_aux(b"MC", Aux::String("3M7N4M")).unwrap();

        let records = transform_record(
            &record,
            SplitOptions {
                mode: SplitMode::Compatibility,
                ..SplitOptions::default()
            },
        )
        .unwrap();

        assert_eq!(records[0].aux(b"MC").unwrap(), Aux::String("3M4S"));
        assert_eq!(records[1].aux(b"MC").unwrap(), Aux::String("3M4S"));
    }

    #[test]
    fn compatibility_mode_keeps_unsplit_record_repairs_mc_and_removes_stale_tags() {
        let mut record = record_with_cigar(vec![Cigar::Match(5)], 0x1 | 0x2 | 0x40);
        record.push_aux(b"NM", Aux::I32(1)).unwrap();
        record.push_aux(b"MD", Aux::String("5")).unwrap();
        record.push_aux(b"NH", Aux::I32(1)).unwrap();
        record.push_aux(b"MC", Aux::String("3M7N4M")).unwrap();

        let records = transform_record(
            &record,
            SplitOptions {
                mode: SplitMode::Compatibility,
                ..SplitOptions::default()
            },
        )
        .unwrap();

        assert_eq!(records.len(), 1);
        assert_eq!(records[0].aux(b"MC").unwrap(), Aux::String("3M4S"));
        assert!(records[0].aux(b"NM").is_err());
        assert!(records[0].aux(b"MD").is_err());
        assert!(records[0].aux(b"NH").is_err());
    }

    #[test]
    fn compatibility_mode_skips_secondary_but_repairs_mc() {
        let mut record = record_with_cigar(
            vec![Cigar::Match(5), Cigar::RefSkip(10), Cigar::Match(5)],
            0x100,
        );
        record.push_aux(b"MC", Aux::String("3M7N4M")).unwrap();

        let records = transform_record(
            &record,
            SplitOptions {
                mode: SplitMode::Compatibility,
                ..SplitOptions::default()
            },
        )
        .unwrap();

        assert_eq!(records.len(), 1);
        assert_eq!(records[0].flags() & 0x100, 0x100);
        assert_eq!(records[0].cigar().to_string(), "5M10N5M");
        assert_eq!(records[0].aux(b"MC").unwrap(), Aux::String("3M4S"));
    }

    #[test]
    fn overhang_fixing_soft_clips_left_mismatching_overhang() {
        let mut record = record_with_cigar(vec![Cigar::SoftClip(2), Cigar::Match(8)], 0);
        record.set_pos(103);
        let splice = Splice {
            tid: 0,
            start: 100,
            end: 104,
            reference: b"TTTTT".to_vec(),
        };
        let mut read = ManagedRead::new(record);

        fix_split_with_options(&mut read, &splice, OverhangOptions::default()).unwrap();

        assert_eq!(read.record.pos(), 105);
        assert_eq!(read.record.cigar().to_string(), "4S6M");
    }

    #[test]
    fn overhang_fixing_soft_clips_right_mismatching_overhang() {
        let mut record = record_with_cigar(vec![Cigar::Match(8), Cigar::SoftClip(2)], 0);
        record.set_pos(95);
        let splice = Splice {
            tid: 0,
            start: 100,
            end: 106,
            reference: b"TTTTTTT".to_vec(),
        };
        let mut read = ManagedRead::new(record);

        fix_split_with_options(&mut read, &splice, OverhangOptions::default()).unwrap();

        assert_eq!(read.record.pos(), 95);
        assert_eq!(read.record.cigar().to_string(), "5M5S");
    }

    #[test]
    fn overhang_fixing_preserves_opposite_end_soft_clip() {
        let mut record = record_with_cigar(
            vec![Cigar::SoftClip(1), Cigar::Match(72), Cigar::SoftClip(27)],
            0,
        );
        record.set_pos(1_717_961 - 1);

        soft_clip_by_read_coordinates(&mut record, 73, 99).unwrap();

        assert_eq!(record.pos(), 1_717_961 - 1);
        assert_eq!(record.cigar().to_string(), "1S72M27S");
    }

    #[test]
    fn overhang_fixing_left_clip_across_deletion_preserves_read_length() {
        let mut record =
            record_with_cigar(vec![Cigar::Match(30), Cigar::Del(1), Cigar::Match(55)], 0);
        record.set_pos(3_576_525 - 1);

        soft_clip_by_read_coordinates(&mut record, 0, 33).unwrap();

        assert_eq!(record.pos(), 3_576_560 - 1);
        assert_eq!(record.cigar().to_string(), "34S51M");
        assert_eq!(
            record.cigar().iter().copied().map(read_len).sum::<u32>() as usize,
            record.seq_len()
        );
    }

    #[test]
    fn overhang_fixing_left_clip_preserves_trailing_soft_clip_and_insertions() {
        let mut record = record_with_cigar(
            vec![
                Cigar::SoftClip(4),
                Cigar::Match(56),
                Cigar::Ins(2),
                Cigar::Match(21),
                Cigar::SoftClip(17),
            ],
            0,
        );
        record.set_pos(3_577_698 - 1);

        soft_clip_by_read_coordinates(&mut record, 0, 18).unwrap();

        assert_eq!(record.pos(), 3_577_713 - 1);
        assert_eq!(record.cigar().to_string(), "19S41M2I21M17S");
        assert_eq!(
            record.cigar().iter().copied().map(read_len).sum::<u32>() as usize,
            record.seq_len()
        );
    }

    #[test]
    fn overhang_fixing_keeps_matching_overhang() {
        let mut record = record_with_cigar(vec![Cigar::SoftClip(2), Cigar::Match(8)], 0);
        record.set_pos(103);
        let splice = Splice {
            tid: 0,
            start: 100,
            end: 104,
            reference: b"AAAAA".to_vec(),
        };
        let mut read = ManagedRead::new(record);

        fix_split_with_options(&mut read, &splice, OverhangOptions::default()).unwrap();

        assert_eq!(read.record.pos(), 103);
        assert_eq!(read.record.cigar().to_string(), "2S8M");
    }

    #[test]
    fn transform_skips_secondary_by_default() {
        let record = record_with_cigar(
            vec![Cigar::Match(5), Cigar::RefSkip(10), Cigar::Match(5)],
            0x100,
        );

        let records = transform_record(&record, SplitOptions::default()).unwrap();

        assert_eq!(records.len(), 1);
        assert_eq!(records[0].mapq(), 60);
        assert_eq!(records[0].cigar().to_string(), "5M10N5M");
    }

    #[test]
    fn transform_can_process_secondary_alignments() {
        let record = record_with_cigar(
            vec![Cigar::Match(5), Cigar::RefSkip(10), Cigar::Match(5)],
            0x100,
        );
        let records = transform_record(
            &record,
            SplitOptions {
                process_secondary_alignments: true,
                ..SplitOptions::default()
            },
        )
        .unwrap();

        assert_eq!(records.len(), 2);
        assert_eq!(records[0].flags() & 0x100, 0x100);
        assert_eq!(records[1].flags() & 0x100, 0x100);
        assert_eq!(records[1].flags() & SUPPLEMENTARY_FLAG, SUPPLEMENTARY_FLAG);
    }

    #[test]
    fn transform_preserves_supplementary_and_unmapped_mate_flags() {
        let record = record_with_cigar(
            vec![Cigar::Match(5), Cigar::RefSkip(10), Cigar::Match(5)],
            SUPPLEMENTARY_FLAG | 0x1 | 0x8,
        );

        let records = transform_record(&record, SplitOptions::default()).unwrap();

        assert_eq!(records.len(), 2);
        assert_eq!(records[0].flags() & SUPPLEMENTARY_FLAG, SUPPLEMENTARY_FLAG);
        assert_eq!(records[1].flags() & SUPPLEMENTARY_FLAG, SUPPLEMENTARY_FLAG);
        assert_eq!(records[0].flags() & 0x8, 0x8);
        assert_eq!(records[1].flags() & 0x8, 0x8);
    }
}
