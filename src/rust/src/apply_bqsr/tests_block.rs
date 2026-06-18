#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::{Cigar, CigarString};

    #[test]
    fn synthetic_quality_recalibration_preserves_low_qualities() {
        let model = RecalibrationModel::from_str(
            &synthetic_table(
                &["rg1  M  30.0000  30.0000  1000000  1000.00"],
                &["rg1  30  M  35.0000  1000000  316.00"],
                &[],
            ),
            false,
        )
        .unwrap();

        let recalibrated = model
            .recalibrate_qualities("rg1", b"AAC", &[30, 30, 5], 0, false)
            .unwrap();
        assert_eq!(recalibrated, vec![34, 34, 5]);
    }

    #[test]
    fn synthetic_cycle_and_context_covariates_are_applied() {
        let model = RecalibrationModel::from_str(
            &synthetic_table(
                &["rg1  M  30.0000  30.0000  1000000  1000.00"],
                &["rg1  30  M  30.0000  1000000  1000.00"],
                &[
                    "rg1  30  1   Cycle    M  40.0000  1000000  100.00",
                    "rg1  30  AA  Context  M  40.0000  1000000  100.00",
                ],
            ),
            false,
        )
        .unwrap();

        let recalibrated = model
            .recalibrate_qualities("rg1", b"AA", &[30, 30], 0, false)
            .unwrap();
        assert_eq!(recalibrated, vec![37, 37]);
    }

    #[test]
    fn header_read_group_identifier_prefers_platform_unit() {
        let header = b"@HD\tVN:1.6\n@RG\tID:rg1\tSM:s1\tPU:unit1\n@RG\tID:rg2\tSM:s1\n";
        let read_groups = read_group_identifiers_from_header(header).unwrap();
        assert_eq!(read_groups.get("rg1").unwrap(), "unit1");
        assert_eq!(read_groups.get("rg2").unwrap(), "rg2");
    }

    #[test]
    fn transform_record_rewrites_qualities_and_preserves_alignment_fields() {
        let model = RecalibrationModel::from_str(
            &synthetic_table(
                &["rg1  M  30.0000  30.0000  1000000  1000.00"],
                &["rg1  30  M  35.0000  1000000  316.00"],
                &[],
            ),
            false,
        )
        .unwrap();
        let read_groups = HashMap::from([("rg1".to_owned(), "rg1".to_owned())]);
        let cigar = CigarString(vec![Cigar::Match(3)]);
        let mut record = Record::new();
        record.set(b"read1", Some(&cigar), b"AAC", &[10, 10, 10]);
        record.set_tid(2);
        record.set_pos(123);
        record.set_mtid(2);
        record.set_mpos(456);
        record.set_insert_size(333);
        record.set_flags(0x41);
        record.set_mapq(60);
        record.push_aux(b"RG", Aux::String("rg1")).unwrap();
        record.push_aux(b"OQ", Aux::String("???")).unwrap();
        record.push_aux(b"NM", Aux::I32(1)).unwrap();
        record.push_aux(b"BI", Aux::String("old")).unwrap();
        record.push_aux(b"BD", Aux::String("old")).unwrap();

        transform_record(&mut record, &model, &read_groups, true, false).unwrap();

        assert_eq!(record.qname(), b"read1");
        assert_eq!(record.seq().as_bytes(), b"AAC");
        assert_eq!(record.qual(), &[34, 34, 34]);
        assert_eq!(record.tid(), 2);
        assert_eq!(record.pos(), 123);
        assert_eq!(record.mtid(), 2);
        assert_eq!(record.mpos(), 456);
        assert_eq!(record.insert_size(), 333);
        assert_eq!(record.flags(), 0x41);
        assert_eq!(record.mapq(), 60);
        assert!(matches!(record.aux(b"RG").unwrap(), Aux::String("rg1")));
        assert!(matches!(record.aux(b"OQ").unwrap(), Aux::String("???")));
        assert!(matches!(record.aux(b"NM").unwrap(), Aux::I32(1)));
        assert!(record.aux(b"BI").is_err());
        assert!(record.aux(b"BD").is_err());
    }

    #[test]
    fn apply_bqsr_read_filter_rejects_cigar_n_records() {
        let mut record = Record::new();
        record.set(
            b"read1",
            Some(&CigarString(vec![
                Cigar::Match(2),
                Cigar::RefSkip(10),
                Cigar::Match(2),
            ])),
            b"AACC",
            &[30, 30, 30, 30],
        );

        assert!(!passes_apply_bqsr_read_filters(&record));
    }

    fn synthetic_table(
        read_group_rows: &[&str],
        quality_rows: &[&str],
        covariate_rows: &[&str],
    ) -> String {
        let mut text = String::new();
        text.push_str("#:GATKReport.v1.1:5\n");
        text.push_str("#:GATKTable:2:17:%s:%s:;\n");
        text.push_str(
            "#:GATKTable:Arguments:Recalibration argument collection values used in this run\n",
        );
        text.push_str("Argument                    Value\n");
        text.push_str("covariate                   ReadGroupCovariate,QualityScoreCovariate,ContextCovariate,CycleCovariate\n");
        text.push_str("low_quality_tail            2\n");
        text.push_str("maximum_cycle_value         500\n");
        text.push_str("mismatches_context_size     2\n");
        text.push_str("no_standard_covs            false\n\n");

        text.push_str("#:GATKTable:3:94:%d:%d:%d:;\n");
        text.push_str("#:GATKTable:Quantized:Quality quantization map\n");
        text.push_str("QualityScore  Count      QuantizedScore\n");
        for quality in 0..=MAX_RECALIBRATED_Q_SCORE {
            text.push_str(&format!("{quality}  0  {quality}\n"));
        }
        text.push('\n');

        text.push_str("#:GATKTable:6:1:%s:%s:%.4f:%.4f:%d:%.2f:;\n");
        text.push_str("#:GATKTable:RecalTable0:\n");
        text.push_str("ReadGroup        EventType  EmpiricalQuality  EstimatedQReported  Observations  Errors\n");
        for row in read_group_rows {
            text.push_str(row);
            text.push('\n');
        }
        text.push('\n');

        text.push_str("#:GATKTable:6:1:%s:%d:%s:%.4f:%d:%.2f:;\n");
        text.push_str("#:GATKTable:RecalTable1:\n");
        text.push_str(
            "ReadGroup        QualityScore  EventType  EmpiricalQuality  Observations  Errors\n",
        );
        for row in quality_rows {
            text.push_str(row);
            text.push('\n');
        }
        text.push('\n');

        text.push_str("#:GATKTable:8:1:%s:%d:%s:%s:%s:%.4f:%d:%.2f:;\n");
        text.push_str("#:GATKTable:RecalTable2:\n");
        text.push_str("ReadGroup        QualityScore  CovariateValue  CovariateName  EventType  EmpiricalQuality  Observations  Errors\n");
        for row in covariate_rows {
            text.push_str(row);
            text.push('\n');
        }
        text
    }
}
