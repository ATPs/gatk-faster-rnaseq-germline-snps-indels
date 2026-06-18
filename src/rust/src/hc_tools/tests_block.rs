#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

    #[test]
    fn vcf_compare_counts_private_shared_and_gt_diff() {
        let dir = tempdir().unwrap();
        let a = dir.path().join("a.vcf");
        let b = dir.path().join("b.vcf");
        fs::write(
            &a,
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\nchr1\t10\t.\tA\tG\t50\tPASS\tDP=10\tGT:DP\t0/1:10\nchr1\t20\t.\tAT\tA\t50\tPASS\tDP=5\tGT:DP\t1/1:5\n",
        )
        .unwrap();
        fs::write(
            &b,
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\nchr1\t10\t.\tA\tG\t50\tPASS\tDP=10\tGT:DP\t1/1:10\nchr1\t30\t.\tC\tT\t50\tPASS\tDP=6\tGT:DP\t0/1:6\n",
        )
        .unwrap();
        let comparison = compare_vcfs(&a, &b, "a", "b").unwrap();
        assert_eq!(comparison.pass_records.a_count, 2);
        assert_eq!(comparison.pass_records.b_count, 2);
        assert_eq!(comparison.pass_records.shared, 1);
        assert_eq!(comparison.pass_records.a_private, 1);
        assert_eq!(comparison.pass_records.b_private, 1);
        assert_eq!(comparison.pass_records.gt_diff, 1);
        assert_eq!(
            map_get(&comparison.pass_records.a_private_types, "INDEL_OR_COMPLEX"),
            1
        );
    }

    #[test]
    fn vcf_compare_splits_multiallelic_alt_keys() {
        let dir = tempdir().unwrap();
        let a = dir.path().join("a.vcf");
        let b = dir.path().join("b.vcf");
        fs::write(
            &a,
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\nchr1\t10\t.\tA\tC,G\t50\tPASS\t.\tGT\t1/2\n",
        )
        .unwrap();
        fs::write(
            &b,
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\nchr1\t10\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n",
        )
        .unwrap();
        let comparison = compare_vcfs(&a, &b, "a", "b").unwrap();
        assert_eq!(comparison.pass_records.a_count, 2);
        assert_eq!(comparison.pass_records.b_count, 1);
        assert_eq!(comparison.pass_records.shared, 1);
        assert_eq!(comparison.pass_records.a_private, 1);
        assert_eq!(comparison.pass_records.b_private, 0);
        assert_eq!(map_get(&comparison.pass_records.shared_types, "SNP"), 1);
        assert_eq!(map_get(&comparison.pass_records.a_private_types, "SNP"), 1);
    }

    #[test]
    fn region_selection_splits_categories() {
        let dir = tempdir().unwrap();
        let a = dir.path().join("a.vcf");
        let b = dir.path().join("b.vcf");
        fs::write(
            &a,
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\nchr1\t10\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\nchr1\t20\t.\tAT\tA\t50\tPASS\t.\tGT\t1/1\n",
        )
        .unwrap();
        fs::write(
            &b,
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\nchr1\t10\t.\tA\tG\t50\tPASS\t.\tGT\t1/1\nchr1\t30\t.\tC\tT\t50\tPASS\t.\tGT\t0/1\n",
        )
        .unwrap();
        let rows = select_regions(&a, &b, 5, 10, true).unwrap();
        let categories: BTreeSet<_> = rows.iter().map(|row| row.category.as_str()).collect();
        assert!(categories.contains("a_private_indel"));
        assert!(categories.contains("b_private_snp"));
        assert!(categories.contains("shared_gt_diff"));
    }

    #[test]
    fn stage_diff_accepts_numeric_tolerance() {
        let dir = tempdir().unwrap();
        let java = dir.path().join("java.tsv");
        let rust = dir.path().join("rust.tsv");
        fs::write(&java, "region\tread\tscore\nr1\tread1\t-10.0001\n").unwrap();
        fs::write(&rust, "region\tread\tscore\nr1\tread1\t-10.0002\n").unwrap();
        let config = StageDiffConfig {
            java_path: java,
            rust_path: rust,
            key_columns: vec!["region".to_string(), "read".to_string()],
            numeric_tolerance: 0.001,
            output_prefix: dir.path().join("diff"),
            stage_name: "pairhmm".to_string(),
        };
        let summary = run_stage_diff(&config).unwrap();
        assert_eq!(summary.shared_rows, 1);
        assert_eq!(summary.field_diffs, 0);
    }

    #[test]
    fn vcf_genotype_table_extracts_shared_fields() {
        let dir = tempdir().unwrap();
        let vcf = dir.path().join("calls.vcf");
        let output = dir.path().join("calls.tsv");
        fs::write(
            &vcf,
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\nchr1\t10\t.\tA\tG\t42.5\tPASS\tDP=12;FS=3.0;QD=4.2;DB\tGT:AD:DP:GQ:PL\t0/1:7,5:12:99:100,0,90\n",
        )
        .unwrap();
        let rows = write_vcf_genotype_table(&vcf, &output).unwrap();
        assert_eq!(rows, 1);
        let table = fs::read_to_string(output).unwrap();
        assert!(table
            .contains("chr1\t10\tA\tG\t42.5\tPASS\t0/1\t99\t12\t7\t5\t3.0\t4.2\t100,0,90\ttrue"));
    }
}
