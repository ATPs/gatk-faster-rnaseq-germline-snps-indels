#!/usr/bin/env zsh

set -euo pipefail

source /data/p/anaconda3/etc/profile.d/conda.sh
conda activate base
export PATH=/data/p/bin:$PATH
export LIBCLANG_PATH=/data/p/anaconda3/lib
export LD_LIBRARY_PATH=/data/p/anaconda3/lib:${LD_LIBRARY_PATH:-}
export BINDGEN_EXTRA_CLANG_ARGS='-I/usr/lib/gcc/x86_64-redhat-linux/11/include -I/usr/local/include -I/usr/include'

cd /data/p/gatk/gatk-faster-rnaseq-germline-snps-indels

GOOD_DIR=/XCLabServer002_fastIO/gatk-faster-rnaseq/SRR949115_broad_hg38_run2/module6_good_compare_20260618_round12
RUST_RAW="${GOOD_DIR}/SRR949115.rust_hc.raw.vcf.gz"
RUST_FILTERED="${GOOD_DIR}/SRR949115.rust_hc.filtered.vcf.gz"
JAVA_FILTERED=/XCLabServer002_fastIO/gatk-faster-rnaseq/SRR949115_broad_hg38_run2/module6_gatk_hc_rerun_20260616/SRR949115.gatk_hc_rerun.filtered.vcf.gz

mkdir -p "${GOOD_DIR}"

python src/step_haplotype_caller.py \
  --backend rust \
  --rust-bin src/rust/target/release/rust_haplotype_caller \
  --ref /data1/pub/gatk/broad_hg38/Homo_sapiens_assembly38.fasta \
  --input-bam /XCLabServer002_fastIO/gatk-faster-rnaseq/SRR949115_broad_hg38_run2/baseline/SRR949115.recal.bam \
  --input-interval-list /XCLabServer002_fastIO/gatk-faster-rnaseq/SRR949115_broad_hg38_run2/baseline/reference/exons.interval_list \
  --output-vcf "${RUST_RAW}" \
  --dbsnp /data1/pub/gatk/broad_hg38/Homo_sapiens_assembly38.dbsnp138.vcf.gz \
  --threads 40 \
  --memory-gb 128 \
  --pair-hmm-threads 8 \
  --pair-hmm-implementation native

/data/p/gatk/gatk-4.6.2.0/gatk --java-options "-Xmx4g" VariantFiltration \
  -R /data1/pub/gatk/broad_hg38/Homo_sapiens_assembly38.fasta \
  -V "${RUST_RAW}" \
  --window 35 \
  --cluster 3 \
  --filter-name FS \
  --filter "FS > 30.0" \
  --filter-name QD \
  --filter "QD < 2.0" \
  -O "${RUST_FILTERED}"

bcftools norm -m -any -Oz -o "${GOOD_DIR}/java.split.vcf.gz" "${JAVA_FILTERED}"
bcftools norm -m -any -Oz -o "${GOOD_DIR}/rust.split.vcf.gz" "${RUST_FILTERED}"
bcftools index -f -t "${GOOD_DIR}/java.split.vcf.gz"
bcftools index -f -t "${GOOD_DIR}/rust.split.vcf.gz"

bcftools view -f PASS -i 'FORMAT/DP[0]>=10 && FORMAT/AD[0:1]>=3' -Oz -o "${GOOD_DIR}/java.good.vcf.gz" "${GOOD_DIR}/java.split.vcf.gz"
bcftools view -f PASS -i 'FORMAT/DP[0]>=10 && FORMAT/AD[0:1]>=3' -Oz -o "${GOOD_DIR}/rust.good.vcf.gz" "${GOOD_DIR}/rust.split.vcf.gz"
bcftools index -f -t "${GOOD_DIR}/java.good.vcf.gz"
bcftools index -f -t "${GOOD_DIR}/rust.good.vcf.gz"

src/rust/target/release/rust_hc_vcf_compare \
  --a-vcf "${GOOD_DIR}/java.good.vcf.gz" \
  --b-vcf "${GOOD_DIR}/rust.good.vcf.gz" \
  --a-label java_good \
  --b-label rust_good \
  --output-prefix "${GOOD_DIR}/java_vs_rust.good"

./development/2026-06-17.good-compare-postprocess.zsh "${GOOD_DIR}" round12 20
