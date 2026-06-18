#!/usr/bin/env zsh

set -euo pipefail

source /data/p/anaconda3/bin/activate base
export PATH=/data/p/bin:$PATH

cd /data/p/gatk/gatk-faster-rnaseq-germline-snps-indels

GOOD_DIR=/XCLabServer002_fastIO/gatk-faster-rnaseq/SRR949115_broad_hg38_run2/module6_good_compare_20260617_round6

./development/2026-06-17.good-set-sample-annotate.zsh \
  "${GOOD_DIR}/round6.java_only.sample20.tsv" \
  "${GOOD_DIR}/rust.split.vcf.gz" \
  "${GOOD_DIR}/round6.java_only.annotated.tsv"

./development/2026-06-17.good-set-sample-annotate.zsh \
  "${GOOD_DIR}/round6.rust_only.sample20.tsv" \
  "${GOOD_DIR}/java.split.vcf.gz" \
  "${GOOD_DIR}/round6.rust_only.annotated.tsv"

echo
echo "round6 summary:"
awk -F'\t' 'NR == 1 || $1 == "pass"' "${GOOD_DIR}/java_vs_rust.good.summary.tsv"

echo
echo "round6 java_only gate counts:"
awk -F'\t' 'NR>1 {c[$12]++} END {for (k in c) print k, c[k]}' "${GOOD_DIR}/round6.java_only.annotated.tsv" | sort

echo
echo "round6 rust_only gate counts:"
awk -F'\t' 'NR>1 {c[$12]++} END {for (k in c) print k, c[k]}' "${GOOD_DIR}/round6.rust_only.annotated.tsv" | sort

echo
for region in \
  chr4:176328436-176328436 \
  chr18:23482513-23482513 \
  chr6:21596269-21596269 \
  chr3:31532901-31532901 \
  chr3:152463139-152463139 \
  chr1:192716486-192716486 \
  chr2:203304551-203304551 \
  chr16:2232255-2232255 \
  chr18:13745168-13745168 \
  chrX:153796101-153796101
do
  echo "=== ${region} RUST ==="
  bcftools query -r "${region}" -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER[\t%GT\t%AD\t%DP\t%GQ]\n' "${GOOD_DIR}/rust.split.vcf.gz"
  echo "=== ${region} JAVA ==="
  bcftools query -r "${region}" -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER[\t%GT\t%AD\t%DP\t%GQ]\n' "${GOOD_DIR}/java.split.vcf.gz"
done
