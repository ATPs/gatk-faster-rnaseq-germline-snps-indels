#!/usr/bin/env zsh

set -euo pipefail

source /data/p/anaconda3/etc/profile.d/conda.sh
conda activate base
export PATH=/data/p/bin:$PATH
export LIBCLANG_PATH=/data/p/anaconda3/lib
export LD_LIBRARY_PATH=/data/p/anaconda3/lib:${LD_LIBRARY_PATH:-}
export BINDGEN_EXTRA_CLANG_ARGS='-I/usr/lib/gcc/x86_64-redhat-linux/11/include -I/usr/local/include -I/usr/include'

cd /data/p/gatk/gatk-faster-rnaseq-germline-snps-indels

ROUND_TAG=${ROUND_TAG:-round28}
OUT_DIR=${OUT_DIR:-/XCLabServer002_fastIO/gatk-faster-rnaseq/SRR949115_broad_hg38_run2/module6_good_compare_20260618_round28_targeted_n_active}
PREV_DIR=/XCLabServer002_fastIO/gatk-faster-rnaseq/SRR949115_broad_hg38_run2/module6_good_compare_20260618_round27_targeted_n_segments
JAVA_FILTERED=/XCLabServer002_fastIO/gatk-faster-rnaseq/SRR949115_broad_hg38_run2/module6_gatk_hc_rerun_20260616/SRR949115.gatk_hc_rerun.filtered.vcf.gz
RUST_RAW="${OUT_DIR}/${ROUND_TAG}.rust_hc.raw.vcf.gz"
RUST_FILTERED="${OUT_DIR}/${ROUND_TAG}.rust_hc.filtered.vcf.gz"

mkdir -p "${OUT_DIR}"
cp "${PREV_DIR}/round27.targets.tsv" "${OUT_DIR}/${ROUND_TAG}.targets.tsv"
cp "${PREV_DIR}/round27.targets.plus250.interval_list" "${OUT_DIR}/${ROUND_TAG}.targets.plus250.interval_list"

cargo build --manifest-path src/rust/Cargo.toml --release --bin rust_haplotype_caller

python src/step_haplotype_caller.py \
  --backend rust \
  --rust-bin src/rust/target/release/rust_haplotype_caller \
  --ref /data1/pub/gatk/broad_hg38/Homo_sapiens_assembly38.fasta \
  --input-bam /XCLabServer002_fastIO/gatk-faster-rnaseq/SRR949115_broad_hg38_run2/baseline/SRR949115.recal.bam \
  --input-interval-list "${OUT_DIR}/${ROUND_TAG}.targets.plus250.interval_list" \
  --output-vcf "${RUST_RAW}" \
  --dbsnp /data1/pub/gatk/broad_hg38/Homo_sapiens_assembly38.dbsnp138.vcf.gz \
  --threads 16 \
  --memory-gb 64 \
  --pair-hmm-threads 4 \
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

bcftools norm -m -any -Oz -o "${OUT_DIR}/java.filtered.split.vcf.gz" "${JAVA_FILTERED}"
bcftools norm -m -any -Oz -o "${OUT_DIR}/${ROUND_TAG}.rust_hc.filtered.split.vcf.gz" "${RUST_FILTERED}"
bcftools norm -m -any -Oz -o "${OUT_DIR}/${ROUND_TAG}.rust_hc.raw.split.vcf.gz" "${RUST_RAW}"
bcftools index -f -t "${OUT_DIR}/java.filtered.split.vcf.gz"
bcftools index -f -t "${OUT_DIR}/${ROUND_TAG}.rust_hc.filtered.split.vcf.gz"
bcftools index -f -t "${OUT_DIR}/${ROUND_TAG}.rust_hc.raw.split.vcf.gz"

OUT_DIR="${OUT_DIR}" ROUND_TAG="${ROUND_TAG}" python - <<'PY'
from pathlib import Path
import os
import subprocess

out_dir = Path(os.environ["OUT_DIR"])
round_tag = os.environ["ROUND_TAG"]
targets_path = out_dir / f"{round_tag}.targets.tsv"
filtered_vcf = out_dir / f"{round_tag}.rust_hc.filtered.split.vcf.gz"
raw_vcf = out_dir / f"{round_tag}.rust_hc.raw.split.vcf.gz"
status_path = out_dir / f"{round_tag}.target_exact_status.tsv"


def load_split_vcf(path):
    fmt = "%CHROM\t%POS\t%REF\t%ALT\t%FILTER[\t%GT\t%AD\t%DP]\n"
    result = subprocess.run(
        ["bcftools", "query", "-f", fmt, str(path)],
        check=True,
        text=True,
        capture_output=True,
    )
    records = {}
    for line in result.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) < 8:
            continue
        chrom, pos, ref, alt, filt, gt, ad, dp = fields[:8]
        records[(chrom, int(pos), ref, alt)] = {
            "filter": filt,
            "gt": gt,
            "ad": ad,
            "dp": dp,
        }
    return records


def alt_ad(ad):
    parts = ad.split(",")
    if len(parts) < 2 or parts[1] in {"", "."}:
        return None
    return int(parts[1])


def int_or_none(value):
    if value in {"", "."}:
        return None
    return int(value)


def classify(record):
    if record is None:
        return None
    dp = int_or_none(record["dp"])
    alt_depth = alt_ad(record["ad"])
    if record["filter"] != "PASS":
        return "exact_present_nonpass"
    if dp is None or dp < 10:
        return "exact_present_low_dp"
    if alt_depth is None or alt_depth < 3:
        return "exact_present_low_alt_ad"
    return "exact_present_good"


filtered = load_split_vcf(filtered_vcf)
raw = load_split_vcf(raw_vcf)

header = [
    "chrom",
    "pos",
    "ref",
    "alt",
    "filtered_present",
    "filtered_filter",
    "filtered_gt",
    "filtered_ad",
    "filtered_dp",
    "raw_present",
    "raw_filter",
    "raw_gt",
    "raw_ad",
    "raw_dp",
    "status",
]

with targets_path.open() as targets, status_path.open("w") as out:
    out.write("\t".join(header) + "\n")
    for line in targets:
        chrom, pos_s, ref, alt = line.rstrip("\n").split("\t")
        key = (chrom, int(pos_s), ref, alt)
        filtered_record = filtered.get(key)
        raw_record = raw.get(key)
        status = classify(filtered_record)
        if status is None:
            status = "exact_present_raw_only" if raw_record is not None else "exact_absent"
        row = [chrom, pos_s, ref, alt]
        for record in [filtered_record, raw_record]:
            row.append("1" if record is not None else "0")
            if record is None:
                row.extend(["", "", "", ""])
            else:
                row.extend([record["filter"], record["gt"], record["ad"], record["dp"]])
        row.append(status)
        out.write("\t".join(row) + "\n")
PY

awk 'NR > 1 { count[$NF]++ } END { for (status in count) print status "\t" count[status] }' \
  "${OUT_DIR}/${ROUND_TAG}.target_exact_status.tsv" \
  | sort > "${OUT_DIR}/${ROUND_TAG}.target_exact_status.counts.tsv"
