#!/usr/bin/env zsh

set -euo pipefail

if [[ $# -lt 2 || $# -gt 3 ]]; then
  echo "usage: $0 GOOD_DIR ROUND_LABEL [SAMPLE_SIZE]" >&2
  exit 2
fi

GOOD_DIR=$1
ROUND_LABEL=$2
SAMPLE_SIZE=${3:-20}
DETAILS="${GOOD_DIR}/java_vs_rust.good.details.tsv"
SUMMARY="${GOOD_DIR}/java_vs_rust.good.summary.tsv"
JAVA_SAMPLE="${GOOD_DIR}/${ROUND_LABEL}.java_only.sample${SAMPLE_SIZE}.tsv"
RUST_SAMPLE="${GOOD_DIR}/${ROUND_LABEL}.rust_only.sample${SAMPLE_SIZE}.tsv"

if [[ ! -f "${SUMMARY}" ]]; then
  echo "missing summary file: ${SUMMARY}" >&2
  exit 1
fi

if [[ ! -f "${DETAILS}" ]]; then
  echo "missing detail file: ${DETAILS}" >&2
  exit 1
fi

echo "pass summary:"
awk -F'\t' 'NR == 1 || $1 == "pass"' "${SUMMARY}"

awk -F'\t' -v n="${SAMPLE_SIZE}" 'NR > 1 && $1 == "b_private"' "${DETAILS}" | shuf -n "${SAMPLE_SIZE}" > "${RUST_SAMPLE}"
awk -F'\t' -v n="${SAMPLE_SIZE}" 'NR > 1 && $1 == "a_private"' "${DETAILS}" | shuf -n "${SAMPLE_SIZE}" > "${JAVA_SAMPLE}"

echo
echo "wrote:"
echo "${RUST_SAMPLE}"
echo "${JAVA_SAMPLE}"
