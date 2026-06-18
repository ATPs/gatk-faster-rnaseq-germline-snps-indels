#!/usr/bin/env zsh

set -euo pipefail

if [[ $# -ne 3 ]]; then
  echo "usage: $0 SAMPLE_TSV OTHER_SPLIT_VCF OUTPUT_TSV" >&2
  exit 2
fi

SAMPLE_TSV=$1
OTHER_SPLIT_VCF=$2
OUTPUT_TSV=$3
MIN_DP=10
MIN_ALT_AD=3

extract_alt_ad() {
  local ad=$1
  if [[ -z "${ad}" || "${ad}" == "." ]]; then
    echo ""
    return
  fi
  echo "${ad}" | awk -F',' '{print (NF >= 2 ? $2 : "")}'
}

classify_gate() {
  local filter=$1
  local dp=$2
  local alt_ad=$3

  if [[ -z "${filter}" ]]; then
    echo "exact_absent"
  elif [[ "${filter}" != "PASS" ]]; then
    echo "exact_present_nonpass"
  elif [[ -z "${dp}" || "${dp}" == "." || ${dp} -lt ${MIN_DP} ]]; then
    echo "exact_present_low_dp"
  elif [[ -z "${alt_ad}" || "${alt_ad}" == "." || ${alt_ad} -lt ${MIN_ALT_AD} ]]; then
    echo "exact_present_low_alt_ad"
  else
    echo "exact_present_good"
  fi
}

mkdir -p "$(dirname "${OUTPUT_TSV}")"
{
  printf 'category\tchrom\tpos\tref\talt\ttype\tother_exact_filter\tother_exact_gt\tother_exact_ad\tother_exact_dp\tother_exact_alt_ad\tother_gate\tother_nearby_count\tother_nearby_summary\n'

  while IFS=$'\t' read -r category chrom pos ref alt type a_filter b_filter a_gt b_gt; do
    [[ -z "${category}" ]] && continue
    local_region="${chrom}:${pos}-${pos}"
    exact_record=$(
      bcftools query -r "${local_region}" -f '%REF\t%ALT\t%FILTER[\t%GT\t%AD\t%DP]\n' "${OTHER_SPLIT_VCF}" \
        | awk -F'\t' -v ref="${ref}" -v alt="${alt}" '$1 == ref && $2 == alt {print; exit}'
    )

    other_filter=""
    other_gt=""
    other_ad=""
    other_dp=""
    other_alt_ad=""
    if [[ -n "${exact_record}" ]]; then
      other_filter=$(echo "${exact_record}" | cut -f3)
      other_gt=$(echo "${exact_record}" | cut -f4)
      other_ad=$(echo "${exact_record}" | cut -f5)
      other_dp=$(echo "${exact_record}" | cut -f6)
      other_alt_ad=$(extract_alt_ad "${other_ad}")
    fi

    other_gate=$(classify_gate "${other_filter}" "${other_dp:-0}" "${other_alt_ad:-0}")

    start=$(( pos > 100 ? pos - 100 : 1 ))
    end=$(( pos + 100 ))
    nearby_summary=$(
      bcftools query -r "${chrom}:${start}-${end}" -f '%POS:%REF:%ALT:%FILTER[%GT:%AD:%DP]\n' "${OTHER_SPLIT_VCF}" \
        | paste -sd';' -
    )
    nearby_count=0
    if [[ -n "${nearby_summary}" ]]; then
      nearby_count=$(echo "${nearby_summary}" | awk -F';' '{print NF}')
    fi

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "${category}" "${chrom}" "${pos}" "${ref}" "${alt}" "${type}" \
      "${other_filter}" "${other_gt}" "${other_ad}" "${other_dp}" "${other_alt_ad}" \
      "${other_gate}" "${nearby_count}" "${nearby_summary}"
  done < "${SAMPLE_TSV}"
} > "${OUTPUT_TSV}"
