#!/usr/bin/env zsh

set -euo pipefail

if [[ $# -lt 4 || $# -gt 5 ]]; then
  echo "usage: $0 DETAILS_TSV RAW_SPLIT_VCF FILTERED_SPLIT_VCF OUTPUT_TSV [CATEGORY]" >&2
  exit 2
fi

DETAILS_TSV=$1
RAW_SPLIT_VCF=$2
FILTERED_SPLIT_VCF=$3
OUTPUT_TSV=$4
CATEGORY_FILTER=${5:-a_private}
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

classify_good_gate() {
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

summarize_same_pos() {
  local vcf=$1
  local chrom=$2
  local pos=$3
  bcftools query \
    -r "${chrom}:${pos}-${pos}" \
    -f '%REF\t%ALT\t%QUAL\t%FILTER[\t%GT\t%AD\t%DP\t%GQ]\n' \
    "${vcf}" \
    | awk -F'\t' '
        BEGIN { sep = "" }
        {
          printf "%s%s>%s:%s:%s:%s:%s:%s:%s",
            sep, $1, $2, $3, $4, $5, $6, $7, $8;
          sep = ";";
        }
        END { print "" }
      '
}

summarize_nearby_nonexact() {
  local vcf=$1
  local chrom=$2
  local pos=$3
  local ref=$4
  local alt=$5
  local start=$(( pos > 100 ? pos - 100 : 1 ))
  local end=$(( pos + 100 ))

  bcftools query \
    -r "${chrom}:${start}-${end}" \
    -f '%POS\t%REF\t%ALT\t%QUAL\t%FILTER[\t%GT\t%AD\t%DP\t%GQ]\n' \
    "${vcf}" \
    | awk -F'\t' -v pos="${pos}" -v ref="${ref}" -v alt="${alt}" '
        BEGIN { sep = "" }
        !($1 == pos && $2 == ref && $3 == alt) {
          printf "%s%s:%s>%s:%s:%s:%s:%s:%s:%s",
            sep, $1, $2, $3, $4, $5, $6, $7, $8, $9;
          sep = ";";
        }
        END { print "" }
      '
}

mkdir -p "$(dirname "${OUTPUT_TSV}")"
{
  printf 'category\tchrom\tpos\tref\talt\ttype\tjava_good_filter\tjava_good_gt'
  printf '\trust_raw_exact_present\trust_raw_exact_filter\trust_raw_exact_qual\trust_raw_exact_gt\trust_raw_exact_ad\trust_raw_exact_dp\trust_raw_exact_alt_ad\trust_raw_exact_gq\trust_raw_exact_pl'
  printf '\trust_raw_info_dp\trust_raw_info_ac\trust_raw_info_af\trust_raw_info_an\trust_raw_info_fs\trust_raw_info_qd'
  printf '\trust_filtered_exact_present\trust_filtered_exact_filter\trust_filtered_exact_qual\trust_filtered_exact_gt\trust_filtered_exact_ad\trust_filtered_exact_dp\trust_filtered_exact_alt_ad\trust_filtered_exact_gq\trust_filtered_exact_pl'
  printf '\trust_filtered_info_dp\trust_filtered_info_ac\trust_filtered_info_af\trust_filtered_info_an\trust_filtered_info_fs\trust_filtered_info_qd\trust_filtered_good_gate'
  printf '\trust_raw_same_pos_count\trust_raw_same_pos_summary\trust_raw_nearby_nonexact_count\trust_raw_nearby_nonexact_summary\n'

  while IFS=$'\t' read -r category chrom pos ref alt type a_filter b_filter a_gt b_gt; do
    [[ -z "${category}" || "${category}" == "category" ]] && continue
    [[ "${category}" != "${CATEGORY_FILTER}" ]] && continue

    region="${chrom}:${pos}-${pos}"

    raw_exact=$(
      bcftools query \
        -r "${region}" \
        -f '%REF\t%ALT\t%FILTER\t%QUAL[\t%GT\t%AD\t%DP\t%GQ\t%PL]\t%INFO/DP\t%INFO/AC\t%INFO/AF\t%INFO/AN\t%INFO/FS\t%INFO/QD\n' \
        "${RAW_SPLIT_VCF}" \
        | awk -F'\t' -v ref="${ref}" -v alt="${alt}" '$1 == ref && $2 == alt {print; exit}'
    )

    raw_exact_present=0
    raw_filter=""
    raw_qual=""
    raw_gt=""
    raw_ad=""
    raw_dp=""
    raw_alt_ad=""
    raw_gq=""
    raw_pl=""
    raw_info_dp=""
    raw_info_ac=""
    raw_info_af=""
    raw_info_an=""
    raw_info_fs=""
    raw_info_qd=""
    if [[ -n "${raw_exact}" ]]; then
      raw_exact_present=1
      IFS=$'\t' read -r raw_ref raw_alt raw_filter raw_qual raw_gt raw_ad raw_dp raw_gq raw_pl raw_info_dp raw_info_ac raw_info_af raw_info_an raw_info_fs raw_info_qd <<< "${raw_exact}"
      raw_alt_ad=$(extract_alt_ad "${raw_ad}")
    fi

    filtered_exact=$(
      bcftools query \
        -r "${region}" \
        -f '%REF\t%ALT\t%FILTER\t%QUAL[\t%GT\t%AD\t%DP\t%GQ\t%PL]\t%INFO/DP\t%INFO/AC\t%INFO/AF\t%INFO/AN\t%INFO/FS\t%INFO/QD\n' \
        "${FILTERED_SPLIT_VCF}" \
        | awk -F'\t' -v ref="${ref}" -v alt="${alt}" '$1 == ref && $2 == alt {print; exit}'
    )

    filtered_exact_present=0
    filtered_filter=""
    filtered_qual=""
    filtered_gt=""
    filtered_ad=""
    filtered_dp=""
    filtered_alt_ad=""
    filtered_gq=""
    filtered_pl=""
    filtered_info_dp=""
    filtered_info_ac=""
    filtered_info_af=""
    filtered_info_an=""
    filtered_info_fs=""
    filtered_info_qd=""
    filtered_good_gate="exact_absent"
    if [[ -n "${filtered_exact}" ]]; then
      filtered_exact_present=1
      IFS=$'\t' read -r filtered_ref filtered_alt filtered_filter filtered_qual filtered_gt filtered_ad filtered_dp filtered_gq filtered_pl filtered_info_dp filtered_info_ac filtered_info_af filtered_info_an filtered_info_fs filtered_info_qd <<< "${filtered_exact}"
      filtered_alt_ad=$(extract_alt_ad "${filtered_ad}")
      filtered_good_gate=$(classify_good_gate "${filtered_filter}" "${filtered_dp}" "${filtered_alt_ad}")
    fi

    raw_same_pos_summary=$(summarize_same_pos "${RAW_SPLIT_VCF}" "${chrom}" "${pos}")
    raw_same_pos_count=0
    if [[ -n "${raw_same_pos_summary}" ]]; then
      raw_same_pos_count=$(echo "${raw_same_pos_summary}" | awk -F';' '{print NF}')
    fi

    raw_nearby_nonexact_summary=$(summarize_nearby_nonexact "${RAW_SPLIT_VCF}" "${chrom}" "${pos}" "${ref}" "${alt}")
    raw_nearby_nonexact_count=0
    if [[ -n "${raw_nearby_nonexact_summary}" ]]; then
      raw_nearby_nonexact_count=$(echo "${raw_nearby_nonexact_summary}" | awk -F';' '{print NF}')
    fi

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' \
      "${category}" "${chrom}" "${pos}" "${ref}" "${alt}" "${type}" "${a_filter}" "${a_gt}"
    printf '\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' \
      "${raw_exact_present}" "${raw_filter}" "${raw_qual}" "${raw_gt}" "${raw_ad}" "${raw_dp}" "${raw_alt_ad}" "${raw_gq}" "${raw_pl}"
    printf '\t%s\t%s\t%s\t%s\t%s\t%s' \
      "${raw_info_dp}" "${raw_info_ac}" "${raw_info_af}" "${raw_info_an}" "${raw_info_fs}" "${raw_info_qd}"
    printf '\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' \
      "${filtered_exact_present}" "${filtered_filter}" "${filtered_qual}" "${filtered_gt}" "${filtered_ad}" "${filtered_dp}" "${filtered_alt_ad}" "${filtered_gq}" "${filtered_pl}"
    printf '\t%s\t%s\t%s\t%s\t%s\t%s\t%s' \
      "${filtered_info_dp}" "${filtered_info_ac}" "${filtered_info_af}" "${filtered_info_an}" "${filtered_info_fs}" "${filtered_info_qd}" "${filtered_good_gate}"
    printf '\t%s\t%s\t%s\t%s\n' \
      "${raw_same_pos_count}" "${raw_same_pos_summary}" "${raw_nearby_nonexact_count}" "${raw_nearby_nonexact_summary}"
  done < "${DETAILS_TSV}"
} > "${OUTPUT_TSV}"
