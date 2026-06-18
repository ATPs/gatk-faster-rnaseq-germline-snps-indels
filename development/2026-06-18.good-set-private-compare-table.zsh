#!/usr/bin/env zsh

set -euo pipefail

if [[ $# -ne 8 ]]; then
  echo "usage: $0 DETAILS_TSV PRIMARY_GOOD_SPLIT_VCF OTHER_RAW_SPLIT_VCF OTHER_FILTERED_SPLIT_VCF OUTPUT_TSV CATEGORY PRIMARY_LABEL OTHER_LABEL" >&2
  exit 2
fi

DETAILS_TSV=$1
PRIMARY_GOOD_SPLIT_VCF=$2
OTHER_RAW_SPLIT_VCF=$3
OTHER_FILTERED_SPLIT_VCF=$4
OUTPUT_TSV=$5
CATEGORY_FILTER=$6
PRIMARY_LABEL=$7
OTHER_LABEL=$8
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

query_exact_record() {
  local vcf=$1
  local chrom=$2
  local pos=$3
  local ref=$4
  local alt=$5
  bcftools query \
    -r "${chrom}:${pos}-${pos}" \
    -f '%REF\t%ALT\t%FILTER\t%QUAL[\t%GT\t%AD\t%DP\t%GQ\t%PL]\t%INFO/DP\t%INFO/AC\t%INFO/AF\t%INFO/AN\t%INFO/FS\t%INFO/QD\n' \
    "${vcf}" \
    | awk -F'\t' -v ref="${ref}" -v alt="${alt}" '$1 == ref && $2 == alt {print; exit}'
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
  printf 'category\tchrom\tpos\tref\talt\ttype'
  printf '\t%s_good_exact_present\t%s_good_filter\t%s_good_qual\t%s_good_gt\t%s_good_ad\t%s_good_dp\t%s_good_alt_ad\t%s_good_gq\t%s_good_pl' \
    "${PRIMARY_LABEL}" "${PRIMARY_LABEL}" "${PRIMARY_LABEL}" "${PRIMARY_LABEL}" "${PRIMARY_LABEL}" "${PRIMARY_LABEL}" "${PRIMARY_LABEL}" "${PRIMARY_LABEL}" "${PRIMARY_LABEL}"
  printf '\t%s_good_info_dp\t%s_good_info_ac\t%s_good_info_af\t%s_good_info_an\t%s_good_info_fs\t%s_good_info_qd' \
    "${PRIMARY_LABEL}" "${PRIMARY_LABEL}" "${PRIMARY_LABEL}" "${PRIMARY_LABEL}" "${PRIMARY_LABEL}" "${PRIMARY_LABEL}"
  printf '\t%s_raw_exact_present\t%s_raw_exact_filter\t%s_raw_exact_qual\t%s_raw_exact_gt\t%s_raw_exact_ad\t%s_raw_exact_dp\t%s_raw_exact_alt_ad\t%s_raw_exact_gq\t%s_raw_exact_pl' \
    "${OTHER_LABEL}" "${OTHER_LABEL}" "${OTHER_LABEL}" "${OTHER_LABEL}" "${OTHER_LABEL}" "${OTHER_LABEL}" "${OTHER_LABEL}" "${OTHER_LABEL}" "${OTHER_LABEL}"
  printf '\t%s_raw_info_dp\t%s_raw_info_ac\t%s_raw_info_af\t%s_raw_info_an\t%s_raw_info_fs\t%s_raw_info_qd' \
    "${OTHER_LABEL}" "${OTHER_LABEL}" "${OTHER_LABEL}" "${OTHER_LABEL}" "${OTHER_LABEL}" "${OTHER_LABEL}"
  printf '\t%s_filtered_exact_present\t%s_filtered_exact_filter\t%s_filtered_exact_qual\t%s_filtered_exact_gt\t%s_filtered_exact_ad\t%s_filtered_exact_dp\t%s_filtered_exact_alt_ad\t%s_filtered_exact_gq\t%s_filtered_exact_pl' \
    "${OTHER_LABEL}" "${OTHER_LABEL}" "${OTHER_LABEL}" "${OTHER_LABEL}" "${OTHER_LABEL}" "${OTHER_LABEL}" "${OTHER_LABEL}" "${OTHER_LABEL}" "${OTHER_LABEL}"
  printf '\t%s_filtered_info_dp\t%s_filtered_info_ac\t%s_filtered_info_af\t%s_filtered_info_an\t%s_filtered_info_fs\t%s_filtered_info_qd\t%s_filtered_good_gate' \
    "${OTHER_LABEL}" "${OTHER_LABEL}" "${OTHER_LABEL}" "${OTHER_LABEL}" "${OTHER_LABEL}" "${OTHER_LABEL}" "${OTHER_LABEL}"
  printf '\t%s_raw_same_pos_count\t%s_raw_same_pos_summary\t%s_raw_nearby_nonexact_count\t%s_raw_nearby_nonexact_summary\n' \
    "${OTHER_LABEL}" "${OTHER_LABEL}" "${OTHER_LABEL}" "${OTHER_LABEL}"

  while IFS=$'\t' read -r category chrom pos ref alt type a_filter b_filter a_gt b_gt; do
    [[ -z "${category}" || "${category}" == "category" ]] && continue
    [[ "${category}" != "${CATEGORY_FILTER}" ]] && continue

    primary_exact=$(query_exact_record "${PRIMARY_GOOD_SPLIT_VCF}" "${chrom}" "${pos}" "${ref}" "${alt}")
    primary_exact_present=0
    primary_filter=""
    primary_qual=""
    primary_gt=""
    primary_ad=""
    primary_dp=""
    primary_alt_ad=""
    primary_gq=""
    primary_pl=""
    primary_info_dp=""
    primary_info_ac=""
    primary_info_af=""
    primary_info_an=""
    primary_info_fs=""
    primary_info_qd=""
    if [[ -n "${primary_exact}" ]]; then
      primary_exact_present=1
      IFS=$'\t' read -r primary_ref primary_alt primary_filter primary_qual primary_gt primary_ad primary_dp primary_gq primary_pl primary_info_dp primary_info_ac primary_info_af primary_info_an primary_info_fs primary_info_qd <<< "${primary_exact}"
      primary_alt_ad=$(extract_alt_ad "${primary_ad}")
    fi

    other_raw_exact=$(query_exact_record "${OTHER_RAW_SPLIT_VCF}" "${chrom}" "${pos}" "${ref}" "${alt}")
    other_raw_exact_present=0
    other_raw_filter=""
    other_raw_qual=""
    other_raw_gt=""
    other_raw_ad=""
    other_raw_dp=""
    other_raw_alt_ad=""
    other_raw_gq=""
    other_raw_pl=""
    other_raw_info_dp=""
    other_raw_info_ac=""
    other_raw_info_af=""
    other_raw_info_an=""
    other_raw_info_fs=""
    other_raw_info_qd=""
    if [[ -n "${other_raw_exact}" ]]; then
      other_raw_exact_present=1
      IFS=$'\t' read -r other_raw_ref other_raw_alt other_raw_filter other_raw_qual other_raw_gt other_raw_ad other_raw_dp other_raw_gq other_raw_pl other_raw_info_dp other_raw_info_ac other_raw_info_af other_raw_info_an other_raw_info_fs other_raw_info_qd <<< "${other_raw_exact}"
      other_raw_alt_ad=$(extract_alt_ad "${other_raw_ad}")
    fi

    other_filtered_exact=$(query_exact_record "${OTHER_FILTERED_SPLIT_VCF}" "${chrom}" "${pos}" "${ref}" "${alt}")
    other_filtered_exact_present=0
    other_filtered_filter=""
    other_filtered_qual=""
    other_filtered_gt=""
    other_filtered_ad=""
    other_filtered_dp=""
    other_filtered_alt_ad=""
    other_filtered_gq=""
    other_filtered_pl=""
    other_filtered_info_dp=""
    other_filtered_info_ac=""
    other_filtered_info_af=""
    other_filtered_info_an=""
    other_filtered_info_fs=""
    other_filtered_info_qd=""
    other_filtered_good_gate="exact_absent"
    if [[ -n "${other_filtered_exact}" ]]; then
      other_filtered_exact_present=1
      IFS=$'\t' read -r other_filtered_ref other_filtered_alt other_filtered_filter other_filtered_qual other_filtered_gt other_filtered_ad other_filtered_dp other_filtered_gq other_filtered_pl other_filtered_info_dp other_filtered_info_ac other_filtered_info_af other_filtered_info_an other_filtered_info_fs other_filtered_info_qd <<< "${other_filtered_exact}"
      other_filtered_alt_ad=$(extract_alt_ad "${other_filtered_ad}")
      other_filtered_good_gate=$(classify_good_gate "${other_filtered_filter}" "${other_filtered_dp}" "${other_filtered_alt_ad}")
    fi

    other_raw_same_pos_summary=$(summarize_same_pos "${OTHER_RAW_SPLIT_VCF}" "${chrom}" "${pos}")
    other_raw_same_pos_count=0
    if [[ -n "${other_raw_same_pos_summary}" ]]; then
      other_raw_same_pos_count=$(echo "${other_raw_same_pos_summary}" | awk -F';' '{print NF}')
    fi

    other_raw_nearby_nonexact_summary=$(summarize_nearby_nonexact "${OTHER_RAW_SPLIT_VCF}" "${chrom}" "${pos}" "${ref}" "${alt}")
    other_raw_nearby_nonexact_count=0
    if [[ -n "${other_raw_nearby_nonexact_summary}" ]]; then
      other_raw_nearby_nonexact_count=$(echo "${other_raw_nearby_nonexact_summary}" | awk -F';' '{print NF}')
    fi

    printf '%s\t%s\t%s\t%s\t%s\t%s' \
      "${category}" "${chrom}" "${pos}" "${ref}" "${alt}" "${type}"
    printf '\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' \
      "${primary_exact_present}" "${primary_filter}" "${primary_qual}" "${primary_gt}" "${primary_ad}" "${primary_dp}" "${primary_alt_ad}" "${primary_gq}" "${primary_pl}"
    printf '\t%s\t%s\t%s\t%s\t%s\t%s' \
      "${primary_info_dp}" "${primary_info_ac}" "${primary_info_af}" "${primary_info_an}" "${primary_info_fs}" "${primary_info_qd}"
    printf '\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' \
      "${other_raw_exact_present}" "${other_raw_filter}" "${other_raw_qual}" "${other_raw_gt}" "${other_raw_ad}" "${other_raw_dp}" "${other_raw_alt_ad}" "${other_raw_gq}" "${other_raw_pl}"
    printf '\t%s\t%s\t%s\t%s\t%s\t%s' \
      "${other_raw_info_dp}" "${other_raw_info_ac}" "${other_raw_info_af}" "${other_raw_info_an}" "${other_raw_info_fs}" "${other_raw_info_qd}"
    printf '\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' \
      "${other_filtered_exact_present}" "${other_filtered_filter}" "${other_filtered_qual}" "${other_filtered_gt}" "${other_filtered_ad}" "${other_filtered_dp}" "${other_filtered_alt_ad}" "${other_filtered_gq}" "${other_filtered_pl}"
    printf '\t%s\t%s\t%s\t%s\t%s\t%s\t%s' \
      "${other_filtered_info_dp}" "${other_filtered_info_ac}" "${other_filtered_info_af}" "${other_filtered_info_an}" "${other_filtered_info_fs}" "${other_filtered_info_qd}" "${other_filtered_good_gate}"
    printf '\t%s\t%s\t%s\t%s\n' \
      "${other_raw_same_pos_count}" "${other_raw_same_pos_summary}" "${other_raw_nearby_nonexact_count}" "${other_raw_nearby_nonexact_summary}"
  done < "${DETAILS_TSV}"
} > "${OUTPUT_TSV}"
