#!/usr/bin/env bash
# Join the AoU 100k person list to Dragen sex ploidy and v9 CRAM paths.
#
# Inputs (defaults match a Verily Workbench Jupyter terminal):
#   aou2_100k.tsv.gz
#     person_id
#   .../v9/wgs/short_read/snpindel/aux/qc/genomics_metrics.tsv
#     research_id, dragen_sex_ploidy (XX/XY, not survey sex_at_birth)
#   .../v9/wgs/cram/manifest.csv
#     person_id,cram_uri,cram_index_uri
#
# Output TSV columns:
#   person_id, dragen_sex_ploidy, sex, cram_uri, cram_index_uri
# sex is female/male for XX/XY and blank for any other ploidy.
# Every 100k id is kept. Missing metrics or CRAMs are blank.
#
#   ./scripts/join_aou2_100k_cram_sex.sh
#   ./scripts/join_aou2_100k_cram_sex.sh ~/aou2_100k.tsv.gz /path/metrics.tsv /path/manifest.csv ~/out.tsv

set -euo pipefail

IDS="${1:-${HOME}/aou2_100k.tsv.gz}"
METRICS="${2:-${HOME}/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/short_read/snpindel/aux/qc/genomics_metrics.tsv}"
MANIFEST="${3:-${HOME}/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/cram/manifest.csv}"
OUT="${4:-${HOME}/aou2_100k_cram_sex.tsv}"

for path in "${IDS}" "${METRICS}" "${MANIFEST}"; do
  if [[ ! -f "${path}" ]]; then
    echo "missing input: ${path}" >&2
    exit 1
  fi
done

tmpdir="$(mktemp -d)"
trap 'rm -rf "${tmpdir}"' EXIT

export LC_ALL=C

gzip -dc "${IDS}" | awk -F'\t' 'NR > 1 && $1 != "" { print $1 }' | sort -u > "${tmpdir}/ids"

# research_id and dragen_sex_ploidy by header name, so column order can move.
awk -F'\t' '
  NR == 1 {
    for (i = 1; i <= NF; i++) {
      if ($i == "research_id") id = i
      if ($i == "dragen_sex_ploidy") sex = i
    }
    if (!id || !sex) {
      print "genomics_metrics.tsv is missing research_id or dragen_sex_ploidy" > "/dev/stderr"
      exit 1
    }
    next
  }
  $id != "" { print $id "\t" $sex }
' "${METRICS}" | sort -k1,1 > "${tmpdir}/metrics"

awk -F',' 'NR > 1 && $1 != "" { print $1 "\t" $2 "\t" $3 }' "${MANIFEST}" | sort -k1,1 > "${tmpdir}/cram"

# One row per id. Later duplicate metric or CRAM rows are dropped.
awk -F'\t' 'seen[$1]++ { next } { print }' "${tmpdir}/metrics" > "${tmpdir}/metrics.uniq"
awk -F'\t' 'seen[$1]++ { next } { print }' "${tmpdir}/cram" > "${tmpdir}/cram.uniq"

{
  printf 'person_id\tdragen_sex_ploidy\tsex\tcram_uri\tcram_index_uri\n'
  join -t $'\t' -a 1 -e '' -o '0,2.2' "${tmpdir}/ids" "${tmpdir}/metrics.uniq" \
    | sort -k1,1 \
    | join -t $'\t' -a 1 -e '' -o '0,1.2,2.2,2.3' - "${tmpdir}/cram.uniq" \
    | awk -F'\t' '
        {
          sex = ""
          if ($2 == "XX") sex = "female"
          else if ($2 == "XY") sex = "male"
          print $1 "\t" $2 "\t" sex "\t" $3 "\t" $4
        }
      '
} > "${OUT}"

awk -F'\t' '
  NR == 1 { next }
  { n++ }
  $2 == "" { no_ploidy++ }
  $2 != "" && $2 != "XX" && $2 != "XY" { other++ }
  $4 == "" || $5 == "" { no_cram++ }
  END {
    printf "wrote %s\n", out
    printf "ids %d\n", n
    printf "missing dragen_sex_ploidy %d\n", no_ploidy + 0
    printf "ploidy other than XX/XY %d\n", other + 0
    printf "missing cram or crai %d\n", no_cram + 0
  }
' out="${OUT}" "${OUT}"
