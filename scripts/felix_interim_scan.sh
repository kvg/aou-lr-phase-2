#!/usr/bin/env bash
# Merge the per-shard result TSVs a FelixGenome run has produced so far, without
# waiting for Concat and Summarize.
#
# Step-2 shards are written as they finish, so a run that is still going — or one
# that failed before its collection tasks — already holds usable results. This
# reports per-phenotype chromosome coverage and merges whichever phenotypes are
# complete, streaming through `gsutil cat` so nothing but the merged TSVs is
# downloaded.
#
# The step-2 scatters are nested, call-<Name>/shard-<chrom>/shard-<pheno>/, and
# the inner index is the line number in selected_phenotypes. The outer index is
# not needed: the chromosome is in the CHROM column of the data.
#
# Writes into <OUT>/:
#   coverage.tsv                   phenotype  index  chromosomes_present  complete
#   <phenotype>.merged.felix.tsv.gz   one per complete phenotype
#
# --run is repeatable (and accepts a comma-separated list), because the shards
# of one logical run can be spread across several submissions. A resubmission
# that call-cached a task does not necessarily re-create that task's directory,
# so the outputs stay under whichever run actually executed them. Give the
# original run and the rescue; shards are unioned on their path below
# call-<Name>/, and where both hold the same shard the LAST --run wins.
#
# Usage:
#   ./scripts/felix_interim_scan.sh --run gs://.../FelixGenome/<wid> \
#       --phenotypes gs://BUCKET/.../selected_phenotypes.txt --arm main --out interim_main
#   ... --run gs://.../<original-wid> --run gs://.../<rescue-wid>   # union both
#   ... --arm ru          # repeat-dosage arm (RuStep2)
#   ... --arm ru_biallelic # carrier-coded comparator (BaStep2)
#   ... --min-chroms 22   # merge only phenotypes with this many chromosomes (default 22)
#   ... --coverage-only   # report coverage, merge nothing
#
# Set GOOGLE_PROJECT to bill a requester-pays bucket (adds gsutil -u).
#
# Requires: bash, gsutil, awk.

set -euo pipefail

ROOTS=()
PHENOTYPES=""
ARM="main"
OUT=""
MIN_CHROMS=22
COVERAGE_ONLY=0
LISTING_DIR=""

usage() { sed -n '2,31p' "$0" | sed 's/^# \?//'; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    -r|--run)
      # Repeatable, and accepts a comma-separated list. Results are spread
      # across runs: a resubmission that call-cached a task does not
      # necessarily re-create that task's directory, so the shards stay under
      # whichever run actually executed them.
      _old_ifs="${IFS}"; IFS=','
      for _r in $2; do [[ -n "${_r}" ]] && ROOTS+=("${_r%/}"); done
      IFS="${_old_ifs}"; shift ;;
    -p|--phenotypes) PHENOTYPES="$2"; shift ;;
    -a|--arm) ARM="$2"; shift ;;
    -o|--out) OUT="$2"; shift ;;
    -m|--min-chroms) MIN_CHROMS="$2"; shift ;;
    --coverage-only) COVERAGE_ONLY=1 ;;
    -l|--listing-dir) LISTING_DIR="${2%/}"; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 1 ;;
  esac
  shift
done

[[ ${#ROOTS[@]} -gt 0 ]] || { echo "--run is required" >&2; exit 1; }
[[ -n "${PHENOTYPES}" ]] || { echo "--phenotypes is required" >&2; exit 1; }

case "${ARM}" in
  main)          CALL="Step2";   FNAME="results.felix.tsv" ;;
  ru)            CALL="RuStep2"; FNAME="results.ru.felix.tsv" ;;
  ru_biallelic)  CALL="BaStep2"; FNAME="results.ru.felix.tsv" ;;
  *) echo "--arm must be main, ru or ru_biallelic" >&2; exit 1 ;;
esac
OUT="${OUT:-interim_${ARM}}"

gs() {
  if [[ -n "${GOOGLE_PROJECT:-}" ]]; then gsutil -u "${GOOGLE_PROJECT}" "$@"; else gsutil "$@"; fi
}

mkdir -p "${OUT}"

# Phenotype names in scatter order; the inner shard index is the 0-based line.
if [[ "${PHENOTYPES}" == gs://* ]]; then
  gs cp "${PHENOTYPES}" "${OUT}/selected_phenotypes.txt" >/dev/null 2>&1
  PHENO_FILE="${OUT}/selected_phenotypes.txt"
else
  PHENO_FILE="${PHENOTYPES}"
fi
NPHENO="$(awk 'NF' "${PHENO_FILE}" | wc -l | tr -d ' ')"
echo "${NPHENO} phenotypes; arm ${ARM} (call-${CALL}/${FNAME})" >&2

echo "Listing ${FNAME} for ${CALL} in ${#ROOTS[@]} run root(s) ..." >&2
# List from the WORKFLOW ROOT, not from ${ROOT}/call-${CALL}/.
#
# Cromwell does not always put a call at the top level. An inner scatter can be
# compiled to an anonymous node, and then its calls live under
# call-ScatterAt<line>_<col>/shard-<outer>/call-<Name>/shard-<inner>/. In
# FelixGenome exactly the three per-phenotype step-2 scatters land that way --
# Step2, RuStep2 and BaStep2 -- while Pack, RuLoci, RuVcf and BaVcf stay at the
# top level. Listing ${ROOT}/call-${CALL}/ therefore matched nothing at all for
# precisely the three tasks this script reads.
#
# The call segment is matched as (call-)?<Name>/ anchored on a slash, so
# /RuStep2/ cannot satisfy a request for Step2, and RuStep2 and BaStep2 stay
# apart even though both write results.ru.felix.tsv.
#
# A full recursive listing is reused across arms when --listing-dir is given:
# one listing per root instead of one per root per arm.
: > "${OUT}/shards.txt"
RANK=0
for ROOT in "${ROOTS[@]}"; do
  RANK=$((RANK + 1))
  KEY="$(printf '%s' "${ROOT}" | shasum | cut -c1-12)"
  CACHE="${LISTING_DIR:-${OUT}}/listing-${KEY}.txt"
  mkdir -p "$(dirname "${CACHE}")"
  if [[ -s "${CACHE}" ]]; then
    echo "  ${ROOT}" >&2
    echo "    reusing $(wc -l < "${CACHE}" | tr -d ' ') listed objects from ${CACHE}" >&2
  else
    : > "${OUT}/ls.err"
    gs ls -r "${ROOT}/" 2>"${OUT}/ls.err" > "${CACHE}.part" || true
    mv "${CACHE}.part" "${CACHE}"
    echo "  ${ROOT}" >&2
    echo "    listed $(wc -l < "${CACHE}" | tr -d ' ') objects" >&2
    if [[ ! -s "${CACHE}" && -s "${OUT}/ls.err" ]]; then
      sed 's/^/      /' "${OUT}/ls.err" >&2
    fi
  fi
  awk -v fn="${FNAME}" -v call="${CALL}" -v r="${RANK}" -v root="${ROOT}/" '
    $0 ~ fn "$" && $0 ~ ("/(call-)?" call "/") {
      rel = $0; sub("^" root, "", rel); print r "\t" rel "\t" $0 }' "${CACHE}" \
    >> "${OUT}/shards.txt" || true
  echo "    $(awk -F'\t' -v r="${RANK}" '$1 == r' "${OUT}/shards.txt" | wc -l | tr -d ' ') ${FNAME} objects for ${CALL}" >&2
done
sort -o "${OUT}/shards.txt" "${OUT}/shards.txt"

if [[ ! -s "${OUT}/shards.txt" ]]; then
  echo "No ${FNAME} objects for ${CALL} in any run root given." >&2
  for ROOT in "${ROOTS[@]}"; do
    KEY="$(printf '%s' "${ROOT}" | shasum | cut -c1-12)"
    CACHE="${LISTING_DIR:-${OUT}}/listing-${KEY}.txt"
    echo "  call directories present under ${ROOT}:" >&2
    sed "s|^${ROOT}/||" "${CACHE}" 2>/dev/null \
      | grep -o '^call-[^/]*' | sort -u | head -20 | sed 's/^/    /' >&2
  done
  echo "  (if ${CALL} is not listed, the run has not reached it; the scan looks" >&2
  echo "   at any depth, so a nested ScatterAt node is already handled)" >&2
  exit 1
fi

# relpath -> inner shard index (the LAST shard-N before the filename) and a
# per-attempt dedupe: a retried shard has attempt-K, keep the highest.
awk -F'\t' '
  { rank = $1; rel = $2; uri = $3
    n = split(rel, p, "/")
    inner = -1; attempt = 1
    for (i = 1; i <= n; i++) {
      if (p[i] ~ /^shard-/) inner = substr(p[i], 7)
      else if (p[i] ~ /^attempt-/) attempt = substr(p[i], 9)
    }
    key = rel; sub("/attempt-[0-9]+/", "/", key)
    # Later root wins over an earlier one; within a root, the highest attempt.
    score = rank * 1000 + attempt
    if (score + 0 >= best[key] + 0) { best[key] = score; uriof[key] = uri; idxof[key] = inner }
  }
  END { for (k in uriof) print idxof[k] "\t" uriof[k] }' "${OUT}/shards.txt" \
  | sort -k1,1n -k2,2 > "${OUT}/by_pheno.tsv"

{
  printf 'phenotype\tindex\tchromosomes_present\tcomplete\n'
  awk -F'\t' -v nf="${PHENO_FILE}" -v minc="${MIN_CHROMS}" '
    BEGIN { i = 0; while ((getline line < nf) > 0) if (line ~ /[^ \t]/) name[i++] = line }
    { count[$1]++ }
    END { for (j = 0; j < i; j++)
            printf "%s\t%d\t%d\t%s\n", name[j], j, count[j] + 0,
                   (count[j] + 0 >= minc ? "yes" : "no") }' "${OUT}/by_pheno.tsv"
} > "${OUT}/coverage.tsv"

echo "" >&2
awk -F'\t' 'NR>1 {printf "  %-12s idx %-3s %2s/%s chromosomes  %s\n", $1, $2, $3, "'"${MIN_CHROMS}"'", ($4=="yes"?"COMPLETE":"")}' \
  "${OUT}/coverage.tsv" >&2
NCOMPLETE="$(awk -F'\t' 'NR>1 && $4=="yes"' "${OUT}/coverage.tsv" | wc -l | tr -d ' ')"
echo "" >&2
echo "${NCOMPLETE} of ${NPHENO} phenotypes have >= ${MIN_CHROMS} chromosomes." >&2
echo "Wrote ${OUT}/coverage.tsv" >&2

if [[ "${COVERAGE_ONLY}" -eq 1 ]]; then
  echo "--coverage-only: nothing merged." >&2
  exit 0
fi
if [[ "${NCOMPLETE}" -eq 0 ]]; then
  echo "Nothing to merge. Lower --min-chroms to merge partial phenotypes (label them as partial)." >&2
  exit 0
fi

# Merge: header from the first shard, then every shard's body. The shards of one
# phenotype are disjoint chromosomes, so order does not matter for lambdaGC or a
# hit table; Manhattan plotting sorts by CHROM/POS anyway.
awk -F'\t' 'NR>1 && $4=="yes" {print $2 "\t" $1}' "${OUT}/coverage.tsv" \
| while IFS=$'\t' read -r idx pheno; do
    dest="${OUT}/${pheno}.merged.felix.tsv.gz"
    if [[ -s "${dest}" ]]; then echo "  ${pheno}: exists, skipping" >&2; continue; fi
    mapfile -t uris < <(awk -F'\t' -v i="${idx}" '$1 == i {print $2}' "${OUT}/by_pheno.tsv") \
      2>/dev/null || uris=($(awk -F'\t' -v i="${idx}" '$1 == i {print $2}' "${OUT}/by_pheno.tsv"))
    echo "  ${pheno}: merging ${#uris[@]} shards" >&2
    first=1
    : > "${dest}.part"
    for u in ${uris[@]+"${uris[@]}"}; do
      if [[ "${first}" -eq 1 ]]; then
        gs cat "${u}" >> "${dest}.part"; first=0
      else
        gs cat "${u}" | tail -n +2 >> "${dest}.part"
      fi
    done
    mv "${dest}.part" "${dest}"
    echo "    $(wc -l < "${dest}" | tr -d ' ') lines" >&2
  done

echo "" >&2
echo "Merged TSVs in ${OUT}/. Summarise them with the same script the workflow uses:" >&2
echo "  python3 felix/scripts/summarize_felix_results.py \\" >&2
echo "    --results ${OUT}/*.merged.felix.tsv.gz \\" >&2
echo "    --phenotypes \$(awk -F'\\t' 'NR>1 && \$4==\"yes\" {print \$1}' ${OUT}/coverage.tsv) \\" >&2
echo "    --pheno-cov <pheno_cov.tsv> --out-dir ${OUT}/summary" >&2
echo "(summarize_felix_results.py takes no --arm; the WDL uses that only to name its bundle.)" >&2