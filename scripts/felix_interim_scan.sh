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
#   <phenotype>.merged.felix.tsv   one per complete phenotype
#
# Usage:
#   ./scripts/felix_interim_scan.sh --run gs://.../FelixGenome/<wid> \
#       --phenotypes gs://BUCKET/.../selected_phenotypes.txt --arm main --out interim_main
#   ... --arm ru          # repeat-dosage arm (RuStep2)
#   ... --arm ru_biallelic # carrier-coded comparator (BaStep2)
#   ... --min-chroms 22   # merge only phenotypes with this many chromosomes (default 22)
#   ... --coverage-only   # report coverage, merge nothing
#
# Set GOOGLE_PROJECT to bill a requester-pays bucket (adds gsutil -u).
#
# Requires: bash, gsutil, awk.

set -euo pipefail

RUN_ROOT=""
PHENOTYPES=""
ARM="main"
OUT=""
MIN_CHROMS=22
COVERAGE_ONLY=0

usage() { sed -n '2,31p' "$0" | sed 's/^# \?//'; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    -r|--run) RUN_ROOT="$2"; shift ;;
    -p|--phenotypes) PHENOTYPES="$2"; shift ;;
    -a|--arm) ARM="$2"; shift ;;
    -o|--out) OUT="$2"; shift ;;
    -m|--min-chroms) MIN_CHROMS="$2"; shift ;;
    --coverage-only) COVERAGE_ONLY=1 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 1 ;;
  esac
  shift
done

[[ -n "${RUN_ROOT}" ]]   || { echo "--run is required" >&2; exit 1; }
[[ -n "${PHENOTYPES}" ]] || { echo "--phenotypes is required" >&2; exit 1; }
RUN_ROOT="${RUN_ROOT%/}"

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

echo "Listing ${FNAME} under call-${CALL} ..." >&2
# Filter on the call segment as well as the filename: RuStep2 and BaStep2 write
# the same results.ru.felix.tsv, so a listing that came back wider than asked
# for would silently merge both arms together.
gs ls "${RUN_ROOT}/call-${CALL}/**/${FNAME}" 2>/dev/null \
  | awk -v fn="${FNAME}" -v call="/call-${CALL}/" '$0 ~ fn "$" && index($0, call)' \
  | sort > "${OUT}/shards.txt" || true

if [[ ! -s "${OUT}/shards.txt" ]]; then
  echo "No ${FNAME} objects under ${RUN_ROOT}/call-${CALL}/" >&2
  exit 1
fi

# relpath -> inner shard index (the LAST shard-N before the filename) and a
# per-attempt dedupe: a retried shard has attempt-K, keep the highest.
awk -v root="${RUN_ROOT}/" '
  { uri = $0; rel = uri; sub("^" root, "", rel)
    n = split(rel, p, "/")
    inner = -1; attempt = 1
    for (i = 1; i <= n; i++) {
      if (p[i] ~ /^shard-/) inner = substr(p[i], 7)
      else if (p[i] ~ /^attempt-/) attempt = substr(p[i], 9)
    }
    key = rel; sub("/attempt-[0-9]+/", "/", key)
    if (attempt + 0 >= best[key] + 0) { best[key] = attempt; uriof[key] = uri; idxof[key] = inner }
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
    dest="${OUT}/${pheno}.merged.felix.tsv"
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
echo "    --results ${OUT}/*.merged.felix.tsv \\" >&2
echo "    --phenotypes \$(awk -F'\\t' 'NR>1 && \$4==\"yes\" {print \$1}' ${OUT}/coverage.tsv) \\" >&2
echo "    --pheno-cov <pheno_cov.tsv> --out-dir ${OUT}/summary" >&2
echo "(summarize_felix_results.py takes no --arm; the WDL uses that only to name its bundle.)" >&2