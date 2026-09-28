#!/usr/bin/env bash
# Cloud Batch / dsub: Nearline-efficient minicram for ExpansionHunter.
# CRAM must be a gs:// URI (--env), not a localized --input.
set -euo pipefail
date

: "${SAMPLE_ID:?}"
: "${CRAM:?}"
: "${CRAI:?}"
: "${REF_FA:?}"
: "${REF_FAI:?}"
: "${CATALOG:?}"
: "${MINICRAM:?}"
: "${MINICRAI:?}"
: "${TRANSFER_STATS:?}"

WINDOW_SIZE="${WINDOW_SIZE:-1000}"
MERGE_REGIONS_DISTANCE="${MERGE_REGIONS_DISTANCE:-1000}"
MAX_RETRY="${MAX_RETRY:-3}"
WAIT_TIME="${WAIT_TIME:-30}"

python3 -c "import json; json.load(open('${CATALOG}'))" \
  || { echo "catalog is not valid JSON: ${CATALOG}" >&2; exit 1; }

ln -sf "${REF_FA}" reference.fa
ln -sf "${REF_FAI}" reference.fa.fai

PROJ="${GCLOUD_PROJECT:-${GOOGLE_CLOUD_PROJECT:-}}"
if [[ -n "${PROJ}" ]]; then
  echo "requester-pays project: ${PROJ}"
fi

tmp_cram="${TMPDIR:-/tmp}/${SAMPLE_ID}.minicram.cram"

run_minicram() {
  local -a cmd=(
    python3 -u -m str_analysis.make_minicram_for_expansion_hunter
    -R reference.fa
    -c "${CATALOG}"
    -i "${CRAI}"
    -o "${tmp_cram}"
    -w "${WINDOW_SIZE}"
    -d "${MERGE_REGIONS_DISTANCE}"
    --verbose
    --output-data-transfer-stats
  )
  if [[ -n "${PROJ}" ]]; then
    cmd+=(--gcloud-project "${PROJ}")
  fi
  cmd+=("${CRAM}")
  "${cmd[@]}"
}

CURR=0
until (( CURR == MAX_RETRY )); do
  if run_minicram; then
    break
  fi
  echo "make_minicram_for_expansion_hunter failed; retry $((CURR + 1))/${MAX_RETRY}" >&2
  rm -f "${tmp_cram}" "${tmp_cram}.crai" *.data_transfer_stats.tsv
  CURR=$((CURR + 1))
  if (( CURR == MAX_RETRY )); then
    echo "make_minicram_for_expansion_hunter failed after ${MAX_RETRY} attempts" >&2
    exit 1
  fi
  sleep "${WAIT_TIME}"
done

if [[ ! -s "${tmp_cram}" || ! -s "${tmp_cram}.crai" ]]; then
  echo "make_minicram wrote no CRAM/CRAI" >&2
  exit 1
fi

mkdir -p "$(dirname "${MINICRAM}")" "$(dirname "${MINICRAI}")" "$(dirname "${TRANSFER_STATS}")"
cp -f "${tmp_cram}" "${MINICRAM}"
cp -f "${tmp_cram}.crai" "${MINICRAI}"
stats=$(ls -1 *.data_transfer_stats.tsv | head -n 1)
cp -f "${stats}" "${TRANSFER_STATS}"
ls -lh "${MINICRAM}" "${MINICRAI}" "${TRANSFER_STATS}"
date
