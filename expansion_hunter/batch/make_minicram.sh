#!/usr/bin/env bash
# Cloud Batch / dsub: Nearline-efficient minicram for ExpansionHunter.
# CRAM and CRAI must be gs:// URIs (--env), not localized --input.
# The CRAI lives in the same requester-pays bucket; localizing it is what
# failed MakeMinicram under both Cromwell and dsub.
set -euo pipefail

# dsub --logging never lands in GCS here (PET cannot read Cloud Logging).
# Re-exec once, capturing stdout/stderr, then upload that file ourselves.
upload_worker_log() {
  local local_log="$1"
  [[ -n "${WORKER_LOG_GCS:-}" && -f "${local_log}" ]] || return 0
  python3 - "${local_log}" "${WORKER_LOG_GCS}" <<'PY'
import sys
from google.cloud import storage

local, dest = sys.argv[1], sys.argv[2]
if not dest.startswith("gs://"):
    raise SystemExit(f"WORKER_LOG_GCS is not gs://: {dest}")
bucket_name, blob_name = dest[5:].split("/", 1)
storage.Client().bucket(bucket_name).blob(blob_name).upload_from_filename(local)
print(f"uploaded worker log to {dest}", flush=True)
PY
}

if [[ -z "${WORKER_LOG_INNER:-}" ]]; then
  export WORKER_LOG_INNER=1
  status=0
  bash "$0" "$@" >"${PWD}/worker.log" 2>&1 || status=$?
  upload_worker_log "${PWD}/worker.log" || true
  exit "${status}"
fi

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

echo "SAMPLE_ID=${SAMPLE_ID}"
echo "CRAM=${CRAM}"
echo "CRAI=${CRAI}"
echo "REF_FA=${REF_FA}"
echo "CATALOG=${CATALOG}"
python3 -c "import str_analysis; print('str_analysis', str_analysis.__file__)"

python3 -c "import json; json.load(open('${CATALOG}'))" \
  || { echo "catalog is not valid JSON: ${CATALOG}" >&2; exit 1; }

ln -sf "${REF_FA}" reference.fa
ln -sf "${REF_FAI}" reference.fa.fai

PROJ="${GCLOUD_PROJECT:-${GOOGLE_CLOUD_PROJECT:-}}"
if [[ -z "${PROJ}" ]]; then
  echo "GCLOUD_PROJECT / GOOGLE_CLOUD_PROJECT empty; requester-pays CRAM reads will 403" >&2
  exit 1
fi
echo "requester-pays project: ${PROJ}"

if [[ "${CRAM}" != gs://* || "${CRAI}" != gs://* ]]; then
  echo "CRAM and CRAI must be gs:// URIs (not localized files): CRAM=${CRAM} CRAI=${CRAI}" >&2
  exit 1
fi

# Write on the dsub data disk, not /tmp (boot disk).
tmp_cram="${PWD}/${SAMPLE_ID}.minicram.cram"

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
    --gcloud-project "${PROJ}"
    "${CRAM}"
  )
  printf '+'
  printf ' %q' "${cmd[@]}"
  printf '\n'
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
