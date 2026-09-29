#!/usr/bin/env bash
# Cloud Batch / dsub: Nearline-efficient minicram for Locityper (print_reads).
# CRAM and CRAI must be gs:// URIs (--env), not localized --input.
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
: "${BED:?}"
: "${MINICRAM:?}"
: "${MINICRAI:?}"
: "${TRANSFER_STATS:?}"

WINDOW_GRAB="${WINDOW_GRAB:-3000}"
# GRCh38 Locityper background interval (Isaac). Empty skips the extra grab.
BG_REGION_BED="${BG_REGION_BED:-chr17	72062001	76562000}"
MAX_RETRY="${MAX_RETRY:-3}"
WAIT_TIME="${WAIT_TIME:-30}"

echo "SAMPLE_ID=${SAMPLE_ID}"
echo "CRAM=${CRAM}"
echo "CRAI=${CRAI}"
echo "REF_FA=${REF_FA}"
echo "BED=${BED}"
python3 -c "import str_analysis; print('str_analysis', str_analysis.__file__)"

awk -v w="${WINDOW_GRAB}" 'BEGIN { OFS="\t" } {
  s = $2 - w
  if (s < 0) s = 0
  print $1, s, $3 + w
}' "${BED}" > intervals.bed

if [[ -n "${BG_REGION_BED}" ]]; then
  echo -e "${BG_REGION_BED}" | awk -v w="${WINDOW_GRAB}" 'BEGIN { OFS="\t" } {
    s = $2 - w
    if (s < 0) s = 0
    print $1, s, $3 + w
  }' >> intervals.bed
fi

echo "minicram intervals:"
cat intervals.bed
wc -l intervals.bed

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

tmp_cram="${PWD}/${SAMPLE_ID}.minicram.cram"

run_print_reads() {
  local -a cmd=(
    python3 -m str_analysis.print_reads
    -R "${REF_FA}"
    --read-index "${CRAI}"
    -L intervals.bed
    --padding 0
    -o "${tmp_cram}"
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
  if run_print_reads; then
    break
  fi
  echo "print_reads failed; retry $((CURR + 1))/${MAX_RETRY}" >&2
  rm -f "${tmp_cram}" "${tmp_cram}.crai" *.data_transfer_stats.tsv
  CURR=$((CURR + 1))
  if (( CURR == MAX_RETRY )); then
    echo "print_reads failed after ${MAX_RETRY} attempts" >&2
    exit 1
  fi
  sleep "${WAIT_TIME}"
done

if [[ ! -s "${tmp_cram}" || ! -s "${tmp_cram}.crai" ]]; then
  echo "print_reads wrote no CRAM/CRAI" >&2
  exit 1
fi

mkdir -p "$(dirname "${MINICRAM}")" "$(dirname "${MINICRAI}")" "$(dirname "${TRANSFER_STATS}")"
cp -f "${tmp_cram}" "${MINICRAM}"
cp -f "${tmp_cram}.crai" "${MINICRAI}"
stats=$(ls -1 *.data_transfer_stats.tsv | head -n 1)
cp -f "${stats}" "${TRANSFER_STATS}"
ls -lh "${MINICRAM}" "${MINICRAI}" "${TRANSFER_STATS}"
date
