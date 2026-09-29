#!/usr/bin/env bash
# Cloud Batch / dsub: Nearline-efficient minicram for ExpansionHunter.
# CRAM and CRAI must be gs:// URIs (--env), not localized --input.
set -euo pipefail

# PET cannot read Cloud Logging; dsub --logging therefore never appears in GCS.
# Capture the worker transcript and upload it ourselves (success or failure).
# Do not re-exec $0 — dsub does not launch --script that way.
upload_worker_log() {
  local local_log="$1"
  if [[ -z "${WORKER_LOG_GCS:-}" ]]; then
    echo "WORKER_LOG_GCS is unset; cannot upload ${local_log}" >&2
    return 1
  fi
  if [[ ! -f "${local_log}" ]]; then
    echo "no local worker log at ${local_log}" >&2
    return 1
  fi
  # Batch copies --input/--output from the agent. ADC inside the user
  # container is often missing; use the GCE metadata token instead.
  python3 - "${local_log}" "${WORKER_LOG_GCS}" <<'PY'
import json, sys, urllib.error, urllib.parse, urllib.request

local, dest = sys.argv[1], sys.argv[2]
if not dest.startswith("gs://"):
    raise SystemExit(f"WORKER_LOG_GCS is not gs://: {dest}")
bucket, _, blob = dest[5:].partition("/")
blob_q = urllib.parse.quote(blob, safe="")

def metadata(path: str) -> bytes:
    last = None
    for host in ("metadata.google.internal", "169.254.169.254"):
        req = urllib.request.Request(
            f"http://{host}/computeMetadata/v1/{path}",
            headers={"Metadata-Flavor": "Google"},
        )
        try:
            return urllib.request.urlopen(req, timeout=10).read()
        except Exception as exc:
            last = exc
    raise RuntimeError(f"metadata token failed: {last}")

token = json.loads(metadata("instance/service-accounts/default/token"))["access_token"]
url = (
    "https://storage.googleapis.com/upload/storage/v1/b/"
    f"{urllib.parse.quote(bucket, safe='')}/o?uploadType=media&name={blob_q}"
)
body = open(local, "rb").read()
req = urllib.request.Request(url, data=body, method="POST")
req.add_header("Authorization", f"Bearer {token}")
req.add_header("Content-Type", "text/plain")
try:
    urllib.request.urlopen(req, timeout=60)
except urllib.error.HTTPError as exc:
    sys.stderr.write(exc.read().decode("utf-8", "replace")[:4000] + "\n")
    raise
print(f"uploaded worker log to {dest} ({len(body)} bytes)", flush=True)
PY
}

LOG_FILE="${PWD}/worker.log"
echo "minicram start $(date -u) SAMPLE_ID=${SAMPLE_ID:-?} WORKER_LOG_GCS=${WORKER_LOG_GCS:-UNSET}" > "${LOG_FILE}"
upload_worker_log "${LOG_FILE}" || echo "heartbeat log upload failed" >&2

run() {
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
  python3 -c "import pkgutil; print('has make_minicram', bool(pkgutil.find_loader('str_analysis.make_minicram_for_expansion_hunter')))"

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
}

status=0
run >>"${LOG_FILE}" 2>&1 || status=$?
upload_worker_log "${LOG_FILE}" || echo "final log upload failed" >&2
exit "${status}"
