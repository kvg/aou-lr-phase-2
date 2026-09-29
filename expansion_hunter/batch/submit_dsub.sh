#!/usr/bin/env bash
# Submit ExpansionHunter as two Cloud Batch jobs via dsub (VWB).
#
# Stage 1: make_minicram (print_reads image). CRAM is --env gs:// (not localized).
# Stage 2: ExpansionHunter on the minicram (--after stage 1).
#
# Usage (in a VWB Jupyter app; `dsub_activate` in the terminal):
#   ./submit_dsub.sh --csv ../../expansion_hunter/configs/batch.header.csv
#   ./submit_dsub.sh --csv smoke.csv --stage minicram
#   ./submit_dsub.sh --csv smoke.csv --stage genotype --after JOBID

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
# shellcheck source=../../scripts/vwb_dsub.sh
source "${REPO_ROOT}/scripts/vwb_dsub.sh"

CSV=""
STAGE="all"
AFTER_JOB=""
OUT_PREFIX="${EH_OUT_PREFIX:-gs://${GOOGLE_CLOUD_PROJECT:-aou-lr-phase2-resources}/batchRuns/expansion_hunter}"
# Prefer the workspace resources bucket when set by the notebook.
if [[ -n "${OUTPUT_BUCKET_GS:-}" ]]; then
  OUT_PREFIX="${OUTPUT_BUCKET_GS%/}/batchRuns/expansion_hunter"
fi
PRINT_READS_DOCKER="${PRINT_READS_DOCKER:-us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/aou-locityper-print-reads:0.1.1}"
EH_DOCKER="${EH_DOCKER:-us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/aou-expansion-hunter:0.1.0}"
WAIT=0
TASKS_DIR="${SCRIPT_DIR}/.tasks"

usage() {
  cat <<EOF
Submit ExpansionHunter minicram + genotype on Google Cloud Batch (dsub).

Options:
  --csv PATH           Per-sample CSV (see configs/batch.header.csv)
  --stage STAGE        minicram | genotype | all (default: all)
  --after JOB_ID       For --stage genotype: wait on this minicram job
  --out-prefix gs://   Output prefix (default: bucket/batchRuns/expansion_hunter)
  --wait               dsub --wait on the last submitted job
  -h, --help
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --csv) CSV="$2"; shift 2 ;;
    --stage) STAGE="$2"; shift 2 ;;
    --after) AFTER_JOB="$2"; shift 2 ;;
    --out-prefix) OUT_PREFIX="$2"; shift 2 ;;
    --wait) WAIT=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 1 ;;
  esac
done

if [[ -z "${CSV}" || ! -f "${CSV}" ]]; then
  echo "--csv is required and must exist" >&2
  exit 1
fi

vwb_dsub_require_env
mapfile -t BASE < <(vwb_dsub_base_args)
LOG_ROOT="${OUT_PREFIX%/}/logs"

python3 "${REPO_ROOT}/scripts/csv_to_dsub_tasks.py" \
  --pipeline expansion_hunter \
  --csv "${CSV}" \
  --out-prefix "${OUT_PREFIX}" \
  --outdir "${TASKS_DIR}"

# Set THIS_WAIT=1 before a call to pass dsub --wait (last job only).
THIS_WAIT=0
run_dsub() {
  local name="$1" image="$2" script="$3" tasks="$4" ram="$5" cores="$6" disk="$7" boot="$8"
  shift 8
  local extra=("$@")
  vwb_dsub_export_cloud_sdk_image
  local -a cmd=(
    dsub
    "${BASE[@]}"
    --name "${name}"
    --image "${image}"
    --script "${script}"
    --tasks "${tasks}"
    --logging "${LOG_ROOT}"
    --min-ram "${ram}"
    --min-cores "${cores}"
    --disk-size "${disk}"
    --boot-disk-size "${boot}"
  )
  if [[ ${#extra[@]} -gt 0 ]]; then
    cmd+=("${extra[@]}")
  fi
  # dsub: --retries is only legal together with --wait.
  if [[ "${THIS_WAIT}" -eq 1 ]]; then
    cmd+=(--wait --retries 2)
  fi
  printf '+' >&2
  printf ' %q' "${cmd[@]}" >&2
  printf '\n' >&2
  mkdir -p "${TASKS_DIR}"
  local log="${TASKS_DIR}/${name}.dsub.log"
  "${cmd[@]}" 2>&1 | tee "${log}"
  local id
  id=$(vwb_dsub_parse_job_id < "${log}")
  if [[ -z "${id}" ]]; then
    echo "failed to parse 'Launched job-id:' from ${log}" >&2
    exit 1
  fi
  printf '%s\n' "${id}" > "${TASKS_DIR}/${name}.job_id"
  echo "${name} job-id: ${id}"
  vwb_dsub_print_log_uris "${LOG_ROOT}" "${id}"
}

MINI_ID=""
if [[ "${STAGE}" == "minicram" || "${STAGE}" == "all" ]]; then
  THIS_WAIT=0
  [[ "${WAIT}" -eq 1 && "${STAGE}" == "minicram" ]] && THIS_WAIT=1
  run_dsub eh-minicram "${PRINT_READS_DOCKER}" \
    "${SCRIPT_DIR}/make_minicram.sh" "${TASKS_DIR}/minicram.tasks.tsv" \
    16 4 50 20 \
    --env "WINDOW_SIZE=${WINDOW_SIZE:-1000}" \
    --env "MERGE_REGIONS_DISTANCE=${MERGE_REGIONS_DISTANCE:-1000}" \
    --env "MAX_RETRY=${MAX_RETRY:-3}" \
    --env "WAIT_TIME=${WAIT_TIME:-30}"
  MINI_ID=$(<"${TASKS_DIR}/eh-minicram.job_id")
  echo "MINICRAM_JOB_ID=${MINI_ID}"
fi

if [[ "${STAGE}" == "genotype" || "${STAGE}" == "all" ]]; then
  after=()
  if [[ -n "${AFTER_JOB}" ]]; then
    after+=(--after "${AFTER_JOB}")
  elif [[ -n "${MINI_ID}" ]]; then
    after+=(--after "${MINI_ID}")
  else
    echo "genotype needs --after JOB_ID (or --stage all)" >&2
    exit 1
  fi
  THIS_WAIT=0
  [[ "${WAIT}" -eq 1 ]] && THIS_WAIT=1
  # --after blocks until minicram finishes, then submits genotype.
  run_dsub eh-genotype "${EH_DOCKER}" \
    "${SCRIPT_DIR}/genotype.sh" "${TASKS_DIR}/genotype.tasks.tsv" \
    8 4 50 30 \
    --preemptible 2 \
    "${after[@]}"
  echo "GENOTYPE_JOB_ID=$(<"${TASKS_DIR}/eh-genotype.job_id")"
fi
