#!/usr/bin/env bash
# Native Google Batch minicram (no dsub). Same PET SA + private VPC as the
# successful eh-ctr-smoke python GCS copy.
#
# In a VWB Jupyter terminal:
#   ./submit_batch.sh --csv ../configs/batch.header.csv --stage minicram

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
# shellcheck source=../../scripts/vwb_dsub.sh
source "${REPO_ROOT}/scripts/vwb_dsub.sh"

CSV=""
STAGE="minicram"
OUT_PREFIX="${EH_OUT_PREFIX:-gs://${GOOGLE_CLOUD_PROJECT:-aou-lr-phase2-resources}/batchRuns/expansion_hunter}"
if [[ -n "${OUTPUT_BUCKET_GS:-}" ]]; then
  OUT_PREFIX="${OUTPUT_BUCKET_GS%/}/batchRuns/expansion_hunter"
fi
PRINT_READS_DOCKER="${PRINT_READS_DOCKER:-us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/aou-locityper-print-reads:0.1.2}"

usage() {
  cat <<EOF
Submit ExpansionHunter minicram on native Google Batch (no dsub).

Options:
  --csv PATH           Per-sample CSV (see configs/batch.header.csv)
  --stage STAGE        minicram (default). genotype is not wired yet.
  --out-prefix gs://   Output prefix
  -h, --help
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --csv) CSV="$2"; shift 2 ;;
    --stage) STAGE="$2"; shift 2 ;;
    --out-prefix) OUT_PREFIX="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 1 ;;
  esac
done

if [[ -z "${CSV}" || ! -f "${CSV}" ]]; then
  echo "--csv is required and must exist" >&2
  exit 1
fi
if [[ "${STAGE}" != "minicram" ]]; then
  echo "only --stage minicram is implemented (drop dsub first, then genotype)" >&2
  exit 1
fi

vwb_dsub_require_env
REGION="${DSUB_REGION:-us-central1}"
JOB_ID="eh-mini-$(date -u +%y%m%d-%H%M%S)"
CFG="$(mktemp)"
trap 'rm -f "${CFG}"' EXIT

python3 "${SCRIPT_DIR}/submit_batch.py" \
  --csv "${CSV}" \
  --worker "${SCRIPT_DIR}/make_minicram.py" \
  --image "${PRINT_READS_DOCKER}" \
  --out-prefix "${OUT_PREFIX}" \
  --project "${GOOGLE_CLOUD_PROJECT}" \
  --region "${REGION}" \
  --sa "${PET_SA_EMAIL}" \
  --config "${CFG}"

echo "JOB_ID=${JOB_ID}"
echo "IMAGE=${PRINT_READS_DOCKER}"
echo "OUT_PREFIX=${OUT_PREFIX}"

gcloud batch jobs submit "${JOB_ID}" \
  --location="${REGION}" \
  --project="${GOOGLE_CLOUD_PROJECT}" \
  --config="${CFG}"

echo "Watch:"
echo "  gcloud batch jobs describe ${JOB_ID} --location=${REGION} --format='yaml(status.state,status.statusEvents)'"
echo "When SUCCEEDED, for each sample:"
echo "  gcloud storage ls ${OUT_PREFIX%/}/SAMPLE_ID/"
