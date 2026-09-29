#!/usr/bin/env bash
# Native Google Batch minicram / genotype (no dsub).
#
# In a VWB Jupyter terminal:
#   ./submit_batch.sh --csv ../configs/batch.header.csv --stage minicram
#   ./submit_batch.sh --csv ../configs/batch.header.csv --stage genotype

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
EH_DOCKER="${EH_DOCKER:-us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/aou-expansion-hunter:0.1.0}"

usage() {
  cat <<EOF
Submit ExpansionHunter on native Google Batch (no dsub).

Options:
  --csv PATH           Per-sample CSV (see configs/batch.header.csv)
  --stage STAGE        minicram | genotype (default: minicram)
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
if [[ "${STAGE}" != "minicram" && "${STAGE}" != "genotype" ]]; then
  echo "--stage must be minicram or genotype" >&2
  exit 1
fi

vwb_dsub_require_env
REGION="${DSUB_REGION:-us-central1}"
if [[ "${STAGE}" == "minicram" ]]; then
  JOB_ID="eh-mini-$(date -u +%y%m%d-%H%M%S)"
  IMAGE="${PRINT_READS_DOCKER}"
  WORKER="${SCRIPT_DIR}/make_minicram.py"
else
  JOB_ID="eh-gt-$(date -u +%y%m%d-%H%M%S)"
  IMAGE="${EH_DOCKER}"
  WORKER="${SCRIPT_DIR}/genotype.py"
fi
CFG="$(mktemp)"
trap 'rm -f "${CFG}"' EXIT

python3 "${SCRIPT_DIR}/submit_batch.py" \
  --csv "${CSV}" \
  --worker "${WORKER}" \
  --image "${IMAGE}" \
  --out-prefix "${OUT_PREFIX}" \
  --project "${GOOGLE_CLOUD_PROJECT}" \
  --region "${REGION}" \
  --sa "${PET_SA_EMAIL}" \
  --stage "${STAGE}" \
  --config "${CFG}"

echo "JOB_ID=${JOB_ID}"
echo "STAGE=${STAGE}"
echo "IMAGE=${IMAGE}"
echo "OUT_PREFIX=${OUT_PREFIX}"

gcloud batch jobs submit "${JOB_ID}" \
  --location="${REGION}" \
  --project="${GOOGLE_CLOUD_PROJECT}" \
  --config="${CFG}"

echo "Watch:"
echo "  gcloud batch jobs describe ${JOB_ID} --location=${REGION} --format='yaml(status.state,status.statusEvents)'"
echo "When SUCCEEDED:"
echo "  gcloud storage ls ${OUT_PREFIX%/}/SAMPLE_ID/"
