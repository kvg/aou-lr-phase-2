#!/usr/bin/env bash
# Perimeter-safe version of the VWB dsub example:
# https://support.workbench.verily.com/docs/guides/workflows/dsub/
#
# Same Batch flags (PET SA, private VPC, --use-private-address, --wait).
# The doc's quay.io samtools image and gs://genomics-public-data BAM are
# outside AoU VPC-SC, so this copies a tiny object already in the workspace
# bucket using the print-reads image (gcloud + restricted-VIP DNS).
#
# In a VWB Jupyter terminal:
#   dsub_activate
#   ./scripts/vwb_dsub_smoke.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
# shellcheck source=vwb_dsub.sh
source "${REPO_ROOT}/scripts/vwb_dsub.sh"

vwb_dsub_require_env
vwb_dsub_export_cloud_sdk_image
mapfile -t BASE < <(vwb_dsub_base_args)

BUCKET="${VWB_DSUB_SMOKE_BUCKET:-gs://aou-lr-phase2-resources}"
WORKDIR="${BUCKET%/}/dsub_example"
IMAGE="${DSUB_CLOUD_SDK_IMAGE}"
IN_URI="${VWB_DSUB_SMOKE_IN:-${BUCKET%/}/expansion_hunter/smoke.catalog.json}"

echo "WORKDIR=${WORKDIR}"
echo "IN_URI=${IN_URI}"
echo "IMAGE=${IMAGE}"

dsub \
  "${BASE[@]}" \
  --name vwb-dsub-smoke \
  --logging "${WORKDIR}/logs" \
  --image "${IMAGE}" \
  --input IN="${IN_URI}" \
  --output OUT="${WORKDIR}/smoke.catalog.copy.json" \
  --command 'cp "${IN}" "${OUT}"' \
  --wait

echo "SUCCESS. Objects:"
gcloud storage ls "${WORKDIR}/**"
echo "Logs:"
gcloud storage ls "${WORKDIR}/logs/**" || true
