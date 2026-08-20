#!/usr/bin/env bash
# Upload Tractor-Mix / SAIGE WDL scripts to a Terra workspace bucket.
#
# Usage:
#   WORKSPACE_BUCKET=gs://fc-secure-... ./scripts/stage_tractor_scripts.sh
#   ./scripts/stage_tractor_scripts.sh gs://fc-secure-...
#
# Requires: gsutil, repo scripts/ directory.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SCRIPTS="${REPO_ROOT}/scripts"
BUCKET="${1:-${WORKSPACE_BUCKET:-}}"

if [[ -z "${BUCKET}" ]]; then
  echo "Set WORKSPACE_BUCKET or pass gs://... as the first argument." >&2
  exit 1
fi
BUCKET="${BUCKET%/}"
DEST="${BUCKET}/scripts"

REQUIRED=(
  fit_null.R
  fit_saige_null.R
  reorder_dosages.py
  sparsify_grm.R
  make_plink_keep.py
  run_saige_step2.R
  build_saige_plink_and_grm.sh
  compare_calibration.py
  compare_tractor_runs.py
  plot_tractor_results.py
  resolve_flare_uris.py
  select_phenotypes.py
  terra_notebook.py
  workspace_paths.py
)

echo "Uploading WDL scripts to ${DEST}/"
for name in "${REQUIRED[@]}"; do
  src="${SCRIPTS}/${name}"
  if [[ ! -f "${src}" ]]; then
    echo "Missing required script: ${src}" >&2
    exit 1
  fi
  gsutil cp "${src}" "${DEST}/${name}"
done

echo ""
echo "Verified objects:"
gsutil ls -l "${DEST}/fit_null.R" "${DEST}/sparsify_grm.R"
echo ""
echo "Tractor-Mix WDL (TractorMixPilot.wdl):"
echo "  fit_null_script -> ${DEST}/fit_null.R"
echo "  Score task uses tractor-mix-score --threads 8 (docker 0.4.2)"
echo "  Shared workflow docker: us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/tractor-mix-pilot:0.4.2"
echo ""
echo "Legacy fit_null_and_score.R is not staged; do not use it for new runs."
