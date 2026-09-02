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
  fit_felix_null.R
  export_felix_null.R
  reorder_dosages.py
  sparsify_grm.R
  make_plink_keep.py
  run_saige_step2.R
  run_felix_step2.R
  build_saige_plink_and_grm.sh
  summarize_tractor_genome_results.py
  summarize_felix_results.py
  annotate_repeat_units.py
  compare_calibration.py
  compare_felix_tractor_calibration.py
  compare_repeat_encodings.py
  simulate_repeat_dosage.py
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
gsutil ls -l "${DEST}/fit_null.R" "${DEST}/fit_saige_null.R" "${DEST}/build_saige_plink_and_grm.sh"
echo ""
echo "Shared workflow docker: us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/tractor-mix-pilot:0.4.2"
echo ""
echo "Tractor-Mix WDL (TractorMixPilot.wdl):"
echo "  fit_null_script          -> ${DEST}/fit_null.R"
echo "  make_plink_keep_script   -> ${DEST}/make_plink_keep.py"
echo "  Score uses tractor-mix-score --threads 8"
echo ""
echo "SAIGE WDL (SaigePilot.wdl):"
echo "  build_saige_grm_script   -> ${DEST}/build_saige_plink_and_grm.sh"
echo "  make_plink_keep_script   -> ${DEST}/make_plink_keep.py"
echo "  fit_saige_null_script    -> ${DEST}/fit_saige_null.R"
echo "  run_saige_step2_script   -> ${DEST}/run_saige_step2.R"
echo ""
echo "Genome-wide Tractor-Mix (TractorMixGenome.wdl):"
echo "  same scripts as TractorMixPilot; supply chroms + flare_vcfs arrays"
echo "  resolve URIs: python3 scripts/resolve_flare_uris.py --from-firecloud --autosomes"
echo "  summarize_script          -> ${DEST}/summarize_tractor_genome_results.py"
echo ""
echo "FELIX WDL (FelixPilot.wdl):"
echo "  image                    -> us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/felix-pilot:0.1.0"
echo "  build_saige_grm_script   -> ${DEST}/build_saige_plink_and_grm.sh"
echo "  make_plink_keep_script   -> ${DEST}/make_plink_keep.py"
echo "  fit_felix_null_script    -> ${DEST}/fit_felix_null.R"
echo "  export_felix_null_script -> ${DEST}/export_felix_null.R"
echo "  run_felix_step2_script   -> ${DEST}/run_felix_step2.R"
echo "  summarize_script         -> ${DEST}/summarize_felix_results.py"
echo ""
echo "Genome-wide FELIX (FelixGenome.wdl):"
echo "  same scripts as FelixPilot; supply chroms + phase_vcfs + flare_vcfs (or joint_vcfs)"
echo "  optional RU_TEST: simple_repeat_bed + annotate_repeat_units_script"
echo "  resolve URIs: python3 scripts/resolve_flare_uris.py --from-firecloud --autosomes"
echo "  eval gates: tractor_mix/eval/README.md"
echo ""
echo "Legacy fit_null_and_score.R is not staged; do not use it for new runs."
