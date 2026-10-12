#!/usr/bin/env bash
# Upload FELIX WDL scripts to a Terra workspace bucket.
#
# Usage:
#   WORKSPACE_BUCKET=gs://fc-secure-... ./scripts/stage_felix_scripts.sh
#   ./scripts/stage_felix_scripts.sh gs://fc-secure-...
#
# Shared helpers (GRM build, repeat-unit annotator) stay in scripts/ — use
# ./scripts/stage_tractor_scripts.sh for those.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
FELIX_SCRIPTS="${REPO_ROOT}/felix/scripts"
BUCKET="${1:-${WORKSPACE_BUCKET:-}}"

if [[ -z "${BUCKET}" ]]; then
  echo "Set WORKSPACE_BUCKET or pass gs://... as the first argument." >&2
  exit 1
fi
BUCKET="${BUCKET%/}"
DEST="${BUCKET}/felix/scripts"

REQUIRED=(
  fit_felix_null.R
  export_felix_null.R
  run_felix_step2.R
  summarize_felix_results.py
  compare_felix_tractor_calibration.py
  compare_repeat_encodings.py
  compare_ru_arms.py
  simulate_repeat_dosage.py
)

echo "Uploading FELIX scripts to ${DEST}/"
for name in "${REQUIRED[@]}"; do
  src="${FELIX_SCRIPTS}/${name}"
  if [[ ! -f "${src}" ]]; then
    echo "Missing required script: ${src}" >&2
    exit 1
  fi
  gsutil cp "${src}" "${DEST}/${name}"
done

# Run-table configuration. felix_03_run_table.ipynb clones a row from the
# entity TSV, so the VM needs a copy; the inputs JSON goes up beside it so the
# method config can be re-uploaded from the bucket after a column change.
CONFIG_SRC="${REPO_ROOT}/felix/configs/terra"
CONFIG_DEST="${BUCKET}/felix/configs"
CONFIGS=(
  felix_runs.entity.tsv
  felix_runs.membership.tsv
  FelixGenomeTable.inputs.json
)

echo ""
echo "Uploading run-table configuration to ${CONFIG_DEST}/"
for name in "${CONFIGS[@]}"; do
  src="${CONFIG_SRC}/${name}"
  if [[ ! -f "${src}" ]]; then
    echo "Missing required config: ${src}" >&2
    exit 1
  fi
  gsutil cp "${src}" "${CONFIG_DEST}/${name}"
done

echo ""
echo "Verified objects:"
gsutil ls -l "${DEST}/fit_felix_null.R" "${DEST}/run_felix_step2.R" "${DEST}/summarize_felix_results.py"
echo ""
echo "FELIX docker: us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/felix-pilot:0.2.0"
echo ""
echo "FelixPilot.wdl / FelixGenome.wdl:"
echo "  fit_felix_null_script    -> ${DEST}/fit_felix_null.R"
echo "  export_null_script       -> ${DEST}/export_felix_null.R   (optional; overrides the copy baked into the image)"
echo "  run_felix_step2_script   -> ${DEST}/run_felix_step2.R"
echo "  compare_ru_arms_script   -> ${DEST}/compare_ru_arms.py   (FelixGenome ru_biallelic_arm 2x2 table)"
echo "  summarize_script         -> FelixGenome: ${DEST}/summarize_felix_results.py"
echo "                              FelixPilot: ${BUCKET}/scripts/summarize_tractor_genome_results.py (stage_tractor_scripts.sh)"
echo ""
echo "Shared scripts (stage separately):"
echo "  build_saige_grm_script   -> ${BUCKET}/scripts/build_saige_plink_and_grm.sh"
echo "  make_plink_keep_script   -> ${BUCKET}/scripts/make_plink_keep.py"
echo "  annotate_repeat_units    -> ${BUCKET}/scripts/annotate_repeat_units.py"
echo ""
echo "Eval gates: felix/eval/README.md"
