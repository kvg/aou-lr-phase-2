#!/usr/bin/env bash
# Stage fixed AF-divergent LAI eval panels under $WORKSPACE_BUCKET/refs/flare/eval_panels/.
#
# Requires: bcftools, bgzip/tabix, WORKSPACE_BUCKET, local or gs:// paths to
# gnomAD LAI ref VCF + aou_1000genomes.refmap for the chromosomes of interest.
#
# Usage (Terra notebook / shell)::
#
#   export WORKSPACE_BUCKET=gs://...
#   export REF_VCF_CHR22=gs://.../chr22.gnomad_lai_90.vcf.bgz   # or local
#   export REF_PANEL=gs://.../aou_1000genomes.refmap
#   bash flare/scripts/stage_eval_af_panels.sh
#
set -euo pipefail

: "${WORKSPACE_BUCKET:?set WORKSPACE_BUCKET}"
: "${REF_PANEL:?set REF_PANEL to aou_1000genomes.refmap}"

SCRIPTS_DIR="${SCRIPTS_DIR:-$(cd "$(dirname "$0")/../../scripts" && pwd)}"
OUT_LOCAL="${OUT_LOCAL:-/tmp/flare_eval_panels}"
DEST="${WORKSPACE_BUCKET%/}/refs/flare/eval_panels"
mkdir -p "$OUT_LOCAL"

build_one() {
  local label="$1"
  local region="$2"
  local ref_vcf="$3"
  local prefix="$OUT_LOCAL/${label}"
  echo "=== building ${label} (${region}) ==="
  python3 "$SCRIPTS_DIR/flare_build_af_panel.py" \
    --ref-vcf "$ref_vcf" \
    --ref-panel "$REF_PANEL" \
    --region "$region" \
    --top-n "${TOP_N:-5000}" \
    --min-mac "${MIN_MAC:-50}" \
    --min-af-range "${MIN_AF_RANGE:-0.05}" \
    --out-prefix "$prefix"
}

# Chr22 method window (10 Mb)
if [[ -n "${REF_VCF_CHR22:-}" ]]; then
  build_one "chr22_10mb" "chr22:26897597-36897597" "$REF_VCF_CHR22"
else
  echo "skip chr22: set REF_VCF_CHR22" >&2
fi

# Full chr20 (or set CHR20_REGION to a narrower eval span)
if [[ -n "${REF_VCF_CHR20:-}" ]]; then
  build_one "chr20_full" "${CHR20_REGION:-chr20}" "$REF_VCF_CHR20"
else
  echo "skip chr20: set REF_VCF_CHR20" >&2
fi

echo "=== staging to ${DEST}/ ==="
gsutil -m rsync -r "$OUT_LOCAL/" "${DEST}/"
echo "staged panels under ${DEST}/"
ls -lh "$OUT_LOCAL" || true
