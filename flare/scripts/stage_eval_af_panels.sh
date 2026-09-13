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

# Resolve against this file's location (not the caller's cwd). Override with SCRIPTS_DIR
# if running from a tree that only has flare/ (e.g. point at the synced repo clone).
_HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -z "${SCRIPTS_DIR:-}" ]]; then
  if [[ -f "${_HERE}/../../scripts/flare_build_af_panel.py" ]]; then
    SCRIPTS_DIR="$(cd "${_HERE}/../../scripts" && pwd)"
  else
    echo "Cannot find repo scripts/ next to flare/scripts (${_HERE}/../../scripts)." >&2
    echo "Run from the synced git clone root, or set SCRIPTS_DIR, e.g.:" >&2
    echo "  export SCRIPTS_DIR=\$HOME/AoU_DRC_LongReads_PhaseTwo_Storage/edit/aou-lr-phase-2/scripts" >&2
    exit 1
  fi
fi
if [[ ! -f "${SCRIPTS_DIR}/flare_build_af_panel.py" ]]; then
  echo "SCRIPTS_DIR=${SCRIPTS_DIR} missing flare_build_af_panel.py" >&2
  exit 1
fi
OUT_LOCAL="${OUT_LOCAL:-/tmp/flare_eval_panels}"
CACHE_LOCAL="${CACHE_LOCAL:-/tmp/flare_eval_refs}"
DEST="${WORKSPACE_BUCKET%/}/refs/flare/eval_panels"
mkdir -p "$OUT_LOCAL" "$CACHE_LOCAL"

# pathlib collapses gs:// → gs:/; localize GCS inputs before calling Python/bcftools.
localize() {
  local uri="$1"
  local dest="$2"
  if [[ "$uri" == gs://* ]]; then
    if [[ ! -s "$dest" ]]; then
      echo "localize $uri -> $dest" >&2
      gsutil -q cp "$uri" "$dest"
    fi
    # sibling index when present (VCF)
    if [[ "$uri" == *.vcf.gz || "$uri" == *.vcf.bgz ]]; then
      local idx_dest="${dest}.tbi"
      if [[ ! -s "$idx_dest" ]]; then
        gsutil -q cp "${uri}.tbi" "$idx_dest" 2>/dev/null || true
      fi
    elif [[ "$uri" == *.bcf ]]; then
      local idx_dest="${dest}.csi"
      if [[ ! -s "$idx_dest" ]]; then
        gsutil -q cp "${uri}.csi" "$idx_dest" 2>/dev/null || true
      fi
    fi
    printf '%s' "$dest"
  else
    printf '%s' "$uri"
  fi
}

REF_PANEL_LOCAL="$(localize "$REF_PANEL" "$CACHE_LOCAL/$(basename "$REF_PANEL")")"

build_one() {
  local label="$1"
  local region="$2"
  local ref_vcf="$3"
  local prefix="$OUT_LOCAL/${label}"
  local vcf_local
  vcf_local="$(localize "$ref_vcf" "$CACHE_LOCAL/$(basename "$ref_vcf")")"
  echo "=== building ${label} (${region}) ==="
  python3 "$SCRIPTS_DIR/flare_build_af_panel.py" \
    --ref-vcf "$vcf_local" \
    --ref-panel "$REF_PANEL_LOCAL" \
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
