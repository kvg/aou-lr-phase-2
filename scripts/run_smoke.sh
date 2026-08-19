#!/usr/bin/env bash
# Local smoke: unit tests + optional miniwdl/womtool syntax check.
set -euo pipefail
SCRIPTS_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "${SCRIPTS_DIR}/.." && pwd)"
PKG="${REPO}/sv_annotation"
cd "$PKG"

if command -v python3 >/dev/null; then
  PY=python3
fi
if [[ -x "${REPO}/venv/bin/python" ]]; then
  PY="${REPO}/venv/bin/python"
fi

"$PY" -m pytest tests/ -q

MINIWDL="$(dirname "$PY")/miniwdl"
if [[ -x "$MINIWDL" ]]; then
  "$MINIWDL" check wdl/AnnotateSvCallset.wdl
  echo "miniwdl check OK"
else
  echo "miniwdl not installed; skipped WDL check"
fi

if [[ -n "${SV_ANNOTATION_IMAGE:-}" ]]; then
  docker run --rm --platform linux/amd64 "$SV_ANNOTATION_IMAGE" bash -lc \
    'caddsv --help >/dev/null && bcftools --version >/dev/null && python /opt/aou_sv/scripts/run_caddsv.py --help >/dev/null'
  echo "Docker image smoke OK"
fi
