#!/usr/bin/env bash
set -euo pipefail

# Run R-oracle vs tractor-mix-score parity inside the pilot image (GMMAT + binaries).
# Does NOT require cargo in the image.
#
# Usage:
#   ./scripts/run_oracle_parity_in_docker.sh
#
# Optional:
#   IMAGE=.../tractor-mix-pilot:0.4.2 ./scripts/run_oracle_parity_in_docker.sh
#   RUST_BIN=/path/to/local/tractor-mix-score ./scripts/run_oracle_parity_in_docker.sh

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CRATE_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
IMAGE="${IMAGE:-us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/tractor-mix-pilot:0.4.2}"
PLATFORM="${PLATFORM:-linux/amd64}"
AC_THRESHOLD="${AC_THRESHOLD:-10}"
CHUNK_SIZE="${CHUNK_SIZE:-16}"

echo "Running oracle parity in Docker image: ${IMAGE} (platform=${PLATFORM})"

docker run --rm --platform "${PLATFORM}" \
  -v "${CRATE_DIR}:/work" \
  -w /work \
  -e AC_THRESHOLD="${AC_THRESHOLD}" \
  -e CHUNK_SIZE="${CHUNK_SIZE}" \
  -e RUST_BIN="${RUST_BIN:-/usr/local/bin/tractor-mix-score}" \
  "${IMAGE}" \
  bash -lc '
set -euo pipefail
WORKDIR=/tmp/oracle_parity
FIX="${WORKDIR}/fixture"
mkdir -p "${WORKDIR}"

echo "=== GMMAT check ==="
Rscript -e "suppressPackageStartupMessages(library(GMMAT)); cat(\"GMMAT OK\n\")"

echo "=== Generate synthetic fixture ==="
Rscript testdata/gen_oracle.R "${FIX}"

DOSAGES=()
for i in 00 01 02 03 04; do
  DOSAGES+=("${FIX}/dosages/anc_${i}.dosage.txt.gz")
done

echo "=== R oracle (pinned TractorMix.score.R) ==="
Rscript testdata/run_oracle.R \
  --null-rds "${FIX}/null.rds" \
  --out "${WORKDIR}/r.tsv" \
  --ac-threshold "${AC_THRESHOLD}" \
  --n-core 1 \
  --chunk-size "${CHUNK_SIZE}" \
  --oracle-r /work/oracle/TractorMix.score.R \
  --dosage-files "${DOSAGES[@]}"

echo "=== Rust scorer (${RUST_BIN}) ==="
"${RUST_BIN}" --help >/dev/null
"${RUST_BIN}" \
  --null-export "${FIX}/null_export" \
  --dosage-files "${DOSAGES[@]}" \
  --out "${WORKDIR}/rust.tsv" \
  --ac-threshold "${AC_THRESHOLD}" \
  --chunk-size "${CHUNK_SIZE}" \
  --threads 1

echo "=== Compare TSVs (round/signif parity rules) ==="
python3 testdata/compare_oracle_tsv.py "${WORKDIR}/r.tsv" "${WORKDIR}/rust.tsv"
'

echo "Oracle parity passed."
