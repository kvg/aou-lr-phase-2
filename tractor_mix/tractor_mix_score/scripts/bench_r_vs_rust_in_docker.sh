#!/usr/bin/env bash
set -euo pipefail

# Benchmark R TractorMix.score vs tractor-mix-score inside the pilot Docker image.
# Reports wall-clock time and parity summary (same TSV formatting rules as R).
#
# Usage:
#   ./scripts/bench_r_vs_rust_in_docker.sh
#
# Optional env:
#   IMAGE=.../tractor-mix-pilot:0.4.2
#   N_SAMPLES=800 N_SITES=2000 N_ANC=5
#   AC_THRESHOLD=50 R_CORES=4 RUST_THREADS=8
#   SKIP_SMALL_PARITY=1          # skip 48-row byte-identical gate
#   RUST_BIN=/usr/local/bin/tractor-mix-score

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CRATE_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
IMAGE="${IMAGE:-us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/tractor-mix-pilot:0.4.2}"
R_CORES="${R_CORES:-4}"
PLATFORM="${PLATFORM:-linux/amd64}"
N_SAMPLES="${N_SAMPLES:-800}"
N_SITES="${N_SITES:-2000}"
N_ANC="${N_ANC:-5}"
AC_THRESHOLD="${AC_THRESHOLD:-50}"
RUST_THREADS="${RUST_THREADS:-8}"
R_CHUNK="${R_CHUNK:-256}"
RUST_CHUNK="${RUST_CHUNK:-2048}"
SKIP_SMALL_PARITY="${SKIP_SMALL_PARITY:-0}"

echo "Benchmark R vs Rust in Docker image: ${IMAGE} (platform=${PLATFORM})"
echo "Scale: n=${N_SAMPLES}, sites=${N_SITES}, ancestries=${N_ANC}, ac_threshold=${AC_THRESHOLD}"

docker run --rm --platform "${PLATFORM}" --cpus="${DOCKER_CPUS:-8}" \
  -v "${CRATE_DIR}:/work" \
  -w /work \
  -e SKIP_SMALL_PARITY="${SKIP_SMALL_PARITY}" \
  -e N_SAMPLES="${N_SAMPLES}" \
  -e N_SITES="${N_SITES}" \
  -e N_ANC="${N_ANC}" \
  -e AC_THRESHOLD="${AC_THRESHOLD}" \
  -e R_CORES="${R_CORES}" \
  -e RUST_THREADS="${RUST_THREADS}" \
  -e R_CHUNK="${R_CHUNK}" \
  -e RUST_CHUNK="${RUST_CHUNK}" \
  -e RUST_BIN="${RUST_BIN:-/usr/local/bin/tractor-mix-score}" \
  "${IMAGE}" \
  bash -lc '
set -euo pipefail
WORKDIR=/tmp/r_vs_rust_bench
mkdir -p "${WORKDIR}"

time_score() {
  local label="$1"
  shift
  local timing_file="/tmp/time.$$.${RANDOM}"
  set +e
  /usr/bin/time -f "%e" -o "${timing_file}" "$@" >/dev/null
  local status=$?
  set -e
  local sec
  sec=$(grep -E "^[0-9]+(\\.[0-9]+)?\$" "${timing_file}" || true)
  rm -f "${timing_file}"
  if [[ ${status} -ne 0 || -z "${sec}" ]]; then
    echo "${label}: command failed (exit ${status})" >&2
    exit 1
  fi
  echo "${label}=${sec}"
}

if [[ "${SKIP_SMALL_PARITY}" != "1" ]]; then
  echo ""
  echo "=== Small parity gate (48 variants, byte-identical expected) ==="
  FIX="${WORKDIR}/small"
  Rscript testdata/gen_oracle.R "${FIX}"
  SDOS=()
  for i in 00 01 02 03 04; do SDOS+=("${FIX}/dosages/anc_${i}.dosage.txt.gz"); done
  Rscript testdata/run_oracle.R \
    --null-rds "${FIX}/null.rds" --out "${WORKDIR}/small_r.tsv" \
    --ac-threshold 10 --n-core 1 --chunk-size 16 \
    --oracle-r /work/oracle/TractorMix.score.R --dosage-files "${SDOS[@]}"
  "${RUST_BIN}" \
    --null-export "${FIX}/null_export" --dosage-files "${SDOS[@]}" \
    --out "${WORKDIR}/small_rust.tsv" \
    --ac-threshold 10 --chunk-size 16 --threads 1
  python3 testdata/compare_oracle_tsv.py "${WORKDIR}/small_r.tsv" "${WORKDIR}/small_rust.tsv"
  if cmp -s "${WORKDIR}/small_r.tsv" "${WORKDIR}/small_rust.tsv"; then
    echo "Small fixture: raw TSV byte-identical"
  fi
fi

echo ""
echo "=== Large benchmark fixture ==="
LARGE="${WORKDIR}/large"
Rscript testdata/bench_large_fixture.R "${LARGE}" "${N_SAMPLES}" "${N_SITES}" "${N_ANC}" 2>/dev/null | tail -1

LDOS=()
for ((i = 0; i < N_ANC; i++)); do
  LDOS+=("${LARGE}/dosages/anc_$(printf '%02d' "${i}").dosage.txt.gz")
done
if [[ ${#LDOS[@]} -ne "${N_ANC}" ]]; then
  echo "Expected ${N_ANC} dosage files, found ${#LDOS[@]}" >&2
  exit 1
fi

echo ""
echo "=== Timing (wall clock, seconds) ==="
R_TIME=$(time_score "R" Rscript testdata/run_oracle.R \
  --null-rds "${LARGE}/null.rds" --out "${LARGE}/r.tsv" \
  --ac-threshold "${AC_THRESHOLD}" --n-core "${R_CORES}" --chunk-size "${R_CHUNK}" \
  --oracle-r /work/oracle/TractorMix.score.R --dosage-files "${LDOS[@]}")
RUST_TIME=$(time_score "Rust" "${RUST_BIN}" \
  --null-export "${LARGE}/null_export" --dosage-files "${LDOS[@]}" \
  --out "${LARGE}/rust.tsv" \
  --ac-threshold "${AC_THRESHOLD}" --chunk-size "${RUST_CHUNK}" --threads "${RUST_THREADS}")

echo "${R_TIME}"
echo "${RUST_TIME}"

python3 - <<PY
r = float("${R_TIME#R=}")
rs = float("${RUST_TIME#Rust=}")
sites = int("${N_SITES}")
print(f"Speedup (R ${R_CORES} core vs Rust ${RUST_THREADS} thread): {r/rs:.1f}x")
print(f"Rust throughput: {sites/rs:.0f} variants/sec")
PY

echo ""
echo "=== Parity report (large benchmark) ==="
python3 testdata/parity_report.py "${LARGE}/r.tsv" "${LARGE}/rust.tsv"
'

echo ""
echo "Benchmark finished."
