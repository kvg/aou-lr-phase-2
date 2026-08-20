#!/usr/bin/env bash
set -euo pipefail

# Realistic-scale R vs Rust benchmark inside the pilot Docker image.
# Defaults approximate AoU pilot scale (thousands of samples, tens of thousands
# of variants). Tune upward toward full chr22 before Terra re-submission.
#
# Usage:
#   ./scripts/bench_realistic_in_docker.sh
#
# Env:
#   N_SAMPLES=4000 N_SITES=30000 N_ANC=5
#   R_CORES=4          # matches legacy TractorMix.score n_core
#   RUST_THREADS=8     # Terra Score task cpus
#   IMAGE=.../tractor-mix-pilot:0.4.2

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CRATE_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
IMAGE="${IMAGE:-us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/tractor-mix-pilot:0.4.2}"
PLATFORM="${PLATFORM:-linux/amd64}"
N_SAMPLES="${N_SAMPLES:-2000}"
N_SITES="${N_SITES:-20000}"
N_ANC="${N_ANC:-5}"
AC_THRESHOLD="${AC_THRESHOLD:-50}"
R_CORES="${R_CORES:-4}"
RUST_THREADS="${RUST_THREADS:-8}"
DOCKER_CPUS="${DOCKER_CPUS:-8}"
# On laptops Docker may expose fewer CPUs; cap R n_core below detectCores().
if [[ "${DOCKER_CPUS}" =~ ^[0-9]+$ ]] && [[ "${DOCKER_CPUS}" -lt 4 ]]; then
  R_CORES="${R_CORES:-1}"
  RUST_THREADS="${RUST_THREADS:-${DOCKER_CPUS}}"
fi

echo "Realistic benchmark in ${IMAGE}"
echo "Scale: n=${N_SAMPLES}, sites=${N_SITES}, ancestries=${N_ANC}"
echo "Compare: R n_core=${R_CORES} vs Rust threads=${RUST_THREADS} (docker --cpus ${DOCKER_CPUS})"

docker run --rm --platform "${PLATFORM}" --cpus="${DOCKER_CPUS}" \
  -v "${CRATE_DIR}:/work" \
  -w /work \
  -e N_SAMPLES="${N_SAMPLES}" \
  -e N_SITES="${N_SITES}" \
  -e N_ANC="${N_ANC}" \
  -e AC_THRESHOLD="${AC_THRESHOLD}" \
  -e R_CORES="${R_CORES}" \
  -e RUST_THREADS="${RUST_THREADS}" \
  "${IMAGE}" \
  bash -lc '
set -euo pipefail
WORKDIR=/tmp/realistic_bench
FIX="${WORKDIR}/fixture"
mkdir -p "${WORKDIR}"

echo "=== Build fixture (glmmkin null + dosages) ==="
Rscript testdata/bench_large_fixture.R "${FIX}" "${N_SAMPLES}" "${N_SITES}" "${N_ANC}" 2>/dev/null | tail -1

LDOS=()
for ((i = 0; i < N_ANC; i++)); do
  LDOS+=("${FIX}/dosages/anc_$(printf "%02d" "${i}").dosage.txt.gz")
done

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
    echo "${label}: FAILED (exit ${status})" >&2
    exit 1
  fi
  echo "${label}=${sec}"
}

echo ""
echo "=== Score timing (production-like configs) ==="
T_R=$(time_score "R" Rscript testdata/run_oracle.R \
  --null-rds "${FIX}/null.rds" --out "${WORKDIR}/r.tsv" \
  --ac-threshold "${AC_THRESHOLD}" --n-core "${R_CORES}" --chunk-size 256 \
  --oracle-r /work/oracle/TractorMix.score.R --dosage-files "${LDOS[@]}")
T_RS1=$(time_score "Rust1" /usr/local/bin/tractor-mix-score \
  --null-export "${FIX}/null_export" --dosage-files "${LDOS[@]}" \
  --out "${WORKDIR}/rust1.tsv" --ac-threshold "${AC_THRESHOLD}" \
  --chunk-size 2048 --threads 1)
T_RS8=$(time_score "Rust8" /usr/local/bin/tractor-mix-score \
  --null-export "${FIX}/null_export" --dosage-files "${LDOS[@]}" \
  --out "${WORKDIR}/rust8.tsv" --ac-threshold "${AC_THRESHOLD}" \
  --chunk-size 2048 --threads "${RUST_THREADS}")

echo "${T_R}"
echo "${T_RS1}"
echo "${T_RS8}"

python3 - <<PY
import os
r = float("${T_R#R=}")
rs1 = float("${T_RS1#Rust1=}")
rs8 = float("${T_RS8#Rust8=}")
sites = int("${N_SITES}")
r_cores = int(os.environ["R_CORES"])
rust_threads = int(os.environ["RUST_THREADS"])
print(f"Speedup vs R ({r_cores} core): {r/rs1:.1f}x (Rust 1 thread), {r/rs8:.1f}x (Rust {rust_threads} thread)")
print(f"Throughput Rust {rust_threads}-thread: {sites/rs8:.0f} variants/sec")
if rs8 >= r:
    print("FAIL: Rust is not faster than parallel R at this scale — do not submit Terra yet.")
    raise SystemExit(1)
print("PASS: Rust beats parallel R on this fixture.")
PY

echo ""
echo "=== Parity (Rust 8-thread vs R) ==="
python3 testdata/parity_report.py "${WORKDIR}/r.tsv" "${WORKDIR}/rust8.tsv"

echo ""
echo "=== Parallel determinism (Rust 1 vs 8 thread) ==="
/usr/local/bin/tractor-mix-score \
  --null-export "${FIX}/null_export" --dosage-files "${LDOS[@]}" \
  --out "${WORKDIR}/det1.tsv" --ac-threshold "${AC_THRESHOLD}" \
  --chunk-size 2048 --threads 1 >/dev/null
/usr/local/bin/tractor-mix-score \
  --null-export "${FIX}/null_export" --dosage-files "${LDOS[@]}" \
  --out "${WORKDIR}/det8.tsv" --ac-threshold "${AC_THRESHOLD}" \
  --chunk-size 2048 --threads "${RUST_THREADS}" >/dev/null
cmp -s "${WORKDIR}/det1.tsv" "${WORKDIR}/det8.tsv"
echo "OK: 1-thread and ${RUST_THREADS}-thread outputs byte-identical"
'

echo "Realistic benchmark passed."
