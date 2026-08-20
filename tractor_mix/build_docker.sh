#!/usr/bin/env bash
# Build (and optionally push) the Tractor-Mix pilot Docker image.
#
# Defaults (from gcloud config / Artifact Registry):
#   project  broad-dsp-lrma
#   repo     us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr
#   image    .../aou-lr/tractor-mix-pilot:0.4.2
#
# 0.4.2 implements parallel variant scoring (--threads) + progress logging.
# 0.4.1 fixes CSC colptr export in fit_null.R (Sigma_i@p is already 0-based).
# 0.4.0 adds Alpine-static tractor-mix-score (Rust TractorMix.score sparse GRM).
# 0.3.2 Alpine-static extract-tracts-flare (verified to exec on SAIGE base).
# 0.3.1 bookworm-cross musl was often dynamically linked (ENOENT on SAIGE).
# 0.3.0 gnu/bookworm needed GLIBC 2.34.
# 0.2.0 added pinned SAIGE (+ deps) alongside Tractor-Mix / GMMAT.
#
# Usage:
#   ./build_docker.sh                     # Cloud Build → push default IMAGE:TAG
#   ./build_docker.sh --local             # local docker build (no push)
#   ./build_docker.sh --local --push      # local build + docker push
#   ./build_docker.sh --tag 0.2.1         # bump tag only
#
# Env overrides: IMAGE, TAG, PROJECT, REGION, MACHINE_TYPE, DISK_SIZE_GB, TIMEOUT

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DOCKER_DIR="${SCRIPT_DIR}/docker"
CONTEXT_DIR="${SCRIPT_DIR}"

# Discovered via: gcloud config get-value project
#               + gcloud artifacts repositories list (aou-lr @ us-central1)
DEFAULT_PROJECT="broad-dsp-lrma"
DEFAULT_REGION="us-central1"
DEFAULT_IMAGE="us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/tractor-mix-pilot"
DEFAULT_TAG="0.4.2"

IMAGE="${IMAGE:-${DEFAULT_IMAGE}}"
TAG="${TAG:-${DEFAULT_TAG}}"
PROJECT="${PROJECT:-${DEFAULT_PROJECT}}"
REGION="${REGION:-${DEFAULT_REGION}}"
MACHINE_TYPE="${MACHINE_TYPE:-e2-highcpu-32}"
DISK_SIZE_GB="${DISK_SIZE_GB:-200}"
TIMEOUT="${TIMEOUT:-7200s}"

MODE="cloudbuild"   # cloudbuild | local
DO_PUSH=0
DRY_RUN=0

usage() {
  cat <<EOF
Build the Tractor-Mix pilot image.

Defaults:
  project  ${DEFAULT_PROJECT}
  region   ${DEFAULT_REGION}
  image    ${DEFAULT_IMAGE}:${DEFAULT_TAG}

Options:
  --local              Build with local docker (default: Google Cloud Build)
  --push               After --local build, docker push IMAGE:TAG
  --cloudbuild         Force Cloud Build (default)
  --image NAME         Image repository path without tag (default: ${IMAGE})
  --tag TAG            Image tag (default: ${TAG})
  --project PROJECT    GCP project for Cloud Build (default: ${PROJECT})
  --region REGION      Artifact Registry region (default: ${REGION})
  --machine-type TYPE  Cloud Build machine type (default: ${MACHINE_TYPE})
  --disk-size GB       Cloud Build disk size GB (default: ${DISK_SIZE_GB})
  --timeout DURATION   Cloud Build timeout (default: ${TIMEOUT})
  --dry-run            Print commands only
  -h, --help           Show this help

Examples:
  ./build_docker.sh                      # Cloud Build → default Artifact Registry image
  ./build_docker.sh --tag 0.2.1          # same repo, new tag
  ./build_docker.sh --local              # local docker build only
  ./build_docker.sh --local --push       # local build + push to default image
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --local) MODE="local"; shift ;;
    --cloudbuild) MODE="cloudbuild"; shift ;;
    --push) DO_PUSH=1; shift ;;
    --image) IMAGE="$2"; shift 2 ;;
    --tag) TAG="$2"; shift 2 ;;
    --project) PROJECT="$2"; shift 2 ;;
    --region) REGION="$2"; shift 2 ;;
    --machine-type) MACHINE_TYPE="$2"; shift 2 ;;
    --disk-size) DISK_SIZE_GB="$2"; shift 2 ;;
    --timeout) TIMEOUT="$2"; shift 2 ;;
    --dry-run) DRY_RUN=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 1 ;;
  esac
done

FULL_IMAGE="${IMAGE}:${TAG}"

run() {
  if [[ "${DRY_RUN}" -eq 1 ]]; then
    printf '+'
    printf ' %q' "$@"
    printf '\n'
  else
    "$@"
  fi
}

if [[ ! -f "${DOCKER_DIR}/Dockerfile" ]]; then
  echo "Dockerfile not found at ${DOCKER_DIR}/Dockerfile" >&2
  exit 1
fi

REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
SHARED_SCRIPTS="${REPO_ROOT}/scripts"
STAGE="${SCRIPT_DIR}/.docker_scripts"
if [[ ! -d "${SHARED_SCRIPTS}" ]]; then
  echo "Shared scripts not found at ${SHARED_SCRIPTS}" >&2
  exit 1
fi
rm -rf "${STAGE}"
mkdir -p "${STAGE}"
cp -a "${SHARED_SCRIPTS}/." "${STAGE}/"
echo "Staged scripts: ${SHARED_SCRIPTS} -> ${STAGE}"

echo "Mode:        ${MODE}"
echo "Project:     ${PROJECT}"
echo "Image:       ${FULL_IMAGE}"
echo "Context:     ${CONTEXT_DIR}"
echo "Dockerfile:  ${DOCKER_DIR}/Dockerfile"

case "${MODE}" in
  local)
    if ! command -v docker >/dev/null 2>&1; then
      echo "docker not found; use Cloud Build instead (omit --local)" >&2
      exit 1
    fi
    run docker build \
      -f "${DOCKER_DIR}/Dockerfile" \
      -t "${FULL_IMAGE}" \
      "${CONTEXT_DIR}"
    if [[ "${DO_PUSH}" -eq 1 ]]; then
      run docker push "${FULL_IMAGE}"
    else
      echo "Built locally. Push with: docker push ${FULL_IMAGE}"
      echo "Or re-run: $0 --local --push --image ${IMAGE} --tag ${TAG}"
    fi
    ;;

  cloudbuild)
    if ! command -v gcloud >/dev/null 2>&1; then
      echo "gcloud not found; install Google Cloud SDK or use --local" >&2
      exit 1
    fi

    CB_ARGS=(
      builds submit
      "${CONTEXT_DIR}"
      --config="${DOCKER_DIR}/cloudbuild.yaml"
      --substitutions="_IMAGE=${IMAGE},_TAG=${TAG}"
      --timeout="${TIMEOUT}"
      --machine-type="${MACHINE_TYPE}"
      --disk-size="${DISK_SIZE_GB}"
      --project="${PROJECT}"
    )

    echo "Submitting Cloud Build (this image is large; expect 30–90+ minutes)..."
    run gcloud "${CB_ARGS[@]}"
    echo "Cloud Build finished. Image: ${FULL_IMAGE}"
    echo "WDL runtime.docker default already matches this image path."
    ;;

  *)
    echo "Invalid MODE=${MODE}" >&2
    exit 1
    ;;
esac
