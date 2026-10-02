#!/usr/bin/env bash
# Build (and optionally push) the FLARE bcftools/htslib image.
#
# Defaults:
#   image  us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/aou-flare-bcftools:1.24
#
# Usage:
#   ./build_docker.sh
#   ./build_docker.sh --local
#   ./build_docker.sh --local --push
#   ./build_docker.sh --tag 1.24
#   ./build_docker.sh --flare2       # FLARE 0.6 + FLARE2 clustering image (aou-flare2)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DOCKER_DIR="${SCRIPT_DIR}/docker"
CONTEXT_DIR="${SCRIPT_DIR}"

DEFAULT_PROJECT="broad-dsp-lrma"
DEFAULT_REGION="us-central1"
DEFAULT_IMAGE="us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/aou-flare-bcftools"
DEFAULT_TAG="1.24"
FLARE2_IMAGE="us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/aou-flare2"
FLARE2_TAG="0.6.0-87573be"
DOCKERFILE_NAME="Dockerfile"
FLARE2=0

IMAGE="${IMAGE:-}"
TAG="${TAG:-}"
PROJECT="${PROJECT:-${DEFAULT_PROJECT}}"
REGION="${REGION:-${DEFAULT_REGION}}"
MACHINE_TYPE="${MACHINE_TYPE:-e2-highcpu-8}"
DISK_SIZE_GB="${DISK_SIZE_GB:-100}"
TIMEOUT="${TIMEOUT:-3600s}"

MODE="cloudbuild"
DO_PUSH=0
DRY_RUN=0

usage() {
  cat <<EOF
Build the FLARE bcftools image (htslib libcurl retry, PR 1987).

Defaults:
  project  ${DEFAULT_PROJECT}
  region   ${DEFAULT_REGION}
  image    ${DEFAULT_IMAGE}:${DEFAULT_TAG}

Options:
  --flare2             Build docker/Dockerfile.flare2 (${FLARE2_IMAGE}:${FLARE2_TAG})
  --local              Build with local docker (default: Google Cloud Build)
  --push               After --local build, docker push IMAGE:TAG
  --cloudbuild         Force Cloud Build (default)
  --image NAME         Image repository path without tag
  --tag TAG            Image tag (default: ${DEFAULT_TAG})
  --project PROJECT    GCP project for Cloud Build
  --region REGION      Artifact Registry region
  --machine-type TYPE  Cloud Build machine type
  --disk-size GB       Cloud Build disk size GB
  --timeout DURATION   Cloud Build timeout
  --dry-run            Print commands only
  -h, --help           Show this help
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --local) MODE="local"; shift ;;
    --flare2) FLARE2=1; shift ;;
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

if [[ "${FLARE2}" -eq 1 ]]; then
  DOCKERFILE_NAME="Dockerfile.flare2"
  IMAGE="${IMAGE:-${FLARE2_IMAGE}}"
  TAG="${TAG:-${FLARE2_TAG}}"
else
  IMAGE="${IMAGE:-${DEFAULT_IMAGE}}"
  TAG="${TAG:-${DEFAULT_TAG}}"
fi
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

if [[ ! -f "${DOCKER_DIR}/${DOCKERFILE_NAME}" ]]; then
  echo "Dockerfile not found at ${DOCKER_DIR}/${DOCKERFILE_NAME}" >&2
  exit 1
fi

echo "Mode:        ${MODE}"
echo "Project:     ${PROJECT}"
echo "Image:       ${FULL_IMAGE}"
echo "Context:     ${CONTEXT_DIR}"
echo "Dockerfile:  ${DOCKER_DIR}/${DOCKERFILE_NAME}"

case "${MODE}" in
  local)
    if ! command -v docker >/dev/null 2>&1; then
      echo "docker not found; use Cloud Build instead (omit --local)" >&2
      exit 1
    fi
    run docker build \
      --platform linux/amd64 \
      -f "${DOCKER_DIR}/${DOCKERFILE_NAME}" \
      -t "${FULL_IMAGE}" \
      "${CONTEXT_DIR}"
    if [[ "${DO_PUSH}" -eq 1 ]]; then
      run gcloud auth configure-docker "${REGION}-docker.pkg.dev" --quiet
      run docker push "${FULL_IMAGE}"
      echo "Pushed ${FULL_IMAGE}"
    else
      echo "Built locally: ${FULL_IMAGE}"
    fi
    ;;

  cloudbuild)
    if ! command -v gcloud >/dev/null 2>&1; then
      echo "gcloud not found; install Google Cloud SDK or use --local" >&2
      exit 1
    fi
    run gcloud builds submit \
      "${CONTEXT_DIR}" \
      --config="${DOCKER_DIR}/cloudbuild.yaml" \
      --substitutions="_IMAGE=${IMAGE},_TAG=${TAG},_DOCKERFILE=docker/${DOCKERFILE_NAME}" \
      --timeout="${TIMEOUT}" \
      --machine-type="${MACHINE_TYPE}" \
      --disk-size="${DISK_SIZE_GB}" \
      --project="${PROJECT}"
    echo "Cloud Build finished. Image: ${FULL_IMAGE}"
    if [[ "${FLARE2}" -eq 1 ]]; then
      echo "Point FlareByPopulation.flare2_docker at ${FULL_IMAGE}."
    else
      echo "Point FlareByPopulation.bcftools_docker at ${FULL_IMAGE}."
    fi
    ;;

  *)
    echo "Invalid MODE=${MODE}" >&2
    exit 1
    ;;
esac
