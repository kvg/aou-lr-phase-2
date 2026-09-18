#!/usr/bin/env bash
# Build (and optionally push) the bw2 ExpansionHunter image.
#
# Defaults:
#   image  us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/aou-expansion-hunter:0.1.0
#
# Usage:
#   ./build_docker.sh
#   ./build_docker.sh --local
#   ./build_docker.sh --local --push
#   ./build_docker.sh --tag 0.1.0

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DOCKER_DIR="${SCRIPT_DIR}/docker"
CONTEXT_DIR="${SCRIPT_DIR}"

DEFAULT_PROJECT="broad-dsp-lrma"
DEFAULT_REGION="us-central1"
DEFAULT_IMAGE="us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/aou-expansion-hunter"
DEFAULT_TAG="0.1.0"

IMAGE="${IMAGE:-${DEFAULT_IMAGE}}"
TAG="${TAG:-${DEFAULT_TAG}}"
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
Build the ExpansionHunter image (bw2 fork; no baked FASTA).

Defaults:
  project  ${DEFAULT_PROJECT}
  region   ${DEFAULT_REGION}
  image    ${DEFAULT_IMAGE}:${DEFAULT_TAG}

Options:
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
      --platform linux/amd64 \
      -f "${DOCKER_DIR}/Dockerfile" \
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
      --substitutions="_IMAGE=${IMAGE},_TAG=${TAG}" \
      --timeout="${TIMEOUT}" \
      --machine-type="${MACHINE_TYPE}" \
      --disk-size="${DISK_SIZE_GB}" \
      --project="${PROJECT}"
    echo "Cloud Build finished. Image: ${FULL_IMAGE}"
    echo "Point ExpansionHunterMinicram.eh_docker at ${FULL_IMAGE}."
    ;;

  *)
    echo "Invalid MODE=${MODE}" >&2
    exit 1
    ;;
esac
