#!/usr/bin/env bash
# Build (and optionally push) the TRGT haplotype-phasing Docker image.
#
# Used by wdl/TrgtHaplotypePhase.wdl: whatshap + TRGT + samtools/bcftools.
# Nothing from the repo is baked in — the WDL passes phase_trgt_by_haplotag.py
# as a File input from the workspace bucket — so the build context is just
# repeat_loci/ and no scripts are staged.
#
# Defaults:
#   project  broad-dsp-lrma
#   image    us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/aou-trgt-phase:0.1.0
#
# Usage:
#   repeat_loci/build_docker.sh                 # Cloud Build -> push default IMAGE:TAG
#   repeat_loci/build_docker.sh --local         # local docker build (no push)
#   repeat_loci/build_docker.sh --local --push  # local build + docker push
#   repeat_loci/build_docker.sh --tag 0.1.1     # bump tag only
#   repeat_loci/build_docker.sh --dry-run       # print the command and stop
#
# Env overrides: IMAGE, TAG, PROJECT, REGION, MACHINE_TYPE, DISK_SIZE_GB, TIMEOUT
#
# After a push, pin TrgtHaplotypePhase.docker to the @sha256 digest rather than
# the tag for any run whose outputs are kept: a re-pushed tag silently changes
# what ran. The digest is printed at the end of a --local --push build and is in
# the Cloud Build output otherwise.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DOCKER_DIR="${SCRIPT_DIR}/docker"
CONTEXT_DIR="${SCRIPT_DIR}"

DEFAULT_PROJECT="broad-dsp-lrma"
DEFAULT_REGION="us-central1"
DEFAULT_IMAGE="us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/aou-trgt-phase"
DEFAULT_TAG="0.1.0"

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
Build the TRGT haplotype-phasing image (whatshap + TRGT + samtools/bcftools).

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
    run docker tag "${FULL_IMAGE}" "aou-trgt-phase:${TAG}"
    if [[ "${DO_PUSH}" -eq 1 ]]; then
      run gcloud auth configure-docker "${REGION}-docker.pkg.dev" --quiet
      run docker push "${FULL_IMAGE}"
      echo "Pushed ${FULL_IMAGE}"
      echo "Digest:"
      run docker inspect --format='{{index .RepoDigests 0}}' "${FULL_IMAGE}"
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

    echo "Submitting Cloud Build..."
    run gcloud "${CB_ARGS[@]}"
    echo "Cloud Build finished. Image: ${FULL_IMAGE}"
    echo "Point TrgtHaplotypePhase.docker at ${FULL_IMAGE} (or its @sha256 digest)."
    ;;

  *)
    echo "Invalid MODE=${MODE}" >&2
    exit 1
    ;;
esac
