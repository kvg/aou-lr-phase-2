#!/usr/bin/env bash
# Shared Verily Workbench / Cloud Batch flags for dsub.
# Source from submit scripts. Requires GOOGLE_CLOUD_PROJECT and PET_SA_EMAIL
# (both set in VWB Jupyter apps).
#
# VWB flags (PET SA, private VPC, --use-private-address) match
# https://support.workbench.verily.com/docs/guides/workflows/dsub/
#
# dsub google-batch copies --input/--output/--logging from a sidecar
# (DSUB_CLOUD_SDK_IMAGE), not from --image. Default gcr.io cloud-sdk is
# outside AoU VPC-SC, so we use the print_reads Broad AR tag (python3 +
# gcloud + str-analysis).

# Same Broad AR tag as minicram --image: python3 + gcloud (dsub prepare /
# localize / --logging) and str-analysis. Google's gcr.io / pkg.dev cloud-sdk
# images are outside the AoU VPC-SC perimeter.
DEFAULT_DSUB_CLOUD_SDK_IMAGE="us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/aou-locityper-print-reads:0.1.2"

vwb_dsub_require_env() {
  if [[ -z "${GOOGLE_CLOUD_PROJECT:-}" ]]; then
    echo "GOOGLE_CLOUD_PROJECT is empty. Open this in a VWB Jupyter app." >&2
    exit 1
  fi
  if [[ -z "${PET_SA_EMAIL:-}" ]]; then
    echo "PET_SA_EMAIL is empty. VWB sets this in cloud apps." >&2
    exit 1
  fi
  # Must run in the submitting shell. mapfile < <(vwb_dsub_base_args) is a
  # subshell, so an export there never reaches dsub.
  vwb_dsub_export_cloud_sdk_image
}

vwb_dsub_export_cloud_sdk_image() {
  export DSUB_CLOUD_SDK_IMAGE="${DSUB_CLOUD_SDK_IMAGE:-${DEFAULT_DSUB_CLOUD_SDK_IMAGE}}"
  echo "DSUB_CLOUD_SDK_IMAGE=${DSUB_CLOUD_SDK_IMAGE}" >&2
}

# Prints base dsub args (one per line) for `mapfile` / xargs.
vwb_dsub_base_args() {
  vwb_dsub_require_env
  local region="${DSUB_REGION:-us-central1}"
  printf '%s\n' \
    --provider google-batch \
    --project "${GOOGLE_CLOUD_PROJECT}" \
    --regions "${region}" \
    --service-account "${PET_SA_EMAIL}" \
    --network "projects/${GOOGLE_CLOUD_PROJECT}/global/networks/network" \
    --subnetwork "projects/${GOOGLE_CLOUD_PROJECT}/regions/${region}/subnetworks/subnetwork" \
    --use-private-address \
    --env "GOOGLE_CLOUD_PROJECT=${GOOGLE_CLOUD_PROJECT}"
}

vwb_dsub_parse_job_id() {
  sed -n 's/.*Launched job-id:[[:space:]]*//p' | tail -n 1
}

# dsub --logging DIR writes DIR/{job-id}.{task-id}.log plus -stdout/-stderr.
vwb_dsub_print_log_uris() {
  local log_root="$1" job_id="$2" task="${3:-1}"
  local prefix="${log_root%/}/${job_id}.${task}"
  echo "dsub --logging objects (not {sample}.minicram.worker.log):" >&2
  printf '  %s\n' "${prefix}.log" "${prefix}-stdout.log" "${prefix}-stderr.log" >&2
  echo "Batch job id: ${job_id}-${task}-0" >&2
  echo "  gcloud batch jobs describe ${job_id}-${task}-0 --project=${GOOGLE_CLOUD_PROJECT} --location=${DSUB_REGION:-us-central1}" >&2
}
