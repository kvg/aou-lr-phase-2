#!/usr/bin/env bash
# Shared Verily Workbench / Cloud Batch flags for dsub.
# Source from submit scripts. Requires GOOGLE_CLOUD_PROJECT and PET_SA_EMAIL
# (both set in VWB Jupyter apps).

vwb_dsub_require_env() {
  if [[ -z "${GOOGLE_CLOUD_PROJECT:-}" ]]; then
    echo "GOOGLE_CLOUD_PROJECT is empty. Open this in a VWB Jupyter app." >&2
    exit 1
  fi
  if [[ -z "${PET_SA_EMAIL:-}" ]]; then
    echo "PET_SA_EMAIL is empty. VWB sets this in cloud apps." >&2
    exit 1
  fi
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
    --use-private-address
}

vwb_dsub_parse_job_id() {
  sed -n 's/.*Launched job-id:[[:space:]]*//p' | tail -n 1
}
