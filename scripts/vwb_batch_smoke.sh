#!/usr/bin/env bash
# Native Google Batch smoke (no dsub sidecars).
#
# hello-world-3 proved Batch VMs start in this project. This job uses the
# same PET SA + private VPC as dsub, one print-reads container, and copies a
# tiny in-perimeter GCS object. If this SUCCEEDS, drop dsub and submit
# minicram the same way.
#
# In a VWB Jupyter terminal:
#   ./scripts/vwb_batch_smoke.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
# shellcheck source=vwb_dsub.sh
source "${REPO_ROOT}/scripts/vwb_dsub.sh"

vwb_dsub_require_env

REGION="${DSUB_REGION:-us-central1}"
IMAGE="${PRINT_READS_DOCKER:-${DEFAULT_DSUB_CLOUD_SDK_IMAGE}}"
BUCKET="${VWB_BATCH_SMOKE_BUCKET:-gs://aou-lr-phase2-resources}"
IN_URI="${VWB_BATCH_SMOKE_IN:-${BUCKET%/}/expansion_hunter/smoke.catalog.json}"
OUT_URI="${BUCKET%/}/dsub_example/batch_smoke.catalog.copy.json"
JOB_ID="eh-batch-smoke-$(date -u +%y%m%d-%H%M%S)"
CFG="$(mktemp)"
trap 'rm -f "${CFG}"' EXIT

export BATCH_SMOKE_IMAGE="${IMAGE}"
export BATCH_SMOKE_IN="${IN_URI}"
export BATCH_SMOKE_OUT="${OUT_URI}"
export DSUB_REGION="${REGION}"

python3 - "${CFG}" <<'PY'
import json, os, sys

cfg_path = sys.argv[1]
project = os.environ["GOOGLE_CLOUD_PROJECT"]
region = os.environ["DSUB_REGION"]
image = os.environ["BATCH_SMOKE_IMAGE"]
in_uri = os.environ["BATCH_SMOKE_IN"]
out_uri = os.environ["BATCH_SMOKE_OUT"]
sa = os.environ["PET_SA_EMAIL"]
# 0.1.2 has python + google-cloud-storage, not bash. gcloud's shebang is bash
# → /usr/bin/gcloud is exit 127. Copy with the GCS client instead.
py = (
    "from google.cloud import storage\n"
    f"in_uri, out_uri = {in_uri!r}, {out_uri!r}\n"
    "def split(u):\n"
    "    b, _, p = u[5:].partition('/')\n"
    "    return b, p\n"
    "c = storage.Client()\n"
    "sb, so = split(in_uri)\n"
    "db, do = split(out_uri)\n"
    "c.bucket(db).blob(do).upload_from_string("
    "c.bucket(sb).blob(so).download_as_bytes())\n"
    "print('copied', in_uri, '->', out_uri)\n"
)
job = {
    "taskGroups": [
        {
            "taskCount": 1,
            "taskSpec": {
                "computeResource": {"cpuMilli": 2000, "memoryMib": 4096},
                "maxRetryCount": 0,
                "runnables": [
                    {
                        "container": {
                            "imageUri": image,
                            "entrypoint": "python3",
                            "commands": ["-c", py],
                        }
                    }
                ],
            },
        }
    ],
    "allocationPolicy": {
        "location": {"allowedLocations": [f"regions/{region}"]},
        "serviceAccount": {"email": sa},
        "network": {
            "networkInterfaces": [
                {
                    "network": f"projects/{project}/global/networks/network",
                    "subnetwork": f"projects/{project}/regions/{region}/subnetworks/subnetwork",
                    "noExternalIpAddress": True,
                }
            ]
        },
    },
    "logsPolicy": {"destination": "CLOUD_LOGGING"},
}
json.dump(job, open(cfg_path, "w"), indent=2)
print(open(cfg_path).read())
PY

echo "JOB_ID=${JOB_ID}"
echo "IMAGE=${IMAGE}"
echo "IN=${IN_URI}"
echo "OUT=${OUT_URI}"

gcloud batch jobs submit "${JOB_ID}" \
  --location="${REGION}" \
  --project="${GOOGLE_CLOUD_PROJECT}" \
  --config="${CFG}"

echo "Watch:"
echo "  gcloud batch jobs describe ${JOB_ID} --location=${REGION} --format='yaml(status.state,status.statusEvents)'"
echo "When SUCCEEDED:"
echo "  gcloud storage ls ${OUT_URI}"
