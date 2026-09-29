"""VWB PET SA + private VPC allocation for native Google Batch."""

from __future__ import annotations

import os
import re
from datetime import datetime, timezone


def require_vwb_env() -> tuple[str, str]:
    project = os.environ.get("GOOGLE_CLOUD_PROJECT", "").strip()
    sa = os.environ.get("PET_SA_EMAIL", "").strip()
    if not project:
        raise SystemExit("GOOGLE_CLOUD_PROJECT is empty. Open this in a VWB Jupyter app.")
    if not sa:
        raise SystemExit("PET_SA_EMAIL is empty. VWB sets this in cloud apps.")
    return project, sa


def region() -> str:
    return os.environ.get("DSUB_REGION", os.environ.get("BATCH_REGION", "us-central1")).strip()


def allocation_policy(project: str, region_name: str, sa: str) -> dict:
    return {
        "location": {"allowedLocations": [f"regions/{region_name}"]},
        "serviceAccount": {"email": sa},
        "network": {
            "networkInterfaces": [
                {
                    "network": f"projects/{project}/global/networks/network",
                    "subnetwork": f"projects/{project}/regions/{region_name}/subnetworks/subnetwork",
                    "noExternalIpAddress": True,
                }
            ]
        },
    }


def sanitize_label(value: str, *, max_len: int = 63) -> str:
    cleaned = re.sub(r"[^a-z0-9_-]+", "-", value.lower()).strip("-")
    return (cleaned or "x")[:max_len]


def sanitize_job_id(value: str) -> str:
    cleaned = re.sub(r"[^a-z0-9-]+", "-", value.lower()).strip("-")
    if not cleaned or not cleaned[0].isalpha():
        cleaned = f"j-{cleaned}"
    return cleaned[:63].rstrip("-")


def job_id(*, pipeline: str, stage: str, run_id: str, shard: int) -> str:
    shorts = {
        "expansion_hunter": "eh",
        "expansion-hunter": "eh",
        "locityper": "lt",
        "minicram": "mc",
        "genotype": "gt",
        "print_reads": "pr",
    }
    stamp = datetime.now(timezone.utc).strftime("%H%M%S")
    pipe = shorts.get(pipeline, pipeline[:2])
    st = shorts.get(stage, stage[:2])
    run = sanitize_label(run_id, max_len=20)
    return sanitize_job_id(f"{pipe}-{st}-{run}-s{shard:03d}-{stamp}")


def job_labels(*, pipeline: str, run_id: str, stage: str, shard: int) -> dict[str, str]:
    return {
        "pipeline": sanitize_label(pipeline.replace("_", "-")),
        "run-id": sanitize_label(run_id),
        "stage": sanitize_label(stage.replace("_", "-")),
        "shard": str(int(shard)),
    }
