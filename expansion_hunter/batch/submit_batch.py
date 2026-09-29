#!/usr/bin/env python3
"""Build a native Google Batch job JSON for ExpansionHunter minicram.

One print-reads container per sample. No dsub sidecars. CRAM stays gs://.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path


def _gs(prefix: str, sample_id: str, name: str) -> str:
    return f"{prefix.rstrip('/')}/{sample_id}/{name}"


def _allocation(project: str, region: str, sa: str) -> dict:
    return {
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
    }


def _rows(csv_path: Path, out_prefix: str, project: str) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with csv_path.open(newline="") as fh:
        for rec in csv.DictReader(fh):
            sid = rec["sample_id"]
            rows.append(
                {
                    "SAMPLE_ID": sid,
                    "CRAM": rec["cram"],
                    "CRAI": rec["crai"],
                    "REF_FA": rec["ref_fa"],
                    "REF_FAI": rec["ref_fai"],
                    "CATALOG": rec["catalog"],
                    "GCLOUD_PROJECT": rec.get("gcloud_project")
                    or rec.get("GCLOUD_PROJECT")
                    or project,
                    "WORKER_LOG_GCS": _gs(out_prefix, sid, f"{sid}.minicram.worker.log"),
                    "MINICRAM": _gs(out_prefix, sid, f"{sid}.minicram.cram"),
                    "MINICRAI": _gs(out_prefix, sid, f"{sid}.minicram.cram.crai"),
                    "TRANSFER_STATS": _gs(out_prefix, sid, f"{sid}.data_transfer_stats.tsv"),
                    "WINDOW_SIZE": os.environ.get("WINDOW_SIZE", "1000"),
                    "MERGE_REGIONS_DISTANCE": os.environ.get("MERGE_REGIONS_DISTANCE", "1000"),
                    "MAX_RETRY": os.environ.get("MAX_RETRY", "3"),
                    "WAIT_TIME": os.environ.get("WAIT_TIME", "30"),
                }
            )
    if not rows:
        raise SystemExit(f"no rows in {csv_path}")
    return rows


def build_minicram_job(
    *,
    csv_path: Path,
    worker_path: Path,
    image: str,
    out_prefix: str,
    project: str,
    region: str,
    sa: str,
) -> dict:
    worker = worker_path.read_text(encoding="utf-8")
    rows = _rows(csv_path, out_prefix, project)
    common = {
        "WINDOW_SIZE": os.environ.get("WINDOW_SIZE", "1000"),
        "MERGE_REGIONS_DISTANCE": os.environ.get("MERGE_REGIONS_DISTANCE", "1000"),
        "MAX_RETRY": os.environ.get("MAX_RETRY", "3"),
        "WAIT_TIME": os.environ.get("WAIT_TIME", "30"),
    }
    spec: dict = {
        "computeResource": {
            "cpuMilli": 4000,
            "memoryMib": 16384,
            "bootDiskMib": 51200,
        },
        "maxRetryCount": 0,
        "maxRunDuration": "14400s",
        "environment": {"variables": common},
        "runnables": [
            {
                "container": {
                    "imageUri": image,
                    "entrypoint": "python3",
                    "commands": ["-c", worker],
                }
            }
        ],
    }
    group: dict = {"taskCount": len(rows), "taskSpec": spec}
    if len(rows) == 1:
        spec["environment"]["variables"].update(rows[0])
    else:
        group["taskEnvironments"] = [{"variables": row} for row in rows]
    return {
        "taskGroups": [group],
        "allocationPolicy": _allocation(project, region, sa),
        "logsPolicy": {"destination": "CLOUD_LOGGING"},
    }


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--csv", type=Path, required=True)
    p.add_argument("--worker", type=Path, required=True)
    p.add_argument("--image", required=True)
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--project", required=True)
    p.add_argument("--region", default="us-central1")
    p.add_argument("--sa", required=True)
    p.add_argument("--config", type=Path, required=True)
    args = p.parse_args()
    job = build_minicram_job(
        csv_path=args.csv,
        worker_path=args.worker,
        image=args.image,
        out_prefix=args.out_prefix,
        project=args.project,
        region=args.region,
        sa=args.sa,
    )
    args.config.write_text(json.dumps(job, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {args.config} ({len(job['taskGroups'][0].get('taskEnvironments') or [1])} task(s))")


if __name__ == "__main__":
    main()
