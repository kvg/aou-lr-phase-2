#!/usr/bin/env python3
"""Build native Google Batch job JSON for ExpansionHunter minicram / genotype.

Minicram: one print-reads container; CRAM stays gs://.
Genotype: host gcloud downloads the minicram + fasta, EH container runs
ExpansionHunter, host gcloud uploads JSON/VCF. The EH 0.1.0 image has no
GCS client.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path

WORKDIR = "/mnt/disks/eh"

DOWNLOAD_SCRIPT = r"""#!/bin/bash
set -eu
WORKDIR=/mnt/disks/eh
mkdir -p "${WORKDIR}"
chmod 777 "${WORKDIR}"
gcloud storage cp --billing-project="${GCLOUD_PROJECT}" "${MINICRAM}" "${WORKDIR}/reads.cram"
gcloud storage cp --billing-project="${GCLOUD_PROJECT}" "${MINICRAI}" "${WORKDIR}/reads.cram.crai"
gcloud storage cp --billing-project="${GCLOUD_PROJECT}" "${REF_FA}" "${WORKDIR}/reference.fa"
gcloud storage cp --billing-project="${GCLOUD_PROJECT}" "${REF_FAI}" "${WORKDIR}/reference.fa.fai"
gcloud storage cp --billing-project="${GCLOUD_PROJECT}" "${CATALOG}" "${WORKDIR}/catalog.json"
ls -lh "${WORKDIR}"
"""

UPLOAD_SCRIPT = r"""#!/bin/bash
set -eu
WORKDIR=/mnt/disks/eh
json="${WORKDIR}/${SAMPLE_ID}.EH.json"
vcf="${WORKDIR}/${SAMPLE_ID}.EH.vcf"
log="${WORKDIR}/worker.log"
if [ -s "${json}" ]; then
  gcloud storage cp --billing-project="${GCLOUD_PROJECT}" "${json}" "${EH_JSON}"
fi
if [ -s "${vcf}" ]; then
  gcloud storage cp --billing-project="${GCLOUD_PROJECT}" "${vcf}" "${EH_VCF}"
fi
if [ -s "${log}" ]; then
  gcloud storage cp --billing-project="${GCLOUD_PROJECT}" "${log}" "${WORKER_LOG_GCS}"
fi
ls -lh "${json}" "${vcf}" "${log}" 2>/dev/null || true
"""


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


def _project_for_row(rec: dict, project: str) -> str:
    return rec.get("gcloud_project") or rec.get("GCLOUD_PROJECT") or project


def minicram_rows(csv_path: Path, out_prefix: str, project: str) -> list[dict[str, str]]:
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
                    "GCLOUD_PROJECT": _project_for_row(rec, project),
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


def genotype_rows(csv_path: Path, out_prefix: str, project: str) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with csv_path.open(newline="") as fh:
        for rec in csv.DictReader(fh):
            sid = rec["sample_id"]
            rows.append(
                {
                    "SAMPLE_ID": sid,
                    "SEX": rec.get("sex") or "female",
                    "MINICRAM": _gs(out_prefix, sid, f"{sid}.minicram.cram"),
                    "MINICRAI": _gs(out_prefix, sid, f"{sid}.minicram.cram.crai"),
                    "REF_FA": rec["ref_fa"],
                    "REF_FAI": rec["ref_fai"],
                    "CATALOG": rec["catalog"],
                    "GCLOUD_PROJECT": _project_for_row(rec, project),
                    "WORKER_LOG_GCS": _gs(out_prefix, sid, f"{sid}.EH.worker.log"),
                    "EH_JSON": _gs(out_prefix, sid, f"{sid}.EH.json"),
                    "EH_VCF": _gs(out_prefix, sid, f"{sid}.EH.vcf"),
                }
            )
    if not rows:
        raise SystemExit(f"no rows in {csv_path}")
    return rows


def _attach_rows(group: dict, spec: dict, rows: list[dict[str, str]]) -> None:
    group["taskCount"] = len(rows)
    if len(rows) == 1:
        spec["environment"]["variables"].update(rows[0])
    else:
        group["taskEnvironments"] = [{"variables": row} for row in rows]


def _job(group: dict, project: str, region: str, sa: str) -> dict:
    return {
        "taskGroups": [group],
        "allocationPolicy": _allocation(project, region, sa),
        "logsPolicy": {"destination": "CLOUD_LOGGING"},
    }


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
    rows = minicram_rows(csv_path, out_prefix, project)
    spec: dict = {
        "computeResource": {
            "cpuMilli": 4000,
            "memoryMib": 16384,
            "bootDiskMib": 51200,
        },
        "maxRetryCount": 0,
        "maxRunDuration": "14400s",
        "environment": {"variables": {}},
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
    group: dict = {"taskSpec": spec}
    _attach_rows(group, spec, rows)
    return _job(group, project, region, sa)


def build_genotype_job(
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
    rows = genotype_rows(csv_path, out_prefix, project)
    spec: dict = {
        "computeResource": {
            "cpuMilli": 4000,
            "memoryMib": 8192,
            "bootDiskMib": 51200,
        },
        "maxRetryCount": 0,
        "maxRunDuration": "14400s",
        "environment": {"variables": {}},
        "runnables": [
            {"script": {"text": DOWNLOAD_SCRIPT}},
            {
                "container": {
                    "imageUri": image,
                    "entrypoint": "python3",
                    "commands": ["-c", worker],
                    "volumes": [f"{WORKDIR}:/work"],
                }
            },
            {"script": {"text": UPLOAD_SCRIPT}, "alwaysRun": True},
        ],
    }
    group: dict = {"taskSpec": spec}
    _attach_rows(group, spec, rows)
    return _job(group, project, region, sa)


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
    p.add_argument("--stage", choices=("minicram", "genotype"), default="minicram")
    args = p.parse_args()
    builder = build_minicram_job if args.stage == "minicram" else build_genotype_job
    job = builder(
        csv_path=args.csv,
        worker_path=args.worker,
        image=args.image,
        out_prefix=args.out_prefix,
        project=args.project,
        region=args.region,
        sa=args.sa,
    )
    args.config.write_text(json.dumps(job, indent=2) + "\n", encoding="utf-8")
    n = len(job["taskGroups"][0].get("taskEnvironments") or [1])
    print(f"wrote {args.config} ({args.stage}, {n} task(s))")


if __name__ == "__main__":
    main()
