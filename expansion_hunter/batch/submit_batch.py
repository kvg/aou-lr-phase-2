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

# Host Batch scripts: if TASKS_TSV is set, pull row BATCH_TASK_INDEX into
# /mnt/disks/eh/task.env + task.json so later runnables see SAMPLE_ID etc.
LOAD_TASK_ENV = r"""
if [ -n "${TASKS_TSV:-}" ]; then
  python3 - <<'PY'
import csv, json, os, subprocess, sys
uri = os.environ["TASKS_TSV"]
idx = int(os.environ.get("BATCH_TASK_INDEX", "0"))
project = os.environ.get("GCLOUD_PROJECT") or os.environ.get("GOOGLE_CLOUD_PROJECT") or ""
workdir = "/mnt/disks/eh"
os.makedirs(workdir, exist_ok=True)
dest = os.path.join(workdir, "tasks.tsv")
cmd = ["gcloud", "storage", "cp"]
if project:
    cmd.append("--billing-project=" + project)
cmd.extend([uri, dest])
print("+", " ".join(cmd), flush=True)
subprocess.check_call(cmd)
with open(dest, newline="", encoding="utf-8") as fh:
    rows = list(csv.DictReader(fh, delimiter="\t"))
if idx < 0 or idx >= len(rows):
    sys.exit(f"BATCH_TASK_INDEX {idx} out of range n={len(rows)}")
row = rows[idx]
with open(os.path.join(workdir, "task.json"), "w", encoding="utf-8") as fh:
    json.dump(row, fh)
with open(os.path.join(workdir, "task.env"), "w", encoding="utf-8") as fh:
    for key, value in row.items():
        if not key:
            continue
        fh.write(f"export {key}={json.dumps('' if value is None else str(value))}\n")
print(f"TASKS_TSV row {idx} SAMPLE_ID={row.get('SAMPLE_ID')}", flush=True)
PY
  # shellcheck disable=SC1091
  . /mnt/disks/eh/task.env
fi
"""

DOWNLOAD_SCRIPT = r"""#!/bin/bash
set -eu
WORKDIR=/mnt/disks/eh
mkdir -p "${WORKDIR}"
chmod 777 "${WORKDIR}"
""" + LOAD_TASK_ENV + r"""
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
if [ -f "${WORKDIR}/task.env" ]; then
  # shellcheck disable=SC1091
  . "${WORKDIR}/task.env"
fi
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
if [ -n "${RESOURCE_STATS:-}" ] && [ -s "${WORKDIR}/resource_stats.tsv" ]; then
  gcloud storage cp --billing-project="${GCLOUD_PROJECT}" "${WORKDIR}/resource_stats.tsv" "${RESOURCE_STATS}"
fi
ls -lh "${json}" "${vcf}" "${log}" "${WORKDIR}/resource_stats.tsv" 2>/dev/null || true
"""

DEFAULT_COMPUTE = {
    "minicram": {"cpuMilli": 4000, "memoryMib": 16384, "bootDiskMib": 51200},
    "genotype": {"cpuMilli": 4000, "memoryMib": 8192, "bootDiskMib": 51200},
}


def _worker_source(worker_path: Path) -> str:
    monitor = Path(__file__).resolve().parent / "resource_monitor.py"
    worker = worker_path.read_text(encoding="utf-8")
    worker = worker.replace("from __future__ import annotations\n", "", 1)
    return monitor.read_text(encoding="utf-8") + "\n" + worker


def _compute(stage: str, override: dict | None) -> dict:
    spec = dict(DEFAULT_COMPUTE[stage])
    if override:
        spec.update({k: int(override[k]) for k in spec if k in override})
    return spec


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


def minicram_rows_from_records(
    records: list[dict[str, str]], out_prefix: str, project: str
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for rec in records:
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
                "RESOURCE_STATS": _gs(out_prefix, sid, f"{sid}.minicram.resources.tsv"),
                "WINDOW_SIZE": os.environ.get("WINDOW_SIZE", "1000"),
                "MERGE_REGIONS_DISTANCE": os.environ.get("MERGE_REGIONS_DISTANCE", "1000"),
                "MAX_RETRY": os.environ.get("MAX_RETRY", "3"),
                "WAIT_TIME": os.environ.get("WAIT_TIME", "30"),
            }
        )
    if not rows:
        raise SystemExit("no minicram records")
    return rows


def genotype_rows_from_records(
    records: list[dict[str, str]], out_prefix: str, project: str
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for rec in records:
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
                "RESOURCE_STATS": _gs(out_prefix, sid, f"{sid}.EH.resources.tsv"),
            }
        )
    if not rows:
        raise SystemExit("no genotype records")
    return rows


def minicram_rows(csv_path: Path, out_prefix: str, project: str) -> list[dict[str, str]]:
    with csv_path.open(newline="") as fh:
        records = list(csv.DictReader(fh))
    rows = minicram_rows_from_records(records, out_prefix, project)
    if not rows:
        raise SystemExit(f"no rows in {csv_path}")
    return rows


def genotype_rows(csv_path: Path, out_prefix: str, project: str) -> list[dict[str, str]]:
    with csv_path.open(newline="") as fh:
        records = list(csv.DictReader(fh))
    rows = genotype_rows_from_records(records, out_prefix, project)
    if not rows:
        raise SystemExit(f"no rows in {csv_path}")
    return rows


def _attach_rows(
    group: dict,
    spec: dict,
    rows: list[dict[str, str]],
    *,
    tasks_tsv_uri: str | None = None,
    parallelism: int | None = None,
) -> None:
    group["taskCount"] = len(rows)
    if parallelism is not None:
        group["parallelism"] = max(1, min(int(parallelism), len(rows)))
    if tasks_tsv_uri:
        spec["environment"]["variables"]["TASKS_TSV"] = tasks_tsv_uri
        if rows:
            spec["environment"]["variables"].setdefault("GCLOUD_PROJECT", rows[0]["GCLOUD_PROJECT"])
            for key in ("WINDOW_SIZE", "MERGE_REGIONS_DISTANCE", "MAX_RETRY", "WAIT_TIME"):
                if key in rows[0]:
                    spec["environment"]["variables"].setdefault(key, rows[0][key])
        return
    if len(rows) == 1:
        spec["environment"]["variables"].update(rows[0])
    else:
        group["taskEnvironments"] = [{"variables": row} for row in rows]


def _job(group: dict, project: str, region: str, sa: str, labels: dict[str, str] | None = None) -> dict:
    job = {
        "taskGroups": [group],
        "allocationPolicy": _allocation(project, region, sa),
        "logsPolicy": {"destination": "CLOUD_LOGGING"},
    }
    if labels:
        job["labels"] = labels
    return job


def build_minicram_job_from_rows(
    *,
    rows: list[dict[str, str]],
    worker_path: Path,
    image: str,
    out_prefix: str,
    project: str,
    region: str,
    sa: str,
    max_run_duration: str = "14400s",
    tasks_tsv_uri: str | None = None,
    parallelism: int | None = None,
    labels: dict[str, str] | None = None,
    compute: dict | None = None,
) -> dict:
    del out_prefix  # rows already contain output URIs
    worker = _worker_source(worker_path)
    spec: dict = {
        "computeResource": _compute("minicram", compute),
        "maxRetryCount": 0,
        "maxRunDuration": max_run_duration,
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
    _attach_rows(group, spec, rows, tasks_tsv_uri=tasks_tsv_uri, parallelism=parallelism)
    return _job(group, project, region, sa, labels=labels)


def build_genotype_job_from_rows(
    *,
    rows: list[dict[str, str]],
    worker_path: Path,
    image: str,
    out_prefix: str,
    project: str,
    region: str,
    sa: str,
    max_run_duration: str = "14400s",
    tasks_tsv_uri: str | None = None,
    parallelism: int | None = None,
    labels: dict[str, str] | None = None,
    compute: dict | None = None,
) -> dict:
    del out_prefix
    worker = _worker_source(worker_path)
    spec: dict = {
        "computeResource": _compute("genotype", compute),
        "maxRetryCount": 0,
        "maxRunDuration": max_run_duration,
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
    _attach_rows(group, spec, rows, tasks_tsv_uri=tasks_tsv_uri, parallelism=parallelism)
    return _job(group, project, region, sa, labels=labels)


def build_minicram_job(
    *,
    csv_path: Path,
    worker_path: Path,
    image: str,
    out_prefix: str,
    project: str,
    region: str,
    sa: str,
    max_run_duration: str = "14400s",
    tasks_tsv_uri: str | None = None,
    parallelism: int | None = None,
    labels: dict[str, str] | None = None,
    compute: dict | None = None,
) -> dict:
    rows = minicram_rows(csv_path, out_prefix, project)
    return build_minicram_job_from_rows(
        rows=rows,
        worker_path=worker_path,
        image=image,
        out_prefix=out_prefix,
        project=project,
        region=region,
        sa=sa,
        max_run_duration=max_run_duration,
        tasks_tsv_uri=tasks_tsv_uri,
        parallelism=parallelism,
        labels=labels,
        compute=compute,
    )


def build_genotype_job(
    *,
    csv_path: Path,
    worker_path: Path,
    image: str,
    out_prefix: str,
    project: str,
    region: str,
    sa: str,
    max_run_duration: str = "14400s",
    tasks_tsv_uri: str | None = None,
    parallelism: int | None = None,
    labels: dict[str, str] | None = None,
    compute: dict | None = None,
) -> dict:
    rows = genotype_rows(csv_path, out_prefix, project)
    return build_genotype_job_from_rows(
        rows=rows,
        worker_path=worker_path,
        image=image,
        out_prefix=out_prefix,
        project=project,
        region=region,
        sa=sa,
        max_run_duration=max_run_duration,
        tasks_tsv_uri=tasks_tsv_uri,
        parallelism=parallelism,
        labels=labels,
        compute=compute,
    )


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
    n = job["taskGroups"][0].get("taskCount") or len(
        job["taskGroups"][0].get("taskEnvironments") or [1]
    )
    print(f"wrote {args.config} ({args.stage}, {n} task(s))")


if __name__ == "__main__":
    main()
