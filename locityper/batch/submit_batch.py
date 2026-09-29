#!/usr/bin/env python3
"""Build native Google Batch job JSON for Locityper minicram / genotype.

Minicram: one print-reads container; CRAM stays gs://.
Genotype: host gcloud downloads the minicram + fasta + jellyfish + BED + db,
Locityper container runs preproc/genotype, host gcloud uploads tar/CSV.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path

WORKDIR = "/mnt/disks/lt"


def require_artifact_registry(image: str, *, role: str = "LOCITYPER_DOCKER") -> None:
    """VWB Batch VMs cannot pull Docker Hub under VPC-SC."""
    if "pkg.dev/" not in (image or ""):
        raise SystemExit(
            f"{role}={image!r} is not Artifact Registry. "
            "VPC-SC cannot pull Docker Hub (eichlerlab/locityper:1.4.5.0). "
            "Mirror the image into a workspace-readable AR repo and set LOCITYPER_DOCKER."
        )

LOAD_TASK_ENV = r"""
if [ -n "${TASKS_TSV:-}" ]; then
  python3 - <<'PY'
import csv, json, os, subprocess, sys
uri = os.environ["TASKS_TSV"]
idx = int(os.environ.get("BATCH_TASK_INDEX", "0"))
project = os.environ.get("GCLOUD_PROJECT") or os.environ.get("GOOGLE_CLOUD_PROJECT") or ""
workdir = "/mnt/disks/lt"
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
  . /mnt/disks/lt/task.env
fi
"""

DOWNLOAD_SCRIPT = r"""#!/bin/sh
set -eu
WORKDIR=/mnt/disks/lt
mkdir -p "${WORKDIR}"
chmod 777 "${WORKDIR}"
HOST_LOG="${WORKDIR}/host.log"
{
  echo "=== download start $(date -u +%Y-%m-%dT%H:%M:%SZ) ==="
  echo "MINICRAM=${MINICRAM:-}"
  echo "COUNTS_JF=${COUNTS_JF:-}"
  echo "DB_TAR=${DB_TAR:-}"
} >> "${HOST_LOG}"
""" + LOAD_TASK_ENV + r"""
gcloud storage cp --billing-project="${GCLOUD_PROJECT}" "${MINICRAM}" "${WORKDIR}/reads.cram"
gcloud storage cp --billing-project="${GCLOUD_PROJECT}" "${MINICRAI}" "${WORKDIR}/reads.cram.crai"
gcloud storage cp --billing-project="${GCLOUD_PROJECT}" "${REF_FA}" "${WORKDIR}/reference.fa"
gcloud storage cp --billing-project="${GCLOUD_PROJECT}" "${REF_FAI}" "${WORKDIR}/reference.fa.fai"
gcloud storage cp --billing-project="${GCLOUD_PROJECT}" "${COUNTS_JF}" "${WORKDIR}/counts.jf"
gcloud storage cp --billing-project="${GCLOUD_PROJECT}" "${BED}" "${WORKDIR}/loci.bed"
gcloud storage cp --billing-project="${GCLOUD_PROJECT}" "${DB_TAR}" "${WORKDIR}/vcf_db.tar.gz"
ls -lh "${WORKDIR}"
ls -lh "${WORKDIR}" >> "${WORKDIR}/host.log" || true
echo "=== download done $(date -u +%Y-%m-%dT%H:%M:%SZ) ===" >> "${WORKDIR}/host.log"
"""

UPLOAD_SCRIPT = r"""#!/bin/sh
set -eu
WORKDIR=/mnt/disks/lt
if [ -f "${WORKDIR}/task.env" ]; then
  # shellcheck disable=SC1091
  . "${WORKDIR}/task.env"
fi
HOST_LOG="${WORKDIR}/host.log"
{
  echo "=== upload start $(date -u +%Y-%m-%dT%H:%M:%SZ) ==="
  echo "SAMPLE_ID=${SAMPLE_ID:-}"
  ls -lh "${WORKDIR}" || true
} >> "${HOST_LOG}"
csv="${WORKDIR}/gts.filtered.csv"
tgz="${WORKDIR}/${SAMPLE_ID}.locityper.tar.gz"
log="${WORKDIR}/worker.log"
if [ -s "${csv}" ]; then
  gcloud storage cp --billing-project="${GCLOUD_PROJECT}" "${csv}" "${SUMMARY_CSV}"
fi
if [ -s "${tgz}" ]; then
  gcloud storage cp --billing-project="${GCLOUD_PROJECT}" "${tgz}" "${RESULTS_TAR}"
fi
if [ -s "${log}" ]; then
  gcloud storage cp --billing-project="${GCLOUD_PROJECT}" "${log}" "${WORKER_LOG_GCS}"
fi
if [ -n "${RESOURCE_STATS:-}" ] && [ -s "${WORKDIR}/resource_stats.tsv" ]; then
  gcloud storage cp --billing-project="${GCLOUD_PROJECT}" "${WORKDIR}/resource_stats.tsv" "${RESOURCE_STATS}"
fi
if [ -n "${HOST_LOG_GCS:-}" ] && [ -s "${HOST_LOG}" ]; then
  echo "=== upload done $(date -u +%Y-%m-%dT%H:%M:%SZ) ===" >> "${HOST_LOG}"
  gcloud storage cp --billing-project="${GCLOUD_PROJECT}" "${HOST_LOG}" "${HOST_LOG_GCS}"
fi
ls -lh "${csv}" "${tgz}" "${log}" "${WORKDIR}/resource_stats.tsv" "${HOST_LOG}" 2>/dev/null || true
"""

# dsub used 4 vCPU / 16 GiB minicram and 2 vCPU / (4*n+6) GiB genotype.
DEFAULT_COMPUTE = {
    "minicram": {"cpuMilli": 4000, "memoryMib": 16384, "bootDiskMib": 51200},
    "genotype": {"cpuMilli": 2000, "memoryMib": 14336, "bootDiskMib": 81920},
}


def _summarize_script() -> str:
    py = Path(__file__).with_name("genotype.py").read_text(encoding="utf-8")
    return (
        "#!/bin/sh\n"
        "set -eu\n"
        "WORKDIR=/mnt/disks/lt\n"
        'if [ -f "${WORKDIR}/task.env" ]; then\n'
        "  # shellcheck disable=SC1091\n"
        '  . "${WORKDIR}/task.env"\n'
        "fi\n"
        'export LT_WORKDIR="${WORKDIR}"\n'
        "python3 - <<'PY'\n"
        f"{py}"
        "PY\n"
    )


def _worker_source(worker_path: Path) -> str:
    monitor = Path(__file__).resolve().parents[2] / "expansion_hunter" / "batch" / "resource_monitor.py"
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


def _ref_fa(rec: dict) -> str:
    return rec.get("ref_fa") or rec.get("ref_fa_uncompressed") or ""


def _ref_fai(rec: dict) -> str:
    return rec.get("ref_fai") or rec.get("ref_fai_uncompressed") or ""


def _db_tar(rec: dict) -> str:
    return rec.get("db_tar") or rec.get("locityper_db_tar_gz") or ""


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
                "REF_FA": _ref_fa(rec),
                "REF_FAI": _ref_fai(rec),
                "BED": rec["bed"],
                "GCLOUD_PROJECT": _project_for_row(rec, project),
                "WORKER_LOG_GCS": _gs(out_prefix, sid, f"{sid}.minicram.worker.log"),
                "MINICRAM": _gs(out_prefix, sid, f"{sid}.minicram.cram"),
                "MINICRAI": _gs(out_prefix, sid, f"{sid}.minicram.cram.crai"),
                "TRANSFER_STATS": _gs(out_prefix, sid, f"{sid}.data_transfer_stats.tsv"),
                "RESOURCE_STATS": _gs(out_prefix, sid, f"{sid}.minicram.resources.tsv"),
                "WINDOW_GRAB": os.environ.get("WINDOW_GRAB", "3000"),
                "BG_REGION_BED": os.environ.get("BG_REGION_BED", "chr17\t72062001\t76562000"),
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
    n_cpu = os.environ.get("LOCITYPER_N_CPU", "2")
    for rec in records:
        sid = rec["sample_id"]
        rows.append(
            {
                "SAMPLE_ID": sid,
                "SEX": rec.get("sex") or "female",
                "MINICRAM": _gs(out_prefix, sid, f"{sid}.minicram.cram"),
                "MINICRAI": _gs(out_prefix, sid, f"{sid}.minicram.cram.crai"),
                "REF_FA": _ref_fa(rec),
                "REF_FAI": _ref_fai(rec),
                "COUNTS_JF": rec["counts_jf"],
                "BED": rec["bed"],
                "DB_TAR": _db_tar(rec),
                "GCLOUD_PROJECT": _project_for_row(rec, project),
                "WORKER_LOG_GCS": _gs(out_prefix, sid, f"{sid}.locityper.worker.log"),
                "HOST_LOG_GCS": _gs(out_prefix, sid, f"{sid}.locityper.host.log"),
                "SUMMARY_CSV": _gs(out_prefix, sid, f"{sid}.gts.filtered.csv"),
                "RESULTS_TAR": _gs(out_prefix, sid, f"{sid}.locityper.tar.gz"),
                "RESOURCE_STATS": _gs(out_prefix, sid, f"{sid}.locityper.resources.tsv"),
                "LOCITYPER_N_CPU": n_cpu,
                "TECHNOLOGY": os.environ.get("TECHNOLOGY", "illumina"),
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
            for key in (
                "WINDOW_GRAB",
                "BG_REGION_BED",
                "MAX_RETRY",
                "WAIT_TIME",
                "LOCITYPER_N_CPU",
                "TECHNOLOGY",
            ):
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
    del out_prefix
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
    if worker_path.suffix != ".sh":
        raise SystemExit(
            f"genotype worker must be genotype.sh (Locityper image has no python3), got {worker_path}"
        )
    shell = worker_path.read_text(encoding="utf-8")
    cr = _compute("genotype", compute)
    spec: dict = {
        "computeResource": cr,
        "maxRetryCount": 0,
        "maxRunDuration": max_run_duration,
        "environment": {
            "variables": {"LOCITYPER_N_CPU": str(max(int(cr["cpuMilli"]) // 1000, 1))}
        },
        "runnables": [
            {"script": {"text": DOWNLOAD_SCRIPT}},
            {
                "container": {
                    "imageUri": image,
                    "entrypoint": "/bin/bash",
                    "commands": ["-c", shell],
                    "volumes": [f"{WORKDIR}:/work"],
                }
            },
            {"script": {"text": _summarize_script()}},
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
