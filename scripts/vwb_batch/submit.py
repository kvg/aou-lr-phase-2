"""Submit / describe / list native Google Batch jobs via gcloud."""

from __future__ import annotations

import json
import subprocess
import tempfile
from pathlib import Path
from typing import Any

MAX_JOB_JSON_BYTES = 900_000  # Batch limit is 1 MiB; leave headroom


def _run(cmd: list[str], *, check: bool = True, quiet: bool = False) -> subprocess.CompletedProcess[str]:
    if not quiet:
        print("$", " ".join(cmd), flush=True)
    proc = subprocess.run(cmd, check=False, text=True, capture_output=True)
    if not quiet:
        if proc.stdout:
            print(proc.stdout, end="" if proc.stdout.endswith("\n") else "\n")
        if proc.stderr:
            print(proc.stderr, end="" if proc.stderr.endswith("\n") else "\n")
    if check and proc.returncode != 0:
        err = (proc.stderr or proc.stdout or "").strip()
        raise RuntimeError(f"command failed ({proc.returncode}): {' '.join(cmd)}" + (f"\n{err}" if err else ""))
    return proc


def job_json_size(job: dict[str, Any]) -> int:
    return len(json.dumps(job, separators=(",", ":")).encode("utf-8"))


def assert_job_json_ok(job: dict[str, Any]) -> int:
    size = job_json_size(job)
    if size > MAX_JOB_JSON_BYTES:
        raise SystemExit(
            f"job JSON is {size} bytes (limit {MAX_JOB_JSON_BYTES}). "
            "Use TASKS_TSV + BATCH_TASK_INDEX instead of taskEnvironments."
        )
    return size


def submit_job(
    *,
    job_id: str,
    job: dict[str, Any],
    project: str,
    region: str,
) -> str:
    size = assert_job_json_ok(job)
    print(f"job JSON {size} bytes, id={job_id}", flush=True)
    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as fh:
        json.dump(job, fh, indent=2)
        fh.write("\n")
        cfg = Path(fh.name)
    try:
        _run(
            [
                "gcloud",
                "batch",
                "jobs",
                "submit",
                job_id,
                f"--location={region}",
                f"--project={project}",
                f"--config={cfg}",
            ]
        )
    finally:
        cfg.unlink(missing_ok=True)
    return job_id


def describe_job(*, job_id: str, project: str, region: str, quiet: bool = False) -> dict[str, Any]:
    proc = _run(
        [
            "gcloud",
            "batch",
            "jobs",
            "describe",
            job_id,
            f"--location={region}",
            f"--project={project}",
            "--format=json",
        ],
        check=False,
        quiet=quiet,
    )
    if proc.returncode != 0:
        err = (proc.stderr or proc.stdout or "").strip()
        raise RuntimeError(f"describe {job_id} failed: {err}")
    return json.loads(proc.stdout)


def list_jobs(*, project: str, region: str, filter_expr: str = "") -> list[dict[str, Any]]:
    cmd = [
        "gcloud",
        "batch",
        "jobs",
        "list",
        f"--location={region}",
        f"--project={project}",
        "--format=json",
    ]
    if filter_expr:
        cmd.append(f"--filter={filter_expr}")
    proc = _run(cmd, check=False, quiet=True)
    if proc.returncode != 0:
        err = (proc.stderr or proc.stdout or "").strip()
        raise RuntimeError(f"jobs list failed: {err}")
    data = json.loads(proc.stdout or "[]")
    if isinstance(data, dict):
        return data.get("jobs") or []
    return data


def list_tasks(*, job_id: str, project: str, region: str, quiet: bool = False) -> list[dict[str, Any]]:
    proc = _run(
        [
            "gcloud",
            "batch",
            "tasks",
            "list",
            f"--job={job_id}",
            f"--location={region}",
            f"--project={project}",
            "--format=json",
        ],
        check=False,
        quiet=quiet,
    )
    if proc.returncode != 0:
        err = (proc.stderr or proc.stdout or "").strip()
        raise RuntimeError(f"tasks list {job_id} failed: {err}")
    data = json.loads(proc.stdout or "[]")
    if isinstance(data, dict):
        return data.get("tasks") or []
    return data
