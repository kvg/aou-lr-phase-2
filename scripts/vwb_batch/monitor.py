"""Read-only status for native Google Batch jobs and GCS run ledgers.

Does not submit. PET cannot read Cloud Logging; worker logs are GCS objects.
"""

from __future__ import annotations

import re
from collections import Counter
from datetime import datetime, timezone
from typing import Any

from . import poll, submit

ACTIVE = frozenset({"QUEUED", "SCHEDULED", "RUNNING"})
DISPLAY_COLUMNS = (
    "job_id",
    "state",
    "age",
    "pipeline",
    "stage",
    "run_id",
    "shard",
    "tasks",
    "ok",
    "fail",
    "run",
    "vcpu",
    "machine",
)


def job_id_of(job: dict[str, Any]) -> str:
    name = str(job.get("name") or "")
    if "/jobs/" in name:
        return name.rsplit("/jobs/", 1)[-1]
    return str(job.get("uid") or name)


def job_state(job: dict[str, Any]) -> str:
    status = job.get("status") or {}
    return str(status.get("state") or job.get("state") or "UNKNOWN")


def labels_of(job: dict[str, Any]) -> dict[str, str]:
    raw = job.get("labels") or {}
    return {str(k): str(v) for k, v in raw.items()}


def task_counts(job: dict[str, Any]) -> dict[str, int]:
    status = job.get("status") or {}
    groups = status.get("taskGroups") or {}
    if isinstance(groups, dict):
        iterable = groups.values()
    elif isinstance(groups, list):
        iterable = groups
    else:
        iterable = []
    out: dict[str, int] = {}
    for group in iterable:
        if not isinstance(group, dict):
            continue
        counts = group.get("counts") or group.get("countsByStatus") or {}
        for key, value in counts.items():
            try:
                n = int(value)
            except (TypeError, ValueError):
                continue
            out[str(key).upper()] = out.get(str(key).upper(), 0) + n
    return out


def _first_task_group(job: dict[str, Any]) -> dict[str, Any]:
    groups = job.get("taskGroups") or []
    if isinstance(groups, list) and groups:
        return groups[0] if isinstance(groups[0], dict) else {}
    if isinstance(groups, dict):
        for value in groups.values():
            if isinstance(value, dict):
                return value
    return {}


def requested_tasks(job: dict[str, Any]) -> int:
    group = _first_task_group(job)
    try:
        return int(group.get("taskCount") or 0)
    except (TypeError, ValueError):
        return 0


def parallelism_of(job: dict[str, Any]) -> int:
    group = _first_task_group(job)
    try:
        return int(group.get("parallelism") or 0)
    except (TypeError, ValueError):
        return 0


def machine_type(job: dict[str, Any]) -> str:
    alloc = job.get("allocationPolicy") or {}
    instances = alloc.get("instances") or []
    if instances and isinstance(instances[0], dict):
        policy = instances[0].get("policy") or instances[0]
        mt = str(policy.get("machineType") or "").strip()
        if mt:
            return mt
    group = _first_task_group(job)
    cr = (group.get("taskSpec") or {}).get("computeResource") or {}
    try:
        milli = int(cr.get("cpuMilli") or 0)
    except (TypeError, ValueError):
        milli = 0
    if milli:
        return f"cpuMilli={milli}"
    return ""


def vcpu_of(job: dict[str, Any]) -> int:
    mt = machine_type(job)
    parts = mt.replace("cpuMilli=", "").split("-")
    if parts and parts[-1].isdigit() and "cpuMilli" not in mt:
        return int(parts[-1])
    group = _first_task_group(job)
    cr = (group.get("taskSpec") or {}).get("computeResource") or {}
    try:
        milli = int(cr.get("cpuMilli") or 0)
    except (TypeError, ValueError):
        milli = 0
    if milli:
        return max(milli // 1000, 1)
    return 0


def parse_run_duration(text: str | None) -> float | None:
    if not text:
        return None
    raw = str(text).strip()
    if raw.endswith("s"):
        raw = raw[:-1]
    try:
        return float(raw)
    except ValueError:
        return None


def age_seconds(job: dict[str, Any], *, now: datetime | None = None) -> float | None:
    status = job.get("status") or {}
    duration = parse_run_duration(status.get("runDuration"))
    if duration is not None:
        return duration
    created = poll.parse_rfc3339(job.get("createTime") or job.get("create_time"))
    if created is None:
        return None
    clock = now or datetime.now(timezone.utc)
    if created.tzinfo is None:
        created = created.replace(tzinfo=timezone.utc)
    return max((clock - created).total_seconds(), 0.0)


def fmt_seconds(seconds: float | None) -> str:
    if seconds is None:
        return ""
    total = int(round(seconds))
    hours, rem = divmod(total, 3600)
    minutes, secs = divmod(rem, 60)
    if hours:
        return f"{hours}h{minutes:02d}m"
    if minutes:
        return f"{minutes}m{secs:02d}s"
    return f"{secs}s"


def inflight_tasks(row: dict[str, Any]) -> int:
    if row.get("state") not in ACTIVE:
        return 0
    requested = int(row.get("task_count") or 0)
    ok = int(row.get("succeeded") or 0)
    fail = int(row.get("failed") or 0)
    remaining = max(requested - ok - fail, 0)
    running = int(row.get("running") or 0)
    if remaining:
        par = int(row.get("parallelism") or remaining)
        return min(par, remaining) if par else remaining
    if running:
        return running
    return int(row.get("parallelism") or 1)


def inflight_vcpu(row: dict[str, Any]) -> int:
    n = inflight_tasks(row)
    vcpu = int(row.get("vcpu") or 0)
    if not n or not vcpu:
        return 0
    return n * vcpu


def summarize_job(job: dict[str, Any], *, now: datetime | None = None) -> dict[str, Any]:
    labels = labels_of(job)
    counts = task_counts(job)
    state = job_state(job)
    requested = requested_tasks(job)
    if not requested:
        requested = sum(counts.values())
    row = {
        "job_id": job_id_of(job),
        "state": state,
        "create_time": job.get("createTime") or job.get("create_time") or "",
        "age": fmt_seconds(age_seconds(job, now=now)),
        "age_seconds": age_seconds(job, now=now),
        "pipeline": labels.get("pipeline", ""),
        "stage": labels.get("stage", ""),
        "run_id": labels.get("run-id", ""),
        "shard": labels.get("shard", ""),
        "task_count": requested,
        "succeeded": counts.get("SUCCEEDED", 0),
        "failed": counts.get("FAILED", 0),
        "running": counts.get("RUNNING", 0) + counts.get("ASSIGNED", 0),
        "parallelism": parallelism_of(job),
        "vcpu": vcpu_of(job),
        "machine": machine_type(job),
        "labels": labels,
        "counts": counts,
    }
    row["tasks"] = str(requested) if requested else ""
    row["ok"] = str(row["succeeded"]) if counts else ""
    row["fail"] = str(row["failed"]) if counts else ""
    row["run"] = str(row["running"]) if counts else ""
    row["inflight_vcpu"] = inflight_vcpu(row)
    return row


def filter_jobs(
    jobs: list[dict[str, Any]],
    *,
    name_prefix: str = "",
    run_id: str = "",
    pipeline: str = "",
    active_only: bool = False,
) -> list[dict[str, Any]]:
    prefix = name_prefix.strip().lower()
    wanted_run = re.sub(r"[^a-z0-9_-]+", "-", run_id.lower()).strip("-") if run_id else ""
    wanted_pipe = re.sub(r"[^a-z0-9_-]+", "-", pipeline.lower()).strip("-") if pipeline else ""
    out: list[dict[str, Any]] = []
    for job in jobs:
        jid = job_id_of(job)
        if prefix and not jid.lower().startswith(prefix):
            continue
        labels = labels_of(job)
        if wanted_run and labels.get("run-id", "") != wanted_run:
            if wanted_run not in jid.lower():
                continue
        if wanted_pipe and labels.get("pipeline", "") != wanted_pipe:
            continue
        if active_only and job_state(job) not in ACTIVE:
            continue
        out.append(job)
    return out


def count_by_state(rows: list[dict[str, Any]]) -> dict[str, int]:
    return dict(Counter(str(r.get("state") or "UNKNOWN") for r in rows))


def total_inflight_vcpu(rows: list[dict[str, Any]]) -> int:
    return sum(int(r.get("inflight_vcpu") or 0) for r in rows)


def format_table(rows: list[dict[str, Any]], columns: tuple[str, ...] = DISPLAY_COLUMNS) -> str:
    if not rows:
        return "(no jobs)"
    widths = {col: len(col) for col in columns}
    rendered: list[dict[str, str]] = []
    for row in rows:
        line = {col: str(row.get(col, "") or "") for col in columns}
        rendered.append(line)
        for col, value in line.items():
            widths[col] = max(widths[col], len(value))
    header = "  ".join(col.ljust(widths[col]) for col in columns)
    sep = "  ".join("-" * widths[col] for col in columns)
    body = ["  ".join(line[col].ljust(widths[col]) for col in columns) for line in rendered]
    return "\n".join([header, sep, *body])


def list_summaries(
    *,
    project: str,
    region: str,
    name_prefix: str = "",
    run_id: str = "",
    pipeline: str = "",
    active_only: bool = False,
    enrich_active: bool = True,
) -> list[dict[str, Any]]:
    jobs = submit.list_jobs(project=project, region=region)
    jobs = filter_jobs(
        jobs,
        name_prefix=name_prefix,
        run_id=run_id,
        pipeline=pipeline,
        active_only=active_only,
    )
    rows: list[dict[str, Any]] = []
    for job in jobs:
        state = job_state(job)
        detail = job
        if enrich_active and state in ACTIVE:
            try:
                detail = submit.describe_job(
                    job_id=job_id_of(job),
                    project=project,
                    region=region,
                    quiet=True,
                )
            except RuntimeError as exc:
                row = summarize_job(job)
                row["error"] = str(exc)
                rows.append(row)
                continue
        rows.append(summarize_job(detail))
    rows.sort(key=lambda r: (r.get("create_time") or "", r.get("job_id") or ""), reverse=True)
    return rows


def task_breakdown(*, job_id: str, project: str, region: str) -> dict[str, Any]:
    job = submit.describe_job(job_id=job_id, project=project, region=region, quiet=True)
    tasks = submit.list_tasks(job_id=job_id, project=project, region=region, quiet=True)
    parsed = []
    for task in tasks:
        parsed.append(
            {
                "index": poll.task_index(task),
                "state": poll.task_state(task),
                "seconds": poll.task_seconds(task),
            }
        )
    parsed.sort(key=lambda row: row["index"])
    by_state = dict(Counter(t["state"] for t in parsed))
    return {
        "job": summarize_job(job),
        "by_state": by_state,
        "tasks": parsed,
    }
