"""Poll Batch jobs and reconcile per-sample success against GCS objects."""

from __future__ import annotations

from datetime import datetime, timezone
from typing import Any, Callable

from . import gcs, submit

TERMINAL = frozenset({"SUCCEEDED", "FAILED", "DELETION_IN_PROGRESS", "CANCELLED"})


def parse_rfc3339(value: str | None) -> datetime | None:
    if not value:
        return None
    text = value.strip()
    if text.endswith("Z"):
        text = text[:-1] + "+00:00"
    try:
        return datetime.fromisoformat(text)
    except ValueError:
        return None


def duration_seconds(start: str | None, end: str | None) -> float | None:
    a = parse_rfc3339(start)
    b = parse_rfc3339(end)
    if a is None or b is None:
        return None
    if a.tzinfo is None:
        a = a.replace(tzinfo=timezone.utc)
    if b.tzinfo is None:
        b = b.replace(tzinfo=timezone.utc)
    return max((b - a).total_seconds(), 0.0)


def task_index(task: dict[str, Any]) -> int:
    name = str(task.get("name") or "")
    if "/tasks/" in name:
        return int(name.rsplit("/tasks/", 1)[-1])
    return int(task.get("task", 0) or 0)


def task_state(task: dict[str, Any]) -> str:
    status = task.get("status") or {}
    return str(status.get("state") or task.get("state") or "UNKNOWN")


def _event_time(task: dict[str, Any], wanted: str) -> str | None:
    status = task.get("status") or {}
    for event in status.get("statusEvents") or []:
        kind = str(event.get("type") or event.get("description") or "").upper()
        if wanted in kind:
            return event.get("eventTime") or event.get("event_time")
    return None


def task_seconds(task: dict[str, Any]) -> float | None:
    start = _event_time(task, "RUNNING") or _event_time(task, "ASSIGNED")
    end = _event_time(task, "SUCCEEDED") or _event_time(task, "FAILED")
    return duration_seconds(start, end)


def objects_ok(uris: list[str], listing: set[str] | None = None) -> bool:
    if listing is None:
        return all(gcs.exists(uri) for uri in uris)
    return all(uri in listing for uri in uris)


def poll_job(*, job_id: str, project: str, region: str) -> dict[str, Any]:
    job = submit.describe_job(job_id=job_id, project=project, region=region)
    status = job.get("status") or {}
    state = str(status.get("state") or "UNKNOWN")
    tasks = submit.list_tasks(job_id=job_id, project=project, region=region)
    parsed = []
    for task in tasks:
        parsed.append(
            {
                "index": task_index(task),
                "state": task_state(task),
                "seconds": task_seconds(task),
                "raw": task,
            }
        )
    parsed.sort(key=lambda row: row["index"])
    return {"job_id": job_id, "state": state, "job": job, "tasks": parsed}


def reconcile_shard(
    *,
    sample_ids: list[str],
    stage: str,
    job_id: str,
    shard: int,
    poll: dict[str, Any],
    expected: Callable[[str, str], list[str]],
    now: str,
    listing: set[str] | None = None,
) -> list[dict[str, Any]]:
    by_index = {t["index"]: t for t in poll["tasks"]}
    rows: list[dict[str, Any]] = []
    for idx, sample_id in enumerate(sample_ids):
        task = by_index.get(idx)
        state = task["state"] if task else poll["state"]
        uris = expected(sample_id, stage)
        ok = objects_ok(uris, listing) if uris else False
        if state == "SUCCEEDED" and not ok:
            state = "MISSING_OUTPUTS"
        elif ok and state != "SUCCEEDED":
            # Objects exist from a previous attempt; treat as done.
            state = "SUCCEEDED"
        rows.append(
            {
                "sample_id": sample_id,
                "stage": stage,
                "state": "SUCCEEDED" if ok else state,
                "objects_ok": ok,
                "job_id": job_id,
                "shard": shard,
                "task_index": idx,
                "seconds": (task or {}).get("seconds") or "",
                "network_bytes": "",
                "updated_at": now,
            }
        )
    return rows


def needs_work(
    status_rows: list[dict[str, Any]],
    *,
    stage: str,
    force: bool = False,
) -> list[str]:
    out: list[str] = []
    for row in status_rows:
        if row.get("stage") != stage:
            continue
        ok = row.get("objects_ok") in {True, "True", "true", "1"}
        if ok and not force:
            continue
        out.append(str(row["sample_id"]))
    return out


def pending_records(
    records: list[dict[str, str]],
    *,
    stage: str,
    out_prefix: str,
    expected,
    force: bool = False,
    listing: set[str] | None = None,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    pending: list[dict[str, str]] = []
    done: list[dict[str, str]] = []
    for rec in records:
        uris = expected(rec["sample_id"], stage, out_prefix)
        ok = objects_ok(uris, listing)
        (done if ok else pending).append(rec)
    if force:
        return list(records), []
    return pending, done


def inflight_sample_ids(events: list[dict[str, Any]], stage: str) -> set[str]:
    submitted: dict[str, set[str]] = {}
    finished: set[str] = set()
    for ev in events:
        if ev.get("stage") != stage:
            continue
        job_id = ev.get("job_id")
        if not job_id:
            continue
        if ev.get("event") == "job_submitted":
            submitted[str(job_id)] = set(ev.get("sample_ids") or [])
        if ev.get("event") == "job_polled" and ev.get("state") in TERMINAL:
            finished.add(str(job_id))
    out: set[str] = set()
    for job_id, sample_ids in submitted.items():
        if job_id not in finished:
            out.update(str(s) for s in sample_ids)
    return out
