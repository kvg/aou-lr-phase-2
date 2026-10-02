"""Choose the next ExpansionHunter submit wave from the run ledger.

First attempts (never submitted, or stranded before the task started) go out
before any sample whose last task actually failed. Samples with expected
objects, and samples in a job that is still queued or running, are skipped.
"""

from __future__ import annotations

from typing import Any

from .poll import TERMINAL

CATEGORIES = ("done", "running", "failed", "not_started")

# Terminal job, but the task never started. Batch drops these when a sibling fails.
_NEVER_RAN = frozenset(
    {
        "",
        "UNKNOWN",
        "STATE_UNSPECIFIED",
        "PENDING",
        "QUEUED",
        "SCHEDULED",
        "UNEXECUTED",
    }
)


def latest_submissions(events: list[dict[str, Any]], stage: str) -> dict[str, tuple[str, int]]:
    """sample_id -> (job_id, task_index) for the latest job_submitted of this stage."""
    out: dict[str, tuple[str, int]] = {}
    for ev in events:
        if ev.get("event") != "job_submitted" or ev.get("stage") != stage:
            continue
        job_id = ev.get("job_id")
        if not job_id:
            continue
        for idx, sample_id in enumerate(ev.get("sample_ids") or []):
            out[str(sample_id)] = (str(job_id), idx)
    return out


def _task_state(polled: dict[str, Any], index: int) -> str:
    for task in polled.get("tasks") or []:
        try:
            task_index = int(task.get("index"))
        except (TypeError, ValueError):
            continue
        if task_index == index:
            return str(task.get("state") or "")
    return ""


def classify_stage(
    sample_ids: list[str],
    *,
    submissions: dict[str, tuple[str, int]],
    polls: dict[str, dict[str, Any]],
    done_ids: set[str],
) -> dict[str, str]:
    """Map each sample to done, running, failed, or not_started.

    ``polls`` is job_id -> {"state", "tasks": [{"index", "state"}]}. A missing
    poll for a submitted job counts as running so a later wave does not submit
    a duplicate.
    """
    out: dict[str, str] = {}
    for sample_id in sample_ids:
        if sample_id in done_ids:
            out[sample_id] = "done"
            continue
        submission = submissions.get(sample_id)
        if submission is None:
            out[sample_id] = "not_started"
            continue
        job_id, index = submission
        polled = polls.get(job_id)
        if polled is None:
            out[sample_id] = "running"
            continue
        state = str(polled.get("state") or "UNKNOWN")
        if state not in TERMINAL or state == "DELETION_IN_PROGRESS":
            out[sample_id] = "running"
            continue
        task_state = _task_state(polled, index)
        if task_state in _NEVER_RAN:
            out[sample_id] = "not_started"
            continue
        out[sample_id] = "failed"
    return out


def count_categories(categories: dict[str, str]) -> dict[str, int]:
    counts = {key: 0 for key in CATEGORIES}
    for value in categories.values():
        if value in counts:
            counts[value] += 1
    return counts


def choose_to_submit(sample_ids: list[str], categories: dict[str, str]) -> list[str]:
    """First attempts in list order, or retries once none of those remain."""
    first = [sid for sid in sample_ids if categories.get(sid) == "not_started"]
    if first:
        return first
    return [sid for sid in sample_ids if categories.get(sid) == "failed"]
