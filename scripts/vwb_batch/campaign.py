"""Run a large cohort through a one-task-per-sample Batch pipeline in waves.

Each run of the submit notebook calls ``sync`` (one ``gcloud batch jobs list``
plus a task listing for each job that finished since last time), then ``plan``
(what to submit next, or why not), then ``submit``. Every rule that protects
the budget lives in ``plan``:

- Canary: the first submission for a new ``code_version`` is ``canary_n``
  samples, and nothing else goes out until that canary job has finished.
- Failure stop: if more than ``max_failure_rate`` of recent first attempts on
  the current code failed, submission halts. Pushing a code fix changes the
  code version, which starts a new canary.
- Budget: estimated spend so far, plus in-flight samples, plus the new wave,
  must stay under ``budget_usd``.
- Window: at most ``wave_size`` samples in flight. A new wave goes out once
  in-flight drops below ``refill_fraction`` of that.
- Retries: failed samples are retried only after every sample has had a first
  attempt, and at most ``max_attempts`` times in total. Spot preemptions are
  retried by Batch itself and do not count.

State lives next to the run ledger in GCS:

  run.json                 settings (see ``Settings``)
  samples.csv              samples in submission order
  jobs.json                every job: samples, state, per-task outcome, cost estimate
  jobs/{job_id}.tasks.json task timings and VMs, written once when the job finishes
  tasks/w{wave}/{job}.tsv  the TASKS_TSV each job reads
  ledger.jsonl             a few events per wave
  lock                     held by the notebook that is submitting
"""

from __future__ import annotations

import csv
import io
import json
import math
import re
import socket
import time
from dataclasses import asdict, dataclass, field, fields
from datetime import datetime, timezone
from typing import Any, Callable

from . import alloc, cost, gcs, monitor, poll, registry, submit
from .registry import RunPaths

TERMINAL = frozenset({"SUCCEEDED", "FAILED", "CANCELLED", "LOST"})
FIRST_ATTEMPT_KINDS = ("canary", "wave")
NEVER_RAN = frozenset({"", "UNKNOWN", "STATE_UNSPECIFIED", "PENDING", "QUEUED", "SCHEDULED", "ASSIGNED", "UNEXECUTED"})
ACTIVE_POINTER = "ACTIVE_RUN"
SUBMITTING_GRACE_S = 900


@dataclass
class Settings:
    out_prefix: str
    budget_usd: float = 3000.0
    wave_size: int = 10_000
    refill_fraction: float = 0.5
    canary_n: int = 200
    job_size: int = 1000
    max_attempts: int = 3
    max_failure_rate: float = 0.05
    min_finished_for_rate: int = 50
    failure_window: int = 1000
    prior_unit_usd: float = cost.PRIOR_UNIT_USD
    spot: bool = True
    machine_type: str = ""
    compute: dict = field(default_factory=lambda: {"cpuMilli": 1000, "memoryMib": 2048, "bootDiskMib": 20480})
    region: str = "us-central1"
    max_run_duration: str = "10800s"
    print_reads_image: str = ""
    eh_image: str = ""
    catalog: str = ""
    ref_fa: str = ""
    ref_fai: str = ""
    billing_table: str = ""

    @classmethod
    def from_run_json(cls, meta: dict[str, Any]) -> "Settings":
        raw = meta.get("campaign") or {}
        known = {f.name for f in fields(cls)}
        return cls(**{k: v for k, v in raw.items() if k in known})

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


# --------------------------------------------------------------------------
# Paths and state I/O
# --------------------------------------------------------------------------


def active_pointer_uri(output_bucket: str, pipeline: str) -> str:
    return f"{output_bucket.rstrip('/')}/batchRuns/{pipeline}/runs/{ACTIVE_POINTER}"


def read_active_run(output_bucket: str, pipeline: str) -> str:
    uri = active_pointer_uri(output_bucket, pipeline)
    if not gcs.exists(uri):
        return ""
    return gcs.cat(uri).strip()


def set_active_run(output_bucket: str, pipeline: str, run_id: str) -> None:
    gcs.upload_text(active_pointer_uri(output_bucket, pipeline), run_id + "\n")


def samples_uri(paths: RunPaths) -> str:
    return f"{paths.prefix}/samples.csv"


def jobs_uri(paths: RunPaths) -> str:
    return f"{paths.prefix}/jobs.json"


def lock_uri(paths: RunPaths) -> str:
    return f"{paths.prefix}/lock"


def task_detail_uri(paths: RunPaths, job_id: str) -> str:
    return f"{paths.prefix}/jobs/{job_id}.tasks.json"


def tasks_tsv_uri(paths: RunPaths, wave: int, job_id: str) -> str:
    return f"{paths.prefix}/tasks/w{wave:03d}/{job_id}.tsv"


SAMPLE_FIELDS = ("sample_id", "cram", "crai", "sex", "sex_raw")


def write_samples(paths: RunPaths, records: list[dict[str, str]]) -> None:
    buf = io.StringIO()
    writer = csv.DictWriter(buf, fieldnames=list(SAMPLE_FIELDS), extrasaction="ignore")
    writer.writeheader()
    for rec in records:
        writer.writerow({k: rec.get(k, "") for k in SAMPLE_FIELDS})
    gcs.upload_text(samples_uri(paths), buf.getvalue())


def load_samples(paths: RunPaths, settings: Settings) -> list[dict[str, str]]:
    text = gcs.cat(samples_uri(paths))
    out = []
    for rec in csv.DictReader(io.StringIO(text)):
        rec = {k: (v or "").strip() for k, v in rec.items()}
        rec.setdefault("ref_fa", settings.ref_fa)
        rec.setdefault("ref_fai", settings.ref_fai)
        rec.setdefault("catalog", settings.catalog)
        out.append(rec)
    return out


def empty_state() -> dict[str, Any]:
    return {"version": 1, "jobs": {}}


def load_state(paths: RunPaths) -> dict[str, Any]:
    uri = jobs_uri(paths)
    if not gcs.exists(uri):
        return empty_state()
    return json.loads(gcs.cat(uri))


def save_state(paths: RunPaths, state: dict[str, Any]) -> None:
    state["saved_at"] = registry.utcnow()
    gcs.upload_text(jobs_uri(paths), json.dumps(state, separators=(",", ":")) + "\n")


# --------------------------------------------------------------------------
# Lock
# --------------------------------------------------------------------------


class LockHeld(RuntimeError):
    pass


def acquire_lock(paths: RunPaths, *, stale_after_s: float = 1800.0) -> str:
    """Take the run lock. Returns a holder token. A lock older than stale_after_s is broken."""
    holder = f"{socket.gethostname()}:{time.time():.0f}"
    body = json.dumps({"holder": holder, "ts": time.time(), "at": registry.utcnow()})
    for _ in range(2):
        if gcs.upload_text_if_absent(lock_uri(paths), body):
            return holder
        try:
            cur = json.loads(gcs.cat(lock_uri(paths)))
        except Exception:
            cur = {}
        age = time.time() - float(cur.get("ts") or 0)
        if age < stale_after_s:
            raise LockHeld(
                f"another submit notebook holds {lock_uri(paths)} ({cur.get('holder')}, "
                f"{age / 60:.0f} min old). If that notebook is gone, wait "
                f"{(stale_after_s - age) / 60:.0f} min or delete the lock object."
            )
        print(f"breaking stale lock from {cur.get('holder')} ({age / 60:.0f} min old)")
        gcs.delete(lock_uri(paths))
    raise LockHeld(f"could not take {lock_uri(paths)}")


def refresh_lock(paths: RunPaths, holder: str) -> None:
    gcs.upload_text(lock_uri(paths), json.dumps({"holder": holder, "ts": time.time(), "at": registry.utcnow()}))


def release_lock(paths: RunPaths) -> None:
    gcs.delete(lock_uri(paths))


# --------------------------------------------------------------------------
# Sync with Batch
# --------------------------------------------------------------------------

_INSTANCE = re.compile(r"instances/(\S+)")


def parse_task(task: dict[str, Any]) -> dict[str, Any]:
    """index, state, start, end, instance from a ``gcloud batch tasks list`` row."""
    status = task.get("status") or {}
    start = end = instance = None
    for ev in status.get("statusEvents") or []:
        kind = str(ev.get("taskState") or ev.get("type") or "").upper()
        desc = str(ev.get("description") or "")
        m = _INSTANCE.search(desc)
        if m:
            instance = m.group(1)
        if kind == "RUNNING" and start is None:
            start = ev.get("eventTime")
        if kind in {"SUCCEEDED", "FAILED"}:
            end = ev.get("eventTime")
    return {
        "index": poll.task_index(task),
        "state": poll.task_state(task),
        "start": start,
        "end": end,
        "instance": instance,
    }


def _instance_info(job: dict[str, Any]) -> dict[str, Any]:
    groups = ((job.get("status") or {}).get("taskGroups") or {})
    for group in (groups.values() if isinstance(groups, dict) else groups):
        inst = (group.get("instances") or [{}])[0] if isinstance(group, dict) else {}
        if inst:
            return {
                "machine": inst.get("machineType") or "",
                "provisioning": inst.get("provisioningModel") or "",
                "task_pack": int(inst.get("taskPack") or 1),
                "boot_disk_gb": float((inst.get("bootDisk") or {}).get("sizeGb") or 53),
            }
    return {}


def _tasks_tsv_of(job: dict[str, Any]) -> str:
    group = monitor._first_task_group(job)
    env = ((group.get("taskSpec") or {}).get("environment") or {}).get("variables") or {}
    return str(env.get("TASKS_TSV") or "")


def _sample_ids_from_tsv(uri: str) -> list[str]:
    text = gcs.cat(uri)
    return [row["SAMPLE_ID"] for row in csv.DictReader(io.StringIO(text), delimiter="\t")]


def _age_s(ts: str | None) -> float:
    t = poll.parse_rfc3339(ts)
    if t is None:
        return 0.0
    if t.tzinfo is None:
        t = t.replace(tzinfo=timezone.utc)
    return (datetime.now(timezone.utc) - t).total_seconds()


def sync(
    paths: RunPaths,
    state: dict[str, Any],
    *,
    project: str,
    region: str,
    save_details: bool = True,
    log: Callable[[str], None] = print,
) -> dict[str, Any]:
    """Refresh job states from Batch. Lists tasks once per job when it finishes."""
    run_label = alloc.sanitize_label(paths.run_id)
    live = [j for j in submit.list_jobs(project=project, region=region) if monitor.labels_of(j).get("run-id") == run_label]
    live_by_id = {monitor.job_id_of(j): j for j in live}
    jobs = state.setdefault("jobs", {})

    for job_id, job in live_by_id.items():
        if job_id in jobs:
            continue
        tsv = _tasks_tsv_of(job)
        if not tsv or monitor.labels_of(job).get("stage") != "fused":
            continue
        labels = monitor.labels_of(job)
        log(f"recovering job {job_id} found in Batch but missing from jobs.json")
        jobs[job_id] = {
            "job_id": job_id,
            "wave": int(labels.get("wave") or 0),
            "kind": labels.get("kind") or "wave",
            "code_version": labels.get("code-version") or "",
            "submitted_at": (job.get("createTime") or "")[:19] + "Z",
            "tasks_uri": tsv,
            "sample_ids": _sample_ids_from_tsv(tsv),
            "state": "QUEUED",
        }

    for job_id, entry in jobs.items():
        if entry.get("state") in TERMINAL and entry.get("task_states") is not None:
            continue
        live_job = live_by_id.get(job_id)
        if live_job is None:
            if entry.get("state") == "SUBMITTING" and _age_s(entry.get("submitted_at")) < SUBMITTING_GRACE_S:
                continue
            log(f"job {job_id} is not in Batch any more; its samples go back in the queue")
            entry["state"] = "LOST"
            entry["task_states"] = []
            continue
        entry["state"] = monitor.job_state(live_job)
        entry["counts"] = monitor.task_counts(live_job)
        entry.update(_instance_info(live_job))
        if entry["state"] not in TERMINAL:
            continue
        tasks = [parse_task(t) for t in submit.list_tasks(job_id=job_id, project=project, region=region, quiet=True)]
        by_index = {t["index"]: t for t in tasks}
        entry["task_states"] = [str((by_index.get(i) or {}).get("state") or "") for i in range(len(entry["sample_ids"]))]
        spot = (entry.get("provisioning") or "").upper() == "SPOT"
        entry["cost"] = cost.job_cost_usd(
            tasks,
            machine=entry.get("machine") or "",
            spot=spot,
            boot_disk_gb=float(entry.get("boot_disk_gb") or 53),
        )
        entry["finished_at"] = registry.utcnow()
        if save_details:
            gcs.upload_text(
                task_detail_uri(paths, job_id),
                json.dumps({"job_id": job_id, "machine": entry.get("machine"), "spot": spot, "tasks": tasks}) + "\n",
            )
    return state


# --------------------------------------------------------------------------
# Pure planning logic
# --------------------------------------------------------------------------


@dataclass
class SampleState:
    status: str = "not_started"  # done, running, failed, gave_up, not_started
    attempts: int = 0
    last_job: str = ""
    last_index: int = -1


def _jobs_in_order(state: dict[str, Any]) -> list[dict[str, Any]]:
    return sorted(state.get("jobs", {}).values(), key=lambda j: (j.get("submitted_at") or "", j.get("job_id") or ""))


def sample_states(sample_ids: list[str], state: dict[str, Any], *, max_attempts: int) -> dict[str, SampleState]:
    out = {sid: SampleState() for sid in sample_ids}
    for job in _jobs_in_order(state):
        terminal = job.get("state") in TERMINAL and job.get("task_states") is not None
        task_states = job.get("task_states") or []
        for idx, sid in enumerate(job.get("sample_ids") or []):
            st = out.get(sid)
            if st is None or st.status == "done":
                continue
            if not terminal:
                st.status, st.last_job, st.last_index = "running", job["job_id"], idx
                continue
            ts = task_states[idx] if idx < len(task_states) else ""
            if ts == "SUCCEEDED":
                st.status, st.last_job, st.last_index = "done", job["job_id"], idx
            elif ts in NEVER_RAN:
                if st.status == "running":
                    st.status = "failed" if st.attempts else "not_started"
            else:
                st.attempts += 1
                st.status, st.last_job, st.last_index = "failed", job["job_id"], idx
    for st in out.values():
        if st.status == "failed" and st.attempts >= max_attempts:
            st.status = "gave_up"
    return out


def count_status(states: dict[str, SampleState]) -> dict[str, int]:
    counts = {k: 0 for k in ("done", "running", "failed", "gave_up", "not_started")}
    for st in states.values():
        counts[st.status] += 1
    return counts


def failure_window(state: dict[str, Any], *, code_version: str, window: int) -> tuple[int, int]:
    """(failed, finished) over the most recent first-attempt tasks on this code version."""
    failed = finished = 0
    for job in reversed(_jobs_in_order(state)):
        if job.get("code_version") != code_version or job.get("kind") not in FIRST_ATTEMPT_KINDS:
            continue
        if job.get("task_states") is not None:
            f = sum(1 for s in job["task_states"] if s not in NEVER_RAN and s != "SUCCEEDED")
            ok = sum(1 for s in job["task_states"] if s == "SUCCEEDED")
        else:
            counts = job.get("counts") or {}
            f, ok = int(counts.get("FAILED", 0)), int(counts.get("SUCCEEDED", 0))
        failed += f
        finished += f + ok
        if finished >= window:
            break
    return failed, finished


def spend(state: dict[str, Any]) -> dict[str, float]:
    vm = near = 0.0
    n_ran = 0.0
    for job in state.get("jobs", {}).values():
        c = job.get("cost") or {}
        vm += float(c.get("vm_usd") or 0)
        near += float(c.get("nearline_usd") or 0)
        n_ran += float(c.get("n_ran") or 0)
    return {"vm_usd": vm, "nearline_usd": near, "total_usd": vm + near, "n_ran": n_ran}


def unit_cost(state: dict[str, Any], settings: Settings, *, min_attempts: int = 500, margin: float = 1.10) -> tuple[float, str]:
    s = spend(state)
    if s["n_ran"] >= min_attempts:
        return s["total_usd"] / s["n_ran"] * margin, f"this run's mean over {s['n_ran']:.0f} attempts x {margin:.2f}"
    return settings.prior_unit_usd, "prior (Sep 29 on-demand + headroom) until 500 attempts finish"


@dataclass
class Plan:
    action: str  # submit, wait, halt, finished
    reason: str
    kind: str = ""
    sample_ids: list[str] = field(default_factory=list)
    numbers: dict[str, Any] = field(default_factory=dict)

    def explain(self) -> str:
        if self.action == "submit":
            n = self.numbers
            return (
                f"SUBMIT {len(self.sample_ids):,} samples ({self.kind}): {self.reason}. "
                f"Est. ${len(self.sample_ids) * n['unit_usd']:,.2f} at ${n['unit_usd']:.4f}/sample; "
                f"committed after this ${n['committed_usd'] + len(self.sample_ids) * n['unit_usd']:,.2f} "
                f"of ${n['budget_usd']:,.0f}."
            )
        return f"{self.action.upper()}: {self.reason}"


def plan(
    settings: Settings,
    sample_ids: list[str],
    state: dict[str, Any],
    *,
    code_version: str,
) -> Plan:
    states = sample_states(sample_ids, state, max_attempts=settings.max_attempts)
    counts = count_status(states)
    s = spend(state)
    unit, unit_note = unit_cost(state, settings)
    inflight = counts["running"]
    committed = s["total_usd"] + inflight * unit
    failed, finished = failure_window(state, code_version=code_version, window=settings.failure_window)
    numbers = {
        **counts,
        "spent_usd": s["total_usd"],
        "unit_usd": unit,
        "unit_note": unit_note,
        "committed_usd": committed,
        "budget_usd": settings.budget_usd,
        "code_version": code_version,
        "window_failed": failed,
        "window_finished": finished,
    }
    first = [sid for sid in sample_ids if states[sid].status == "not_started"]
    retry = [sid for sid in sample_ids if states[sid].status == "failed"]

    if not first and not retry:
        if inflight:
            return Plan("wait", f"last {inflight:,} samples in flight; nothing left to submit", numbers=numbers)
        return Plan(
            "finished",
            f"{counts['done']:,} done, {counts['gave_up']:,} gave up after {settings.max_attempts} attempts",
            numbers=numbers,
        )

    canaries = [j for j in _jobs_in_order(state) if j.get("code_version") == code_version and j.get("kind") == "canary"]
    canary_open = [j for j in canaries if j.get("task_states") is None]
    canary_ran = [
        j for j in canaries
        if j.get("task_states") is not None
        and sum(s not in NEVER_RAN for s in j["task_states"]) >= min(settings.min_finished_for_rate, len(j["sample_ids"]))
    ]
    if canary_open:
        return Plan("wait", f"canary {canary_open[0]['job_id']} still running; waves start after it passes", numbers=numbers)
    if not canary_ran:
        kind, pool = "canary", (first or retry)
        why = f"first submission for code version {code_version}" if not canaries else "the last canary never ran; new canary"
        room = settings.canary_n
    else:
        if finished >= settings.min_finished_for_rate and failed / finished > settings.max_failure_rate:
            return Plan(
                "halt",
                f"{failed:,} of the last {finished:,} first attempts on code {code_version} failed "
                f"({failed / finished:.0%} > {settings.max_failure_rate:.0%}). Check the failure logs, push a fix "
                f"(that starts a new canary), or raise max_failure_rate in expansion_hunter_03_setup if the failures are expected",
                numbers=numbers,
            )
        room = settings.wave_size - inflight
        threshold = math.ceil(settings.wave_size * settings.refill_fraction)
        if inflight >= threshold:
            return Plan(
                "wait",
                f"{inflight:,} samples in flight; next wave when fewer than {threshold:,}",
                numbers=numbers,
            )
        if first:
            kind, pool, why = "wave", first, f"{len(first):,} samples never attempted"
        else:
            kind, pool, why = "retry", retry, f"{len(retry):,} failed samples with attempts left"

    affordable = int(max(0.0, settings.budget_usd - committed) / unit + 1e-9) if unit > 0 else 0
    n = min(room, len(pool), affordable)
    if n <= 0:
        return Plan(
            "halt",
            f"budget: committed ${committed:,.2f} of ${settings.budget_usd:,.0f} leaves room for "
            f"{affordable} samples at ${unit:.4f}. Raise budget_usd in expansion_hunter_03_setup to continue",
            numbers=numbers,
        )
    if n < min(room, len(pool)):
        why += f"; capped at {n:,} by the budget"
    return Plan("submit", why, kind=kind, sample_ids=pool[:n], numbers=numbers)


# --------------------------------------------------------------------------
# Submit
# --------------------------------------------------------------------------


def make_job_id(pipeline_short: str, run_id: str, wave: int, n: int) -> str:
    stamp = datetime.now(timezone.utc).strftime("%H%M%S")
    run = alloc.sanitize_label(run_id, max_len=24)
    return alloc.sanitize_job_id(f"{pipeline_short}-{run}-w{wave:03d}-{n:03d}-{stamp}")


def submit_plan(
    paths: RunPaths,
    state: dict[str, Any],
    the_plan: Plan,
    *,
    settings: Settings,
    records_by_id: dict[str, dict[str, str]],
    code_version: str,
    project: str,
    region: str,
    build_rows: Callable[[list[dict[str, str]]], list[dict[str, str]]],
    build_job: Callable[..., dict[str, Any]],
    pipeline_short: str,
    dry_run: bool = False,
    log: Callable[[str], None] = print,
) -> list[str]:
    if the_plan.action != "submit" or not the_plan.sample_ids:
        return []
    jobs = state.setdefault("jobs", {})
    wave = 1 + max([int(j.get("wave") or 0) for j in jobs.values()] or [0])
    size = min(settings.job_size, settings.canary_n) if the_plan.kind == "canary" else settings.job_size
    ids = the_plan.sample_ids
    submitted: list[str] = []
    for n, start in enumerate(range(0, len(ids), size)):
        chunk = ids[start : start + size]
        job_id = make_job_id(pipeline_short, paths.run_id, wave, n)
        rows = build_rows([records_by_id[sid] for sid in chunk])
        tsv_uri = tasks_tsv_uri(paths, wave, job_id)
        labels = {
            "pipeline": alloc.sanitize_label(paths.pipeline.replace("_", "-")),
            "run-id": alloc.sanitize_label(paths.run_id),
            "stage": "fused",
            "kind": the_plan.kind,
            "wave": str(wave),
            "code-version": alloc.sanitize_label(code_version),
        }
        job = build_job(rows=rows, tasks_tsv_uri=tsv_uri, labels=labels, parallelism=len(chunk))
        nbytes = submit.assert_job_json_ok(job)
        if dry_run:
            log(f"dry run: {job_id} {len(chunk)} tasks, job JSON {nbytes:,} bytes")
            submitted.append(job_id)
            continue
        buf = io.StringIO()
        writer = csv.DictWriter(buf, fieldnames=list(rows[0].keys()), delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
        gcs.upload_text(tsv_uri, buf.getvalue())
        gcs.upload_text(f"{paths.prefix}/tasks/w{wave:03d}/{job_id}.job.json", json.dumps(job, indent=1) + "\n")
        jobs[job_id] = {
            "job_id": job_id,
            "wave": wave,
            "kind": the_plan.kind,
            "code_version": code_version,
            "submitted_at": registry.utcnow(),
            "tasks_uri": tsv_uri,
            "sample_ids": chunk,
            "state": "SUBMITTING",
        }
        save_state(paths, state)
        try:
            submit.submit_job(job_id=job_id, job=job, project=project, region=region, quiet=True)
        except Exception:
            jobs.pop(job_id, None)
            save_state(paths, state)
            raise
        jobs[job_id]["state"] = "QUEUED"
        submitted.append(job_id)
        log(f"submitted {job_id}: {len(chunk):,} samples")
    if not dry_run:
        save_state(paths, state)
        registry.append_event(
            paths,
            "wave_submitted",
            wave=wave,
            kind=the_plan.kind,
            code_version=code_version,
            n_samples=len(ids),
            job_ids=submitted,
            unit_usd=the_plan.numbers.get("unit_usd"),
            committed_usd=the_plan.numbers.get("committed_usd"),
        )
    return submitted


# --------------------------------------------------------------------------
# Reporting
# --------------------------------------------------------------------------


def failed_samples(
    sample_ids: list[str], state: dict[str, Any], *, max_attempts: int, limit: int = 50
) -> list[tuple[str, SampleState]]:
    states = sample_states(sample_ids, state, max_attempts=max_attempts)
    rows = [(sid, st) for sid, st in states.items() if st.status in {"failed", "gave_up"}]
    rows.sort(key=lambda r: r[1].last_job, reverse=True)
    return rows[:limit]


def projection(the_plan: Plan, settings: Settings) -> dict[str, float]:
    n = the_plan.numbers
    remaining = n["not_started"] + n["running"] + n["failed"]
    return {
        "remaining_samples": float(remaining),
        "projected_total_usd": n["spent_usd"] + remaining * n["unit_usd"],
    }


def summary_lines(the_plan: Plan, settings: Settings, *, total: int) -> list[str]:
    n = the_plan.numbers
    proj = projection(the_plan, settings)
    rate = (n["window_failed"] / n["window_finished"]) if n["window_finished"] else 0.0
    return [
        f"samples {total:,}: done {n['done']:,} ({n['done'] / max(total, 1):.1%})  running {n['running']:,}  "
        f"failed {n['failed']:,}  gave up {n['gave_up']:,}  not started {n['not_started']:,}",
        f"cost: est. spent ${n['spent_usd']:,.2f}  committed ${n['committed_usd']:,.2f}  "
        f"budget ${n['budget_usd']:,.0f}  projected total ${proj['projected_total_usd']:,.0f}",
        f"unit ${n['unit_usd']:.4f}/sample ({n['unit_note']})",
        f"code {n['code_version']}: recent first-attempt failures {n['window_failed']:,}/{n['window_finished']:,} ({rate:.1%})",
        the_plan.explain(),
    ]


# --------------------------------------------------------------------------
# One notebook run
# --------------------------------------------------------------------------


def tick(
    paths: RunPaths,
    *,
    project: str,
    code_version: str,
    pipeline_short: str,
    make_builders: Callable[[Settings], tuple[Callable[..., Any], Callable[..., Any]]],
    dry_run: bool = False,
    log: Callable[[str], None] = print,
) -> tuple[Plan, list[str], dict[str, Any], Settings]:
    """Sync, plan, and submit at most one wave. Safe to run again at any time."""
    cfg = Settings.from_run_json(registry.load_run_json(paths))
    if not dry_run:
        acquire_lock(paths)
    try:
        records = load_samples(paths, cfg)
        ids = [r["sample_id"] for r in records]
        state = load_state(paths)
        sync(paths, state, project=project, region=cfg.region, save_details=not dry_run, log=log)
        if not dry_run:
            save_state(paths, state)
        the_plan = plan(cfg, ids, state, code_version=code_version)
        for line in summary_lines(the_plan, cfg, total=len(ids)):
            log(line)
        build_rows, build_job = make_builders(cfg)
        submitted = submit_plan(
            paths,
            state,
            the_plan,
            settings=cfg,
            records_by_id={r["sample_id"]: r for r in records},
            code_version=code_version,
            project=project,
            region=cfg.region,
            build_rows=build_rows,
            build_job=build_job,
            pipeline_short=pipeline_short,
            dry_run=dry_run,
            log=log,
        )
        if the_plan.action == "halt" and not dry_run:
            registry.append_event(paths, "halted", reason=the_plan.reason, code_version=code_version)
    finally:
        if not dry_run:
            release_lock(paths)
    return the_plan, submitted, state, cfg
