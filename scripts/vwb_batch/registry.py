"""GCS run ledger: run.json, append-only JSONL, status TSV, keep/pilot CSVs."""

from __future__ import annotations

import json
from dataclasses import dataclass
from datetime import datetime, timezone
from typing import Any
from uuid import uuid4

from . import gcs


def utcnow() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


@dataclass(frozen=True)
class RunPaths:
    prefix: str
    pipeline: str
    run_id: str

    @property
    def run_json(self) -> str:
        return f"{self.prefix}/run.json"

    @property
    def ledger(self) -> str:
        return f"{self.prefix}/ledger.jsonl"

    @property
    def status_tsv(self) -> str:
        return f"{self.prefix}/status.tsv"

    @property
    def cohort_csv(self) -> str:
        return f"{self.prefix}/cohort.csv"

    @property
    def pilot_csv(self) -> str:
        return f"{self.prefix}/pilot.csv"

    @property
    def keep_csv(self) -> str:
        return f"{self.prefix}/keep.csv"

    @property
    def cost_pilot_tsv(self) -> str:
        return f"{self.prefix}/cost_pilot.tsv"

    def shard_prefix(self, stage: str, shard: int) -> str:
        return f"{self.prefix}/shards/{stage}/shard-{shard:04d}"

    def tasks_tsv(self, stage: str, shard: int) -> str:
        return f"{self.shard_prefix(stage, shard)}/tasks.tsv"

    def job_json(self, stage: str, shard: int) -> str:
        return f"{self.shard_prefix(stage, shard)}/job.json"


def list_runs(*, output_bucket: str, pipeline: str) -> list[str]:
    prefix = f"{output_bucket.rstrip('/')}/batchRuns/{pipeline}/runs/"
    ids: list[str] = []
    for uri in gcs.ls(prefix):
        rid = uri.rstrip("/").rsplit("/", 1)[-1]
        if rid and rid != "runs":
            ids.append(rid)
    return sorted(set(ids))


def run_paths(*, output_bucket: str, pipeline: str, run_id: str) -> RunPaths:
    prefix = f"{output_bucket.rstrip('/')}/batchRuns/{pipeline}/runs/{run_id}"
    return RunPaths(prefix=prefix, pipeline=pipeline, run_id=run_id)


def append_event(paths: RunPaths, event: str, **fields: Any) -> dict[str, Any]:
    record = {
        "ts": utcnow(),
        "event": event,
        "pipeline": paths.pipeline,
        "run_id": paths.run_id,
        "uuid": uuid4().hex[:12],
        **fields,
    }
    existing = ""
    if gcs.exists(paths.ledger):
        existing = gcs.cat(paths.ledger)
        if existing and not existing.endswith("\n"):
            existing += "\n"
    gcs.upload_text(paths.ledger, existing + json.dumps(record, default=str) + "\n")
    return record


def load_events(paths: RunPaths) -> list[dict[str, Any]]:
    if not gcs.exists(paths.ledger):
        return []
    events: list[dict[str, Any]] = []
    for line in gcs.cat(paths.ledger).splitlines():
        line = line.strip()
        if not line:
            continue
        events.append(json.loads(line))
    return events


def write_run_json(paths: RunPaths, payload: dict[str, Any]) -> None:
    body = dict(payload)
    body.setdefault("pipeline", paths.pipeline)
    body.setdefault("run_id", paths.run_id)
    body.setdefault("updated_at", utcnow())
    gcs.upload_text(paths.run_json, json.dumps(body, indent=2, default=str) + "\n")


def load_run_json(paths: RunPaths) -> dict[str, Any]:
    if not gcs.exists(paths.run_json):
        return {}
    return json.loads(gcs.cat(paths.run_json))


STATUS_FIELDS = (
    "sample_id",
    "stage",
    "state",
    "objects_ok",
    "job_id",
    "shard",
    "task_index",
    "seconds",
    "network_bytes",
    "updated_at",
)


def format_status_tsv(rows: list[dict[str, Any]]) -> str:
    lines = ["\t".join(STATUS_FIELDS)]
    for row in rows:
        lines.append("\t".join(str(row.get(k, "")) for k in STATUS_FIELDS))
    return "\n".join(lines) + "\n"


def write_status(paths: RunPaths, rows: list[dict[str, Any]]) -> None:
    gcs.upload_text(paths.status_tsv, format_status_tsv(rows))


def load_status(paths: RunPaths) -> list[dict[str, str]]:
    if not gcs.exists(paths.status_tsv):
        return []
    lines = [ln for ln in gcs.cat(paths.status_tsv).splitlines() if ln.strip()]
    if not lines:
        return []
    header = lines[0].split("\t")
    return [dict(zip(header, ln.split("\t"))) for ln in lines[1:]]
