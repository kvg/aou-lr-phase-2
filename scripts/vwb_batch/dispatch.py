"""Write a shard's TASKS_TSV, build job JSON, optionally submit."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Callable

from . import alloc, cohort, gcs, registry, submit
from .registry import RunPaths


def prepare_shard(
    *,
    paths: RunPaths,
    stage: str,
    shard: int,
    records: list[dict[str, str]],
    env_rows: list[dict[str, str]],
    build_job: Callable[..., dict[str, Any]],
    image: str,
    out_prefix: str,
    project: str,
    region: str,
    sa: str,
    parallelism: int,
    max_run_duration: str,
    scratch: Path,
) -> dict[str, Any]:
    if not records:
        raise ValueError("empty shard")
    local_dir = scratch / paths.run_id / stage / f"shard-{shard:04d}"
    local_dir.mkdir(parents=True, exist_ok=True)
    local_tsv = local_dir / "tasks.tsv"
    local_job = local_dir / "job.json"
    cohort.write_tasks_tsv(local_tsv, env_rows)
    tasks_uri = paths.tasks_tsv(stage, shard)
    gcs.upload_file(local_tsv, tasks_uri)
    job_id = alloc.job_id(pipeline=paths.pipeline, stage=stage, run_id=paths.run_id, shard=shard)
    labels = alloc.job_labels(pipeline=paths.pipeline, run_id=paths.run_id, stage=stage, shard=shard)
    job = build_job(
        stage=stage,
        rows=env_rows,
        tasks_tsv_uri=tasks_uri,
        out_prefix=out_prefix,
        project=project,
        region=region,
        sa=sa,
        image=image,
        parallelism=min(int(parallelism), len(env_rows)),
        labels=labels,
        max_run_duration=max_run_duration,
    )
    size = submit.assert_job_json_ok(job)
    local_job.write_text(json.dumps(job, indent=2) + "\n", encoding="utf-8")
    gcs.upload_file(local_job, paths.job_json(stage, shard))
    sample_ids = [r["sample_id"] for r in records]
    return {
        "job_id": job_id,
        "job": job,
        "job_json_bytes": size,
        "tasks_uri": tasks_uri,
        "n_tasks": len(env_rows),
        "sample_ids": sample_ids,
        "local_job": str(local_job),
    }


def maybe_submit(
    *,
    prepared: dict[str, Any],
    paths: RunPaths,
    stage: str,
    shard: int,
    project: str,
    region: str,
    do_submit: bool,
) -> str | None:
    if not do_submit:
        print(
            f"{stage} shard {shard}: dry-run {prepared['n_tasks']} tasks "
            f"job_id={prepared['job_id']} json={prepared['job_json_bytes']} bytes"
        )
        return None
    submit.submit_job(
        job_id=prepared["job_id"],
        job=prepared["job"],
        project=project,
        region=region,
    )
    registry.append_event(
        paths,
        "job_submitted",
        stage=stage,
        shard=shard,
        job_id=prepared["job_id"],
        n_tasks=prepared["n_tasks"],
        sample_ids=prepared["sample_ids"],
        tasks_uri=prepared["tasks_uri"],
    )
    return prepared["job_id"]
