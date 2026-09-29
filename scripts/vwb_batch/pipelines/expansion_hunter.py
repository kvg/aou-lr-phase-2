"""ExpansionHunter native Batch pipeline (minicram → genotype)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from typing import Any

NAME = "expansion_hunter"
STAGES = ("minicram", "genotype")

PRINT_READS_DOCKER = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/aou-locityper-print-reads:0.1.2"
EH_DOCKER = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/aou-expansion-hunter:0.1.0"


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def _load_submit_batch():
    path = _repo_root() / "expansion_hunter" / "batch" / "submit_batch.py"
    name = "eh_submit_batch"
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise FileNotFoundError(path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def _gs(prefix: str, sample_id: str, name: str) -> str:
    return f"{prefix.rstrip('/')}/{sample_id}/{name}"


def expected_objects(sample_id: str, stage: str, out_prefix: str) -> list[str]:
    if stage == "minicram":
        return [
            _gs(out_prefix, sample_id, f"{sample_id}.minicram.cram"),
            _gs(out_prefix, sample_id, f"{sample_id}.minicram.cram.crai"),
            _gs(out_prefix, sample_id, f"{sample_id}.data_transfer_stats.tsv"),
        ]
    if stage == "genotype":
        return [
            _gs(out_prefix, sample_id, f"{sample_id}.EH.json"),
            _gs(out_prefix, sample_id, f"{sample_id}.EH.vcf"),
        ]
    raise ValueError(f"unknown stage {stage!r}")


def predecessor(stage: str) -> str | None:
    if stage == "genotype":
        return "minicram"
    return None


def transfer_stats_uri(sample_id: str, out_prefix: str) -> str:
    return _gs(out_prefix, sample_id, f"{sample_id}.data_transfer_stats.tsv")


def resource_stats_uri(sample_id: str, stage: str, out_prefix: str) -> str:
    if stage == "minicram":
        return _gs(out_prefix, sample_id, f"{sample_id}.minicram.resources.tsv")
    if stage == "genotype":
        return _gs(out_prefix, sample_id, f"{sample_id}.EH.resources.tsv")
    raise ValueError(f"unknown stage {stage!r}")


def worker_path(stage: str) -> Path:
    batch = _repo_root() / "expansion_hunter" / "batch"
    if stage == "minicram":
        return batch / "make_minicram.py"
    if stage == "genotype":
        return batch / "genotype.py"
    raise ValueError(f"unknown stage {stage!r}")


def image_for(stage: str, *, print_reads: str = PRINT_READS_DOCKER, eh: str = EH_DOCKER) -> str:
    if stage == "minicram":
        return print_reads
    if stage == "genotype":
        return eh
    raise ValueError(f"unknown stage {stage!r}")


def env_rows(
    records: list[dict[str, str]],
    *,
    stage: str,
    out_prefix: str,
    project: str,
) -> list[dict[str, str]]:
    eh = _load_submit_batch()
    if stage == "minicram":
        return eh.minicram_rows_from_records(records, out_prefix, project)
    if stage == "genotype":
        return eh.genotype_rows_from_records(records, out_prefix, project)
    raise ValueError(f"unknown stage {stage!r}")


def build_job(
    *,
    stage: str,
    rows: list[dict[str, str]],
    tasks_tsv_uri: str,
    out_prefix: str,
    project: str,
    region: str,
    sa: str,
    image: str,
    parallelism: int,
    labels: dict[str, str],
    max_run_duration: str = "28800s",
    compute: dict | None = None,
) -> dict[str, Any]:
    eh = _load_submit_batch()
    worker = worker_path(stage)
    common = dict(
        rows=rows,
        worker_path=worker,
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
    if stage == "minicram":
        return eh.build_minicram_job_from_rows(**common)
    if stage == "genotype":
        return eh.build_genotype_job_from_rows(**common)
    raise ValueError(f"unknown stage {stage!r}")
