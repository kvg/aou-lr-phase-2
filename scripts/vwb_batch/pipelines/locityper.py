"""Locityper Batch adapter (ledger-ready; native Batch job builder later).

`locityper_02_run_batch` still uses dsub. When that path is ported, implement
`build_job` the same way as expansion_hunter.py and reuse registry/poll/cost.
"""

from __future__ import annotations

NAME = "locityper"
STAGES = ("minicram", "genotype")


def expected_objects(sample_id: str, stage: str, out_prefix: str) -> list[str]:
    prefix = f"{out_prefix.rstrip('/')}/{sample_id}"
    if stage == "minicram":
        return [
            f"{prefix}/{sample_id}.minicram.cram",
            f"{prefix}/{sample_id}.minicram.cram.crai",
        ]
    if stage == "genotype":
        return [f"{prefix}/{sample_id}.locityper.tar.gz"]
    raise ValueError(f"unknown stage {stage!r}")


def predecessor(stage: str) -> str | None:
    if stage == "genotype":
        return "minicram"
    return None


def build_job(**_kwargs):
    raise NotImplementedError(
        "Locityper native Batch is not wired yet; locityper_02 still uses dsub."
    )
