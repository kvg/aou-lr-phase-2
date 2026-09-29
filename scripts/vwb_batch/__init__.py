"""Shared Verily Workbench / native Google Batch helpers.

Notebooks and pipelines use this for PET/VPC allocation, a GCS run ledger,
job submit/poll, and VM+Nearline cost estimates. ExpansionHunter is the first
pipeline; Locityper should plug in the same registry later.

Success is Batch task SUCCEEDED **and** the expected GCS objects existing.
Resubmit FAILED / missing-output samples with a new job id. Never retry a
sample that already has objects unless FORCE_RERUN is set.
"""

from __future__ import annotations

__all__ = [
    "alloc",
    "cohort",
    "cost",
    "gcs",
    "poll",
    "registry",
    "submit",
]
