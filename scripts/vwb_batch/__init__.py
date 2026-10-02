"""Shared Verily Workbench / native Google Batch helpers.

Notebooks and pipelines use this for PET/VPC allocation, a GCS run ledger,
job submit/poll, a read-only job monitor, and VM+Nearline cost estimates. ExpansionHunter is the first
pipeline; Locityper should plug in the same registry later.

Cohort runs go through ``campaign``: one task per sample, submitted in
budget-capped waves after a canary, with a failure-rate stop. ``cost`` holds
the VM + Nearline model calibrated against billing, and ``billing`` reads the
BigQuery billing export when it is available.
"""

from __future__ import annotations

__all__ = [
    "alloc",
    "billing",
    "campaign",
    "cohort",
    "cost",
    "gcs",
    "monitor",
    "poll",
    "registry",
    "submit",
]
