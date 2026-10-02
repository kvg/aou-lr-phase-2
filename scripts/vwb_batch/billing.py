"""Actual cost from the Cloud Billing BigQuery export (optional; lags about a day).

Needs the *detailed usage cost* export enabled by a billing admin, and BigQuery
read access on that table for whoever runs the notebook. The export is not
retroactive. Batch VMs carry a ``batch-job-id`` label, so VM cost joins to
jobs exactly. Cloud Storage (Nearline reads) is billed to the bucket, so it is
summed by SKU over the run's time window and is shared with anything else that
read Nearline in this project at the same time.
"""

from __future__ import annotations

import json
import re
import subprocess
from typing import Any

_TABLE = re.compile(r"^[A-Za-z0-9_\-]+[.:][A-Za-z0-9_]+\.[A-Za-z0-9_]+$")


def _bq(sql: str) -> list[dict[str, Any]]:
    proc = subprocess.run(
        ["bq", "query", "--use_legacy_sql=false", "--format=json", "--max_rows=100000", sql],
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        raise RuntimeError((proc.stderr or proc.stdout).strip()[:600])
    out = proc.stdout.strip()
    return json.loads(out) if out else []


def _check(table: str, project: str) -> str:
    table = table.replace(":", ".", 1)
    if not _TABLE.match(table):
        raise ValueError(f"billing_table should look like project.dataset.table, got {table!r}")
    if not re.match(r"^[a-z0-9\-]+$", project):
        raise ValueError(f"bad project id {project!r}")
    return table


def vm_cost_by_job(*, table: str, project: str, start: str, end: str, job_prefix: str) -> dict[str, float]:
    """batch-job-id -> net cost (after credits) for jobs whose id starts with job_prefix."""
    table = _check(table, project)
    prefix = re.sub(r"[^a-z0-9\-]", "", job_prefix.lower())
    rows = _bq(
        f"""
        SELECT l.value AS job_id,
               SUM(cost) + SUM(IFNULL((SELECT SUM(c.amount) FROM UNNEST(credits) c), 0)) AS usd
        FROM `{table}`, UNNEST(labels) l
        WHERE l.key = 'batch-job-id' AND STARTS_WITH(l.value, '{prefix}')
          AND project.id = '{project}'
          AND usage_start_time BETWEEN TIMESTAMP('{start}') AND TIMESTAMP('{end}')
        GROUP BY job_id"""
    )
    return {r["job_id"]: float(r["usd"]) for r in rows}


def storage_cost_by_sku(*, table: str, project: str, start: str, end: str) -> list[dict[str, Any]]:
    table = _check(table, project)
    rows = _bq(
        f"""
        SELECT sku.description AS sku, SUM(cost) AS usd,
               SUM(usage.amount_in_pricing_units) AS qty, ANY_VALUE(usage.pricing_unit) AS unit
        FROM `{table}`
        WHERE service.description = 'Cloud Storage' AND project.id = '{project}'
          AND usage_start_time BETWEEN TIMESTAMP('{start}') AND TIMESTAMP('{end}')
        GROUP BY sku ORDER BY usd DESC"""
    )
    return [{**r, "usd": float(r["usd"])} for r in rows]


CONSOLE_HINT = (
    "Actual cost: billing export not configured (billing_table in run.json is empty).\n"
    "  Console: Billing > Reports, project {project}, Group by 'Label key: batch-job-id' for VM cost\n"
    "  per job, and Group by SKU with service Cloud Storage for Nearline. Days are Pacific time;\n"
    "  costs appear about a day later."
)
