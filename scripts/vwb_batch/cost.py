"""VM + Nearline cost estimates. Billing export lags; PET often cannot see it.

Rates are us-central1 on-demand list prices used as planning numbers, not
invoices. Prefer p90/p95 of completed pilot samples over the mean.
"""

from __future__ import annotations

import csv
import io
import math
from typing import Any

# us-central1 on-demand, approximate. Batch maps cpuMilli/memoryMib onto a
# machine; e2-standard-4 matches 4 vCPU / 16 GiB (minicram). e2-standard-2
# matches 2 vCPU / 8 GiB; genotype requests 4 vCPU / 8 GiB so we still price
# e2-standard-4 as a conservative default.
VM_HOURLY_USD = {
    "e2-standard-2": 0.067014,
    "e2-standard-4": 0.134028,
    "n2-standard-4": 0.194236,
}
PD_STANDARD_GB_MONTH_USD = 0.04
NEARLINE_RETRIEVAL_GB_USD = 0.01
HOURS_PER_MONTH = 24 * 30


def quantile(values: list[float], q: float) -> float:
    if not values:
        raise ValueError("quantile of empty list")
    if not 0 <= q <= 1:
        raise ValueError("q must be in [0, 1]")
    xs = sorted(values)
    if len(xs) == 1:
        return xs[0]
    pos = q * (len(xs) - 1)
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        return xs[lo]
    frac = pos - lo
    return xs[lo] * (1 - frac) + xs[hi] * frac


def vm_cost_usd(
    *,
    seconds: float,
    machine: str = "e2-standard-4",
    boot_disk_gib: float = 50.0,
) -> float:
    hours = max(seconds, 0.0) / 3600.0
    vm = VM_HOURLY_USD[machine] * hours
    disk = PD_STANDARD_GB_MONTH_USD * boot_disk_gib * hours / HOURS_PER_MONTH
    return vm + disk


def nearline_cost_usd(network_bytes: float) -> float:
    return max(network_bytes, 0.0) / (1024**3) * NEARLINE_RETRIEVAL_GB_USD


def parse_transfer_stats(text: str) -> int:
    """Return total_network_bytes from a str-analysis transfer-stats TSV."""
    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    total = 0
    found = False
    for rec in reader:
        raw = rec.get("total_network_bytes")
        if raw is None or raw == "":
            continue
        total += int(float(raw))
        found = True
    if not found:
        raise ValueError("no total_network_bytes in transfer stats")
    return total


def sample_cost_usd(
    *,
    minicram_seconds: float,
    genotype_seconds: float,
    network_bytes: float,
    minicram_machine: str = "e2-standard-4",
    genotype_machine: str = "e2-standard-4",
    boot_disk_gib: float = 50.0,
) -> dict[str, float]:
    mini = vm_cost_usd(
        seconds=minicram_seconds, machine=minicram_machine, boot_disk_gib=boot_disk_gib
    )
    gt = vm_cost_usd(
        seconds=genotype_seconds, machine=genotype_machine, boot_disk_gib=boot_disk_gib
    )
    near = nearline_cost_usd(network_bytes)
    return {
        "minicram_vm_usd": mini,
        "genotype_vm_usd": gt,
        "nearline_usd": near,
        "total_usd": mini + gt + near,
        "network_bytes": float(network_bytes),
        "minicram_seconds": float(minicram_seconds),
        "genotype_seconds": float(genotype_seconds),
    }


def summarize_pilot(costs: list[dict[str, float]], *, q: float = 0.90) -> dict[str, Any]:
    totals = [c["total_usd"] for c in costs]
    if not totals:
        raise ValueError("no completed pilot costs")
    return {
        "n": len(totals),
        "mean_usd": sum(totals) / len(totals),
        "p50_usd": quantile(totals, 0.50),
        "p90_usd": quantile(totals, 0.90),
        "p95_usd": quantile(totals, 0.95),
        "max_usd": max(totals),
        "planning_usd": quantile(totals, q),
        "quantile": q,
        "mean_network_gib": sum(c["network_bytes"] for c in costs) / len(costs) / (1024**3),
    }


def n_keep(*, budget_usd: float, planning_usd: float, retry_margin: float, n_available: int) -> int:
    unit = planning_usd * retry_margin
    if unit <= 0:
        raise ValueError("planning cost must be > 0")
    return max(0, min(n_available, int(budget_usd // unit)))
