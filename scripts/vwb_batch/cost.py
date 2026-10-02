"""VM + Nearline cost estimates. Billing export lags; PET often cannot see it.

Rates are us-central1 on-demand list prices used as planning numbers, not
invoices. Prefer p90/p95 of completed pilot samples over the mean.

The calibrated model (``job_cost_usd``) prices each VM from its tasks' start
and end times plus a fixed startup allowance, which reproduced the Sep 29 2026
EH bill (modeled / billed VM cost 1.000 over nine labeled shards). The older
``sample_cost_usd`` helpers price every task as its own e2-standard-4 and
overestimate about 1.6x; they stay for the pilot notebooks.
"""

from __future__ import annotations

import csv
import io
import math
import re
from collections import defaultdict
from datetime import datetime
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


DEFAULT_COMPUTE = {
    "minicram": {"cpuMilli": 4000, "memoryMib": 16384, "bootDiskMib": 51200},
    "genotype": {"cpuMilli": 4000, "memoryMib": 8192, "bootDiskMib": 51200},
}


def parse_resource_stats(text: str) -> dict[str, float]:
    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    rec = next(reader, None)
    if rec is None:
        raise ValueError("empty resource stats")
    return {
        "peak_rss_bytes": float(rec["peak_rss_bytes"]),
        "peak_disk_bytes": float(rec["peak_disk_bytes"]),
        "wall_sec": float(rec["wall_sec"]),
        "cpu_sec": float(rec["cpu_sec"]),
        "cpu_cores_avg": float(rec["cpu_cores_avg"]),
        "peak_nproc": float(rec.get("peak_nproc") or 0),
        "n_samples": float(rec.get("n_samples") or 0),
        "n_cpus_os": float(rec.get("n_cpus_os") or 0),
    }


def machine_for(*, cpu_milli: int, memory_mib: int) -> str:
    if cpu_milli <= 2000 and memory_mib <= 8192:
        return "e2-standard-2"
    return "e2-standard-4"


def recommend_compute(
    rows: list[dict[str, float]],
    *,
    stage: str,
    q: float = 0.90,
    mem_headroom: float = 1.4,
    min_memory_mib: int = 2048,
    min_disk_mib: int = 20480,
    defaults: dict[str, int] | None = None,
) -> dict[str, Any]:
    """Pilot p90 RSS/CPU/disk → Batch computeResource. Never larger than the default."""
    if not rows:
        raise ValueError(f"no resource stats for {stage}")
    default = dict(defaults or DEFAULT_COMPUTE[stage])
    rss_mib = [r["peak_rss_bytes"] / (1024**2) for r in rows]
    disk_mib = [r["peak_disk_bytes"] / (1024**2) for r in rows]
    cores = [r["cpu_cores_avg"] for r in rows]
    p90_rss = quantile(rss_mib, q)
    p90_disk = quantile(disk_mib, q)
    p90_cores = quantile(cores, q)
    memory_mib = int(math.ceil(p90_rss * mem_headroom / 1024.0) * 1024)
    memory_mib = max(memory_mib, min_memory_mib)
    memory_mib = min(memory_mib, default["memoryMib"])
    if p90_cores < 0.75:
        cpu_milli = 1000
    elif p90_cores < 1.75:
        cpu_milli = 2000
    else:
        cpu_milli = 4000
    cpu_milli = min(cpu_milli, default["cpuMilli"])
    boot = int(math.ceil(p90_disk * 2.0 / 1024.0) * 1024)
    boot = max(boot, min_disk_mib)
    boot = min(boot, default["bootDiskMib"])
    return {
        "stage": stage,
        "n": len(rows),
        "quantile": q,
        "p90_rss_mib": p90_rss,
        "p90_disk_mib": p90_disk,
        "p90_cpu_cores": p90_cores,
        "cpuMilli": cpu_milli,
        "memoryMib": memory_mib,
        "bootDiskMib": boot,
        "machine": machine_for(cpu_milli=cpu_milli, memory_mib=memory_mib),
        "default": default,
    }


# --- Calibrated model (Sep 29 2026 billing, us-central1 on-demand list prices) ---
CORE_HOUR_USD = {"e2": 0.021813, "n2": 0.031611, "n1": 0.031611}
GIB_HOUR_USD = {"e2": 0.002926, "n2": 0.004237, "n1": 0.004237}
PD_BALANCED_GIB_MONTH_USD = 0.10
HOURS_PER_BILLING_MONTH = 730
# ASSUMPTION until a Spot run is billed: Spot VM price as a fraction of on-demand.
SPOT_FRACTION = 0.35
# Billed VM time no task timestamp covers (boot, image pull, teardown), fitted to billing.
VM_STARTUP_S = 85.0
# Nearline per minicram sample: retrieval 0.92 GiB ($0.0092), Class B reads ($0.0040),
# other Cloud Storage ($0.0002). Billed to the bucket, so it cannot be split by job.
NEARLINE_USD_PER_SAMPLE = 0.0134
# Tasks shorter than this failed before reading the CRAM and are not charged Nearline.
NEARLINE_MIN_RUN_S = 120.0
# Planning cost per attempt before a run has its own numbers: Sep 29 on-demand 2.3c plus headroom.
PRIOR_UNIT_USD = 0.025

_SHAPE = re.compile(r"^(e2|n1|n2)-(standard|highcpu|highmem)-(\d+)$")
_CUSTOM = re.compile(r"^(?:(e2|n1|n2)-)?custom-(\d+)-(\d+)$")
_GIB_PER_VCPU = {"standard": 4.0, "highcpu": 1.0, "highmem": 8.0}


def machine_shape(machine: str) -> tuple[int, float, str]:
    """(vCPU, GiB, family) for a Compute Engine machine type. Unknown types price as e2-highcpu-2."""
    m = _SHAPE.match(machine or "")
    if m:
        fam, kind, n = m.group(1), m.group(2), int(m.group(3))
        gib = n * _GIB_PER_VCPU[kind]
        if fam == "n1" and kind == "highcpu":
            gib = n * 0.9
        if fam == "n1" and kind == "highmem":
            gib = n * 6.5
        if fam == "n1" and kind == "standard":
            gib = n * 3.75
        return n, gib, fam
    m = _CUSTOM.match(machine or "")
    if m:
        return int(m.group(2)), int(m.group(3)) / 1024.0, m.group(1) or "n1"
    if machine == "e2-medium":
        return 2, 4.0, "e2"
    return 2, 2.0, "e2"


def vm_hourly_usd(machine: str, *, spot: bool = False, boot_disk_gb: float = 53.0) -> float:
    vcpu, gib, fam = machine_shape(machine)
    rate = vcpu * CORE_HOUR_USD[fam] + gib * GIB_HOUR_USD[fam]
    if spot:
        rate *= SPOT_FRACTION
    return rate + boot_disk_gb * PD_BALANCED_GIB_MONTH_USD / HOURS_PER_BILLING_MONTH


def _ts(value: str | None) -> datetime | None:
    if not value:
        return None
    text = str(value).strip().replace("Z", "+00:00")
    # Batch timestamps carry nanoseconds; fromisoformat takes at most microseconds.
    text = re.sub(r"(\.\d{6})\d+", r"\1", text)
    try:
        return datetime.fromisoformat(text)
    except ValueError:
        return None


def job_cost_usd(
    tasks: list[dict[str, Any]],
    *,
    machine: str,
    spot: bool,
    boot_disk_gb: float = 53.0,
    startup_s: float = VM_STARTUP_S,
    nearline_usd_per_sample: float = NEARLINE_USD_PER_SAMPLE,
    nearline_min_run_s: float = NEARLINE_MIN_RUN_S,
) -> dict[str, float]:
    """Estimate one finished job's cost from its tasks.

    Each task needs ``start``, ``end`` (RFC 3339) and ``instance`` (VM id).
    A VM is billed from its first task start to its last task end plus
    ``startup_s``. Nearline is charged per task that ran long enough to read
    the CRAM.
    """
    spans: dict[str, list[datetime]] = defaultdict(list)
    n_ran = 0
    n_nearline = 0
    for task in tasks:
        start, end = _ts(task.get("start")), _ts(task.get("end"))
        if start is None or end is None:
            continue
        n_ran += 1
        spans[str(task.get("instance") or f"solo-{id(task)}")].extend([start, end])
        if (end - start).total_seconds() >= nearline_min_run_s:
            n_nearline += 1
    rate = vm_hourly_usd(machine, spot=spot, boot_disk_gb=boot_disk_gb)
    vm_hours = sum(
        (max(times) - min(times)).total_seconds() / 3600.0 + startup_s / 3600.0
        for times in spans.values()
    )
    vm = vm_hours * rate
    near = n_nearline * nearline_usd_per_sample
    return {
        "vm_usd": vm,
        "nearline_usd": near,
        "total_usd": vm + near,
        "vm_hours": vm_hours,
        "n_vms": float(len(spans)),
        "n_ran": float(n_ran),
        "n_nearline": float(n_nearline),
    }
