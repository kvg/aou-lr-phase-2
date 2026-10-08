#!/usr/bin/env python3
"""
Figure 2C inputs: allele rarefaction at repeat loci and PLVI-ranked example loci.

Reads allele tables for both phases (``trgt_allele_counts.py``, or
``aggregate_repeat_loci.py``; one row per locus x distinct dosage, no sample
ids) and the PLVI table of ``trgt_plvi.py merge``.

Locus set: catalog loci with a length-changing record in either phase and,
in each phase, at least ``--min-call-rate`` of that phase's haplotypes called.
A locus with no record in one phase is REF-length for every haplotype there.

Rarefaction is exact: the expected number of distinct alleles among n
haplotypes drawn without replacement from a locus with N called haplotypes
and allele counts N_i is ``sum_i 1 - C(N - N_i, n) / C(N, n)``.

Outputs (no sample ids)::

  --out-rarefaction  phase n_hap mean q25 q75 n_loci
  --out-examples     rank locus_id chrom start end motif period plvi plvi_pct_in_group
                     phase dosage_bp dosage_units n_hap_allele in_phase1
  --out-summary      locus-set size, alleles per locus, share of Phase 2 alleles
                     absent from Phase 1, selection settings
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
from pathlib import Path
from typing import Iterable, Optional, TextIO

import numpy as np


def open_text(path: str, mode: str = "rt") -> TextIO:
    return gzip.open(path, mode) if path.endswith((".gz", ".bgz")) else open(path, mode)


class Locus:
    __slots__ = ("chrom", "start", "end", "motif", "period", "n_hap", "n_total", "bp", "cnt")

    def __init__(self, r: dict):
        self.chrom, self.start, self.end = r["chrom"], int(r["start"]), int(r["end"])
        self.motif, self.period = r["motif"], int(r["period"])
        self.n_hap = int(r["n_hap"])
        self.n_total = self.n_hap + int(r["n_hap_missing"])
        self.bp: list[int] = []
        self.cnt: list[int] = []


def load_alleles(paths: Iterable[str]) -> tuple[dict[str, Locus], int]:
    loci: dict[str, Locus] = {}
    total = 0
    for path in paths:
        with open_text(path) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                L = loci.get(r["locus_id"])
                if L is None:
                    L = loci[r["locus_id"]] = Locus(r)
                    total = max(total, L.n_total)
                L.bp.append(int(r["dosage_bp"]))
                L.cnt.append(int(r["n_hap_allele"]))
    return loci, total


def log_factorials(n: int) -> np.ndarray:
    return np.concatenate([[0.0], np.cumsum(np.log(np.arange(1, n + 1, dtype=np.float64)))])


def expected_alleles(counts: list[np.ndarray], grid: np.ndarray, lf: np.ndarray) -> np.ndarray:
    """(n_loci, len(grid)) expected distinct alleles; each locus has N = sum of its counts >= max(grid)."""
    sizes = np.array([c.size for c in counts])
    ni = np.concatenate(counts).astype(np.int64)
    big_n = np.repeat(np.array([c.sum() for c in counts], np.int64), sizes)
    owner = np.repeat(np.arange(len(counts)), sizes)
    out = np.empty((len(counts), grid.size))
    for j, n in enumerate(grid):
        rest = big_n - ni
        ok = rest >= n
        log_p = np.full(ni.size, -np.inf)
        log_p[ok] = lf[rest[ok]] - lf[rest[ok] - n] - lf[big_n[ok]] + lf[big_n[ok] - n]
        out[:, j] = np.bincount(owner, weights=1.0 - np.exp(log_p), minlength=len(counts))
    return out


def locus_counts(L: Optional[Locus], n_total: int) -> np.ndarray:
    if L is None:
        return np.array([n_total], np.int64)
    return np.array(L.cnt, np.int64)


def select_examples(plvi_path: str, keep: set[str], p1: dict[str, Locus], p2: dict[str, Locus], *,
                    n: int, min_alleles: int, max_motif_len: int, distinct_groups: bool) -> list[dict]:
    chosen: list[dict] = []
    groups: set[str] = set()
    with open_text(plvi_path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            lid = r["trid"]
            if lid not in keep or int(r["motif_len"]) > max_motif_len:
                continue
            L1, L2 = p1.get(lid), p2.get(lid)
            if L1 is None or L2 is None or len(L1.bp) < 2 or len(L2.bp) < min_alleles:
                continue
            if distinct_groups and r["motif_group"] in groups:
                continue
            groups.add(r["motif_group"])
            chosen.append(r)
            if len(chosen) == n:
                break
    return chosen


def main(argv: Optional[list[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--p1-alleles", nargs="+", required=True)
    p.add_argument("--p2-alleles", nargs="+", required=True)
    p.add_argument("--plvi", required=True, help="trgt_plvi.py merge output (sorted by within-group percentile)")
    p.add_argument("--min-call-rate", type=float, default=0.9)
    p.add_argument("--grid-points", type=int, default=60)
    p.add_argument("--grid-min", type=int, default=50)
    p.add_argument("--n-examples", type=int, default=3)
    p.add_argument("--min-alleles", type=int, default=5, help="Distinct Phase 2 alleles an example locus needs")
    p.add_argument("--max-motif-len", type=int, default=6)
    p.add_argument("--any-groups", action="store_true", help="Allow two examples from one motif-length group")
    p.add_argument("--out-rarefaction", required=True)
    p.add_argument("--out-examples", required=True)
    p.add_argument("--out-summary", required=True)
    args = p.parse_args(argv)

    p1, h1 = load_alleles(args.p1_alleles)
    p2, h2 = load_alleles(args.p2_alleles)
    min1, min2 = math.ceil(args.min_call_rate * h1), math.ceil(args.min_call_rate * h2)
    keep = sorted(
        lid for lid in set(p1) | set(p2)
        if (lid not in p1 or p1[lid].n_hap >= min1) and (lid not in p2 or p2[lid].n_hap >= min2)
    )
    keep_set = set(keep)

    rows = []
    summary: dict = {
        "haplotypes": {"phase1": h1, "phase2": h2},
        "min_call_rate": args.min_call_rate,
        "loci": {"phase1_tables": len(p1), "phase2_tables": len(p2), "rarefaction_set": len(keep)},
    }
    for phase, loci, total, cap in (("phase1", p1, h1, min1), ("phase2", p2, h2, min2)):
        counts = [locus_counts(loci.get(lid), total) for lid in keep]
        grid = np.unique(np.geomspace(args.grid_min, cap, args.grid_points).astype(np.int64))
        if phase == "phase2" and h1 <= cap:
            grid = np.unique(np.r_[grid, h1])
        exp = expected_alleles(counts, grid, log_factorials(max(total, 1)))
        for j, n in enumerate(grid):
            col = exp[:, j]
            rows.append((phase, int(n), col.mean(), np.percentile(col, 25), np.percentile(col, 75), len(keep)))
        full = np.array([c.size for c in counts])
        summary[f"{phase}_alleles_per_locus"] = {"mean": float(full.mean()), "median": float(np.median(full))}
        if phase == "phase2" and h1 <= cap:
            summary["phase2_at_phase1_size_mean"] = float(exp[:, list(grid).index(h1)].mean())
    with open_text(args.out_rarefaction, "wt") as fh:
        fh.write("phase\tn_hap\tmean\tq25\tq75\tn_loci\n")
        for phase, n, m, a, b, k in rows:
            fh.write(f"{phase}\t{n}\t{m:.4f}\t{a:.4f}\t{b:.4f}\t{k}\n")

    new = total2 = 0
    for lid in keep:
        a2 = set(p2[lid].bp) if lid in p2 else {0}
        a1 = set(p1[lid].bp) if lid in p1 else {0}
        new += len(a2 - a1)
        total2 += len(a2)
    summary["phase2_alleles_absent_from_phase1"] = {"alleles": new, "of": total2, "frac": new / max(total2, 1)}

    chosen = select_examples(args.plvi, keep_set, p1, p2, n=args.n_examples, min_alleles=args.min_alleles,
                             max_motif_len=args.max_motif_len, distinct_groups=not args.any_groups)
    with open_text(args.out_examples, "wt") as fh:
        fh.write("rank\tlocus_id\tchrom\tstart\tend\tmotif\tperiod\tplvi\tplvi_pct_in_group\tphase\t"
                 "dosage_bp\tdosage_units\tn_hap_allele\tin_phase1\n")
        for rank, r in enumerate(chosen, 1):
            lid = r["trid"]
            seen = set(p1[lid].bp)
            for phase, L in (("phase1", p1[lid]), ("phase2", p2[lid])):
                for bp, c in sorted(zip(L.bp, L.cnt)):
                    fh.write(f"{rank}\t{lid}\t{L.chrom}\t{L.start}\t{L.end}\t{L.motif}\t{L.period}\t{r['plvi']}\t"
                             f"{r['plvi_pct_in_group']}\t{phase}\t{bp}\t{bp / L.period:.4g}\t{c}\t{int(bp in seen)}\n")
    summary["examples"] = [
        {"rank": i, "locus_id": r["trid"], "motif": r["motif"], "plvi": float(r["plvi"]),
         "plvi_pct_in_group": float(r["plvi_pct_in_group"]),
         "alleles": {"phase1": len(p1[r["trid"]].bp), "phase2": len(p2[r["trid"]].bp)}}
        for i, r in enumerate(chosen, 1)
    ]
    summary["selection"] = {"min_alleles": args.min_alleles, "max_motif_len": args.max_motif_len,
                            "distinct_groups": not args.any_groups}
    Path(args.out_summary).write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({k: summary[k] for k in ("loci", "phase2_alleles_absent_from_phase1")}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
