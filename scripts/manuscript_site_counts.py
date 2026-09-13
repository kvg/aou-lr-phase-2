#!/usr/bin/env python3
"""Emit manuscript SV site-count table rows for one or both phases."""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from sv_site_utils import iter_sites, open_text  # noqa: E402

RESOLVED = ("DEL", "DUP", "INS", "INV")
METRICS = (
    ("deletions", "DEL"),
    ("duplications", "DUP"),
    ("insertions", "INS"),
    ("inversions", "INV"),
)
# GIAB-style sequence-context strata for DEL/INS (parent class row is ALL).
# US is the complement of RM∪SD∪SR; RM/SD/SR may overlap; CMRG is independent.
CONTEXT_CLASSES = ("DEL", "INS")
STRATA = ("us", "rm", "sd", "sr", "cmrg")


def site_strata(site: dict) -> tuple[str, ...]:
    rm = bool(site.get("hit_rmsk"))
    sd = bool(site.get("hit_genomicSuperDups"))
    sr = bool(site.get("hit_simpleRepeat"))
    labels = []
    if not (rm or sd or sr):
        labels.append("us")
    if rm:
        labels.append("rm")
    if sd:
        labels.append("sd")
    if sr:
        labels.append("sr")
    if site.get("hit_cmrg"):
        labels.append("cmrg")
    return tuple(labels)


def _empty_resolved() -> dict:
    return {sv: None for _, sv in METRICS}


def _empty_context() -> dict:
    return {sv: {k: None for k in STRATA} for sv in CONTEXT_CLASSES}


def empty_phase() -> dict:
    return {
        "available": False,
        "resolved_ge50": _empty_resolved(),
        "resolved_ge20": _empty_resolved(),
        "context_ge50": _empty_context(),
        "context_ge20": _empty_context(),
        "breakends": None,
        "large_events_gt10kb": None,
    }


def _pack(
    c50: Counter,
    c20: Counter,
    c50_ctx: Counter,
    c20_ctx: Counter,
    n_bnd: int,
    n_large: int,
) -> dict:
    def ctx(counter: Counter, sv: str) -> dict:
        return {k: counter[(sv, k)] for k in STRATA}

    return {
        "available": True,
        "resolved_ge50": {sv: c50[sv] for _, sv in METRICS},
        "resolved_ge20": {sv: c20[sv] for _, sv in METRICS},
        "context_ge50": {sv: ctx(c50_ctx, sv) for sv in CONTEXT_CLASSES},
        "context_ge20": {sv: ctx(c20_ctx, sv) for sv in CONTEXT_CLASSES},
        "breakends": n_bnd,
        "large_events_gt10kb": n_large,
    }


def count_phase(path: str) -> dict:
    c50: Counter = Counter()
    c20: Counter = Counter()
    c50_ctx: Counter = Counter()
    c20_ctx: Counter = Counter()
    n_bnd = 0
    n_large = 0
    for s in iter_sites(path):
        src = s["source_vcf"]
        if src == "bnd":
            n_bnd += 1
            continue
        if src == "large":
            n_large += 1
            continue
        if src != "main" or s["svtype"] not in RESOLVED:
            continue
        labels = site_strata(s) if s["svtype"] in CONTEXT_CLASSES else ()
        if s.get("size_bin_ge50"):
            c50[s["svtype"]] += 1
            for lab in labels:
                c50_ctx[(s["svtype"], lab)] += 1
        if s.get("size_bin_ge20"):
            c20[s["svtype"]] += 1
            for lab in labels:
                c20_ctx[(s["svtype"], lab)] += 1
    return _pack(c50, c20, c50_ctx, c20_ctx, n_bnd, n_large)


def display_pair(ge50, ge20, available: bool) -> str:
    if not available or ge50 is None:
        return "—"
    return f"{ge50}; {ge20}"


def pct_display(phase: dict, sv: str) -> str:
    """Percent of ≥50 bp sites in each GIAB stratum (US/RM/SD/SR/CMRG)."""
    if not phase["available"]:
        return "—"
    total = phase["resolved_ge50"][sv]
    if not total:
        return "—"
    ctx = phase["context_ge50"][sv]
    return " / ".join(f"{100.0 * ctx[k] / total:.1f}" for k in STRATA)


def _resolved_row(metric: str, sv: str, p1: dict, p2: dict) -> dict:
    return {
        "metric": metric,
        "phase1_display": display_pair(
            p1["resolved_ge50"][sv], p1["resolved_ge20"][sv], p1["available"]
        ),
        "phase2_display": display_pair(
            p2["resolved_ge50"][sv], p2["resolved_ge20"][sv], p2["available"]
        ),
        "phase1_ge50": p1["resolved_ge50"][sv],
        "phase1_ge20": p1["resolved_ge20"][sv],
        "phase2_ge50": p2["resolved_ge50"][sv],
        "phase2_ge20": p2["resolved_ge20"][sv],
    }


def _context_row(metric: str, sv: str, stratum: str, p1: dict, p2: dict) -> dict:
    return {
        "metric": metric,
        "phase1_display": display_pair(
            p1["context_ge50"][sv][stratum],
            p1["context_ge20"][sv][stratum],
            p1["available"],
        ),
        "phase2_display": display_pair(
            p2["context_ge50"][sv][stratum],
            p2["context_ge20"][sv][stratum],
            p2["available"],
        ),
        "phase1_ge50": p1["context_ge50"][sv][stratum],
        "phase1_ge20": p1["context_ge20"][sv][stratum],
        "phase2_ge50": p2["context_ge50"][sv][stratum],
        "phase2_ge20": p2["context_ge20"][sv][stratum],
    }


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--phase1-sites")
    p.add_argument("--phase2-sites")
    p.add_argument("--out-tsv", required=True)
    p.add_argument("--out-json", required=True)
    args = p.parse_args()

    if not args.phase1_sites and not args.phase2_sites:
        raise SystemExit("Provide --phase1-sites and/or --phase2-sites")

    empty = empty_phase()
    p1 = count_phase(args.phase1_sites) if args.phase1_sites else empty
    p2 = count_phase(args.phase2_sites) if args.phase2_sites else empty

    rows = []
    for metric, sv in METRICS:
        rows.append(_resolved_row(metric, sv, p1, p2))
        if sv in CONTEXT_CLASSES:
            for stratum in STRATA:
                rows.append(_context_row(f"{metric}_{stratum}", sv, stratum, p1, p2))
            rows.append(
                {
                    "metric": f"{metric}_context_pct",
                    "phase1_display": pct_display(p1, sv),
                    "phase2_display": pct_display(p2, sv),
                    "phase1_ge50": None,
                    "phase1_ge20": None,
                    "phase2_ge50": None,
                    "phase2_ge20": None,
                }
            )

    rows.append(
        {
            "metric": "breakends",
            "phase1_display": "—" if not p1["available"] else str(p1["breakends"]),
            "phase2_display": "—" if not p2["available"] else str(p2["breakends"]),
            "phase1_ge50": p1["breakends"],
            "phase1_ge20": p1["breakends"],
            "phase2_ge50": p2["breakends"],
            "phase2_ge20": p2["breakends"],
        }
    )
    rows.append(
        {
            "metric": "large_events_gt10kb",
            "phase1_display": "—" if not p1["available"] else str(p1["large_events_gt10kb"]),
            "phase2_display": "—" if not p2["available"] else str(p2["large_events_gt10kb"]),
            "phase1_ge50": p1["large_events_gt10kb"],
            "phase1_ge20": p1["large_events_gt10kb"],
            "phase2_ge50": p2["large_events_gt10kb"],
            "phase2_ge20": p2["large_events_gt10kb"],
        }
    )

    # Phase 1 never has bnd/large partitions in this project
    if p1["available"]:
        for r in rows:
            if r["metric"] in ("breakends", "large_events_gt10kb"):
                r["phase1_display"] = "—"
                r["phase1_ge50"] = None
                r["phase1_ge20"] = None

    Path(args.out_tsv).parent.mkdir(parents=True, exist_ok=True)
    with open_text(args.out_tsv, "wt") as out:
        out.write(
            "metric\tphase1_display\tphase2_display\t"
            "phase1_ge50\tphase1_ge20\tphase2_ge50\tphase2_ge20\n"
        )
        for r in rows:
            out.write(
                f"{r['metric']}\t{r['phase1_display']}\t{r['phase2_display']}\t"
                f"{'' if r['phase1_ge50'] is None else r['phase1_ge50']}\t"
                f"{'' if r['phase1_ge20'] is None else r['phase1_ge20']}\t"
                f"{'' if r['phase2_ge50'] is None else r['phase2_ge50']}\t"
                f"{'' if r['phase2_ge20'] is None else r['phase2_ge20']}\n"
            )

    Path(args.out_json).write_text(
        json.dumps({"rows": rows, "phase1": p1, "phase2": p2}, indent=2) + "\n"
    )
    print(f"Wrote {args.out_tsv}")


if __name__ == "__main__":
    main()
