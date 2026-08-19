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


def count_resolved(sites, size_flag: str) -> Counter:
    c: Counter = Counter()
    for s in sites:
        if s["source_vcf"] != "main":
            continue
        if s["svtype"] not in RESOLVED:
            continue
        if not s.get(size_flag):
            continue
        c[s["svtype"]] += 1
    return c


def count_phase(path: str) -> dict:
    c50: Counter = Counter()
    c20: Counter = Counter()
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
        if s.get("size_bin_ge50"):
            c50[s["svtype"]] += 1
        if s.get("size_bin_ge20"):
            c20[s["svtype"]] += 1
    return {
        "available": True,
        "resolved_ge50": {sv: c50[sv] for _, sv in METRICS},
        "resolved_ge20": {sv: c20[sv] for _, sv in METRICS},
        "breakends": n_bnd,
        "large_events_gt10kb": n_large,
    }


def phase_block(sites: list[dict] | None, available: bool) -> dict:
    if not available or sites is None:
        return {
            "available": False,
            "resolved_ge50": {k: None for _, k in METRICS},
            "resolved_ge20": {k: None for _, k in METRICS},
            "breakends": None,
            "large_events_gt10kb": None,
        }
    c50 = count_resolved(sites, "size_bin_ge50")
    c20 = count_resolved(sites, "size_bin_ge20")
    return {
        "available": True,
        "resolved_ge50": {sv: c50[sv] for _, sv in METRICS},
        "resolved_ge20": {sv: c20[sv] for _, sv in METRICS},
        "breakends": sum(1 for s in sites if s["source_vcf"] == "bnd"),
        "large_events_gt10kb": sum(1 for s in sites if s["source_vcf"] == "large"),
    }


def display_pair(ge50, ge20, available: bool) -> str:
    if not available or ge50 is None:
        return "—"
    return f"{ge50}; {ge20}"


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--phase1-sites")
    p.add_argument("--phase2-sites")
    p.add_argument("--out-tsv", required=True)
    p.add_argument("--out-json", required=True)
    args = p.parse_args()

    if not args.phase1_sites and not args.phase2_sites:
        raise SystemExit("Provide --phase1-sites and/or --phase2-sites")

    empty = {
        "available": False,
        "resolved_ge50": {k: None for _, k in METRICS},
        "resolved_ge20": {k: None for _, k in METRICS},
        "breakends": None,
        "large_events_gt10kb": None,
    }
    p1 = count_phase(args.phase1_sites) if args.phase1_sites else empty
    p2 = count_phase(args.phase2_sites) if args.phase2_sites else empty

    rows = []
    for metric, sv in METRICS:
        rows.append(
            {
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
