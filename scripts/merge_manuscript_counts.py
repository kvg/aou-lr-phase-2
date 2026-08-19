#!/usr/bin/env python3
"""Merge phase1-only and phase2-only manuscript count TSVs into one comparison table."""

from __future__ import annotations

import argparse
from pathlib import Path


def load(path: str) -> dict[str, dict[str, str]]:
    rows = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            row = dict(zip(header, parts))
            rows[row["metric"]] = row
    return rows


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--phase1-tsv", required=True)
    p.add_argument("--phase2-tsv", required=True)
    p.add_argument("--out-tsv", required=True)
    args = p.parse_args()

    p1 = load(args.phase1_tsv)
    p2 = load(args.phase2_tsv)
    metrics = list(dict.fromkeys(list(p1) + list(p2)))

    Path(args.out_tsv).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out_tsv, "w") as out:
        out.write(
            "metric\tphase1_display\tphase2_display\t"
            "phase1_ge50\tphase1_ge20\tphase2_ge50\tphase2_ge20\n"
        )
        for m in metrics:
            a = p1.get(m, {})
            b = p2.get(m, {})
            out.write(
                f"{m}\t{a.get('phase1_display', '—')}\t{b.get('phase2_display', '—')}\t"
                f"{a.get('phase1_ge50', '')}\t{a.get('phase1_ge20', '')}\t"
                f"{b.get('phase2_ge50', '')}\t{b.get('phase2_ge20', '')}\n"
            )


if __name__ == "__main__":
    main()
