#!/usr/bin/env python3
"""Prepare CADD-SV v2 coordinate-scoring BED (DEL/DUP/INS/INV, >=50 bp)."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from sv_site_utils import iter_sites, open_text  # noqa: E402


SUPPORTED = {"DEL", "DUP", "INS", "INV"}


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sites", required=True)
    p.add_argument(
        "--out-bed",
        required=True,
        help="CADD-SV BED: chrom, 0-based start, 1-based end, SVTYPE",
    )
    p.add_argument(
        "--max-dup-inv-bp",
        type=int,
        default=1_000_000,
        help="Omit DUP/INV longer than this (CADD-SV fetches the full span)",
    )
    args = p.parse_args()

    Path(args.out_bed).parent.mkdir(parents=True, exist_ok=True)
    n = 0
    n_skip = 0
    with open_text(args.out_bed, "wt") as out:
        for s in iter_sites(args.sites):
            if s["svtype"] not in SUPPORTED:
                continue
            if s.get("svlen") in (None, ""):
                continue
            if abs(int(s["svlen"])) < 50:
                continue
            start0 = int(s["pos"]) - 1
            end0 = int(s["end"]) if s.get("end") not in (None, "") else start0 + abs(int(s["svlen"]))
            if end0 <= start0:
                end0 = start0 + max(1, abs(int(s["svlen"])))
            span = end0 - start0
            if s["svtype"] in {"DUP", "INV"} and span > args.max_dup_inv_bp:
                n_skip += 1
                continue
            out.write(f"{s['chrom']}\t{start0}\t{end0}\t{s['svtype']}\n")
            n += 1
    print(
        f"Wrote {n} CADD-SV candidate intervals to {args.out_bed} "
        f"(skipped {n_skip} DUP/INV > {args.max_dup_inv_bp} bp)",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
