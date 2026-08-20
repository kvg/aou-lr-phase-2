#!/usr/bin/env python3
"""Compare R-oracle vs Rust TractorMix.score TSV outputs (post round/signif rules)."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from typing import Dict, List


def r_round5(s: str) -> str:
    if s in ("", "NA"):
        return "NA"
    x = float(s)
    if math.isnan(x):
        return "NA"
    return str(round(x, 5))


def r_signif5(s: str) -> str:
    if s in ("", "NA"):
        return "NA"
    x = float(s)
    if math.isnan(x):
        return "NA"
    if x == 0.0:
        return "0"
    exp = math.floor(math.log10(abs(x)))
    factor = 10 ** (5 - 1 - int(exp))
    return str(round(x * factor) / factor)


def read_tsv(path: str) -> tuple[List[str], List[Dict[str, str]]]:
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        header = reader.fieldnames or []
        rows = [dict(row) for row in reader]
    return list(header), rows


def compare(r_path: str, rust_path: str) -> int:
    h1, r_rows = read_tsv(r_path)
    h2, rs_rows = read_tsv(rust_path)
    if h1 != h2:
        print("header mismatch", file=sys.stderr)
        print("R:", h1, file=sys.stderr)
        print("Rust:", h2, file=sys.stderr)
        return 1
    if len(r_rows) != len(rs_rows):
        print(f"row count mismatch: R={len(r_rows)} Rust={len(rs_rows)}", file=sys.stderr)
        return 1

    round5 = lambda c: c == "Chi2" or c.startswith("Eff_anc") or c.startswith("SE_anc")
    signif5 = lambda c: c == "P" or c.startswith("Pval_anc")

    mismatches = 0
    for i, (r, rs) in enumerate(zip(r_rows, rs_rows)):
        for col in h1:
            a = r.get(col, "")
            b = rs.get(col, "")
            if round5(col):
                aa, bb = r_round5(a), r_round5(b)
            elif signif5(col):
                aa, bb = r_signif5(a), r_signif5(b)
            else:
                aa, bb = a, b
            if aa != bb:
                mismatches += 1
                print(f"row {i} col {col}: R={a!r} Rust={b!r} (cmp {aa!r} vs {bb!r})", file=sys.stderr)
    if mismatches:
        print(f"FAIL: {mismatches} column mismatches", file=sys.stderr)
        return 1
    print(f"OK: {len(r_rows)} rows match between {r_path} and {rust_path}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("r_tsv", help="R oracle output TSV")
    ap.add_argument("rust_tsv", help="Rust tractor-mix-score output TSV")
    args = ap.parse_args()
    return compare(args.r_tsv, args.rust_tsv)


if __name__ == "__main__":
    raise SystemExit(main())
