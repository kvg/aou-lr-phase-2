#!/usr/bin/env python3
"""Summarize R-oracle vs Rust TSV comparison (same rules as compare_oracle_tsv.py)."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from collections import Counter
from typing import Dict, List, Tuple


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


def read_tsv(path: str) -> Tuple[List[str], List[Dict[str, str]]]:
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        header = reader.fieldnames or []
        rows = [dict(row) for row in reader]
    return list(header), rows


def column_kind(col: str) -> str:
    if col in ("CHR", "POS", "ID", "REF", "ALT"):
        return "variant_id"
    if col == "Chi2" or col.startswith("Eff_anc") or col.startswith("SE_anc"):
        return "round5"
    if col == "P" or col.startswith("Pval_anc"):
        return "signif5"
    if col.startswith("AC_count") or col.startswith("include_anc"):
        return "integer"
    return "other"


def compare_rows(
    header: List[str], r_rows: List[Dict[str, str]], rs_rows: List[Dict[str, str]]
) -> Tuple[int, Counter[str], List[str]]:
    mismatches = 0
    by_col: Counter[str] = Counter()
    examples: List[str] = []

    round5_cols = lambda c: column_kind(c) == "round5"
    signif5_cols = lambda c: column_kind(c) == "signif5"

    for i, (r, rs) in enumerate(zip(r_rows, rs_rows)):
        for col in header:
            a = r.get(col, "")
            b = rs.get(col, "")
            if round5_cols(col):
                aa, bb = r_round5(a), r_round5(b)
            elif signif5_cols(col):
                aa, bb = r_signif5(a), r_signif5(b)
            else:
                aa, bb = a, b
            if aa != bb:
                mismatches += 1
                by_col[col] += 1
                if len(examples) < 8:
                    examples.append(
                        f"  row {i} {col}: R={a!r} Rust={b!r} (cmp {aa!r} vs {bb!r})"
                    )
    return mismatches, by_col, examples


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("r_tsv", help="R oracle output TSV")
    ap.add_argument("rust_tsv", help="Rust tractor-mix-score output TSV")
    ap.add_argument(
        "--strict",
        action="store_true",
        help="Exit 1 if any column mismatch (same as compare_oracle_tsv.py)",
    )
    args = ap.parse_args()

    h1, r_rows = read_tsv(args.r_tsv)
    h2, rs_rows = read_tsv(args.rust_tsv)
    if h1 != h2:
        print("header mismatch", file=sys.stderr)
        return 1
    if len(r_rows) != len(rs_rows):
        print(
            f"row count mismatch: R={len(r_rows)} Rust={len(rs_rows)}",
            file=sys.stderr,
        )
        return 1

    scored = sum(1 for r in r_rows if r.get("Chi2", "NA") not in ("", "NA"))
    mismatches, by_col, examples = compare_rows(h1, r_rows, rs_rows)

    try:
        import filecmp

        byte_identical = filecmp.cmp(args.r_tsv, args.rust_tsv, shallow=False)
    except OSError:
        byte_identical = False

    chi2_cols = [c for c in by_col if c == "Chi2" or c.startswith("Eff_anc") or c.startswith("SE_anc")]
    p_cols = [c for c in by_col if c == "P" or c.startswith("Pval_anc")]

    print(f"Rows: {len(r_rows)} total, {scored} scored (Chi2 present)")
    print(f"Byte-identical TSV: {'yes' if byte_identical else 'no'}")
    print(f"Column mismatches (R formatting rules): {mismatches}")
    if mismatches == 0:
        print(f"OK: all rows match between {args.r_tsv} and {args.rust_tsv}")
    else:
        print("Mismatch breakdown (top columns):")
        for col, count in by_col.most_common(10):
            print(f"  {col}: {count}")
        if chi2_cols:
            print("  includes Chi2/Eff/SE (round-5) columns")
        if p_cols and not chi2_cols:
            print("  only P / Pval (signif-5) columns — likely floating-point CDF tie at display precision")
        print("Examples:")
        for line in examples:
            print(line)

    if args.strict and mismatches:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
