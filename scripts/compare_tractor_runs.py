#!/usr/bin/env python3
"""Compare two Tractor-Mix result directories (e.g. legacy R score vs Rust score).

Uses the same round(5) / signif(5) rules as the R-oracle parity tests.
"""

from __future__ import annotations

import argparse
import csv
import filecmp
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional


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


def read_tsv(path: Path) -> tuple[List[str], List[Dict[str, str]]]:
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        header = reader.fieldnames or []
        rows = [dict(row) for row in reader]
    return list(header), rows


def phenotype_key(path: Path) -> str:
    for suffix in (".tractor_mix.tsv", ".tractor_mix.txt", ".tsv"):
        if path.name.endswith(suffix):
            return path.name[: -len(suffix)]
    return path.stem


def discover(dir_path: Path) -> dict[str, Path]:
    out: dict[str, Path] = {}
    for p in sorted(dir_path.glob("*.tsv")):
        out[phenotype_key(p)] = p
    return out


@dataclass
class CompareResult:
    phenotype: str
    old_path: Path
    new_path: Path
    n_rows: int
    n_scored: int
    byte_identical: bool
    mismatches: int
    chi2_mismatches: int
    p_only_mismatches: int
    ok: bool
    note: str = ""


def compare_pair(old_path: Path, new_path: Path) -> CompareResult:
    pheno = phenotype_key(old_path)
    try:
        h1, old_rows = read_tsv(old_path)
        h2, new_rows = read_tsv(new_path)
    except OSError as exc:
        return CompareResult(
            pheno, old_path, new_path, 0, 0, False, 0, 0, 0, False, str(exc)
        )

    if h1 != h2:
        return CompareResult(
            pheno, old_path, new_path, 0, 0, False, 0, 0, 0, False, "header mismatch"
        )
    if len(old_rows) != len(new_rows):
        return CompareResult(
            pheno,
            old_path,
            new_path,
            len(old_rows),
            0,
            False,
            -1,
            0,
            0,
            False,
            f"row count {len(old_rows)} vs {len(new_rows)}",
        )

    byte_identical = filecmp.cmp(old_path, new_path, shallow=False)
    round5 = lambda c: c == "Chi2" or c.startswith("Eff_anc") or c.startswith("SE_anc")
    signif5 = lambda c: c == "P" or c.startswith("Pval_anc")

    mismatches = chi2_mm = p_only = 0
    for r, n in zip(old_rows, new_rows):
        for col in h1:
            a, b = r.get(col, ""), n.get(col, "")
            if round5(col):
                aa, bb = r_round5(a), r_round5(b)
            elif signif5(col):
                aa, bb = r_signif5(a), r_signif5(b)
            else:
                aa, bb = a, b
            if aa != bb:
                mismatches += 1
                if col == "Chi2" or col.startswith("Eff_anc") or col.startswith("SE_anc"):
                    chi2_mm += 1
                elif col == "P" or col.startswith("Pval_anc"):
                    p_only += 1

    scored = sum(1 for r in old_rows if r.get("Chi2", "NA") not in ("", "NA"))
    ok = mismatches == 0
    return CompareResult(
        pheno,
        old_path,
        new_path,
        len(old_rows),
        scored,
        byte_identical,
        mismatches,
        chi2_mm,
        p_only,
        ok,
        "",
    )


def write_summary(path: Path, results: List[CompareResult]) -> None:
    cols = [
        "phenotype",
        "n_rows",
        "n_scored",
        "byte_identical",
        "mismatches",
        "chi2_mismatches",
        "p_only_mismatches",
        "ok",
        "note",
        "old_path",
        "new_path",
    ]
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in results:
            w.writerow(
                {
                    "phenotype": r.phenotype,
                    "n_rows": r.n_rows,
                    "n_scored": r.n_scored,
                    "byte_identical": int(r.byte_identical),
                    "mismatches": r.mismatches,
                    "chi2_mismatches": r.chi2_mismatches,
                    "p_only_mismatches": r.p_only_mismatches,
                    "ok": int(r.ok),
                    "note": r.note,
                    "old_path": str(r.old_path),
                    "new_path": str(r.new_path),
                }
            )


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("old_dir", type=Path, help="Legacy R-score result directory")
    ap.add_argument("new_dir", type=Path, help="Rust-score result directory")
    ap.add_argument(
        "--out",
        type=Path,
        default=Path("tractor_run_compare"),
        help="Write run_compare_summary.tsv here",
    )
    ap.add_argument(
        "--strict",
        action="store_true",
        help="Exit 1 if any phenotype fails parity",
    )
    ap.add_argument(
        "--verbose",
        action="store_true",
        help="Print per-phenotype status",
    )
    args = ap.parse_args()

    old_map = discover(args.old_dir)
    new_map = discover(args.new_dir)
    if not old_map:
        print(f"No *.tsv in {args.old_dir}", file=sys.stderr)
        return 1
    if not new_map:
        print(f"No *.tsv in {args.new_dir}", file=sys.stderr)
        return 1

    phenotypes = sorted(set(old_map) & set(new_map))
    only_old = sorted(set(old_map) - set(new_map))
    only_new = sorted(set(new_map) - set(old_map))
    if only_old:
        print(f"Only in old: {only_old}", file=sys.stderr)
    if only_new:
        print(f"Only in new: {only_new}", file=sys.stderr)
    if not phenotypes:
        print("No overlapping phenotypes", file=sys.stderr)
        return 1

    results = [compare_pair(old_map[p], new_map[p]) for p in phenotypes]
    args.out.mkdir(parents=True, exist_ok=True)
    summary_path = args.out / "run_compare_summary.tsv"
    write_summary(summary_path, results)

    n_ok = sum(1 for r in results if r.ok)
    n_byte = sum(1 for r in results if r.byte_identical)
    print(f"Compared {len(results)} phenotypes: {n_ok} OK, {n_byte} byte-identical")
    print(f"Wrote {summary_path}")

    for r in results:
        if args.verbose or not r.ok:
            status = "OK" if r.ok else "FAIL"
            extra = f" byte_identical={r.byte_identical}" if r.byte_identical else ""
            mm = "" if r.ok else f" mismatches={r.mismatches} (chi2={r.chi2_mismatches} p_only={r.p_only_mismatches})"
            note = f" {r.note}" if r.note else ""
            print(f"  [{status}] {r.phenotype}{extra}{mm}{note}")

    if args.strict and n_ok != len(results):
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
