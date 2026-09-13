#!/usr/bin/env python3
"""Part 4.0 helpers: FLARE log retained markers + VCF variant-class breakdown.

FLARE (browning-lab) writes under ``Statistics``::

    markers           :  N

That N is markers retained after ref∩target identity matching and FLARE's
own min-maf / min-mac / array / excludemarkers filters — not the raw subset
VCF site count.

Example::

    python3 scripts/flare_site_stats.py parse-log shard.log
    python3 scripts/flare_site_stats.py vcf-classes shard.vcf.gz
    python3 scripts/flare_site_stats.py summarize \\
      --pop AFR --log shard.log --vcf shard.vcf.gz --out site_stats.tsv
"""

from __future__ import annotations

import argparse
import csv
import re
import subprocess
import sys
from pathlib import Path
from typing import Optional

# FLARE AdmixMain.statistics(): "  markers           :  N"
MARKERS_RE = re.compile(
    r"^\s*markers\s*:\s*([0-9][0-9,]*)\s*$", re.IGNORECASE | re.MULTILINE
)
# Fallbacks seen in related Beagle-family tools / older builds.
MARKERS_FALLBACK_RES = [
    re.compile(r"^\s*Number of markers\s*:\s*([0-9][0-9,]*)\s*$", re.I | re.M),
    re.compile(r"^\s*Retained markers\s*:\s*([0-9][0-9,]*)\s*$", re.I | re.M),
]


def _parse_int(raw: str) -> int:
    return int(raw.replace(",", "").strip())


def parse_flare_retained_markers(log_text: str) -> Optional[int]:
    """Return FLARE's retained marker count from a ``.log`` body, or None."""
    m = MARKERS_RE.search(log_text)
    if m:
        return _parse_int(m.group(1))
    for rx in MARKERS_FALLBACK_RES:
        m = rx.search(log_text)
        if m:
            return _parse_int(m.group(1))
    return None


def parse_flare_log_file(path: Path | str) -> Optional[int]:
    return parse_flare_retained_markers(Path(path).read_text(errors="replace"))


def _bcftools_count(vcf: Path, *extra: str) -> int:
    cmd = ["bcftools", "view", "-H", *extra, str(vcf)]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        raise RuntimeError(
            f"bcftools view failed ({proc.returncode}): {proc.stderr.strip()[:500]}"
        )
    # Count non-empty lines without loading huge stdout into a list twice.
    n = 0
    for line in proc.stdout.splitlines():
        if line.strip():
            n += 1
    return n


def vcf_index_n(vcf: Path) -> Optional[int]:
    proc = subprocess.run(
        ["bcftools", "index", "-n", str(vcf)],
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        return None
    text = proc.stdout.strip()
    if not text:
        return None
    try:
        return int(text.split()[0])
    except ValueError:
        return None


def vcf_variant_classes(vcf: Path) -> dict[str, Optional[float]]:
    """Breakdown for Part 4.0: total, biallelic SNVs, indels, multiallelic."""
    total = vcf_index_n(vcf)
    if total is None:
        # Slow path if no index.
        total = _bcftools_count(vcf)
    n_biallelic_snv = _bcftools_count(vcf, "-v", "snps", "-m2", "-M2")
    n_snps_any = _bcftools_count(vcf, "-v", "snps")
    n_indels = _bcftools_count(vcf, "-v", "indels")
    # Multiallelic ≈ anything with >2 alleles; approximate via total − biallelic records.
    # Prefer explicit: records that are not -m2 -M2.
    n_biallelic_any = _bcftools_count(vcf, "-m2", "-M2")
    n_multiallelic = max(0, int(total) - n_biallelic_any) if total is not None else None
    frac = (n_biallelic_snv / total) if total else None
    return {
        "n_sites_input": float(total) if total is not None else None,
        "n_biallelic_snv": float(n_biallelic_snv),
        "n_snps_any": float(n_snps_any),
        "n_indels": float(n_indels),
        "n_biallelic_any": float(n_biallelic_any),
        "n_multiallelic": float(n_multiallelic) if n_multiallelic is not None else None,
        "frac_biallelic_snv": frac,
    }


def summarize_row(
    *,
    population: str,
    log_path: Optional[Path] = None,
    vcf_path: Optional[Path] = None,
    n_sites_flare_retained: Optional[int] = None,
) -> dict[str, object]:
    retained = n_sites_flare_retained
    if retained is None and log_path is not None:
        retained = parse_flare_log_file(log_path)
    classes: dict[str, Optional[float]] = {}
    if vcf_path is not None:
        classes = vcf_variant_classes(Path(vcf_path))
    n_in = classes.get("n_sites_input")
    frac = classes.get("frac_biallelic_snv")
    return {
        "population": population,
        "n_sites_input": int(n_in) if n_in is not None else "",
        "n_sites_flare_retained": retained if retained is not None else "",
        "frac_biallelic_snv": f"{frac:.6g}" if frac is not None else "",
        "n_biallelic_snv": int(classes["n_biallelic_snv"]) if classes.get("n_biallelic_snv") is not None else "",
        "n_indels": int(classes["n_indels"]) if classes.get("n_indels") is not None else "",
        "n_multiallelic": int(classes["n_multiallelic"]) if classes.get("n_multiallelic") is not None else "",
    }


SITE_STATS_FIELDS = [
    "population",
    "n_sites_input",
    "n_sites_flare_retained",
    "frac_biallelic_snv",
    "n_biallelic_snv",
    "n_indels",
    "n_multiallelic",
]


def write_site_stats_tsv(rows: list[dict[str, object]], dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    with dest.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=SITE_STATS_FIELDS, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for row in rows:
            w.writerow(row)


def main(argv: Optional[list[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="cmd", required=True)

    pl = sub.add_parser("parse-log", help="Print retained marker count from a FLARE .log")
    pl.add_argument("log", type=Path)

    vc = sub.add_parser("vcf-classes", help="Print variant-class breakdown for a VCF")
    vc.add_argument("vcf", type=Path)

    sm = sub.add_parser("summarize", help="Write one-row site_stats.tsv")
    sm.add_argument("--pop", required=True)
    sm.add_argument("--log", type=Path)
    sm.add_argument("--vcf", type=Path)
    sm.add_argument("--out", type=Path, required=True)

    mg = sub.add_parser("merge", help="Concatenate per-pop site_stats TSVs")
    mg.add_argument("--tables", nargs="+", type=Path, required=True)
    mg.add_argument("--out", type=Path, required=True)

    args = p.parse_args(argv)
    if args.cmd == "parse-log":
        n = parse_flare_log_file(args.log)
        if n is None:
            raise SystemExit(f"no markers line in {args.log}")
        print(n)
        return 0
    if args.cmd == "vcf-classes":
        for k, v in vcf_variant_classes(args.vcf).items():
            print(f"{k}\t{v}")
        return 0
    if args.cmd == "summarize":
        row = summarize_row(population=args.pop, log_path=args.log, vcf_path=args.vcf)
        write_site_stats_tsv([row], args.out)
        print("wrote", args.out, row)
        return 0
    if args.cmd == "merge":
        rows: list[dict[str, object]] = []
        for path in args.tables:
            with path.open() as fh:
                r = csv.DictReader(fh, delimiter="\t")
                rows.extend(r)
        write_site_stats_tsv(rows, args.out)
        print("wrote", args.out, "n=", len(rows))
        return 0
    raise SystemExit(f"unknown cmd {args.cmd}")


if __name__ == "__main__":
    raise SystemExit(main())
