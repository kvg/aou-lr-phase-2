#!/usr/bin/env python3
"""Rewrite a chr:start-end region to match VCF contig names (chr22 vs 22)."""

from __future__ import annotations

import argparse
import re
import sys
from typing import Iterable, Optional


def contig_aliases(chrom: str) -> list[str]:
    chrom = chrom.strip()
    out = [chrom]
    if chrom.startswith("chr"):
        out.append(chrom[3:])
    elif chrom:
        out.append("chr" + chrom)
    # unique, preserve order
    seen: set[str] = set()
    uniq: list[str] = []
    for name in out:
        if name not in seen:
            seen.add(name)
            uniq.append(name)
    return uniq


def match_region_chrom(region: str, contigs: Iterable[str]) -> tuple[str, str]:
    """Return (resolved_region, matched_chrom)."""
    text = region.strip()
    if not text:
        raise ValueError("empty region")
    if ":" not in text:
        chrom, rest = text, ""
    else:
        chrom, rest = text.split(":", 1)
    contig_set = set(contigs)
    for alias in contig_aliases(chrom):
        if alias in contig_set:
            resolved = alias if not rest else f"{alias}:{rest}"
            return resolved, alias
    sample = ", ".join(list(contigs)[:12])
    raise SystemExit(
        f"region chrom {chrom!r} is not in the VCF contig list. "
        f"Have: {sample}{' …' if len(contig_set) > 12 else ''}. "
        f"Point gt_vcf/ref_vcf/map_file at the same chromosome as region, "
        f"or set region to a window on that contig."
    )


def parse_contigs_from_header(text: str) -> list[str]:
    contigs: list[str] = []
    for line in text.splitlines():
        if not line.startswith("##contig=<"):
            continue
        match = re.search(r"ID=([^,>]+)", line)
        if match:
            contigs.append(match.group(1))
    return contigs


def main(argv: Optional[list[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--region", required=True)
    p.add_argument("--header", help="VCF header text file (bcftools view -h)")
    p.add_argument("--contigs", help="Comma-separated contig names (tests)")
    p.add_argument("--out", help="Write resolved region here")
    args = p.parse_args(argv)
    if args.header:
        contigs = parse_contigs_from_header(open(args.header).read())
    elif args.contigs:
        contigs = [c for c in args.contigs.split(",") if c]
    else:
        raise SystemExit("need --header or --contigs")
    if not contigs:
        raise SystemExit("no ##contig=<ID=...> lines in header")
    resolved, matched = match_region_chrom(args.region, contigs)
    if args.out:
        open(args.out, "w").write(resolved + "\n")
    print(f"region {args.region} -> {resolved} (contig {matched})", file=sys.stderr)
    print(resolved)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
