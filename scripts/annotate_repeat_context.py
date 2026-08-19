#!/usr/bin/env python3
"""
Annotate sites with repetitive-region union overlaps.

Inputs are BED files (possibly gzipped) with at least chrom, start, end
(0-based, half-open). Track names map to hit_* columns.
"""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from sv_site_utils import iter_sites, open_text, site_row, write_site_header  # noqa: E402


def load_intervals(bed_path: str) -> dict[str, list[tuple[int, int]]]:
    intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with open_text(bed_path, "rt") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#") or line.startswith("track"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            if end < start:
                start, end = end, start
            intervals[chrom].append((start, end))
    for chrom in intervals:
        intervals[chrom].sort()
    return intervals


def overlaps(intervals: list[tuple[int, int]], start0: int, end0: int) -> bool:
    """True if [start0, end0) overlaps any interval (binary search over starts)."""
    if end0 <= start0:
        end0 = start0 + 1
    lo, hi = 0, len(intervals)
    # find first interval with end > start0
    while lo < hi:
        mid = (lo + hi) // 2
        if intervals[mid][1] <= start0:
            lo = mid + 1
        else:
            hi = mid
    i = lo
    while i < len(intervals) and intervals[i][0] < end0:
        a, b = intervals[i]
        if a < end0 and b > start0:
            return True
        i += 1
    return False


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sites", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--rmsk-bed", required=True)
    p.add_argument("--simple-repeat-bed", required=True)
    p.add_argument("--segdup-bed", required=True)
    args = p.parse_args()

    tracks = {
        "hit_rmsk": load_intervals(args.rmsk_bed),
        "hit_simpleRepeat": load_intervals(args.simple_repeat_bed),
        "hit_genomicSuperDups": load_intervals(args.segdup_bed),
    }

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with open_text(args.out, "wt") as out:
        write_site_header(out)
        for s in iter_sites(args.sites):
            n += 1
            if n % 500000 == 0:
                print(f"[annotate_repeat] {n:,} sites", file=sys.stderr, flush=True)
            start0 = int(s["pos"]) - 1
            end0 = int(s["end"]) if s.get("end") not in (None, "") else start0 + 1
            if end0 <= start0:
                end0 = start0 + 1
            chrom = s["chrom"]
            any_hit = False
            for col, ivals in tracks.items():
                hit = overlaps(ivals.get(chrom, []), start0, end0)
                s[col] = hit
                any_hit = any_hit or hit
            s["region_class"] = "repetitive" if any_hit else "non_repetitive"
            out.write(site_row(s))
    print(f"[annotate_repeat] wrote {n:,} sites", file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
