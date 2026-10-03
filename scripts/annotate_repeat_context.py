#!/usr/bin/env python3
"""
Assign each site one sequence-context class: US, RM, SD, or SR.

This is the breakpoint rule from Xuefang Zhao's
``annotate_genomic_context.sh`` (gatk-sv ``xz_fixes_3`` / Zhao et al., AJHG
2021), not a whole-interval union. SR overrides SD, which overrides RM.

* Either breakpoint can assign the class. Breakpoints are the 0-based BED
  start (VCF POS - 1) and the BED end (VCF END, 1-based inclusive).
* DEL, DUP, and CNV whose reference span is greater than 5 kb ignore
  breakpoints. A class applies only when merged-track coverage of the body
  is above 0.5, with the same priority. Otherwise the site is US.
* Insertions stay on the breakpoint rule. When the site table filled ``end``
  from SVLEN (``end - pos == |SVLEN|``), both breakpoints are the anchor at
  POS, not a span of the inserted length.
* ``hit_rmsk`` / ``hit_simpleRepeat`` / ``hit_genomicSuperDups`` record the
  evidence before that priority collapse, so more than one can be true.
* ``hit_cmrg`` is still any overlap of the reference span with the CMRG BED.
  It is not part of ``region_class``.

Track BEDs should be the merged GRCh38 files shipped with that script
(``hg38.RM/SD/SR.sorted.merged.bed.gz``). Intervals are merged here so an
unmerged track cannot inflate the body fraction.
"""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from sv_site_utils import iter_sites, open_text, site_row, write_site_header  # noqa: E402

BODY_SVTYPES = frozenset({"DEL", "DUP", "CNV"})
INSERTION_SVTYPES = frozenset({"INS", "ALU", "LINE1", "SVA", "MEI"})
BODY_SPAN_GT = 5000

TRACK_FIELDS = (
    ("RM", "hit_rmsk"),
    ("SD", "hit_genomicSuperDups"),
    ("SR", "hit_simpleRepeat"),
)


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
            if end == start:
                continue
            intervals[chrom].append((start, end))
    return {chrom: _merge_intervals(iv) for chrom, iv in intervals.items()}


def _merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if not intervals:
        return []
    ordered = sorted(intervals)
    merged = [ordered[0]]
    for start, end in ordered[1:]:
        prev_start, prev_end = merged[-1]
        if start <= prev_end:
            merged[-1] = (prev_start, max(prev_end, end))
        else:
            merged.append((start, end))
    return merged


def _first_overlapping(intervals: list[tuple[int, int]], start0: int) -> int:
    """Index of the first interval with end > start0."""
    lo, hi = 0, len(intervals)
    while lo < hi:
        mid = (lo + hi) // 2
        if intervals[mid][1] <= start0:
            lo = mid + 1
        else:
            hi = mid
    return lo


def contains_point(intervals: list[tuple[int, int]], x: int) -> bool:
    """True if some half-open interval contains the 0-based coordinate x."""
    i = _first_overlapping(intervals, x)
    if i >= len(intervals):
        return False
    start, end = intervals[i]
    return start <= x < end


def covered_bases(intervals: list[tuple[int, int]], start0: int, end0: int) -> int:
    """Unique bases of [start0, end0) covered by merged intervals."""
    if end0 <= start0:
        return 0
    i = _first_overlapping(intervals, start0)
    covered = 0
    while i < len(intervals) and intervals[i][0] < end0:
        start, end = intervals[i]
        covered += max(0, min(end, end0) - max(start, start0))
        i += 1
    return covered


def overlaps(intervals: list[tuple[int, int]], start0: int, end0: int) -> bool:
    if end0 <= start0:
        end0 = start0 + 1
    return covered_bases(intervals, start0, end0) > 0


def reference_interval(site: dict) -> tuple[int, int]:
    """0-based half-open reference span used for breakpoints and body coverage.

    Returns (start, end) with end > start. Insertions whose ``end`` was filled
    from SVLEN collapse to the 1 bp anchor at POS.
    """
    start0 = int(site["pos"]) - 1
    end_field = site.get("end")
    if end_field in (None, ""):
        end0 = start0
    else:
        end0 = int(end_field)
    if end0 < start0:
        start0, end0 = end0, start0

    svtype = str(site.get("svtype") or "").upper()
    svlen = site.get("svlen")
    if svtype in INSERTION_SVTYPES and svlen not in (None, ""):
        try:
            alen = abs(int(svlen))
        except (TypeError, ValueError):
            alen = 0
        if alen > 1 and end_field not in (None, "") and (int(end_field) - int(site["pos"])) == alen:
            end0 = start0

    if end0 <= start0:
        end0 = start0 + 1
    return start0, end0


def context_label(
    svtype: str,
    start0: int,
    end0: int,
    tracks: dict[str, list[tuple[int, int]]],
) -> tuple[str, dict[str, bool]]:
    """Return (region_class, per-track evidence) for one reference span.

    ``tracks`` keys are RM, SD, and SR. Evidence is a breakpoint hit, or body
    coverage above 0.5 for a long DEL/DUP/CNV. The class is the last evidence
    flag in RM, SD, SR order.
    """
    span = end0 - start0
    use_body = str(svtype or "").upper() in BODY_SVTYPES and span > BODY_SPAN_GT
    hits: dict[str, bool] = {}
    if use_body:
        for name, _field in TRACK_FIELDS:
            frac = covered_bases(tracks.get(name, []), start0, end0) / span
            hits[name] = frac > 0.5
    else:
        for name, _field in TRACK_FIELDS:
            iv = tracks.get(name, [])
            hits[name] = contains_point(iv, start0) or contains_point(iv, end0)
    label = "US"
    for name, _field in TRACK_FIELDS:
        if hits[name]:
            label = name
    return label, hits


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sites", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--rmsk-bed", required=True)
    p.add_argument("--simple-repeat-bed", required=True)
    p.add_argument("--segdup-bed", required=True)
    p.add_argument("--cmrg-bed", required=True)
    args = p.parse_args()

    tracks = {
        "RM": load_intervals(args.rmsk_bed),
        "SD": load_intervals(args.segdup_bed),
        "SR": load_intervals(args.simple_repeat_bed),
    }
    cmrg = load_intervals(args.cmrg_bed)

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with open_text(args.out, "wt") as out:
        write_site_header(out)
        for s in iter_sites(args.sites):
            n += 1
            if n % 500000 == 0:
                print(f"[annotate_repeat] {n:,} sites", file=sys.stderr, flush=True)
            start0, end0 = reference_interval(s)
            chrom = s["chrom"]
            chrom_tracks = {name: tracks[name].get(chrom, []) for name in tracks}
            label, hits = context_label(s.get("svtype") or "", start0, end0, chrom_tracks)
            s["region_class"] = label
            for name, field in TRACK_FIELDS:
                s[field] = hits[name]
            s["hit_cmrg"] = overlaps(cmrg.get(chrom, []), start0, end0)
            out.write(site_row(s))
    print(f"[annotate_repeat] wrote {n:,} sites", file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
