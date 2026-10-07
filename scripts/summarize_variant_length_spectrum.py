#!/usr/bin/env python3
"""Count alleles by signed length (1 bp) and sequence context for Figure 2.

Uses the same sources and rules as Table 2 (tab:callset):

  * ``--glnexus-vcf`` (one DeepVariant + GLnexus shard per chromosome,
    FILTER PASS or ``.``): SNVs (signed length 0) and indels with
    |L| < ``--small-max-bp`` (20). Equal-length multi-base substitutions
    (MNPs), spanning ``*`` and symbolic alleles are counted in the summary only.
  * ``--main-vcf`` (v3 integrated SV main partition): resolved DEL / INS with
    |SVLEN| >= ``--sv-min-bp`` (20), i.e. ``size_bin_ge20`` in the site table.
    SVTYPE / SVLEN / END are derived as in ``extract_sites_from_vcf.py``.
  * ``--large-vcf``: ultralong companion, every record.
  * ``--bnd-vcf``: breakends; no length (``signed_len`` 0). The class comes from
    either end, the record's own position or its mate's (ALT bracket notation,
    else INFO CHR2 / END), with the same SR > SD > RM precedence.

GLnexus is per chromosome; the companions are genome-wide and are scanned once
(by MergeVariantLengthSpectrum), optionally restricted with ``--contigs``.

Every row carries ``region_class`` (US / RM / SD / SR) from
``annotate_repeat_context.context_label`` with the merged RM / SR / SD tracks:
either breakpoint assigns the class (SR > SD > RM), DEL/DUP/CNV > 5 kb use body
coverage > 0.5, insertions collapse to the anchor. Small variants use the same
breakpoint rule on their REF span.

Output rows: ``chrom  partition  signed_len  region_class  n_sites`` with
partition in {small, sv, ultralong, bnd}.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from bisect import bisect_right
from collections import Counter
from pathlib import Path
from typing import Iterable, Iterator, Optional

sys.path.insert(0, str(Path(__file__).resolve().parent))

from annotate_repeat_context import (  # noqa: E402
    TRACK_FIELDS,
    context_label,
    load_intervals,
    reference_interval,
)
from sv_site_utils import (  # noqa: E402
    bcftools_bin,
    bcftools_stdout,
    infer_svlen,
    infer_svtype,
    iter_site_records,
    open_text,
    parse_info,
    size_flags,
)

RESOLVED = ("DEL", "DUP", "INS", "INV")
LENGTH_SVTYPES = frozenset({"DEL", "INS"})
TRACK_NAMES = tuple(name for name, _field in TRACK_FIELDS)
SKIP_SMALL_ALTS = frozenset({"*", "<*>", "<NON_REF>", "."})
BINS_HEADER = "chrom\tpartition\tsigned_len\tregion_class\tn_sites\n"


def chrom_contig(chrom_label: str) -> str:
    """Terra entity id ``…v1.chr1`` → contig ``chr1``; plain ``chr1`` unchanged."""
    label = chrom_label.strip().strip("\"'")
    if "." in label:
        return label.rsplit(".", 1)[-1]
    return label


def chrom_matches(record_chrom: str, want_contig: str) -> bool:
    a = record_chrom.strip()
    b = want_contig.strip()
    if a == b:
        return True
    a_core = a[3:] if a.lower().startswith("chr") else a
    b_core = b[3:] if b.lower().startswith("chr") else b
    return a_core == b_core and a_core != ""


class Tracks:
    """Merged RM / SD / SR intervals with a fast point lookup for small variants."""

    def __init__(self, beds: dict[str, str]) -> None:
        self.intervals = {name: load_intervals(path) for name, path in beds.items()}
        self._points: dict[str, list[tuple[list[int], list[int]]]] = {}
        self._chrom: dict[str, dict[str, list[tuple[int, int]]]] = {}

    def for_chrom(self, chrom: str) -> dict[str, list[tuple[int, int]]]:
        hit = self._chrom.get(chrom)
        if hit is None:
            hit = {name: self.intervals[name].get(chrom, []) for name in TRACK_NAMES}
            self._chrom[chrom] = hit
        return hit

    def points(self, chrom: str) -> list[tuple[list[int], list[int]]]:
        """Per track (RM, SD, SR order): (starts, ends) of merged intervals."""
        hit = self._points.get(chrom)
        if hit is None:
            hit = []
            for name in TRACK_NAMES:
                iv = self.intervals[name].get(chrom, [])
                hit.append(([s for s, _e in iv], [e for _s, e in iv]))
            self._points[chrom] = hit
        return hit


def breakpoint_class(points: list[tuple[list[int], list[int]]], start0: int, end0: int) -> str:
    """Same result as ``context_label`` on the breakpoint path, via bisect."""
    label = "US"
    for name, (starts, ends) in zip(TRACK_NAMES, points):
        i = bisect_right(starts, start0) - 1
        if i >= 0 and start0 < ends[i]:
            label = name
            continue
        i = bisect_right(starts, end0) - 1
        if i >= 0 and end0 < ends[i]:
            label = name
    return label


def small_allele_length(ref: str, alt: str) -> tuple[str, Optional[int]]:
    """Classify one GLnexus ALT: (kind, signed_len). kind in snv/indel/mnp/ref/other."""
    if alt in SKIP_SMALL_ALTS or alt.startswith("<") or "[" in alt or "]" in alt:
        return "other", None
    if len(alt) == len(ref):
        diff = sum(1 for a, b in zip(ref, alt) if a != b)
        if diff == 0:
            return "ref", None
        if diff == 1:
            return "snv", 0
        return "mnp", None
    return "indel", len(alt) - len(ref)


def iter_small_records(path: str) -> Iterator[tuple[str, str, str, str, str]]:
    """Yield (chrom, pos, ref, alt_field, filter) without INFO or genotypes."""
    if bcftools_bin():
        fmt = r"%CHROM\t%POS\t%REF\t%ALT\t%FILTER\n"
        with bcftools_stdout(["query", "-f", fmt], path) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) == 5:
                    yield parts[0], parts[1], parts[2], parts[3], parts[4]
        return
    with open_text(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t", 7)
            if len(parts) >= 7:
                yield parts[0], parts[1], parts[3], parts[4], parts[6]


def summarize_glnexus(
    path: str,
    *,
    contig: Optional[str],
    contigs: Optional[set[str]],
    tracks: Tracks,
    small_max_bp: int,
    filters: frozenset[str],
    counts: Counter,
    totals: Counter,
) -> None:
    points_by_chrom: dict[str, list] = {}
    n = 0
    for chrom, pos_s, ref, alt_field, filt in iter_small_records(path):
        n += 1
        if n % 2_000_000 == 0:
            print(f"[length-spectrum] {path}: {n:,} records", file=sys.stderr, flush=True)
        off = (contig is not None and not chrom_matches(chrom, contig)) or (contigs is not None and chrom not in contigs)
        if off:
            totals["glnexus_records_off_contig"] += 1
            continue
        totals["glnexus_records"] += 1
        if filt not in filters:
            totals["glnexus_records_filtered"] += 1
            continue
        totals["glnexus_records_pass"] += 1
        points = points_by_chrom.get(chrom)
        if points is None:
            points = tracks.points(chrom)
            points_by_chrom[chrom] = points
        start0 = int(pos_s) - 1
        end0 = start0 + len(ref)
        if end0 <= start0:
            end0 = start0 + 1
        label: Optional[str] = None
        has_snv = has_indel = False
        for alt in alt_field.split(","):
            kind, slen = small_allele_length(ref, alt)
            totals[f"glnexus_alleles_{kind}"] += 1
            if slen is None:
                continue
            if kind == "snv":
                has_snv = True
            else:
                has_indel = True
                if abs(slen) >= small_max_bp:
                    totals["glnexus_indel_ge_small_max_dropped"] += 1
                    continue
                totals["ins_lt_small_max" if slen > 0 else "del_lt_small_max"] += 1
            if label is None:
                label = breakpoint_class(points, start0, end0)
            counts[("small", slen, label)] += 1
        if has_snv:
            totals["snv_records"] += 1
        if has_indel:
            totals["indel_records"] += 1


def sv_site_fields(
    pos: int, ref: str, alt: str, info: dict[str, str], source_vcf: str
) -> tuple[str, Optional[int], int]:
    """(svtype, svlen, end) exactly as ``extract_sites_from_vcf.py`` writes them."""
    svtype = infer_svtype(alt, info, source_vcf)
    svlen = infer_svlen(ref, alt, info)
    end: Optional[int] = None
    if svlen is None and "END" in info:
        try:
            end = int(info["END"])
            svlen = end - pos
            if svtype == "DEL":
                svlen = -abs(svlen)
        except ValueError:
            end = None
    if end is None and "END" in info:
        try:
            end = int(info["END"])
        except ValueError:
            end = pos
    if end is None:
        end = pos + abs(svlen or 0)
    return svtype, svlen, end


def signed_sv_length(svtype: str, svlen: int) -> int:
    if svtype == "DEL":
        return -abs(svlen)
    if svtype in ("INS", "DUP"):
        return abs(svlen)
    return int(svlen)


BND_MATE = re.compile(r"[\[\]]([^\[\]:]+):(\d+)[\[\]]")
CLASS_RANK = {"US": 0, **{name: i + 1 for i, name in enumerate(TRACK_NAMES)}}


def iter_bnd_records(path: str) -> Iterator[tuple[str, int, str, str, dict[str, str]]]:
    """Yield (chrom, pos, id, alt, info); unlike ``iter_site_records`` this keeps ALT."""
    if bcftools_bin():
        with bcftools_stdout(["query", "-f", r"%CHROM\t%POS\t%ID\t%ALT\t%INFO\n"], path) as fh:
            for line in fh:
                chrom, pos, vid, alt, info = line.rstrip("\n").split("\t", 4)
                yield chrom, int(pos), vid, alt, parse_info(info)
        return
    with open_text(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t", 8)
            if len(parts) >= 8:
                yield parts[0], int(parts[1]), parts[2], parts[4], parse_info(parts[7])


def bnd_mate(alt: str, info: dict[str, str]) -> Optional[tuple[str, int]]:
    m = BND_MATE.search(alt)
    if m:
        return m.group(1), int(m.group(2))
    if "CHR2" in info and "END" in info:
        try:
            return info["CHR2"], int(info["END"])
        except ValueError:
            return None
    return None


def summarize_bnd(path: str, *, contigs: Optional[set[str]], tracks: Tracks, counts: Counter, totals: Counter) -> None:
    """One count per breakend event: a reciprocal pair of records (A's mate is B, B's mate is A) counts once."""
    ids: set[str] = set()
    mate_ids: list[str] = []
    recs: list[tuple[tuple[str, int], Optional[tuple[str, int]], str]] = []
    for chrom, pos, vid, alt, info in iter_bnd_records(path):
        if contigs is not None and chrom not in contigs:
            totals["bnd_records_off_contig"] += 1
            continue
        totals["bnd_records"] += 1
        ends = [(chrom, pos)]
        mate = bnd_mate(alt, info)
        if mate is None:
            totals["bnd_mate_unparsed"] += 1
        else:
            ends.append(mate)
        label = "US"
        for c, p in ends:
            p0 = p - 1
            cls = breakpoint_class(tracks.points(c), p0, p0)
            if CLASS_RANK[cls] > CLASS_RANK[label]:
                label = cls
        recs.append(((chrom, pos), mate, label))
        if vid and vid != ".":
            ids.add(vid)
        if info.get("MATEID"):
            mate_ids.append(info["MATEID"])
    links = {(own, mate) for own, mate, _label in recs if mate is not None}
    for own, mate, label in recs:
        if mate is not None and mate != own and (mate, own) in links and own > mate:
            totals["bnd_reciprocal_pairs"] += 1
            continue
        totals["bnd_events"] += 1
        counts[("bnd", 0, label)] += 1
    totals["bnd_records_with_mateid"] = len(mate_ids)
    totals["bnd_records_mate_in_file"] = sum(1 for m in mate_ids if m in ids)


def summarize_companion(
    path: str,
    *,
    source_vcf: str,
    contigs: Optional[set[str]],
    tracks: Tracks,
    sv_min_bp: int,
    counts: Counter,
    totals: Counter,
) -> None:
    n = 0
    for chrom, pos_s, _vid, ref, alt, _filt, info_s in iter_site_records(path):
        n += 1
        if n % 500_000 == 0:
            print(f"[length-spectrum] {path}: {n:,} records", file=sys.stderr, flush=True)
        if contigs is not None and chrom not in contigs:
            totals[f"{source_vcf}_records_off_contig"] += 1
            continue
        totals[f"{source_vcf}_records"] += 1
        pos = int(pos_s)
        info = parse_info(info_s)
        svtype, svlen, end = sv_site_fields(pos, ref, alt, info, source_vcf)
        svtype_u = svtype.upper()
        if svlen is None:
            totals[f"{source_vcf}_unknown_length"] += 1
            continue
        if source_vcf == "main":
            if svtype not in RESOLVED:
                totals[f"main_unresolved_{svtype}"] += 1
                continue
            ge20, ge50 = size_flags(svlen, svtype)
            if ge20:
                totals[f"main_ge20_{svtype}"] += 1
            if ge50:
                totals[f"main_ge50_{svtype}"] += 1
            if svtype_u not in LENGTH_SVTYPES:
                continue
            if abs(svlen) < sv_min_bp:
                totals["main_lt_sv_min_dropped"] += 1
                continue
            partition = "sv"
        else:
            partition = "ultralong"
        site = {"pos": pos, "end": end, "svtype": svtype, "svlen": svlen}
        start0, end0 = reference_interval(site)
        label, _hits = context_label(svtype, start0, end0, tracks.for_chrom(chrom))
        counts[(partition, signed_sv_length(svtype_u, svlen), label)] += 1


def write_bins(path: Path, *, chrom: str, counts: Counter) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wt", encoding="utf-8") as out:
        out.write(BINS_HEADER)
        for (part, slen, label), n in sorted(counts.items()):
            out.write(f"{chrom}\t{part}\t{slen}\t{label}\t{n}\n")


def main(argv: Optional[Iterable[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--glnexus-vcf", help="One DeepVariant + GLnexus chromosome shard")
    p.add_argument("--chrom", help="Shard label (Terra GL_INTERVAL_set id), or 'all' for a whole-genome VCF; required with --glnexus-vcf")
    p.add_argument("--main-vcf", help="v3 integrated SV main partition (genome-wide)")
    p.add_argument("--large-vcf", help="v3 ultralong partition (genome-wide)")
    p.add_argument("--bnd-vcf", help="v3 breakend partition (genome-wide)")
    p.add_argument("--contigs", help="Comma-separated contigs to keep from all inputs (default: all)")
    p.add_argument("--rmsk-bed", required=True)
    p.add_argument("--simple-repeat-bed", required=True)
    p.add_argument("--segdup-bed", required=True)
    p.add_argument("--small-max-bp", type=int, default=20, help="GLnexus indels kept when |L| < this")
    p.add_argument("--sv-min-bp", type=int, default=20, help="Main SVs kept when |SVLEN| >= this")
    p.add_argument("--apply-filters", default="PASS,.", help="GLnexus FILTER values kept (bcftools stats -f)")
    p.add_argument("--out-bins", required=True)
    p.add_argument("--out-summary", required=True)
    args = p.parse_args(list(argv) if argv is not None else None)

    companions = [(s, getattr(args, f"{s}_vcf")) for s in ("main", "large", "bnd")]
    companions = [(s, v) for s, v in companions if v]
    if not args.glnexus_vcf and not companions:
        p.error("give --glnexus-vcf and/or at least one of --main-vcf / --large-vcf / --bnd-vcf")
    if args.glnexus_vcf and not args.chrom:
        p.error("--chrom is required with --glnexus-vcf")

    tracks = Tracks({"RM": args.rmsk_bed, "SD": args.segdup_bed, "SR": args.simple_repeat_bed})
    counts: Counter = Counter()
    totals: Counter = Counter()
    contig = chrom_contig(args.chrom) if args.chrom and args.chrom.lower() != "all" else None
    keep = {c.strip() for c in args.contigs.split(",") if c.strip()} if args.contigs else None

    if args.glnexus_vcf:
        summarize_glnexus(
            args.glnexus_vcf,
            contig=contig,
            contigs=keep,
            tracks=tracks,
            small_max_bp=args.small_max_bp,
            filters=frozenset(f.strip() for f in args.apply_filters.split(",")),
            counts=counts,
            totals=totals,
        )
    for source_vcf, path in companions:
        if source_vcf == "bnd":
            summarize_bnd(path, contigs=keep, tracks=tracks, counts=counts, totals=totals)
            continue
        summarize_companion(
            path,
            source_vcf=source_vcf,
            contigs=keep,
            tracks=tracks,
            sv_min_bp=args.sv_min_bp,
            counts=counts,
            totals=totals,
        )

    label = args.chrom if args.glnexus_vcf else "companions"
    write_bins(Path(args.out_bins), chrom=label, counts=counts)
    summary = {
        "chrom": label,
        "contig": contig,
        "binning": "1bp_signed_len",
        "region_rule": "annotate_repeat_context.context_label",
        "sources": {
            "glnexus": args.glnexus_vcf,
            **{s: v for s, v in companions},
        },
        "small_max_bp": args.small_max_bp,
        "sv_min_bp": args.sv_min_bp,
        "apply_filters": args.apply_filters,
        "n_length_rows": len(counts),
        "counts": {k: int(v) for k, v in sorted(totals.items())},
    }
    Path(args.out_summary).write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2), file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
