#!/usr/bin/env python3
"""
Pure Length Variability Index (PLVI) per tandem-repeat locus from per-sample TRGT VCFs.

Danzi et al. (bioRxiv 2025.01.06.631535): the longest pure segment (LPS) of an
allele is its longest uninterrupted run of the locus motif, and PLVI is the
standard deviation of LPS length across haplotypes, ranked within motif-length
groups. Here LPS is the longest exact run of whole motif copies in the allele
sequence (any rotation of the motif; only the motif as written for motifs
longer than ``MAX_ROTATE_BP``), in bp. Loci with more than one catalog motif are
skipped.

Subcommands::

    catalog  one TRGT VCF  -> BED: chrom start end TRID motifs
    stats    TRGT VCFs     -> per-locus n, sum and sum of squares of LPS (one shard)
    merge    shard stats   -> PLVI table with within-group percentile ranks

Shards come from VCFs on the same catalog, so their rows line up by TRID.
"""

from __future__ import annotations

import argparse
import gzip
import math
import re
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
from functools import lru_cache
from pathlib import Path
from typing import Iterator, Optional, TextIO

import numpy as np

MAX_ROTATE_BP = 12
QUERY = "%CHROM\t%POS\t%INFO/END\t%INFO/TRID\t%INFO/MOTIFS\t%REF\t%ALT[\t%GT]\n"
STATS_COLUMNS = ("trid", "chrom", "start", "end", "motif", "n_hap", "lps_sum", "lps_sumsq")
PLVI_COLUMNS = (
    "trid", "chrom", "start", "end", "motif", "motif_len", "motif_group", "n_hap",
    "lps_mean", "plvi", "plvi_pct_in_group", "n_in_group",
)


def open_text(path: str, mode: str = "rt") -> TextIO:
    return gzip.open(path, mode) if path.endswith((".gz", ".bgz")) else open(path, mode)


@lru_cache(maxsize=100_000)
def motif_patterns(motif: str) -> tuple[re.Pattern, ...]:
    m = motif.upper()
    rots = {m[i:] + m[:i] for i in range(len(m))} if len(m) <= MAX_ROTATE_BP else {m}
    return tuple(re.compile(f"(?:{re.escape(r)})+") for r in sorted(rots))


def lps_bp(seq: str, motif: str) -> int:
    """Longest run of whole motif copies in seq, in bp."""
    if not motif or not seq:
        return 0
    s = seq.upper()
    best = 0
    for pat in motif_patterns(motif):
        for m in pat.finditer(s):
            n = m.end() - m.start()
            if n > best:
                best = n
    return best


def motif_group(k: int) -> str:
    if k <= 6:
        return str(k)
    if k <= 12:
        return "7-12"
    if k <= 24:
        return "13-24"
    return "25+"


def query_lines(vcf: str) -> Iterator[list[str]]:
    proc = subprocess.Popen(["bcftools", "query", "-f", QUERY, vcf], stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    assert proc.stdout is not None
    for line in proc.stdout:
        yield line.rstrip("\n").split("\t")
    proc.stdout.close()
    if proc.wait() != 0:
        raise SystemExit(f"bcftools query failed on {vcf}")


def write_catalog(vcf: str, out: str) -> int:
    n = 0
    with open_text(out, "wt") as fh:
        for f in query_lines(vcf):
            fh.write(f"{f[0]}\t{int(f[1]) - 1}\t{f[2]}\t{f[3]}\t{f[4]}\n")
            n += 1
    return n


def sample_lps(vcf: str, catalog_trids: list[str]) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    """(n_hap, sum, sumsq) per catalog row for one sample, plus rows looked up off-order."""
    n = len(catalog_trids)
    cnt = np.zeros(n, np.int32)
    s = np.zeros(n, np.float64)
    ss = np.zeros(n, np.float64)
    index: Optional[dict[str, int]] = None
    off_order = 0
    for i, f in enumerate(query_lines(vcf)):
        trid, motifs, ref, alt, gt = f[3], f[4], f[5], f[6], f[7]
        if "," in motifs or gt in (".", "./.", ".|."):
            continue
        j = i
        if i >= n or catalog_trids[i] != trid:
            if index is None:
                index = {t: k for k, t in enumerate(catalog_trids)}
            j = index.get(trid, -1)
            off_order += 1
            if j < 0:
                continue
        alleles = [ref] + (alt.split(",") if alt not in (".", "") else [])
        for a in re.split(r"[/|]", gt):
            if not a.isdigit():
                continue
            k = int(a)
            if k >= len(alleles):
                continue
            v = lps_bp(alleles[k], motifs)
            cnt[j] += 1
            s[j] += v
            ss[j] += v * v
    return cnt, s, ss, off_order


_CATALOG: list[str] = []


def _init(trids: list[str]) -> None:
    global _CATALOG
    _CATALOG = trids


def _one(vcf: str):
    return vcf, sample_lps(vcf, _CATALOG)


def read_catalog(path: str) -> list[tuple[str, str, str, str, str]]:
    rows = []
    with open_text(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            rows.append((f[3], f[0], f[1], f[2], f[4]))
    return rows


def cmd_stats(args: argparse.Namespace) -> int:
    vcfs = [ln.strip() for ln in Path(args.vcf_list).read_text().splitlines() if ln.strip()] if args.vcf_list else []
    vcfs += args.vcf or []
    catalog = read_catalog(args.catalog_bed)
    trids = [r[0] for r in catalog]
    n = len(trids)
    cnt = np.zeros(n, np.int64)
    s = np.zeros(n, np.float64)
    ss = np.zeros(n, np.float64)
    done = 0
    with ProcessPoolExecutor(max_workers=args.threads, initializer=_init, initargs=(trids,)) as pool:
        for vcf, (c1, s1, ss1, off) in pool.map(_one, vcfs):
            cnt += c1
            s += s1
            ss += ss1
            done += 1
            print(f"[plvi] {done}/{len(vcfs)} {Path(vcf).name} off_order={off}", file=sys.stderr, flush=True)
    with open_text(args.out, "wt") as fh:
        fh.write("\t".join(STATS_COLUMNS) + "\n")
        for (trid, chrom, start, end, motif), c, a, b in zip(catalog, cnt, s, ss):
            fh.write(f"{trid}\t{chrom}\t{start}\t{end}\t{motif}\t{int(c)}\t{a:.0f}\t{b:.0f}\n")
    print(f"[plvi] wrote {n:,} loci from {done} samples -> {args.out}", file=sys.stderr)
    return 0


def cmd_merge(args: argparse.Namespace) -> int:
    handles = [open_text(p) for p in args.stats]
    heads = [h.readline().rstrip("\n").split("\t") for h in handles]
    if any(hd != list(STATS_COLUMNS) for hd in heads):
        raise SystemExit("unexpected stats header")
    rows = []
    for lines in zip(*handles):
        parts = [ln.rstrip("\n").split("\t") for ln in lines]
        trid = parts[0][0]
        if any(p[0] != trid for p in parts):
            raise SystemExit(f"shards are out of order at {trid}; build them from one catalog")
        c = sum(int(p[5]) for p in parts)
        a = sum(float(p[6]) for p in parts)
        b = sum(float(p[7]) for p in parts)
        rows.append((parts[0][:5], c, a, b))
    for h in handles:
        h.close()
    max_hap = max((c for _r, c, _a, _b in rows), default=0)
    min_hap = math.ceil(args.min_call_rate * max_hap)
    out_rows = []
    for (trid, chrom, start, end, motif), c, a, b in rows:
        if "," in motif or c < max(min_hap, 2):
            continue
        mean = a / c
        var = max(b / c - mean * mean, 0.0)
        k = len(motif)
        out_rows.append([trid, chrom, start, end, motif, k, motif_group(k), c, mean, math.sqrt(var)])
    groups: dict[str, list[int]] = {}
    for i, r in enumerate(out_rows):
        groups.setdefault(r[6], []).append(i)
    for idx in groups.values():
        vals = np.array([out_rows[i][9] for i in idx])
        order = np.argsort(vals, kind="mergesort")
        ranks = np.empty(len(idx))
        ranks[order] = np.arange(1, len(idx) + 1)
        # Ties share the highest rank so equal PLVI gives equal percentile.
        uniq, inv = np.unique(vals, return_inverse=True)
        top = np.zeros(uniq.size)
        np.maximum.at(top, inv, ranks)
        pct = 100.0 * top[inv] / len(idx)
        for i, p in zip(idx, pct):
            out_rows[i].extend([p, len(idx)])
    out_rows.sort(key=lambda r: (-r[10], -r[9]))
    with open_text(args.out, "wt") as fh:
        fh.write("\t".join(PLVI_COLUMNS) + "\n")
        for r in out_rows:
            fh.write("\t".join([*map(str, r[:8]), f"{r[8]:.4g}", f"{r[9]:.4g}", f"{r[10]:.4f}", str(r[11])]) + "\n")
    print(f"[plvi] {len(out_rows):,} ranked loci (min haplotypes {min_hap:,}) -> {args.out}", file=sys.stderr)
    return 0


def main(argv: Optional[list[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("catalog")
    c.add_argument("--vcf", required=True)
    c.add_argument("--out", required=True)
    s = sub.add_parser("stats")
    s.add_argument("--catalog-bed", required=True)
    s.add_argument("--vcf", action="append")
    s.add_argument("--vcf-list", default=None)
    s.add_argument("--threads", type=int, default=1)
    s.add_argument("--out", required=True)
    m = sub.add_parser("merge")
    m.add_argument("--stats", nargs="+", required=True)
    m.add_argument("--min-call-rate", type=float, default=0.8,
                   help="Rank only loci with at least this fraction of the best-called locus's haplotypes")
    m.add_argument("--out", required=True)
    args = p.parse_args(argv)
    if args.cmd == "catalog":
        n = write_catalog(args.vcf, args.out)
        print(f"[plvi] catalog: {n:,} loci -> {args.out}", file=sys.stderr)
        return 0
    return cmd_stats(args) if args.cmd == "stats" else cmd_merge(args)


if __name__ == "__main__":
    raise SystemExit(main())
