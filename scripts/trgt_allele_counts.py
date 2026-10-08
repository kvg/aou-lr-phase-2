#!/usr/bin/env python3
"""
Distinct TRGT allele lengths per catalog locus, for Figure 2C.

A haplotype's dosage is the TRGT allele length minus the reference span,
``AL - (END - POS + 1)``, in bp. Repeat units are that change divided by the
first motif's length. Sequences of the same length are one allele, including
a length that is not a whole number of motif copies. ``./.`` is missing, not
reference. A catalog locus absent from a VCF counts as two missing haplotypes.
A hemizygous call counts as one haplotype.

The allele table matches ``aggregate_repeat_loci.ALLELE_COLUMNS`` (``locus_id``
is the TRID) so ``panel_c_repeat_alleles.py`` can read it. No sample ids.

Subcommands::

    counts   one shard of per-sample TRGT VCFs  -> counts.npz
    merge    shard npz files + the catalog BED  -> alleles.tsv.gz and summary.json

Both phases use one TRExplorer catalog. Run the workflow once on each VCF list.
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Optional, TextIO

import numpy as np

QUERY = "%POS\t%INFO/END\t%INFO/TRID\t[%AL]\t[%GT]\n"
ALLELE_COLUMNS = (
    "locus_id", "chrom", "start", "end", "period", "motif", "n_motifs", "n_records",
    "n_hap", "n_hap_missing", "dosage_bp", "dosage_units", "n_hap_allele",
)
GT_SPLIT = re.compile(r"[/|]")


def open_text(path: str, mode: str = "rt") -> TextIO:
    return gzip.open(path, mode) if path.endswith((".gz", ".bgz")) else open(path, mode)


def fmt_units(v: float) -> str:
    return f"{v:.4g}" if v != int(v) else str(int(v))


def read_catalog(path: str) -> list[tuple[str, str, str, str, str]]:
    rows = []
    with open_text(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                raise SystemExit(f"catalog row needs chrom, start, end, TRID, motifs: {path}")
            rows.append((f[3], f[0], f[1], f[2], f[4]))
    if not rows:
        raise SystemExit(f"empty catalog: {path}")
    return rows


def haplotype_dosages(al: str, gt: str, ref_len: int) -> list[Optional[int]]:
    """Dosage in bp for each haplotype, or None when that haplotype is missing.

    TRGT ``AL`` is already one length per haplotype, in genotype order. ``GT``
    only says which of those haplotypes were called.
    """
    if gt in (".", "./.", ".|."):
        return [None] * (1 if gt == "." else 2)
    gts = GT_SPLIT.split(gt)
    als = al.split(",")
    if len(als) != len(gts):
        return [None] * len(gts)
    out: list[Optional[int]] = []
    for length, allele in zip(als, gts):
        length = length.strip()
        if allele == "." or length in (".", ""):
            out.append(None)
            continue
        try:
            out.append(int(length) - ref_len)
        except ValueError:
            out.append(None)
    return out


def query_lines(vcf: str):
    proc = subprocess.Popen(
        ["bcftools", "query", "-f", QUERY, vcf], stdout=subprocess.PIPE, text=True, bufsize=1 << 20,
    )
    assert proc.stdout is not None
    for line in proc.stdout:
        fields = line.rstrip("\n").split("\t")
        if len(fields) != 5:
            proc.kill()
            raise SystemExit(
                f"{vcf} produced {len(fields)} query columns; TrgtAlleleCounts expects one-sample TRGT VCFs"
            )
        yield fields
    proc.stdout.close()
    if proc.wait() != 0:
        raise SystemExit(f"bcftools query failed on {vcf}")


def reduce_sparse(idx: np.ndarray, dosage: np.ndarray, count: np.ndarray):
    if idx.size == 0:
        return (
            np.array([], np.int32),
            np.array([], np.int32),
            np.array([], np.int64),
        )
    order = np.lexsort((dosage, idx))
    idx = np.asarray(idx, np.int32)[order]
    dosage = np.asarray(dosage, np.int32)[order]
    count = np.asarray(count, np.int64)[order]
    change = np.empty(idx.size, dtype=bool)
    change[0] = True
    change[1:] = (idx[1:] != idx[:-1]) | (dosage[1:] != dosage[:-1])
    starts = np.flatnonzero(change)
    return idx[starts], dosage[starts], np.add.reduceat(count, starts)


def accumulate(acc, idx: np.ndarray, dosage: np.ndarray, count: np.ndarray):
    if acc is None or acc[0].size == 0:
        return reduce_sparse(idx, dosage, count)
    if idx.size == 0:
        return acc
    return reduce_sparse(
        np.concatenate([acc[0], np.asarray(idx, np.int32)]),
        np.concatenate([acc[1], np.asarray(dosage, np.int32)]),
        np.concatenate([acc[2], np.asarray(count, np.int64)]),
    )


def sample_counts(vcf: str, catalog_trids: list[str]):
    """Called and missing haplotypes, plus non-reference (locus, dosage) counts, for one VCF."""
    n = len(catalog_trids)
    n_hap = np.zeros(n, np.int32)
    n_miss = np.zeros(n, np.int32)
    seen = np.zeros(n, dtype=bool)
    obs_i: list[int] = []
    obs_d: list[int] = []
    index: Optional[dict[str, int]] = None
    off_order = off_catalog = dup = 0
    for i, (pos, end, trid, al, gt) in enumerate(query_lines(vcf)):
        j = i
        if i >= n or catalog_trids[i] != trid:
            if index is None:
                index = {t: k for k, t in enumerate(catalog_trids)}
            j = index.get(trid, -1)
            off_order += 1
        if j < 0:
            off_catalog += 1
            continue
        if seen[j]:
            dup += 1
            continue
        seen[j] = True
        try:
            ref_len = int(end) - int(pos) + 1
        except ValueError:
            ref_len = 0
        if ref_len <= 0:
            n_miss[j] += 2
            continue
        for dosage in haplotype_dosages(al, gt, ref_len):
            if dosage is None:
                n_miss[j] += 1
            else:
                n_hap[j] += 1
                if dosage != 0:
                    obs_i.append(j)
                    obs_d.append(dosage)
    # A locus this VCF does not contain was not genotyped. Diploid missing keeps
    # it from being treated as reference. Hemizygous loci that are present are
    # already counted as one haplotype above.
    n_miss[~seen] += 2
    idx = np.asarray(obs_i, np.int32)
    dosage = np.asarray(obs_d, np.int32)
    if idx.size:
        idx, dosage, count = reduce_sparse(idx, dosage, np.ones(idx.size, np.int64))
    else:
        count = np.array([], np.int64)
    return n_hap, n_miss, idx, dosage, count, off_order, off_catalog, dup


_CATALOG: list[str] = []


def _init(trids: list[str]) -> None:
    global _CATALOG
    _CATALOG = trids


def _one(vcf: str):
    return vcf, sample_counts(vcf, _CATALOG)


def write_counts(path: str, n_hap, n_miss, acc, n_vcf: int, n_loci: int) -> None:
    idx, dosage, count = acc if acc is not None else (
        np.array([], np.int32), np.array([], np.int32), np.array([], np.int64),
    )
    np.savez_compressed(
        path, n_hap=np.asarray(n_hap, np.int64), n_miss=np.asarray(n_miss, np.int64),
        idx=np.asarray(idx, np.int32), dosage=np.asarray(dosage, np.int32),
        count=np.asarray(count, np.int64), n_vcf=np.int32(n_vcf), n_loci=np.int32(n_loci),
    )


def vcf_paths_from_lines(lines) -> list[str]:
    """One VCF path per line, or ``sample_id<TAB>path`` as in ``trgt_table.txt``."""
    out = []
    for line in lines:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if "\t" in line:
            sample, path = line.split("\t", 1)
            path = path.strip()
            if sample in ("sample_id", "sample") or not path:
                continue
            line = path
        out.append(line)
    return out


def cmd_counts(args: argparse.Namespace) -> int:
    vcfs = vcf_paths_from_lines(Path(args.vcf_list).read_text().splitlines()) if args.vcf_list else []
    vcfs += args.vcf or []
    if not vcfs:
        raise SystemExit("give --vcf-list or --vcf")
    catalog = read_catalog(args.catalog_bed)
    trids = [r[0] for r in catalog]
    n = len(trids)
    n_hap = np.zeros(n, np.int64)
    n_miss = np.zeros(n, np.int64)
    acc = None
    done = 0
    with ProcessPoolExecutor(max_workers=args.threads, initializer=_init, initargs=(trids,)) as pool:
        for vcf, (c_hap, c_miss, idx, dosage, count, off, off_cat, dup) in pool.map(_one, vcfs):
            n_hap += c_hap
            n_miss += c_miss
            acc = accumulate(acc, idx, dosage, count)
            done += 1
            print(
                f"[trgt-alleles] {done}/{len(vcfs)} {Path(vcf).name} "
                f"off_order={off} off_catalog={off_cat} dup_trid={dup}",
                file=sys.stderr, flush=True,
            )
    write_counts(args.out, n_hap, n_miss, acc, done, n)
    print(f"[trgt-alleles] {n:,} loci from {done} samples -> {args.out}", file=sys.stderr)
    return 0


def _saved(z, key: str, dtype) -> np.ndarray:
    # Copy out of the compressed archive before it closes. A view would dangle.
    return np.array(z[key], dtype=dtype, copy=True)


def load_shard(path: str, n_loci: int):
    with np.load(path) as z:
        if int(z["n_loci"]) != n_loci:
            raise SystemExit(f"{path} has {int(z['n_loci'])} loci; catalog has {n_loci}. Use one catalog for every shard.")
        n_vcf = int(z["n_vcf"])
        return (
            _saved(z, "n_hap", np.int64), _saved(z, "n_miss", np.int64),
            _saved(z, "idx", np.int32), _saved(z, "dosage", np.int32), _saved(z, "count", np.int64),
            n_vcf,
        )


def locus_slices(idx: np.ndarray):
    if idx.size == 0:
        return np.array([], np.int64), np.array([], np.int64)
    change = np.empty(idx.size, dtype=bool)
    change[0] = True
    change[1:] = idx[1:] != idx[:-1]
    starts = np.flatnonzero(change)
    ends = np.empty_like(starts)
    ends[:-1] = starts[1:]
    ends[-1] = idx.size
    return starts, ends


def write_alleles(catalog, n_hap, n_miss, idx, dosage, count, out_path: str) -> dict:
    n = len(catalog)
    sums = np.zeros(n, np.int64)
    if idx.size:
        np.add.at(sums, idx, count)
    ref = n_hap - sums
    if np.any(ref < 0):
        bad = int(np.flatnonzero(ref < 0)[0])
        raise SystemExit(f"allele counts exceed called haplotypes at {catalog[bad][0]}")
    starts, ends = locus_slices(idx)
    g = 0
    summary = {
        "n_loci_catalog": n,
        "hap_called": int(n_hap.sum()),
        "hap_missing": int(n_miss.sum()),
        "loci_with_calls": int(np.count_nonzero(n_hap)),
        "loci_polymorphic": 0,
        "alleles": 0,
        "alleles_nonref": 0,
    }
    with open_text(out_path, "wt") as fh:
        fh.write("\t".join(ALLELE_COLUMNS) + "\n")
        for i, (trid, chrom, start, end, motif) in enumerate(catalog):
            nh, nm = int(n_hap[i]), int(n_miss[i])
            if nh == 0 and nm == 0:
                continue
            first = motif.split(",")[0] if motif else ""
            period = max(len(first), 1)
            n_motifs = motif.count(",") + 1 if motif else 0
            head = f"{trid}\t{chrom}\t{start}\t{end}\t{period}\t{motif or '.'}\t{n_motifs}\t1\t{nh}\t{nm}"
            alleles: list[tuple[int, int]] = []
            if int(ref[i]) > 0:
                alleles.append((0, int(ref[i])))
            if g < starts.size and int(idx[starts[g]]) == i:
                for d, c in zip(dosage[starts[g]:ends[g]], count[starts[g]:ends[g]]):
                    alleles.append((int(d), int(c)))
                g += 1
            elif g < starts.size and int(idx[starts[g]]) < i:
                raise SystemExit(f"sparse counts skipped locus index {int(idx[starts[g]])}")
            if nh == 0:
                alleles = [(0, 0)]
            alleles.sort()
            called = [(d, c) for d, c in alleles if c > 0]
            if len(called) > 1:
                summary["loci_polymorphic"] += 1
            for d, c in alleles:
                if c == 0 and nh > 0:
                    continue
                summary["alleles"] += 1
                if d != 0 and c > 0:
                    summary["alleles_nonref"] += 1
                fh.write(f"{head}\t{d}\t{fmt_units(d / period)}\t{c}\n")
    if g != starts.size:
        raise SystemExit("sparse counts left over after the catalog")
    return summary


def cmd_merge(args: argparse.Namespace) -> int:
    catalog = read_catalog(args.catalog_bed)
    n = len(catalog)
    n_hap = np.zeros(n, np.int64)
    n_miss = np.zeros(n, np.int64)
    acc = None
    n_vcf = 0
    for path in args.counts:
        c_hap, c_miss, idx, dosage, count, n_one = load_shard(path, n)
        n_hap += c_hap
        n_miss += c_miss
        acc = accumulate(acc, idx, dosage, count)
        n_vcf += n_one
        print(f"[trgt-alleles] merged {path} ({n_one} samples)", file=sys.stderr, flush=True)
    idx, dosage, count = acc if acc is not None else (
        np.array([], np.int32), np.array([], np.int32), np.array([], np.int64),
    )
    summary = write_alleles(catalog, n_hap, n_miss, idx, dosage, count, args.out)
    summary["n_vcf"] = n_vcf
    Path(args.out_summary).write_text(json.dumps(summary, indent=2) + "\n")
    print(
        f"[trgt-alleles] {summary['loci_with_calls']:,} loci with calls, "
        f"{summary['loci_polymorphic']:,} polymorphic, {n_vcf} samples -> {args.out}",
        file=sys.stderr,
    )
    return 0


def main(argv: Optional[list[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("counts")
    s.add_argument("--catalog-bed", required=True)
    s.add_argument("--vcf", action="append")
    s.add_argument("--vcf-list", default=None)
    s.add_argument("--threads", type=int, default=1)
    s.add_argument("--out", required=True, help="counts.npz for this shard")
    m = sub.add_parser("merge")
    m.add_argument("--catalog-bed", required=True)
    m.add_argument("--counts", nargs="+", required=True)
    m.add_argument("--out", required=True, help="allele table (.tsv or .tsv.gz)")
    m.add_argument("--out-summary", required=True)
    args = p.parse_args(argv)
    return cmd_counts(args) if args.cmd == "counts" else cmd_merge(args)


if __name__ == "__main__":
    raise SystemExit(main())
