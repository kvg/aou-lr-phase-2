#!/usr/bin/env python3
"""
Compare integrated-callset locus dosages with TRGT at the same catalog loci.

Integrated side: the locus VCF from ``aggregate_repeat_loci.py --out-vcf`` (GT
over ``<TR:+Nbp>`` alleles). Catalog loci absent from it had no length change
in the integrated callset, so every sample is REF-length there. TRGT side: one
per-sample TRGT VCF per sample; a haplotype's change is ``len(allele) -
len(REF)`` from the allele sequences, so padding conventions cancel.

Genotypes are compared as unordered pairs (TRGT calls are unphased). A pair
is ``exact`` when both changes match, ``within_unit`` when both differ by at
most one motif length, and ``overcount`` when the integrated total exceeds the
TRGT total in magnitude by >= 50 bp (the signature of an indel and an SV
record describing the same event).

Outputs hold no sample ids: ``--out-loci`` (one row per locus) and
``--out-summary`` (totals overall and by integrated record count and TRGT
change size).
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from typing import Optional, TextIO

OVERCOUNT_BP = 50
LOCI_COLUMNS = (
    "locus_id", "period", "n_records", "n_pairs", "n_exact", "n_within_unit", "n_overcount",
    "mean_abs_sum_diff_bp",
)


def open_text(path: str) -> TextIO:
    return gzip.open(path, "rt") if path.endswith((".gz", ".bgz")) else open(path)


def load_catalog(path: str, chrom: str, lo: Optional[int] = None, hi: Optional[int] = None) -> dict[str, int]:
    """locus id -> period (first motif length). lo/hi are 1-based positions on chrom."""
    out = {}
    with open_text(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0] != chrom:
                continue
            if lo is not None and not (lo <= int(f[1]) + 1 <= hi):
                continue
            out[f[3]] = max(len(f[4].split(",")[0]), 1) if len(f) > 4 else 1
    return out


def load_locus_vcf(path: str, keep: set[str]):
    """{locus: (n_records, {sample: (d1, d2) or None})} and the sample list."""
    loci: dict[str, tuple[int, dict[str, Optional[tuple[int, int]]]]] = {}
    samples: list[str] = []
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            f = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                samples = f[9:]
                continue
            bp = [0] + [int(re.match(r"<TR:([+-]?\d+)bp>", a).group(1)) for a in f[4].split(",")]
            info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
            gi = f[8].split(":").index("GT")
            calls: dict[str, Optional[tuple[int, int]]] = {}
            for s, cell in zip(samples, f[9:]):
                if s not in keep:
                    continue
                al = re.split(r"[/|]", cell.split(":")[gi])
                if len(al) != 2 or not all(a.isdigit() for a in al):
                    calls[s] = None
                    continue
                calls[s] = (bp[int(al[0])], bp[int(al[1])])
            loci[f[2]] = (int(info.get("LOCUS_RECORDS", "1")), calls)
    return loci, samples


def trgt_changes(vcf: str, region: str) -> dict[str, tuple[int, int]]:
    cmd = ["bcftools", "query", "-r", region, "-f", "%INFO/TRID\t%REF\t%ALT[\t%GT]\n", vcf]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    assert proc.stdout is not None
    out = {}
    for line in proc.stdout:
        trid, ref, alt, gt = line.rstrip("\n").split("\t")
        al = re.split(r"[/|]", gt)
        if len(al) != 2 or not all(a.isdigit() for a in al):
            continue
        seqs = [ref] + (alt.split(",") if alt not in (".", "") else [])
        out[trid] = (len(seqs[int(al[0])]) - len(ref), len(seqs[int(al[1])]) - len(ref))
    proc.stdout.close()
    if proc.wait() != 0:
        raise SystemExit(f"bcftools query failed on {vcf}")
    return out


def size_class(bp: int) -> str:
    a = abs(bp)
    if a == 0:
        return "0"
    if a < 20:
        return "1-19"
    if a < 50:
        return "20-49"
    return ">=50"


def run(args: argparse.Namespace) -> dict:
    region = args.region or args.chrom
    if not region:
        raise SystemExit("--chrom or --region is required")
    if ":" in region:
        chrom, span = region.rsplit(":", 1)
        lo, hi = (int(x) for x in span.split("-"))
    else:
        chrom, lo, hi = region, None, None
    catalog = load_catalog(args.catalog_bed, chrom, lo, hi)
    pairs = [ln.split("\t") for ln in Path(args.trgt_tsv).read_text().splitlines() if ln.strip()]
    trgt_vcfs = {s.strip(): p.strip() for s, p in pairs}
    loci, vcf_samples = load_locus_vcf(args.locus_vcf, set(trgt_vcfs))
    samples = [s for s in vcf_samples if s in trgt_vcfs]
    if not samples:
        raise SystemExit("no sample is in both the locus VCF and the TRGT table")

    per_locus: dict[str, list[float]] = defaultdict(lambda: [0, 0, 0, 0, 0.0])
    strata: dict[str, dict[str, list[int]]] = {"n_records": defaultdict(lambda: [0, 0, 0, 0]), "trgt_size": defaultdict(lambda: [0, 0, 0, 0])}
    xs: list[int] = []
    ys: list[int] = []
    n_trgt_only = 0
    for i, s in enumerate(samples):
        tr = trgt_changes(trgt_vcfs[s], region)
        print(f"[concordance] {i + 1}/{len(samples)} {len(tr):,} TRGT loci", file=sys.stderr, flush=True)
        for lid, (d1, d2) in tr.items():
            if lid not in catalog:
                n_trgt_only += 1
                continue
            period = catalog[lid]
            nrec, calls = loci.get(lid, (0, {}))
            integ = calls.get(s, (0, 0))
            if integ is None:
                continue
            a, b = sorted(integ), sorted((d1, d2))
            exact = a == b
            within = abs(a[0] - b[0]) <= period and abs(a[1] - b[1]) <= period
            over = abs(sum(a)) - abs(sum(b)) >= OVERCOUNT_BP
            row = per_locus[lid]
            row[0] += 1
            row[1] += exact
            row[2] += within
            row[3] += over
            row[4] += abs(sum(a) - sum(b))
            for key, label in (("n_records", "3+" if nrec >= 3 else str(nrec)), ("trgt_size", size_class(max(b, key=abs)))):
                st = strata[key][label]
                st[0] += 1
                st[1] += exact
                st[2] += within
                st[3] += over
            xs.append(sum(a))
            ys.append(sum(b))

    out = gzip.open(args.out_loci, "wt") if args.out_loci.endswith(".gz") else open(args.out_loci, "w")
    with out as fh:
        fh.write("\t".join(LOCI_COLUMNS) + "\n")
        for lid, (n, ex, wi, ov, absdiff) in per_locus.items():
            nrec = loci[lid][0] if lid in loci else 0
            fh.write(f"{lid}\t{catalog[lid]}\t{nrec}\t{n}\t{ex}\t{wi}\t{ov}\t{absdiff / n:.3f}\n")

    def frac(t: list[int]) -> dict:
        n = max(t[0], 1)
        return {"pairs": t[0], "exact": t[1] / n, "within_unit": t[2] / n, "overcount": t[3] / n}

    total = [sum(r[k] for r in per_locus.values()) for k in range(4)]
    r = None
    if len(xs) > 2:
        import numpy as np

        x, y = np.array(xs, float), np.array(ys, float)
        if x.std() > 0 and y.std() > 0:
            r = float(np.corrcoef(x, y)[0, 1])
    return {
        "chrom": args.chrom,
        "n_samples": len(samples),
        "n_loci_compared": len(per_locus),
        "overall": frac(total),
        "by_integrated_records": {k: frac(v) for k, v in sorted(strata["n_records"].items())},
        "by_trgt_change": {k: frac(v) for k, v in sorted(strata["trgt_size"].items())},
        "pearson_r_total_change": r,
        "trgt_loci_not_in_catalog": n_trgt_only,
        "overcount_bp": OVERCOUNT_BP,
    }


def main(argv: Optional[list[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--locus-vcf", required=True, help="aggregate_repeat_loci.py --out-vcf output")
    p.add_argument("--catalog-bed", required=True)
    p.add_argument("--trgt-tsv", required=True, help="sample<TAB>TRGT VCF path, one per line")
    p.add_argument("--chrom", default=None)
    p.add_argument("--region", default=None, help="chrom or chrom:start-end; limits both callsets to this window")
    p.add_argument("--out-loci", required=True)
    p.add_argument("--out-summary", required=True)
    args = p.parse_args(argv)
    summary = run(args)
    Path(args.out_summary).write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary["overall"]), file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
