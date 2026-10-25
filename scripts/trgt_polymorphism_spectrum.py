#!/usr/bin/env python3
"""Per-locus repeat-length polymorphism from a sample of production TRGT calls.

Genotyping all 4,439,813 TRExplorer loci costs 209 of the 268 minutes of a
whole-genome TrgtHaplotypePhase task -- 78% of the wall clock and, at 8 cores,
27.9 core-hours per sample. Across 12,261 participants that is ~342,000
core-hours, and it is the only part of the run whose cost is not fixed: the
locus count scales it linearly. Nothing else in the design moves the bill.

So the question is how many loci carry testable length variation. A locus whose
repeat length never varies cannot contribute to a dosage test no matter how
many participants are genotyped, and one that is invariant across a few hundred
unrelated people is very unlikely to vary usefully across the cohort.

This reads a RANDOM SAMPLE of the production TRGT VCFs -- which already exist
for every participant, so this costs an egress read rather than a re-genotyping
run -- and accumulates, per locus, enough to describe the length spectrum:

    n_called, n_alleles, n_het, al_min, al_max, al_sum, al_sumsq

Those are sufficient statistics for the mean, SD, range and length-heterozygosity
rate, and they fit in fixed-size arrays, so memory does not grow with the number
of samples read. Holding the full genotype matrix would not fit: 4.44M loci x
200 samples is 888M calls.

WHAT THIS DOES NOT MEASURE. Length-heterozygosity rate in a sample of
participants is an estimate of the rate in the cohort, not of statistical power.
A locus can be polymorphic and still untestable -- if the variation is confined
to one ancestry, if it sits in a segmental duplication, or if the calls are
driven by coverage rather than biology. Treat the retained-locus counts as an
upper bound on a useful catalog, and note that the sampled participants inherit
whatever ancestry composition the sample has: a 200-participant draw will
underrepresent amr and sas, whose haplotype counts already bound the per-ancestry
test (1,641 and 3,200 respectively).

Locus ORDER is assumed identical across VCFs -- true when every file was called
against one catalog, which is how the production batches were run -- and that
assumption is CHECKED, not trusted: the TRID column is hashed per file and the
digests must agree. A mismatch means the files disagree about the catalog, and
every per-locus number here would be silently misaligned.
"""

from __future__ import annotations

import argparse
import hashlib
import random
import subprocess
import sys
from typing import Dict

import numpy as np
import pandas as pd

MISSING_MAX = np.iinfo(np.int32).max
FIELDS = ("n_called", "n_alleles", "n_het", "al_min", "al_max", "al_sum", "al_sumsq")


def new_accumulator(n_loci: int) -> Dict[str, np.ndarray]:
    """Fixed-size per-locus accumulators; memory is independent of sample count."""
    return {
        "n_called": np.zeros(n_loci, np.int32),
        "n_alleles": np.zeros(n_loci, np.int32),
        "n_het": np.zeros(n_loci, np.int32),
        "al_min": np.full(n_loci, MISSING_MAX, np.int32),
        "al_max": np.full(n_loci, -1, np.int32),
        "al_sum": np.zeros(n_loci, np.int64),
        "al_sumsq": np.zeros(n_loci, np.int64),
    }


def merge(into: Dict[str, np.ndarray], other: Dict[str, np.ndarray]) -> Dict[str, np.ndarray]:
    for key in ("n_called", "n_alleles", "n_het", "al_sum", "al_sumsq"):
        into[key] += other[key]
    np.minimum(into["al_min"], other["al_min"], out=into["al_min"])
    np.maximum(into["al_max"], other["al_max"], out=into["al_max"])
    return into


def accumulate_frame(acc: Dict[str, np.ndarray], offset: int, al: pd.Series) -> int:
    """Fold one chunk of a single sample's AL column into the accumulators.

    ``al`` holds TRGT's FORMAT/AL as written: "123,130" for a diploid call, a
    single value where the locus is hemizygous, "." or an empty string where
    there is no call. Both alleles count towards the length statistics; only a
    diploid call with two DIFFERENT lengths counts as length-heterozygous, which
    is the quantity that makes a locus orientable and a dosage test informative.
    """
    n = len(al)
    parts = al.str.split(",", n=1, expand=True)
    if parts.shape[1] < 2:                       # whole chunk hemizygous
        parts[1] = None
    a0 = pd.to_numeric(parts[0], errors="coerce").to_numpy(dtype=float)
    a1 = pd.to_numeric(parts[1], errors="coerce").to_numpy(dtype=float)
    ok0, ok1 = ~np.isnan(a0), ~np.isnan(a1)

    i0 = np.where(ok0, a0, 0).astype(np.int64)
    i1 = np.where(ok1, a1, 0).astype(np.int64)
    sl = slice(offset, offset + n)

    acc["n_called"][sl] += (ok0 | ok1).astype(np.int32)
    acc["n_alleles"][sl] += ok0.astype(np.int32) + ok1.astype(np.int32)
    acc["n_het"][sl] += (ok0 & ok1 & (i0 != i1)).astype(np.int32)
    acc["al_sum"][sl] += i0 + i1
    acc["al_sumsq"][sl] += i0 * i0 + i1 * i1

    lo = np.where(ok0, i0, MISSING_MAX)
    np.minimum(lo, np.where(ok1, i1, MISSING_MAX), out=lo)
    hi = np.where(ok0, i0, -1)
    np.maximum(hi, np.where(ok1, i1, -1), out=hi)
    np.minimum(acc["al_min"][sl], lo.astype(np.int32), out=acc["al_min"][sl])
    np.maximum(acc["al_max"][sl], hi.astype(np.int32), out=acc["al_max"][sl])
    return n


QUERY_FORMAT = "%INFO/TRID\t[%AL]\n"


def open_query(vcf: str, bcftools: str, gcs_cat: str):
    """Start a `bcftools query` over ``vcf``, streaming it if it lives on GCS.

    htslib in these images has no GCS backend -- handing it a gs:// path fails
    with "Protocol not supported" -- so a remote object is piped in through
    `gcloud storage cat` instead. bcftools detects BGZF on stdin, so reading
    from "-" is equivalent to opening the file, and nothing is written to disk:
    at ~184 MiB per production VCF, 200 of them would otherwise be 37 GiB.

    Returns (query_process, upstream_process_or_None).
    """
    if vcf.startswith("gs://"):
        cat = subprocess.Popen([*gcs_cat.split(), vcf],
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        query = subprocess.Popen([bcftools, "query", "-f", QUERY_FORMAT, "-"],
                                 stdin=cat.stdout, stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE)
        cat.stdout.close()      # only the child keeps the read end
        return query, cat
    query = subprocess.Popen([bcftools, "query", "-f", QUERY_FORMAT, vcf],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return query, None


def read_one(vcf: str, acc: Dict[str, np.ndarray], n_loci: int,
             bcftools: str = "bcftools", chunk_rows: int = 500_000,
             gcs_cat: str = "gcloud storage cat") -> str:
    """Stream one VCF's TRID/AL columns into ``acc``; return the TRID digest."""
    proc, upstream = open_query(vcf, bcftools, gcs_cat)
    hasher = hashlib.blake2b(digest_size=16)
    offset = 0
    try:
        frames = pd.read_csv(proc.stdout, sep="\t", header=None, names=["trid", "al"],
                             dtype=str, na_filter=False, chunksize=chunk_rows,
                             quoting=3)
        for frame in frames:
            if offset + len(frame) > n_loci:
                raise SystemExit(
                    f"{vcf}: more than --n-loci {n_loci} records; the file was called "
                    f"against a different catalog, so per-locus stats would misalign")
            hasher.update("\n".join(frame["trid"]).encode())
            offset += accumulate_frame(acc, offset, frame["al"])
    finally:
        err = proc.stderr.read().decode(errors="replace") if proc.stderr else ""
        proc.stdout.close()
        rc = proc.wait()
        up_rc, up_err = 0, ""
        if upstream is not None:
            up_err = upstream.stderr.read().decode(errors="replace")
            up_rc = upstream.wait()
    if rc != 0:
        raise SystemExit(f"{vcf}: bcftools query exited {rc}\n{err.strip()}")
    # A failed `gcloud storage cat` can still let bcftools exit 0 on a short
    # stream, so the record-count guard below is not the only thing standing
    # between a truncated download and a quietly wrong answer.
    if up_rc != 0:
        raise SystemExit(f"{vcf}: streaming the object exited {up_rc}\n{up_err.strip()}")
    if offset != n_loci:
        raise SystemExit(
            f"{vcf}: {offset} records but --n-loci is {n_loci}; a short file means a "
            f"different catalog or a truncated object, not a tolerable difference")
    return hasher.hexdigest()


def read_block(paths, n_loci: int, bcftools: str, gcs_cat: str):
    """Fold a block of VCFs into one accumulator. Runs in a worker process.

    Blocks rather than one task per file: the accumulator is ~160 MB at the
    production catalog size, and returning one per FILE would ship that back
    200 times. One per worker ships it once.
    """
    acc = new_accumulator(n_loci)
    digests = set()
    for vcf in paths:
        digests.add(read_one(vcf, acc, n_loci, bcftools=bcftools, gcs_cat=gcs_cat))
        print(f"  read {vcf}", file=sys.stderr, flush=True)
    return acc, digests


def read_all(paths, n_loci: int, bcftools: str, gcs_cat: str, threads: int):
    """Accumulate every VCF, in parallel when asked, and verify catalog identity."""
    if threads <= 1 or len(paths) == 1:
        blocks = [read_block(paths, n_loci, bcftools, gcs_cat)]
    else:
        from concurrent.futures import ProcessPoolExecutor
        from functools import partial
        chunks = [list(c) for c in np.array_split(np.array(paths, dtype=object),
                                                  min(threads, len(paths)))]
        worker = partial(read_block, n_loci=n_loci, bcftools=bcftools, gcs_cat=gcs_cat)
        try:
            with ProcessPoolExecutor(max_workers=len(chunks)) as pool:
                blocks = list(pool.map(worker, chunks))
        except (PermissionError, OSError) as exc:
            # Some sandboxes refuse to fork. Falling back is slower by roughly
            # `threads`x but still correct, which beats losing the run.
            print(f"worker pool unavailable ({exc}); falling back to one process",
                  file=sys.stderr, flush=True)
            blocks = [read_block(paths, n_loci, bcftools, gcs_cat)]
    total = new_accumulator(n_loci)
    digests: set = set()
    for acc, block_digests in blocks:
        merge(total, acc)
        digests |= block_digests
    if len(digests) > 1:
        raise SystemExit(
            "the sampled VCFs do not share a TRID column: they disagree about the "
            "catalog or its order, so accumulating them by locus ordinal would "
            f"silently mix different loci together ({len(digests)} distinct layouts)")
    return total


def summarise(acc: Dict[str, np.ndarray], thresholds, cohort_n: int) -> pd.DataFrame:
    called = acc["n_called"]
    seen = called > 0
    het_rate = np.divide(acc["n_het"], called, out=np.zeros(len(called)), where=seen)
    spread = np.where(seen, acc["al_max"].astype(np.int64) - acc["al_min"], -1)
    rows = [
        ("loci_total", len(called), ""),
        ("loci_with_any_call", int(seen.sum()), f"{seen.mean():.4f} of total"),
        ("loci_monomorphic", int((seen & (spread == 0)).sum()),
         "no length variation among sampled participants"),
        ("loci_polymorphic", int((seen & (spread > 0)).sum()), ""),
    ]
    for thr in thresholds:
        keep = int((het_rate >= thr).sum())
        rows.append((f"loci_het_rate_ge_{thr:g}", keep,
                     f"{keep / len(called):.4f} of catalog; "
                     f"~{int(round(thr * cohort_n)):,} het participants at the threshold"))
    return pd.DataFrame(rows, columns=["metric", "value", "note"])


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcf-list", required=True,
                    help="file of production TRGT VCF paths, one per line")
    ap.add_argument("--max-samples", type=int, default=200,
                    help="how many to draw at random (default 200)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--n-loci", type=int, default=4439813,
                    help="records expected per VCF; a mismatch is a hard error")
    ap.add_argument("--cohort-n", type=int, default=12261,
                    help="cohort size, for translating rates into participant counts")
    ap.add_argument("--thresholds", default="0.001,0.005,0.01,0.02,0.05,0.10",
                    help="length-heterozygosity rates to report retained loci for")
    ap.add_argument("--threads", type=int, default=4,
                    help="worker processes; each streams its own block of VCFs")
    ap.add_argument("--bcftools", default="bcftools")
    ap.add_argument("--gcs-cat", default="gcloud storage cat",
                    help="command used to stream gs:// objects into bcftools")
    ap.add_argument("--out-summary", default=None, help="summary TSV (default stdout)")
    ap.add_argument("--out-loci", default=None,
                    help="optional per-locus TSV -- one row per locus, so large")
    args = ap.parse_args()

    with open(args.vcf_list) as handle:
        paths = [line.strip() for line in handle if line.strip()]
    if not paths:
        raise SystemExit(f"{args.vcf_list}: no paths")
    if len(paths) > args.max_samples:
        paths = random.Random(args.seed).sample(paths, args.max_samples)

    print(f"reading {len(paths)} VCFs with {args.threads} worker(s)",
          file=sys.stderr, flush=True)
    acc = read_all(paths, args.n_loci, args.bcftools, args.gcs_cat, args.threads)

    thresholds = [float(x) for x in args.thresholds.split(",") if x.strip()]
    summary = summarise(acc, thresholds, args.cohort_n)
    summary.insert(0, "n_samples_read", len(paths))
    out = args.out_summary or sys.stdout
    summary.to_csv(out, sep="\t", index=False)

    if args.out_loci:
        called = acc["n_called"]
        seen = called > 0
        mean = np.divide(acc["al_sum"], acc["n_alleles"], out=np.zeros(len(called)),
                         where=acc["n_alleles"] > 0)
        var = np.divide(acc["al_sumsq"], acc["n_alleles"], out=np.zeros(len(called)),
                        where=acc["n_alleles"] > 0) - mean ** 2
        pd.DataFrame({
            "locus_ordinal": np.arange(len(called)),
            "n_called": called,
            "n_het": acc["n_het"],
            "het_rate": np.divide(acc["n_het"], called, out=np.zeros(len(called)), where=seen),
            "al_min": np.where(seen, acc["al_min"], -1),
            "al_max": acc["al_max"],
            "al_mean": np.round(mean, 3),
            "al_sd": np.round(np.sqrt(np.maximum(var, 0)), 3),
        }).to_csv(args.out_loci, sep="\t", index=False)
    return 0


if __name__ == "__main__":
    sys.exit(main())
