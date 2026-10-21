#!/usr/bin/env python3
"""Expected per-ancestry repeat-dosage exposure across the cohort.

A per-ancestry repeat *length* test only carries information at a locus where
the two haplotypes differ in length AND differ in local ancestry. If both
haplotypes are ancestry k, swapping their lengths leaves DS{k} unchanged,
because DS{k} is a sum over one bin -- so such a locus contributes nothing that
distinguishes a per-ancestry effect from a pooled one.

Counting those person-loci exactly needs local ancestry at every TR locus for
every participant, which means streaming a ~29 GiB annotated VCF per
chromosome. This script estimates it instead, from two cheap inputs:

  * global ancestry proportions per participant (FLARE's global.anc.gz, ~126
    KiB), giving P(the two haplotypes differ in ancestry) = 1 - sum_k p_k^2;
  * the length-heterozygosity rate measured directly from TRGT calls.

The independence assumption in 1 - sum_k p_k^2 treats a person's two haplotypes
as independent draws from their own global proportions. Local ancestry is
correlated along a chromosome and between homologues in a recently admixed
population, so this is an UPPER bound on discordance, and the exposure figures
here should be read as optimistic. The point is the order of magnitude and the
shape of the distribution across participants, not a precise count.

It also reports per-ancestry haplotype totals, because those bound power
independently: DS{k} has no variance among participants carrying none of
ancestry k, so an ancestry present at a few per cent cohort-wide supports a far
smaller test than its locus count suggests.
"""

from __future__ import annotations

import argparse
import gzip
import sys

import numpy as np
import pandas as pd


def read_global_anc(path: str) -> pd.DataFrame:
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        df = pd.read_csv(fh, sep="\t")
    first = df.columns[0]
    df = df.rename(columns={first: "sample"})
    df["sample"] = df["sample"].astype(str)
    pops = [c for c in df.columns if c != "sample"]
    assert pops, "no ancestry columns found"
    return df, pops


def discordance(props: np.ndarray) -> np.ndarray:
    """P(two independent haplotype draws differ in ancestry), per row."""
    tot = props.sum(axis=1, keepdims=True)
    p = np.divide(props, tot, out=np.zeros_like(props), where=tot > 0)
    return 1.0 - (p ** 2).sum(axis=1)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--global-anc", required=True,
                    help="FLARE global.anc.gz (SAMPLE + one column per population)")
    ap.add_argument("--length-het-loci", type=int, default=2431,
                    help="TR loci length-heterozygous in the measured window "
                         "(default 2431, sample 1000234, chr20:10-20Mb)")
    ap.add_argument("--window-loci", type=int, default=14618,
                    help="TR loci in that window (default 14618)")
    ap.add_argument("--catalog-loci", type=int, default=4439813,
                    help="loci in the full catalog (default 4,439,813)")
    ap.add_argument("--min-prop", type=float, default=0.05,
                    help="proportion above which a participant counts as carrying "
                         "an ancestry (default 0.05)")
    ap.add_argument("--out", default=None, help="optional per-participant TSV")
    args = ap.parse_args()

    df, pops = read_global_anc(args.global_anc)
    props = df[pops].to_numpy(dtype=float)
    het_rate = args.length_het_loci / args.window_loci

    df["p_anc_discordant"] = discordance(props)
    df["exposed_loci"] = df.p_anc_discordant * het_rate * args.catalog_loci

    print(f"participants:              {len(df):,}")
    print(f"populations:               {', '.join(pops)}")
    print(f"length-het rate (measured): {het_rate:.4f}"
          f"  ({args.length_het_loci:,}/{args.window_loci:,})")
    print(f"catalog loci:              {args.catalog_loci:,}")

    print("\nP(haplotypes differ in ancestry), per participant:")
    q = df.p_anc_discordant.describe(percentiles=[0.05, 0.25, 0.5, 0.75, 0.95])
    for k in ("mean", "5%", "25%", "50%", "75%", "95%", "max"):
        print(f"  {k:>5}  {q[k]:.4f}")
    for thr in (0.01, 0.05, 0.10):
        n = int((df.p_anc_discordant < thr).sum())
        print(f"  below {thr:>4}: {n:,} participants ({n / len(df):.1%})"
              " -- effectively no per-ancestry information")

    print("\nexpected exposed loci per participant (genome-wide):")
    q = df.exposed_loci.describe(percentiles=[0.05, 0.5, 0.95])
    for k in ("mean", "5%", "50%", "95%", "max"):
        print(f"  {k:>5}  {q[k]:,.0f}")
    print(f"\ncohort total exposed person-loci: {df.exposed_loci.sum():,.0f}")

    print("\nper-ancestry haplotype totals (bounds power independently of loci):")
    tot = props.sum(axis=1, keepdims=True)
    p = np.divide(props, tot, out=np.zeros_like(props), where=tot > 0)
    rows = []
    for j, pop in enumerate(pops):
        rows.append({
            "ancestry": pop,
            "haplotypes": 2.0 * p[:, j].sum(),
            "cohort_frac": p[:, j].mean(),
            f"carriers_gt_{args.min_prop}": int((p[:, j] > args.min_prop).sum()),
        })
    anc = pd.DataFrame(rows).sort_values("haplotypes", ascending=False)
    anc["haplotypes"] = anc.haplotypes.round(0).astype(int)
    anc["cohort_frac"] = anc.cohort_frac.round(4)
    print(anc.to_string(index=False))

    if args.out:
        # Not %g: exposed_loci runs to six figures, and %.6g would silently
        # drop everything after the decimal point on exactly those rows.
        df[["sample"] + pops + ["p_anc_discordant", "exposed_loci"]].to_csv(
            args.out, sep="\t", index=False, float_format="%.6f")
        print(f"\nwrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
