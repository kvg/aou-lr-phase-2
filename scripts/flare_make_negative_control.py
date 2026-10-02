#!/usr/bin/env python3
"""Negative-control LAI: permute ancestry tracks between samples of one population.

Renames the sample columns of a recipe ``anc.vcf.gz`` with a within-population
derangement (``bcftools reheader -s``), so sample X carries sample Y's AN1/AN2
tracks. Tract lengths, switch rates and population-level ancestry composition
are unchanged; the link between a person's alleles (and relatives) and their
ancestry calls is broken.

A recipe-selection metric that scores this control as well as the real recipe
cannot tell good LAI from bad, so the control must fail every gate.

Scorers read genotypes from ``gt_vcf`` and ancestry from ``anc.vcf``, so the
permuted GT column in the control is never used.

Example::

    python3 scripts/flare_make_negative_control.py \\
      --anc-vcf chr20_flat_props_pin_t.anc.vcf.gz \\
      --covariates covariates.source_rebuilt.csv.gz \\
      --out negctl_chr20_flat_props_pin_t.anc.vcf.gz --seed 7
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import random
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Optional, Sequence


def load_populations(path: Path, id_column: str, pop_column: str) -> dict[str, str]:
    opener = gzip.open if str(path).endswith(".gz") else open
    out: dict[str, str] = {}
    with opener(path, "rt", newline="") as fh:
        sample = fh.read(4096)
        fh.seek(0)
        delim = "\t" if sample.count("\t") > sample.count(",") else ","
        for row in csv.DictReader(fh, delimiter=delim):
            sid = (row.get(id_column) or "").strip()
            pop = (row.get(pop_column) or "").strip().upper()
            if sid:
                out[sid] = pop or "NA"
    return out


def derangement(items: Sequence[str], rng: random.Random) -> list[str]:
    """Random permutation with no fixed points (identity when len < 2)."""
    items = list(items)
    if len(items) < 2:
        return items
    # Sattolo's algorithm: uniform over cyclic permutations, never a fixed point.
    perm = items[:]
    for i in range(len(perm) - 1, 0, -1):
        j = rng.randrange(i)
        perm[i], perm[j] = perm[j], perm[i]
    return perm


def within_population_mapping(
    samples: Sequence[str], pops: dict[str, str], seed: int
) -> dict[str, str]:
    """old column name -> new column name; each sample keeps its population."""
    rng = random.Random(seed)
    groups: dict[str, list[str]] = {}
    for s in samples:
        groups.setdefault(pops.get(s, "NA"), []).append(s)
    mapping: dict[str, str] = {}
    for pop in sorted(groups):
        members = groups[pop]
        for old, new in zip(members, derangement(members, rng)):
            mapping[old] = new
    return mapping


def vcf_samples(vcf: str) -> list[str]:
    out = subprocess.run(["bcftools", "query", "-l", vcf], check=True, capture_output=True, text=True)
    return [ln for ln in out.stdout.splitlines() if ln.strip()]


def main(argv: Optional[list[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--anc-vcf", required=True)
    p.add_argument("--covariates", required=True, type=Path)
    p.add_argument("--id-column", default="research_id")
    p.add_argument("--pop-column", default="ancestry_pred")
    p.add_argument("--seed", type=int, default=7)
    p.add_argument("--out", required=True, type=Path)
    args = p.parse_args(argv)

    samples = vcf_samples(args.anc_vcf)
    pops = load_populations(args.covariates, args.id_column, args.pop_column)
    mapping = within_population_mapping(samples, pops, args.seed)
    n_moved = sum(1 for k, v in mapping.items() if k != v)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", suffix=".rename.txt", delete=False) as fh:
        for s in samples:
            fh.write(f"{s} {mapping[s]}\n")
        rename = fh.name
    subprocess.run(
        ["bcftools", "reheader", "-s", rename, "-o", str(args.out), args.anc_vcf], check=True
    )
    subprocess.run(["bcftools", "index", "-f", "-t", str(args.out)], check=True)
    meta = {
        "source_anc_vcf": args.anc_vcf,
        "out": str(args.out),
        "seed": args.seed,
        "pop_column": args.pop_column,
        "n_samples": len(samples),
        "n_permuted": n_moved,
        "n_unlabeled": sum(1 for s in samples if s not in pops),
    }
    Path(str(args.out) + ".negctl.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(json.dumps(meta), file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
