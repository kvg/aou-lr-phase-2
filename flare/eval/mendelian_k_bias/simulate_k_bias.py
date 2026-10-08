#!/usr/bin/env python3
"""Is `violations_per_informative_locus` comparable across models with different nanc?

Selection v2 ranks LAI recipes by the pedigree Mendelian violation rate
(`scripts/flare_score_mendelian_lai.py`). For `sel_chr20_flare2_nanc5` vs
`sel_chr20_flare2_nanc6` that is a comparison across a different number of
ancestry labels, so the metric has to be K-free to be usable as a ranking.

A locus counts as a hard violation when the child's two haplotype labels cannot
be drawn one from each parent: neither `c1 in F and c2 in M` nor the reverse.
Compatibility is preserved under *coarsening* — if two labels are merged, every
locus that was compatible stays compatible and some violations disappear — so
the rate is expected to rise with K at fixed LAI accuracy, independent of
whether the painting is any better.

Two regimes, both scored with the production `score_trio_path`:

  uniform   K equiprobable ancestries, misclassification rate eps per
            haplotype-locus. Upper bound on the effect.
  split     The true demography is fixed at 5 ancestries. The K=6 model is the
            same truth with one ancestry split into two labels at a given
            purity, assigned independently per haplotype-locus — a noisy
            refinement, which is what an extra GMM cluster looks like when it
            is not capturing real structure. Isolates the label partition from
            the demography.

Also scores the proposed fix: project every label onto its dominant reference
panel (`flare2_build_model.py` already records `dominant_panel` per cluster)
and score on the projection, so both candidates are compared at one resolution.

Usage::

    python3 flare/eval/mendelian_k_bias/simulate_k_bias.py --out-dir flare/eval/mendelian_k_bias
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "scripts"))
from flare_score_mendelian_lai import score_trio_path  # noqa: E402

# chr20: ~100 cM, and the selection rows use gen=10, so ~10 ancestry switches
# per founder haplotype and ~1 meiotic crossover per transmitted haplotype.
CHROM_MORGANS = 1.0
T_GEN = 10


def founder_mosaic(rng: np.random.Generator, n_hap: int, n_loci: int, pi: np.ndarray) -> np.ndarray:
    """(n_hap, n_loci) ancestry labels: Markov admixture mosaic, rate T per Morgan."""
    switch = rng.random((n_hap, n_loci)) < (T_GEN * CHROM_MORGANS / n_loci)
    switch[:, 0] = True
    draws = rng.choice(len(pi), size=(n_hap, n_loci), p=pi)
    idx = np.maximum.accumulate(np.where(switch, np.arange(n_loci)[None, :], 0), axis=1)
    return np.take_along_axis(draws, idx, axis=1)


def transmit(rng: np.random.Generator, h1: np.ndarray, h2: np.ndarray, n_loci: int) -> np.ndarray:
    """Meiosis: mosaic of a parent's two haplotypes, ~1 crossover per Morgan."""
    xo = rng.random((h1.shape[0], n_loci)) < (CHROM_MORGANS / n_loci)
    which = np.cumsum(xo, axis=1) % 2
    start = rng.integers(0, 2, size=(h1.shape[0], 1))
    which = (which + start) % 2
    return np.where(which == 0, h1, h2)


def misclassify(rng: np.random.Generator, truth: np.ndarray, eps: float, k: int) -> np.ndarray:
    """With probability eps, replace the label with a different one (uniform)."""
    if eps <= 0:
        return truth
    hit = rng.random(truth.shape) < eps
    shift = rng.integers(1, k, size=truth.shape)  # never 0, so always a real error
    return np.where(hit, (truth + shift) % k, truth)


def score(c1, c2, f1, f2, m1, m2) -> tuple[int, int]:
    """Pooled (informative, hard violations) over trios, via the production scorer."""
    tot_inf = tot_hard = 0
    for i in range(c1.shape[0]):
        loci = list(zip(c1[i].tolist(), c2[i].tolist(), f1[i].tolist(),
                        f2[i].tolist(), m1[i].tolist(), m2[i].tolist()))
        s = score_trio_path(loci)
        tot_inf += s["n_informative"]
        tot_hard += s["n_hard_violations"]
    return tot_inf, tot_hard


def simulate_uniform(rng, n_trios, n_loci, k, eps):
    pi = np.full(k, 1.0 / k)
    f1 = founder_mosaic(rng, n_trios, n_loci, pi)
    f2 = founder_mosaic(rng, n_trios, n_loci, pi)
    m1 = founder_mosaic(rng, n_trios, n_loci, pi)
    m2 = founder_mosaic(rng, n_trios, n_loci, pi)
    c1 = transmit(rng, f1, f2, n_loci)
    c2 = transmit(rng, m1, m2, n_loci)
    obs = [misclassify(rng, x, eps, k) for x in (c1, c2, f1, f2, m1, m2)]
    return score(*obs)


def simulate_split(rng, n_trios, n_loci, k, eps, purity, pi5):
    """True demography is always 5 ancestries; K=6 splits ancestry 4 in two.

    The split is noisy: a haplotype of true ancestry 4 is labelled 4 with
    probability `purity` and 5 otherwise, drawn independently per locus. K=5
    scores the same truth with no split.
    """
    f1 = founder_mosaic(rng, n_trios, n_loci, pi5)
    f2 = founder_mosaic(rng, n_trios, n_loci, pi5)
    m1 = founder_mosaic(rng, n_trios, n_loci, pi5)
    m2 = founder_mosaic(rng, n_trios, n_loci, pi5)
    c1 = transmit(rng, f1, f2, n_loci)
    c2 = transmit(rng, m1, m2, n_loci)
    truth = [c1, c2, f1, f2, m1, m2]
    if k == 6:
        truth = [np.where((x == 4) & (rng.random(x.shape) >= purity), 5, x) for x in truth]
    obs = [misclassify(rng, x, eps, k) for x in truth]
    inf_raw, hard_raw = score(*obs)
    # Projection fix: collapse the split label back to its dominant panel.
    proj = [np.where(x == 5, 4, x) for x in obs] if k == 6 else obs
    inf_p, hard_p = score(*proj)
    return (inf_raw, hard_raw), (inf_p, hard_p)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--n-trios", type=int, default=101)   # the Phase 2 pedigree
    ap.add_argument("--n-loci", type=int, default=3000)
    ap.add_argument("--reps", type=int, default=3)
    ap.add_argument("--seed", type=int, default=20261008)
    a = ap.parse_args()
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    eps_grid = [0.0, 0.005, 0.01, 0.02, 0.05, 0.10]
    rows = []

    for k in (3, 4, 5, 6, 7, 8):
        for eps in eps_grid:
            for rep in range(a.reps):
                rng = np.random.default_rng(a.seed + 1000 * k + int(eps * 1e4) + rep)
                inf, hard = simulate_uniform(rng, a.n_trios, a.n_loci, k, eps)
                rows.append({"regime": "uniform", "K": k, "eps": eps, "rep": rep,
                             "scoring": "raw", "n_informative": inf, "n_violations": hard,
                             "rate": hard / inf if inf else np.nan})

    # Admixed-cohort-like proportions; ancestry 4 is the one that gets split.
    pi5 = np.array([0.42, 0.10, 0.06, 0.04, 0.38])
    for purity in (0.85, 0.95):
        for k in (5, 6):
            for eps in eps_grid:
                for rep in range(a.reps):
                    rng = np.random.default_rng(a.seed + 7000 * k + int(eps * 1e4)
                                                + rep + int(purity * 100))
                    (inf, hard), (infp, hardp) = simulate_split(
                        rng, a.n_trios, a.n_loci, k, eps, purity, pi5)
                    rows.append({"regime": f"split_purity{purity}", "K": k, "eps": eps,
                                 "rep": rep, "scoring": "raw",
                                 "n_informative": inf, "n_violations": hard,
                                 "rate": hard / inf if inf else np.nan})
                    rows.append({"regime": f"split_purity{purity}", "K": k, "eps": eps,
                                 "rep": rep, "scoring": "projected",
                                 "n_informative": infp, "n_violations": hardp,
                                 "rate": hardp / infp if infp else np.nan})

    df = pd.DataFrame(rows)
    df.to_csv(out / "k_bias_raw.tsv", sep="\t", index=False)
    agg = (df.groupby(["regime", "scoring", "K", "eps"])["rate"]
             .agg(rate_mean="mean", rate_sd="std").reset_index())
    agg.to_csv(out / "k_bias_summary.tsv", sep="\t", index=False)

    u = agg[(agg.regime == "uniform") & (agg.scoring == "raw")]
    piv = u.pivot(index="eps", columns="K", values="rate_mean")

    def pick(regime, k, scoring, eps=0.02):
        m = agg[(agg.regime == regime) & (agg.K == k) & (agg.scoring == scoring)
                & (np.isclose(agg.eps, eps))]
        return float(m.rate_mean.iloc[0])

    summary = {
        "n_trios": a.n_trios, "n_loci": a.n_loci, "reps": a.reps,
        "uniform_rate_by_K": {f"eps={e}": {f"K={k}": round(float(piv.loc[e, k]), 6)
                                           for k in piv.columns} for e in piv.index},
    }
    for purity in (0.85, 0.95):
        reg = f"split_purity{purity}"
        summary[reg] = {
            "at_eps": 0.02,
            "K5_raw": round(pick(reg, 5, "raw"), 6),
            "K6_raw": round(pick(reg, 6, "raw"), 6),
            "K6_projected": round(pick(reg, 6, "projected"), 6),
        }
    (out / "k_bias_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
