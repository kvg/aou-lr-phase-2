"""Best screened P per run and arm: the per-phenotype readout behind the Figure 3 PheWAS panel.

felix_08's sweep keeps only the hits (P < 5e-8), so a phenotype with no hit has no P value in any sweep table. This module
reads the same harvest cache (``harvest_cache/<run>/<arm>/**/top_hits.tsv``, written by felix_phenotype_report.harvest, no
new download) and reports, for every run and arm, the strongest test that survives the artifact screen:

    spa_fallback      SPA attempted and failed, so P is the uncorrected normal approximation
    no_case_carriers  no case carries the allele
    impossible_p      P is exactly 0 (or missing)
    excess_het        controls more heterozygous than Hardy-Weinberg allows (main and ru_biallelic arms)

Same four exclusions as felix_phenotype_report.hit_report, imported rather than re-implemented, so the PheWAS panel and
the hit table agree.

Limits, carried into the columns. top_hits.tsv holds only the strongest N loci per arm (N = n_rows), ordered by the arm's own
joint P column, which is not always the column screened here (P_cct_admixed_c). The best screened P is therefore the best
among those rows. When every row fails the screen (``censored``) the true best screened P is unknown; ``p_bound`` is the
weakest P in the file and the test cannot be placed more precisely than that.

Terra:  gsutil cp $WORKSPACE_BUCKET/felix/scripts/sweep_best_p.py .   (after 00_sync_repo), then see fig3_terra_export_phewas.txt.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import felix_phenotype_report as fpr  # noqa: E402

COLUMNS = ["run", "phenotype", "arm", "n_rows", "n_screened_in", "best_p", "MarkerID", "CHR", "POS", "case_carriers",
           "censored", "p_bound"]


def best_p_report(bundle: dict, num_ancs: int = 5, p_col: str = "P_cct_admixed_c") -> pd.DataFrame:
    """One row per arm that has a top_hits table with ``p_col``."""
    rows = []
    for arm in fpr.ARMS:
        a = bundle["arms"].get(arm, {})
        th = a.get("top_hits")
        if th is None or p_col not in th.columns or not len(th):
            continue
        p = pd.to_numeric(th[p_col], errors="coerce")
        spa, _ = fpr.spa_fallback(th, num_ancs)
        cc = fpr.case_carriers(th)
        het, _, _ = fpr.excess_het(th, arm)
        ok = ~(spa | (cc == 0) | (p == 0) | p.isna() | het)
        positive = p[p > 0]
        row = dict(arm=arm, phenotype=a.get("phenotype"), n_rows=int(len(th)), n_screened_in=int(ok.sum()),
                   p_bound=float(positive.max()) if len(positive) else np.nan)
        if ok.any():
            i = p[ok].idxmin()
            row.update(best_p=float(p[i]), MarkerID=th.at[i, "MarkerID"], CHR=th.at[i, "CHR"], POS=th.at[i, "POS"],
                       case_carriers=int(cc[i]), censored=False)
        else:
            row.update(best_p=np.nan, MarkerID=None, CHR=None, POS=None, case_carriers=None, censored=True)
        rows.append(row)
    return pd.DataFrame(rows)


def sweep_best_p(store: str = "harvest_cache", num_ancs: int = 5, p_col: str = "P_cct_admixed_c") -> pd.DataFrame:
    """best_p_report over every harvested run in ``store`` (directories holding a .done marker)."""
    out = []
    for d in sorted(Path(store).iterdir()):
        if not (d / ".done").exists():
            continue
        r = best_p_report(fpr.read_harvest(d), num_ancs, p_col)
        if len(r):
            r.insert(0, "run", d.name)
            out.append(r)
    if not out:
        return pd.DataFrame(columns=COLUMNS)
    df = pd.concat(out, ignore_index=True)
    df["phenotype"] = df.groupby("run")["phenotype"].transform(lambda s: s.dropna().iloc[0] if s.notna().any() else None)
    return df[COLUMNS]
