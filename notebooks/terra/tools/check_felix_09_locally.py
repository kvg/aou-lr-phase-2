#!/usr/bin/env python3
"""Run felix_09's cell logic locally against the real tables, standing in for the Terra-only calls.

What is NOT exercised here: the FISS entity call (cell 2), the gsutil module copy (cell 1) and the gsutil cat stream
(cell 7's extract_many, which has its own unit tests with a faked gsutil). Everything else -- the locus reduction, the
lead match and its sanity check, the phenotype-context merge, the region compaction and the checks dry run -- runs on
the real 329-run sweep tables and the five region extracts in the manuscript repo's data/.

The region extracts stand in for top_hits.tsv: they carry the same 161 columns, so lead_rows cannot tell them apart.
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "felix" / "scripts"))
import felix_hit_dossier as fhd  # noqa: E402

DATA = REPO.parent / "aou-lr-phase-2-manuscript" / "data"
ARMS = fhd.ARMS


def main() -> int:
    hits = pd.read_csv(DATA / "felix_sweep_all_hits.tsv", sep="\t")
    cal = pd.read_csv(DATA / "felix_sweep_calibration.tsv", sep="\t")

    # ---- cell 3: screen, then reduce to loci
    loci = fhd.dossier_loci(hits)
    kept = fhd.retained(hits)
    print(f"[3] {len(hits)} reported -> {int((hits.verdict == 'keep').sum())} retained "
          f"-> {len(kept)} with >= {fhd.MIN_CASE_CARRIERS} carriers -> {len(loci)} loci "
          f"across {loci.phenotype.nunique()} phenotypes")
    assert len(loci) and loci.locus_id.is_unique

    # ---- cell 4: lead rows, with the region extracts standing in for top_hits
    top = {}
    for path in sorted(DATA.glob("fig3_locus_*.tsv")):
        stem = path.name.removeprefix("fig3_locus_").removesuffix(".tsv")
        tag, arm = stem.split(".", 1)
        top[(f"felix.prod.{tag}", arm)] = pd.read_csv(path, sep="\t", low_memory=False)
    print(f"[4] {len(top)} stand-in run x arm tables for "
          f"{len({k[0] for k in top})} runs")
    leads = fhd.lead_rows(loci, top)
    found = leads[leads.found.astype(bool)]
    print(f"[4] {len(leads)} locus x arm rows; {len(found)} matched, {len(leads) - len(found)} did not")

    chk = found.merge(loci[["locus_id", "lead_arm", "p"]], on="locus_id")
    chk = chk[chk.arm == chk.lead_arm]
    delta = (pd.to_numeric(chk["P_cct_admixed_c"], errors="coerce") / chk["p"] - 1).abs()
    print(f"[4] lead p agrees with the sweep for {int((delta < 1e-6).sum())}/{len(chk)} loci "
          f"(largest relative difference {delta.max():.2e})")
    assert len(chk) and (delta < 1e-6).all(), "a dossier lead disagrees with the sweep's p"

    # ---- cell 5: phenotype context
    wide = cal.pivot_table(index="run", columns="arm", values=["lambda_gc", "n_cases", "N_eff"], aggfunc="first")
    wide.columns = [f"{a}_{b}" for a, b in wide.columns]
    wide = wide.reset_index()
    ctx_cols = wide[["run"] + [c for c in wide.columns if c.startswith("lambda_gc_")]].copy()
    ctx_cols["n_cases"] = wide[[c for c in wide.columns if c.startswith("n_cases_")]].max(axis=1)
    ctx_cols["N_eff"] = wide[[c for c in wide.columns if c.startswith("N_eff_")]].max(axis=1)
    before = len(loci)
    loci = loci.merge(ctx_cols, on="run", how="left")
    assert len(loci) == before, "the calibration merge changed the row count"
    lead_lambda = loci.apply(lambda r: r.get(f"lambda_gc_{r.lead_arm}"), axis=1)
    print(f"[5] lead-arm lambda: median {lead_lambda.median():.3f}, "
          f"range {lead_lambda.min():.3f}-{lead_lambda.max():.3f}; "
          f"{int((lead_lambda < 0.90).sum())} loci below 0.90; "
          f"{int(loci.n_cases.isna().sum())} loci with no case count")
    assert loci.n_cases.notna().all(), "a locus has no case count"

    # ---- cell 7: compaction (the gsutil stream itself is unit-tested separately)
    tmp = Path(tempfile.mkdtemp(prefix="dossier_regions_"))
    by_id = {r.locus_id: r for _, r in loci.iterrows()}
    made = []
    for (run, arm), raw in top.items():
        for lid, L in by_id.items():
            if L.run != run:
                continue
            src = DATA / f"fig3_locus_{run.removeprefix('felix.prod.')}.{arm}.tsv"
            out = fhd.region_path(tmp, lid, arm)
            n = fhd.compact_region(src, out, L)
            made.append(dict(locus_id=lid, arm=arm, kept_rows=n, kb=round(out.stat().st_size / 1024, 1)))
    m = pd.DataFrame(made)
    covered = m[m.kept_rows > 0]
    print(f"[7] {len(m)} compact files, {len(covered)} non-empty, "
          f"{m.kb.sum() / 1024:.2f} MB total, {m.kb.max():.0f} KB largest")
    print(f"[7] projected for 90 loci x 3 arms at the mean non-empty size: "
          f"{covered.kb.mean() * 90 * 3 / 1024:.0f} MB")
    one = pd.read_csv(fhd.region_path(tmp, "BI_174__chr9_8991058", "main"), sep="\t")
    assert list(one.columns) == [c for c in fhd.REGION_COLS if c in one.columns]
    assert fhd.P_COL in one.columns and len(one) > 0
    print(f"[7] BI_174 main region file: {len(one)} rows, {len(one.columns)} columns, "
          f"{(one[fhd.P_COL] < fhd.REGION_KEEP_BELOW).sum()} markers below p 1e-3")

    # ---- cell 8: the checks dry run over every locus
    context = pd.read_csv(DATA / fhd.CONTEXT_FILE, sep="\t")
    rows = []
    for _, L in loci.iterrows():
        lr = leads[(leads.locus_id == L.locus_id) & (leads.arm == L.lead_arm) & leads.found.astype(bool)]
        lead = lr.iloc[0] if len(lr) else None
        cx = context[context.locus_id == L.locus_id]
        rp = fhd.region_path(tmp, L.locus_id, L.lead_arm)
        region = pd.read_csv(rp, sep="\t") if rp.exists() else None
        flags = fhd.checks(L, lead, cx.iloc[0] if len(cx) else None, region)
        rows.append(dict(locus_id=L.locus_id, phenotype=L.phenotype, p=L.p, verdict=fhd.verdict(flags),
                         **{f["check"]: f["status"] for f in flags}))
    review = pd.DataFrame(rows)
    print("[8] verdicts:", review.verdict.value_counts().to_dict())
    for c in [c for c in review.columns if c not in ("locus_id", "phenotype", "p", "verdict")]:
        print(f"      {c:24s} " + ", ".join(f"{k} {v}" for k, v in review[c].value_counts().items()))
    assert len(review) == len(loci)
    out = Path(tempfile.gettempdir()) / "felix_09_local_review.tsv"
    review.to_csv(out, sep="\t", index=False)
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
