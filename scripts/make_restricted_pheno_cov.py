#!/usr/bin/env python3
"""Write a sex-restricted copy of pheno_cov.tsv for one phenotype.

A phecode PhecodeX labels "Both" but whose cases are all one sex leaves `sex` in
the covariate set, where it separates the cases perfectly. SAIGE step 1 then
diverges and still writes a .rda (observed on CA_142.2: converged FALSE,
theta 2.0e13, intercept -1.7e15, every fitted value 2.22e-16), and every step-2
shard for that phenotype aborts with a boost gamma_q_inv overflow.

The fix is to make the excluded sex MISSING for that phenotype, so SAIGE drops
those participants and sex becomes constant among the rest. Only that phenotype's
column is touched, and the file is written per phenotype, so the shared
pheno_cov.tsv keeps its content hash and the other 947 nulls stay cache-hit.

Do NOT subset analysis_samples to achieve this instead: that file is an input to
PackFelixla and MakeGRM, so changing it re-packs every chromosome.

Usage:
  make_restricted_pheno_cov.py --pheno-cov pheno_cov.tsv --covariates covariates.v8.csv.gz \
      --phenotype CA_142_2 --restrict-to Female --out pheno_cov.CA_142_2.tsv

The phenotype name is the one in selected_phenotypes.txt (underscore dialect).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

SEX_FEMALE = {"F", "f", "female", "Female", "2", "XX", "xx"}
SEX_MALE = {"M", "m", "male", "Male", "1", "XY", "xy"}


def encode_sex(s: pd.Series) -> pd.Series:
    a = s.astype(str).str.strip()
    out = pd.Series(pd.NA, index=s.index, dtype="object")
    out[a.isin(SEX_FEMALE)] = "Female"
    out[a.isin(SEX_MALE)] = "Male"
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pheno-cov", type=Path, required=True)
    ap.add_argument("--covariates", type=Path, required=True,
                    help="table carrying sex_at_birth / inferred_sex")
    ap.add_argument("--phenotype", required=True, help="column to restrict")
    ap.add_argument("--restrict-to", required=True, choices=("Female", "Male"))
    ap.add_argument("--sex-col", default=None,
                    help="sex column IN pheno_cov to restrict on (default: 'sex' if present). "
                         "Pass an empty string to force use of the covariate table instead.")
    ap.add_argument("--min-sex-purity", type=float, default=0.99,
                    help="required purity when labelling a pheno_cov sex level [0.99]")
    ap.add_argument("--id-col", default=None, help="id column in pheno_cov (auto-detected)")
    ap.add_argument("--cov-id-col", default="research_id")
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)

    pc = pd.read_csv(a.pheno_cov, sep="\t", dtype=str, low_memory=False)
    if a.phenotype not in pc.columns:
        raise SystemExit(f"{a.phenotype!r} is not a column in {a.pheno_cov}; "
                         f"first columns are {list(pc.columns[:6])}")
    id_col = a.id_col or next(
        (c for c in ("ID", "IID", "person_id", "research_id", "sample_id", "s") if c in pc.columns),
        pc.columns[0])

    cov = pd.read_csv(a.covariates, dtype=str, low_memory=False)
    sex_src = next((c for c in ("sex_at_birth", "inferred_sex") if c in cov.columns), None)
    if sex_src is None:
        raise SystemExit(f"{a.covariates} has neither sex_at_birth nor inferred_sex")
    cov["_sex"] = encode_sex(cov[sex_src])
    fallback = cov["_sex"].isna()
    if "inferred_sex" in cov.columns and fallback.any():
        cov.loc[fallback, "_sex"] = encode_sex(cov.loc[fallback, "inferred_sex"])

    sex_by_id = dict(zip(cov[a.cov_id_col].astype(str).str.strip(), cov["_sex"]))
    ids = pc[id_col].astype(str).str.strip()
    cov_sex = ids.map(sex_by_id)

    # Prefer pheno_cov's OWN sex column when it has one. That column is what
    # goes into the design matrix and what separates the cases, so it is the
    # one that has to become constant. On the production file it is coded 0/1
    # and does NOT equal sex_at_birth: 222 participants whose sex_at_birth is
    # "PMI: Skip" / "I prefer not to answer" / "Intersex" / "None of these"
    # still carry a 0/1 here (presumably from inferred sex). Restricting on
    # sex_at_birth would set those to NA even though the model calls them
    # female, costing 11 cases on CA_142_2 alone.
    #
    # The covariate table is used only to LABEL the levels, never to decide who
    # is dropped: each level must map to one label essentially purely, or we
    # refuse rather than guess.
    sex_col = (a.sex_col or None) if a.sex_col is not None else (
        "sex" if "sex" in pc.columns else None)
    if sex_col is not None and sex_col not in pc.columns:
        raise SystemExit(f"--sex-col {sex_col!r} is not a column in {a.pheno_cov}")

    if sex_col is None:
        sex, sex_desc = cov_sex, f"{a.covariates.name}:{sex_src}"
    else:
        raw = pc[sex_col].astype(str).str.strip()
        mapping, report = {}, []
        for level, grp in raw.groupby(raw):
            counts = cov_sex[grp.index].value_counts()
            if counts.empty:
                raise SystemExit(f"sex level {level!r} has no labelled participants")
            label, n = counts.index[0], int(counts.iloc[0])
            purity = n / int(counts.sum())
            report.append(f"{level!r}->{label} ({purity:.4%} of {int(counts.sum())} labelled)")
            if purity < a.min_sex_purity:
                raise SystemExit(
                    f"sex level {level!r} in {a.pheno_cov} is only {purity:.2%} {label} "
                    f"in {sex_src}; refusing to guess. Pass --sex-col '' to restrict on "
                    f"the covariate table instead, or raise --min-sex-purity deliberately.")
            mapping[level] = label
        if len(set(mapping.values())) < 2:
            raise SystemExit(f"every level of {sex_col!r} maps to {set(mapping.values())}; "
                             f"the two sexes are not distinguishable")
        sex = raw.map(mapping)
        sex_desc = f"{a.pheno_cov.name}:{sex_col} [{'; '.join(report)}]"

    # The file is read as str so nothing is reformatted on write, but case
    # status must be counted numerically: the production file stores "1.0"/
    # "0.0", not "1"/"0", and a string comparison silently counts zero cases
    # and then trips the no-cases-left guard.
    def case_control(series):
        num = pd.to_numeric(series, errors="coerce")
        return int((num == 1).sum()), int((num == 0).sum())

    before_cases, before_ctrl = case_control(pc[a.phenotype])
    drop = sex.ne(a.restrict_to)          # NA sex is also dropped: it cannot be verified
    pc.loc[drop, a.phenotype] = "NA"
    after_cases, after_ctrl = case_control(pc[a.phenotype])

    pc.to_csv(a.out, sep="\t", index=False, na_rep="NA")

    print(f"{a.phenotype}: restricted to {a.restrict_to}")
    print(f"  id column      : {id_col}")
    print(f"  sex from       : {sex_desc}")
    print(f"  cases          : {before_cases} -> {after_cases}")
    print(f"  controls       : {before_ctrl} -> {after_ctrl}")
    print(f"  set to NA      : {int(drop.sum())} participants")
    print(f"  other columns  : untouched ({pc.shape[1] - 1} besides {id_col})")
    print(f"  wrote {a.out}")
    if after_cases != before_cases:
        print(f"  WARNING: {before_cases - after_cases} case(s) were in the excluded sex "
              f"and have been dropped. For a phecode that is genuinely both-sex "
              f"(e.g. male breast cancer) decide deliberately whether that is right.",
              file=sys.stderr)
    if after_cases == 0:
        raise SystemExit("no cases left after restriction — wrong --restrict-to?")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
