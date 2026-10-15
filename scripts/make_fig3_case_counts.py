#!/usr/bin/env python3
"""Write the case-count table behind Figure 3 panel A from a select_phenotypes.py run.

Reads the selector's own cohort (``analysis_samples.txt`` and the ``sex`` column of ``pheno_cov.tsv``), so the
counts describe exactly the participants the FELIX scan was fitted on, and applies the selector's own phenotype
encoding and ``--sex-specific restrict`` rule (a sex-specific phecode is counted only in the matching sex). Every
candidate phecode is written, not only the selected ones, because the cutoff curve needs counts below the floor.

Outputs
  * ``--out-counts``: one row per phecode: n_cases, n_eligible, cases_<group> for the seven participant
    groups of Figs. 1-2, and hap_<ancestry>, the expected case haplotypes per local ancestry (2 x the sum of a
    case's chr1 FLARE global ancestry proportion).
  * ``--out-summary``: cohort size, sex split, participants per group and expected haplotypes per local ancestry.

``--legacy-sex-at-birth`` instead keeps only participants with a defined sex_at_birth and restricts on that field.
That is the cohort the first version of these tables used (n = 9,375); running it must reproduce the old files
exactly, which is the check that this script and the old producer agree.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import select_phenotypes as sp  # noqa: E402

GROUPS = ["AFR", "AMR", "EAS", "EUR", "MID", "SAS", "OTH"]
LOCAL = ["AFR", "AMR", "EAS", "EUR", "SAS"]
ANC_COLS = {"eas": "EAS", "amr": "AMR", "eur": "EUR", "afr": "AFR", "sas": "SAS"}  # FLARE panel order eas,amr,eur,afr,sas


def truthy(x) -> bool:
    return str(x).strip().lower() in ("true", "1", "1.0", "t", "yes")


def _num(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return None if v != v else v


def load_covariates(path: Path) -> dict[str, dict]:
    """research_id -> participant group and sex_at_birth, for releasable PacBio participants."""
    out = {}
    with gzip.open(path, "rt", newline="") as fh:
        for r in csv.DictReader(fh):
            if not truthy(r.get("final_releasable_v9", "")) or r.get("technology") != "PacBio":
                continue
            a = (r.get("ancestry_pred_other") or "oth").upper()
            out[str(r["research_id"]).strip()] = {
                "group": a if a in GROUPS[:-1] else "OTH",
                "sex_at_birth": (r.get("sex_at_birth") or "").strip(),
                "complete": all(_num(r.get(c)) is not None for c in ("age", "coverage", *[f"PC{i}" for i in range(1, 11)])),
            }
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--selection-dir", required=True, type=Path, help="select_phenotypes.py --out-dir")
    ap.add_argument("--phenotype-csv", required=True, type=Path)
    ap.add_argument("--covariates-csv", required=True, type=Path)
    ap.add_argument("--phecode-sex-tsv", required=True, type=Path)
    ap.add_argument("--global-anc", required=True, type=Path, help="chr1 FLARE global ancestry (SAMPLE eas amr eur afr sas)")
    ap.add_argument("--out-counts", required=True, type=Path)
    ap.add_argument("--out-summary", required=True, type=Path)
    ap.add_argument("--legacy-sex-at-birth", action="store_true")
    args = ap.parse_args()

    cov = load_covariates(args.covariates_csv)
    pheno = pd.read_csv(args.phenotype_csv, dtype=str, low_memory=False)
    pheno = pheno.drop(columns=[c for c in pheno.columns if c.startswith("Unnamed") or c == ""], errors="ignore")
    pheno["person_id"] = pheno["person_id"].astype(str)
    pheno = pheno.set_index("person_id")

    if args.legacy_sex_at_birth:
        # first-version cohort: releasable PacBio in the phenotype table with a defined sex_at_birth and complete
        # age, coverage and PC1-10 (no inferred_sex requirement). Sex-specific phecodes restrict on sex_at_birth.
        ids = sorted(i for i, c in cov.items() if i in pheno.index and c["sex_at_birth"] in ("Female", "Male") and c["complete"])
        sex = pd.Series({i: (0.0 if cov[i]["sex_at_birth"] == "Female" else 1.0) for i in ids})
    else:
        ids = [s.strip() for s in (args.selection_dir / "analysis_samples.txt").read_text().splitlines() if s.strip()]
        sel_sex = pd.read_csv(args.selection_dir / "pheno_cov.tsv", sep="\t", usecols=["ID", "sex"], dtype={"ID": str}).set_index("ID")["sex"]
        missing = [i for i in ids if i not in cov]
        assert not missing, f"{len(missing)} analysis samples are not releasable PacBio in the covariates, e.g. {missing[:3]}"
        sex = sel_sex.reindex(ids)
        assert sex.notna().all(), "selector sex must be defined for every analysis sample"
    n = len(ids)

    pheno = pheno.loc[ids]
    cols = [c for c in pheno.columns if c != "sex"]  # the selector's own candidate set: every column but the id and sex

    phecode_sex = sp.read_phecode_sex(args.phecode_sex_tsv)
    anc = pd.read_csv(args.global_anc, sep="\t", dtype={"SAMPLE": str}).set_index("SAMPLE").rename(columns=ANC_COLS)
    prop = anc.loc[ids, LOCAL].to_numpy(float)
    prop = prop / prop.sum(axis=1, keepdims=True)  # the file rounds to 3 decimals; each participant's proportions sum to 1
    group = np.array([cov[i]["group"] for i in ids])
    gmat = np.stack([(group == g) for g in GROUPS], axis=1).astype(float)  # n x 7

    rows = []
    for col in cols:
        code = sp.phecode_of(col)
        enc = sp.encode_binary(pheno[col])
        side = phecode_sex.get(code, "")
        assert side in ("Male", "Female", "Both"), f"{code} is not in the phecode definitions table"
        if side in ("Male", "Female"):
            enc = enc.mask(sex.reindex(enc.index) == (0.0 if side == "Male" else 1.0), pd.NA)
        case = (enc == 1.0).fillna(False).astype(bool).to_numpy(float)  # NA (the other sex) is neither case nor control
        rows.append((code, int(case.sum()), int(enc.notna().sum()), case @ gmat, 2.0 * (case @ prop)))

    args.out_counts.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out_counts, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["phecode", "n_cases", "n_eligible"] + [f"cases_{g}" for g in GROUPS] + [f"hap_{k}" for k in LOCAL])
        for code, nc, ne, cg, hp in rows:
            assert int(cg.sum()) == nc, f"{code}: group cases {cg.sum()} != pooled {nc}"
            w.writerow([code, nc, ne] + [int(x) for x in cg] + [f"{x:.1f}" for x in hp])

    sex_arr = sex.to_numpy()
    summary = {
        "n": n,
        "n_female": int((sex_arr == 0.0).sum()),
        "n_male": int((sex_arr == 1.0).sum()),
        "n_group": {g: int((group == g).sum()) for g in GROUPS},
        "hap_total": {k: round(float(2.0 * prop[:, j].sum()), 1) for j, k in enumerate(LOCAL)},
        "source": (
            "select_phenotypes.py cohort (releasable PacBio x phenotype table x full-model covariates; "
            + ("defined sex_at_birth, legacy" if args.legacy_sex_at_birth else "selector sex")
            + ") x chr1 FLARE global ancestry; counts after the sex-specific restriction"
        ),
    }
    args.out_summary.write_text(json.dumps(summary, indent=1) + "\n")
    print(f"{n:,} participants, {len(rows):,} phecodes -> {args.out_counts}")


if __name__ == "__main__":
    main()
