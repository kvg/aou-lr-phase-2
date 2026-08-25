#!/usr/bin/env python3
"""Merge long-read within-population PCs into covariates.

Adds:
  - lr_pop_PC1–lr_pop_PC32 from population_pcs.tsv
  - has_lr_pop_pcs flag
  - lr_pop_population (continental subset used for that sample's within-pop PCA)

Does not overwrite short-read pop_PC* / has_pop_pcs or covariates ``population``.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

LR_POP_PC_COLS = [f"lr_pop_PC{i}" for i in range(1, 33)]
NEW_COLS = [*LR_POP_PC_COLS, "has_lr_pop_pcs", "lr_pop_population"]


def order_columns(columns: list[str]) -> list[str]:
    cols = [c for c in columns if c not in NEW_COLS]

    if "lr_PC32" in cols:
        idx = cols.index("lr_PC32") + 1
        for i, c in enumerate(LR_POP_PC_COLS):
            cols.insert(idx + i, c)
    elif "PC32" in cols:
        idx = cols.index("PC32") + 1
        for i, c in enumerate(LR_POP_PC_COLS):
            cols.insert(idx + i, c)
    else:
        cols.extend(LR_POP_PC_COLS)

    if "has_lr_pcs" in cols:
        idx = cols.index("has_lr_pcs") + 1
        cols.insert(idx, "has_lr_pop_pcs")
    elif "has_pop_pcs" in cols:
        idx = cols.index("has_pop_pcs") + 1
        cols.insert(idx, "has_lr_pop_pcs")
    else:
        cols.append("has_lr_pop_pcs")

    if "population" in cols:
        idx = cols.index("population") + 1
        cols.insert(idx, "lr_pop_population")
    else:
        cols.append("lr_pop_population")

    return cols


def update_data_dictionary(dd: pd.DataFrame, final_columns: list[str]) -> pd.DataFrame:
    by_col = {r["column"]: dict(r) for _, r in dd.iterrows()}
    for i, c in enumerate(LR_POP_PC_COLS, start=1):
        by_col[c] = {
            "column": c,
            "source": "population_pcs.tsv (DeepVariant long-read within-population PCA)",
            "notes": (
                f"long-read within-population PC{i}; computed within ancestry_pred_other "
                "subsets; does not replace short-read pop_PC*"
            ),
        }
    by_col["has_lr_pop_pcs"] = {
        "column": "has_lr_pop_pcs",
        "source": "derived",
        "notes": "True if matched to a long-read within-population PCA row in population_pcs.tsv",
    }
    by_col["lr_pop_population"] = {
        "column": "lr_pop_population",
        "source": "population_pcs.tsv population",
        "notes": (
            "continental subset used for this sample's long-read within-pop PCA "
            "(lowercase afr/amr/eas/eur/sas/mid/oth); usually equals ancestry_pred_other"
        ),
    }
    rows = []
    for c in final_columns:
        if c in by_col:
            rows.append(by_col[c])
        else:
            rows.append({"column": c, "source": "", "notes": ""})
    return pd.DataFrame(rows, columns=["column", "source", "notes"])


def merge_lr_pop_pcs(cov: pd.DataFrame, pcs: pd.DataFrame) -> pd.DataFrame:
    pcs = pcs.copy()
    pcs["research_id"] = pcs["research_id"].astype(str)
    need = ["research_id", "population", *LR_POP_PC_COLS]
    missing = [c for c in need if c not in pcs.columns]
    if missing:
        raise ValueError(f"population_pcs.tsv missing columns: {missing}")
    if pcs["research_id"].duplicated().any():
        raise ValueError("population_pcs.tsv has duplicate research_id values")

    drop_cols = [c for c in cov.columns if c in NEW_COLS]
    out = cov.drop(columns=drop_cols, errors="ignore").copy()
    out["research_id"] = out["research_id"].astype(str)

    pcs_join = pcs[["research_id", "population", *LR_POP_PC_COLS]].rename(
        columns={"population": "lr_pop_population"}
    )
    pcs_join["lr_pop_population"] = pcs_join["lr_pop_population"].astype(str).str.lower()

    out = out.merge(pcs_join, on="research_id", how="left", validate="one_to_one")
    out["has_lr_pop_pcs"] = out[LR_POP_PC_COLS[0]].notna()
    return out[order_columns(list(out.columns))]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--covariates",
        type=Path,
        default=Path("tractor_mix/covariates.source_rebuilt.csv.gz"),
    )
    parser.add_argument(
        "--population-pcs",
        type=Path,
        default=Path("tractor_mix/pca/deepvariant_lr_v1/population_pcs.tsv"),
    )
    parser.add_argument(
        "--data-dictionary",
        type=Path,
        default=Path("tractor_mix/covariates.source_rebuilt.data_dictionary.tsv"),
    )
    parser.add_argument("--out-covariates", type=Path, default=None)
    parser.add_argument("--out-data-dictionary", type=Path, default=None)
    args = parser.parse_args()
    out_cov = args.out_covariates or args.covariates
    out_dd = args.out_data_dictionary or args.data_dictionary

    cov = pd.read_csv(args.covariates, dtype={"research_id": str}, low_memory=False)
    pcs = pd.read_csv(args.population_pcs, sep="\t", dtype={"research_id": str})
    dd = pd.read_csv(args.data_dictionary, sep="\t")

    n_before = len(cov)
    out = merge_lr_pop_pcs(cov, pcs)
    dd_out = update_data_dictionary(dd, list(out.columns))
    assert list(dd_out["column"]) == list(out.columns)
    assert out["research_id"].is_unique

    out.to_csv(out_cov, index=False, compression="gzip")
    dd_out.to_csv(out_dd, sep="\t", index=False)

    n_hit = int(out["has_lr_pop_pcs"].sum())
    print(f"Wrote {out_cov}")
    print(f"Wrote {out_dd}")
    print(f"rows: {n_before:,} (unchanged {len(out):,})")
    print(f"has_lr_pop_pcs: {n_hit:,}")
    print("lr_pop_population counts:")
    print(out.loc[out["has_lr_pop_pcs"], "lr_pop_population"].value_counts().to_string())


if __name__ == "__main__":
    main()
