#!/usr/bin/env python3
"""Merge long-read global PCs and HPRC/HGSVC3 control rows into covariates.

Adds:
  - lr_PC1–lr_PC32 for every sample in global_pcs.tsv (AoU + controls)
  - has_lr_pcs / is_reference_control flags
  - HG/NA reference-control rows from control_sample_metadata.tsv (ancestry,
    sex, Phase-1 SV counts where available; AoU phenotype / CDR fields left
    missing or False). Controls in metadata but not in global_pcs.tsv are
    still appended, with lr_PC* missing and has_lr_pcs=False.

Does not overwrite existing short-read PC1–PC32.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

LR_PC_COLS = [f"lr_PC{i}" for i in range(1, 33)]

# Boolean / flag columns that should be False for non-AoU reference controls
# unless we explicitly fill them.
CONTROL_FALSE_FLAGS = [
    "lr_meet_qc",
    "lr_releasable_v9",
    "rna_meet_qc",
    "has_rna",
    "proteomics_meet_qc",
    "has_proteomics",
    "has_exposomics",
    "has_srWGS",
    "has_demographics",
    "has_extraction_method",
    "has_trgt",
    "in_pedigree",
    "pedigree_is_founder",
    "has_global_pcs",  # short-read global PCs
    "has_pop_pcs",
    "in_cdr_v7",
    "in_cdr_v8",
    "in_cdr_v9",
    "has_methylation",
    "in_phenotype_table",
    "is_AIAN",
    "withdrawn",
    "unavailable",
    "final_releasable_v9",
    "in_final_set",
    "has_lr_tech",
    "has_ONT",
    "has_PacBio",
    "has_asm_metrics",
    "has_phase1_aux_metrics",
    "has_phase1_asm_metrics",
    "has_phase2_asm_metrics",
]


def _is_hg_na(series: pd.Series) -> pd.Series:
    return series.astype(str).str.match(r"^(HG|NA)\d", na=False)


def merge_lr_pcs(cov: pd.DataFrame, pcs: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    pcs = pcs.copy()
    pcs["research_id"] = pcs["research_id"].astype(str)
    missing = [c for c in ["research_id", *LR_PC_COLS] if c not in pcs.columns]
    if missing:
        raise ValueError(f"global_pcs.tsv missing columns: {missing}")

    # Drop any prior lr_PC / flag columns so re-runs are idempotent.
    drop_cols = [c for c in cov.columns if c in LR_PC_COLS or c in {"has_lr_pcs", "is_reference_control"}]
    out = cov.drop(columns=drop_cols, errors="ignore").copy()
    out["research_id"] = out["research_id"].astype(str)

    # Remove previously appended control rows before re-adding.
    out = out.loc[~_is_hg_na(out["research_id"])].copy()

    overlap = out.merge(pcs[["research_id", *LR_PC_COLS]], on="research_id", how="left", validate="one_to_one")
    return overlap, pcs


def build_control_rows(
    template_columns: list[str],
    pcs: pd.DataFrame,
    meta: pd.DataFrame,
) -> pd.DataFrame:
    meta = meta.copy()
    meta["research_id"] = meta["research_id"].astype(str)
    pcs_controls = pcs.loc[_is_hg_na(pcs["research_id"])].copy()
    pcs_ids = set(pcs_controls["research_id"])
    meta_ids = set(meta["research_id"])
    only_pcs = sorted(pcs_ids - meta_ids)
    if only_pcs:
        raise ValueError(
            "global_pcs HG/NA IDs missing from control metadata: "
            f"{only_pcs[:10]}"
        )

    rows = meta.merge(pcs_controls, on="research_id", how="left", validate="one_to_one")
    assert rows["ancestry_pred"].notna().all()

    out = pd.DataFrame({c: pd.NA for c in template_columns}, index=range(len(rows)))
    out["research_id"] = rows["research_id"].to_numpy()

    for col in LR_PC_COLS:
        out[col] = rows[col].to_numpy()

    out["is_reference_control"] = True
    out["has_lr_pcs"] = rows[LR_PC_COLS[0]].notna().to_numpy()
    out["ancestry_pred"] = rows["ancestry_pred"].to_numpy()
    out["ancestry_pred_other"] = rows["ancestry_pred_other"].to_numpy()
    out["has_ancestry_annotation"] = True
    out["sex_at_birth"] = rows["sex_at_birth"].to_numpy()
    out["inferred_sex"] = rows["inferred_sex"].to_numpy()
    # Continental codes only (AFR/AMR/EAS/EUR/SAS/OTH/MID), matching AoU
    # `population` and within-population PCA subset keys. Fine-grained 1KG/HPRC
    # codes stay in control_sample_metadata.tsv as population_code.
    if "population" in out.columns:
        # Prefer ancestry_pred_other (permits oth) then ancestry_pred.
        continental = rows["ancestry_pred_other"].fillna(rows["ancestry_pred"]).astype(str).str.upper()
        out["population"] = continental.to_numpy()

    out["has_sv_counts"] = rows["has_sv_counts"].fillna(False).astype(bool).to_numpy()
    out["n_sv_del"] = rows["n_sv_del"].to_numpy()
    out["n_sv_ins"] = rows["n_sv_ins"].to_numpy()
    out["n_trgt_vcfs"] = 0
    out["pedigree_n_parents"] = 0
    out["n_final_rows"] = 0
    out["n_lr_technical_rows"] = 0

    for col in CONTROL_FALSE_FLAGS:
        if col in out.columns:
            out[col] = False

    # Nullable bool in the AoU table.
    if "has_ehr_data" in out.columns:
        out["has_ehr_data"] = pd.NA

    return out


def order_columns(columns: list[str]) -> list[str]:
    """Insert new columns next to related existing fields."""
    cols = list(columns)
    for c in LR_PC_COLS + ["has_lr_pcs", "is_reference_control"]:
        if c in cols:
            cols.remove(c)

    if "PC32" in cols:
        idx = cols.index("PC32") + 1
        for i, c in enumerate(LR_PC_COLS):
            cols.insert(idx + i, c)
    else:
        cols.extend(LR_PC_COLS)

    if "has_global_pcs" in cols:
        idx = cols.index("has_global_pcs") + 1
        cols.insert(idx, "has_lr_pcs")
    else:
        cols.append("has_lr_pcs")

    if "has_ancestry_annotation" in cols:
        idx = cols.index("has_ancestry_annotation") + 1
        cols.insert(idx, "is_reference_control")
    else:
        cols.append("is_reference_control")

    return cols


def update_data_dictionary(dd: pd.DataFrame, final_columns: list[str]) -> pd.DataFrame:
    by_col = {r["column"]: r for _, r in dd.iterrows()}

    for i, c in enumerate(LR_PC_COLS, start=1):
        by_col[c] = {
            "column": c,
            "source": "tractor_mix/pca/deepvariant_lr_v1/global_pcs.tsv (DeepVariant long-read global PCA)",
            "notes": f"long-read global PC{i}; includes HPRC/HGSVC3 HG/NA controls; does not replace short-read PC{i}",
        }
    by_col["has_lr_pcs"] = {
        "column": "has_lr_pcs",
        "source": "derived",
        "notes": "True if matched to a long-read global PCA row in global_pcs.tsv",
    }
    by_col["is_reference_control"] = {
        "column": "is_reference_control",
        "source": "derived",
        "notes": "True for HPRC/HGSVC3/GIAB reference controls with HG*/NA* sample IDs",
    }

    # Clarify research_id now includes HG/NA controls.
    if "research_id" in by_col:
        by_col["research_id"] = {
            "column": "research_id",
            "source": "multiomics research_id (+ HG/NA reference-control sample IDs)",
            "notes": "primary key; AoU numeric research IDs plus HPRC/HGSVC3 HG*/NA* controls",
        }

    rows = []
    for c in final_columns:
        if c in by_col:
            rows.append(dict(by_col[c]))
        else:
            rows.append({"column": c, "source": "", "notes": ""})
    return pd.DataFrame(rows, columns=["column", "source", "notes"])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--covariates",
        type=Path,
        default=Path("tractor_mix/covariates.source_rebuilt.csv.gz"),
    )
    parser.add_argument(
        "--global-pcs",
        type=Path,
        default=Path("tractor_mix/pca/deepvariant_lr_v1/global_pcs.tsv"),
    )
    parser.add_argument(
        "--control-metadata",
        type=Path,
        default=Path("tractor_mix/reference_controls/control_sample_metadata.tsv"),
    )
    parser.add_argument(
        "--data-dictionary",
        type=Path,
        default=Path("tractor_mix/covariates.source_rebuilt.data_dictionary.tsv"),
    )
    parser.add_argument(
        "--out-covariates",
        type=Path,
        default=None,
        help="Defaults to --covariates (in-place update).",
    )
    parser.add_argument(
        "--out-data-dictionary",
        type=Path,
        default=None,
        help="Defaults to --data-dictionary (in-place update).",
    )
    args = parser.parse_args()
    out_cov = args.out_covariates or args.covariates
    out_dd = args.out_data_dictionary or args.data_dictionary

    cov = pd.read_csv(args.covariates, low_memory=False)
    pcs = pd.read_csv(args.global_pcs, sep="\t")
    meta = pd.read_csv(args.control_metadata, sep="\t")
    dd = pd.read_csv(args.data_dictionary, sep="\t")

    n_before = len(cov)
    merged, pcs = merge_lr_pcs(cov, pcs)

    # Template columns after lr_PC insertion.
    template_cols = order_columns(list(merged.columns) + LR_PC_COLS + ["has_lr_pcs", "is_reference_control"])
    # Ensure merged has the new columns.
    for c in LR_PC_COLS:
        if c not in merged.columns:
            raise RuntimeError(f"missing {c} after merge")
    merged["has_lr_pcs"] = merged[LR_PC_COLS[0]].notna()
    merged["is_reference_control"] = False

    controls = build_control_rows(template_cols, pcs, meta)
    # Align control frame to merged columns then concat.
    for c in merged.columns:
        if c not in controls.columns:
            controls[c] = pd.NA
    controls = controls[merged.columns]

    out = pd.concat([merged, controls], ignore_index=True)
    out = out[order_columns(list(out.columns))]
    assert out["research_id"].is_unique
    assert _is_hg_na(out["research_id"]).sum() == len(meta)

    dd_out = update_data_dictionary(dd, list(out.columns))
    assert list(dd_out["column"]) == list(out.columns)

    out_cov.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_cov, index=False, compression="gzip")
    dd_out.to_csv(out_dd, sep="\t", index=False)

    n_lr = int(out["has_lr_pcs"].sum())
    n_ctrl = int(out["is_reference_control"].sum())
    n_ctrl_anc = int(out.loc[out["is_reference_control"], "has_ancestry_annotation"].sum())
    n_ctrl_pcs = int(out.loc[out["is_reference_control"], "has_lr_pcs"].sum())
    print(f"Wrote {out_cov}")
    print(f"Wrote {out_dd}")
    print(f"rows: {n_before:,} -> {len(out):,} (+{len(out) - n_before})")
    print(f"has_lr_pcs: {n_lr:,}")
    print(f"reference controls: {n_ctrl} (with lr_PCs: {n_ctrl_pcs}; ancestry filled: {n_ctrl_anc})")
    print(f"control ancestry: {out.loc[out['is_reference_control'], 'ancestry_pred'].value_counts().to_dict()}")


if __name__ == "__main__":
    main()
