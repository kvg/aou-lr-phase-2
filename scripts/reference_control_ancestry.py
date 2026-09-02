#!/usr/bin/env python3
"""Curated ancestry labels for HPRC / GIAB / 1KG reference controls in joint VCFs."""

from __future__ import annotations

from pathlib import Path

# HPRC year-1 assembly metadata (Superpopulation column), mapped to AoU-style labels.
# Source: human-pangenomics/HPP_Year1_Assemblies sample_metadata.
_HG_SUPERPOP: dict[str, str] = {
    "HG01123": "amr",
    "HG01258": "amr",
    "HG01358": "amr",
    "HG01361": "amr",
    "HG01891": "afr",
    "HG02257": "afr",
    "HG02486": "afr",
    "HG02559": "afr",
    "HG02572": "afr",
    "HG03516": "afr",
    "HG00438": "eas",
    "HG00621": "eas",
    "HG00673": "eas",
    "HG00735": "amr",
    "HG00741": "amr",
    "HG01071": "amr",
    "HG01106": "amr",
    "HG01175": "amr",
    "HG01928": "amr",
    "HG01952": "amr",
    "HG01978": "amr",
    "HG02148": "amr",
    "HG02622": "afr",
    "HG02630": "afr",
    "HG02717": "afr",
    "HG02886": "afr",
    "HG03453": "afr",
    "HG03540": "afr",
    "HG03579": "afr",
    "HG002": "oth",  # GIAB Ashkenazi trio
    "HG003": "oth",
    "HG004": "oth",
    "HG005": "eas",  # GIAB Chinese trio
    "HG006": "eas",
    "HG007": "eas",
    "HG00733": "amr",
    "HG01109": "amr",
    "HG01243": "amr",
    "HG02080": "eas",
    "HG02109": "afr",
    "HG02145": "afr",
    "HG02723": "afr",
    "HG02818": "afr",
    "HG03486": "afr",
    "HG03492": "sas",
    "HG02055": "afr",
    "HG03098": "afr",
    "NA18906": "afr",
    "NA19240": "afr",
    "NA20129": "afr",
    "NA21309": "afr",
}

# Common Coriell aliases for GIAB / 1KG controls.
_NA_ALIASES: dict[str, str] = {
    "NA24385": "HG002",
    "NA24149": "HG003",
    "NA24143": "HG004",
    "NA12878": "eur",  # CEU benchmark; not in HPRC table above
}

REFERENCE_CONTROL_PREFIXES = ("HG", "NA")

# Expanded Phase-2 control table (HG/NA IDs) written beside covariates.
_METADATA_CANDIDATES = (
    Path(__file__).resolve().parent.parent / "tractor_mix/reference_controls/control_sample_metadata.tsv",
    Path("tractor_mix/reference_controls/control_sample_metadata.tsv"),
)


def _load_control_metadata() -> dict[str, tuple[str, str]]:
    for path in _METADATA_CANDIDATES:
        if not path.is_file():
            continue
        import pandas as pd

        df = pd.read_csv(path, sep="\t")
        out: dict[str, tuple[str, str]] = {}
        for _, row in df.iterrows():
            sid = str(row["research_id"]).strip()
            pred = row.get("ancestry_pred")
            other = row.get("ancestry_pred_other", pred)
            if pd.isna(pred):
                continue
            out[sid] = (str(pred), str(other) if not pd.isna(other) else str(pred))
        return out
    return {}


_CONTROL_METADATA_ANCESTRY = _load_control_metadata()


def is_reference_control(sample_id: str) -> bool:
    token = str(sample_id).strip()
    return token.startswith(REFERENCE_CONTROL_PREFIXES)


def lookup_reference_ancestry(sample_id: str) -> tuple[str, str] | None:
    """Return (ancestry_pred, ancestry_pred_other) for a known reference control."""
    token = str(sample_id).strip()
    if token in _CONTROL_METADATA_ANCESTRY:
        return _CONTROL_METADATA_ANCESTRY[token]
    if token in _HG_SUPERPOP:
        label = _HG_SUPERPOP[token]
        return label, label
    if token in _NA_ALIASES:
        mapped = _NA_ALIASES[token]
        if mapped in _HG_SUPERPOP:
            label = _HG_SUPERPOP[mapped]
            return label, label
        return mapped, mapped
    return None


def fill_reference_control_ancestry(
    df,
    *,
    id_column: str = "s",
    ancestry_pred_column: str = "ancestry_pred",
    ancestry_pred_other_column: str = "ancestry_pred_other",
):
    """Fill missing ancestry fields for HG/NA controls from curated public metadata."""
    import pandas as pd

    out = df.copy()
    if id_column not in out.columns:
        raise ValueError(f"missing {id_column}")

    for col in (ancestry_pred_column, ancestry_pred_other_column):
        if col not in out.columns:
            out[col] = pd.NA

    missing = out[ancestry_pred_other_column].isna() | (
        out[ancestry_pred_other_column].astype(str) == "NA"
    )
    is_control = out[id_column].map(is_reference_control).fillna(False)
    for idx in out.index[missing & is_control]:
        sample_id = str(out.at[idx, id_column])
        labels = lookup_reference_ancestry(sample_id)
        if labels is None:
            continue
        pred, pred_other = labels
        if pd.isna(out.at[idx, ancestry_pred_column]) or str(out.at[idx, ancestry_pred_column]) == "NA":
            out.at[idx, ancestry_pred_column] = pred
        if pd.isna(out.at[idx, ancestry_pred_other_column]) or str(out.at[idx, ancestry_pred_other_column]) == "NA":
            out.at[idx, ancestry_pred_other_column] = pred_other
    return out
