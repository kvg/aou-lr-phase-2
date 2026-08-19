#!/usr/bin/env python3
"""Build a PLINK2 --keep file from analysis_samples ∩ PLINK FAM/PSAM IIDs.

Prints clear diagnostics when overlap is empty (the usual MakeGRM failure mode).
Writes a ``#FID IID`` keep file with FID=IID (required for ``plink2 --double-id``
VCF imports; a lone ``#IID`` column often yields ``--keep: 0 samples remaining``).
"""

from __future__ import annotations

import argparse
from pathlib import Path


def read_ids(path: Path) -> list[str]:
    ids: list[str] = []
    text = path.read_text(encoding="utf-8-sig")  # strip BOM if present
    for line in text.splitlines():
        s = line.strip().replace("\r", "").strip('"').strip("'")
        if not s or s.startswith("#"):
            continue
        # allow FID IID or single-column ID lists
        ids.append(s.split()[0].strip('"').strip("'"))
    return ids


def read_fam_iids(path: Path) -> list[str]:
    """Return IID column from .fam or .psam."""
    ids: list[str] = []
    text = path.read_text(encoding="utf-8-sig")
    lines = text.splitlines()
    if not lines:
        return ids
    if path.suffix == ".psam" or lines[0].startswith("#"):
        header = lines[0].lstrip("#").split()
        # IID may be col 0 or 1 depending on whether FID is present
        try:
            iid_idx = header.index("IID")
        except ValueError:
            iid_idx = 1 if "FID" in header else 0
        for line in lines[1:]:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) > iid_idx:
                ids.append(parts[iid_idx])
        return ids

    # Classic FAM: FID IID ...
    for line in lines:
        parts = line.split()
        if len(parts) >= 2:
            ids.append(parts[1])
        elif len(parts) == 1:
            ids.append(parts[0])
    return ids


def read_id_list_file(path: Path) -> list[str]:
    """Single-column sample IDs (e.g. bcftools query -l), or FAM/PSAM IIDs."""
    if path.suffix in {".fam", ".psam"} or (
        path.exists() and path.read_text(encoding="utf-8-sig", errors="replace")[:1] == "#"
    ):
        return read_fam_iids(path)
    # Heuristic: if every non-comment line has >=2 columns, treat as FAM; else col1
    ids: list[str] = []
    multi = 0
    single = 0
    rows: list[list[str]] = []
    text = path.read_text(encoding="utf-8-sig")
    for line in text.splitlines():
        s = line.strip().replace("\r", "")
        if not s or s.startswith("#"):
            continue
        parts = s.split()
        rows.append(parts)
        if len(parts) >= 2:
            multi += 1
        else:
            single += 1
    if multi > single:
        return [r[1] for r in rows if len(r) >= 2]
    return [r[0] for r in rows]


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--analysis-samples", required=True, type=Path)
    p.add_argument(
        "--fam",
        required=True,
        type=Path,
        help="PLINK .fam/.psam, or a single-column sample ID list (bcftools query -l).",
    )
    p.add_argument("--out-keep", required=True, type=Path)
    p.add_argument(
        "--out-samples",
        type=Path,
        default=None,
        help="Optional rewritten analysis sample list (intersection, FAM/VCF order).",
    )
    args = p.parse_args()

    as_path = args.analysis_samples
    as_size = as_path.stat().st_size if as_path.exists() else -1
    print(f"analysis_samples path: {as_path} (bytes={as_size})")
    if as_size >= 0:
        raw_head = as_path.read_bytes()[:200]
        print(f"analysis_samples raw head: {raw_head!r}")

    wanted = read_ids(as_path)
    fam_iids = read_id_list_file(args.fam)
    fam_set = set(fam_iids)
    wanted_set = set(wanted)

    # Preserve FAM order for stable GRM orientation
    keep = [s for s in fam_iids if s in wanted_set]
    missing_in_fam = [s for s in wanted if s not in fam_set]
    extra_in_fam = len(fam_iids) - len(keep)

    print(f"analysis_samples: {len(wanted):,}")
    print(f"plink fam/psam:   {len(fam_iids):,}")
    print(f"intersection:     {len(keep):,}")
    print(f"analysis not in fam: {len(missing_in_fam):,}")
    print(f"fam not in analysis: {extra_in_fam:,}")
    print("analysis_samples examples:", wanted[:3])
    print("fam IID examples:        ", fam_iids[:3])
    if missing_in_fam:
        print("missing-in-fam examples:", missing_in_fam[:5])
    if fam_iids and wanted and not keep:
        fam0 = fam_iids[0]
        want0 = wanted[0]
        print(
            "ID mismatch detail: "
            f"fam[0]={fam0!r} (len={len(fam0)}), "
            f"analysis[0]={want0!r} (len={len(want0)})"
        )

    if not wanted:
        raise SystemExit(
            "ERROR: analysis_samples file is empty (0 IDs after parsing). "
            "This is not a VCF/FAM ID mismatch — the Terra input File has no sample IDs. "
            "Check gsutil cat/wc -l on the analysis_samples URI in your inputs JSON, "
            "re-run notebooks/tractor_01_prepare_inputs.ipynb, and re-upload "
            "tractor_mix_pilot/analysis_samples.txt."
        )

    if not keep:
        raise SystemExit(
            "ERROR: 0 samples remain after intersecting analysis_samples with PLINK FAM. "
            "GRM VCF sample IDs must match analysis_samples.txt (from the prep notebook / "
            "chr22 FLARE VCF). Re-run prep against the same VCF sample IDs used in grm_vcfs, "
            "or fix the analysis_samples path in the Terra inputs JSON."
        )

    # FID=IID pairs for plink2 --double-id (VCF sample ID becomes both FID and IID)
    args.out_keep.parent.mkdir(parents=True, exist_ok=True)
    body = "\n".join(f"{s}\t{s}" for s in keep)
    args.out_keep.write_text("#FID\tIID\n" + body + "\n", encoding="utf-8")
    print(f"Wrote {args.out_keep} ({len(keep)} samples, #FID IID / double-id)")

    if args.out_samples is not None:
        args.out_samples.write_text("\n".join(keep) + "\n", encoding="utf-8")
        print(f"Wrote {args.out_samples}")


if __name__ == "__main__":
    main()
