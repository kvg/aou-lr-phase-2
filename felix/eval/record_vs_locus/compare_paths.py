#!/usr/bin/env python3
"""Record-level vs locus-level repeat dosage on the defect fixture.

Runs both candidate association inputs over
``scripts/testdata/repeat_locus_assoc/`` and compares the per-ancestry dosage
``D_ik`` each one hands to FELIX against the fixture's hand-written truth:

  record path  annotate_repeat_units.py --simple-repeat-bed
                 -> write_admixed_dosage_vcf.py        (one test per VCF record)
  locus path   aggregate_repeat_loci.py --out-vcf --with-ancestry
                 -> write_admixed_dosage_vcf.py        (one test per catalog locus)

Outputs (``--out-dir``):
  expected_haplotype_dosage.tsv   the fixture truth, per locus/sample/haplotype
  per_sample_ancestry_dosage.tsv  D_ik: truth vs locus path vs summed record path
  defect_summary.tsv              per-defect verdict and counts
  record_vs_locus.png             produced dosage against expected dosage
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

SAMPLES = ["S1", "S2", "S3", "S4", "S5", "S6"]
# AN1, AN2 per sample; constant across records in the fixture.
ANCESTRY = {"S1": (0, 1), "S2": (0, 0), "S3": (1, 1), "S4": (0, 1), "S5": (1, 0), "S6": (0, 0)}
PERIOD = 3
NUM_ANCS = 2

# Fixture truth. Per locus: the defect it exercises, and the signed bp change
# each haplotype carries, summed over the locus's records.
#   locus -> (defect, {sample: (hap1_bp, hap2_bp)})
TRUTH = {
    "TR_L1": (
        "1_split_across_records",
        {"S1": (66, 0), "S2": (0, 6), "S3": (0, 60), "S4": (6, 6), "S5": (60, 60), "S6": (0, 0)},
    ),
    "TR_L2": (
        "2_impure_allele",
        {"S1": (6, 0), "S2": (7, 0), "S3": (6, 7), "S4": (0, 0), "S5": (0, 0), "S6": (0, 0)},
    ),
    "TR_L3": (
        "3_seqres_deletion_no_svtype",
        {"S1": (0, 0), "S2": (-6, 0), "S3": (0, -6), "S4": (-6, -6), "S5": (0, 0), "S6": (0, 0)},
    ),
    "TR_L4": (
        "4_mixed_symbolic_record",
        {"S1": (-9, 0), "S2": (9, 0), "S3": (-9, 9), "S4": (0, 0), "S5": (0, 0), "S6": (0, 0)},
    ),
    "TR_L5": (
        "4b_mixed_seqresolved_record",
        {"S1": (-6, 0), "S2": (6, 0), "S3": (-6, 6), "S4": (0, 0), "S5": (0, 0), "S6": (0, 0)},
    ),
}
# Catalog interval each locus occupies, to map a record-path site back to a locus.
LOCUS_SPAN = {
    "TR_L1": (1000, 1060),
    "TR_L2": (2000, 2060),
    "TR_L3": (3000, 3060),
    "TR_L4": (4000, 4060),
    "TR_L5": (5000, 5060),
}


def sh(cmd: list[str], **kw) -> subprocess.CompletedProcess:
    r = subprocess.run(cmd, capture_output=True, text=True, **kw)
    if r.returncode != 0:
        raise SystemExit(f"command failed ({r.returncode}): {' '.join(cmd)}\n{r.stderr}")
    return r


def expected_frames() -> tuple[pd.DataFrame, pd.DataFrame]:
    """(per-haplotype truth, per-sample-per-ancestry truth) in repeat units."""
    hap_rows = []
    for locus, (defect, per_sample) in TRUTH.items():
        for s, (h1, h2) in per_sample.items():
            for hap, (bp, anc) in enumerate(zip((h1, h2), ANCESTRY[s]), start=1):
                hap_rows.append(
                    {"locus": locus, "defect": defect, "sample": s, "hap": hap,
                     "ancestry": anc, "dosage_bp": bp, "expected_units": bp / PERIOD}
                )
    hap = pd.DataFrame(hap_rows)
    dik = (
        hap.groupby(["locus", "defect", "sample", "ancestry"], as_index=False)["expected_units"]
        .sum()
        .rename(columns={"expected_units": "expected"})
    )
    return hap, dik


def read_dosage_vcf(path: Path) -> pd.DataFrame:
    """Per-sample per-ancestry dosage from a FELIX admixed dosage VCF, unscaled."""
    rows = []
    names: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith("##"):
            continue
        if line.startswith("#CHROM"):
            names = line.split("\t")[9:]
            continue
        f = line.split("\t")
        info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
        scale = float(info["RU_SCALE"])
        keys = f[8].split(":")
        for name, cell in zip(names, f[9:]):
            vals = dict(zip(keys, cell.split(":")))
            for k in range(1, NUM_ANCS + 1):
                rows.append({"pos": int(f[1]), "site": f[2], "sample": name, "ancestry": k - 1,
                             "produced": float(vals[f"DS{k}"]) * scale,
                             "anc_count": int(vals[f"ANC{k}"])})
    return pd.DataFrame(rows)


def locus_of_pos(pos: int) -> str | None:
    for locus, (lo, hi) in LOCUS_SPAN.items():
        if lo <= pos <= hi:
            return locus
    return None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fixture-dir", required=True)
    ap.add_argument("--scripts-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--work-dir", required=True)
    a = ap.parse_args()
    fx, sc = Path(a.fixture_dir), Path(a.scripts_dir)
    out, work = Path(a.out_dir), Path(a.work_dir)
    out.mkdir(parents=True, exist_ok=True)
    work.mkdir(parents=True, exist_ok=True)

    gz = work / "joint.vcf.gz"
    sh(["bash", "-c", f"bgzip -c {fx / 'joint.vcf'} > {gz} && bcftools index -t -f {gz}"])

    # ---- record path -------------------------------------------------------
    ru_vcf = work / "record.ru.vcf"
    sh([sys.executable, str(sc / "annotate_repeat_units.py"),
        "--vcf", str(gz), "--simple-repeat-bed", str(fx / "simple_repeat.bed"),
        "--out", str(ru_vcf)])
    rec_dose = work / "record.dosage.vcf"
    rec_err = ""
    try:
        sh([sys.executable, str(sc / "write_admixed_dosage_vcf.py"),
            "--vcf", str(ru_vcf), "--num-ancs", str(NUM_ANCS), "--out", str(rec_dose)])
    except SystemExit as e:
        rec_err = str(e)

    # ---- locus path --------------------------------------------------------
    loc_vcf = work / "locus.vcf"
    sh([sys.executable, str(sc / "aggregate_repeat_loci.py"),
        "--vcf", str(gz), "--catalog-bed", str(fx / "catalog.bed"), "--chrom", "chr22",
        "--with-ancestry", "--out-vcf", str(loc_vcf),
        "--out-alleles", str(work / "locus.alleles.tsv"),
        "--out-summary", str(work / "locus.summary.json")])
    loc_dose = work / "locus.dosage.vcf"
    sh([sys.executable, str(sc / "write_admixed_dosage_vcf.py"),
        "--vcf", str(loc_vcf), "--num-ancs", str(NUM_ANCS), "--out", str(loc_dose)])

    # ---- compare -----------------------------------------------------------
    hap, dik = expected_frames()
    hap.to_csv(out / "expected_haplotype_dosage.tsv", sep="\t", index=False)

    loc = read_dosage_vcf(loc_dose)
    loc["locus"] = loc["site"]
    loc_dik = loc.groupby(["locus", "sample", "ancestry"], as_index=False)["produced"].sum()
    loc_dik = loc_dik.rename(columns={"produced": "locus_path"})
    n_loc_tests = loc.groupby("locus")["site"].nunique()

    if rec_err:
        rec_dik = pd.DataFrame(columns=["locus", "sample", "ancestry", "record_path"])
        n_rec_tests = pd.Series(dtype=int)
        rec_sites: list[int] = []
    else:
        rec = read_dosage_vcf(rec_dose)
        rec["locus"] = rec["pos"].map(locus_of_pos)
        rec_sites = sorted(rec["pos"].unique())
        rec_dik = (
            rec.groupby(["locus", "sample", "ancestry"], as_index=False)["produced"]
            .sum().rename(columns={"produced": "record_path"})
        )
        n_rec_tests = rec.groupby("locus")["pos"].nunique()

    cmp = dik.merge(loc_dik, on=["locus", "sample", "ancestry"], how="left")
    cmp = cmp.merge(rec_dik, on=["locus", "sample", "ancestry"], how="left")
    cmp["locus_path"] = cmp["locus_path"].fillna(0.0)
    cmp["locus_ok"] = np.isclose(cmp["locus_path"], cmp["expected"], atol=1e-6)
    cmp["record_tested"] = cmp["locus"].map(lambda l: int(n_rec_tests.get(l, 0)) > 0)
    cmp["record_ok"] = np.where(
        cmp["record_tested"], np.isclose(cmp["record_path"].fillna(0.0), cmp["expected"], atol=1e-6), False
    )
    cmp.to_csv(out / "per_sample_ancestry_dosage.tsv", sep="\t", index=False)

    # Defect-level verdict.
    srows = []
    for locus, (defect, _) in TRUTH.items():
        sub = cmp[cmp["locus"] == locus]
        nz = sub[sub["expected"] != 0]
        n_rec = int(n_rec_tests.get(locus, 0))
        srows.append({
            "locus": locus,
            "defect": defect,
            "tests_locus_path": int(n_loc_tests.get(locus, 0)),
            "tests_record_path": n_rec,
            "nonzero_dik": int(len(nz)),
            "locus_path_correct": int(nz["locus_ok"].sum()),
            "record_path_correct": int(nz["record_ok"].sum()),
            "verdict": (
                "record path drops the locus entirely" if n_rec == 0
                else "record path splits the locus into separate tests" if n_rec > 1
                else "record path agrees with the locus path"
                if bool(nz["record_ok"].all())
                else "record path mis-signs the dosage"
            ),
        })
    summary = pd.DataFrame(srows).sort_values("defect")
    summary.to_csv(out / "defect_summary.tsv", sep="\t", index=False)

    with open(out / "run_notes.txt", "w") as fh:
        fh.write(f"record_path_sites_emitted\t{len(rec_sites)}\t{rec_sites}\n")
        fh.write(f"locus_path_sites_emitted\t{int(loc['site'].nunique())}\n")
        fh.write(f"record_path_writer_error\t{rec_err}\n")

    print(summary.to_string(index=False))
    print("\nrecord-path sites:", rec_sites)
    print("writer error:", rec_err or "(none)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
