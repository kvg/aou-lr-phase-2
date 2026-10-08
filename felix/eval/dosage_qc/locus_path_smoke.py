#!/usr/bin/env python3
"""End-to-end smoke for the locus-level repeat-dosage input, plus A-1 / A-2.

Drives the production chain on a synthetic cohort, up to the point FELIX step 2
would take over:

    joint phased VCF (GT:AN1:AN2)
      -> aggregate_repeat_loci.py --out-vcf --with-ancestry   (one record per locus)
      -> write_admixed_dosage_vcf.py                          (DS{k} / ANC{k})
      -> bcftools view -Oz + index                            (VCF is well formed)

and re-establishes the two invariants from felix/REPEAT_DOSAGE.md §4.1 at the
writer boundary, where they can be checked without running FELIX:

  A-1  equivalence. At a locus holding exactly one record with one whole-unit
       ALT and no missing calls, the locus path and the record-level path must
       hand FELIX the same DS{k} and ANC{k}. This is the case where the two
       codings coincide, so any difference is a bug in the new path rather than
       a property of locus aggregation.
  A-2  scale invariance. Multiplying every input dosage at a locus by a
       constant must leave DS{k} unchanged and scale INFO RU_SCALE by that
       constant. DS is what enters the score test, so this is the
       writer-boundary form of "score tests do not change under a linear
       rescaling of the dosage".

Rerun this after moving to a newer FELIX or editing either script.

Usage::

    python3 felix/eval/dosage_qc/locus_path_smoke.py --out-dir /tmp/smoke
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2].parent
SCRIPTS = ROOT / "scripts"
PERIOD = 3
MOTIF = "CAG"
NUM_ANCS = 3


def sh(cmd: list[str]) -> str:
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit(f"command failed ({r.returncode}): {' '.join(map(str, cmd))}\n{r.stderr}")
    return r.stdout


def build_cohort(out: Path, n_samples: int, n_loci: int, seed: int) -> dict:
    """Synthetic joint VCF + catalogs. Returns the per-locus layout."""
    rng = np.random.default_rng(seed)
    samples = [f"s{i:05d}" for i in range(n_samples)]

    header = [
        "##fileformat=VCFv4.2",
        "##contig=<ID=chr1,length=250000000>",
        '##FILTER=<ID=PASS,Description="All filters passed">',
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">',
        '##INFO=<ID=SVLEN,Number=A,Type=Integer,Description="Length">',
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=AN1,Number=1,Type=Integer,Description="Ancestry hap 1">',
        '##FORMAT=<ID=AN2,Number=1,Type=Integer,Description="Ancestry hap 2">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples),
    ]

    # Local ancestry is a property of the haplotype at the locus, so it is drawn
    # once per locus per sample and shared by that locus's records.
    rows: list[str] = []
    catalog: list[str] = []
    simple: list[str] = []
    layout: dict[str, dict] = {}

    for l in range(n_loci):
        base = 100_000 + l * 1_000
        lid = f"TR{l:04d}"
        catalog.append(f"chr1\t{base - 1}\t{base + 60}\t{lid}\t{MOTIF}")
        simple.append(f"chr1\t{base - 1}\t{base + 60}\t{PERIOD}\t20.0\t{MOTIF}")

        an1 = rng.integers(0, NUM_ANCS, n_samples)
        an2 = rng.integers(0, NUM_ANCS, n_samples)

        # A third of loci are the A-1 case: one record, one whole-unit ALT, no
        # missing calls. The rest mix record counts, impure alleles and dropouts.
        kind = "a1" if l % 3 == 0 else ("multi" if l % 3 == 1 else "impure")
        n_rec = 1 if kind in ("a1", "impure") else int(rng.integers(2, 4))
        allow_missing = kind != "a1"

        recs = []
        for r in range(n_rec):
            pos = base + 5 + r * 15
            if kind == "impure" and r == 0:
                delta = int(rng.choice([-7, 5, 7, -5]))       # not a multiple of 3
            else:
                delta = int(PERIOD * rng.choice([-3, -2, -1, 1, 2, 4]))
            if delta > 0:
                ref = "A"
                alt = "A" + (MOTIF * ((delta // len(MOTIF)) + 1))[:delta]
            else:
                ref = "A" + (MOTIF * ((-delta // len(MOTIF)) + 1))[:-delta]
                alt = "A"
            gt1 = rng.integers(0, 2, n_samples)
            gt2 = rng.integers(0, 2, n_samples)
            miss = rng.random(n_samples) < (0.02 if allow_missing else 0.0)
            cells = []
            for i in range(n_samples):
                if miss[i]:
                    cells.append(f".|.:{an1[i]}:{an2[i]}")
                else:
                    cells.append(f"{gt1[i]}|{gt2[i]}:{an1[i]}:{an2[i]}")
            rows.append(
                "\t".join(["chr1", str(pos), f"{lid}_r{r}", ref, alt, ".", "PASS", ".",
                           "GT:AN1:AN2", *cells])
            )
            recs.append({"pos": pos, "delta": delta})
        layout[lid] = {"kind": kind, "n_rec": n_rec, "records": recs}

    rows.sort(key=lambda s: int(s.split("\t")[1]))
    plain = out / "joint.vcf"
    plain.write_text("\n".join(header + rows) + "\n")
    gz = out / "joint.vcf.gz"
    sh(["bcftools", "view", "-Oz", "-o", str(gz), str(plain)])
    sh(["bcftools", "index", "-t", str(gz)])
    (out / "catalog.bed").write_text("\n".join(catalog) + "\n")
    (out / "simple_repeat.bed").write_text("\n".join(simple) + "\n")
    (out / "samples.txt").write_text("\n".join(samples) + "\n")
    (out / "layout.json").write_text(json.dumps(layout, indent=2) + "\n")
    return layout


def read_dosage(path: Path) -> tuple[dict[str, dict], list[str]]:
    """{site_id: {scale, info, ds: (n, K) array, anc: (n, K) array}}, sample names."""
    sites: dict[str, dict] = {}
    names: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith("##"):
            continue
        f = line.split("\t")
        if line.startswith("#CHROM"):
            names = f[9:]
            continue
        info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
        keys = f[8].split(":")
        ds = np.zeros((len(f) - 9, NUM_ANCS))
        anc = np.zeros((len(f) - 9, NUM_ANCS))
        for i, cell in enumerate(f[9:]):
            v = dict(zip(keys, cell.split(":")))
            for k in range(NUM_ANCS):
                ds[i, k] = float(v[f"DS{k + 1}"])
                anc[i, k] = float(v[f"ANC{k + 1}"])
        sites[f[2]] = {"scale": float(info["RU_SCALE"]), "info": info, "ds": ds, "anc": anc,
                       "pos": int(f[1])}
    return sites, names


def run_locus_path(work: Path, tag: str, *extra: str) -> Path:
    sh([sys.executable, str(SCRIPTS / "aggregate_repeat_loci.py"),
        "--vcf", str(work / "joint.vcf.gz"), "--catalog-bed", str(work / "catalog.bed"),
        "--chrom", "chr1", "--samples", str(work / "samples.txt"),
        "--with-ancestry", "--num-ancs", str(NUM_ANCS),
        "--out-vcf", str(work / f"{tag}.loci.vcf"),
        "--out-alleles", str(work / f"{tag}.alleles.tsv.gz"),
        "--out-summary", str(work / f"{tag}.summary.json")])
    dose = work / f"{tag}.dosage.vcf"
    sh([sys.executable, str(SCRIPTS / "write_admixed_dosage_vcf.py"),
        "--vcf", str(work / f"{tag}.loci.vcf"), "--num-ancs", str(NUM_ANCS),
        "--samples", str(work / "samples.txt"), "--max-hap-missing", "0.2",
        "--out", str(dose), "--stats-json", str(work / f"{tag}.dosage.stats.json"), *extra])
    return dose


def scale_locus_vcf(src: Path, dst: Path, factor: float) -> None:
    """Multiply every RU_DOSAGE (and SVLEN) at every locus by `factor`."""
    out = []
    for line in src.read_text().splitlines():
        if line.startswith("#"):
            out.append(line)
            continue
        f = line.split("\t")
        parts = []
        for kv in f[7].split(";"):
            if kv.startswith("RU_DOSAGE="):
                vals = [float(v) * factor for v in kv.split("=", 1)[1].split(",")]
                parts.append("RU_DOSAGE=" + ",".join(f"{v:.17g}" for v in vals))
            else:
                parts.append(kv)
        f[7] = ";".join(parts)
        out.append("\t".join(f))
    dst.write_text("\n".join(out) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--n-samples", type=int, default=500)
    ap.add_argument("--n-loci", type=int, default=210)
    ap.add_argument("--seed", type=int, default=20261007)
    a = ap.parse_args()
    work = Path(a.out_dir)
    work.mkdir(parents=True, exist_ok=True)

    layout = build_cohort(work, a.n_samples, a.n_loci, a.seed)
    report: dict[str, object] = {"n_samples": a.n_samples, "n_loci": a.n_loci, "seed": a.seed}

    # ---- smoke: the chain runs and produces a VCF bcftools accepts ---------
    dose = run_locus_path(work, "locus")
    gz = work / "locus.dosage.vcf.gz"
    sh(["bcftools", "view", "-Oz", "-o", str(gz), str(dose)])
    sh(["bcftools", "index", "-c", str(gz)])
    n_out = int(sh(["bcftools", "index", "-n", str(gz)]).strip())
    agg_summary = json.loads((work / "locus.summary.json").read_text())
    writer_stats = json.loads((work / "locus.dosage.stats.json").read_text())
    report["loci_written_by_aggregator"] = agg_summary["counts"]["loci_written"]
    report["records_in_dosage_vcf"] = n_out
    report["writer_stats"] = writer_stats
    report["aggregator_hap_missing_frac"] = agg_summary["hap_missing_frac"]
    report["aggregator_hap_ancestry_inconsistent"] = agg_summary["counts"]["hap_ancestry_inconsistent"]
    assert n_out == writer_stats["written"] > 0, (n_out, writer_stats)

    sites, names = read_dosage(dose)
    # |DS{k}| <= ANC{k} <= 2 everywhere, the invariant FELIX's readers assume.
    worst = max(float(np.max(np.abs(s["ds"]) - s["anc"])) for s in sites.values())
    report["max_ds_minus_anc"] = worst
    assert worst <= 1e-9, f"|DS| exceeds ANC by {worst}"
    assert all(float(np.max(s["anc"])) <= 2.0 for s in sites.values())

    # ---- A-1: single-record whole-unit loci agree with the record path -----
    ru_vcf = work / "record.ru.vcf"
    sh([sys.executable, str(SCRIPTS / "annotate_repeat_units.py"),
        "--vcf", str(work / "joint.vcf.gz"),
        "--simple-repeat-bed", str(work / "simple_repeat.bed"), "--out", str(ru_vcf)])
    rec_dose = work / "record.dosage.vcf"
    sh([sys.executable, str(SCRIPTS / "write_admixed_dosage_vcf.py"),
        "--vcf", str(ru_vcf), "--num-ancs", str(NUM_ANCS),
        "--samples", str(work / "samples.txt"), "--max-hap-missing", "0.2",
        "--out", str(rec_dose), "--stats-json", str(work / "record.dosage.stats.json")])
    rec_sites, rec_names = read_dosage(rec_dose)
    assert names == rec_names

    a1_ids = [lid for lid, v in layout.items() if v["kind"] == "a1"]
    rec_by_pos = {s["pos"]: s for s in rec_sites.values()}
    checked, max_ds_diff, max_anc_diff = 0, 0.0, 0.0
    for lid in a1_ids:
        if lid not in sites:
            continue
        pos = layout[lid]["records"][0]["pos"]
        if pos not in rec_by_pos:
            continue
        loc, rec = sites[lid], rec_by_pos[pos]
        # Compare unscaled dosages: the two paths pick their own per-locus scale.
        d = np.max(np.abs(loc["ds"] * loc["scale"] - rec["ds"] * rec["scale"]))
        n = np.max(np.abs(loc["anc"] - rec["anc"]))
        max_ds_diff = max(max_ds_diff, float(d))
        max_anc_diff = max(max_anc_diff, float(n))
        checked += 1
    report["a1_loci_checked"] = checked
    report["a1_max_abs_ds_diff"] = max_ds_diff
    report["a1_max_abs_anc_diff"] = max_anc_diff
    assert checked >= 0.5 * len(a1_ids), f"only {checked} of {len(a1_ids)} A-1 loci compared"
    assert max_ds_diff <= 1e-9 and max_anc_diff == 0.0, (max_ds_diff, max_anc_diff)

    # ---- A-2: DS invariant under a constant rescaling of the dosage -------
    factor = 7.5
    scale_locus_vcf(work / "locus.loci.vcf", work / "scaled.loci.vcf", factor)
    scaled_dose = work / "scaled.dosage.vcf"
    sh([sys.executable, str(SCRIPTS / "write_admixed_dosage_vcf.py"),
        "--vcf", str(work / "scaled.loci.vcf"), "--num-ancs", str(NUM_ANCS),
        "--samples", str(work / "samples.txt"), "--max-hap-missing", "0.2",
        "--out", str(scaled_dose)])
    scaled, _ = read_dosage(scaled_dose)
    assert set(scaled) == set(sites), "rescaling changed which loci are testable"
    ds_diff = max(float(np.max(np.abs(scaled[k]["ds"] - sites[k]["ds"]))) for k in sites)
    scale_err = max(abs(scaled[k]["scale"] / sites[k]["scale"] - factor) for k in sites)
    report["a2_factor"] = factor
    report["a2_max_abs_ds_diff"] = ds_diff
    report["a2_max_scale_ratio_err"] = scale_err
    assert ds_diff <= 1e-7, ds_diff
    assert scale_err <= 1e-6, scale_err

    (work / "locus_path_smoke.report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
