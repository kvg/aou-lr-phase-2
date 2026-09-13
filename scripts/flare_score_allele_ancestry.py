#!/usr/bin/env python3
"""Allele–ancestry concordance for FLARE recipes (association screen, Part 2).

For each fixed-panel marker, read alleles from the unfiltered phased ``gt_vcf``
and map ancestry from the recipe ``anc.vcf`` with FELIXla / propagate covering
semantics (first FLARE POS >= query, else last on contig). Score each haplotype
allele under the panel AF for the inferred ancestry.

Primary gate metric: ``mean_ll`` (mean per-haplotype log-likelihood).
Secondary: ``mean_brier`` (mean (1 - p_correct)^2).

Example::

    python3 scripts/flare_score_allele_ancestry.py \\
      --anc-vcf recipe.anc.vcf.gz \\
      --gt-vcf chr22.phased.vcf.gz \\
      --panel eval_panels/chr22_10mb.markers.tsv \\
      --samples analysis_samples.txt \\
      --out allele_ancestry_score.json
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
from pathlib import Path
from typing import Any, Iterable, Optional, Sequence

# Match scripts/flare_switch_qc.py
ANCESTRY = {0: "eas", 1: "amr", 2: "eur", 3: "afr", 4: "sas"}
ANCESTRY_TO_INT = {v: k for k, v in ANCESTRY.items()}
PANELS = ("eas", "amr", "eur", "afr", "sas")


def load_panel(path: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    with path.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < len(header):
                continue
            rec = dict(zip(header, parts))
            afs: dict[str, float] = {}
            for p in PANELS:
                key = f"af_{p}"
                raw = rec.get(key, "")
                afs[p] = float("nan") if raw in {"", ".", "nan", "NA"} else float(raw)
            rows.append(
                {
                    "chrom": rec["chrom"],
                    "pos": int(rec["pos"]),
                    "ref": rec["ref"],
                    "alt": rec["alt"],
                    "afs": afs,
                }
            )
    return rows


def load_sample_list(path: Optional[Path]) -> Optional[list[str]]:
    if path is None:
        return None
    return [ln.strip() for ln in path.read_text().splitlines() if ln.strip()]


def parse_gt_alleles(gt: str) -> Optional[tuple[int, int]]:
    g = gt.replace("/", "|")
    if "|" not in g:
        return None
    a, b = g.split("|", 1)
    if a in {".", ""} or b in {".", ""}:
        return None
    try:
        return int(a), int(b)
    except ValueError:
        return None


def parse_an(val: str) -> Optional[int]:
    if val in {".", "", "NA"}:
        return None
    try:
        v = int(val)
    except ValueError:
        return None
    return v if v in ANCESTRY else None


def covering_index(lai_positions: Sequence[int], query_pos: int) -> Optional[int]:
    """FELIXla covering: first index with pos >= query, else last."""
    if not lai_positions:
        return None
    lo, hi = 0, len(lai_positions)
    while lo < hi:
        mid = (lo + hi) // 2
        if lai_positions[mid] < query_pos:
            lo = mid + 1
        else:
            hi = mid
    if lo < len(lai_positions):
        return lo
    return len(lai_positions) - 1


def allele_loglik(allele: int, af_alt: float, *, eps: float = 1e-6) -> float:
    """Log P(allele | AF_alt) for a biallelic site."""
    if math.isnan(af_alt):
        return float("nan")
    p = min(max(af_alt, eps), 1.0 - eps)
    if allele == 1:
        return math.log(p)
    if allele == 0:
        return math.log(1.0 - p)
    return float("nan")


def allele_brier(allele: int, af_alt: float) -> float:
    if math.isnan(af_alt):
        return float("nan")
    p_correct = af_alt if allele == 1 else (1.0 - af_alt)
    p_correct = min(max(p_correct, 0.0), 1.0)
    return (1.0 - p_correct) ** 2


def score_haplotype(
    allele: int,
    anc: Optional[int],
    afs: dict[str, float],
) -> tuple[float, float, bool]:
    """Return (ll, brier, scored)."""
    if anc is None or allele not in (0, 1):
        return float("nan"), float("nan"), False
    label = ANCESTRY.get(anc)
    if label is None:
        return float("nan"), float("nan"), False
    af = afs.get(label, float("nan"))
    ll = allele_loglik(allele, af)
    br = allele_brier(allele, af)
    if math.isnan(ll):
        return float("nan"), float("nan"), False
    return ll, br, True


def bcftools_query(args: list[str]) -> Iterable[str]:
    cmd = ["bcftools", "query", *args]
    print("+", " ".join(cmd), file=sys.stderr)
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True, bufsize=1)
    assert proc.stdout is not None
    for line in proc.stdout:
        yield line
    rc = proc.wait()
    if rc != 0:
        raise SystemExit(f"bcftools query failed ({rc}): {' '.join(cmd)}")


def load_anc_sites(
    anc_vcf: str,
    samples: list[str],
    *,
    region: str = "",
) -> tuple[list[int], list[list[tuple[Optional[int], Optional[int]]]]]:
    """Return (positions, per_site list of (AN1, AN2) per sample index)."""
    sample_arg: list[str] = []
    tmp: Optional[Path] = None
    if samples:
        tmp = Path(anc_vcf + ".score_samples.txt")
        # Prefer stdin-less file next to out; use NamedTemp via /tmp
        import tempfile

        tf = tempfile.NamedTemporaryFile("w", delete=False, suffix=".samples.txt")
        tf.write("\n".join(samples) + "\n")
        tf.close()
        tmp = Path(tf.name)
        sample_arg = ["-S", str(tmp)]
    cmd = [*sample_arg, "-f", r"%POS[\t%AN1\t%AN2]\n"]
    if region.strip():
        cmd.extend(["-r", region.strip()])
    cmd.append(anc_vcf)
    positions: list[int] = []
    ancs: list[list[tuple[Optional[int], Optional[int]]]] = []
    n_samp = len(samples)
    for line in bcftools_query(cmd):
        parts = line.rstrip("\n").split("\t")
        if not parts:
            continue
        pos = int(parts[0])
        vals = parts[1:]
        if n_samp and len(vals) != 2 * n_samp:
            continue
        row: list[tuple[Optional[int], Optional[int]]] = []
        for i in range(0, len(vals), 2):
            row.append((parse_an(vals[i]), parse_an(vals[i + 1])))
        positions.append(pos)
        ancs.append(row)
    if tmp is not None:
        tmp.unlink(missing_ok=True)
    return positions, ancs


def load_gt_alleles_at_panel(
    gt_vcf: str,
    samples: list[str],
    panel: list[dict[str, Any]],
    *,
    region: str = "",
) -> dict[tuple[str, int], list[Optional[tuple[int, int]]]]:
    """Map (chrom,pos) -> per-sample (a1,a2) or None if missing / allele mismatch."""
    if not panel:
        return {}
    import tempfile

    # Restrict query to panel positions via regions file when possible
    with tempfile.NamedTemporaryFile("w", delete=False, suffix=".samples.txt") as sf:
        sf.write("\n".join(samples) + "\n")
        s_path = sf.name
    with tempfile.NamedTemporaryFile("w", delete=False, suffix=".regions.txt") as rf:
        for m in panel:
            rf.write(f"{m['chrom']}\t{m['pos']}\n")
        r_path = rf.name

    want = {(m["chrom"], m["pos"], m["ref"], m["alt"]) for m in panel}
    want_pos = {(m["chrom"], m["pos"]) for m in panel}
    out: dict[tuple[str, int], list[Optional[tuple[int, int]]]] = {}

    cmd = [
        "-S",
        s_path,
        "-R",
        r_path,
        "-f",
        r"%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n",
    ]
    if region.strip():
        # -R already limits; optional -r further restricts
        cmd.extend(["-r", region.strip()])
    cmd.append(gt_vcf)
    try:
        for line in bcftools_query(cmd):
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4 + len(samples):
                continue
            chrom, pos_s, ref, alt = parts[0], parts[1], parts[2], parts[3]
            pos = int(pos_s)
            key = (chrom, pos)
            if key not in want_pos:
                continue
            # Accept only exact REF/ALT match to panel (first ALT if multi)
            alt0 = alt.split(",")[0]
            panel_match = (chrom, pos, ref, alt0) in want
            gts = parts[4 : 4 + len(samples)]
            row: list[Optional[tuple[int, int]]] = []
            for gt in gts:
                al = parse_gt_alleles(gt) if panel_match else None
                row.append(al)
            out[key] = row
    finally:
        Path(s_path).unlink(missing_ok=True)
        Path(r_path).unlink(missing_ok=True)
    return out


def score_recipe(
    *,
    panel: list[dict[str, Any]],
    lai_positions: Sequence[int],
    lai_ancs: Sequence[Sequence[tuple[Optional[int], Optional[int]]]],
    gt_by_pos: dict[tuple[str, int], list[Optional[tuple[int, int]]]],
    n_samples: int,
) -> dict[str, Any]:
    sum_ll = 0.0
    sum_brier = 0.0
    n_scored = 0
    n_carried = 0
    markers_used = 0
    markers_missing_gt = 0
    markers_missing_lai = 0

    for m in panel:
        key = (m["chrom"], m["pos"])
        gts = gt_by_pos.get(key)
        if not gts:
            markers_missing_gt += 1
            continue
        idx = covering_index(lai_positions, m["pos"])
        if idx is None:
            markers_missing_lai += 1
            continue
        exact = lai_positions[idx] == m["pos"]
        site_anc = lai_ancs[idx]
        site_scored = 0
        for s_i in range(n_samples):
            al = gts[s_i] if s_i < len(gts) else None
            if al is None:
                continue
            an1, an2 = site_anc[s_i] if s_i < len(site_anc) else (None, None)
            for allele, anc in ((al[0], an1), (al[1], an2)):
                ll, br, ok = score_haplotype(allele, anc, m["afs"])
                if not ok:
                    continue
                sum_ll += ll
                sum_brier += br
                n_scored += 1
                site_scored += 1
                if not exact:
                    n_carried += 1
        if site_scored:
            markers_used += 1

    mean_ll = sum_ll / n_scored if n_scored else float("nan")
    mean_brier = sum_brier / n_scored if n_scored else float("nan")
    return {
        "mean_ll": mean_ll,
        "mean_brier": mean_brier,
        "n_scored_haps": n_scored,
        "n_markers_used": markers_used,
        "n_markers_panel": len(panel),
        "n_markers_missing_gt": markers_missing_gt,
        "n_markers_missing_lai": markers_missing_lai,
        "n_carried_forward_haps": n_carried,
        "frac_carried_forward": (n_carried / n_scored) if n_scored else float("nan"),
        "n_samples": n_samples,
        "n_lai_markers": len(lai_positions),
    }


def vcf_samples(vcf: str) -> list[str]:
    out = subprocess.run(
        ["bcftools", "query", "-l", vcf],
        check=True,
        capture_output=True,
        text=True,
    )
    return [ln for ln in out.stdout.splitlines() if ln.strip()]


def main(argv: Optional[list[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--anc-vcf", required=True, help="Recipe FLARE anc.vcf.gz")
    p.add_argument("--gt-vcf", required=True, help="Unfiltered phased genotypes")
    p.add_argument("--panel", required=True, type=Path, help="markers.tsv from Part 1")
    p.add_argument("--samples", type=Path, default=None, help="Optional keep-list")
    p.add_argument("--region", default="", help="Optional chr:start-end")
    p.add_argument("--out", required=True, type=Path)
    p.add_argument("--experiment", default="", help="Optional recipe id for JSON")
    args = p.parse_args(argv)

    panel = load_panel(args.panel)
    if not panel:
        raise SystemExit(f"empty panel: {args.panel}")

    keep = load_sample_list(args.samples)
    anc_samps = vcf_samples(args.anc_vcf)
    gt_samps = vcf_samples(args.gt_vcf)
    shared = [s for s in anc_samps if s in set(gt_samps)]
    if keep is not None:
        keep_set = set(keep)
        shared = [s for s in shared if s in keep_set]
    if not shared:
        raise SystemExit("no shared samples between anc, gt, and keep-list")

    lai_pos, lai_ancs = load_anc_sites(args.anc_vcf, shared, region=args.region)
    gt_by_pos = load_gt_alleles_at_panel(
        args.gt_vcf, shared, panel, region=args.region
    )
    summary = score_recipe(
        panel=panel,
        lai_positions=lai_pos,
        lai_ancs=lai_ancs,
        gt_by_pos=gt_by_pos,
        n_samples=len(shared),
    )
    summary["experiment"] = args.experiment
    summary["panel"] = str(args.panel)
    summary["anc_vcf"] = args.anc_vcf
    summary["gt_vcf"] = args.gt_vcf
    summary["region"] = args.region

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as fh:
        json.dump(summary, fh, indent=2, allow_nan=True)
        fh.write("\n")
    print(
        f"mean_ll={summary['mean_ll']:.4f} n_scored={summary['n_scored_haps']} "
        f"markers={summary['n_markers_used']}/{summary['n_markers_panel']}",
        file=sys.stderr,
    )
    print("wrote", args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
