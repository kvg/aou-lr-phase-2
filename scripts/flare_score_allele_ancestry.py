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
import os
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


def bcftools_query(args: list[str], *, threads: int = 0):
    """Yield lines from ``bcftools query``.

    ``threads`` is accepted for API compatibility but **ignored**: older Terra
    bcftools builds reject ``query --threads``. Parallelize across recipes
    instead (``ASSOC_JOBS``).
    """
    del threads  # not supported on Terra's bcftools query
    cmd = ["bcftools", "query", *args]
    print("+", " ".join(cmd), file=sys.stderr)
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    assert proc.stdout is not None
    return proc


def _iter_query(proc: subprocess.Popen) -> Iterable[str]:
    assert proc.stdout is not None
    try:
        for line in proc.stdout:
            yield line
    finally:
        if proc.stdout is not None and not proc.stdout.closed:
            proc.stdout.close()
        rc = proc.wait()
        # 141 = SIGPIPE after early close (we got all panel coverings)
        if rc not in (0, 141, -13):
            raise SystemExit(f"bcftools query failed ({rc})")


def _parse_an_row(vals: list[str], n_samp: int) -> Optional[list[tuple[Optional[int], Optional[int]]]]:
    if n_samp and len(vals) != 2 * n_samp:
        return None
    row: list[tuple[Optional[int], Optional[int]]] = []
    for i in range(0, len(vals), 2):
        row.append((parse_an(vals[i]), parse_an(vals[i + 1])))
    return row


def load_anc_sites(
    anc_vcf: str,
    samples: list[str],
    *,
    region: str = "",
    threads: int = 0,
) -> tuple[list[int], list[list[tuple[Optional[int], Optional[int]]]]]:
    """Return (positions, per_site list of (AN1, AN2) per sample index).

    Prefer :func:`load_anc_covering_for_panel` for scoring — this loads *all*
    LAI sites and is much heavier on full chromosomes.
    """
    import tempfile

    tmp: Optional[Path] = None
    sample_arg: list[str] = []
    if samples:
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
    try:
        for line in _iter_query(bcftools_query(cmd, threads=threads)):
            parts = line.rstrip("\n").split("\t")
            if not parts:
                continue
            pos = int(parts[0])
            row = _parse_an_row(parts[1:], n_samp)
            if row is None:
                continue
            positions.append(pos)
            ancs.append(row)
    finally:
        if tmp is not None:
            tmp.unlink(missing_ok=True)
    return positions, ancs


def load_anc_covering_for_panel(
    anc_vcf: str,
    samples: list[str],
    panel_positions: Sequence[int],
    *,
    region: str = "",
    threads: int = 0,
) -> tuple[dict[int, list[tuple[Optional[int], Optional[int]]]], dict[int, int], int]:
    """Stream LAI once; keep covering AN only at panel positions.

    Covering = first LAI POS >= panel POS (else last LAI on contig). Returns
    ``(anc_by_panel_pos, lai_pos_by_panel_pos, n_lai_sites_scanned)``.
    """
    import tempfile

    wanted = sorted({int(p) for p in panel_positions})
    if not wanted:
        return {}, {}, 0

    tmp: Optional[Path] = None
    sample_arg: list[str] = []
    if samples:
        tf = tempfile.NamedTemporaryFile("w", delete=False, suffix=".samples.txt")
        tf.write("\n".join(samples) + "\n")
        tf.close()
        tmp = Path(tf.name)
        sample_arg = ["-S", str(tmp)]
    cmd = [*sample_arg, "-f", r"%POS[\t%AN1\t%AN2]\n"]
    if region.strip():
        cmd.extend(["-r", region.strip()])
    cmd.append(anc_vcf)

    n_samp = len(samples)
    anc_by_pos: dict[int, list[tuple[Optional[int], Optional[int]]]] = {}
    lai_pos_by_panel: dict[int, int] = {}
    pi = 0
    n_lai = 0
    last_row: Optional[list[tuple[Optional[int], Optional[int]]]] = None
    last_pos: Optional[int] = None
    try:
        for line in _iter_query(bcftools_query(cmd, threads=threads)):
            parts = line.rstrip("\n").split("\t")
            if not parts:
                continue
            pos = int(parts[0])
            row = _parse_an_row(parts[1:], n_samp)
            if row is None:
                continue
            n_lai += 1
            last_row, last_pos = row, pos
            while pi < len(wanted) and wanted[pi] <= pos:
                ppos = wanted[pi]
                anc_by_pos[ppos] = row
                lai_pos_by_panel[ppos] = pos
                pi += 1
            if pi >= len(wanted):
                # Still drain? Can break early — remaining LAI not needed.
                break
        # Panel sites past the last LAI site → FELIX "else last"
        if last_row is not None and last_pos is not None and pi < len(wanted):
            while pi < len(wanted):
                ppos = wanted[pi]
                anc_by_pos[ppos] = last_row
                lai_pos_by_panel[ppos] = last_pos
                pi += 1
    finally:
        if tmp is not None:
            tmp.unlink(missing_ok=True)
    return anc_by_pos, lai_pos_by_panel, n_lai


def load_gt_alleles_at_panel(
    gt_vcf: str,
    samples: list[str],
    panel: list[dict[str, Any]],
    *,
    region: str = "",
    threads: int = 0,
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
            # chrom, start, end (1-based inclusive) — portable for bcftools -R
            rf.write(f"{m['chrom']}\t{m['pos']}\t{m['pos']}\n")
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
        for line in _iter_query(bcftools_query(cmd, threads=threads)):
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
    anc_by_panel_pos: dict[int, list[tuple[Optional[int], Optional[int]]]],
    lai_pos_by_panel: dict[int, int],
    gt_by_pos: dict[tuple[str, int], list[Optional[tuple[int, int]]]],
    n_samples: int,
    n_lai_markers: int,
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
        site_anc = anc_by_panel_pos.get(m["pos"])
        if site_anc is None:
            markers_missing_lai += 1
            continue
        lai_pos = lai_pos_by_panel.get(m["pos"])
        exact = lai_pos == m["pos"]
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
        "n_lai_markers": n_lai_markers,
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
    p.add_argument(
        "--threads",
        type=int,
        default=0,
        help="Ignored (Terra bcftools query has no --threads); kept for CLI compat",
    )
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

    anc_by_pos, lai_pos_by_panel, n_lai = load_anc_covering_for_panel(
        args.anc_vcf,
        shared,
        [m["pos"] for m in panel],
        region=args.region,
        threads=args.threads,
    )
    gt_by_pos = load_gt_alleles_at_panel(
        args.gt_vcf,
        shared,
        panel,
        region=args.region,
        threads=args.threads,
    )
    summary = score_recipe(
        panel=panel,
        anc_by_panel_pos=anc_by_pos,
        lai_pos_by_panel=lai_pos_by_panel,
        gt_by_pos=gt_by_pos,
        n_samples=len(shared),
        n_lai_markers=n_lai,
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
