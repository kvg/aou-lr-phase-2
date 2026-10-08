#!/usr/bin/env python3
"""Pedigree Mendelian screen for FLARE local ancestry (association screen, Part 3).

Loads trios from ``aou_phase2.ped`` (id/father/mother), optionally restricts to
complete Phase 2 families via covariates, and scores child AN against parents
at a locus grid (fixed panel or thinned LAI sites).

Single-locus Mendelian rule (unordered parental transmission):
  (c1 in F and c2 in M) or (c1 in M and c2 in F).

Along the window, assignment flips that restore compatibility are counted as
recombinations; hard failures (neither assignment works) are violations.
Expected crossovers come from the PLINK genetic map span × 2 meioses
(father + mother). Gate metric: ``violations_per_informative_locus`` and
``excess_recomb_over_expected``.

Example::

    python3 scripts/flare_score_mendelian_lai.py \\
      --anc-vcf recipe.anc.vcf.gz \\
      --ped tractor_mix/resources/legacy_covariates/aou_phase2.ped \\
      --map maps/plink.chrchr22.GRCh38.map \\
      --panel eval_panels/chr22_10mb.markers.tsv \\
      --region chr22:26897597-36897597 \\
      --out mendelian_lai_score.json
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Optional, Sequence

ANCESTRY = {0: "eas", 1: "amr", 2: "eur", 3: "afr", 4: "sas"}


def vcf_header_text(vcf: str) -> str:
    return subprocess.run(["bcftools", "view", "-h", vcf], check=True,
                          capture_output=True, text=True).stdout


def parse_ancestry_header(header_text: str) -> Optional[dict[int, str]]:
    """``##ANCESTRY=<eas=0,amr=1,...>`` -> {0: "eas", ...}; None when absent.

    Same convention as ``flare_score_allele_ancestry.py``. FLARE2 cluster models
    write ``anc0_afr`` … here, and with ``nanc`` > 5 they emit ancestry codes the
    five-panel ``ANCESTRY`` default does not contain.
    """
    for line in header_text.splitlines():
        if line.startswith("##ANCESTRY=<") and line.rstrip().endswith(">"):
            body = line.strip()[len("##ANCESTRY=<"):-1]
            out: dict[int, str] = {}
            for part in body.split(","):
                name, _, idx = part.partition("=")
                out[int(idx)] = name.strip()
            return out or None
    return None


def load_projection(labels_tsv: Path) -> dict[int, int]:
    """``index -> projected code`` from a ``flare2_build_model.py`` labels.tsv.

    Each ancestry is mapped to its ``dominant_panel``, and panels are coded by
    position in ``ANCESTRY`` so that every recipe lands in the same alphabet.
    Scoring the Mendelian metric on the projection makes recipes with different
    ``nanc`` comparable: the raw violation rate rises with the number of labels
    at fixed painting accuracy, because label compatibility is preserved under
    coarsening but not under refinement (flare/eval/mendelian_k_bias/).
    """
    panel_code = {v: k for k, v in ANCESTRY.items()}
    out: dict[int, int] = {}
    with labels_tsv.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ii, di = header.index("index"), header.index("dominant_panel")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(ii, di):
                continue
            panel = parts[di].strip().lower()
            if panel not in panel_code:
                raise SystemExit(f"{labels_tsv}: unknown dominant_panel {panel!r}")
            out[int(parts[ii])] = panel_code[panel]
    if not out:
        raise SystemExit(f"{labels_tsv}: no ancestry rows")
    return out


@dataclass(frozen=True)
class Trio:
    child: str
    father: str
    mother: str
    family_id: str = ""


def load_ped(path: Path) -> list[Trio]:
    """Parse aou_phase2.ped CSV (ped,id,father,mother,...)."""
    trios: list[Trio] = []
    with path.open() as fh:
        header = fh.readline().rstrip("\n").split(",")
        cols = {c: i for i, c in enumerate(header)}
        for line in fh:
            parts = line.rstrip("\n").split(",")
            if len(parts) < 4:
                continue
            cid = parts[cols["id"]].strip()
            fa = parts[cols["father"]].strip()
            mo = parts[cols["mother"]].strip()
            fam = parts[cols["ped"]].strip() if "ped" in cols else ""
            if fa in {"", "0"} or mo in {"", "0"}:
                continue
            trios.append(Trio(child=cid, father=fa, mother=mo, family_id=fam))
    return trios


def filter_complete_trios(
    trios: list[Trio],
    *,
    covariates: Optional[Path] = None,
    sample_set: Optional[set[str]] = None,
) -> list[Trio]:
    """Keep trios where child+parents are present; optionally complete families.

    When covariates are provided, keep only families where every recorded member
    with that pedigree_family_id is in the Phase 2 PacBio discovery set
    (``long_read_phase == phase_2`` and ``in_pedigree``) — same spirit as
    tableS complete pedigrees. Without covariates, ``sample_set`` alone is used.
    """
    if covariates is not None:
        import pandas as pd

        df = pd.read_csv(covariates, dtype=str, compression="infer")
        if "research_id" not in df.columns:
            raise SystemExit("covariates need research_id")
        df["research_id"] = df["research_id"].astype(str)
        if "long_read_phase" in df.columns:
            df = df[df["long_read_phase"].astype(str).str.lower().eq("phase_2")]
        if "in_pedigree" in df.columns:
            ped_mask = df["in_pedigree"].astype(str).str.lower().isin(
                {"true", "1", "yes"}
            )
            df = df[ped_mask]
        present = set(df["research_id"])
        if "pedigree_family_id" in df.columns and "pedigree_family_size" in df.columns:
            # complete: n people with this family id == pedigree_family_size
            sizes = (
                df.groupby("pedigree_family_id")["research_id"]
                .size()
                .rename("n_obs")
                .reset_index()
            )
            meta = (
                df[["pedigree_family_id", "pedigree_family_size"]]
                .drop_duplicates("pedigree_family_id")
            )
            meta["pedigree_family_size"] = pd.to_numeric(
                meta["pedigree_family_size"], errors="coerce"
            )
            merged = sizes.merge(meta, on="pedigree_family_id")
            complete_fams = set(
                merged.loc[
                    merged["n_obs"] == merged["pedigree_family_size"],
                    "pedigree_family_id",
                ].astype(str)
            )
            trios = [t for t in trios if t.family_id in complete_fams]
        sample_set = present if sample_set is None else (sample_set & present)

    if sample_set is not None:
        trios = [
            t
            for t in trios
            if t.child in sample_set and t.father in sample_set and t.mother in sample_set
        ]
    return trios


def _covering_index(lai_positions: Sequence[int], query_pos: int) -> Optional[int]:
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


def locus_compatible(
    c1: Optional[int],
    c2: Optional[int],
    f1: Optional[int],
    f2: Optional[int],
    m1: Optional[int],
    m2: Optional[int],
) -> bool:
    if None in (c1, c2, f1, f2, m1, m2):
        return False
    F, M = {f1, f2}, {m1, m2}
    return (c1 in F and c2 in M) or (c1 in M and c2 in F)


def assignment_ok(
    c1: int,
    c2: int,
    f_set: set[int],
    m_set: set[int],
    *,
    child_first_from_father: bool,
) -> bool:
    if child_first_from_father:
        return c1 in f_set and c2 in m_set
    return c1 in m_set and c2 in f_set


def score_trio_path(
    loci: Sequence[tuple[Optional[int], Optional[int], Optional[int], Optional[int], Optional[int], Optional[int]]],
) -> dict[str, int]:
    """Walk loci; count informative, hard violations, and assignment flips (recombs)."""
    n_inf = 0
    n_hard = 0
    n_recomb = 0
    assignment: Optional[bool] = None  # True => c1 from father
    for c1, c2, f1, f2, m1, m2 in loci:
        if None in (c1, c2, f1, f2, m1, m2):
            continue
        assert c1 is not None and c2 is not None and f1 is not None and f2 is not None
        assert m1 is not None and m2 is not None
        F, M = {f1, f2}, {m1, m2}
        ok_fm = assignment_ok(c1, c2, F, M, child_first_from_father=True)
        ok_mf = assignment_ok(c1, c2, F, M, child_first_from_father=False)
        if not ok_fm and not ok_mf:
            n_inf += 1
            n_hard += 1
            continue
        n_inf += 1
        if assignment is None:
            assignment = True if ok_fm else False
            continue
        if assignment and ok_fm:
            continue
        if (not assignment) and ok_mf:
            continue
        # Prefer flip if the other assignment works
        if assignment and (not ok_fm) and ok_mf:
            assignment = False
            n_recomb += 1
            continue
        if (not assignment) and (not ok_mf) and ok_fm:
            assignment = True
            n_recomb += 1
            continue
        # Current fails but both somehow inconsistent — treat as hard
        n_hard += 1
    return {
        "n_informative": n_inf,
        "n_hard_violations": n_hard,
        "n_recombinations": n_recomb,
    }


def load_map_cm(path: Path, chrom: str) -> list[tuple[int, float]]:
    """PLINK map: chrom pos(cm?)/ or chrom id pos cm — handle common layouts.

    Supports:
      chrom  rs  cM  bp
      chrom  bp  cM
    """
    pts: list[tuple[int, float]] = []
    with path.open() as fh:
        for line in fh:
            parts = line.strip().split()
            if len(parts) < 3:
                continue
            c = parts[0]
            if c != chrom and c.replace("chr", "") != chrom.replace("chr", ""):
                # also accept chrchr22 style vs chr22
                if c.replace("chr", "") != chrom.replace("chr", ""):
                    continue
            if len(parts) >= 4:
                # chrom id cM bp  OR chrom id bp cM — detect by magnitude
                try:
                    a, b = float(parts[2]), float(parts[3])
                except ValueError:
                    continue
                if a > 1000 and b < 500:  # bp, cM
                    bp, cm = int(a), b
                elif b > 1000 and a < 500:  # cM, bp
                    cm, bp = a, int(b)
                else:
                    # default plink: chrom snp cM bp
                    cm, bp = a, int(b)
                pts.append((bp, cm))
            else:
                try:
                    bp, cm = int(parts[1]), float(parts[2])
                except ValueError:
                    continue
                pts.append((bp, cm))
    pts.sort()
    return pts


def map_cm_at(pts: Sequence[tuple[int, float]], pos: int) -> Optional[float]:
    if not pts:
        return None
    if pos <= pts[0][0]:
        return pts[0][1]
    if pos >= pts[-1][0]:
        return pts[-1][1]
    lo, hi = 0, len(pts) - 1
    while lo + 1 < hi:
        mid = (lo + hi) // 2
        if pts[mid][0] <= pos:
            lo = mid
        else:
            hi = mid
    p0, c0 = pts[lo]
    p1, c1 = pts[hi]
    if p1 == p0:
        return c0
    t = (pos - p0) / (p1 - p0)
    return c0 + t * (c1 - c0)


def expected_crossovers(
    map_pts: Sequence[tuple[int, float]],
    start: int,
    end: int,
    *,
    n_meioses: int = 2,
) -> float:
    """Expected number of crossovers across n_meioses in [start,end]."""
    c0 = map_cm_at(map_pts, start)
    c1 = map_cm_at(map_pts, end)
    if c0 is None or c1 is None:
        return float("nan")
    morgan = abs(c1 - c0) / 100.0
    return n_meioses * morgan


WHOLE_CONTIG_END = 10**10


def parse_region(spec: str) -> tuple[str, int, int]:
    """``chr:start-end``, or a bare contig (``chr20``) for the whole chromosome."""
    spec = spec.strip()
    if not spec:
        raise ValueError("empty region")
    if ":" not in spec:
        return spec, 1, WHOLE_CONTIG_END
    chrom, rest = spec.split(":", 1)
    start_s, end_s = rest.split("-", 1)
    return chrom, int(start_s), int(end_s)


def bcftools_samples(vcf: str) -> list[str]:
    out = subprocess.run(
        ["bcftools", "query", "-l", vcf], check=True, capture_output=True, text=True
    )
    return [ln for ln in out.stdout.splitlines() if ln.strip()]


def load_trio_anc_grid(
    anc_vcf: str,
    samples: list[str],
    positions: Sequence[int],
    *,
    region: str,
    allowed: Optional[set[int]] = None,
    projection: Optional[dict[int, int]] = None,
    stats: Optional[dict[str, int]] = None,
) -> dict[int, dict[str, tuple[Optional[int], Optional[int]]]]:
    """pos -> sample -> (AN1, AN2) for requested samples at (approx) grid positions.

    Queries all LAI sites in region then subsets to covering indices for each
    requested panel/grid position (FELIXla covering).
    """
    with tempfile.NamedTemporaryFile("w", delete=False, suffix=".samples.txt") as sf:
        sf.write("\n".join(samples) + "\n")
        s_path = sf.name
    cmd = [
        "bcftools",
        "query",
        "-S",
        s_path,
        "-f",
        r"%POS[\t%AN1\t%AN2]\n",
    ]
    if region.strip():
        cmd.extend(["-r", region.strip()])
    cmd.append(anc_vcf)
    print("+", " ".join(cmd), file=sys.stderr)
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
    assert proc.stdout is not None
    lai_pos: list[int] = []
    lai_rows: list[list[tuple[Optional[int], Optional[int]]]] = []
    for line in proc.stdout:
        parts = line.rstrip("\n").split("\t")
        pos = int(parts[0])
        vals = parts[1:]
        row: list[tuple[Optional[int], Optional[int]]] = []
        for i in range(0, len(vals), 2):
            def _p(x: str) -> Optional[int]:
                if x in {".", ""}:
                    if stats is not None:
                        stats["n_missing_calls"] = stats.get("n_missing_calls", 0) + 1
                    return None
                try:
                    v = int(x)
                except ValueError:
                    if stats is not None:
                        stats["n_unparsable_calls"] = stats.get("n_unparsable_calls", 0) + 1
                    return None
                if allowed is not None and v not in allowed:
                    # Previously this silently returned None, so every locus
                    # touching an out-of-range code was dropped from the
                    # denominator. Count it; main() refuses to report a score.
                    if stats is not None:
                        stats["n_out_of_range_calls"] = stats.get("n_out_of_range_calls", 0) + 1
                        stats["max_code_seen"] = max(stats.get("max_code_seen", -1), v)
                    return None
                return projection.get(v, v) if projection else v

            row.append((_p(vals[i]), _p(vals[i + 1])))
        lai_pos.append(pos)
        lai_rows.append(row)
    rc = proc.wait()
    Path(s_path).unlink(missing_ok=True)
    if rc != 0:
        raise SystemExit(f"bcftools query failed ({rc})")

    out: dict[int, dict[str, tuple[Optional[int], Optional[int]]]] = {}
    for q in positions:
        idx = _covering_index(lai_pos, q)
        if idx is None:
            continue
        row = lai_rows[idx]
        out[q] = {samples[i]: row[i] for i in range(len(samples))}
    return out


def load_panel_positions(panel: Path) -> list[tuple[str, int]]:
    rows: list[tuple[str, int]] = []
    with panel.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ci, pi = header.index("chrom"), header.index("pos")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            rows.append((parts[ci], int(parts[pi])))
    return rows


def main(argv: Optional[list[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--anc-vcf", required=True)
    p.add_argument("--ped", required=True, type=Path)
    p.add_argument("--map", required=True, type=Path, help="PLINK genetic map")
    p.add_argument("--panel", type=Path, default=None, help="Optional fixed marker panel")
    p.add_argument("--region", required=True, help="chr:start-end or bare contig (whole chromosome)")
    p.add_argument("--covariates", type=Path, default=None)
    p.add_argument("--samples", type=Path, default=None, help="Optional analysis keep-list")
    p.add_argument("--thin-every", type=int, default=1, help="Keep every Nth panel/grid locus")
    p.add_argument("--out", required=True, type=Path)
    p.add_argument("--experiment", default="")
    p.add_argument(
        "--num-ancs", type=int, default=0,
        help="Number of ancestry codes (0..N-1). Default: read the ##ANCESTRY "
             "header, else the five reference panels",
    )
    p.add_argument(
        "--project-labels", type=Path, default=None,
        help="labels.tsv from flare2_build_model.py; score on each ancestry's "
             "dominant_panel so recipes with different nanc are comparable",
    )
    p.add_argument(
        "--allow-out-of-range", action="store_true",
        help="Score even when ancestry codes fall outside the alphabet (they "
             "are treated as missing, which biases the denominator)",
    )
    args = p.parse_args(argv)

    chrom, start, end = parse_region(args.region)
    trios = load_ped(args.ped)
    sample_set: Optional[set[str]] = None
    if args.samples:
        sample_set = {
            ln.strip() for ln in args.samples.read_text().splitlines() if ln.strip()
        }
    vcf_samps = set(bcftools_samples(args.anc_vcf))
    sample_set = vcf_samps if sample_set is None else (sample_set & vcf_samps)
    trios = filter_complete_trios(
        trios, covariates=args.covariates, sample_set=sample_set
    )
    if not trios:
        raise SystemExit("no scorable trios after filters")

    needed = sorted({s for t in trios for s in (t.child, t.father, t.mother)})
    if args.panel:
        pos_rows = [(c, pos) for c, pos in load_panel_positions(args.panel) if c == chrom]
        positions = [pos for _, pos in pos_rows if start <= pos <= end]
    else:
        # thin LAI sites via bcftools
        cmd = ["bcftools", "query", "-f", "%POS\n", "-r", args.region, args.anc_vcf]
        positions = [
            int(x)
            for x in subprocess.run(cmd, check=True, capture_output=True, text=True).stdout.split()
        ]
    if args.thin_every > 1:
        positions = positions[:: args.thin_every]
    if not positions:
        raise SystemExit("no grid positions in region")

    # Ancestry alphabet comes from the VCF's own ##ANCESTRY header when present
    # (FLARE2 cluster models with nanc > 5 emit codes the five-panel default
    # does not contain). Without this, every locus touching such a code was
    # read as missing and silently left out of the denominator.
    anc_names = parse_ancestry_header(vcf_header_text(args.anc_vcf))
    if args.num_ancs:
        allowed = set(range(args.num_ancs))
        source = f"--num-ancs={args.num_ancs}"
    elif anc_names:
        allowed = set(anc_names)
        source = "##ANCESTRY header"
    else:
        allowed = set(ANCESTRY)
        source = "five-panel default"
    projection = load_projection(args.project_labels) if args.project_labels else None
    if projection is not None:
        missing = sorted(allowed - set(projection))
        if missing:
            raise SystemExit(
                f"--project-labels does not cover ancestry codes {missing}"
            )
    grid_stats: dict[str, int] = {}
    print(f"ancestry alphabet: {sorted(allowed)} (from {source})"
          + (f"; projected to {sorted(set(projection.values()))}" if projection else ""),
          file=sys.stderr)
    grid = load_trio_anc_grid(
        args.anc_vcf, needed, positions, region=args.region,
        allowed=allowed, projection=projection, stats=grid_stats,
    )
    n_oor = grid_stats.get("n_out_of_range_calls", 0)
    if n_oor and not args.allow_out_of_range:
        raise SystemExit(
            f"{n_oor} ancestry calls outside {sorted(allowed)} "
            f"(max code seen {grid_stats.get('max_code_seen')}). These would be "
            "dropped as missing, shrinking the denominator and making the "
            "violation rate incomparable across recipes. Pass --num-ancs for "
            "this model, or --allow-out-of-range to score anyway."
        )
    map_pts = load_map_cm(args.map, chrom)
    # try chrchrN naming if empty
    if not map_pts and not chrom.startswith("chrchr"):
        map_pts = load_map_cm(args.map, "chr" + chrom if not chrom.startswith("chr") else chrom)
    # Span of the scored grid (not the requested region), so a bare-contig
    # region does not count map distance outside the LAI markers.
    exp_xo = expected_crossovers(map_pts, positions[0], positions[-1], n_meioses=2)

    tot_inf = 0
    tot_hard = 0
    tot_recomb = 0
    per_trio: list[dict[str, Any]] = []
    for t in trios:
        loci = []
        for pos in positions:
            site = grid.get(pos)
            if not site:
                loci.append((None, None, None, None, None, None))
                continue
            c = site.get(t.child, (None, None))
            f = site.get(t.father, (None, None))
            m = site.get(t.mother, (None, None))
            loci.append((c[0], c[1], f[0], f[1], m[0], m[1]))
        sc = score_trio_path(loci)
        tot_inf += sc["n_informative"]
        tot_hard += sc["n_hard_violations"]
        tot_recomb += sc["n_recombinations"]
        per_trio.append(
            {
                "child": t.child,
                "father": t.father,
                "mother": t.mother,
                "family_id": t.family_id,
                **sc,
                "excess_recomb": sc["n_recombinations"]
                - (0.0 if math.isnan(exp_xo) else exp_xo),
            }
        )

    n_trios = len(trios)
    summary = {
        "experiment": args.experiment,
        "region": args.region,
        # The violation rate is only comparable across recipes scored on the
        # same alphabet; record it so a ranking cannot silently mix resolutions.
        "ancestry_alphabet": sorted(allowed),
        "ancestry_alphabet_source": source,
        "n_ancestries": len(allowed),
        "projected_to": sorted(set(projection.values())) if projection else None,
        "grid_call_stats": grid_stats,
        "n_trios_scored": n_trios,
        "n_loci_grid": len(positions),
        "n_informative_locus_calls": tot_inf,
        "n_hard_violations": tot_hard,
        "violations_per_informative_locus": (
            tot_hard / tot_inf if tot_inf else float("nan")
        ),
        "n_recombinations": tot_recomb,
        "mean_recombinations_per_trio": tot_recomb / n_trios if n_trios else float("nan"),
        "expected_crossovers_per_trio": exp_xo,
        "excess_recomb_over_expected": (
            (tot_recomb / n_trios - exp_xo)
            if n_trios and not math.isnan(exp_xo)
            else float("nan")
        ),
        "ped": str(args.ped),
        "map": str(args.map),
        "panel": str(args.panel) if args.panel else "",
        "trios": per_trio,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as fh:
        json.dump(summary, fh, indent=2, allow_nan=True)
        fh.write("\n")
    print(
        f"trios={n_trios} viol_rate={summary['violations_per_informative_locus']:.4g} "
        f"excess_xo={summary['excess_recomb_over_expected']:.3g}",
        file=sys.stderr,
    )
    print("wrote", args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
