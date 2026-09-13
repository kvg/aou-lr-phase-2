#!/usr/bin/env python3
"""Build a fixed AF-divergent marker panel for LAI recipe evaluation (Part 1).

Computes ancestry-specific allele frequencies from the FLARE reference VCF +
refmap, ranks biallelic SNVs by cross-panel AF divergence, and writes a fixed
marker list used by every recipe's allele–ancestry concordance score.

Writes:
  * ``{out_prefix}.markers.tsv`` — chrom pos ref alt + af_<panel> + af_range
  * ``{out_prefix}.bed.gz`` (+ ``.tbi`` when bgzip/tabix available)

Example::

    python3 scripts/flare_build_af_panel.py \\
      --ref-vcf chr22.gnomad_lai_90.vcf.bgz \\
      --ref-panel aou_1000genomes.refmap \\
      --region chr22:26897597-36897597 \\
      --top-n 5000 \\
      --out-prefix eval_panels/chr22_10mb
"""

from __future__ import annotations

import argparse
import math
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

PANELS = ("eas", "amr", "eur", "afr", "sas")


def run(cmd: list[str], **kwargs) -> subprocess.CompletedProcess[str]:
    print("+", " ".join(map(str, cmd)), file=sys.stderr)
    return subprocess.run(cmd, check=True, text=True, **kwargs)


def load_refmap(path: Path) -> dict[str, str]:
    """sample -> panel (lowercase)."""
    out: dict[str, str] = {}
    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            sample, panel = parts[0], parts[1].strip().lower()
            out[sample] = panel
    return out


def af_range(afs: dict[str, float]) -> float:
    vals = [afs[p] for p in PANELS if p in afs and not math.isnan(afs[p])]
    if len(vals) < 2:
        return 0.0
    return max(vals) - min(vals)


def site_af_from_genotypes(
    gts: list[str],
    sample_panels: list[str],
) -> tuple[dict[str, float], int]:
    """Return (af_by_panel, total_allele_count_alt_plus_ref_called)."""
    alt = defaultdict(int)
    tot = defaultdict(int)
    for gt, panel in zip(gts, sample_panels):
        if panel not in PANELS:
            continue
        g = gt.replace("|", "/").split("/")
        if len(g) != 2 or g[0] in {".", ""} or g[1] in {".", ""}:
            continue
        try:
            a0, a1 = int(g[0]), int(g[1])
        except ValueError:
            continue
        tot[panel] += 2
        alt[panel] += int(a0 == 1) + int(a1 == 1)
    afs: dict[str, float] = {}
    mac_total = 0
    alleles_total = 0
    for p in PANELS:
        if tot[p] == 0:
            afs[p] = float("nan")
        else:
            afs[p] = alt[p] / tot[p]
            mac_total += min(alt[p], tot[p] - alt[p])
            alleles_total += tot[p]
    return afs, alleles_total


def select_markers(
    rows: list[dict],
    *,
    top_n: int,
    min_mac: int,
    min_af_range: float,
) -> list[dict]:
    """Filter and rank marker dicts that already have af_* and af_range."""
    kept = []
    for r in rows:
        if int(r.get("mac", 0)) < min_mac:
            continue
        if float(r.get("af_range", 0.0)) < min_af_range:
            continue
        # require at least 2 panels with finite AF
        n_fin = sum(1 for p in PANELS if not math.isnan(float(r.get(f"af_{p}", "nan"))))
        if n_fin < 2:
            continue
        kept.append(r)
    kept.sort(key=lambda r: (-float(r["af_range"]), int(r["pos"])))
    return kept[:top_n]


def bgzip_tabix_bed(bed: Path) -> Path:
    gz = Path(str(bed) + ".gz") if not str(bed).endswith(".gz") else bed
    if not str(bed).endswith(".gz"):
        run(["bgzip", "-f", str(bed)])
        bed = gz
    if shutil.which("tabix"):
        run(["tabix", "-f", "-p", "bed", str(bed)])
    return bed


def stream_ref_sites(
    ref_vcf: str,
    *,
    region: str,
    samples: list[str],
    sample_panels: list[str],
    sample_file: Path,
) -> list[dict]:
    sample_file.write_text("\n".join(samples) + "\n")
    cmd = [
        "bcftools",
        "query",
        "-S",
        str(sample_file),
        "-m2",
        "-M2",
        "-v",
        "snps",
        "-f",
        "%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n",
    ]
    if region.strip():
        cmd.extend(["-r", region.strip()])
    cmd.append(ref_vcf)
    print("+", " ".join(cmd), file=sys.stderr)
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True, bufsize=1)
    assert proc.stdout is not None
    rows: list[dict] = []
    n = 0
    for line in proc.stdout:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 4 + len(samples):
            continue
        chrom, pos_s, ref, alt = parts[0], parts[1], parts[2], parts[3]
        if len(ref) != 1 or len(alt) != 1:
            continue
        gts = parts[4 : 4 + len(samples)]
        afs, _tot_all = site_af_from_genotypes(gts, sample_panels)
        alt_c = 0
        tot_c = 0
        for gt, panel in zip(gts, sample_panels):
            if panel not in PANELS:
                continue
            g = gt.replace("|", "/").split("/")
            if len(g) != 2 or "." in g:
                continue
            try:
                a0, a1 = int(g[0]), int(g[1])
            except ValueError:
                continue
            tot_c += 2
            alt_c += (a0 == 1) + (a1 == 1)
        mac = min(alt_c, tot_c - alt_c) if tot_c else 0
        row = {
            "chrom": chrom,
            "pos": int(pos_s),
            "ref": ref,
            "alt": alt,
            "mac": mac,
            "n_alleles": tot_c,
            "af_range": af_range(afs),
        }
        for p in PANELS:
            row[f"af_{p}"] = afs.get(p, float("nan"))
        rows.append(row)
        n += 1
        if n % 100000 == 0:
            print(f"... scanned {n} sites", file=sys.stderr)
    rc = proc.wait()
    if rc != 0:
        raise SystemExit(f"bcftools query failed ({rc})")
    return rows


def write_outputs(rows: list[dict], out_prefix: Path) -> None:
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    tsv = Path(str(out_prefix) + ".markers.tsv")
    bed = Path(str(out_prefix) + ".bed")
    cols = ["chrom", "pos", "ref", "alt", "mac", "af_range"] + [f"af_{p}" for p in PANELS]
    with tsv.open("w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            vals = []
            for c in cols:
                v = r[c]
                if isinstance(v, float):
                    vals.append("" if math.isnan(v) else f"{v:.6g}")
                else:
                    vals.append(str(v))
            fh.write("\t".join(vals) + "\n")
    with bed.open("w") as fh:
        for r in rows:
            fh.write(f"{r['chrom']}\t{r['pos'] - 1}\t{r['pos']}\n")
    if shutil.which("bgzip"):
        bgzip_tabix_bed(bed)
    print(f"wrote {tsv} ({len(rows)} markers)", file=sys.stderr)


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--ref-vcf", required=True, help="FLARE reference VCF (gnomAD LAI)")
    p.add_argument("--ref-panel", required=True, type=Path, help="sample\\tpanel refmap")
    p.add_argument("--region", default="", help="chr or chr:start-end")
    p.add_argument("--top-n", type=int, default=5000)
    p.add_argument("--min-mac", type=int, default=50, help="minor allele count across ref")
    p.add_argument(
        "--min-af-range",
        type=float,
        default=0.05,
        help="minimum max_panel_AF - min_panel_AF",
    )
    p.add_argument("--out-prefix", required=True)
    args = p.parse_args(argv)

    if not shutil.which("bcftools"):
        raise SystemExit("bcftools is required")

    refmap = load_refmap(args.ref_panel)
    vcf_samples = run(["bcftools", "query", "-l", args.ref_vcf], capture_output=True).stdout.splitlines()
    samples = [s for s in vcf_samples if s in refmap and refmap[s] in PANELS]
    if not samples:
        raise SystemExit("no overlapping samples between VCF and refmap panels")
    sample_panels = [refmap[s] for s in samples]
    print(
        f"ref samples used: {len(samples)} "
        f"({ {p: sample_panels.count(p) for p in PANELS} })",
        file=sys.stderr,
    )

    sample_file = Path(str(args.out_prefix) + ".samples.txt")
    Path(args.out_prefix).parent.mkdir(parents=True, exist_ok=True)
    all_rows = stream_ref_sites(
        args.ref_vcf,
        region=args.region,
        samples=samples,
        sample_panels=sample_panels,
        sample_file=sample_file,
    )
    selected = select_markers(
        all_rows,
        top_n=args.top_n,
        min_mac=args.min_mac,
        min_af_range=args.min_af_range,
    )
    write_outputs(selected, Path(args.out_prefix))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
