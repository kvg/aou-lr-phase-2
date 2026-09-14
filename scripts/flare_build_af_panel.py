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


def materialize(path_str: str, dest_dir: Path) -> Path:
    """Return a local Path; gsutil-copy gs:// inputs (Path collapses gs:// → gs:/)."""
    raw = str(path_str)
    if raw.startswith("gs:/") and not raw.startswith("gs://"):
        # pathlib / shell mangling
        raw = "gs://" + raw[len("gs:/") :]
    if raw.startswith("gs://"):
        dest_dir.mkdir(parents=True, exist_ok=True)
        dest = dest_dir / Path(raw).name
        if not dest.is_file() or dest.stat().st_size == 0:
            run(["gsutil", "-q", "cp", raw, str(dest)])
        return dest
    p = Path(raw)
    if not p.is_file():
        raise FileNotFoundError(raw)
    return p


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
    """Stream biallelic SNPs via ``bcftools view | query`` (Terra-safe).

    Older bcftools (common on Terra) reject ``query -m2/-M2/-v``; those filters
    belong on ``view``.
    """
    sample_file.write_text("\n".join(samples) + "\n")
    view_cmd = [
        "bcftools",
        "view",
        "-m2",
        "-M2",
        "-v",
        "snps",
        "-S",
        str(sample_file),
        "--force-samples",
        "-Ou",
    ]
    if region.strip():
        view_cmd.extend(["-r", region.strip()])
    view_cmd.append(ref_vcf)
    query_cmd = [
        "bcftools",
        "query",
        "-f",
        "%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n",
    ]
    print("+", " ".join(view_cmd), "|", " ".join(query_cmd), file=sys.stderr)
    view = subprocess.Popen(view_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    proc = subprocess.Popen(
        query_cmd,
        stdin=view.stdout,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1,
    )
    if view.stdout is not None:
        view.stdout.close()  # allow view to receive SIGPIPE if query exits
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
    rc_q = proc.wait()
    view_err = view.stderr.read().decode() if view.stderr is not None else ""
    rc_v = view.wait()
    if rc_v != 0:
        raise SystemExit(f"bcftools view failed ({rc_v}): {view_err.strip()}")
    if rc_q != 0:
        q_err = proc.stderr.read() if proc.stderr is not None else ""
        raise SystemExit(f"bcftools query failed ({rc_q}): {q_err.strip()}")
    return rows


def write_outputs(rows: list[dict], out_prefix: Path) -> None:
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    tsv = Path(str(out_prefix) + ".markers.tsv")
    bed = Path(str(out_prefix) + ".bed")
    cols = ["chrom", "pos", "ref", "alt", "mac", "af_range"] + [f"af_{p}" for p in PANELS]
    # TSV keeps AF-range ranking (panel priority order).
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
    # BED must be chrom/pos sorted for tabix (ranking order is not genomic).
    bed_rows = sorted(rows, key=lambda r: (str(r["chrom"]), int(r["pos"])))
    with bed.open("w") as fh:
        for r in bed_rows:
            fh.write(f"{r['chrom']}\t{r['pos'] - 1}\t{r['pos']}\n")
    if shutil.which("bgzip"):
        try:
            bgzip_tabix_bed(bed)
        except subprocess.CalledProcessError as exc:
            print(f"warning: bgzip/tabix failed ({exc}); markers.tsv is still usable", file=sys.stderr)
    print(f"wrote {tsv} ({len(rows)} markers)", file=sys.stderr)


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--ref-vcf", required=True, help="FLARE reference VCF (gnomAD LAI); local or gs://")
    p.add_argument("--ref-panel", required=True, help="sample\\tpanel refmap; local or gs://")
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

    cache = Path(str(args.out_prefix)).parent / ".cache"
    ref_panel = materialize(str(args.ref_panel), cache)
    ref_vcf = materialize(str(args.ref_vcf), cache)
    # pull VCF index when remote
    raw_vcf = str(args.ref_vcf)
    if raw_vcf.startswith("gs:/") and not raw_vcf.startswith("gs://"):
        raw_vcf = "gs://" + raw_vcf[len("gs:/") :]
    if raw_vcf.startswith("gs://"):
        idx_uri = raw_vcf + (".tbi" if raw_vcf.endswith((".vcf.gz", ".vcf.bgz")) else ".csi")
        try:
            materialize(idx_uri, cache)
        except Exception:
            print(f"warning: no index at {idx_uri}", file=sys.stderr)

    refmap = load_refmap(ref_panel)
    vcf_samples = run(["bcftools", "query", "-l", str(ref_vcf)], capture_output=True).stdout.splitlines()
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
        str(ref_vcf),
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
