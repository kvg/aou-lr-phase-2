#!/usr/bin/env python3
"""Pull per-sample FORMAT fields at a handful of named sites out of any indexed VCF, without localising it.

Written for the hit dossier's call-rate question, but it is deliberately general: give it a site list and a VCF URI and
it returns one row per site and sample. The dossier needs 90 lead positions out of callsets that are hundreds of GB, so
it range-reads through the index -- one ``bcftools query -R`` per VCF, one index seek per site -- rather than streaming
or merging anything.

    python3 pull_site_fields.py --sites leads.tsv --vcf gs://.../chr9.vcf.gz --out chr9.fields.tsv.gz

``--sites`` is a TSV with ``locus_id, chrom, pos`` and optionally ``ref, alt``. Supply ref and alt: several FELIX
markers can share one POS, and a POS-only match has already picked the wrong marker once in this project.

On Terra, htslib range-reads ``gs://`` directly if the environment carries credentials::

    export GCS_OAUTH_TOKEN=$(gcloud auth print-access-token)
    export GCS_REQUESTER_PAYS_PROJECT=$GOOGLE_PROJECT

Why this and not a propagate/annotate workflow: the information wanted is a few hundred numbers, the targets are known
in advance, and every indexed callset -- pre-phasing DeepVariant/GLnexus, the SV callset, the phased ``gt_vcf``, the
FLARE anc VCF -- answers the same query with the same code.
"""

from __future__ import annotations

import argparse
import gzip
import math
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

DEFAULT_FIELDS = ("GT", "DP", "GQ", "RNC")
NOCALL = re.compile(r"^[.](?:[/|][.])*$")  # '.', './.', '.|.'
SITE_COLS = ("locus_id", "chrom", "pos")


def _chrom_key(v: str) -> str:
    return re.sub(r"^chr", "", str(v).strip(), flags=re.IGNORECASE)


def read_sites(path) -> pd.DataFrame:
    """The site list: locus_id, chrom, pos, and ref/alt when given."""
    df = pd.read_csv(path, sep="\t")
    missing = [c for c in SITE_COLS if c not in df.columns]
    if missing:
        raise KeyError(f"{path}: site list lacks column(s) {missing}; has {list(df.columns)}")
    df = df.copy()
    df["pos"] = df["pos"].astype(int)
    for c in ("ref", "alt"):
        if c not in df.columns:
            df[c] = None
    return df


def probe(vcf: str, gcs_project: str | None = None) -> str:
    """One-line verdict on whether bcftools can open ``vcf`` at all, with htslib's own message when it cannot.

    Worth running once before a loop over chromosomes: a gs:// URI needs an htslib built with libcurl plus
    GCS_OAUTH_TOKEN in the environment, and when either is missing every query fails the same way.
    """
    out = subprocess.run(_bcftools(["view", "-h", vcf], gcs_project), capture_output=True, text=True)
    if out.returncode == 0:
        n = len([l for l in out.stdout.splitlines() if l.startswith("##")])
        return f"OK: header read ({n} header lines)"
    env = [k for k in ("GCS_OAUTH_TOKEN", "GCS_REQUESTER_PAYS_PROJECT") if os.environ.get(k)]
    return (f"FAILED (exit {out.returncode}). bcftools said:\n    "
            + "\n    ".join(out.stderr.strip().splitlines()[-6:] or ["<no stderr>"])
            + f"\n  environment carrying credentials: {env or 'none of GCS_OAUTH_TOKEN, GCS_REQUESTER_PAYS_PROJECT'}")


def available_format_tags(vcf: str, gcs_project: str | None = None) -> set[str]:
    """FORMAT tags the VCF header declares. Asking bcftools for an undeclared tag is an error, so this gates the query."""
    out = subprocess.run(_bcftools(["view", "-h", vcf], gcs_project), capture_output=True, text=True)
    if out.returncode != 0:
        raise RuntimeError(f"bcftools view -h failed for {vcf}\n  {out.stderr.strip()[-700:] or '<no stderr>'}")
    return set(re.findall(r"^##FORMAT=<ID=([^,>]+)", out.stdout, flags=re.MULTILINE))


def samples(vcf: str, gcs_project: str | None = None) -> list[str]:
    """Sample names in file order, which is the order the per-sample fields come back in."""
    out = subprocess.run(_bcftools(["query", "-l", vcf], gcs_project), capture_output=True, text=True)
    if out.returncode != 0:
        raise RuntimeError(f"bcftools query -l failed for {vcf}:\n{out.stderr.strip()[:500]}")
    return out.stdout.split()


def _bcftools(args: list[str], gcs_project: str | None) -> list[str]:
    exe = os.environ.get("BCFTOOLS", "bcftools")
    return [exe] + args


def write_regions(sites: pd.DataFrame, path: Path, chrom_prefix: str = "") -> int:
    """A bcftools -R regions file, one line per site. Returns the number written."""
    with open(path, "w") as fh:
        for _, s in sites.sort_values(["chrom", "pos"]).iterrows():
            fh.write(f"{chrom_prefix}{s['chrom']}\t{int(s['pos'])}\t{int(s['pos'])}\n")
    return len(sites)


def query_sites(vcf: str, sites: pd.DataFrame, fields: tuple[str, ...] = DEFAULT_FIELDS,
                gcs_project: str | None = None, chrom_prefix: str = "chr") -> pd.DataFrame:
    """One row per (site, sample) with the requested FORMAT fields, from one indexed range-read of ``vcf``.

    Only tags the header declares are requested; the rest come back as a column of NaN so the caller can see which were
    unavailable rather than having the query fail. Rows are matched to the site list on (chrom, pos) and, when the site
    list gives ref/alt, on those too -- several records can share a position.
    """
    have = available_format_tags(vcf, gcs_project)
    want = [f for f in fields if f in have]
    absent = [f for f in fields if f not in have]
    names = samples(vcf, gcs_project)
    if not names:
        raise RuntimeError(f"{vcf} has no samples")

    with tempfile.TemporaryDirectory() as tmp:
        reg = Path(tmp) / "regions.tsv"
        write_regions(sites, reg, chrom_prefix=chrom_prefix)
        per_sample = "".join("%" + f for f in want)
        fmt = "%CHROM\\t%POS\\t%REF\\t%ALT[\\t" + ":".join("%" + f for f in want) + "]\\n" if want else \
              "%CHROM\\t%POS\\t%REF\\t%ALT\\n"
        cmd = _bcftools(["query", "-R", str(reg), "-f", fmt, vcf], gcs_project)
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(f"bcftools query failed for {vcf}:\n{proc.stderr.strip()[:800]}")
        text = proc.stdout
    del per_sample

    rows = []
    for line in text.splitlines():
        if not line:
            continue
        f = line.split("\t")
        chrom, pos, ref, alt = _chrom_key(f[0]), int(f[1]), f[2], f[3]
        cells = f[4:]
        if len(cells) != len(names):
            raise RuntimeError(f"{vcf}: {len(cells)} sample columns for {len(names)} samples at {chrom}:{pos}")
        for name, cell in zip(names, cells):
            vals = cell.split(":")
            rec = dict(chrom=chrom, pos=pos, ref=ref, alt=alt, sample=name)
            for i, tag in enumerate(want):
                rec[tag] = vals[i] if i < len(vals) else "."
            rows.append(rec)
    out = pd.DataFrame(rows)
    for tag in absent:
        out[tag] = math.nan
    if out.empty:
        return out
    return _match_to_sites(out, sites)


def _match_to_sites(rows: pd.DataFrame, sites: pd.DataFrame) -> pd.DataFrame:
    """Attach locus_id, matching on (chrom, pos) and additionally on ref/alt where the site list gives them."""
    s = sites.copy()
    s["chrom"] = s["chrom"].map(_chrom_key)
    keyed = s[s.ref.notna() & s.alt.notna()]
    loose = s[~(s.ref.notna() & s.alt.notna())]
    out = []
    if len(keyed):
        out.append(rows.merge(keyed[["locus_id", "chrom", "pos", "ref", "alt"]],
                              on=["chrom", "pos", "ref", "alt"], how="inner"))
    if len(loose):
        out.append(rows.merge(loose[["locus_id", "chrom", "pos"]], on=["chrom", "pos"], how="inner"))
    return pd.concat(out, ignore_index=True) if out else rows.iloc[0:0]


def is_nocall(gt) -> bool:
    return bool(NOCALL.match(str(gt).strip())) if gt is not None else True


def call_rate_by_group(rows: pd.DataFrame, groups: dict[str, str], dp_gq: bool = True) -> pd.DataFrame:
    """Per locus, the call rate in cases and in controls, with DP and GQ medians.

    ``groups`` maps sample name to 'case' or 'control'; samples absent from it are ignored, which is how the analysis
    cohort and the phecode's own eligible set are applied. Returns one row per locus_id with the two call rates, their
    difference and a two-sided Fisher p on the 2x2 of called against no-call.
    """
    if rows.empty:
        return pd.DataFrame()
    r = rows.copy()
    r["group"] = r["sample"].map(groups)
    r = r[r.group.isin(("case", "control"))]
    if r.empty:
        return pd.DataFrame()
    r["nocall"] = r["GT"].map(is_nocall) if "GT" in r.columns else True
    out = []
    for lid, g in r.groupby("locus_id"):
        rec = dict(locus_id=lid)
        for who in ("case", "control"):
            sub = g[g.group == who]
            n = len(sub)
            miss = int(sub.nocall.sum())
            rec[f"n_{who}"] = n
            rec[f"n_{who}_nocall"] = miss
            rec[f"call_rate_{who}"] = (n - miss) / n if n else math.nan
            if dp_gq:
                for tag in ("DP", "GQ"):
                    v = pd.to_numeric(sub.get(tag), errors="coerce") if tag in sub.columns else None
                    rec[f"median_{tag}_{who}"] = float(v.median()) if v is not None and v.notna().any() else math.nan
        rec["call_rate_diff"] = rec["call_rate_case"] - rec["call_rate_control"]
        rec["fisher_p"] = _fisher(rec["n_case"] - rec["n_case_nocall"], rec["n_case_nocall"],
                                  rec["n_control"] - rec["n_control_nocall"], rec["n_control_nocall"])
        out.append(rec)
    return pd.DataFrame(out)


def imputed_fraction_by_group(pre: pd.DataFrame, post: pd.DataFrame | None, groups: dict[str, str]) -> pd.DataFrame:
    """Per locus, the share of the genotypes the scan USED that phasing filled in rather than observed.

    ``pre`` is the pre-phasing callset's rows. A genotype counted as imputed is a no-call in ``pre`` that is called
    afterwards. Split by case and control, because a difference between the two is the failure mode a plain call-rate
    comparison cannot see: if the scan leaned on filled-in genotypes more for cases than for controls, the association
    is partly about phasing.

    ``post`` may be None, which asserts that phasing filled every no-call so every pre-phasing no-call is an imputed
    genotype. That is the measured situation here -- MissingRate is exactly 0 at every record of FELIX's input -- and
    it makes the phased query unnecessary. Pass the phased rows instead to verify it rather than assume it.
    """
    if pre.empty:
        return pd.DataFrame()
    a = pre[["locus_id", "sample", "GT"]].rename(columns={"GT": "gt_pre"})
    if post is None or len(post) == 0:
        m = a.copy()
        m["gt_post"] = "0/0"  # called, whatever the call is: phasing leaves no no-calls behind
    else:
        b = post[["locus_id", "sample", "GT"]].rename(columns={"GT": "gt_post"})
        m = a.merge(b, on=["locus_id", "sample"], how="inner")
    m["group"] = m["sample"].map(groups)
    m = m[m.group.isin(("case", "control"))]
    if m.empty:
        return pd.DataFrame()
    m["imputed"] = m.gt_pre.map(is_nocall) & ~m.gt_post.map(is_nocall)
    out = []
    for lid, g in m.groupby("locus_id"):
        rec = dict(locus_id=lid, n_compared=len(g))
        for who in ("case", "control"):
            sub = g[g.group == who]
            rec[f"n_{who}_compared"] = len(sub)
            rec[f"imputed_frac_{who}"] = float(sub.imputed.mean()) if len(sub) else math.nan
        rec["imputed_frac_diff"] = rec["imputed_frac_case"] - rec["imputed_frac_control"]
        out.append(rec)
    return pd.DataFrame(out)


def _fisher(a: int, b: int, c: int, d: int) -> float:
    """Two-sided Fisher exact p on [[a, b], [c, d]]; scipy when available, else nan (the counts still stand)."""
    try:
        from scipy.stats import fisher_exact

        return float(fisher_exact([[a, b], [c, d]])[1])
    except ImportError:
        return math.nan


def phenotype_column(df: pd.DataFrame, names) -> str:
    """The first of ``names`` that is a column of ``df``; raises KeyError naming all of them if none is.

    pheno_cov.tsv names its phenotype columns with the pipeline's underscore form (``BI_164_6``), while the PhecodeX
    tables and the figure code use the dotted form (``BI_164.6``). They coincide only for codes without a sub-code, so
    looking up the dotted form alone silently finds 82 of 193 dossier phenotypes and reports the rest as absent.
    """
    if isinstance(names, str):
        names = [names]
    for n in names:
        if n in df.columns:
            return n
    raise KeyError(f"pheno_cov has none of the column name(s) {list(names)}")


def case_control_map(pheno_cov, phecode, analysis_samples=None, id_col: str = "person_id") -> dict[str, str]:
    """{sample: 'case'|'control'} for one phenotype from pheno_cov.tsv, restricted to the analysis sample list.

    ``phecode`` may be one name or several candidates, tried in order -- pass both the underscore and the dotted form
    and whichever pheno_cov uses is found (see ``phenotype_column``).

    The phenotype columns are 0/1; anything else (blank, NA) is neither a case nor a control and is left out, which is
    how a phenotype's own eligible set differs from the cohort.
    """
    df = pheno_cov if isinstance(pheno_cov, pd.DataFrame) else pd.read_csv(pheno_cov, sep="\t", low_memory=False)
    if id_col not in df.columns:
        raise KeyError(f"pheno_cov lacks {id_col}; has {list(df.columns)[:10]}...")
    phecode = phenotype_column(df, phecode)
    keep = None
    if analysis_samples is not None:
        if isinstance(analysis_samples, (str, Path)):
            keep = {l.split()[0] for l in Path(analysis_samples).read_text().split("\n") if l.strip()}
        else:
            keep = set(analysis_samples)
    v = pd.to_numeric(df[phecode], errors="coerce")
    out = {}
    for sid, val in zip(df[id_col].astype(str), v):
        if keep is not None and sid not in keep:
            continue
        if val == 1:
            out[sid] = "case"
        elif val == 0:
            out[sid] = "control"
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sites", required=True, help="TSV with locus_id, chrom, pos and optionally ref, alt")
    ap.add_argument("--vcf", required=True, action="append",
                    help="indexed VCF/BCF (local or gs://); repeat for several, results are concatenated with a label")
    ap.add_argument("--label", action="append", default=None, help="label per --vcf (default: the file name)")
    ap.add_argument("--fields", default=",".join(DEFAULT_FIELDS))
    ap.add_argument("--chrom-prefix", default="chr", help="prefix the VCF's contig names use (default chr)")
    ap.add_argument("--gcs-project", default=os.environ.get("GOOGLE_PROJECT", ""))
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    sites = read_sites(a.sites)
    labels = a.label or [Path(v).name for v in a.vcf]
    if len(labels) != len(a.vcf):
        ap.error("--label must be given as often as --vcf")
    frames = []
    for vcf, label in zip(a.vcf, labels):
        chrom_sites = sites
        got = query_sites(vcf, chrom_sites, tuple(a.fields.split(",")), a.gcs_project or None, a.chrom_prefix)
        got["source"] = label
        print(f"{label}: {len(got)} site x sample rows for {got.locus_id.nunique() if len(got) else 0} loci",
              file=sys.stderr)
        frames.append(got)
    out = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    out.to_csv(a.out, sep="\t", index=False, compression="gzip" if str(a.out).endswith(".gz") else None)
    print(f"wrote {a.out}: {len(out)} rows", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
