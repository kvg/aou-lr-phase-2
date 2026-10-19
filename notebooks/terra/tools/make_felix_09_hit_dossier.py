#!/usr/bin/env python3
"""Generate ../felix_09_hit_dossier.ipynb (notebooks/terra/felix_09_hit_dossier.ipynb).

The notebook is written from here rather than edited by hand so the cell order and the markdown stay reviewable in a
diff. Re-run it to regenerate; ids are derived from the cell index so the file is stable between runs.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

NB = Path(__file__).resolve().parents[1] / "felix_09_hit_dossier.ipynb"

CELLS: list[tuple[str, str]] = []


def md(text: str) -> None:
    CELLS.append(("markdown", text.strip("\n")))


def code(text: str) -> None:
    CELLS.append(("code", text.strip("\n")))


md("""
# FELIX hit dossier — one page of diagnostics per significant locus

Assembles everything needed to judge whether each retained hit is real, short of the read pileups, and exports it for
the review plate (`manuscript_figures/figS_hit_dossier.py`, one page per locus).

**What a locus is here.** A hit the screen retained (`verdict == 'keep'`: none of spa_fallback, no_case_carriers,
impossible_p, excess_het fired) with at least 10 case carriers. Hits of one phenotype within 500 kb are one locus.
That is the same rule the Figure 3 table and PheWAS panel count, enforced by a test, so the plate reviews exactly the
hits the figure reports.

**What this notebook writes**

| file | one row per | cost |
|---|---|---|
| `felix_dossier_loci.tsv` | locus: identity, lead, arms, phenotype calibration, window | free |
| `felix_dossier_leads.tsv` | locus x arm: the lead marker's full per-ancestry block | free (harvest cache) |
| `felix_dossier_context.tsv` | locus: cytoband, segmental duplications, tandem repeat, genes | UCSC API |
| `felix_dossier_genes.tsv` | gene in a locus window, for the regional panel's gene track | UCSC API |
| `felix_dossier_region.<locus>.tsv.gz` | marker near the locus, thinned, all arms | **one stream per run** |

Everything except the region files comes out of the harvest cache that `felix_08_hit_sweep` already built, so sections
1-6 take about a minute. Section 7 is the expensive one: it streams each phenotype's `merged.felix.tsv` (0.5-27 GB)
once for all of that phenotype's loci.

**Run order.** `00_sync_repo` (so the bucket copies of `felix_phenotype_report.py`, `felix_hit_dossier.py` and
`extract_region.py` are current), then `felix_08_hit_sweep` (builds `harvest_cache/`), then this notebook.
""")

md("## 1 · Setup")

code("""
import os, json, subprocess, sys, warnings, importlib
from pathlib import Path
import numpy as np, pandas as pd
from IPython.display import display, Markdown
warnings.filterwarnings("ignore", category=FutureWarning)
pd.set_option("display.width", 220); pd.set_option("display.max_columns", 80)

ETYPE    = "flare_lai_prod_set"
ID_COL   = f"entity:{ETYPE}_id"
TAR_COL  = "run_summary_tar"
NUM_ANCS = 5
STORE    = "harvest_cache"      # built by felix_08; shared, not rebuilt here unless it is missing
REGION_DIR = Path("dossier_regions")   # compact per-locus region files land here

WS_NS   = os.environ.get("WORKSPACE_NAMESPACE")
WS_NAME = os.environ.get("WORKSPACE_NAME")
BUCKET  = os.environ.get("WORKSPACE_BUCKET")
print("workspace:", WS_NS, "/", WS_NAME)
assert WS_NS and WS_NAME, "set WS_NS / WS_NAME by hand if the env vars are absent"

# ALWAYS re-copy the modules. Several notebooks write these filenames into this working directory, and an older sync
# leaves an older copy, so a "fetch only if absent" guard silently imports stale code.
MODULES = ["felix_phenotype_report.py", "felix_hit_dossier.py", "extract_region.py",
           "pull_site_fields.py"]
for m in MODULES:
    subprocess.run(["gsutil", "cp", f"{BUCKET}/felix/scripts/{m}", m], check=True)

importlib.invalidate_caches()
for m in ("felix_phenotype_report", "felix_hit_dossier", "extract_region", "pull_site_fields"):
    sys.modules.pop(m, None)
import felix_phenotype_report as fpr
import felix_hit_dossier as fhd
import extract_region as er
import pull_site_fields as ps

REQUIRED = {
    "felix_phenotype_report": (fpr, ["sweep", "harvest", "read_harvest", "hit_report", "calibration_report",
                                     "excess_het", "corroborated", "recurrence"]),
    "felix_hit_dossier":      (fhd, ["dossier_loci", "lead_rows", "annotate_context", "region_plan",
                                     "compact_region", "checks", "verdict"]),
    "extract_region":         (er,  ["extract_many", "read_windows", "split_windows"]),
    "pull_site_fields":       (ps,  ["query_sites", "call_rate_by_group", "imputed_fraction_by_group",
                                     "case_control_map", "available_format_tags"]),
}
for name, (mod, funcs) in REQUIRED.items():
    missing = [f for f in funcs if not hasattr(mod, f)]
    if missing:
        raise RuntimeError(
            f"{name}.py is missing {missing}.\\n"
            f"The copy at {BUCKET}/felix/scripts/{name}.py predates the dossier work. "
            f"Re-run 00_sync_repo (it rsyncs felix/scripts/ to the bucket from origin/main), then re-run this cell.")
print("modules ok:", ", ".join(REQUIRED))
""")

md("## 2 · Run table")

code("""
from firecloud import api as fapi

def entities_to_df(ns, name, etype):
    r = fapi.get_entities(ns, name, etype); r.raise_for_status()
    rows = []
    for e in r.json():
        rec = {f"entity:{etype}_id": e["name"]}
        for k, v in e.get("attributes", {}).items():
            if isinstance(v, dict) and "items" in v:
                v = list(v["items"])
            rec[k] = v
        rows.append(rec)
    return pd.DataFrame(rows)

runs = entities_to_df(WS_NS, WS_NAME, ETYPE).sort_values(ID_COL).reset_index(drop=True)
done = runs[runs[TAR_COL].notna()].copy()
print(f"{len(runs)} rows in {ETYPE}; {len(done)} have {TAR_COL}")
""")

md("""
## 3 · Screen every run, then reduce the hits to loci

Re-runs the same screen as `felix_08` off the shared harvest cache (seconds when the cache is warm, a few minutes
cold), then applies the locus rule. The counts printed here are the ones the plate will show.
""")

code("""
%%time
cal, hits, failures = fpr.sweep(done, ID_COL, TAR_COL, store=STORE, num_ancs=NUM_ANCS, workers=8)
print(f"{len(cal)} run x arm calibration rows; {len(hits)} genome-wide hits; {len(failures)} harvest failures")
if len(failures):
    display(failures)

loci = fhd.dossier_loci(hits)
kept = fhd.retained(hits)
print(f"\\n{len(hits)} reported hits -> {int((hits.verdict == 'keep').sum())} retained by the screen "
      f"-> {len(kept)} with at least {fhd.MIN_CASE_CARRIERS} case carriers -> {len(loci)} loci "
      f"across {loci.phenotype.nunique()} phenotypes")
print("\\nby lead arm:", loci.lead_arm.value_counts().to_dict())
print("arms per locus:", loci.arms.value_counts().to_dict())
display(loci.head(10))
""")

md("""
## 4 · Lead statistics, with every per-ancestry column

The lead of each arm that retained the locus, matched inside the locus window on the arm's own `top_hits.tsv`. These
are the columns the effect forest, the genotype block and the joint-test coherence check are computed from, and they
are already in the harvest cache -- no streaming.
""")

code("""
%%time
top = {}
for run_dir in sorted(Path(STORE).iterdir()):
    if not run_dir.is_dir() or run_dir.name.startswith("_"):
        continue
    b = fpr.read_harvest(run_dir)
    for arm in fpr.ARMS:
        th = b["arms"].get(arm, {}).get("top_hits")
        if th is not None and len(th):
            top[(run_dir.name, arm)] = th
print(f"{len(top)} run x arm top_hits tables in the cache")

leads = fhd.lead_rows(loci, top)
found = leads[leads.found.astype(bool)]
print(f"{len(leads)} locus x arm rows; {len(found)} matched a top_hits row, {len(leads) - len(found)} did not")
if len(leads) > len(found):
    display(leads[~leads.found.astype(bool)][["locus_id", "arm", "note"]])

# Sanity: the lead p a dossier row reports must be the p the sweep reported for that locus and arm.
chk = found.merge(loci[["locus_id", "lead_arm", "p"]], on="locus_id")
chk = chk[chk.arm == chk.lead_arm]
delta = (pd.to_numeric(chk["P_cct_admixed_c"], errors="coerce") / chk["p"] - 1).abs()
print(f"lead p agrees with the sweep for {int((delta < 1e-6).sum())}/{len(chk)} loci "
      f"(largest relative difference {delta.max():.2e})")
display(found[["locus_id", "arm", "MarkerID", "P_cct_admixed_c", "P_het_admixed_c", "P_hom_admixed_c"]].head(10))
""")

md("""
## 5 · Phenotype context

Each locus carries its phenotype's case count and per-arm genomic-control lambda, so the page can say whether the
scan that produced the hit was itself calibrated. `lambda_gc` below about 0.9 means the arm is deflated and no hit
screen will change that.
""")

code("""
wide = cal.pivot_table(index="run", columns="arm", values=["lambda_gc", "n_cases", "N_eff"], aggfunc="first")
wide.columns = [f"{a}_{b}" for a, b in wide.columns]
wide = wide.reset_index()
keep_cols = ["run"] + [c for c in wide.columns if c.startswith("lambda_gc_")]
ctx_cols = wide[keep_cols].copy()
ctx_cols["n_cases"] = wide[[c for c in wide.columns if c.startswith("n_cases_")]].max(axis=1)
ctx_cols["N_eff"] = wide[[c for c in wide.columns if c.startswith("N_eff_")]].max(axis=1)

loci = loci.merge(ctx_cols, on="run", how="left")
print("lambda of the lead arm, across the dossier loci:")
lead_lambda = loci.apply(lambda r: r.get(f"lambda_gc_{r.lead_arm}"), axis=1)
print(lead_lambda.describe().round(3).to_string())
print(f"\\nloci whose lead arm has lambda below 0.90: {int((lead_lambda < 0.90).sum())}")
display(loci[["locus_id", "phenotype", "n_cases", "lambda_gc_main", "lambda_gc_ru", "lambda_gc_ru_biallelic"]].head(10))
""")

md("""
## 6 · Genomic context

Cytoband, segmental duplications, tandem repeats and genes from the UCSC hg38 API — the questions "is this
pericentromeric", "is it in a duplication" and "what gene is it in" that a reviewer asks before anything else. About
four requests per locus.
""")

code("""
%%time
context, genes, tracks = fhd.annotate_context(loci)
bad = context[context.error != ""]
print(f"{len(context)} loci annotated; {len(bad)} lookups failed")
if len(bad):
    display(bad[["locus_id", "error"]])
print(f"pericentromeric: {int(context.pericentromeric.sum())}; "
      f"inside a segmental duplication: {int(context.in_segdup.sum())}; "
      f"inside a tandem repeat: {int(context.in_tandem_repeat.sum())}; "
      f"{len(genes)} gene records and {len(tracks)} repeat/duplication features in the windows")
display(context[["locus_id", "cytoband", "in_segdup", "segdup_frac", "in_tandem_repeat", "gene",
                 "nearest_gene", "nearest_gene_bp"]].head(10))
""")

md("""
## 7 · Region extracts — the expensive step

The only diagnostic that needs the full results: does the peak have neighbours? One `gsutil cat` pass per run and arm
covers every locus of that phenotype, and each window is immediately reduced to the ~26 columns the regional panel
reads, keeping every marker below p 1e-3 plus a 3,000-row background sample. The raw cut is deleted as soon as it is
compacted, so peak disk stays in the hundreds of MB rather than ~1 GB.

Re-running is cheap for loci already done: set `REDO = True` to force them.
""")

code("""
%%time
REDO = False
REGION_DIR.mkdir(exist_ok=True)
plan = fhd.region_plan(loci)
by_id = {r.locus_id: r for _, r in loci.iterrows()}   # Series, so compact_region can read ["window_start"]
print(f"{sum(len(v) for v in plan.values())} windows across {len(plan)} runs")

written, skipped, problems = [], [], []
for i, (run, spans) in enumerate(sorted(plan.items()), 1):
    ids = [lid for s in spans for lid in s["locus_id"].split(",")]
    if not REDO and all(fhd.region_path(REGION_DIR, lid, arm).exists() for lid in ids for arm in fpr.ARMS):
        skipped += ids
        continue
    row = done[done[ID_COL] == run]
    if len(row) != 1:
        problems.append(dict(run=run, error=f"{len(row)} rows in the run table"))
        continue
    root = row[TAR_COL].iloc[0].split("/call-Bundle")[0]
    print(f"[{i}/{len(plan)}] {run}: {len(spans)} window(s)", flush=True)
    for arm in fpr.ARMS:
        objs = er.list_call(root, er.ARM_CALL[arm])
        uri = er.find_merged(objs)
        if not uri:
            problems.append(dict(run=run, arm=arm, error=f"no merged table under call-{er.ARM_CALL[arm]}"))
            continue
        targets = [(str(s["chrom"]), s["start"], s["end"], Path(f"_raw.{arm}.{j}.tsv"))
                   for j, s in enumerate(spans)]
        try:
            counts = er.extract_many(uri, targets)
        except Exception as exc:
            problems.append(dict(run=run, arm=arm, error=f"{type(exc).__name__}: {exc}"))
            continue
        for (span, (_, _, _, raw), n) in zip(spans, targets, counts):
            for lid in span["locus_id"].split(","):
                out = fhd.region_path(REGION_DIR, lid, arm)
                try:
                    k = fhd.compact_region(raw, out, by_id[lid])
                    written.append(dict(locus_id=lid, arm=arm, raw_rows=n, kept_rows=k,
                                        kb=round(out.stat().st_size / 1024, 1)))
                except Exception as exc:
                    problems.append(dict(run=run, arm=arm, locus_id=lid, error=f"{type(exc).__name__}: {exc}"))
            raw.unlink(missing_ok=True)

w = pd.DataFrame(written)
print(f"\\nwrote {len(w)} compact region files ({w.kb.sum() / 1024:.1f} MB total)" if len(w) else "\\nnothing written")
print(f"skipped {len(set(skipped))} loci already done; {len(problems)} problems")
if problems:
    display(pd.DataFrame(problems))
if len(w):
    display(w.head(10))
""")

md("""
## 8 · Call rate at the lead sites

The one diagnostic FELIX's own output cannot give. `phase_vcfs` is the ligated phased callset, where phasing has filled
every genotype, so `MissingRate` is 0 at every record — not a clean call rate, an absent one.

This pulls the genotypes back from the indexed callsets directly: one `bcftools query -R` per VCF, one index seek per
lead, nothing localised and nothing merged. Two sources per locus — the **pre-phasing** DeepVariant/GLnexus VCF and the
**phased** `gt_vcf` the scan actually read — which give two numbers: the case/control call rate, and the share of the
genotypes the scan used that phasing *imputed* rather than observed. A difference between cases and controls in the
second is the failure mode a call-rate comparison alone cannot see.

Fill `PRE_PHASING_VCFS` with the per-chromosome DeepVariant/GLnexus URIs (the `GL_INTERVAL_set` table that
`propagate_annotations/wdl/AnnotateFlareGqDp.wdl` names as its `deepvariant_vcf`). Leave it empty and the cell reports
the phased call rate only, which is 1.0 everywhere by construction and says so.
""")

code("""
%%time
# htslib range-reads gs:// when the environment carries credentials; without these every query fails on auth.
os.environ.setdefault("GCS_OAUTH_TOKEN", subprocess.run(["gcloud", "auth", "print-access-token"],
                                                        capture_output=True, text=True).stdout.strip())
os.environ.setdefault("GCS_REQUESTER_PAYS_PROJECT", os.environ.get("GOOGLE_PROJECT", ""))

# The phased callset the scan read, per chromosome, straight off the entity table.
chrom_rows = entities_to_df(WS_NS, WS_NAME, "flare_lai_prod")
PHASED_VCFS = {str(r["chrom"]): r["gt_vcf"] for _, r in chrom_rows.iterrows() if r.get("gt_vcf")}
print(f"phased gt_vcf for {len(PHASED_VCFS)} chromosomes")

PRE_PHASING_VCFS = {}   # {'chr9': 'gs://.../chr9.glnexus.vcf.gz', ...}  <- fill this
if not PRE_PHASING_VCFS:
    print("PRE_PHASING_VCFS is empty: call rate will be reported from the phased callset only, "
          "which is 1.0 by construction. Fill it to get a real call rate.")

sites = fhd.lead_sites(loci)
sites["chrom_name"] = "chr" + sites.chrom.astype(str)
# pheno_cov is an attribute of the run set itself (the config passes it as ${this.pheno_cov}), so take it from there
# rather than asking; every run in the table should name the same file.
uris = sorted({u for u in done.get("pheno_cov", pd.Series(dtype=str)).dropna().unique() if str(u).strip()})
if len(uris) != 1:
    raise RuntimeError(f"expected one pheno_cov across the run table, found {len(uris)}: {uris[:3]}. "
                       f"Set PHENO_COV by hand and re-run this cell.")
PHENO_COV = uris[0]
if PHENO_COV.startswith("gs://"):
    subprocess.run(["gsutil", "cp", PHENO_COV, "pheno_cov.tsv"], check=True)
    PHENO_COV = "pheno_cov.tsv"
pheno = pd.read_csv(PHENO_COV, sep="\t", low_memory=False)
print(f"pheno_cov ({uris[0]}): {len(pheno)} rows, {len(pheno.columns)} columns")
missing_pheno = sorted(set(sites.phenotype) - set(pheno.columns))
if missing_pheno:
    print(f"{len(missing_pheno)} dossier phenotypes have no column in pheno_cov, e.g. {missing_pheno[:5]}")

pulled = {}
for label, mapping in (("pre", PRE_PHASING_VCFS), ("post", PHASED_VCFS)):
    frames = []
    for chrom_name, uri in sorted(mapping.items()):
        want = sites[sites.chrom_name == chrom_name]
        if want.empty:
            continue
        try:
            got = ps.query_sites(uri, want, fields=("GT", "DP", "GQ", "RNC"),
                                 gcs_project=os.environ.get("GOOGLE_PROJECT"))
        except Exception as exc:
            print(f"  {label} {chrom_name}: {type(exc).__name__}: {str(exc)[:200]}")
            continue
        print(f"  {label} {chrom_name}: {len(got)} rows for {got.locus_id.nunique() if len(got) else 0} loci",
              flush=True)
        frames.append(got)
    pulled[label] = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()

# Each locus is its own phenotype, so the case/control split is per locus, not global.
rows = []
for _, S in sites.iterrows():
    try:
        groups = ps.case_control_map(pheno, S.phenotype)
    except KeyError:
        rows.append(dict(locus_id=S.locus_id, note=f"no {S.phenotype} column in pheno_cov"))
        continue
    rec = dict(locus_id=S.locus_id, phenotype=S.phenotype, note="")
    for label in ("pre", "post"):
        df = pulled[label]
        sub = df[df.locus_id == S.locus_id] if len(df) else pd.DataFrame()
        cr = ps.call_rate_by_group(sub, groups)
        if len(cr):
            for k, v in cr.iloc[0].items():
                if k != "locus_id":
                    rec[f"{label}_{k}"] = v
    if len(pulled["pre"]) and len(pulled["post"]):
        imp = ps.imputed_fraction_by_group(pulled["pre"][pulled["pre"].locus_id == S.locus_id],
                                           pulled["post"][pulled["post"].locus_id == S.locus_id], groups)
        if len(imp):
            for k, v in imp.iloc[0].items():
                if k != "locus_id":
                    rec[k] = v
    rows.append(rec)

callrate = pd.DataFrame(rows)
callrate.to_csv(fhd.CALLRATE_FILE, sep="\t", index=False)
print(f"\\nwrote {fhd.CALLRATE_FILE}: {len(callrate)} loci")
cols = [c for c in ("locus_id", "pre_call_rate_case", "pre_call_rate_control", "pre_fisher_p",
                    "imputed_frac_case", "imputed_frac_control", "imputed_frac_diff") if c in callrate.columns]
if len(cols) > 1:
    display(callrate[cols].head(12))
""")

md("## 8 · Export")

code("""
loci.to_csv(fhd.LOCI_FILE, sep="\\t", index=False)
leads.to_csv(fhd.LEADS_FILE, sep="\\t", index=False)
context.to_csv(fhd.CONTEXT_FILE, sep="\\t", index=False)
genes.to_csv(fhd.GENES_FILE, sep="\\t", index=False)
tracks.to_csv(fhd.TRACKS_FILE, sep="\\t", index=False)
print("wrote:", sorted(p.name for p in Path(".").glob("felix_dossier_*.tsv")))
print(f"region files: {len(list(REGION_DIR.glob('*.tsv.gz')))}")

# A dry run of the checks, so the notebook says what the plate will say before anyone renders it.
rows = []
for _, L in loci.iterrows():
    lr = leads[(leads.locus_id == L.locus_id) & (leads.arm == L.lead_arm) & leads.found.astype(bool)]
    lead = lr.iloc[0] if len(lr) else None
    cx = context[context.locus_id == L.locus_id]
    reg_path = fhd.region_path(REGION_DIR, L.locus_id, L.lead_arm)
    region = pd.read_csv(reg_path, sep="\\t") if reg_path.exists() else None
    flags = fhd.checks(L, lead, cx.iloc[0] if len(cx) else None, region)
    rows.append(dict(locus_id=L.locus_id, phenotype=L.phenotype, p=L.p, verdict=fhd.verdict(flags),
                     **{f["check"]: f["status"] for f in flags}))
review = pd.DataFrame(rows)
review.to_csv("felix_dossier_review.tsv", sep="\\t", index=False)
display(Markdown("### Verdict over the dossier loci"))
print(review.verdict.value_counts().to_string())
for c in [c for c in review.columns if c not in ("locus_id", "phenotype", "p", "verdict")]:
    print(f"  {c:24s} " + ", ".join(f"{k} {v}" for k, v in review[c].value_counts().items()))
display(review[review.verdict.isin(["fail", "warn"])].head(20))
""")

md("## 9 · Bundle for download")

code("""
!rm -f felix_dossier.tar.gz felix_dossier.tar
!tar -cf felix_dossier.tar felix_dossier_*.tsv dossier_regions/
!gzip felix_dossier.tar
!gsutil cp felix_dossier.tar.gz $WORKSPACE_BUCKET/scratch/kvg/
!ls -la felix_dossier.tar.gz && tar -tzf felix_dossier.tar.gz | head -20
""")

md("""
## Caveats

- **The screen comes first.** Every locus here already passed `hit_report`, so the plate reviews survivors, not the
  5,774 excluded hits. The exclusion counts are in `felix_sweep_all_hits.tsv`; the plate's header reports them so the
  reader can see what fraction of the scan reached this document.
- **Only the top-50 rows of each bundle are screenable.** A hit outside a run's `top_hits.tsv` is invisible to both the
  screen and this dossier, and those concentrate in the lowest-case-count runs.
- **No read pileups.** Everything here is derived from the association tables. A locus that passes every check still
  has not been looked at directly, and for the handful that survive, the pileup is the next step, not the last one.
- **The repeat-dosage (`ru`) arm has no genotype counts.** Its `N_*_het/hom` columns hold dosage sums, so
  Hardy-Weinberg, carrier fractions and the fixed-heterozygote check do not apply to it and report `n/a`.
- **Call rate cannot come from FELIX's own input, but it is not out of reach.** `phase_vcfs` is
  `flare_lai_prod.gt_vcf`, the ligated phased callset, and statistical phasing fills every genotype: `MissingRate` is
  exactly 0 in all 331,518 cells of the five region extracts, so the page reports it as not assessed rather than as a
  pass. The call rate, and the fraction of genotypes that phasing **imputed** rather than observed, live in the
  pre-phasing DeepVariant/GLnexus callset, which `propagate_annotations/wdl/AnnotateFlareGqDp.wdl` already reaches
  (it copies `FORMAT/GQ`, `DP` and `RNC` onto FLARE sites, region-indexed). 88 of the 90 leads are SNV or indel and so
  are answerable that way; 1 SV lead needs the SV callset and 1 repeat lead the repeat branch.
- **Context is hg38 UCSC annotation**, not a mappability or coverage track: "unique sequence" here means no overlapping
  segmental duplication record, which is weaker than a mappability estimate from the actual reads.
""")


def check_parses(src: str, index: int) -> None:
    """Fail here, not on Terra. Cell magics and ! shell lines are stripped before parsing.

    This exists because the code cells above are ordinary triple-quoted strings, so a backslash escape meant for the
    notebook is interpreted once by this file first: `print(f"\\n...")` here is what puts `\n` in the notebook, and
    writing one backslash instead put a real newline inside an f-string, which only showed up as a SyntaxError when the
    cell was run in Terra.
    """
    import ast

    body = src
    if body.lstrip().startswith("%%"):
        body = body.split("\n", 1)[1] if "\n" in body else ""
    body = "\n".join(l for l in body.split("\n") if not l.lstrip().startswith("!"))
    if not body.strip():
        return
    try:
        ast.parse(body)
    except SyntaxError as exc:
        lines = body.split("\n")
        near = lines[max(0, (exc.lineno or 1) - 2):(exc.lineno or 1) + 1]
        raise SystemExit(f"cell {index} will not compile: {exc.msg} at line {exc.lineno}\n  "
                         + "\n  ".join(repr(l) for l in near)) from exc


def build() -> dict:
    cells = []
    for i, (kind, src) in enumerate(CELLS):
        if kind == "code":
            check_parses(src, i)
        cid = hashlib.sha1(f"felix_09:{i}".encode()).hexdigest()[:8]
        cell = {"cell_type": kind, "id": cid, "metadata": {},
                "source": [l + "\n" for l in src.split("\n")[:-1]] + [src.split("\n")[-1]]}
        if kind == "code":
            cell["execution_count"] = None
            cell["outputs"] = []
        cells.append(cell)
    return {
        "cells": cells,
        "metadata": {
            "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
            "language_info": {"name": "python", "version": "3.10"},
        },
        "nbformat": 4,
        "nbformat_minor": 5,
    }


if __name__ == "__main__":
    NB.write_text(json.dumps(build(), indent=1, ensure_ascii=False) + "\n")
    print(f"wrote {NB} ({len(CELLS)} cells)")
