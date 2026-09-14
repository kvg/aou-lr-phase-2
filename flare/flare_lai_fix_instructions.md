# FLARE LAI fix instructions

Parts **0–3** (model I/O, props/µ blocks, switch-QC metrics, blocked chr20
prep) are implemented in this repo and remain prerequisites for any grid
rerun. See the plan file and `flare/README.md`. **Do not reorder:** Part 0
(model parser) and Part 1.1 (`out_model` glob) must land before any Terra
rerun that interprets models or pin-µ rows.

Parts **4–5** below are the site-restriction addendum.

---

## Part 4.0: measure what FLARE is already discarding, before adding filters

Nothing in the current workflow restricts variant class. `SubsetOnePopulation`
runs `bcftools view -S <keep> -r <region>` and `CutRefRegion` runs
`bcftools view -r <region>`, so no `-m2 -M2`, no `-v snps`, and no mask. The
only site-level thinning is FLARE's own `min-maf=0.005` and `min-mac=50`.

FLARE may drop multiallelic markers internally. Determine this from data
rather than assuming it either way:

1. Read the per-shard FLARE log (already captured as `per_pop_log`) and extract
   the marker count FLARE reports as retained (`Statistics` → `markers : N`).
2. Compare against `bcftools index -n` on that shard's subset VCF, and against
   a breakdown of that VCF by variant class:
   `bcftools view -H -v snps -m2 -M2 | wc -l` versus indels versus
   multiallelic.
3. Write the result into the experiment summary as columns
   `n_sites_input`, `n_sites_flare_retained`, `frac_biallelic_snv`.

If FLARE already retains only biallelic SNVs, the `biallelic_snvs_only` input
below is a no-op and the remaining leverage is entirely in the context mask and
the concordance list. That changes what is worth running, so do this step
first and report the numbers before implementing 4.1.

**Code:** `scripts/flare_site_stats.py` (log parse + optional VCF class
breakdown). Workflow emits `per_pop_site_stats` / merged `site_stats.tsv`.

**Source note (pre-Terra):** FLARE's Java `Marker` accepts multi-allelic
records (`nAlleles` up to a large max); common markers require identical
CHROM/POS/REF/ALT. So `biallelic_snvs_only` is **not** assumed a no-op until
Part 4.0 numbers from a real shard confirm otherwise. Live counts need
workspace credentials (local GCS was AccessDenied).

**Part 4.0 status:** tooling is in place (`flare_site_stats.py`, per-shard
`site_classes.tsv`, merged `site_stats.tsv` with `n_sites_input` /
`n_sites_flare_retained` / `frac_biallelic_snv`). **Report live numbers from
the next Terra submission** (or a notebook pass over existing `per_pop_log` +
subset VCF) before treating step 2 of the filter ladder as optional. Until
those numbers exist, keep `biallelic_snvs_only` on the ladder.

---

## Part 4.1: site restriction inputs

Add to the workflow and thread into both `SubsetOnePopulation` and
`CutRefRegion`:

```wdl
Boolean biallelic_snvs_only = false   # bcftools view -m2 -M2 -v snps
File? include_sites          # BED or tabixed sites file: keep only these
File? exclude_regions        # BED: segdups, repeats, low-complexity, indel flanks
Int   indel_flank_bp = 0     # pad around indels in the target when > 0
```

Requirements:

- **Apply the identical site treatment to the target and the reference.** Fail
  the workflow if `include_sites` or `exclude_regions` is set and the two tasks
  did not receive the same value (checksum gate).
- Order the bcftools operations region, then sites, then samples, using
  `-R`/`-T` against the tabix index. Keep `bcftools_view_retry` + auth refresh.
- Extend `counts.tsv` to
  `population, n_samples, n_sites_pre_filter, n_sites_post_filter` (and the
  same for the reference slice).
- Hard-fail on zero retained sites; also fail when post-filter density falls
  below `min_sites_per_mb` (default 1000).

Call-rate filter: target is documented as phased with no missing GT — verify
on one shard; if complete, document and skip (do not add a silent no-op filter).

Concordance `include_sites` list is out of scope to build in the WDL; the file
input is the drop-in hook. Stage both resources on Terra with
`notebooks/terra/flare_03_stage_site_filters.ipynb` (context mask always;
concordance when HiFi + srWGS VCFs are set).

---

## Part 4.2: per-population frequency and HWE

FLARE's `min-maf` / `min-mac` operate on the **reference** panel. Target-side
frequency / HWE filtering is separate and must run **per shard after sample
subsetting** (`bcftools +fill-tags` + expression exclude inside
`SubsetOnePopulation`). Record how many sites each criterion removed.

---

## Part 5: metric changes the filters force

Thinning sites reduces switch opportunity, so `switches_per_hap_per_mb` falls
mechanically. Add to the comparison table:

- `markers_per_mb`
- `switches_per_hap_per_marker`
- filter provenance: `biallelic_snvs_only`, `include_sites` basename,
  `exclude_regions` basename

Require claims of improvement on **both** per-Mb and per-marker rates. Guard
against comparing rows whose filter provenance tuples differ.

### Rerun ladder (after 4.0 numbers)

Single pinned T on chr20, props held at Part 3 choice:

1. current site handling (baseline)
2. plus biallelic SNVs only, if 4.0 showed this is not already the case
3. plus the context mask
4. plus the concordance list (deferred; not in `lai_exp.tsv` for now)

Theory: per-marker switch rate falls across steps 2–4; if not, platform
mismatch is wrong and the cause is phasing or the panel.

---

## Part 6: FLARE2 (clustering for poorly-matched panels)

FLARE2 ([Browning et al. 2025](https://doi.org/10.1101/2025.10.13.681993);
PLOS Genetics) is **not** a new LAI engine. It is a three-step recipe on
FLARE 0.6 that learns a copying matrix \(P\) when ancestries and reference
panels are not one-to-one:

1. `panel-probs=true` — ancestry-agnostic Li–Stephens; write `.panels`
2. `create-model-file.py nanc panels-file out-prefix` — GMM cluster → `.model`
   (initial \(P\), ρ; flat μ; default θ). Check **min lag-1 ancestry
   autocorrelation ≥ ~0.25**; below that, `nanc` is likely spurious.
3. FLARE again with `model=` that file and `update-p=true` (already our default)

**Where it helps us:** AMR (Indigenous poorly matched by Asian/Amerind proxies)
and MID (Mozabite-like North African / ME). **Where it does not:** AFR/EUR with
good gnomAD/1KG matches (original FLARE can be slightly better); HiFi genotype
error / short tracts (same θ emission model). Prefer **STEAM** \(g\) over FLARE
\(T\) when quoting admixture time from FLARE2 calls.

### Experiment rows (`lai_exp.tsv`)

| Row | Intent |
|---|---|
| `flare2_amr_nanc2` | AMR-only; `nanc=2`; vs `pin_amr_only_gen12` |
| `flare2_amr_nanc3` | AMR-only; `nanc=3`; drop if autocorrelation &lt; 0.25 |
| `flare2_mid_nanc2` | MID-only; `nanc=2` (primary FLARE2 use case) |
| `flare2_afr_nanc2` | AFR-only control; expect little/no gain vs `pin_afr_only_gen8` |

All on the chr22 10 Mb window first. Score with the same switch-QC + global
props as other rows; record clustering autocorrelation and STEAM \(g\) when
available.

### WDL work (required before launching those rows)

`FlareByPopulation` today runs a single FLARE pass per pop. Add a `flare2`
mode (table column `flare2_nanc`, 0 = off):

1. Per-pop (or pooled) `flare_task` with `panel-probs=true`
2. Task wrapping upstream `create-model-file.py` (`nanc=flare2_nanc`)
3. Second `flare_task` with that model + `update-p=true` (reuse existing path)

Until that lands, FLARE2 rows are **blocked** (do not launch; notes say so).
Software already ships in FLARE 0.6 (`panel-probs`, `create-model-file.py`).

---

## Part 7: GQ / DP / RNC call-quality filter

Hypothesis: short tracts / excess switches are enriched at poor DeepVariant /
GLnexus calls (low GQ, low DP, RNC=`I`). Evaluate in **two layers** — do not
conflate them.

### 7.1 Post-hoc (no FLARE rerun) — effect on scored annotations

Requires `AnnotateFlareGqDp` output (FORMAT GQ/DP/RNC on FLARE sites).

`flare_switch_qc.py` already stratifies switches by quality. Extend scoring to
report rates **after dropping** switches that fail call QC:

- `GQ < gq_threshold` (default 20) **or**
- `DP < dp_threshold` (default 10) **or**
- `RNC` contains `I` (incomplete gVCF / no-call reason)

Compare `switches_per_hap_per_mb` / per-marker **all** vs **call_qc_pass**.
If most of the excess vanishes after dropping poor calls, the LAI model is
less at fault than genotype quality — prefer masking those sites for
downstream SV ancestry annotation rather than chasing T/μ.

Notebook: `flare_01_switch_gq_dp.ipynb` (annotate) → `flare_02` / switch-QC
with the filtered block in `summary.json`.

### 7.2 Pre-FLARE site filter (rerun ladder step)

Build an `include_sites` BED from the joint DeepVariant/GLnexus VCF
(`scripts/flare_build_call_qc_sites.py`): keep biallelic SNVs where the
fraction of samples failing GQ/DP/RNC is below a cap (default 5%). Stage via
`flare_03_stage_site_filters.ipynb` Part C. Apply identically to target and
reference (Part 4.1 digest gate).

Ladder (chr20 or chr22 slice; pinned T held fixed):

1. baseline site handling
2. biallelic SNVs only
3. + context mask
4. + **call-QC `include_sites`** (`chr22_filter_call_qc` / `chr20_filter_call_qc`)
5. concordance list (still deferred)

Score both per-Mb and per-marker rates; require improvement on both. Also run
7.1 on the annotated outputs so pre-filter and post-hoc effects are separable.

**Caveat:** FLARE forbids missing GT — we filter **sites**, not genotypes.
Phased `gt_vcf` may lack GQ/DP; the BED must come from the unphased joint VCF.

---

## Part 8: association-facing recipe evaluation (two-stage)

Switch / tract metrics (Parts 3–5, 7) are **diagnostics only**. Recipe choice
for Tractor / FELIX uses a two-stage association-facing rule:

```
fixed AF panel ──► allele–ancestry concordance (gate)
                └► pedigree Mendelian LAI (gate)
                         │
                    clears both?
                    ├── no  → eliminate
                    └── yes → Tractor null-λ (decisive) → winner closest to λ=1
```

Screens do **not** rank. Among survivors, pick the recipe with mean \|λ−1\|
closest to zero. If bootstrap CIs on \|λ−1\| overlap, tie-break with Part 2
`mean_ll` **only inside that overlapping set**.

### Part 8.1 Fixed AF-divergent marker panel

`scripts/flare_build_af_panel.py` ranks biallelic SNVs in the gnomAD LAI ref
by cross-panel AF range (`eas,amr,eur,afr,sas`), keeps top ~5k / window with a
MAC floor, writes `*.markers.tsv` + `*.bed.gz`.

Stage once (chr22 method window + chr20) via Terra notebook (preferred):

```text
notebooks/terra/flare_04_stage_eval_panels.ipynb
# → $WORKSPACE_BUCKET/refs/flare/eval_panels/chr22_10mb.markers.tsv
# → $WORKSPACE_BUCKET/refs/flare/eval_panels/chr20_full.markers.tsv
```

Or from a shell (same builder):

```bash
export WORKSPACE_BUCKET=gs://...
export REF_PANEL=.../aou_1000genomes.refmap
export REF_VCF_CHR22=.../chr22....vcf.bgz
export REF_VCF_CHR20=.../chr20....vcf.bgz
bash flare/scripts/stage_eval_af_panels.sh
```

Windows: `chr22:26897597-36897597` and full `chr20` (override `CHR20_REGION`).
Score finished `lai_exp` rows on their **native** `region`.

### Part 8.2 Allele–ancestry concordance (gate)

`scripts/flare_score_allele_ancestry.py`

- Alleles from the **unfiltered** phased `gt_vcf` for the window.
- Ancestry from recipe `anc.vcf` with FELIXla / `propagate_flare_ancestry`
  covering (first LAI POS ≥ query, else last).
- Gate metric: `mean_ll` (mean per-haplotype log-likelihood under panel AF for
  inferred ancestry). `mean_brier` is secondary, not the gate.
- Output: `allele_ancestry_score.json`.

### Part 8.3 Pedigree Mendelian LAI (gate)

`scripts/flare_score_mendelian_lai.py` + `tractor_mix/resources/legacy_covariates/aou_phase2.ped`

- Restrict to complete Phase 2 families when covariates are provided (same
  spirit as tableS pedigrees).
- Hard violations = child ancestries incompatible with parents at a locus.
- Assignment flips along the window count as recombinations; compare to map
  expectation (`maps/plink.chrchr{N}.GRCh38.map`, 2 meioses).
- Gate metrics: `violations_per_informative_locus`,
  `excess_recomb_over_expected`.

### Part 8.4 Tractor null-λ (decisive)

`scripts/flare_lai_null_lambda.py`

1. `bcftools view -r WINDOW` on recipe `anc.vcf`
2. `extract-tracts-flare --num-ancs 5 --samples analysis_samples`
3. K=50 phenotype permutations (shuffle an existing binary column in
   `pheno_cov.tsv`, fixed seed) → `fit_null.R` (reuse GRM) →
   `tractor-mix-score` with lowered `ac_threshold` (default 20) → λGC
4. Aggregate mean \|λ−1\| + bootstrap CI; flag `unstable_lambda` if
   `<100` tested sites after AC filter (widen window / lower AC; do not
   silently compare).

### Part 8.5 Positive controls (optional corroboration)

Literature check in parallel — not blocking. If a validated in-window hit
exists later, run Tractor on that trait for survivors and attach effect
recovery beside λ. Do **not** plant phenotypes from the LAI recipe.

### Part 8.6 Notebook + gates

`notebooks/terra/flare_02_lai_exp_compare.ipynb` — section **Part 8** after
switch-QC diagnostics:

1. Build/load fixed panels  
2. Score finished rows (8.2–8.3)  
3. Gate table + survivors (`flare_lai_exp.apply_eval_gates`)  
4. Null-λ on survivors → declare winner (`flare_lai_null_lambda.pick_winner`)

Write empirical thresholds into
`$WORKSPACE_BUCKET/refs/flare/eval_panels/eval_gates.json` after the first
pass (template: `flare/configs/eval_gates.json.example`). No hardcoded
cutoffs in v1 scorers.

Helpers: `flare_lai_exp.enrich_association_scores`,
`flare_lai_exp.apply_eval_gates`.

```bash
python3 scripts/test_flare_lai_association_metrics.py
```
