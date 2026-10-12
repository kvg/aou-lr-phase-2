# FELIX LAI GWAS (Terra)

Chr22 pilot and genome-wide workflows for FELIX local-ancestry-aware association on
All of Us / Terra. Consumes the published `lhu1/felix:latest` image (do not compile
or fork FELIX). Shares cohort inputs with Tractor-Mix / SAIGE pilots under
`tractor_mix_pilot/`.

**Plans, status, and agent coordination:** see [`PLAN.md`](PLAN.md).

## Layout

```
felix/
  build_docker.sh              # Cloud Build or local docker (repo-root context)
  docker/Dockerfile            # FROM lhu1/felix:latest + PLINK2/bcftools + Rust CLIs
  docker/cloudbuild.yaml
  wdl/FelixPilot.wdl           # chr22: mtx GRM + felixla + admixed step2
  wdl/FelixGenome.wdl          # autosomes: shared null, per-chr pack/step2, optional RU_TEST scorer
  wdl/FitFelixNulls.wdl        # standalone step-1 nulls, one per phenotype, reusable across runs
  configs/felix*.json.example
  scripts/                     # FELIX-specific CLIs (staged to gs://BUCKET/felix/scripts/)
  eval/README.md               # M1/M4 calibration gates
scripts/                       # shared: GRM build, make_plink_keep, annotate_repeat_units, …
tractor_mix/
  extract_tracts_flare/        # RU dosage extraction (used by FelixGenome RU_TEST branch)
  tractor_mix_score/           # Rust scorer (--mode felix)
```

## Design

FELIX replaces Tractor-Mix for SNV/indel + unique-SV association (`felixla` pack +
`step2_SPAtests` with `--is_admixed=TRUE`). Repeat-mediated SVs use the optional
repeat-dosage branch, which tests **one catalog repeat locus at a time**:
`scripts/aggregate_repeat_loci.py --out-vcf --with-ancestry` (every
length-changing record inside a TRExplorer interval contributes its signed
length change to that haplotype's locus dosage) →
`scripts/write_admixed_dosage_vcf.py` (REF-relative repeat units per ancestry,
`DS{k}` / `ANC{k}`) → FELIX step 2 on that VCF with carrier-count QC, using the
same FELIX null as the SNV scan (`ru_engine = "felix"`, default).
`ru_engine = "rust"` keeps the old record-level `annotate_repeat_units.py` →
`extract-tracts-flare` → `tractor-mix-score` encodings (no SPA) for comparison
only — it tests each VCF record separately and drops records whose alleles are
not a whole number of repeat units.

**Design, coding rules, status and open items:
[`REPEAT_DOSAGE.md`](REPEAT_DOSAGE.md)** (single source of truth).
Engine-choice history: [`SV_SCORER_DESIGN.md`](SV_SCORER_DESIGN.md).
Fixture evidence for the locus-level switch:
[`eval/record_vs_locus/`](eval/record_vs_locus/).

| Model | Covariates (2×2 calibration) |
|-------|------------------------------|
| FELIX limited | sex, age, PC1–PC10 |
| FELIX full | limited + standardized coverage + GC dummies |

Scan: chr22 FLARE / phased VCFs (pilot) or parallel `chroms` + `phase_vcfs` +
`flare_vcfs` (genome). GRM markers: chr1 + chr22 (same as Tractor-Mix / SAIGE).

## Build image

```bash
felix/build_docker.sh                 # Cloud Build → us-central1-docker.pkg.dev/.../felix-pilot:0.2.0
felix/build_docker.sh --local --push  # or local docker
```

Bump `FelixPilot.docker` / `FelixGenome.docker` in config JSON when retagging.

The image starts from FELIX v0.1 pinned by digest and applies
[`patches/`](patches/) to FELIX's source before reinstalling its R package:

| Patch | Effect |
|---|---|
| `0001-dosage-carrier-qc.patch` | With `FELIX_DOSAGE_QC=carrier` (set by `run_felix_step2.R --dosage-qc carrier`), the admixed dosage-VCF path filters on carrier counts instead of `sum(DS)`. Off by default; FELIXla results are unchanged. |

## Stage scripts

FELIX scripts live under `felix/scripts/`; shared helpers stay in repo `scripts/`.
`notebooks/terra/00_sync_repo.ipynb` stages both (`gs://BUCKET/felix/scripts/`,
`gs://BUCKET/scripts/`). From another shell with bucket access:

```bash
WORKSPACE_BUCKET=gs://fc-secure-... ./scripts/stage_felix_scripts.sh
WORKSPACE_BUCKET=gs://fc-secure-... ./scripts/stage_tractor_scripts.sh   # GRM, annotate_repeat_units, …
```

## Run order

1. Prepare cohort with `notebooks/terra/tractor_01_prepare_inputs.ipynb` (same as Tractor-Mix).
2. Submit **FelixPilot.wdl** with `felix/configs/felix.pilot.inputs.{limited,full}.json.example`.
   Point `gs://BUCKET/...` at your workspace bucket.

   Workflow: **Check** → **MakeGRM** (`build_saige_plink_and_grm.sh`) → **FitFelixNull**
   (`fit_felix_null.R` → FELIX step1 + `export_felix_null.R`) → **Pack** (`felixla`) →
   **Step2** (`run_felix_step2.R` / `step2_SPAtests.R`, `--is_admixed=TRUE`) → summarize
   (FelixPilot: `summarize_tractor_genome_results.py --p-column P_cct_admixed_c`; FelixGenome: `summarize_felix_results.py`; see Notes).

3. Compare λGC vs Tractor-Mix 2×2: `felix/eval/README.md`.

4. Genome-wide: **FelixGenome.wdl** + `felix/configs/felix.genome.inputs.*.json.example`.
   Repeat dosage: set `repeat_catalog_bed` (TRExplorer intervals),
   `aggregate_repeat_loci_script` and `write_admixed_dosage_script`, and pass
   `joint_vcfs` so the branch sees `GT:AN1:AN2`. Set `num_ancs` to the FLARE2
   model's `nanc` — `aggregate_repeat_loci.py --num-ancs` fails the task rather
   than folding an out-of-range ancestry code into another ancestry.

5. Long phenotype lists: fit the step-1 nulls once with **FitFelixNulls.wdl**
   (`felix/configs/fit_felix_nulls.inputs.json.example`) and pass its
   `null_rdas` / `variance_ratios` / `samples_used` / `null_export_tars` arrays
   to FelixGenome as `precomputed_null_*`. A null depends only on the GRM, the
   covariates and the phenotype column, so the same fits serve every
   chromosome, both repeat-coding arms and every resubmission — at P
   phenotypes the scan is P null fits plus 66P step-2 shards (22 chromosomes ×
   SNV + two repeat arms), and only the first term is avoidable.

   Order is the contract: the arrays are indexed against
   `read_lines(selected_phenotypes)`, so give both workflows the same
   phenotype file and check `null_manifest` before reusing the arrays
   elsewhere. `CheckNullArrays` fails the run on a length mismatch before any
   step-2 shard starts. Nulls are cache-keyed on `pheno_cov`,
   `covariate_columns`, `analysis_samples`, the GRM files and
   `fit_felix_null.R`, so settle the covariate arm before fitting hundreds.

   Build the list with `scripts/select_phenotypes.py`. Past the pilot's top-10
   it needs its scale flags: `--n-phenotypes 0` (keep everything passing),
   `--sex-specific restrict` with `--phecode-sex-tsv` so a sex-specific phecode
   does not carry the other sex as controls, `--nesting leaves` or `roots` for
   the parent/child hierarchy, and `--max-r2 0.5` to prune near-duplicate
   phecodes (what All by All used). Counts are of phecodes, not independent
   traits; `phenotype_selection_audit.tsv` records why each candidate is in or
   out, and `--tier-sizes` writes smaller prefixes for a shakedown run.

### Scaling to the full phecode list

`pheno_cov.tsv` does not change shape: `ID`, the covariate columns, then one
column per phecode. Three things about it are already handled and do not need
manual work.

- **Sex-specific phecodes.** `select_phenotypes.py --sex-specific restrict`
  (with `--phecode-sex-tsv`) sets the other sex to `NA` so they are not counted
  as controls.
- **The covariate that goes constant.** On a female-only phecode, `sex` has no
  variance among the complete cases and would make the design singular.
  `fit_felix_null.R` computes each covariate's SD on that phenotype's complete
  cases, drops the zero-variance ones, records them in
  `covariates_dropped_constant`, and fails only if every covariate is constant.
  So one covariate file serves all phenotypes.
- **Missing phenotypes.** `complete.cases` over the phenotype plus covariates
  drops those samples for that phenotype alone.

What does need attention is scale. At 22 chromosomes the call count is
`71P + 92`: `66P` step-2 shards across the three arms, `3P` concats, `P` null
fits, `P` arm comparisons, 88 per-chromosome tasks and 4 collection tasks. That
is 802 calls at P = 10 and **46,242 at P = 650**, the list `--min-cases 100`
actually yields — 58x, in a single Cromwell workflow.

- **Batch the phenotypes across submissions** rather than submitting one
  single run over the whole list. A batch is a `felix_runs.entity.tsv` row
  differing only in
  `selected_phenotypes`; the GRM and all 22 packs call-cache across batches, so
  batching costs nothing but gains a workflow that can be watched, diagnosed and
  restarted. Seven batches of ~95 phenotypes is ~6,800 calls each.
- **`SummarizeFelixResults` is the task to size.** It localizes every
  per-phenotype merged TSV *and* copies each into `results_by_phenotype/`, so
  its disk is twice the whole scan — which is why `disk_gb_multiplier` is 2.0.
  Memory is fine (one phenotype is read at a time and the frame is reassigned),
  but set `FelixGenome.Summarize.preemptible = 0` at scale: a preemption
  restarts a localization measured in hundreds of gigabytes. Measure one
  phenotype's genome-wide TSV from the first run before choosing a batch size.
- **The significance threshold is per test.** `summarize_p_threshold` is 5e-8
  with no correction for the number of phenotypes; at 650 phecodes decide
  whether the headline threshold should be 5e-8/650 and set the column
  accordingly.

Resolve FLARE URIs:

```bash
python3 scripts/resolve_flare_uris.py --from-firecloud --scan-chrom chr22 --grm-chroms chr1 chr22
python3 scripts/resolve_flare_uris.py --from-firecloud --autosomes --grm-chroms chr1 chr22
```

## Task reference (FelixGenome.wdl)

Sixteen tasks. Five run once per run, one runs once per chromosome, the rest
scatter. With C chromosomes and P phenotypes the genome-wide run is P null fits,
C packs, C repeat aggregations, and 3 x C x P step-2 shards — 22 x 10 here, so
220 shards in each of the three arms.

### Guards (fail the run before anything expensive starts)

| Task | What it does |
| --- | --- |
| `CheckChromVcfPairs` | Takes only the array *lengths* and accepts one of three input shapes: `joint_vcfs` alone, `phase_vcfs` + `flare_vcfs`, or all three. Anything else — a length mismatch, a half-supplied pair — fails here rather than at the first `Pack`. |
| `CheckNullArrays` | Checks that the null arrays agree in length with the phenotype list, whether they were fitted here or passed in as `precomputed_null_*`. Emits `n_phenotypes_out`, which is the width every downstream scatter is indexed against. Order is the contract: arrays are positional against `read_lines(selected_phenotypes)`. |

### Shared setup (once per run)

| Task | What it does |
| --- | --- |
| `BuildSaigePlinkAndSparseGRM` (`MakeGRM`) | Runs `build_saige_plink_and_grm.sh` and `make_plink_keep.py` over `grm_vcfs` to produce the PLINK trio and the sparse GRM (`.mtx` + sample ids) at `relatedness_cutoff`. `grm_vcfs` is deliberately pinned to chr1 + chr22, not the whole genome, so the GRM and every null built on it stay cache-valid across resubmissions. |
| `FitFelixNull` (`Null`, scatter over phenotypes) | `fit_felix_null.R` runs FELIX step 1 (`fitNULLGLMM`, `--useSparseGRMtoFitNULL=TRUE`), then `export_felix_null.R` packages the fit. Outputs the null `.rda`, the variance ratio, the sample list actually used, and an export tar. A null depends only on the GRM, the covariates and the phenotype column, so it is shared across all chromosomes and all three arms. Skipped entirely when `precomputed_null_rdas` is supplied. |

### Per chromosome, shared by both lanes

| Task | What it does |
| --- | --- |
| `PackFelixla` (`Pack`) | Normalises the VCF header (`bcftools view -h` / `reheader`), indexes it, then runs `felixla --phase-vcf --flare-vcf --n-ancestries --keep --make-felixla`. **This is where local ancestry meets the phased genotypes:** felixla assigns each phase record to its covering FLARE marker itself, which is why the ordinary scan never needs a pre-propagated VCF. Output is one tarball of the packed binary format. |

### Ordinary scan (SNVs, indels and SVs — one test per variant record)

| Task | What it does |
| --- | --- |
| `RunFelixStep2` (`Step2`, scatter chrom x phenotype) | `run_felix_step2.R --felixla-prefix` against the packed tar and that phenotype's null, with `--min-mac`, `--n-ancestries` and `--pvalcutoff-of-haplotype`. Writes the per-shard results TSV. |
| `ConcatPhenotypeScores` (`Concat`) | Merges the C per-chromosome shards for one phenotype into a single TSV. |
| `SummarizeFelixResults` (`Summarize`, `arm = "main"`) | `summarize_felix_results.py` over the merged TSVs: calibration (lambdaGC on `P_cct_admixed_c`), hit table at `summarize_p_threshold`, `qq_cct.png` and `manhattan_cct.png`, bundled as `main.summary.tar.gz` — manifests and plots only, not the per-phenotype TSVs. |

### Repeat-dosage lane (`ru_engine = "felix"` — one test per catalog locus)

| Task | What it does |
| --- | --- |
| `AggregateRepeatLoci` (`RuLoci`) | The heart of the method. Indexes the source VCF, then `aggregate_repeat_loci.py --catalog-bed --flank --max-bp --with-ancestry` walks the TRExplorer catalog and, for each locus, sums the signed REF-relative length change of every overlapping record on each haplotype. `--with-ancestry` is why this lane — and only this lane — needs the ancestry-propagated VCF: it reads per-record `AN1`/`AN2` to attribute each haplotype's length to an ancestry. Writes `loci.vcf`, a per-allele table, and `loci.summary.json` (the QC numbers the review notebook thresholds). |
| `WriteRuAdmixedVcf` (`RuVcf`, `coding = "ru"`) | `write_admixed_dosage_vcf.py` turns the locus VCF into the per-ancestry dosage VCF FELIX step 2 consumes, applying the missing-haplotype policy (`--missing`, `--max-hap-missing`). `ru.admixed.stats.json` records records in/written and why any were dropped. |
| `WriteRuAdmixedVcf` (`BaVcf`, `coding = "carrier"`) | The comparator arm. Same loci, same null, same carriers — only the coding differs: 0/1 carrier instead of repeat dosage. This is what makes the dosage-vs-biallelic claim a controlled comparison rather than two separate analyses. |
| `RunFelixStep2RuVcf` (`RuStep2` / `BaStep2`) | `run_felix_step2.R --vcf-file` instead of `--felixla-prefix`, plus `--dosage-qc`. Thresholds on `ru_min_carriers` rather than `min_mac`, since a locus's signal is carriers of a length change, not minor-allele count. |
| `ConcatPhenotypeScores` (`RuConcat`, `BaConcat`) | As above, with `.ru` / `.ru_biallelic` suffixes on the phenotype name. |
| `CompareRuArms` (`ArmCompare`) | `compare_ru_arms.py` joins the two arms per phenotype at `ru_arm_p_threshold` and writes the 2x2 table, a summary and a markdown block — the direct input to Figure 3's dosage-vs-biallelic panel. |
| `SummarizeFelixResults` (`SummarizeRu`, `SummarizeBa`) | The same summariser on each repeat arm. Identical joint-P columns, one row per locus instead of per variant, so repeat lambdaGC is directly comparable with the main arm's **from the same null** — which is what makes the calibration gate meaningful. |

### Legacy record-level lane (`ru_engine = "rust"` — not used in this run)

`AnnotateRuTest` → `ExtractRuDosage` → `ScoreRuDosage` test each VCF record
separately against a simple-repeat BED, in three encodings (dosage, collapse,
split). Retained only for the encoding comparison. It drops records whose
alleles are not a whole number of repeat units, which is why the locus-level
path above is the default (`felix/REPEAT_DOSAGE.md` section 6.2).

### Collection

| Task | What it does |
| --- | --- |
| `BundleRunSummary` (`Bundle`) | Merges the per-arm summary tars, the per-chromosome `loci.summary.json` and `ru.admixed.stats.json` (renamed by chromosome), and the ArmCompare outputs into one `felix_run.summary.tar.gz`. This is the single small file to download for `notebooks/terra/felix_02_genome_review.ipynb`; the per-phenotype TSVs stay in the bucket. |


## Notes

- Joint p-value for calibration: `P_cct_admixed_c` (FELIX) vs `P` (Tractor-Mix).
- `fit_null.R --step1-rda` delegates to `felix/scripts/export_felix_null.R` for local dev.
- `export_felix_null.R` normally runs from the copy baked into the image
  (`/opt/felix_scripts`), because Cromwell localizes only declared inputs and the
  bucket copy beside `fit_felix_null.R` is never mounted. To run a staged copy instead
  (script fixes without an image rebuild), set the optional `export_null_script` input
  of `FelixPilot` / `FelixGenome` to `gs://BUCKET/felix/scripts/export_felix_null.R`.
  `felix-pilot:0.2.0` needs this: its copy fails under Matrix 1.7.5 on
  `as(<dsTMatrix>, "dgCMatrix")`.
- `run_felix_step2.R` does not pass the sparse GRM to `step2_SPAtests.R` by default
  (`--use-sparse-grm false`). Step 1 fits the null with `--useSparseGRMtoFitNULL=TRUE`, which in
  FELIX's `fitNULLGLMM` resets `useSparseGRMforVarRatio`, so `varianceRatio.txt` has only `null`
  rows and a step 2 given `--sparseGRMFile` stops in `Get_Variance_Ratio`. The WDLs still
  localize and pass `--sparse-grm`; the wrapper ignores them unless `--use-sparse-grm true`,
  which checks for a `sparse` row first. Tractor-Mix's SAIGE runs use the sparse GRM and its
  sparse variance ratio at step 2, so the two arms of the M1 gate differ here; this was accepted
  on 2026-10-09 in favour of following FELIX's guidance (see the design table in `PLAN.md`).
- `wdl/FelixPilotStep2Resume.wdl` is a temporary workflow that restarts at Step2 from the
  `call-MakeGRM`, `call-Pack` and `call-Null` outputs of an earlier FelixPilot run, so a Step2 fix
  can be tested without the GRM build. Its three downstream tasks are copied verbatim from
  `FelixPilot.wdl`; build its inputs with `scripts/make_resume_inputs.py` (`--limit N` for a quick
  test). Delete it once a full run is cached.
- The two WDLs summarize with different scripts. FelixPilot's Summarize task passes `--p-column`, `--named-suffix`
  and `--report-title` and declares `calibration_summary.md`, `lambda_gc_wide.tsv`, `phewas_genomewide_hits.tsv` and
  `qc/*/qq_joint_acpass.png` as outputs: the CLI and output set of `scripts/summarize_tractor_genome_results.py`, so
  its configs use that as `summarize_script` (both arms of the M1 gate then share the λGC code). FelixGenome's
  Summarize task makes the six-argument call and declares `qq_cct.png` / `manhattan_cct.png`, which is
  `felix/scripts/summarize_felix_results.py`. Pointing FelixPilot at the FELIX script failed the 2026-10-09 run with
  `unrecognized arguments`; pointing FelixGenome at the Tractor script would default `--p-column` to `P`.
- Do not headline encoding comparisons without matched λGC (`felix/eval/README.md`).
