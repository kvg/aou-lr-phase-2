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

Resolve FLARE URIs:

```bash
python3 scripts/resolve_flare_uris.py --from-firecloud --scan-chrom chr22 --grm-chroms chr1 chr22
python3 scripts/resolve_flare_uris.py --from-firecloud --autosomes --grm-chroms chr1 chr22
```

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
