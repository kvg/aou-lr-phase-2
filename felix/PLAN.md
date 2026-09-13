# FELIX association testing — plans and status

Coordination doc for humans and agents working on LAI-informed GWAS in
`aou-lr-phase-2`. Operational how-to lives in [`README.md`](README.md);
evaluation gates in [`eval/README.md`](eval/README.md).

**Do not fork FELIX.** Consume `lhu1/felix:latest` via our `felix-pilot` image.
**Do not put FELIX code back under `tractor_mix/`.** Shared Rust crates stay in
`tractor_mix/` because Tractor-Mix / SAIGE still use them.

---

## Goal

Replace Tractor-Mix as the primary LAI association engine for:

1. **SNV / indel / unique-sequence SVs** — FELIXla pack + FELIX Step 2
   (`--is_admixed=TRUE`), joint p-value `P_cct_admixed_c`.
2. **Repeat-mediated SVs** — annotate `RU_TEST` in the joint VCF →
   `extract-tracts-flare` copy-number × ancestry dosages →
   `tractor-mix-score --mode felix` on the **same** FELIX null export.

Same cohort / covariates / GRM markers as the existing Tractor-Mix + SAIGE 2×2
calibration under `tractor_mix_pilot/`.

```text
Joint phased VCF (GT + AN1/AN2 + SVs)
        │
        ├──► felixla pack ──► FELIX Step 2 (0/1 path)
        │         ▲
        │         │
        │    FELIX Step 1 null ──► also exported as CSC for Rust scorer
        │         │
        └──► RU annotator ──► extract-tracts-flare ──► tractor-mix-score --mode felix
             (RU_TEST only)      (dosage / collapse / split)
```

---

## Status (2026-09)

| Milestone | Code | Terra / data gate |
|-----------|------|-------------------|
| **M1** FELIX chr22 pilot (null + pack + step2 + summarize) | Done | **Not run yet** — needs limited/full 2×2 vs Tractor-Mix λGC |
| **M2** FELIX null → CSC + Rust scorer FELIX mode | Done | Local unit tests pass; GMMAT oracle parity on tractor-mix image |
| **M3** RU annotator + extract dosage encodings | Done | Fixture tests pass; not run on AoU joint VCF |
| **M4** Encoding comparison + simulator | Done (scripts) | Simulator smoke OK; genome-wide encoding compare needs Terra |
| **M5** `FelixGenome.wdl` (autosomes + optional RU branch) | Done | Not submitted |

Image: `us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/felix-pilot:0.1.0`
(built FROM `lhu1/felix:latest` + PLINK2 + bcftools + Rust CLIs + staged scripts).

Local smoke summary: [`eval/overnight_local_smoke.log`](eval/overnight_local_smoke.log).

---

## Immediate next work (Terra)

These are the blocking items before claiming M1 / publishing results.

1. **Stage scripts** to the workspace bucket:
   ```bash
   WORKSPACE_BUCKET=gs://... ./scripts/stage_felix_scripts.sh
   WORKSPACE_BUCKET=gs://... ./scripts/stage_tractor_scripts.sh
   ```
2. **Reuse cohort** from `notebooks/terra/tractor_01_prepare_inputs.ipynb`
   (`analysis_samples`, `pheno_cov`, `selected_phenotypes`,
   `covariate_columns_{limited,full}`).
3. **Submit FelixPilot** twice (limited + full) using
   `configs/felix.pilot.inputs.{limited,full}.json.example`.
   Point `gs://BUCKET/...` at the workspace; keep GRM VCFs chr1+chr22.
4. **M1 gate** — compare λGC to existing Tractor-Mix chr22 limited/full:
   ```bash
   python3 felix/scripts/compare_felix_tractor_calibration.py \
     --felix-limited-dir ... --felix-full-dir ... \
     --tractor-limited-dir ... --tractor-full-dir ... \
     --out-dir eval/chr22_felix_tractor_gate
   ```
   Pass: λGC ≈ 1 on null phenotypes; `comparability_flags.tsv` status `OK`.
5. Only after M1: genome-wide `FelixGenome.wdl`, then optional RU_TEST branch.

---

## Design decisions (do not reopen casually)

| Decision | Choice | Why |
|----------|--------|-----|
| FELIX source | Published image only | No compile / fork; ZhouLabGenetics stays upstream |
| GRM | SAIGE mtx via `build_saige_plink_and_grm.sh` | Same as SaigePilot; FELIX Step 1 with `--useSparseGRMtoFitNULL=TRUE` |
| Primary p-column | `P_cct_admixed_c` | FELIX admixed Cauchy combination (conditional columns mirror unconditional for now) |
| Unique SVs | FELIXla 0/1 split-biallelic | Same path as SNV/indel |
| Repeat SVs | Separate dosage scorer | Length-additive biology; compare vs collapse/split |
| Scorer | Extend `tractor-mix-score`, not SAIGE C++ | f64 dosages, hom/het/CCT, quantitative; SPA out of scope for v1 |
| Layout | Top-level `felix/` | Separated from `tractor_mix/` before any Terra runs |
| Scripts on GCS | `gs://BUCKET/felix/scripts/` for FELIX; `gs://BUCKET/scripts/` for shared | Stage via `stage_felix_scripts.sh` / `stage_tractor_scripts.sh` |

### Out of scope (v1)

- FELIXla native copy-number format
- TRACTOR dosage-VCF through FELIX Step 2
- New LAI / re-running FLARE
- Rust SPA (normal approx only; SPA-retest hits later if λGC inflates)
- Forking ZhouLabGenetics/FELIX
- Backward-compat symlinks under `tractor_mix/`

### Related work (RU dosage motivation)

- **Margoliash, Gymrek et al., *Cell Genomics* 2023** —
  [doi:10.1016/j.xgen.2023.100458](https://doi.org/10.1016/j.xgen.2023.100458)
  ([PMC10726533](https://pmc.ncbi.nlm.nih.gov/articles/PMC10726533/)):
  imputed ~446k STRs into UK Biobank; length-based association with blood/serum
  traits; STRs estimated as **5.2–7.6%** of identifiable causal variants; many
  multi-allelic repeats are imperfectly tagged by biallelic SNPs. Supports our
  M3/M4 premise that **copy-number / length dosage** beats collapse or
  split-allele 0/1 encodings, and that GWAS should not stop at SNPs/indels.
  Complementary to us: they use array-imputed STRs in Europeans; we test
  LR-resolved repeat-mediated SVs with **ancestry-partitioned** dosages in AoU.
  Candidate loci (e.g. APOB CTG, CBL CGG, TAOK1 polyA) and nearby blood traits
  are useful positive-control / phenotype ideas alongside Mendelian expansion
  loci (HTT, FMR1, FXN, C9ORF72, RFC1) when present in the callset.

---

## Ownership map

| Path | Role |
|------|------|
| `felix/wdl/` | Terra workflows (`FelixPilot`, `FelixGenome`) |
| `felix/scripts/` | FELIX-only CLIs (null, step2, summarize, calibration, sim) |
| `felix/docker/` | `felix-pilot` image definition |
| `felix/configs/` | Example Cromwell / Terra inputs |
| `scripts/annotate_repeat_units.py` | Shared RU annotator (also staged into image) |
| `scripts/build_saige_plink_and_grm.sh` | Shared GRM builder |
| `tractor_mix/extract_tracts_flare/` | Rust FLARE extract + RU dosage encodings |
| `tractor_mix/tractor_mix_score/` | Rust scorer (`--mode auto\|legacy\|felix`) |
| `tractor_mix/wdl/` | Tractor-Mix + SAIGE only (do not add FELIX here) |

When changing shared Rust crates, keep **legacy Tractor-Mix** (`--mode legacy` /
GMMAT export) green; FELIX mode is additive.

---

## Workflow shapes

### FelixPilot (chr22)

`Check` → `MakeGRM` → scatter `FitFelixNull` → `PackFelixla` → scatter
`RunFelixStep2` → concat → `SummarizeFelixResults`.

- Pack: `felixla --phase-vcf --flare-vcf --keep --make-felixla`
  (joint `GT:AN1:AN2` VCF may fill both VCF flags).
- Step2: `nThreads=1` for full-chrom FELIXla scans (FELIX only parallelizes with
  `--idstoIncludeFile`).
- `--chrom` must match VCF contig (`chr22`, not `22`).

### FelixGenome (autosomes)

Shared null once; scatter pack/step2 over chroms. Optional RU branch when
`simple_repeat_bed` + `annotate_repeat_units_script` are set:
annotate → extract (dosage + collapse + split) → score × encodings.

**Known WDL note:** miniwdl reports `RuScore` outputs unused at the workflow
level — genome summarize currently focuses on FELIXla concat; wiring RU score
TSVs into summarize / encoding-compare is a follow-up when the RU branch is
turned on.

---

## Calibration and p-value conventions

| Model | Joint p-column | Result suffix |
|-------|----------------|---------------|
| FELIX | `P_cct_admixed_c` | `.felix.tsv` |
| Tractor-Mix | `P` | `.tractor_mix.tsv` |

Covariate arms (matched 2×2):

- **Limited:** sex, age, PC1–PC10
- **Full:** limited + standardized coverage + GC dummies

Do **not** headline encoding “more hits” without matched λGC
(`eval/README.md`).

---

## Agent checklist

Before editing:

- [ ] Read this file + `README.md` for the surface you touch
- [ ] Prefer config / script / WDL changes in `felix/` over `tractor_mix/`
- [ ] Do not compile FELIX; bump `felix-pilot` tag if Docker deps change
- [ ] Stage path: FELIX scripts → `felix/scripts/` on GCS, not repo-root `scripts/`

Before claiming M1 done:

- [ ] Limited + full FelixPilot succeed on Terra
- [ ] `compare_felix_tractor_calibration.py` produces `OK` comparability flags
- [ ] λGC ~1 on null phenotypes for all four models

Before claiming M4 / RU genome done:

- [ ] FelixGenome RU branch produces three encoding TSVs per phenotype/chrom
- [ ] `compare_repeat_encodings.py` + simulator type I / power tables reviewed
- [ ] Positive-control loci (HTT, FMR1, FXN, C9ORF72, RFC1) checked if present

---

## Quick commands

```bash
# Image
felix/build_docker.sh

# Stage
WORKSPACE_BUCKET=gs://... ./scripts/stage_felix_scripts.sh
WORKSPACE_BUCKET=gs://... ./scripts/stage_tractor_scripts.sh

# Local tests (no AoU data)
cd tractor_mix/tractor_mix_score && cargo test
cd tractor_mix/extract_tracts_flare && cargo test
python3 -m pytest sv_annotation/tests/test_annotate_repeat_units.py -q

# FLARE URIs
python3 scripts/resolve_flare_uris.py --from-firecloud --scan-chrom chr22 --grm-chroms chr1 chr22
```
