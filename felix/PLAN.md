# FELIX association testing — plans and status

Coordination doc for humans and agents working on LAI-informed GWAS in
`aou-lr-phase-2`. Operational how-to lives in [`README.md`](README.md);
evaluation gates in [`eval/README.md`](eval/README.md).

Consume `lhu1/felix:latest` via our `felix-pilot` image. A FELIX fork is
**allowed if needed** (decided 2026-10-02) but is the last option for repeat
SVs; see [`SV_SCORER_DESIGN.md`](SV_SCORER_DESIGN.md).
**Do not put FELIX code back under `tractor_mix/`.** Shared Rust crates stay in
`tractor_mix/` because Tractor-Mix / SAIGE still use them.

## Decisions (2026-10-02)

| Decision | Choice |
|---|---|
| LAI engine | **FLARE2** (clustered model trained once, applied to every shard and chromosome) unless selection v2 shows the original-FLARE pin is better |
| Association engine | **FELIX** (replaces Tractor-Mix) |
| Cohort | **Whole cohort** for LAI-informed tests (all populations, including MID) |
| Phasing | SNVs + SVs physically phased with HiPhase, then statistically phased with SHAPEIT4 (same haplotypes) |
| Repeat SVs | Locus-level length dosage × ancestry through FELIX's own tests with SPA, via FELIX's admixed dosage-VCF input plus our carrier-QC patch (`felix/patches/0001-dosage-carrier-qc.patch`); SPA in our Rust scorer is the fallback. **All repeat-dosage detail now lives in [`REPEAT_DOSAGE.md`](REPEAT_DOSAGE.md)**; patch-not-fork is settled and the FELIX `repeat-sv-dosage` branch is shelved. Upstreaming deferred until the pipeline works end to end |
| Recipe selection | Time-boxed selection v2 on chr20 (`flare/README.md` Part 8); no Tractor null-λ |

---

## Goal

Replace Tractor-Mix as the primary LAI association engine for:

1. **SNV / indel / unique-sequence SVs** — FELIXla pack + FELIX Step 2
   (`--is_admixed=TRUE`), joint p-value `P_cct_admixed_c`.
2. **Repeat-mediated SVs** — annotate `RU_TEST` in the joint VCF →
   repeat-length deviation × ancestry dosages → FELIX step 2 (SPA) on the
   **same** FELIX null, through FELIX's admixed dosage-VCF input. Design, spike
   results and engine choice: [`SV_SCORER_DESIGN.md`](SV_SCORER_DESIGN.md). The
   current WDL path (`extract-tracts-flare` → `tractor-mix-score --mode felix`)
   has no SPA; treat its outputs as provisional.

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

## Status (2026-10)

| Milestone | Code | Terra / data gate |
|-----------|------|-------------------|
| **M1** FELIX chr22 pilot (null + pack + step2 + summarize) | Done | **Not run yet** — needs limited/full 2×2 vs Tractor-Mix λGC |
| **M2** FELIX null → CSC + Rust scorer FELIX mode | Done | Local unit tests pass; GMMAT oracle parity on tractor-mix image |
| **M3** RU annotator + extract dosage encodings | Done | Fixture tests pass; not run on AoU joint VCF |
| **M4** Encoding comparison + simulator | Done (scripts) | Simulator smoke OK; genome-wide encoding compare needs Terra |
| **M5** `FelixGenome.wdl` (autosomes + optional RU branch) | Done | Not submitted |
| **M0** LAI recipe (selection v2: FLARE2 vs FLARE pin, chr20) | FLARE2 WDL mode + scorers done | Build `aou-flare2` image; run `sel_chr20_*` rows; flare_02 Part 8 |
| **M6** Repeat-SV scorer with SPA | Done: REF-relative dosage; dosage-VCF writer; FELIX carrier-QC patch (`felix/patches/`); `FelixGenome` RU branch on FELIX step 2. Validated locally on synthetic data; `felix-pilot:0.2.0` pushed | Run the RU branch on chr22 `joint_vcfs` |
| **M6b** Locus-level repeat dosage | Done: RU branch rewired from the record-level annotator to `scripts/aggregate_repeat_loci.py --out-vcf` (one test per catalog locus, signed by allele length, impure alleles kept); missing-genotype policy and multi-digit ancestry in the dosage writer. See [`REPEAT_DOSAGE.md`](REPEAT_DOSAGE.md) | Same chr22 gate as M6 |

Image: `us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/felix-pilot:0.2.0`
(built FROM FELIX v0.1 pinned by digest + `felix/patches/` (carrier QC, off by default) + PLINK2 +
bcftools + Rust CLIs + staged scripts). 0.2.0 pushed 2026-10-02 from commit `f65852f`
(Cloud Build `0219507e-71cd-44ae-b28a-d154059f0bcf`), digest
`sha256:a9b51dad68b6a8881357cbb9ce2317a68a5d00eadf6c6bc65b43612935c5868b`.

Local smoke summary: [`eval/overnight_local_smoke.log`](eval/overnight_local_smoke.log).

---

## Immediate next work (Terra)

These are the blocking items before claiming M1 / publishing results. M1 does
not wait for M0: it compares FELIX with Tractor-Mix on the LAI the Tractor-Mix
pilot already used. Both pilot configs now scan the same FLARE-marker VCF as
that pilot, so the 2×2 compares like with like.

1. **Stage scripts**: run `notebooks/terra/00_sync_repo.ipynb`. It stages
   `scripts/` → `gs://BUCKET/scripts/` and `felix/scripts/` → `gs://BUCKET/felix/scripts/`
   and copies the WDLs to the bucket; register changed WDLs yourself in the
   Terra Methods Repository UI.
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
   Notebook: `notebooks/terra/felix_01_pilot_gate.ipynb`.
5. M1b: rerun FELIX on the full long-read callset (`phase_vcf` + `flare_vcf`).
6. After M0 + M1: genome-wide FLARE2 apply (`flare/configs/flare2.apply.inputs.json.example`),
   then `FelixGenome.wdl`. Set `num_ancs` to the FLARE2 model's `nanc`.
7. M6: build and push `felix-pilot:0.2.0` (`felix/build_docker.sh`; carries the FELIX
   patch), then run the RU branch with `joint_vcfs` on chr22
   ([`SV_SCORER_DESIGN.md`](SV_SCORER_DESIGN.md) §7).

---

## Design decisions (do not reopen casually)

| Decision | Choice | Why |
|----------|--------|-----|
| FELIX source | Published image; fork allowed if needed | Prefer FELIX's own inputs (FELIXla, admixed dosage VCF) over code changes; FELIX is GPL-3.0 |
| GRM | SAIGE mtx via `build_saige_plink_and_grm.sh` | Same as SaigePilot; FELIX Step 1 with `--useSparseGRMtoFitNULL=TRUE` |
| Primary p-column | `P_cct_admixed_c` | FELIX admixed Cauchy combination (conditional columns mirror unconditional for now) |
| Unique SVs | FELIXla 0/1 split-biallelic | Same path as SNV/indel |
| Repeat SVs | Length dosage through FELIX tests with SPA | Length-additive biology; compare vs collapse/split; see `SV_SCORER_DESIGN.md` |
| Scorer | FELIX VCF dosage input first; Rust scorer + SPA as fallback | SPA is required for low-prevalence binary traits across the phenome |
| Layout | Top-level `felix/` | Separated from `tractor_mix/` before any Terra runs |
| Scripts on GCS | `gs://BUCKET/felix/scripts/` for FELIX; `gs://BUCKET/scripts/` for shared | Staged by `00_sync_repo` (`stage_*_scripts.sh` as a fallback) |

### Out of scope (v1)

- FELIXla native copy-number format (only if options A and B in `SV_SCORER_DESIGN.md` fail)
- Backward-compat symlinks under `tractor_mix/`

(Re-running LAI, admixed dosage-VCF input to FELIX step 2, and SPA for repeat
SVs moved into scope on 2026-10-02.)

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
- [ ] Prefer FELIX's existing inputs over a fork; bump `felix-pilot` tag if Docker deps change
- [ ] Stage path: FELIX scripts → `felix/scripts/` on GCS, not repo-root `scripts/`

Before claiming M1 done:

- [ ] Limited + full FelixPilot succeed on Terra
- [ ] `compare_felix_tractor_calibration.py` produces `OK` comparability flags
- [ ] λGC ~1 on null phenotypes for all four models

Before claiming M4 / RU genome done:

- [ ] Repeat dosage coded relative to `CN_REF` and scored with SPA (M6)
- [ ] FelixGenome RU branch produces three encoding TSVs per phenotype/chrom
- [ ] `compare_repeat_encodings.py` + simulator type I / power tables reviewed
- [ ] Positive-control loci (HTT, FMR1, FXN, C9ORF72, RFC1) checked if present

---

## Quick commands

```bash
# Image
felix/build_docker.sh

# Stage: run notebooks/terra/00_sync_repo.ipynb on the Terra VM, or from a shell with bucket access:
WORKSPACE_BUCKET=gs://... ./scripts/stage_felix_scripts.sh
WORKSPACE_BUCKET=gs://... ./scripts/stage_tractor_scripts.sh

# Local tests (no AoU data)
cd tractor_mix/tractor_mix_score && cargo test
cd tractor_mix/extract_tracts_flare && cargo test
python3 -m pytest sv_annotation/tests/test_annotate_repeat_units.py -q

# FLARE URIs
python3 scripts/resolve_flare_uris.py --from-firecloud --scan-chrom chr22 --grm-chroms chr1 chr22
```
