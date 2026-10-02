# Design note: ancestry-aware association for repeat-mediated SVs

Status: **proposal, not implemented** (2026-10-02). Owner: Kiran Garimella.
Related docs: [`PLAN.md`](PLAN.md) (milestones), [`eval/README.md`](eval/README.md)
(calibration gates).

## 1. Problem

We want local-ancestry-aware association tests across all AoU phenotypes for
every long-read variant class, using FELIX. Three classes:

| Class | Encoding | Engine today |
|---|---|---|
| SNVs, indels | haplotype 0/1 × local ancestry | FELIX (FELIXla pack + step 2) |
| Unique-sequence SVs (presence/absence) | haplotype 0/1 × local ancestry | FELIX, split-biallelic like indels |
| **Repeat-mediated SVs** (STR/VNTR length changes, multi-allelic) | **repeat length per haplotype** × local ancestry | Rust `tractor-mix-score --mode felix` (no SPA) |

FELIXla stores one bit per haplotype per ALT allele, so it cannot hold a
repeat length. Splitting a VNTR into one 0/1 test per allele throws away the
length-additive signal and spreads it across many rare alleles
([Margoliash et al. 2023](https://doi.org/10.1016/j.xgen.2023.100458)).

The current repeat-SV path has two problems:

1. **No saddlepoint approximation (SPA).** The Rust scorer uses the normal
   approximation. For low-prevalence binary traits (most of the phenome) that
   inflates p-values for skewed dosages. FELIX's own tests use SPA.
2. **Per-ancestry dosage is mostly a haplotype count.** `extract-tracts-flare`
   adds the **absolute** copy number `C = CN_REF ± RU` for each haplotype
   (`tractor_mix/extract_tracts_flare/src/ru.rs`). So the ancestry-k dosage is
   `CN_REF × (number of ancestry-k haplotypes) + deviations`. The per-ancestry
   and heterogeneity tests then mostly test local-ancestry dosage (admixture
   mapping), not repeat length. The shared-effect test on the sum is unaffected,
   because the offset `2·CN_REF` is a constant. SNVs avoid this because the REF
   allele is coded 0.

## 2. Model

Keep FELIX's model and tests so repeat results read like the SNV results.

For haplotype `h` of sample `i` at a repeat locus:

- `a_h` = local ancestry (from FLARE2, propagated to the SV position)
- `x_h` = repeat-length deviation from the reference allele, in repeat units:
  `x_h = C_h − CN_REF` (0 on a REF haplotype, like an SNV's REF allele)

Per-ancestry dosage: `D_ik = Σ_{h ∈ i, a_h = k} x_h`. Total: `D_i = Σ_k D_ik`.

Tests (same as FELIX `step2_SPAtests.R --is_admixed=TRUE`):

| Column | Test | Meaning |
|---|---|---|
| `P_hom_admixed` | 1-df score test on `D_i` | shared per-repeat-unit effect |
| per-ancestry `p.value_k` | 1-df score test on `D_ik` | ancestry-specific effect |
| `P_het_admixed` | K-df joint test on `(D_i1 … D_iK)` | effects differ by ancestry |
| `P_cct_admixed` | Cauchy combination (SKAT-O, het, hom) | primary p-value |

Null model: the same FELIX step-1 fit used for SNVs (same samples, covariates,
sparse GRM). For binary traits, the per-ancestry and shared tests use SPA, and
the joint test rescales its variance matrix so that it agrees with the
SPA-corrected single tests. FELIX already does this rescaling
(`get_newPhi_scaleFactor_cpp` in `src/Main.cpp`).

**Encodings** (`extract-tracts-flare` already writes dosage / collapse / split
from the same haplotypes):

- **primary:** length deviation `x_h` above
- **comparison:** collapse (any non-REF = 1) and split (one 0/1 test per ALT)
- **sensitivity:** ancestry-centered length, `x_h − mean_k(x)` over ancestry-k
  haplotypes. This separates the repeat-length effect from admixture mapping;
  flag loci where it disagrees with the primary encoding.
- **expansion loci** (HTT, FMR1, FXN, C9orf72, RFC1): a threshold indicator
  (allele ≥ pathogenic repeat count) alongside the length dosage

## 3. Engine options

### Option A — FELIX's VCF dosage input (recommended first; no fork)

FELIX step 2 has a second admixed input path besides FELIXla: a VCF with
per-ancestry FORMAT fields (`genoType == "vcf"` in `SAIGE.Admixed`,
`R/SAIGE_SPATest_Tractor.R`). From FELIX `src/Main.cpp`, it reads:

- `ANC{k}`: the number of ancestry-k haplotypes per sample (0/1/2), used as the
  per-ancestry denominator
- `DS{k}`: the ancestry-k dosage
- `DSALL`: the total dosage (the shared-effect test)

This is the Tractor dosage-VCF layout. The SAIGE VCF reader takes `DS` as a
float, so a repeat-length dosage can go through the same SPA tests as SNVs with
no FELIX code changes.

Scaling: score tests do not change under a linear rescaling of the dosage
(covariates include an intercept). To stay inside dosage-range checks, write
`DS{k} = Σ_{h: a_h=k} x_h / s` with one scale `s` per locus (for example
`max |x_h|`), so `|DS{k}| ≤ ANC{k} ≤ 2`. Report effect sizes per repeat unit by
multiplying `BETA` back by `1/s`.

Unknowns the spike must answer:

1. Does the VCF admixed path accept negative or non-integer `DS{k}`
   (contractions give `x_h < 0`)? If not, shift so the shortest observed allele
   is 0 (an affine change, so p-values are unchanged).
2. Allele frequency and minor-allele count are computed as `Σ DS{k} / Σ ANC{k}`.
   With scaled lengths, MAC no longer counts carriers. Check how this affects the
   MAC filters and the choice of variance ratio. We may need a carrier-count
   filter outside FELIX.
3. Speed of the VCF reader. Repeat loci are about 10⁴–10⁵ genome-wide, so a
   slower reader is acceptable.

### Option B — add SPA to our Rust scorer

Implement SAIGE's binary SPA (the cumulant generating function of a sum of
Bernoullis with fitted `mu` on the covariate-projected dosage) and the
variance rescaling for the joint test in `tractor-mix-score --mode felix`.
`mu` is already in the null export.

- **Pros:** fully under our control; fast; reuses the existing extract → score
  plumbing.
- **Cons:** an independent reimplementation that needs a parity test against
  FELIX on SNVs, and it diverges whenever FELIX changes.

### Option C — fork FELIXla to store float dosages

Extend the packed format with a per-haplotype float plane for flagged loci.
FELIX is GPL-3.0, so a fork is fine to maintain but harder to explain and to
keep in sync. Choose this only if A and B both fail.

**Decision gate:** run the Option A spike first. If A passes the equivalence
check below, use it. Otherwise build B.

## 4. Data flow

```text
joint phased VCF (HiPhase → SHAPEIT4; SNVs + SVs on the same haplotypes)
   │
   ├─ annotate_repeat_units.py  → INFO RU_TEST, RU (units per ALT), CN_REF, motif
   ├─ propagate_flare_ancestry  → AN1/AN2 at every SV site (FLARE2 production recipe)
   │
   └─ write_admixed_dosage_vcf  (new; Rust extract or Python)
         per sample: ANC1..K, DS1..K, DSALL for RU_TEST loci
         │
         └─ FELIX step2 --vcfFile … --is_admixed=TRUE --number_of_ancestry=K
               same null (step-1 .rda + variance ratio) as the SNV run
```

`FelixGenome.wdl`'s RU branch changes from "extract → `tractor-mix-score`" to
"extract → admixed dosage VCF → FELIX step 2". The summary step then handles
repeat results like FELIXla results (today the `RuScore` outputs are not wired
into the summary at all).

## 5. Upstream checks (before any results)

1. **Ancestry at SV loci.** Measure the distance from each repeat SV to its
   covering FLARE marker (`propagate_flare_ancestry --stats-json`). Repeat SVs
   sit in regions the context mask removes, so the production LAI recipe should
   not use the context mask (`flare/README.md`, selection v2).
2. **SV phasing.** SHAPEIT4 statistically re-phased the HiPhase blocks. Measure
   the switch-error rate at SV sites with the trios (child SV haplotype vs
   parents), separately for repeat SVs.
3. **Repeat genotype quality.** Long VNTR alleles have length-calling noise.
   Record per-locus Mendelian length inconsistency, and drop or flag loci above
   a threshold.

## 6. Validation plan

| Step | Check | Pass |
|---|---|---|
| A-1 equivalence | Write SNVs (0/1) as an admixed dosage VCF and run FELIX's VCF path vs the FELIXla path on the same chr22 sites | identical p-values up to rounding |
| A-2 scaling | Same repeat loci with two different scales `s` | identical p-values |
| Null calibration | λGC per encoding on null and low-prevalence binary phenotypes (`compare_repeat_encodings.py`), compared with the SNV λGC from the same null | λGC ≈ 1; no excess at low case counts |
| Simulation | `simulate_repeat_dosage.py` architectures (length-additive, single allele, threshold) | type I error at α; dosage beats split/collapse under length-additive |
| Admixture-mapping check | primary vs ancestry-centered encoding at top loci | disagreements listed, not hidden |
| Positive controls | Expansion loci and published STR–trait pairs (APOB CTG, CBL CGG, TAOK1 polyA; Margoliash 2023) | expected direction where power allows |

Do not headline "more hits" for any encoding without matched λGC
(`eval/README.md`).

## 7. Work items

1. Fix the `x_h` baseline in `extract-tracts-flare` (code relative to `CN_REF`)
   and add the ancestry-centered encoding. Keep the old absolute mode behind a
   flag for comparison.
2. Option A spike: a small writer for the admixed dosage VCF, then run FELIX
   step 2 on the chr22 repeat loci and run checks A-1 and A-2.
3. Depending on the spike: wire A into `FelixGenome.wdl`, or implement SPA in the
   Rust scorer (B) with a parity test against FELIX on SNVs.
4. Run the upstream checks in section 5 on chr22, then genome-wide.

## 8. Open questions

- Which repeat catalog defines `RU_TEST` loci for the callset (TRGT catalog vs
  `simpleRepeat`)? This affects `CN_REF` and the motif definitions.
- Should segmental-duplication-mediated CNVs (gene-family copy number) use the
  same length-dosage path? They are copy-number dosages too, but their phasing
  and ancestry assignment are less reliable.
- How many ancestries does the FLARE2 production model have (`nanc`)? That sets
  `K` for every test (`num_ancs` in the FELIX WDLs).
