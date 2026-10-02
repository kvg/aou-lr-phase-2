# Design note: ancestry-aware association for repeat-mediated SVs

Status: **proposal; Option A spike run locally 2026-10-02 (§3.1)**. Owner: Kiran Garimella.
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

### Option A — FELIX's VCF dosage input

FELIX step 2 has a second admixed input path besides FELIXla: a VCF with
per-ancestry FORMAT fields (`genoType == "vcf"` in `SAIGE.Admixed`,
`R/SAIGE_SPATest_Tractor.R`). From FELIX `src/Main.cpp`, it reads:

- `DS{k}`: the ancestry-k dosage (FORMAT, Float)
- `ANC{k}`: the number of ancestry-k haplotypes per sample (0/1/2), used as the
  per-ancestry denominator

This is the Tractor dosage-VCF layout (FELIX's own fixture:
`tools/felixla/testdata/tiny.tractor_dosage.vcf`). The path is not documented
in the FELIX README, so we tested it (§3.1).

Scaling: score tests do not change under a linear rescaling of the dosage
(covariates include an intercept). Write `DS{k} = Σ_{h: a_h=k} x_h / s` with one
scale `s` per locus (for example `max |x_h|`), so `|DS{k}| ≤ ANC{k} ≤ 2`.
Report effect sizes per repeat unit by multiplying `BETA` back by `1/s`.

**Never shift the coding** (for example "shortest allele = 0"). A per-haplotype
shift `c` adds `c × ANC{k}` to `DS{k}`, which puts the local-ancestry offset
back into every per-ancestry test (the problem in §1). Only the shared-effect
test is unaffected.

#### 3.1 Spike results (local, felix-pilot:0.1.0, 2026-10-02)

Synthetic cohort: 2,000 samples, 2 ancestries, binary traits at 10% and 2%
prevalence and a quantitative trait, null fitted once per trait (FELIX step 1).
Scripts and outputs: `scratch/m6spike/` (git-ignored).

| Check | Result |
|---|---|
| A-1: 300 SNVs (incl. rare) through FELIXla vs the dosage VCF | **Identical** p-values in every test column, all three traits |
| A-2: repeat dosages at scale `s` vs `2s` | Identical p-values (≤ 8×10⁻⁵ in −log₁₀ p, print precision); `BETA` doubles as expected |
| Negative and non-integer `DS{k}` | Read correctly |
| SPA on continuous dosages | Applied (`Is.SPA = TRUE` where the score is large), same rule as SNVs |
| Shifted coding vs REF-relative | Per-ancestry, het and CCT p-values change by up to 1.2 in −log₁₀ p (confirms "never shift") |
| **Allele-count filter** | FELIX computes `AC = Σ DS{k}` and `MAC = min(AC, Σ ANC{k} − AC)` and requires `MAC ≥ minMAC`, with `minMAC ≥ 0.5` enforced. An ancestry whose summed repeat dosage is negative (contraction-heavy) or near its haplotype count is **dropped**. Here 6 of 40 loci lost an ancestry test, and the shared-effect p-value at those loci changed as well (it is computed over the kept ancestries). |
| Null calibration, 2,000 repeat loci, binary 2% and 10% prevalence and quantitative | Calibrated for every test that ran: λ 0.95–1.09; type I error 0.047–0.054 at 0.05, 0.007–0.010 at 0.01, ≤ 0.0022 at 0.001. SPA applied at about 4% of loci for binary traits. |
| Dropped tests in that run | 55 of 2,000 loci gave no output; 270 / 278 per-ancestry tests dropped, 254 of the 270 because the summed dosage was negative. Contractions are common in this simulation, so real rates may be lower. |

**Conclusion:** FELIX's dosage-VCF path reproduces its FELIXla tests exactly and
gives repeat dosages the same SPA, but its allele-count filter assumes 0–2
genotype dosages. Without a change, contraction-heavy loci silently lose
ancestry tests.

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

### Recommendation after the spike

**Option A with a one-function FELIX patch.** In the admixed VCF path, replace
the allele-count filter for dosage loci with a carrier-count filter (number of
ancestry-k haplotypes with `x_h ≠ 0`), and keep everything else (SPA,
variance rescaling, CCT) unchanged. Because the SNV equivalence is exact, the
statistics are FELIX's own; the patch only stops valid loci being dropped. It
is easy to explain ("FELIX's dosage input, with the genotype-count filter
replaced by a carrier-count filter for repeat loci").

- Until the patch exists, Option A can run unpatched with each locus oriented so
  ancestry sums are positive where possible (multiply the locus by −1: p-values
  unchanged, `BETA` flips). This helps little: in the simulation only 51 of 584
  affected loci had all ancestry sums of one sign. Report every locus where an
  ancestry was dropped and treat those results as incomplete.
- Option B (SPA in the Rust scorer) remains the fallback if the FELIX patch
  cannot be maintained.
- Option C (float dosages in FELIXla) is not needed.

## 4. Data flow

```text
joint phased VCF (HiPhase → SHAPEIT4; SNVs + SVs on the same haplotypes)
   │
   ├─ annotate_repeat_units.py  → INFO RU_TEST, RU (units per ALT), CN_REF, motif
   ├─ propagate_flare_ancestry  → AN1/AN2 at every SV site (FLARE2 production recipe)
   │
   └─ write_admixed_dosage_vcf  (new; Rust extract or Python)
         per sample: DS1..K, ANC1..K for RU_TEST loci (REF-relative, scaled, oriented)
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
| A-1 equivalence | Write SNVs (0/1) as an admixed dosage VCF and run FELIX's VCF path vs the FELIXla path on the same sites | identical p-values (passed on synthetic data; repeat on real chr22) |
| A-2 scaling | Same repeat loci with two different scales `s` | identical p-values (passed on synthetic data) |
| Null calibration | λGC per encoding on null and low-prevalence binary phenotypes (`compare_repeat_encodings.py`), compared with the SNV λGC from the same null | λGC ≈ 1; no excess at low case counts |
| Simulation | `simulate_repeat_dosage.py` architectures (length-additive, single allele, threshold) | type I error at α; dosage beats split/collapse under length-additive |
| Admixture-mapping check | primary vs ancestry-centered encoding at top loci | disagreements listed, not hidden |
| Positive controls | Expansion loci and published STR–trait pairs (APOB CTG, CBL CGG, TAOK1 polyA; Margoliash 2023) | expected direction where power allows |

Do not headline "more hits" for any encoding without matched λGC
(`eval/README.md`).

## 7. Work items

1. ~~Fix the `x_h` baseline in `extract-tracts-flare`~~ Done: `--ru-baseline ref`
   is the default (REF = 0); `absolute` keeps the old coding. The
   ancestry-centered encoding is still to do (in the dosage-VCF writer).
2. ~~Option A spike~~ Done locally on synthetic data (§3.1).
3. Patch FELIX's allele-count filter for dosage loci (small fork of the
   admixed VCF path) and rerun the spike to confirm no loci are dropped and
   SNV results are unchanged.
4. Write the admixed dosage-VCF writer (`DS{k}`, `ANC{k}` from the joint
   `GT:AN1:AN2` VCF; REF-relative, scaled, oriented) and wire it into
   `FelixGenome.wdl` in place of the `tractor-mix-score` RU branch.
5. Run the upstream checks in section 5 on chr22, then genome-wide.

## 8. Open questions

- Which repeat catalog defines `RU_TEST` loci for the callset (TRGT catalog vs
  `simpleRepeat`)? This affects `CN_REF` and the motif definitions.
- Should segmental-duplication-mediated CNVs (gene-family copy number) use the
  same length-dosage path? They are copy-number dosages too, but their phasing
  and ancestry assignment are less reliable.
- How many ancestries does the FLARE2 production model have (`nanc`)? That sets
  `K` for every test (`num_ancs` in the FELIX WDLs).
