# Repeat-dosage association (single source of truth)

Local-ancestry-informed association for repeat-mediated SVs, where the test
carries the **actual repeat dosage** instead of a biallelic 0/1 indicator.

This note supersedes the repeat-dosage sections of
[`PLAN.md`](PLAN.md) (M6 row), [`SV_SCORER_DESIGN.md`](SV_SCORER_DESIGN.md)
and [`../repeat_loci/README.md`](../repeat_loci/README.md). Those files keep
their history and point here. Calibration gates stay in
[`eval/README.md`](eval/README.md).

Owner: Kiran Garimella. Last substantive change: 2026-10-07.

---

## 1. Settled decisions

| Decision | Choice | Where it lives |
|---|---|---|
| FELIX source | **Patch, do not fork.** FELIX's own admixed dosage-VCF step-2 path plus one opt-in patch | [`patches/0001-dosage-carrier-qc.patch`](patches/0001-dosage-carrier-qc.patch) |
| Tested unit | **One catalog locus**, not one VCF record | `scripts/aggregate_repeat_loci.py` |
| Dosage coding | REF-relative repeat units, `x_h = (haplotype length − reference length) / period` | §2 |
| Shifting the coding | **Never** | §2.1 |
| Per-locus scaling | `s = max|x_h|`, reported as INFO `RU_SCALE`; multiply `BETA` by `1/s` | §2.2 |
| QC filter | Carrier counts, not `sum(DS)` allele counts (`FELIX_DOSAGE_QC=carrier`) | §4 |
| Null model | The same FELIX step-1 fit as the SNV scan — same samples, covariates, sparse GRM | §3 |
| Headline p-value | `P_cct_admixed_c` | §3 |
| Fallback engine | SPA in the Rust `tractor-mix-score` (Option B) if the patch cannot be carried to a newer FELIX | §5 |

Not pursued: forking FELIX to add a continuous-dosage step-2 entry point
(§5.3).

---

## 2. Model

For haplotype `h` of sample `i` at a catalog repeat locus:

- `a_h` — local ancestry at the locus, from FLARE2, propagated onto the site
- `x_h` — repeat-length deviation from the reference allele, in repeat units:

```
x_h = ( Σ_{records r at the locus} [ len(ALT_r[a_h]) − len(REF_r) ] ) / period
```

The sum is what makes this a locus-level test: at one tandem array the
integrated callset holds several records (GLnexus small indels plus one or
more merged SV records), and a haplotype's total length change is the sum of
the signed changes it carries across all of them. A reference haplotype has
`x_h = 0`, exactly as an SNV's REF allele is coded 0.

Per-ancestry dosage and total:

```
D_ik = Σ_{h ∈ i, a_h = k} x_h          D_i = Σ_k D_ik
```

`period` is the motif length of the catalog locus. For a compound locus
(`n_motifs > 1`) the **first** motif's length is used for every record at that
locus, so `RU_DOSAGE` at such loci is "units of the primary motif", not of
whichever motif a given allele changed. Loci with no motif in the catalog get
`period = 1`, i.e. dosage in bp. Both conventions are recorded in the allele
table (`period`, `n_motifs` columns) so they can be filtered on.

### 2.1 Never shift the coding

A per-haplotype shift `c` (for example "shortest observed allele = 0", or
per-ancestry centring) adds `c · ANC{k}` to `DS{k}`. That puts the
local-ancestry haplotype count straight back into every per-ancestry test, so
the per-ancestry and heterogeneity tests become admixture mapping rather than
length association. Only the shared-effect test on `D_i` is invariant, because
its offset is a constant.

This is not hypothetical: it is the defect the original Rust path had.
`extract-tracts-flare` coded **absolute** copy number, so each per-ancestry
dosage was `CN_REF × (number of ancestry-k haplotypes) + deviations`. Fixed by
`--ru-baseline ref`, now the default (`absolute` keeps the old coding). In the
Option A spike, shifted versus REF-relative coding moved per-ancestry, het and
CCT p-values by up to 1.2 in −log₁₀ p.

The same argument rules out the ancestry-centred sensitivity encoding as
originally specified; see §7.

### 2.2 Scaling is free, shifting is not

Score tests are invariant to a linear rescaling of the dosage when the
covariates include an intercept. One scale `s = max|x_h|` per locus keeps
`|DS{k}| ≤ ANC{k} ≤ 2`, which is the range FELIX's readers and its variance
ratio were built for. Report effects per repeat unit by multiplying `BETA` by
`1/s`. Confirmed empirically: the same loci at `s` and `2s` gave identical
p-values (≤ 8×10⁻⁵ in −log₁₀ p, print precision) with `BETA` doubling.

---

## 3. Tests

Identical to FELIX's own admixed single-variant tests, so repeat results read
like the SNV results:

| Column | Test | Meaning |
|---|---|---|
| `P_hom_admixed` | 1-df score test on `D_i` | shared per-repeat-unit effect |
| per-ancestry `p.value_k` | 1-df score test on `D_ik` | ancestry-specific effect |
| `P_het_admixed` | K-df joint test on `(D_i1 … D_iK)` | effects differ by ancestry |
| `P_cct_admixed` / `P_cct_admixed_c` | Cauchy combination | **headline p-value** is the `_c` (haplotype-conditioned) column |

For binary traits the per-ancestry and shared tests use SPA, and the joint test
rescales its variance matrix to agree with the SPA-corrected single tests
(`get_newPhi_scaleFactor_cpp`, FELIX `src/Main.cpp`). Nothing here is
reimplemented; the dosage simply enters FELIX where a 0/1 genotype would.

---

## 4. Why the patch is needed

FELIX computes `AC = Σ DS{k}` and `MAC = min(AC, Σ ANC{k} − AC)`, then
requires `MAC ≥ minMAC` with `minMAC ≥ 0.5` enforced. That treats the summed
dosage as an allele count. REF-relative repeat dosage is not an allele count:
it is negative at a contraction-heavy locus and can exceed the haplotype count
at an expansion. Stock FELIX therefore **silently drops** those tests.

`0001-dosage-carrier-qc.patch` (59 lines in `mainMarkerAdmixedInCPP`) adds an
opt-in switch: with `FELIX_DOSAGE_QC=carrier` and dosage-VCF input, the
`minMAC` / `minMAF` / `MACCutoffforER` checks count **carriers** — samples with
a non-zero `DS{k}`, versus samples without, among samples with an ancestry-k
haplotype (`ANC{k} > 0`; all samples for the combined test). Under REF-relative
coding "non-zero" means "carries a non-reference repeat allele", so `--min-mac`
is a carrier count in this mode. SPA, Firth, the joint tests and the reported
`AC`/`AF` columns are untouched, and with the switch off nothing changes.

`run_felix_step2.R --vcf-file … --dosage-qc carrier` sets it. The image applies
the patch to FELIX's in-image source and reinstalls its R package.

### 4.1 Validation to date — synthetic only

Synthetic cohort: 2,000 samples, 2 ancestries, binary traits at 10% and 2%
prevalence plus a quantitative trait, null fitted once per trait. Scripts and
outputs in `scratch/m6spike/` (git-ignored); rerun harness in
[`eval/dosage_qc/`](eval/dosage_qc/).

| Check | Result |
|---|---|
| Switch off vs stock FELIX (FELIXla SNVs, dosage-VCF SNVs, repeat loci) | byte-identical output files |
| A-1: 300 SNVs through FELIXla vs the dosage VCF | identical p-values in every column, all three traits |
| A-2: repeat loci at scale `s` vs `2s` | identical p-values; `BETA` doubles |
| Negative and non-integer `DS{k}` | read correctly |
| SPA on continuous dosages | applied where the score is large, same rule as SNVs (~4% of loci, binary traits) |
| Stock FELIX, 2,000 null repeat loci | 55 loci produced no output; 270 of 278 per-ancestry tests dropped, 254 of those because the summed dosage was negative |
| Patched, same 2,000 loci (`minMAC` = 20 carriers) | all loci and all per-ancestry tests kept; λ 0.91–1.08; type I error 0.044–0.053 at 0.05, 0.006–0.011 at 0.01, ≤ 0.002 at 0.001 |

λ for `P_cct_admixed` sits a little below 1 (0.91–0.97); the Cauchy combination
of correlated tests is conservative and stock FELIX shows the same on SNVs.

**Nothing above has been run on real data.** The chr22 gate is §6.

### 4.2 Locus-level path, writer-boundary checks (local, 2026-10-07)

`eval/dosage_qc/locus_path_smoke.py` drives the production chain on a synthetic
cohort (500 samples, 3 ancestries, 210 catalog loci; 1–3 records per locus, a
mix of whole-unit and impure alleles, insertions and deletions, ~2% haplotype
dropout) up to the point FELIX step 2 would take over, and re-establishes A-1
and A-2 at the writer boundary — where they can be checked without running
FELIX.

| Check | Result |
|---|---|
| Chain runs end to end; `bcftools` accepts and indexes the dosage VCF | 210 loci in, 210 records out |
| Missing haplotypes under `--missing ref` | 4,942 of 210,000 (2.35%) absorbed; no locus hit the `--max-hap-missing 0.2` cap |
| `\|DS{k}\| ≤ ANC{k} ≤ 2` at every sample and locus | holds (max of `\|DS\| − ANC` = 0) |
| Ancestry inconsistency across records at a locus | 0 haplotypes |
| **A-1 equivalence.** 70 loci holding exactly one record with one whole-unit ALT and no missing calls, locus path vs record-level path | `DS{k}` and `ANC{k}` **identical** (max abs difference 0) |
| **A-2 scale invariance.** Every input dosage multiplied by 7.5 | `DS{k}` **unchanged** (max abs difference 0); `RU_SCALE` scales by 7.5 to within 1.5×10⁻⁷ (the writer's 8-significant-figure output) |

A-2 is stated here as invariance of `DS` rather than "p-values unchanged at two
scales" (§4.1). The writer chooses `s = max|x_h|` per locus, so a constant
rescaling of the input cancels exactly in `DS` and survives only in `RU_SCALE`.
That is the stronger statement and it does not need FELIX to check.

A-1 is restricted to single-record whole-unit loci deliberately: that is the
only configuration in which the two codings should coincide, so a difference
there is a bug in the new path rather than a property of locus aggregation.
Where they are expected to differ, the fixture in `eval/record_vs_locus/`
quantifies it.

**FELIX step 2 was not run locally** for these checks. The image is amd64 and a
sparse-GRM step-1 fit previously ran over 40 minutes under emulation; the
FELIX-side behaviour of the missingness policy and carrier QC is argued from
the patch and FELIX's source (§4, and the module docstring of
`write_admixed_dosage_vcf.py`), not executed. It is executed by the chr22 job.

---

## 5. Engine options as evaluated

### 5.1 Option A — FELIX's dosage-VCF input (chosen)

FELIX step 2 has a second admixed input path besides FELIXla: a VCF with
per-ancestry FORMAT fields (`genoType == "vcf"` in `SAIGE.Admixed`,
`R/SAIGE_SPATest_Tractor.R`), reading `DS{k}` (ancestry-k dosage, Float) and
`ANC{k}` (ancestry-k haplotype count, 0/1/2). This is the Tractor dosage-VCF
layout; FELIX ships a fixture (`tools/felixla/testdata/tiny.tractor_dosage.vcf`)
but does not document the path, so we tested it (§4.1).

FELIXla itself cannot carry a repeat length — it stores one bit per haplotype
per ALT allele, and `tools/felixla/src/tractor_dosage_vcf_to_hybrid.cpp` aborts
on a fractional `DS`.

### 5.2 Option B — SPA in the Rust scorer (fallback)

Implement SAIGE's binary SPA and the joint-test variance rescaling in
`tractor-mix-score --mode felix`; `mu` is already in the null export. Fully
under our control and fast, but an independent reimplementation needing a
parity test against FELIX on SNVs. Kept as the fallback if the patch cannot be
carried forward to a newer FELIX.

### 5.3 Option C — fork FELIX (not pursued)

The branch `repeat-sv-dosage` in the FELIX checkout holds skeletons for a
continuous-dosage step-2 entry point: `SAIGE.RepeatSV()` (`R/repeat_sv.R`) and
`scoreTestAdmixedDosage()` (`src/repeat_sv_shim.cpp`), both marked
`STATUS: skeleton`, plus `docs/repeat-sv-dosage.md`. Shelved, for three
reasons:

1. It is blocked on refactoring the per-marker body of
   `mainMarkerAdmixedInCPP` (~400 lines) into a reusable
   `scoreOneAdmixedMarker`, or duplicating it and accepting drift.
2. Option A reproduces FELIX's own tests exactly (§4.1, A-1), so the fork buys
   no statistical capability — only a different input format.
3. Its R driver specifies **z-scoring `G` per ancestry**, which is a shift and
   therefore violates §2.1: it re-injects the local-ancestry offset into every
   per-ancestry test. If the fork is ever revived, that must become a
   scale-only transform.

Forking also means maintaining a GPL-3.0 fork of an actively changing package
for the lifetime of the manuscript.

---

## 6. Data flow

```text
joint phased VCF  (HiPhase → SHAPEIT4; SNVs + indels + SVs on the same haplotypes)
   │              FORMAT GT:AN1:AN2  (FLARE2 ancestry propagated onto every site)
   │
   ├─ scripts/aggregate_repeat_loci.py --out-vcf --with-ancestry
   │     one record per polymorphic catalog locus
   │     ALT <TR:+Nbp>, INFO RU_TEST / RU_DOSAGE (signed units per ALT) / PERIOD / MOTIF
   │     FORMAT GT:AN1:AN2        ← haplotype dosages summed across the locus's records
   │
   └─ scripts/write_admixed_dosage_vcf.py --num-ancs K
         FORMAT DS1..DSK:ANC1..ANCK, INFO RU_SCALE
         │
         └─ FELIX step 2  --vcfFile … --vcfField=DS --is_admixed=TRUE --number_of_ancestry=K
               FELIX_DOSAGE_QC=carrier; the same step-1 null and variance ratio as the SNV scan
```

In `FelixGenome.wdl` with `ru_engine = "felix"`: `AggregateRepeatLoci` →
`WriteRuAdmixedVcf` → `RunFelixStep2RuVcf` per phenotype →
`<pheno>.ru.felix.tsv`. The branch reads `GT:AN1:AN2`, so it needs
`joint_vcfs`. RU results are written separately and are not yet merged into
`SummarizeFelixResults`.

### 6.1 Locus definition and allele rules

From `aggregate_repeat_loci.py`, shared with Figure 2C so the figure's loci and
the tested loci are the same objects:

- **Catalog.** TRExplorer v1.0.1 intervals (`--catalog-bed`: chrom, start, end,
  TRID, motifs). This settles `period` and the motif definitions.
- **Assignment.** A record goes to the catalog locus it overlaps most, within
  `--flank` bp (default 10). Records with a length change above `--max-bp`
  (default 100 kb) are ignored.
- **Sign.** Per ALT, `len(ALT) − len(REF)` for sequence-resolved alleles;
  symbolic ALTs fall back to `SVLEN` / `END` / `SVTYPE`, with `DEL` negative and
  `INS`/`DUP`/`CNV` positive. INV and BND contribute 0. Each ALT is signed
  independently, so a record carrying both an insertion and a deletion allele
  codes them in opposite directions.
- **Impure alleles are kept.** The dosage is `bp / period` as a float. There is
  no whole-number-of-units or perfect-tiling requirement, so interrupted
  repeats are not discarded.
- **Size split.** With `--small-vcf` and `--sv-vcf`, an allele counts only from
  the callset that owns its size (`|change| < --split-bp`, default 50, is
  small), so an event called by both DeepVariant and the SV caller counts once.
  Haplotypes carrying both a `<50 bp` and a `≥50 bp` change in the same
  direction are counted in `hap_possible_duplicate`.
- **Missingness.** A haplotype is missing at a locus when any of its
  contributing records is missing. A sample is missing on both haplotypes when
  two or more of its records at the locus are unphased heterozygous length
  changes, because the changes cannot be assigned to haplotypes
  (`hap_unphased_ambiguous`). One unphased het plus any number of homozygous
  records is resolved.
- **Ancestry.** `AN1`/`AN2` are taken from the first record at the locus that
  carries them, which assumes ancestry is constant across a tandem array. FLARE
  tracts are orders of magnitude longer than a tandem array, so this should
  hold; `hap_ancestry_inconsistent` in the summary measures it rather than
  assuming it.

### 6.2 Superseded: the record-level path

`scripts/annotate_repeat_units.py` (INFO `RU_TEST`/`RU`/`PERIOD`/`CN_REF`)
annotates **one VCF record at a time** and is no longer the association input.
It is retained for `ru_engine = "rust"` and the encoding comparison. Four
defects made it unsuitable as the association input, all of them resolved by
locus-level aggregation (§6.1):

1. **A locus is split across records.** Each record got its own dosage and its
   own test, so the multi-allelic repeat locus was never a single tested unit.
2. **Impure alleles dropped the whole record.** `_units_for_period` returns
   `None` unless the length change is a whole number of periods or the inserted
   sequence tiles perfectly, and `try_period` then rejects the record — losing
   every allele at it, not just the impure one.
3. **Sequence-resolved deletions got the wrong sign.** `svtype_sign` in the
   dosage writer treats an allele as a deletion only on `SVTYPE=DEL` or an ALT
   starting `<DEL`, so a DeepVariant deletion (`REF=ACAG, ALT=A`) was coded as
   a gain.
4. **Mixed records got a single sign.** The sign came from the first `SVTYPE`
   token, so a record with both an insertion and a deletion ALT coded both in
   the same direction.

Fixture evidence: `eval/record_vs_locus/`.

---

## 7. Open items

| Item | State |
|---|---|
| chr22 RU branch on real `joint_vcfs` | not run — the blocking gate |
| A-1 / A-2 on real data | not run (synthetic only, §4.1) |
| FLARE2 `nanc` → `num_ancs` / `K` | undecided; sets the width of every test |
| FLARE marker distance at repeat SVs | not measured. Repeat SVs sit in regions the context mask removes, so the production LAI recipe should not use the context mask |
| SHAPEIT4 switch error at repeat sites | not measured. This is the phased-data analogue of the aggregator's unphased-ambiguity rule: trio-based switch error at SV sites, repeat SVs separately |
| Per-locus Mendelian length inconsistency | not measured. Long VNTR alleles have length-calling noise; record it and drop or flag loci above a threshold |
| Ancestry-centred sensitivity encoding | needs a redesign. Centring is a shift (§2.1) **and** makes every haplotype non-zero, so it is incompatible with carrier QC as written |
| Segmental-duplication CNVs through the same path | undecided; they are copy-number dosages too, but phasing and ancestry assignment are less reliable |
| RU results in `SummarizeFelixResults` | not wired; written separately |
| Upstreaming the carrier-QC switch | deferred until the pipeline works end to end on AoU data |

---

## 8. Conventions that must not drift

- Headline p-column: `P_cct_admixed_c` (FELIX) vs `P` (Tractor-Mix).
- Covariate arms for the matched 2×2: **limited** = sex, age, PC1–PC10;
  **full** = limited + standardized coverage + GC dummies.
- Do not headline an encoding's "more hits" without matched λGC
  ([`eval/README.md`](eval/README.md)).
- Compare repeat λGC against the SNV λGC **from the same null**, not against a
  separate fit.
