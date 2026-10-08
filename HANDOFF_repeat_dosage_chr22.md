# Handoff: chr22 repeat-dosage run (M6b)

First real-data run of the locus-level repeat-dosage association path. Design
and conventions: [`felix/REPEAT_DOSAGE.md`](felix/REPEAT_DOSAGE.md). Gates:
[`felix/eval/README.md`](felix/eval/README.md).

Written 2026-10-07. Everything below the "Before you submit" section is
unexecuted.

---

## What changed since the last handoff

The repeat branch of `FelixGenome.wdl` now tests **one catalog repeat locus at a
time** instead of one VCF record at a time.

| Before | After |
|---|---|
| `AnnotateRuTest` (`annotate_repeat_units.py --simple-repeat-bed`) | `AggregateRepeatLoci` (`aggregate_repeat_loci.py --out-vcf --with-ancestry`) |
| One test per VCF record; a repeat array split across a GLnexus indel and a merged SV record became two tests | One test per TRExplorer locus; a haplotype's dosage is the sum of its signed length changes across the locus's records |
| A record was dropped unless every ALT was a whole number of repeat units | Dosage is `bp / period` as a float; impure alleles are kept |
| Mixed insertion/deletion records took one sign from the first `SVTYPE` token, or were dropped as `UNK` | Each ALT is signed independently by `len(ALT) − len(REF)` |
| Any missing or unphased genotype aborted the dosage writer | `--missing ref` (default): the haplotype contributes 0 and leaves `ANC{k}`; `--max-hap-missing` drops a locus above a threshold |

Fixture evidence: the record-level path delivered 12 of 18 per-ancestry dosages
correctly on `scripts/testdata/repeat_locus_assoc/`, the locus path 18 of 18
(`felix/eval/record_vs_locus/defect_summary.tsv`).

Two smaller fixes in the same area:

- `RU_DOSAGE` was written at 4 significant figures by `fmt_units`, which is fine
  for the Figure 2C allele table but put a ~1e-4 relative error into every
  non-integer dosage. `--out-vcf` now uses a full-precision formatter;
  `SVLEN` and `PERIOD` are both present so the ratio stays recomputable.
- `AN1`/`AN2` were read one byte at a time, so a two-digit ancestry code folded
  into its first digit (`12` → `1`). Now parsed as decimal integers on both the
  fast and slow paths, with `--num-ancs` to fail loudly on an out-of-range code.
  **This is defensive, not currently load-bearing:** the FLARE2 candidates are
  `nanc` 5 and 6, so today's codes are single-digit.

---

## Before you submit

1. **Settle `num_ancs`.** The config ships `5`. Selection v2 is still choosing
   between `sel_chr20_flare2_nanc5` and `sel_chr20_flare2_nanc6`
   (`flare/README.md` Part 8 → `selection_v2/selection_decision.json`). `K` sets
   the width of every per-ancestry test, so running with the wrong value
   invalidates the comparison against the SNV scan. If the decision is `nanc 6`,
   set `num_ancs` to 6 in the inputs JSON.
2. **Confirm the catalog is staged.** `repeat_catalog_bed` points at
   `gs://BUCKET/refs/repeat_loci/trexplorer_v1.0.1.catalog.bed.gz`, the same
   object `repeat_loci/configs/*.json.example` uses. Build it with
   `trgt_plvi.py catalog --vcf <one TRGT VCF>` if absent.
3. **Stage the scripts.** `aggregate_repeat_loci.py` is new to this workflow and
   has been added to `scripts/stage_tractor_scripts.sh`:
   ```bash
   WORKSPACE_BUCKET=gs://fc-secure-... ./scripts/stage_tractor_scripts.sh
   WORKSPACE_BUCKET=gs://fc-secure-... ./scripts/stage_felix_scripts.sh
   ```
   Or run `notebooks/terra/00_sync_repo.ipynb` on the Terra VM. Re-register
   `FelixGenome.wdl` in the Terra Methods Repository — its input set changed.
4. **Run the local gates first.** They need only `bcftools`:
   ```bash
   python3 -m pytest scripts/test_aggregate_repeat_loci.py \
     scripts/test_write_admixed_dosage_vcf.py -q
   python3 felix/eval/dosage_qc/locus_path_smoke.py --out-dir /tmp/smoke
   ```
5. **Image.** `felix-pilot:0.2.0`
   (`sha256:a9b51dad68b6a8881357cbb9ce2317a68a5d00eadf6c6bc65b43612935c5868b`)
   already carries the carrier-QC patch. No rebuild is needed for this change —
   only staged scripts and the WDL changed.

---

## Submit

```bash
# felix/configs/felix.chr22_repeat_dosage.inputs.limited.json.example
#   chroms        = ["chr22"]
#   joint_vcfs    = aou_lr_phase2_v1.chr22.anc.vcf.gz   (GT:AN1:AN2)
#   grm_vcfs      = chr1 + chr22, same as the Tractor-Mix / SAIGE pilots
#   ru_engine     = "felix"
```

Replace `BUCKET` with the workspace bucket, set `num_ancs`, and submit
`FelixGenome.wdl`. The repeat branch requires `joint_vcfs` — `phase_vcfs` alone
has no `AN1`/`AN2`.

Workflow path for this run:

```
Check → MakeGRM → FitFelixNull (per phenotype)
      → AggregateRepeatLoci → WriteRuAdmixedVcf → RunFelixStep2RuVcf (per phenotype)
      → RuConcat
      (and the ordinary FELIXla Pack → Step2 → Summarize for SNVs/indels/unique SVs)
```

Outputs to pull back:

| Output | Use |
|---|---|
| `ru_results_tsvs` | per-phenotype repeat results (`P_cct_admixed_c`) |
| `ru_locus_summaries` | locus/haplotype tallies — read these before the p-values |
| `ru_admixed_stats` | writer stats: `hap_missing`, `missingness_skipped`, `written` |
| `ru_locus_alleles` | per-locus allele spectrum (no sample ids) |
| `results_tsvs` | the SNV/indel scan from the **same null**, for the λGC comparison |

---

## Read the summaries before the p-values

`ru_locus_summaries` and `ru_admixed_stats` are the first thing to look at;
several of them test assumptions that have only ever been checked on synthetic
data.

| Field | Where | What a bad value means |
|---|---|---|
| `hap_ancestry_inconsistent_frac` | locus summary | Records at one locus disagree on local ancestry. Should be ~0. A non-trivial value means "ancestry is constant across a tandem array" does not hold, and the first record's call (which is the one used) is not safe |
| `anc_code_max` | locus summary | Must be `num_ancs − 1`. Anything else means the FLARE model and `num_ancs` disagree |
| `hap_possible_duplicate_frac` | locus summary | Haplotypes carrying both a `<50 bp` and a `≥50 bp` change in the same direction at one locus — the DeepVariant/SV double-count. Figure 2C's `TrgtLocusConcordance` measures how often the sum overshoots TRGT |
| `hap_unphased_ambiguous` | locus summary | Samples set missing because two or more unphased heterozygous changes could not be assigned to haplotypes. Should be small on a phased callset; a large value points at phasing, not at the aggregator |
| `missingness_skipped` | writer stats | Loci dropped by `--max-hap-missing 0.1`. If this is large, raise the cap deliberately and record it rather than quietly |
| `hap_missing / hap_total` | writer stats | Overall dropout absorbed by `--missing ref` |

---

## Gate criteria

Pass requires all of:

1. **λGC ≈ 1** on null phenotypes for the repeat results, compared against the
   SNV λGC **from the same null fit** — not a separate fit, and not a published
   value. `felix/scripts/compare_repeat_encodings.py` for the per-encoding
   table.
2. **A-1 on real data.** Code a set of biallelic chr22 indels through both the
   locus path and FELIXla and confirm identical p-values. The synthetic version
   passed at the writer boundary (`REPEAT_DOSAGE.md` §4.2); this is the version
   that counts, because it exercises FELIX's two input readers rather than ours.
3. **A-2 on real data.** The same loci at two per-locus scales give identical
   p-values, with `BETA` scaling. `DS` invariance already holds exactly by
   construction (§4.2), so a failure here is a FELIX-side finding.
4. **Carrier QC did what it claims.** With `FELIX_DOSAGE_QC=carrier`, no locus
   should lose an ancestry test because its summed dosage was negative. Compare
   the number of per-ancestry tests reported against `num_ancs` × loci tested;
   stock FELIX dropped ~14% of them on synthetic contraction-heavy data.
5. Positive-control loci checked if present in the callset: HTT, FMR1, FXN,
   C9orf72, RFC1, plus the published STR–trait pairs (APOB CTG, CBL CGG, TAOK1
   polyA; Margoliash et al. 2023) against nearby blood traits.

Do not headline an encoding's "more hits" without matched λGC.

---

## Still unmeasured, and not blocked by this run

These are the upstream checks from `REPEAT_DOSAGE.md` §7. None of them is
addressed by the code changes above, and the first two bound how much the
results can be trusted.

1. **FLARE marker distance at repeat SVs.** Repeat SVs sit in regions the
   context mask removes, so the production LAI recipe should not use the context
   mask. Measure with `propagate_flare_ancestry --stats-json`.
2. **SHAPEIT4 switch error at repeat sites.** SHAPEIT4 statistically re-phased
   the HiPhase blocks. Measure trio-based switch error at SV sites, repeat SVs
   separately. This is the phased-data analogue of the aggregator's
   unphased-ambiguity rule: the aggregator can detect an unassignable change,
   but it cannot detect a confidently wrong haplotype assignment.
3. **Per-locus Mendelian length inconsistency.** Long VNTR alleles have
   length-calling noise. Record it per locus and drop or flag above a threshold.
4. **Ancestry-centred sensitivity encoding.** Needs a redesign: centring is a
   shift (`REPEAT_DOSAGE.md` §2.1) *and* makes every haplotype non-zero, so it
   is incompatible with carrier QC as written. Without it there is no clean
   separation of repeat-length effect from admixture mapping at top loci.
5. **RU results are not merged into `SummarizeFelixResults`.** They come back as
   `ru_results_tsvs` and need their own QQ/Manhattan pass.
6. **Compound loci.** At a TRExplorer locus with more than one motif, `period`
   is the **first** motif's length for every record, so `RU_DOSAGE` there is
   "units of the primary motif". The allele table carries `n_motifs`, so these
   loci can be filtered or re-examined; decide which before reporting effect
   sizes per repeat unit.
