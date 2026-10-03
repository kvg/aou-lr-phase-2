# Handoff: LAI recipe selection v2 + FELIX (updated 2026-10-02)

**Audience:** the next Claude / Cursor session on this work.
**Repo:** `kvg/aou-lr-phase-2`, branch `main`.
**Companion docs:** [`README.md`](README.md) ("FLARE2 mode", "Recipe selection v2"),
[`../felix/PLAN.md`](../felix/PLAN.md), [`../felix/SV_SCORER_DESIGN.md`](../felix/SV_SCORER_DESIGN.md).

---

## 1. Goal

Local-ancestry-informed association testing of long-read SNVs, indels and SVs
(including repeat-mediated SVs) across as many AoU phenotypes as possible,
for the **whole cohort**.

| Piece | Choice |
|---|---|
| LAI | **FLARE2** (unless selection v2 says the original-FLARE pin is better) |
| Association | **FELIX** (replaces Tractor-Mix) |
| Phasing | HiPhase (physical) → SHAPEIT4 (statistical); SVs and SNVs on the same haplotypes |
| Repeat SVs | Length dosage × ancestry with SPA through FELIX's tests (design note) |
| FELIX fork | Allowed if needed; last resort |

---

## 2. What changed on 2026-10-02 and why

A review of the old Part 8 plan found:

1. **`mean_ll` compared recipes on different sample sets.** `lai_exp.tsv` has
   no `analysis_samples` column, so each recipe was scored on whichever samples
   it painted (AFR-only, AMR-only, AFR+AMR, all populations). The AFR-only
   "win" (−0.4885) mostly reflects the sample set. The old table also mixed a
   10 Mb chr22 window with whole chr20.
2. **The Mendelian gate never ran.** In `flare_02`, `PED` was relative to the
   notebook directory and `MAP_TMPL` pointed at a `maps/` directory that does
   not exist. The code skipped silently. The scorer also rejected bare `chr20`
   regions.
3. **The Mendelian check has little power on the chr22 window.** The pedigree
   has 101 children with both parents (71 families). 10 Mb gives about 0.1
   crossovers per meiosis.
4. **`mean_ll` uses allele frequencies from FLARE's own reference panel,** so it
   can favour recipes that over-fit locally. Use it as a gate only.
5. **Null-λ only measures calibration,** was built around Tractor, dropped PCs
   (not the production model), and was never shown to respond to LAI quality.
6. **The old per-population FLARE2 rows would give inconsistent ancestry labels
   across shards.** Upstream clustering names ancestries `anc_0..` in arbitrary
   order with an unseeded GMM.

So: the old 18-recipe `mean_ll` table and the suggested `−0.53` gate are
**retired**. Do not use them to pick a recipe.

---

## 3. What is built (code done, not run on Terra)

| Item | Where |
|---|---|
| FLARE2 mode in the WDL (pooled training subsample → `panel-probs` → clustering → shared `template_model`, EM on) | `wdl/FlareByPopulation.wdl` (`flare2_nanc`, `flare2_model`, …); miniwdl check clean |
| FLARE2 image (flare.jar at upstream commit `87573be` + scikit-learn + bcftools) | `docker/Dockerfile.flare2`, `./build_docker.sh --flare2` → `aou-flare2:0.6.0-87573be` |
| Stable FLARE2 labels + autocorrelation gate | `scripts/flare2_build_model.py` (+ `test_flare2_build_model.py`) |
| Selection rows | `configs/lai_exp.tsv`: `selection_set=v2` on `chr20_flat_props_pin_t`, `sel_chr20_flare2_nanc5`, `sel_chr20_flare2_nanc6`; new columns `flare2_nanc`, `flare2_em`, `selection_set` |
| Genome-wide FLARE2 apply config | `configs/flare2.apply.inputs.json.example` |
| Negative controls (within-population track permutation) | `scripts/flare_make_negative_control.py` |
| Allele scorer reads `##ANCESTRY` header; `--model` maps FLARE2 clusters to panel AFs via P | `scripts/flare_score_allele_ancestry.py` |
| Mendelian scorer accepts bare contigs; expected crossovers over the scored span | `scripts/flare_score_mendelian_lai.py` |
| Decision rule + trio-bootstrap CI | `scripts/flare_lai_exp.py` (`select_recipe`, `trio_bootstrap_violation_ci`); tests in `test_flare_selection_v2.py` |
| Notebook Part 8 rewritten | `notebooks/terra/flare_02_lai_exp_compare.ipynb` |
| FELIX pilot configs on one variant set; M1 gate notebook | `felix/configs/felix.pilot.inputs.*.json.example`, `notebooks/terra/felix_01_pilot_gate.ipynb` |
| Repeat-SV scorer design | `felix/SV_SCORER_DESIGN.md` |

---

## 4. Next steps (Terra), in order

The two tracks below can run in parallel.

### Track A: LAI recipe (M0)

1. Build and push the FLARE2 image: `cd flare && ./build_docker.sh --flare2`.
2. Run `00_sync_repo.ipynb` (stages the new scripts, notebooks and WDL; upserts
   `flare_lai_exp`). Set the new column types (`flare2_nanc` number, `flare2_em`
   boolean). Add the four `flare2_*` inputs from
   `configs/lai_exp.inputs.json.example` to the method config, and map the
   workflow output `flare2_model_out` to a table column of the same name.
3. Launch `sel_chr20_flare2_nanc5` and `sel_chr20_flare2_nanc6`. Expect a
   training subset stream (~90 min), the panel-probs pass, then the usual
   shards. If training fails the autocorrelation gate, that `nanc` is too high.
   Record it and drop the row.
4. Run `flare_02` top to bottom (switch-QC diagnostics, then Part 8).
   Part 8 writes `selection_v2/selection_decision.json`:
   - `winner` → production recipe.
   - `tie_human_decision` → choose by hand (tie-break: smaller
     |excess_recomb_over_expected|) and write the reason in `README.md`.
   - `metric_invalid_…` → the Mendelian metric cannot beat the negative control;
     stop and rethink before choosing.
5. If FLARE2 wins: apply its trained model genome-wide with
   `configs/flare2.apply.inputs.json.example`.

### Track B: FELIX (M1, M1b, M6)

1. M1: submit FelixPilot limited + full (same FLARE-marker VCF as the
   Tractor-Mix pilot), then `felix_01_pilot_gate.ipynb`. Pass = λGC ≈ 1 on null
   phenotypes and comparability `OK`.
2. M1b: rerun FELIX with the full long-read callset (`phase_vcf` + `flare_vcf`).
3. M6: build and push `felix-pilot:0.2.0` (FELIX v0.1 + carrier-QC patch), then
   run `FelixGenome`'s RU branch on chr22 `joint_vcfs`. Patch, writer and WDL
   wiring are done and validated locally (`SV_SCORER_DESIGN.md` §3.2, §7).

---

## 5. Terra quirks (still true; don't regress)

1. Prefer the git clone's `scripts/` over `$WORKSPACE/edit/scripts` (stale copies).
2. Terra bcftools is old: no `query -m2/-M2/-v` or `--threads`; `-R` + `-r`
   together breaks GT pulls. Scorers use the panel `-R` only.
3. `Path("gs://…")` collapses to `gs:/`; localize before using pathlib.
4. Chr20 rows have a blank `region`; Part 8 infers `chr20`.
5. Laptop `gsutil` cannot reach AoU buckets; stage from a Terra notebook.

## 6. Locations

| Item | Location |
|---|---|
| Workspace bucket | `gs://fc-secure-8f7d6a20-04ce-40d7-8c88-aececeac3e09` |
| Eval panels | `$BUCKET/refs/flare/eval_panels/{chr22_10mb,chr20_full}.markers.tsv` |
| Selection outputs | `$OUT/selection_v2/` (scores, controls, `selection_decision.json`); decision staged to `$BUCKET/refs/flare/selection_v2/` |
| Notebook working dir | `~/AoU_DRC_LongReads_PhaseTwo_Storage/edit/flare_lai_exp/` (`OUT`) |
| Git clone on the VM | `~/AoU_DRC_LongReads_PhaseTwo_Storage/edit/aou-lr-phase-2/` |

## 7. Known open items

- Repeat SVs: repeat dosage is REF-relative; FELIX runs it through its
  dosage-VCF input with our carrier-QC patch (no dropped tests, calibrated on
  synthetic nulls). Nothing has run on AoU data yet, and `felix-pilot:0.2.0`
  must be built before any FELIX workflow that defaults to it.
- FLARE2 `flare2_em=true` re-estimates T per shard. Earlier per-population EM
  runs drifted to high T. Check `models.tsv` `t_gen` for the FLARE2 rows before
  trusting them.
- MID has no reference panel; FLARE2 clustering is the intended fix. Check
  where MID haplotypes land in the cluster labels (`flare2_labels`).
- Call-QC site-filter rows (`*_filter_call_qc`) remain deferred.
