# Long-read covariate fills (provenance)

This document records **what** was added or backfilled into
`tractor_mix/covariates.source_rebuilt.csv.gz` for long-read analyses, **how**
each fill was derived, and **where** the audit trails live.

Target use cases: joint-callset PCA, within-population PCA, and association
models that need ancestry / continental population / LR PCs for AoU samples and
HPRC/HGSVC3 `HG*`/`NA*` controls.

Final table size after fills: **17,518** rows (17,226 AoU + **292** controls),
**198** columns.

---

## Column semantics (do not confuse)

| Column | Meaning |
| --- | --- |
| `ancestry_pred` / `ancestry_pred_other` | Predicted continental ancestry (lowercase). `ancestry_pred_other` allows `oth`. |
| `population` | Continental membership for **short-read** within-pop PCA training (uppercase `AFR`/`AMR`/…). |
| `lr_pop_population` | Subset used for **long-read** within-pop PCA (lowercase; equals `ancestry_pred_other` for the joint callset). |
| `lr_PC*` / `has_lr_pcs` | Global long-read PCs from the DeepVariant joint callset. |
| `lr_pop_PC*` / `has_lr_pop_pcs` | Within-population long-read PCs. |
| `PC*` / `pop_PC*` | Short-read PCs (unchanged by these fills). |
| `is_reference_control` | `True` for HPRC/HGSVC3/GIAB `HG*`/`NA*` rows. |

---

## Pipeline order

Rebuild source covariates, then apply LR fills in this order:

1. `notebooks/tractor_00_cov_rebuild_source.ipynb`
2. `scripts/merge_lr_global_pcs_into_covariates.py` — global `lr_PC*` + control rows
3. `notebooks/tractor_05a_fill_lr_ancestry.ipynb` (or `scripts/fill_lr_ancestry_from_pcs.py`) — ancestry / population gaps
4. `notebooks/tractor_05b_pca_within_population.ipynb` — produce `population_pcs.tsv`
5. `scripts/merge_lr_pop_pcs_into_covariates.py` — `lr_pop_PC*`
6. Soft field backfills documented below (`population` / `sex_at_birth` nits)

Canonical PC inputs:

- `tractor_mix/pca/deepvariant_lr_v1/global_pcs.tsv`
- `tractor_mix/pca/deepvariant_lr_v1/population_pcs.tsv`

---

## 1. Reference controls (292 `HG*` / `NA*`)

**Script:** `scripts/merge_lr_global_pcs_into_covariates.py`  
**Metadata:** `control_sample_metadata.tsv`

| Field | How filled |
| --- | --- |
| Row presence | Every sample ID in `global_pcs.tsv` that is `HG*`/`NA*` is appended if absent from AoU covariates |
| `lr_PC1`–`lr_PC32`, `has_lr_pcs` | Joined from `global_pcs.tsv` |
| `ancestry_pred` / `ancestry_pred_other` | Public 1KG / HPRC / Coriell labels mapped to AoU-style continents (`afr`/`amr`/`eas`/`eur`/`sas`/`oth`) |
| `population` | Uppercase of `ancestry_pred_other` (continental only; **not** fine-grained 1KG codes such as `YRI`) |
| `sex_at_birth` / `inferred_sex` | Public sex metadata |
| `is_reference_control` | Set `True` |
| Phase-1 SV counts | For the 47 Phase-1 controls only (`is_phase1_control`) |

Controls intentionally lack AoU phenotypes, CDR/EHR, tech sheets, releasable flags,
and short-read `PC*` / `pop_PC*`.

Fine-grained codes (when known) stay in metadata column `population_code` only.

---

## 2. Long-read global PCs (joint callset)

**Script:** `scripts/merge_lr_global_pcs_into_covariates.py`  
**Input:** `pca/deepvariant_lr_v1/global_pcs.tsv`

- Adds `lr_PC1`–`lr_PC32` and `has_lr_pcs=True` for **12,553** samples
  (12,261 AoU + 292 controls).
- Does **not** overwrite short-read `PC1`–`PC16` / `has_global_pcs`.
- Releasable samples outside the joint callset (~1,269) correctly remain without
  `lr_PC*` (no genotypes → no LR PCs).

---

## 3. Ancestry backfills (AoU gaps)

**Notebook:** `notebooks/tractor_05a_fill_lr_ancestry.ipynb`  
**CLI:** `scripts/fill_lr_ancestry_from_pcs.py`  
**Audit:** `lr_ancestry_knn_fills.tsv` (**80** applied fills)

Only runs on long-read-flagged rows. Never overwrites existing non-missing
`ancestry_pred` / `ancestry_pred_other` / `population`.

### Rule A — copy from continental `population` (58 fills)

If ancestry is missing but `population` is present:

- `ancestry_pred` ← `population.lower()`
- `ancestry_pred_other` ← same
- `has_ancestry_annotation` ← `True`

Audit `method`: `from_population`.

### Rule B — kNN + centroid on long-read PCs (22 fills)

If **both** ancestry and `population` are missing **and** `has_lr_pcs`:

1. **Reference set:** all `has_lr_pcs` samples with non-missing `ancestry_pred_other`.
2. **Features:** `lr_PC1`–`lr_PC10`, z-scored using the reference mean/SD.
3. **Centroid vote:** Euclidean distance to each ancestry-label centroid; nearest
   label = `centroid`. Margin = distance to 2nd-nearest − nearest.
4. **kNN vote:** k = **15** nearest reference neighbors; majority label =
   `suggested`. Fraction = majority count / 15.
5. **High confidence** if centroid and kNN labels **agree** and
   (`knn15_frac` ≥ **0.8** **or** `centroid_margin` ≥ **0.5**).
6. **Applied label:** kNN majority (`suggested`). Also writes
   `population` ← `suggested.upper()` when population was missing.
7. Default run fills low-confidence calls too (`fill_low_confidence=True`;
   CLI omits `--high-confidence-only`). Of the 22 kNN fills, **17** were
   high-confidence and **5** lower-confidence (still applied; see audit
   `high_confidence` column).

Audit `method`: `lr_pc_knn`.

### Result

- `has_lr_pcs` missing ancestry: **0**
- `final_releasable_v9` missing ancestry: **0**
- ~130 leftover “LR” ancestry gaps remain only among samples flagged via
  `lr_meet_qc` / `lr_phase` but **without** genotypes or `lr_PC*` (mostly
  withdrawn). Not filled; not used in joint-callset PCA/association.

---

## 4. Long-read within-population PCs

**Script:** `scripts/merge_lr_pop_pcs_into_covariates.py`  
**Input:** `pca/deepvariant_lr_v1/population_pcs.tsv` (from `tractor_05b`)

- Adds `lr_pop_PC1`–`lr_pop_PC32`, `has_lr_pop_pcs`, `lr_pop_population`
  for the same **12,553** joint-callset samples.
- Does **not** replace short-read `pop_PC*` / `has_pop_pcs`.
- `lr_pop_population` matches `ancestry_pred_other` for every joint-callset row.

---

## 5. Soft field backfills (joint-callset nits)

**Audit:** `lr_soft_field_fills.tsv` (**4** rows)

Applied after PC/ancestry merges so continental population and usable sex are
complete for every `has_lr_pcs` sample.

| research_id | Field | Before | After | Method |
| --- | --- | --- | --- | --- |
| `2638375` | `population` | missing | `EUR` | Uppercase of `lr_pop_population` (`eur`) |
| `1027261` | `sex_at_birth` | missing | `Male` | Map `inferred_sex` `XY` → `Male` |
| `1690804` | `sex_at_birth` | missing | `Female` | Map `inferred_sex` `XX` → `Female` |
| `1902721` | `sex_at_birth` | missing | `Female` | Map `inferred_sex` `XX` → `Female` |

Notes:

- Did **not** overwrite existing `sex_at_birth` sentinels (e.g. `PMI: Skip` on
  `2638375`).
- Sex mapping is only for truly missing `sex_at_birth` with diploid `XX`/`XY`
  `inferred_sex`.

After this step: joint callset missing `population` = **0**, missing
`sex_at_birth` = **0**.

---

## Intentionally not filled

- Fine-grained 1KG/HPRC subpopulation codes into covariates `population`
- Control tech/coverage/platform, phenotypes, short-read PCs
- HPRC pedigree IDs into AoU pedigree columns
- Ancestry for the ~130 non-genotyped / withdrawn leftovers
- Short-read `PC*` for the ~23 AoU joint-callset samples (and all controls)
  that lack short-read PCA membership

---

## File index

| Path | Role |
| --- | --- |
| `control_sample_metadata.tsv` | Control ancestry/sex/SV provenance |
| `lr_ancestry_knn_fills.tsv` | Ancestry Rule A/B audit (80 rows) |
| `lr_soft_field_fills.tsv` | Soft `population` / `sex_at_birth` audit (4 rows) |
| `../pca/deepvariant_lr_v1/global_pcs.tsv` | Global LR PC source |
| `../pca/deepvariant_lr_v1/population_pcs.tsv` | Within-pop LR PC source |
| `../covariates.source_rebuilt.csv.gz` | Filled covariates table |
| `../covariates.source_rebuilt.data_dictionary.tsv` | Column dictionary |
