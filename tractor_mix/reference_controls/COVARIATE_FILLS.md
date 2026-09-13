# Long-read covariate fills (provenance)

This document records **what** was added or backfilled into
`tractor_mix/covariates.source_rebuilt.csv.gz` for long-read analyses, **how**
each fill was derived, and **where** the audit trails live.

Target use cases: joint-callset PCA, within-population PCA, and association
models that need ancestry / continental population / LR PCs for AoU samples and
HPRC/HGSVC3 `HG*`/`NA*` controls.

Final table size after fills: **17,519** rows (17,226 AoU + **293** controls),
**198** columns.

---

## Column semantics (do not confuse)

| Column | Meaning |
| --- | --- |
| `ancestry_pred` / `ancestry_pred_other` | Predicted continental ancestry (lowercase). Hard `ancestry_pred` is `afr`/`amr`/`eas`/`eur`/`sas`/`mid` (never `oth` after Rule C). `ancestry_pred_other` allows `oth`. |
| `population` | Continental membership for **short-read** within-pop PCA training (uppercase `AFR`/`AMR`/…). |
| `lr_pop_population` | Subset used for **long-read** within-pop PCA (lowercase; equals `ancestry_pred_other` for the joint callset). |
| `lr_PC*` / `has_lr_pcs` | Global long-read PCs from the DeepVariant joint callset. |
| `lr_pop_PC*` / `has_lr_pop_pcs` | Within-population long-read PCs. |
| `PC*` / `pop_PC*` | Short-read PCs (unchanged by these fills). |
| `is_reference_control` | `True` for HPRC/HGSVC3/GIAB `HG*`/`NA*` rows. |

---

## Pipeline order

Rebuild source covariates, then apply LR fills in this order:

1. `notebooks/terra/tractor_00_cov_rebuild_source.ipynb`
2. `scripts/merge_lr_global_pcs_into_covariates.py` — global `lr_PC*` + control rows
3. `notebooks/terra/tractor_06_fill_lr_ancestry.ipynb` (or `scripts/fill_lr_ancestry_from_pcs.py`) — ancestry / population gaps
4. `notebooks/terra/tractor_07_pca_within_population.ipynb` — produce `population_pcs.tsv`
5. `scripts/merge_lr_pop_pcs_into_covariates.py` — `lr_pop_PC*`
6. `scripts/apply_lr_soft_field_fills.py` — soft `population` / `sex_at_birth` nits

Canonical PC inputs:

- `tractor_mix/pca/deepvariant_lr_v1/global_pcs.tsv`
- `tractor_mix/pca/deepvariant_lr_v1/population_pcs.tsv`

---

## 1. Reference controls (293 `HG*` / `NA*`)

**Script:** `scripts/merge_lr_global_pcs_into_covariates.py`  
**Metadata:** `control_sample_metadata.tsv`

| Field | How filled |
| --- | --- |
| Row presence | Every `HG*`/`NA*` ID in `control_sample_metadata.tsv` is appended if absent from AoU covariates. **292** of these are in `global_pcs.tsv`. `HG02015` is in the joint VCF but failed PCA sample QC (call rate `< 0.98`), so it is a metadata-only row (`has_lr_pcs=False`, no `lr_PC*`) |
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
  (12,261 AoU + 292 controls). `HG02015` is a 293rd control row without `lr_PC*`
  because it failed PCA sample QC (call rate `< 0.98`), not because it is absent
  from the joint VCF. Other samples without `lr_PC*` are mostly releasable
  long-read IDs that were never in this joint callset.
- Does **not** overwrite short-read `PC1`–`PC16` / `has_global_pcs`.
- Releasable samples outside the joint callset (~1,269) correctly remain without
  `lr_PC*` (no genotypes → no LR PCs).

---

## 3. Ancestry backfills (AoU gaps)

**Notebook:** `notebooks/terra/tractor_06_fill_lr_ancestry.ipynb`  
**CLI:** `scripts/fill_lr_ancestry_from_pcs.py`  
**Audit:** `lr_ancestry_knn_fills.tsv` (**91** applied fills: 80 historical + 11 hard-pred remediations)

Only runs on long-read-flagged rows. Never overwrites existing non-missing
hard `ancestry_pred`, non-missing `ancestry_pred_other`, or `population`, except
Rule C which remediates soft `ancestry_pred=oth` → a hard continental label.

**Hard vs soft:** `ancestry_pred` is the hard six-class call
(`afr`/`amr`/`eas`/`eur`/`sas`/`mid`) and must **not** be `oth`.
`ancestry_pred_other` (and continental `population`) may still be `oth`/`OTH`.
FLARE population splits that want no OTH shard should use `ancestry_pred`
(uppercased), not `population`.

### Rule A — copy from continental `population` (58 fills)

If ancestry is missing but `population` is present:

- `ancestry_pred_other` ← `population.lower()` (including `oth`)
- `ancestry_pred` ← same **only** when the label is hard (not `oth`)
- `has_ancestry_annotation` ← `True`

`population=OTH` therefore fills soft other only; hard `ancestry_pred` comes
from Rule B/C.

Audit `method`: `from_population`.

### Rule B — kNN + centroid on long-read PCs (22 fills)

If **both** ancestry and `population` are missing **and** `has_lr_pcs`:

1. **Soft reference set:** `has_lr_pcs` samples with non-missing
   `ancestry_pred_other` (may include `oth`) → vote for
   `ancestry_pred_other` / `population`.
2. **Hard reference set:** `has_lr_pcs` samples with a hard continental label
   on `ancestry_pred` (else hard `ancestry_pred_other`) → vote for
   `ancestry_pred` (never `oth`).
3. **Features:** `lr_PC1`–`lr_PC10`, z-scored using each reference mean/SD.
4. **Centroid + kNN (k=15)** as before; high confidence if centroid and kNN
   **agree** and (`knn15_frac` ≥ **0.8** **or** `centroid_margin` ≥ **0.5**).
   Gate on the **hard** vote.
5. Writes `ancestry_pred` ← hard kNN majority; `ancestry_pred_other` ← soft
   kNN majority; `population` ← soft label uppercased when population was
   missing.
6. Default run fills low-confidence calls too (`fill_low_confidence=True`).

Audit `method`: `lr_pc_knn`.

### Rule C — hard `ancestry_pred` remediation (11 fills)

If AoU (non-control) `has_lr_pcs` and `ancestry_pred` is missing or `oth`:

- Re-infer a **hard** label from hard-labeled neighbors (same kNN/centroid
  rules as Rule B hard vote).
- Overwrites `ancestry_pred` only; leaves `ancestry_pred_other` and
  `population` unchanged (often still `oth`/`OTH`).
- Skips curated reference controls (e.g. HG002 may remain `oth`).

Applied 2026-09-08 to **11** AoU rows that earlier Rule A/B incorrectly wrote
`oth` into both ancestry columns (9× `from_population` from `population=OTH`,
2× `lr_pc_knn` with soft majority `oth`). Labels after remediation:

| research_id | hard `ancestry_pred` | high_confidence |
| --- | --- | --- |
| 1179960 | eur | False |
| 1286158 | afr | True |
| 1308393 | afr | True |
| 1340678 | amr | False |
| 1360567 | eur | True |
| 1490799 | eur | True |
| 1865868 | eur | True |
| 1958378 | eur | True |
| 2691108 | amr | True |
| 2911678 | afr | True |
| 7959375 | sas | False |

Audit `method`: `lr_pc_knn_hard_pred`. Re-runs merge into the existing audit
TSV (replace prior Rule C rows for the same IDs) unless `--replace-audit`.

### Result

- `has_lr_pcs` missing ancestry: **0**
- `has_lr_pcs` AoU with `ancestry_pred` missing/`oth`: **0**
- `final_releasable_v9` missing ancestry: **0**
- ~130 leftover “LR” ancestry gaps remain only among samples flagged via
  `lr_meet_qc` / `lr_phase` but **without** genotypes or `lr_PC*` (mostly
  withdrawn). Not filled; not used in joint-callset PCA/association.
- MID remains a valid hard `ancestry_pred` label. FLARE keeps MID shards with
  `allow_unrepresented_pops=true` for now (no MID panel ancestry); a dedicated
  panel or remap may replace that later.

---

## 4. Long-read within-population PCs

**Script:** `scripts/merge_lr_pop_pcs_into_covariates.py`  
**Input:** `pca/deepvariant_lr_v1/population_pcs.tsv` (from `tractor_07`)

- Adds `lr_pop_PC1`–`lr_pop_PC32`, `has_lr_pop_pcs`, `lr_pop_population`
  for the same **12,553** joint-callset samples.
- Does **not** replace short-read `pop_PC*` / `has_pop_pcs`.
- `lr_pop_population` matches `ancestry_pred_other` for every joint-callset row.

---

## 5. Soft field backfills (joint-callset nits)

**Script:** `scripts/apply_lr_soft_field_fills.py`  
**Audit:** `lr_soft_field_fills.tsv` (**4** rows)

```bash
python3 scripts/apply_lr_soft_field_fills.py \
  --audit tractor_mix/reference_controls/lr_soft_field_fills.tsv
# optional after rebuilds:
#   --discover
```

Applied after PC/ancestry merges so continental population and usable sex are
complete for every `has_lr_pcs` sample. Default mode replays the audit; `--discover`
also fills any remaining joint-callset gaps with the same rules.

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
| `lr_ancestry_knn_fills.tsv` | Ancestry Rule A/B/C audit (91 rows) |
| `lr_soft_field_fills.tsv` | Soft `population` / `sex_at_birth` audit (4 rows); replay via `scripts/apply_lr_soft_field_fills.py` |
| `../pca/deepvariant_lr_v1/global_pcs.tsv` | Global LR PC source |
| `../pca/deepvariant_lr_v1/population_pcs.tsv` | Within-pop LR PC source |
| `../covariates.source_rebuilt.csv.gz` | Filled covariates table |
| `../covariates.source_rebuilt.data_dictionary.tsv` | Column dictionary |
