# Reference-control metadata (HPRC / HGSVC3 / GIAB)

`control_sample_metadata.tsv` annotates the **292** `HG*` / `NA*` samples present in
`tractor_mix/pca/deepvariant_lr_v1/global_pcs.tsv` (Phase-2 long-read joint callset controls,
including the original Phase-1 set of 47).

| Field | Notes |
| --- | --- |
| `ancestry_pred` / `ancestry_pred_other` | AoU-style labels (`afr`/`amr`/`eas`/`eur`/`sas`/`oth`) |
| `sex_at_birth` / `inferred_sex` | From 1KG / HPRC / Coriell public metadata |
| `population_code` | Optional fine-grained 1KG/HPRC subpopulation (e.g. `GBR`, `YRI`); **not** written to covariates `population` |
| `is_phase1_control` | True for the 47 Phase-1 controls (also have SV sens07 counts) |
| `ancestry_source` | Provenance tag for the ancestry fill |

In `covariates.source_rebuilt.csv.gz`, control `population` is the **continental**
code (`AFR`/`AMR`/`EAS`/`EUR`/`SAS`/`OTH`), derived from `ancestry_pred_other`, so
controls join the same within-population PCA subsets as AoU samples.

Used by `scripts/merge_lr_global_pcs_into_covariates.py` to append control rows to
`covariates.source_rebuilt.csv.gz`.

## Full fill provenance

See **[COVARIATE_FILLS.md](COVARIATE_FILLS.md)** for the complete record of long-read
covariate fills: control rows, global / within-pop `lr_PC*`, ancestry Rule A
(population copy) and Rule B (centroid + kNN on `lr_PC1`–`lr_PC10`), and soft
`population` / `sex_at_birth` backfills.

| Audit | Contents |
| --- | --- |
| `lr_ancestry_knn_fills.tsv` | 80 ancestry fills (`from_population` / `lr_pc_knn`) |
| `lr_soft_field_fills.tsv` | 4 soft field fills (`population` / `sex_at_birth`) |

Replay soft fills: `scripts/apply_lr_soft_field_fills.py` (add `--discover` to
catch new joint-callset gaps).

Mainline notebook for ancestry: `notebooks/tractor_05a_fill_lr_ancestry.ipynb`
(CLI: `scripts/fill_lr_ancestry_from_pcs.py`).
