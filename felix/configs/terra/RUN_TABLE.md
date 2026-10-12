# Table-driven FelixGenome runs

A run configuration is a **row**, not a JSON file. To change a parameter, add a
row with the new value and submit on it; the old row stays as the record of what
produced the old results.

## Files

| File | Upload as | What it is |
| --- | --- | --- |
| `felix_runs.entity.tsv` | Data table | One row per run configuration, 33 analysis columns + `run_description`. Entity type `flare_lai_prod_set`. |
| `felix_runs.membership.tsv` | Data table | The 22 autosome members of each row, in chromosome order. |
| `FelixGenomeTable.inputs.json` | Method config inputs | Reads `${this.<column>}` for analysis parameters and `${this.flare_lai_prods.<column>}` for the per-chromosome arrays. |

Upload the entity TSV first, then the membership TSV. Both are needed: the first
sets the attributes, the second sets the members.

## What lives where

**In the table (33 columns).** Everything that defines the analysis: cohort and
phenotype files, covariate set, GRM inputs and relatedness cutoff, `num_ancs`,
`min_mac`, the whole `ru_*` group, the summariser thresholds, the docker image,
and every staged script path. Column names are the WDL input names verbatim, so
the mapping is mechanical and a new input is a new column plus one JSON line.

**In the JSON (8 keys).** Resource sizing only — `pack_*`, `step2_*`,
`ru_disk_*`. These change cost and runtime, not results, and putting them in the
table would widen it by a quarter for no scientific benefit. Move one into the
table if you ever want a row to be a complete cost record too.

**Not set at all.** `precomputed_null_*` (fit nulls in-workflow),
`simple_repeat_bed` / `annotate_repeat_units_script` (the legacy record-level
lane, inactive while `ru_engine = "felix"`).

## Rows

| Row | Differs from baseline in |
| --- | --- |
| `felix.wgs.limited.v1` | — baseline; reproduces the configuration submitted 2026-10-12 exactly |
| `felix.wgs.full.v1` | `covariate_columns` only |

`felix.wgs.full.v1` reuses the GRM and every chromosome pack by call caching;
only the 10 null fits and the step-2 shards rerun.

## Writing cells

- **Floats as plain decimals.** `0.00000005`, never `5e-8` — Terra's input editor
  rejects scientific notation for a Float.
- **`grm_vcfs` is a list cell**, written as a JSON array: `["gs://...","gs://..."]`.
  Terra parses that into a list attribute on upload.
- **Check the types once after the first upload.** Terra infers a cell's type
  from its content, so an `Int` column with a stray space becomes a String. The
  method config editor validates inputs against the WDL and flags a type
  mismatch in red before you can submit — look at it after the first upload of a
  new column rather than finding out from a failed submission.

## Adding a configuration

1. Copy a row in `felix_runs.entity.tsv`, rename the id, change the cells.
2. Copy that row's 22 lines in `felix_runs.membership.tsv` under the new id.
3. Upload both; submit the same method config on the new row.

Call caching does the rest: a row that changes only `min_mac` reuses the GRM,
the nulls and every pack, and reruns only step 2.

## Note on the set table

`flare_lai_prod_set` also accumulates rows named after workflow submissions
(`PropagateFlareAncestry_2026-10-10T22-08-23` and similar), auto-created by Terra
when a workflow is launched over multiple entities. Those rows have the config
columns empty. Submit FELIX runs only on the `felix.*` rows.
