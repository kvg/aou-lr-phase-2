# Table-driven FelixGenome runs

A run configuration is a **row**, not a JSON file. To change a parameter, add a
row with the new value and submit on it; the old row stays as the record of what
produced the old results.

## Files

| File | Upload as | What it is |
| --- | --- | --- |
| `felix_runs.entity.tsv` | Data table | One row per run configuration, 24 analysis columns + `run_description`. Entity type `flare_lai_prod_set`. |
| `felix_runs.membership.tsv` | Data table | The 22 autosome members of each row, in chromosome order. |
| `FelixGenomeTable.inputs.json` | Method config inputs | Reads `${this.<column>}` for analysis parameters and `${this.flare_lai_prods.<column>}` for the per-chromosome arrays. |

## Uploading

Set tables take **two files in two different formats**, and they are not
interchangeable:

- The membership file's header is `membership:flare_lai_prod_set_id` and its
  second column must be the **member entity type name**, `flare_lai_prod`. It
  takes exactly two columns — it cannot carry attributes.
- The attribute file's header is `entity:flare_lai_prod_set_id` followed by the
  attribute columns. It cannot carry members.

**Upload membership first, then attributes.** Either order creates the rows, but
this one only relies on behaviour already verified in this workspace: an
attribute TSV carrying the id column plus some columns upserts onto existing
rows without touching the columns it omits (the same pattern as
`flare/configs/lai_prod.ligated_cols.tsv`). The reverse order would depend on a
membership upload preserving attributes, which has not been checked here.

`flare_lai_prod` must already exist, which it does — Terra refuses a set table
whose members it does not yet know.

Re-uploading an id overwrites that row and leaves other rows alone, so adding a
configuration never disturbs the existing ones. **To change which chromosomes a
row covers, use a new id** rather than re-uploading the same one: whether a
second membership upload replaces the member list or unions with it is not
established, and a silently unioned set would scatter over duplicated
chromosomes.

## What lives where

**In the table (24 columns).** Everything that defines the analysis: cohort and
phenotype files, covariate set, GRM inputs and relatedness cutoff, `num_ancs`,
`min_mac`, the whole `ru_*` group, the summariser thresholds and the docker
image. Column names are the WDL input names verbatim, so the mapping is
mechanical and a new input is a new column plus one JSON line.

**In the JSON (17 keys).** Two groups that are not analysis parameters. The nine
staged script paths belong to the workflow definition, not to a run — they
change on a restage, and when they do every row should pick the change up at
once rather than needing 24 cells edited. And eight resource-sizing keys
(`pack_*`, `step2_*`, `ru_disk_*`), which change cost and runtime, not results.

The two files have to move together: deleting a column while the config still
reads `${this.<column>}` leaves the method config referencing an attribute that
no longer exists.

**Not set at all.** `precomputed_null_*` (fit nulls in-workflow),
`simple_repeat_bed` / `annotate_repeat_units_script` (the legacy record-level
lane, inactive while `ru_engine = "felix"`).

## The row

One row ships: **`felix.wgs.full.v1`**, the full covariate set, matching the
genome-wide run submitted 2026-10-12. Clone it for a new configuration.

The limited covariate set was the pilot's comparator. It is not run going
forward and has no row here. Should an arm comparison ever be wanted again, the
two differ only in `covariate_columns`, and `MakeGRM`, the 22 packs and the
repeat aggregation never see a covariate file — so they call-cache, and only the
null fits, which are keyed on it, and the step-2 shards rerun.

Note that `FelixGenomeLimited.inputs.json`, one directory up in this folder, is
the pre-table config and is named for the limited arm. It is superseded by
`FelixGenomeTable.inputs.json`; do not take its name as a description of what
ran.

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

## Scaling to the full phecode list

One `pheno_cov.tsv` carrying **every** phecode column, and one phenotype file
per batch. `fit_felix_null.R` reads only the phenotype column it is given and
`summarize_felix_results.py` reads `pheno_cov` with `usecols`, so the extra
columns cost a batch nothing. Batches differ only in `selected_phenotypes`.

Run on 2026-10-12: 650 phenotypes, 58 of them sex-restricted from the
definitions table and none from the case sex ratio; `analysis_samples.txt`
md5 `1f389dd96aa9f709cfe4845f0e5c3eae`, identical to the pilot's.

### 1. Regenerate the phenotype table

In `notebooks/terra/tractor_01_prepare_inputs.ipynb`, change the selection cell
to the full list and the scale flags:

```python
N_PHENOTYPES = 0        # keep every phenotype that passes, not the top 10
```
```
--sex-specific restrict --phecode-sex-tsv <phecode definitions with a sex column>
--nesting leaves        # or roots; `keep` leaves parent and child as separate traits
--max-r2 0.5            # what All by All used; off by default
```

Without `--phecode-sex-tsv` the sex call falls back to the observed case sex
ratio at `--sex-case-fraction 0.99`, which is a guess. Source the definitions
table first.

**Write the outputs to a new prefix.** Do not overwrite
`tractor_mix_pilot/pheno_cov.tsv` — a running workflow still reads it.

**Then check `analysis_samples.txt` is unchanged.** `select_phenotypes.py`
derives the analysis cohort from covariate completeness *before* it selects
phenotypes, so going from 10 phecodes to the full list should leave it
byte-identical:

```bash
gsutil cat gs://BUCKET/tractor_mix_pilot/analysis_samples.txt | md5sum
md5sum <new out-dir>/analysis_samples.txt
```

If they match, the GRM and all 22 chromosome packs call-cache and only the
nulls and step 2 rerun. If they differ, every task in the run is invalidated —
find out why before submitting.

### 2. Generate the batch rows

```bash
python3 scripts/make_felix_run_table.py \
  --phenotypes <new out-dir>/selected_phenotypes.txt \
  --template felix/configs/terra/felix_runs.entity.tsv \
  --template-row felix.wgs.limited.v1 \
  --row-prefix felix.wgs650.limited \
  --batch-size 95 \
  --bucket gs://BUCKET \
  --phenotype-prefix felix/phenotypes/wgs650 \
  --pheno-cov gs://BUCKET/felix/phenotypes/wgs650/pheno_cov.tsv \
  --out-dir felix/configs/terra/wgs650
```

The realised list is 650 phenotypes (`--min-cases 100` on 9,775 participants),
which at `--batch-size 95` gives 7 rows of ~6,800 calls each instead of one
workflow of 46,242. Rows differ from the template only in the id,
`selected_phenotypes` and `run_description` — asserted in
`scripts/test_make_felix_run_table.py`. Run it again with
`--template-row felix.wgs.full.v1` for the full-covariate arm.

### 3. Upload and submit

Copy the per-batch phenotype files to the bucket (the script prints the
`gsutil` line), then upload membership and entity TSVs as in **Uploading**
above, and submit the same method config on each batch row.

Submit one batch first and let it finish. The GRM, the packs and the repeat
aggregation it builds are then cached for the other nine, so the remaining
batches are step 2 and summarization only.
