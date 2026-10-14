# `felix_prod_set`: one row per phenotype

Reshapes the genome-wide scan from 11 batch rows of ~87 phenotypes into 951 rows
of one, so a failure is isolated to a single phenotype and Terra's **Relaunch
failures** re-runs only the rows that failed. Submit on the *set*, not row by
row: that is one submission containing 951 workflows.

## Why the inputs JSON does not change

`FelixGenome.inputs.json` references only

- `${this.<column>}` for 24 set columns, and
- `${this.flare_lai_prods.<column>}` for the four per-chromosome inputs
  (`chrom`, `gt_vcf`, `anc_vcf`, `annotated_vcf`),

plus 9 literal `gs://` script paths and 8 literal numeric constants. **Nothing in
it names the set type.** The new table keeps every column name and keeps
`flare_lai_prod` as the member type, so every `${this.*}` resolves unchanged and
the file is byte-identical to the live one. The only change is the method
config's **root entity type**: `flare_lai_prod_set` -> `felix_prod_set`.

### RESOLVED: `pack_cpu = 1`, `pack_memory_gb = 6`

`PackFelixla` declares `cpu` and `memory_gb` as task *inputs*, not merely as
runtime attributes, so their values are hashed: get them wrong and all 22 packs
re-run in every one of the 951 workflows.

Read off the Inputs pane of a succeeded `Pack` shard (chr1) from the run whose
cache is being inherited: `cpu 1`, `memory_gb 6`, `disk_gb_floor 100`,
`disk_gb_multiplier 3.5`, `preemptible 0`, `num_ancs 5`, `disk_gb 145`.

That matches the config exported from Terra and the working copy of
`FelixGenomeTable.inputs.json`, and contradicts the same file as committed on
`main`, which still carries `pack_cpu = ${4}` / `pack_memory_gb = ${64}` from
before the resource rework. The project's resubmission note recording the
committed values as what batch01 ran is wrong for this run. Regenerating the
config from `main` instead of the working copy would re-pack 20,922 chromosome
shards.

Every task *command* is byte-identical between `main` and the working copy --
the rework changed resource declarations only -- so nothing else in the
uncommitted diff affects the cache.

### Terra will not accept a separate `felix_prod_set`

The member type is derived from the set name: a `felix_prod_set` must contain
`felix_prod` entities. Uploading a membership TSV whose second column is
`flare_lai_prod` fails with

    Invalid membership TSV. Second column header should be felix_prod

A separate table therefore needs a `felix_prod` member table -- 22 rows
duplicating `flare_lai_prod` -- and every `${this.flare_lai_prods.*}` in the
inputs JSON rewritten to `${this.felix_prods.*}`. The duplicate member rows
would have to carry byte-identical paths or every cached call is lost.

So the rows go into the existing `flare_lai_prod_set` table instead. They are
distinguishable by their `felix.prod.` prefix, the 16 existing rows are
untouched (no id collisions), and the root entity type does not change, so the
EXISTING method config runs them with no edit at all. Drop `--set-type`.

## Generating the table

Clone from `felix.wgs951.full.batch01` -- **not** from `felix.wgs.full.v1`, which
points at `tractor_mix_pilot/` paths. `analysis_samples` is an input to
`PackFelixla` and `MakeGRM`, so a different sample list re-packs 22 chromosomes
in every one of the 951 workflows.

    python3 scripts/make_felix_run_table.py \
      --phenotypes selected_phenotypes.txt \
      --template flare_lai_prod_set_entity.tsv \
      --template-row felix.wgs951.full.batch01 \
      --row-prefix felix.prod \
      --name-rows-by-phenotype \
      --batch-size 1 \
      --bucket "${WORKSPACE_BUCKET}" \
      --phenotype-prefix felix/phenotypes/prod \
      --out-dir prod_table

`--template` is the entity TSV exported from Terra. Only `run_description` and
`selected_phenotypes` differ from the template row; every cache-critical column
(`pheno_cov`, `analysis_samples`, `covariate_columns`, `grm_vcfs`, `docker`, the
`ru_*` parameters) is copied verbatim, which is what keeps `Pack`, `MakeGRM`,
`RuLoci`, `RuVcf`, `BaVcf` and the nulls cache-hitting.

Upload `prod_table/phenotypes/*.txt` to
`${WORKSPACE_BUCKET}/felix/phenotypes/prod/`, then the **membership TSV first**
and the entity TSV second. That order relies only on upsert behaviour already
verified in this workspace (an attribute TSV naming a subset of columns leaves
the others alone); the reverse would depend on a membership upload preserving
attributes, which has not been checked here.

## Verify on one row before submitting 951

Submit `felix.prod.CV_401` alone and read the **Call Caching Result** column:

| call | expected |
|---|---|
| `MakeGRM`, `Pack`, `RuLoci`, `RuVcf`, `BaVcf` | Cache Hit |
| `Null` | Cache Hit (the shared `pheno_cov` is unchanged) |
| `Step2`, `RuStep2`, `BaStep2` | Cache Hit for a phenotype that already completed |

A Cache Miss on `Pack` or `MakeGRM` means a column drifted from the template and
the whole reshape is mispriced -- stop and diff the row against
`felix.wgs951.full.batch01` before going further.

One thing left to confirm on first upload, which cannot be checked from outside
Terra: whether AoU's Terra caps workflows per submission at under 951.

The set-type question is settled -- see the rejection above. These rows are
`flare_lai_prod_set` rows, so `${this.flare_lai_prods.*}` resolves exactly as it
does for `felix.wgs951.full.batch01` and the method config needs no edit.

## Cost of the reshape

Every cache hit still *copies* its outputs into the new workflow's directory
(`cacheCopy/`). Measured: `felixla_packed.tar.gz` is 1.62 GB per chromosome and
the sparse GRM 0.18 GB, so the upstream copy is ~36 GB per workflow before the
repeat VCFs -- roughly 76 TB across 951 rows against 0.9 TB across 11. That is
bucket storage, not compute, and it is reclaimable afterwards.

## The four phenotypes that need a restricted `pheno_cov`

Computed from the local phenotype and covariate tables against the real 951-name
list; all four are in it, and all four are labelled `Both` by PhecodeX:

| row | cases | female | male | separation |
|---|---|---|---|---|
| `felix.prod.CA_142_2` (Abnormal mammogram) | 360 | 360 | 0 | complete |
| `felix.prod.GU_614_55` (Candidiasis of vulva and vagina) | 154 | 154 | 0 | complete |
| `felix.prod.CA_142_21` (Mammographic microcalcification) | 79 | 79 | 0 | complete |
| `felix.prod.CA_105` (Malignant neoplasm of the breast) | 229 | 227 | 2 | quasi |

The three with zero cases in one sex must diverge in step 1. `CA_105` is
different in kind: male breast cancer is real, so the `Both` label is correct and
2 cases simply cannot support a sex term -- it may converge to a strained fit and
return wrong p-values with no crash, which is why crash-as-signal is not enough
to find it.

Generate each restricted file with
`scripts/make_restricted_pheno_cov.py --phenotype <row> --restrict-to Female`
and upload to `felix/phenotypes/prod/pheno_cov.<row>.tsv`, which is where those
four rows already point.

## Fixing one phenotype

Point that row's `pheno_cov` at its own file:

    --set pheno_cov=gs://.../felix/phenotypes/prod/pheno_cov.<phenotype>.tsv

Make that file a **copy of the shared `pheno_cov.tsv`** with the excluded
stratum set to `NA` in that one phenotype's column. `FitFelixNull` takes
`pheno_cov` plus a `phenotype` name, so the file can keep all 951 columns.

Do **not** subset `analysis_samples` to achieve the same thing: SAIGE drops
NA-phenotype rows itself, and changing that file re-runs `Pack` and `MakeGRM`.
