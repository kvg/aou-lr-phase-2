# Verily Workbench notebooks

Notebooks that run on Verily Workbench live here. The Workflows GUI is not
reliable. Prefer native Cloud Batch (`gcloud batch`) from the EH notebooks when
Cromwell is spotty; Cromwell `wb workflow` remains in the `_01` notebooks.
`_02` is the one-sample smoke; `_03` / `_04` are the production ledger.
[`batch_monitor.ipynb`](batch_monitor.ipynb) lists every native Batch job
(read-only). ExpansionHunter `_05` queues the next cohort wave and `_06`
draws a slide figure. Locityper `_02` / `_03` / `_04` use native Batch like EH.

Companion CLIs live in [`../../scripts/`](../../scripts/). Batch production
helpers: [`../../scripts/vwb_batch/`](../../scripts/vwb_batch/). Jobs use the
workspace PET SA and private VPC.

| Notebook | Use |
|---|---|
| [`locityper_00_prep_reference.ipynb`](locityper_00_prep_reference.ipynb) | Copy GATK `Homo_sapiens_assembly38`, Jellyfish 25-mers, one-locus smoke BED + toy `vcf_db`; peek at the v9 CRAM manifest |
| [`locityper_01_run_stream.ipynb`](locityper_01_run_stream.ipynb) | Cromwell: stage / register / submit `LocityperStream` |
| [`locityper_02_run_batch.ipynb`](locityper_02_run_batch.ipynb) | Native Cloud Batch **smoke**: one sample, print_reads minicram then Locityper (no dsub) |
| [`locityper_03_budget.ipynb`](locityper_03_budget.ipynb) | Pilot 20 samples, p90/p95 cost, write `keep.csv` to the GCS run ledger |
| [`locityper_04_dispatch.ipynb`](locityper_04_dispatch.ipynb) | Read keep list, shard, poll, retry missing/failed (native Batch, no DAG) |
| [`expansion_hunter_00_prep_reference.ipynb`](expansion_hunter_00_prep_reference.ipynb) | Stage GATK `Homo_sapiens_assembly38` + EH catalogs on the EH VM (no Jellyfish); peek at the v9 CRAM manifest |
| [`expansion_hunter_01_run.ipynb`](expansion_hunter_01_run.ipynb) | Cromwell: stage / register / submit `ExpansionHunterMinicram` |
| [`expansion_hunter_02_run_batch.ipynb`](expansion_hunter_02_run_batch.ipynb) | Native Cloud Batch **smoke**: one sample, minicram then ExpansionHunter (no dsub) |
| [`expansion_hunter_03_budget.ipynb`](expansion_hunter_03_budget.ipynb) | Pilot 20 samples, p90/p95 cost, write `keep.csv` to the GCS run ledger |
| [`expansion_hunter_04_dispatch.ipynb`](expansion_hunter_04_dispatch.ipynb) | Read keep list, shard, poll, retry missing/failed (native Batch, no DAG) |
| [`expansion_hunter_05_submit.ipynb`](expansion_hunter_05_submit.ipynb) | Run All: queue first attempts, then retries, for the EH keep list |
| [`expansion_hunter_06_progress.ipynb`](expansion_hunter_06_progress.ipynb) | Read-only EH progress figure (PNG/PDF) for a slide |
| [`batch_monitor.ipynb`](batch_monitor.ipynb) | Read-only status of all native Batch jobs + GCS run ledgers |
