# Verily Workbench notebooks

Notebooks that run on Verily Workbench live here. The Workflows GUI is not
reliable. Prefer native Cloud Batch (`gcloud batch`) from the EH `_02` notebook when
Cromwell is spotty; Cromwell `wb workflow` remains in the `_01` notebooks.
Locityper `_02` still uses dsub until that path is ported.

Companion CLIs live in [`../../scripts/`](../../scripts/). Batch jobs use the
workspace PET SA and private VPC.

| Notebook | Use |
|---|---|
| [`locityper_00_prep_reference.ipynb`](locityper_00_prep_reference.ipynb) | Copy GATK `Homo_sapiens_assembly38`, Jellyfish 25-mers, one-locus smoke BED + toy `vcf_db`; peek at the v9 CRAM manifest |
| [`locityper_01_run_stream.ipynb`](locityper_01_run_stream.ipynb) | Cromwell: stage / register / submit `LocityperStream` |
| [`locityper_02_run_batch.ipynb`](locityper_02_run_batch.ipynb) | Cloud Batch / dsub: `print_reads` minicram, then genotype (preferred when Cromwell is spotty) |
| [`expansion_hunter_00_prep_reference.ipynb`](expansion_hunter_00_prep_reference.ipynb) | Stage GATK `Homo_sapiens_assembly38` + EH catalogs on the EH VM (no Jellyfish); peek at the v9 CRAM manifest |
| [`expansion_hunter_01_run.ipynb`](expansion_hunter_01_run.ipynb) | Cromwell: stage / register / submit `ExpansionHunterMinicram` |
| [`expansion_hunter_02_run_batch.ipynb`](expansion_hunter_02_run_batch.ipynb) | Native Cloud Batch: `make_minicram_for_expansion_hunter` (no dsub) |
