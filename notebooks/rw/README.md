# Verily Workbench notebooks

Notebooks that run on Verily Workbench live here. The Workflows GUI is not
reliable; submit WDLs with the `wb workflow` CLI from these notebooks.

Companion CLIs live in [`../../scripts/`](../../scripts/).

| Notebook | Use |
|---|---|
| [`locityper_00_prep_reference.ipynb`](locityper_00_prep_reference.ipynb) | Copy GATK `Homo_sapiens_assembly38`, Jellyfish 25-mers, one-locus smoke BED + toy `vcf_db`; peek at the v9 CRAM manifest |
| [`locityper_01_run_stream.ipynb`](locityper_01_run_stream.ipynb) | Stage / register / submit `LocityperStream` (`print_reads` minicram, then genotype) |
| [`expansion_hunter_01_run.ipynb`](expansion_hunter_01_run.ipynb) | Stage catalogs / register / submit `ExpansionHunterMinicram` (`make_minicram_for_expansion_hunter`, then genotype) |
