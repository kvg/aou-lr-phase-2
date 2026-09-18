Adapted from a batch ExpansionHunter WDL that ran the
[bw2/ExpansionHunter](https://github.com/bw2/ExpansionHunter) fork in
`--analysis-mode seeking` directly against `gs://` WGS CRAMs (htslib +
`GCS_OAUTH_TOKEN`), writing a resume-able JSONL. That layout is expensive on
AoU v9 Nearline: every sample seeks the full CRAM.

The WDL here (`wdl/ExpansionHunterMinicram.wdl`) is a **Verily Workbench**
per-sample adaptation:

- One Cromwell job = one sample (same pattern as `LocityperStream`).
- `MakeMinicram` runs
  [str-analysis `make_minicram_for_expansion_hunter`](https://github.com/broadinstitute/str-analysis/blob/main/str_analysis/make_minicram_for_expansion_hunter.py)
  once: catalog `ReferenceRegion` + `OfftargetRegions`, then a mate-discovery
  pass, downloading each CRAM container **once** (GCS client, not `samtools
  view` / EH seeking).
- `ExpansionHunterGenotype` runs on that local minicram with
  `--analysis-mode optimized-streaming` (bw2 flag; falls back to `seeking` if
  the binary lacks it).
- Reference FASTA is a WDL input (GATK `Homo_sapiens_assembly38` already
  staged by the Locityper prep notebook). It is **not** baked into the image.
- Catalog JSON in `configs/` was converted from Python dict literals
  (`ast.literal_eval`) to real JSON. Smoke catalog is the first two loci.

`MakeMinicram` reuses the Locityper print_reads image
(`aou-locityper-print-reads:0.1.0`), which already has `str-analysis`. The EH
image is a slim debian build of the bw2 fork + htslib 1.22.
