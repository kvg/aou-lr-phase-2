# Optional landscape producers

These are **not** required for the manuscript site-count table or Ebert discovery
plots. Use them for Collins-style MEI / TR / remap panels.

Prefer the upstream Dockstore workflows (pin a commit + container digest):

| Layer | Dockstore / source |
|-------|--------------------|
| TR motifs | `talkowski-lab/lr-annotation/AnnotateIndelTRs:aou_scaling` |
| Remap | `talkowski-lab/lr-annotation/AnnotateRemap:aou_scaling` (`CallInsRemap`) |
| MEI | `talkowski-lab/lr-annotation/AnnotateL1MEAID:aou_scaling` |
| TSV→VCF transfer | `talkowski-lab/lr-annotation/AnnotateVcf:main` |
| External AF | local `scratch/yulia/AnnotateExternalAFs.7.wdl` (refactor before production) |

After each producer emits a TSV, normalize TR rows with:

```bash
python3 ../scripts/format_optional_tr_tsv.py \
  --in-tsv tr_raw.tsv \
  --out-tsv tr_for_annotate.tsv
```

Then transfer INFO fields with `bcftools annotate` using the schema in
`schemas/info_schema.json` (`optional_info_fields`).
