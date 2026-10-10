# Manuscript figures

Figure-generation code for the AoU-LR Phase 2 manuscript. It lives in this repo
(version controlled) but reads inputs from, and writes outputs into, the
manuscript checkout.

Until 2026-10-10 this package sat at `aou-lr-phase-2-manuscript/scripts/mockup_figures/`,
where `.gitignore` excluded `scripts/`, so none of it was under version control.

## Where things go

| | Path |
|---|---|
| Code | `aou-lr-phase-2/manuscript_figures/` (this directory) |
| Manuscript root | sibling `aou-lr-phase-2-manuscript/`, or `$AOU_LR_MANUSCRIPT_ROOT` |
| Panel inputs | `<manuscript>/data/` |
| Draft renders | `<manuscript>/mockup_figures/` |
| Final plates | `<manuscript>/figures/` |

`<manuscript>/data/` holds participant-level files, including research IDs. It
is gitignored in the manuscript repo and must stay that way. Never copy its
contents into this repo.

## Running

From the repo root:

```bash
python -m manuscript_figures              # every figure in the FIGURES registry
python -m manuscript_figures fig1 fig3    # a subset
```

If the two repos are not siblings:

```bash
AOU_LR_MANUSCRIPT_ROOT=/path/to/aou-lr-phase-2-manuscript python -m manuscript_figures
```

## Manuscript plates

`python -m manuscript_figures` does **not** build Figure 2 — `fig2_options` is
absent from the `FIGURES` registry and is reached only through `render_main()`:

```bash
python -c "import manuscript_figures.fig2_options as m; m.render_main()"
```

That renders `option_b(rng=7, ratios=(0.48, 1.62), summary="new",
row_scales="symlog")` and copies the result to `<manuscript>/figures/fig2_catalog.{pdf,png}`,
which is what `aou-lr-phase-2.tex` includes.

| Manuscript figure | Built by |
|---|---|
| Fig. 1 `figures/fig1_cohort` | `fig1_null.render` (registry name `fig1_cohort`) |
| Fig. 2 `figures/fig2_catalog` | `fig2_options.render_main` |
| Figs. 3–6, supplementaries | `FIGURES` registry in `__init__.py` |

## Global ancestry

`fig2_options.GLOBAL_ANC` is the FLARE global-ancestry table, and it feeds both
the ancestry strip and `write_discovery_order()`, which writes
`<manuscript>/data/fig2_discovery_order.tsv` — the participant order the Terra
discovery job consumed. Changing `GLOBAL_ANC` therefore changes the discovery
order, and `data/fig2_panel_b_discovery.tsv` has to be regenerated to match.

## Conventions

Science 2-column: 7.2 in wide, Arial/Helvetica, 7 pt minimum in the file
(≈6.2 pt placed). Counts of 1–19 are disclosed as "<20", aggregated, or blanked.
See `style.py` and `scratch/HANDOFF_fig2.md` in the manuscript repo.
