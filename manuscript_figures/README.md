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

## Merging Figures 1 and 2

`fig12_combined` composes the surviving panels of both onto one plate, for the
reviewer request to spend less page on design and catalog. Fig. 1B (PC scatter)
and Fig. 1F (pedigrees) move to the supplement; everything else is kept.

```bash
python -c "import manuscript_figures.fig12_combined as m; m.render_all()"
python -c "import manuscript_figures.fig12_combined as m; m.render('band', 'shares')"
```

`layout` is `band` (cohort row 1.45 in, four loci in the allelic series),
`tall` (cohort row 1.85 in, three loci), or `compact` (cohort row cut to
1.30 in and the catalog gap tightened, with the space handed to the length
spectrum: 0.96 in of axes against 0.66 in before the rebalance, for a 6.4 in
plate). `center` chooses how sequencing center is shown —
`named`, `omit`, `shares` (equal-length composition bars, alphabetical), or
`coverage` (a median/IQR rail under the coverage histogram). Four cohort panels
leave the UpSet too narrow, so `named` and `shares` swap in
`_panel_molecular_compact`, one bar per assay.

The `coverage` rail carries Phase 1 first, in the Phase 1 color and labelled
`HudsonAlpha (Phase 1)` — 1,023 of its 1,027 participants were sequenced there
— then the four PacBio centers with at least 20 participants, alphabetically.
JHU never appears: it contributed no PacBio. ONT came from Baylor (1,022), JHU
(710), UW (396) and Broad (17), so no single rail can be labelled for it.

Letters: the cohort row, the spectrum and the discovery panel share the figure
margin as their letter column, and `_align_letter_to_legend` drops the spectrum
letter onto the baseline of the context legend above it. The cohort row is
placed with `inline_letters=True`, which sets its letters alongside the top of
the row's ink rather than on a band above it — every legend there is anchored
to its own axes, so the letter column stays clear of them and the row reclaims
the letter band, `_LETTER_GAP + _LETTER_SIZE/72` = 0.165 in. Measured: the
catalog region grows by exactly that (3.189 → 3.354 in), of which 0.112 in
reaches the two panels' axes (spectrum 0.964 → 1.014, discovery 1.200 → 1.262)
and 0.053 in goes into the `hspace` gap between them, which scales with row
height. `cat_gap` per layout
sets the space between the cohort row and the spectrum.

`assays="bar_upset"` (the shipped choice) carries ONT as a platform bar over
the assay UpSet. The whole releasable Phase 2 cohort is a three-segment bar in
the UpSet's bar column, on the UpSet's own x scale, so one axis at the foot of
the panel serves both. Two square brackets sit above the bar — ONT on the upper
level, PacBio on the lower — and they overlap over the 876 participants
sequenced on both, which is the point: the platforms are not a partition. A
curly brace under the bar funnels the PacBio span into the UpSet below, a
different glyph for a different relationship (square = "this span is that
platform", curly = "this span breaks out below"). `brackets=False` gives the
plain bar-over-UpSet without either annotation.

`counts="study"` is the shipped framing, and it drops the platform taxonomy
for three statements: a bracket over the whole bar labelled
`PacBio+ONT (13,530)`, a leader onto the dark segment labelled `Both (876)`,
and a bracket under the PacBio span labelled `This study (12,261)` whose stem
runs down into the UpSet. The message is that 12,261 participants enter this
analysis, some of them with ONT alongside, and the remaining ONT data belongs
to a companion paper — not that the bar is a platform partition. `connector`
is ignored in this mode; the lower bracket is the connector. The header is
wider than the bar, which costs nothing: measured, the panel is the same width
with and without it, because the UpSet's count labels already set the row's
right-hand ink and the gap between the matrix and bar insets absorbs the
overhang on the left. Its stem sits at 0.22 of the span rather than the middle
so it passes clear of the right-aligned overlap label on the line below.

`connector` and `counts` are the two knobs on the platform-labelled element, and they are
coupled: whatever sits directly under the bar is what the connector appears to
attach to, so a counts line between the bar and the brace makes the brace read
as if it brackets the text. `counts="under_connector"` and
`counts="in_brackets"` both avoid that; `counts="under_bar"` only works with
`connector="none"` or `"tint"`. Tried and set aside: `connector="bracket"`
(a square brace below mirrors the platform brackets above and reads as a box
around the bar) and `connector="tint"` (a pale band is too heavy and its
bottom edge competes with the UpSet). `counts="in_brackets"` folds the totals
into two-line platform labels, which is the tidiest option but has no room for
the 876 overlap — every spot between the two bracket levels is crossed by one
of their tips — so the 876 goes in the caption. `tip_label` instead states what
the brace delivers at the point where it arrives: 12,261, the denominator the
four UpSet rows sum to and the same n panel E carries. When it is set the
PacBio bracket drops to a bare `PacBio` rather than repeating that number, and
the one-line label is lifted by one text height so it stays on the same first
line as the two-line `ONT` / `2,145`. Put the tip label to the LEFT of the
point: to the right it runs past the end of the 13,530 bar, and `_place_row`
sizes panels from their tight bboxes, so the row answers by narrowing panel C.

`legend_lift=N` appends N blank legend rows at the bottom. Every row legend on
this plate is bottom-anchored at `_ROW_LEGEND_Y`, so adding rows is what moves
the visible lines up — exact by construction, with no measured font constants —
and `legend_lift=1` puts panel C's two lines on panel B's `Phase 1` and
`Mid-pass` rows. The band that frees underneath is spent by setting the top of
`ylim` inside the top label line, letting it overflow the axes: text is not
clipped, and panel A's four-row legend still sets the row's top, so the plate
does not grow. Clearance above and below the bar is equalised through
`BAR_GAP = DROP - TIP - BAR_H / 2`, with the brace's slot deepened by its
`0.40 * CONN_H` end curl so the two gaps look the same.

Six other ONT placements were built and rejected — fourth UpSet column,
coverage-panel depth overlay, two-tone row under the matrix, two stacked
blocks, pyramid outline, 5mC dropped in favour of PacBio/ONT columns, and a
braced bar on a separate inset axis (`assays="bracket"`). They all still
render, for the record.

`ont_row=True` puts ONT in the assay panel as a single
bar below a rule, outside the UpSet's partition of the PacBio cohort: 876 with
PacBio in the ONT color, 1,269 ONT-only in the pale ONT color, labelled
`876 + 1,269`. Deliberately not a fourth matrix column — the 876 are already
inside the rows above, so a column would double-count them, and splitting
every row by ONT is the seven-row version that was too dense. The caption has
to say that the 876 are also counted above, and which segment is which.

`ont_coverage=True` instead overlays ONT depth on the coverage panel as an
outline in
the ONT color (n = 2,145: 876 PacBio+ONT participants carrying `ont_coverage`
on their PacBio row, plus 1,269 ONT-only participants read from the covariate
table, which `_load_model_rows` excludes). Median 34x, so it sits on the
high-pass PacBio mode; no mean rule is drawn because the two would print as
one, and the 11 samples above 52x go into the last bin. The x label becomes
`Coverage (x)` and the center rail is tagged `PacBio, by center`.

`assays="ont"` instead puts ONT in the UpSet as a fourth column with an
ONT-only row. Tried and rejected: seven rows make the matrix hard to read at
this panel height. It also counts `In QTL analysis` only for participants who
actually have Olink or RNA-seq, where `fig1_null._panel_molecular_data` treats
a missing assay as satisfied and so marks the 5mC-only group as in-analysis —
that difference is still unresolved in Fig. 1E.

Panel C's shading is a release-versus-freeze statement, not QTL membership:
data that arrived after a companion study's freeze is in the release but not
in that study. `assays="upset_freeze"` (the shipped choice) keeps
`fig1_null._panel_molecular_data`'s arithmetic and relabels the key to
`In companion analyses` / `Release only`. The vacuous test for a missing assay
is correct under that reading — nothing of theirs was left out — but it does
make 5mC-only participants dark on the assumption that methylation has no
freeze gap, which no `in_meqtl` list exists to check. In practice the whole
pale area is Olink: 390 of 9,243 Olink participants are absent from the pQTL
list (344 in the all-three row, 46 in 5mC+Olink), and RNA-seq has no gap at
all (8,327 assayed, all in sQTL). `upset_any` and `upset_plain` are the
alternatives that were tried and set aside.

`summary=False` drops the repeat-allele carrier histogram and moves its grey /
blue key onto the allelic series.

Current build for the merged plate:

```bash
python -c "import manuscript_figures.fig12_combined as m; \
  m.render('compact', 'coverage', assays='bar_upset', counts='study', \
    summary=False, stem='fig12_v25')"
```

Plates land in `<manuscript>/mockup_figures/` at 7.2 in wide and 6.6–8.4 in tall
depending on the options.

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
