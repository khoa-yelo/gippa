# M12-05c: structured plot guides in the HTML report

Milestone M12 (1.1.0) · depends on: M12-04 PR 2 (captions, `gippa.report.content.CAPTIONS`),
M12-04 PR 4 (V03 index table, empty states), **M12-05a** (V13 lists every novel ORF candidate by
q; only the V13 guide depends on it) · does not touch: M12-05b (`novel_orfs.tsv` gene columns) ·
estimate ≈ 550 lines over 2 PRs (tests ≈ 260, implementation ≈ 290) · roles: spec-writer (this
file), test-author (`tests/unit/report/test_m12_05c_*.py` and updates to the M12-04 PR 2 caption
tests), implementer (`src/gippa/report/content.py`, `src/gippa/report/build.py`,
`src/gippa/report/compare.py` or wherever the compare slot is built, the two templates and
`report.css`) · source: Khoa's review comment on the 1.1.0rc1 report (2026-11-09): "The
instruction on how to interpret each plot should be more structured (e.g. Instruction:
Question: .... Panel....".

## 1. Current state (main at e6f4aa3c)

- `gippa.report.content.CAPTIONS` maps the 16 plot IDs (V01–V03, V05–V08, V10–V18; V04 and V09
  do not exist) to one string: the GIPPA_SPEC.md §10.2 Question, then one run-on "how to read"
  sentence. `content.caption(pid, stem)` returns it; `INDEX_CAPTIONS` gives the V13 index instance
  its own sentence.
- `report/build.py::_plot_slot` puts `content.caption(pid, inst.stem)` on every drawn instance.
  `report.html.j2` renders `<figcaption class="links" data-stem="…">{{ caption }} SVG PNG PDF
  data TSV</figcaption>`. Every V03 locus view repeats the full V03 caption (up to 50 times).
- `compare.html.j2` renders the V16 figcaption as the bare stem (`V16_compare SVG …`): the compare
  report has **no** guide at all.
- The guide is not in the PNG/SVG/PDF files, and there is no `plots/README`; `OUT/README.txt`
  has one line for `plots/`.

Captions that no longer match the plot (checked against the route A–E SVGs of example run
1945698 and the `gippa.viz` module docstrings), fixed by the guides in § 4:

| Plot | Current caption says | The plot draws |
|---|---|---|
| V01 | "features" ring; "transcription units with promoters and terminators" as one ring | a "regions" ring (prophage, plasmid, BGC, CRISPR, AMR arcs); TU starts and promoter/terminator ticks are two rings; rings are numbered and keyed; rings with nothing to draw are left out |
| V02 | points, LOESS, islands | also a rug of undetected genes, origin/terminus lines, and a detection strip when the run has no abundance |
| V03 | "flagged operon or novel ORF"; no order | novel ORFs are drawn by V13, not V03; the order takes the best locus of each kind in turn; undetected genes are hatched white |
| V07 | "observed minus expected number of detected genes" | the x axis is in **percentage points** of the group's genes |
| V10 | "join-count statistic" | panel (b) is detection: the join statistic as a Pearson residual |
| V12 | "each top site" | at most 15 sites; bars are also stacked by flag |
| V15 | four panels always | (b) and (c) are hidden under a banner when the model has no predictive power |
| V17 | "Panels show …" five items, no letters | five lettered panels (a)–(e) plus a warning list |

## 2. Decisions

1. **One labelled structure for every plot**, in this order:

   | Field (label shown) | Required | Content |
   |---|---|---|
   | **Question** | yes | the GIPPA_SPEC.md §10.2 Question, verbatim |
   | **Panels** | yes | one item per panel: what is drawn, the axes and the encodings. Lettered panels start with the letter as the plot prints it, "(a) ", "(b) ", …; a one-panel plot has one item with no letter |
   | **How to read** | yes | what a pattern means; 1–3 sentences |
   | **Caveats** | no | when not to trust the plot; 1–2 sentences; omitted (no empty `dt`) when the guide has none |
   | **Data** | yes | generated, not written by hand: a link to the instance's `plots/data/<stem>.tsv` (link text = that path), then "figure files:" and the SVG/PNG/PDF links already built in `_plot_slot` |

   Lead decision needed only on the labels (Q1); the field set and order are fixed here.
2. **Markup: a description list inside the `<figcaption>`.** Valid HTML (a `figcaption` holds
   flow content), readable without CSS (browsers indent `dd` by default), and accessible: the
   `figcaption` stays the figure's accessible name and `dt`/`dd` give screen readers the field
   names. Panels are a `<ul>` inside their `dd`. See § 3.
3. **Once per plot, not once per instance.** The first drawn instance of a plot carries the full
   guide with `id="guide-<pid>"`. Every later instance of the same plot (V03 locus views, V13
   locus views) carries a short guide: Question, "How to read" = "As for the first <name> figure
   (see the guide above)." with a link to `#guide-<pid>`, and its own Data. The V13 index instance
   keeps its index sentence as its "How to read" (§ 4, V13).
4. **Same guide in every report.** The guide text is static: `redact=True` and the size-budget
   `linked` mode produce the same `<dl>` bytes as the default report. The compare report's V16
   figure gets the full V16 guide, and its `figcaption` gets `data-stem="V16_compare"`.
5. **Images and READMEs do not carry the guide.** PNG/SVG/PDF files are unchanged (titles and
   panel letters already in them are what the guide refers to), so no pytest-mpl baseline moves.
   No `plots/README` is added (Q3). `OUT/README.txt` is unchanged.
6. **Empty states and skipped plots are unchanged.** A plot shown only as a `.notice` or a
   `.note` has no figure and no guide (Q4).
7. **Code shape (the implementer may rename privately).** `CAPTIONS` and `caption()` are replaced
   by `GUIDES: Mapping[str, Guide]`, `Guide` a frozen dataclass `(question: str, panels:
   tuple[str, ...], read: str, caveats: str | None = None)`; `question(pid)` stays and returns
   `GUIDES[pid].question`. `_plot_slot` passes `guide` (the `Guide`), `guide_full: bool` and the
   links; the template renders them. Wording rules of M12-04 AC 2.1 apply to every field (no
   `ADR-`, `§`, `GIPPA_SPEC`, `docs/methods`, no flag-vocabulary token).

## 3. Markup

Full guide (first instance; whitespace for reading only):

```html
<figure>
  …svg or plotly block…
  <figcaption class="links guide" data-stem="V02_abundance_landscape">
    <dl class="guide" id="guide-V02" data-plot="V02">
      <dt>Question</dt>
      <dd>Is abundance clustered along the chromosome?</dd>
      <dt>Panels</dt>
      <dd><ul><li>One panel per replicon: …</li></ul></dd>
      <dt>How to read</dt>
      <dd>A shaded band is …</dd>
      <dt>Caveats</dt>
      <dd>The trend line is …</dd>
      <dt>Data</dt>
      <dd><a href="../plots/data/V02_abundance_landscape.tsv">plots/data/V02_abundance_landscape.tsv</a>;
          figure files: <a href="../plots/V02_abundance_landscape.svg">SVG</a>
          <a href="../plots/V02_abundance_landscape.png">PNG</a>
          <a href="../plots/V02_abundance_landscape.pdf">PDF</a></dd>
    </dl>
  </figcaption>
</figure>
```

Short guide (later instances):

```html
<figcaption class="links guide" data-stem="V03_locus_bgc_2">
  <dl class="guide guide-short" data-plot="V03">
    <dt>Question</dt><dd>What happens at this prophage, BGC, island, flagged operon or novel ORF?</dd>
    <dt>How to read</dt><dd>As for the first locus view (<a href="#guide-V03">see the guide above</a>).</dd>
    <dt>Data</dt><dd><a href="../plots/data/V03_locus_bgc_2.tsv">plots/data/V03_locus_bgc_2.tsv</a>; figure files: …</dd>
  </dl>
</figcaption>
```

CSS (inline `report.css`, optional polish only): `dl.guide` as a two-column grid
(`dt` bold, right column `dd`), `dd ul` without bullets indent. Without CSS the page must still
read top to bottom as "Question / text / Panels / list / …".

## 4. The guides (verbatim text for `GUIDES`)

`·` separates the items of Panels. Text in *italics* is a note to the implementer, not guide text.

### V01 circos
- **Question:** Where in the genome is the proteome?
- **Panels:** One circle; each replicon is a sector with ticks in kb. Numbered rings, outer to
  inner, as listed in the ring key: regions (prophage, plasmid, BGC, CRISPR and AMR arcs); genes
  on the + strand, then the − strand, coloured by detection status; log2 abundance as a heat ring
  (white: no value or imputed); peptide density per 10-kb bin; transcription-unit starts;
  promoter and terminator ticks; abundance islands; cumulative GC skew; variant sites and
  novel-ORF candidates. Radial lines mark the origin and terminus when they are known.
- **How to read:** Arcs where the gene rings are mostly "detected" and the abundance and peptide
  rings are dense are the parts of the genome the proteome covers. A region arc over undetected
  genes and no peptides is a region the sample did not detect; V07 tests whether that is more
  than detectability predicts.
- **Caveats:** A ring with nothing to draw is left out, and the title names the missing layers.
  The circle is descriptive; the tests are in V02, V07 and the Genome context tables.

### V02 abundance_landscape
- **Question:** Is abundance clustered along the chromosome?
- **Panels:** One panel per replicon: x is the position (kb), y the log2 abundance of each
  detected gene with a measured value (one point per gene); undetected genes are ticks along the
  bottom; the line is a LOESS trend; shaded bands are abundance islands; vertical lines mark the
  origin and terminus when they are known. With no measured abundance, each panel is one strip of
  ticks coloured by detection status.
- **How to read:** A shaded band is a run of neighbouring genes more (or less) abundant than the
  rest of the replicon, beyond chance at a family-wise error of 0.05. The title gives the number
  of islands and the lag-1 Moran's I with its q. A wavy trend line without bands is not evidence
  of clustering.
- **Caveats:** The trend line is a visual guide and is not tested. The island test shuffles whole
  transcription units, so operon partners looking alike do not by themselves make an island.

### V03 locus_<id>
- **Question:** What happens at this prophage, BGC, island, flagged operon or novel ORF?
- **Panels:** One locus per figure, with flanks; x is the position (kb). Gene arrows point in
  the direction of transcription, coloured by log2 abundance on one scale for the whole run;
  undetected genes are white with a grey outline and hatching. Mapped peptides are stacked above
  the arrows. Promoters are bent arrows, terminators lollipops; the region or island is shaded. A
  dashed line marks the origin when the window crosses it.
- **How to read:** The title gives the locus, its detected share against the share detectability
  predicts, and the q. Look at which genes carry peptides and how abundant they are against the
  flanks. The views follow the "Locus views, in display order" table: the best locus of each kind
  in turn, then the second best, and so on, at most 50.
- **Caveats:** Novel-ORF candidates are drawn under Novel ORFs (V13), not here. A gene supported
  by a single peptide is marked, and its abundance is not estimated.

### V05 mapping_sankey
- **Question:** Where did my sequences go?
- **Panels:** One flow diagram, left to right: input class (peptide or protein), then the
  matching route (exact, tolerant, miniprot, tblastn or unmapped), then the outcome (placed,
  contaminant, rejected by FDR or unexplained). Band widths are proportional to the number of
  distinct sequences, and the labels give the counts. When every sequence takes one route to one
  outcome, the diagram is a single bar.
- **How to read:** Most of the flow should end in "placed". A large unmapped, unexplained or
  rejected share means the genome, the search database and the input do not match well; V06
  checks the tolerant route, and the unmapped table lists the sequences.
- **Caveats:** Each sequence counts once: one placed by one route and unmapped by another counts
  as placed.

### V06 mapping_fdr
- **Question:** Can I trust tolerant hits and de novo peptides?
- **Panels:** (a) targets vs. decoys: target and decoy placements per peptide length class and
  number of differences, on a log scale (a zero count sits at the floor). · (b) stratum FDR: the
  estimated false discovery rate per number of differences, with the acceptance level dashed; a
  stratum with no targets is NA. · (c) de novo precision: the share of correct de novo peptides
  per score bin, with its confidence interval; drawn only for de novo input.
- **How to read:** Where decoys come close to targets in (a), those placements are mostly chance;
  (b) turns this into an FDR, and only strata at or below the dashed line are accepted. Bins with
  low precision in (c) add errors to every later count. When every peptide matched exactly, the
  plot is one line saying so.
- **Caveats:** Without reference peptides, de novo precision is not measured and an assumed value
  is used; panel (c) names it.

### V07 obs_vs_exp
- **Question:** Which regions and categories are detected less than detectability predicts?
- **Panels:** One row per group, in blocks: replicon, feature class, COG category, region. The
  point is observed minus expected detected genes in percentage points of the group's genes, the
  bar its 95 % confidence interval; the group size n is in the row label and the q at the right.
- **How to read:** A row left of zero has fewer detected genes than its genes' detectability
  predicts, a row right of zero more. Only rows whose q is at or below the level in the title are
  findings; the title counts them.
- **Caveats:** The expected counts come from the detectability model, so check its reliability in
  V18 (a) first: where the model is off, every group can drift.

### V08 cog_abundance
- **Question:** Which functions dominate abundance?
- **Panels:** One column per COG category, highest median abundance first and "none" (no COG
  annotation) last: each measured gene's log2 abundance is a point, or a violin with its median
  above 50 genes; the label above each column is detected / total genes. With no measured
  abundance, each column shows the detected fraction instead.
- **How to read:** Categories to the left hold the most abundant proteins. Compare the medians,
  not the number of points.
- **Caveats:** Imputed values are not drawn. The plot is descriptive; the tested comparison of
  functions is V14.

### V10 correlogram
- **Question:** Is there structure beyond operons?
- **Panels:** (a) Moran's I of abundance against the lag (genes apart, 1 to 20) per replicon,
  with two shaded 95 % null envelopes: gene shuffles (grey) and transcription-unit block shuffles
  (blue). · (b) the same for detection: the join statistic of neighbouring pairs in which both
  genes are detected, as a Pearson residual, against the same two envelopes.
- **How to read:** A filled marker lies outside the blue envelope with q at or below the level:
  genes that far apart look alike beyond what operons explain. Open markers inside the blue
  envelope at every lag mean operons explain the similarity. The title gives the lag-1 result.
- **Caveats:** The grey envelope ignores operons and is too narrow; judge by the blue one. Lags
  with too few gene pairs are not drawn.

### V11 operon_pairs
- **Question:** Do operon partners look alike, and which other neighbours do too?
- **Panels:** (a) |Δ| by pair class: histograms of the absolute log2 abundance difference of
  adjacent genes, for same-operon (OP), non-operon (NO), convergent (CV) and sandwiched (SW)
  pairs, with n in the legend. · (b) the local FDR of each pair against |Δ|, with the flag level
  dashed. · (c) partners: upstream against downstream log2 abundance of OP and NO pairs, coloured
  by class; flagged pairs are larger outlined markers, one shape per flag.
- **How to read:** OP pairs at smaller |Δ| than NO pairs in (a) mean operon partners are made at
  similar levels, as expected. Pairs below the dashed line in (b) are flagged, for example an
  operon pair whose partners differ or a non-operon pair that behaves like one. In (c), points
  near the diagonal look alike.
- **Caveats:** Pairs with an imputed abundance are left out of (a) and (c). A flag is a lead to
  check in the locus view (V03), not a confirmed promoter or terminator.

### V12 variants
- **Question:** Where do peptides support ALT residues?
- **Panels:** (a) variant sites by call (ALT or mixed), one bar per confidence level, stacked by
  flag (a site with several flags counts in each); the total is printed above each bar. · (b)
  alignments of the top ALT sites (at most 15): the reference protein window, then the covering
  peptide with the ALT residue in bold colour, and the number of covering peptides.
- **How to read:** ALT means the peptides show only the alternative residue; mixed means both
  residues are seen. High-confidence ALT sites without flags are the strongest. In (b), the
  peptide should match the reference everywhere except the highlighted residue.
- **Caveats:** The class error rate is pooled over every site of a class, not per site; a rate of
  exactly 1 (marked *) means it was not measured. I and L have the same mass and are written I/L.

### V13 novel_orfs (text depends on M12-05a)
- **Question:** Which new genes or upstream extensions are supported?
- **Panels:** (a) candidates per class (for example N-terminal extension, intergenic, small ORF),
  stacked by status: called, not called or unestimated. · (b) every candidate, smallest q first
  (unestimated last), each in its genomic neighbourhood: the candidate, the nearby gene arrows,
  its peptides and the promoter and terminator marks, with its class, status and q in the title.
- **How to read:** A called candidate passed its class FDR; its q is the lowest FDR at which it
  would be called. Read (b) from the top: the first candidates are the best supported, and a
  candidate with several peptides, a promoter upstream and a start codon with a ribosome-binding
  site is the most convincing.
- **Caveats:** When a class has too few decoys its q is unestimated and nothing in it is called;
  the candidates are still drawn. A single peptide is weak evidence on its own.
- *Index instance (`V13_novel_orfs_locus`), How to read:* "This index lists every candidate drawn
  in (b), in the same order." *(Replaces `INDEX_CAPTIONS["V13"]`.)*
- *If M12-05c merges before M12-05a, ship this interim (b) and index text, and M12-05a switches
  it:* (b) "the top 20 called candidates by rank, each in its genomic neighbourhood …" (rest as
  above); How to read second sentence "The views follow the index, best rank first."; index "This
  index lists the top called candidates drawn as locus views, by rank." *The V13 guide must match
  whatever panel (b) M12-05a finally draws (AC 7); the M12-05a author owns that edit.*

### V14 enrichment
- **Question:** Which functions are high or low?
- **Panels:** (a) GSEA running sums for the top six gene sets: the running enrichment score along
  the genes ranked by abundance, the set's genes as ticks and the leading edge shaded. · (b) each
  set's normalised enrichment score (NES) against −log10 q, coloured by collection; flagged sets
  are hollow, redundant sets small and grey, and the ten most significant are labelled. · (c)
  detection over-representation: observed over expected detected genes per category on a log2
  axis, coloured by direction, with −log10 q in a strip beside it and the level line.
- **How to read:** In (b), sets to the right are more abundant than average and sets to the left
  less; higher means a smaller q. A running sum in (a) that peaks early is driven by the most
  abundant genes. In (c), categories to the right of 1 are detected more often than expected.
- **Caveats:** A set marked "not calibrated" (dashed outline in (a), hollow in (b)) is clustered
  in the genome or tied to detection, so its q is too small. Panel (c) has no confidence
  intervals.

### V15 determinants
- **Question:** What explains abundance?
- **Panels:** (a) block ΔR²: the drop in cross-validated R² when each feature block
  (transcription, translation, codon, protein, context) is left out, for the elastic net (the
  main model) and LightGBM (a check), with intervals. · (b) partial dependence of abundance on
  the codon adaptation index, the tRNA adaptation index, the OSTIR initiation rate and the
  promoter score. · (c) SHAP values of the top 20 features: one point per gene, coloured by the
  feature's value. · (d) observed against out-of-fold predicted log2 abundance; genes more than
  2 SD off are marked as unexplained.
- **How to read:** A large ΔR² in (a) means that block carries information the others lack. (b)
  and (c) show the direction, for example higher codon adaptation going with higher abundance.
  In (d), points near the diagonal are well predicted; the title gives the cross-validated R².
- **Caveats:** When the model does not predict abundance (R² at or below 0, or its interval
  includes 0), only (a) and (d) are drawn, under a banner. Feature blocks are correlated, and the
  intervals are uncalibrated because the cross-validation repeats share data.

### V16 compare (compare report only)
- **Question:** How do search strategies differ?
- **Panels:** (a) genes detected per exact set of runs: each bar counts the genes detected by
  exactly the runs marked by the dots below it; each run's total is at the left; the 30 largest
  sets are drawn and the rest summed as "other". · (b) pairwise concordance: below the diagonal,
  Spearman ρ of the shared measured log2 abundances with its 95 % CI; above it, the shared values
  of each pair against the identity line.
- **How to read:** A tall bar under every dot means the runs agree on those genes; a bar under a
  single dot counts genes only that run finds. In (b), ρ near 1 means the runs rank abundance
  alike; high ρ with points off the identity line is a difference of scale only.
- **Caveats:** Panel (b) uses only genes measured in both runs. Runs with different versions may
  not be comparable; the Runs table says when.

### V17 qc
- **Question:** Are the inputs sane?
- **Panels:** (a) preflight checks: counts by level (ok, info, warn, error). · (b) contig
  lengths. · (c) distinct peptide lengths (aa), placed and unmapped. · (d) log2 abundance,
  measured and imputed. · (e) missing values (%) per gene feature over complete genes, with the
  reason when a feature was not computed.
- **How to read:** An error in (a) needs action and a warning a look; the list under the panels
  names them. Unmapped peptides much longer or shorter than the placed ones point to a wrong
  digest or database. A large imputed share in (d), or a feature near 100 % missing in (e),
  weakens every later plot that uses it.

### V18 stat_diagnostics
- **Question:** Are the tests calibrated on this dataset?
- **Panels:** (a) detectability reliability: the observed detection rate against the mean
  predicted probability per bin, with the slope, intercept and Brier score. · (b) p-value
  histograms, one per test family, with n; families with fewer than 10 tests are listed, not
  drawn. · (c) null distributions for regions, islands and the determinants model, with the
  observed statistic marked and its p.
- **How to read:** In (a), points on the diagonal (slope near 1, intercept near 0) mean the
  detectability model predicts well, so V07 can be trusted. In (b), a calibrated test is flat away
  from zero; a spike near zero is signal, a pile-up near 1 a conservative test. In (c), an
  observed value far in the tail is unlikely under the null.
- **Caveats:** Gene-shuffle nulls ignore operons and are too narrow; judge by the
  transcription-unit block null.

## 5. Examples (input → expected output)

| # | Input | Expected output |
|---|---|---|
| 1 | report of route A; slot V02, one instance `V02_abundance_landscape` with SVG, PNG, PDF, data | one `dl.guide#guide-V02[data-plot=V02]` inside its `figcaption`; `dt` texts exactly `Question, Panels, How to read, Caveats, Data`; Question `dd` = "Is abundance clustered along the chromosome?"; Panels `ul` has 1 `li`; Data `dd` has 4 links, the first with text `plots/data/V02_abundance_landscape.tsv` |
| 2 | V17 guide (no caveats) | `dt` texts `Question, Panels, How to read, Data`; Panels has 5 `li` starting "(a) " … "(e) " |
| 3 | V03 slot with 3 locus instances `bgc_1`, `bgc_2`, `island_1` | first figcaption has the full guide `#guide-V03`; second and third have `dl.guide.guide-short` with `dt` `Question, How to read, Data`, no `id`, a link `href="#guide-V03"`, and their own data TSV link; exactly one `#guide-V03` on the page |
| 4 | V13 with M12-05a: the index instance plus 8 candidate views (route E: 8 candidates, none called) | full guide on the first drawn V13 figure; the index instance has Question, the index sentence as How to read, Data; candidate views have short guides |
| 5 | `build_report(..., redact=True)` vs. default on the same run | the concatenated outer HTML of every `dl.guide` is byte-identical |
| 6 | the same run in `linked` mode (V13 locus views linked as `<img>`) | each linked figure still has its `figcaption` with the guide (full or short as in 3) |
| 7 | compare report of 2 runs | the V16 `figcaption` has `data-stem="V16_compare"` and the full V16 guide with 2 Panels items "(a) …", "(b) …" |
| 8 | V12 present only as an empty-state notice (routes A, D, E) | `.notice` unchanged; no `dl.guide[data-plot=V12]` |
| 9 | an instance with images but no data TSV | Data `dd` = "no data table for this figure; figure files: SVG PNG PDF" (links on the labels) |
| 10 | figcaption text, whitespace collapsed, of V07 | starts "Question Which regions and categories are detected less than detectability predicts? Panels One row per group" |

## 6. Edge cases

- **Plot with twin (Plotly) plus static SVG in `<details>`:** one figcaption, one guide, after
  both (unchanged position).
- **Plot skipped or file absent (`p.note`):** no guide.
- **Every instance of a plot empty:** notices only, no guide (decision 6). If the first instance
  is empty and the second drawn, the second carries the full guide (the "first drawn instance").
- **Index table plots (V03 index table, `INDEX_TABLES`):** the `<table><caption>` is unchanged;
  only figures get guides. A V03 slot whose only instance is the index table has no guide.
- **Replicon count, strand, circular/linear:** the guide text is fixed per plot and does not vary
  with the data (one panel per replicon is described in the Panels text itself).
- **Unknown plot ID in a slot (a future plot without a `GUIDES` entry):** `build_report` raises
  `KeyError` naming the ID (tests keep `GUIDES` keys = `viz.base.registered()`).
- **Compare report with no drawable V16** (`sim_empty`): the note is unchanged, no guide.
- **Wording:** every field passes M12-04 AC 2.1 (no internal references or raw flag tokens);
  "(a)" letters match the letters printed in the figure.

## 7. Acceptance criteria

- **AC 1** `GUIDES` keys equal the registered plot IDs; every `Guide` has non-empty `question`,
  `read` and at least one `panels` item; `caveats` is `None` or non-empty; `question` equals the
  §10.2 Question (the existing `QUESTIONS` helper of the M12-04 tests).
- **AC 2** In every fixture report, every `figure` whose plot has a drawn instance has exactly one
  `figcaption > dl.guide`; its `dt` texts are, in order, `Question, Panels, How to read,
  [Caveats,] Data` (full) or `Question, How to read, Data` (short); every `dd` has non-empty text.
- **AC 3** Panels `dd` holds one `ul` whose `li` count equals `len(guide.panels)`; for a guide with
  more than one item, the k-th `li` text starts with `(a) `, `(b) `, … in order.
- **AC 4** Each plot ID has exactly one element `#guide-<pid>` per page, on its first drawn
  figure; every short guide links to it.
- **AC 5** The Data `dd` of each instance links to `../plots/data/<stem>.tsv` when the instance
  has data (link text `plots/data/<stem>.tsv`) and to each image file; the link checker of §10.3
  still finds no broken link.
- **AC 6** `redact=True` and `linked` reports have byte-identical `dl.guide` markup to the default
  report of the same run.
- **AC 7** The compare report's V16 figure has the full V16 guide and `data-stem="V16_compare"`.
- **AC 8** M12-04 AC 2.1 (wording) holds for the `dt` and `dd` text of every guide; AC 2.2 is
  restated as: the figcaption text starts with "Question " followed by the §10.2 Question.
- **AC 9** Rendering is pure text: no image file, no `plots/data` TSV and no `OUT/README.txt`
  changes (pytest-mpl baselines and byte-identity snapshots pass untouched).
- **AC 10** (after M12-05a) the V13 Panels (b) item and the index sentence describe the M12-05a
  ordering (every candidate, smallest q first); checked by a test that asserts "every candidate"
  and "smallest q first" appear in the V13 guide.
- **AC 11** The verbatim V02 and V07 guides of § 4 are asserted word for word (the two worked
  examples); the other guides are asserted only by AC 1–3.

## 8. Snapshot-test impact

- `tests/snapshot/report_structure.json`: **no change.** The extractor records
  `["figure", data-stem]` from each `figcaption`; report figcaptions keep `data-stem`, and the
  compare figcaption gains `data-stem="V16_compare"`, which equals the first word it records
  today. `dl`, `dt`, `dd` are not recorded. If the extractor is later extended to `dt`, the
  baseline is regenerated with Khoa's approval.
- pytest-mpl baselines and `test_snapshot_byte_identity.py`: no change (AC 9).
- Tests to update (test-author, PR 1): `tests/unit/report/test_m12_04_pr2_captions.py`,
  `tests/unit/report/test_m12_04_pr2_review_r1.py` (they import `CAPTIONS` and compare the whole
  figcaption text) and `tests/helpers/m12_04_pr2.py`; the toy page in
  `tests/snapshot/test_report_structure.py::test_rich_report_inlines…` uses a hand-written
  figcaption and may stay.
- The weekly e2e and integration tests that grep report text: the lead greps
  `tests/integration` for "Are the inputs sane? Panels show" before PR 2 merges.

## 9. Work split

| PR | Title | Role | Lines | Content |
|---|---|---|---|---|
| 1 | `test(report): M12-05c structured plot guides` | test-author | ≈ 260 | AC 1–9, 11 on the M12-04 sim fixtures plus a 3-locus V03 fixture and a compare fixture; update the 3 M12-04 caption test files; AC 10 as `xfail(strict=True)` until M12-05a merges |
| 2 | `feat(report): structured plot guides (M12-05c)` | implementer | ≈ 290 | `Guide`, `GUIDES` (§ 4 text, ≈ 190 lines), `_plot_slot` and the compare slot, both templates, `report.css`; `changelog.d/<PR>.feature.md` ("Each figure in the report now has a labelled guide: Question, Panels, How to read, Caveats, Data.") |

If PR 2 exceeds ~400 lines, split it as 2a (mechanism, templates, V01–V08) and 2b (V10–V18),
with PR 1's per-plot tests parametrised so 2a's missing entries are skipped by ID until 2b.
Order: PR 1 → PR 2; M12-05a's author flips the AC 10 xfail and edits the V13 entry, or, if
M12-05a lands first, PR 2 ships the M12-05a text directly. Snapshot baselines do not change, so
no baseline hand-off to Khoa is needed; the lead shows Khoa a rendered route A report before
PR 2 merges (review outputs first).

## 10. Questions for the lead / Khoa (defaults are used if unanswered)

- **Q1 labels.** "Question / Panels / How to read / Caveats / Data", or Khoa's "Instruction:"
  as an outer heading? *Default: the five labels, no "Instruction" heading (the `dl` is
  already the instruction block).*
- **Q2 open or collapsed.** Show the guide always, or inside `<details>` ("How to read this
  plot")? *Default: always shown, so it prints and is searchable.*
- **Q3 outside the report.** Also write the guides to `plots/README.txt` or to
  docs/user-guide/plots-and-report.md (generated from `GUIDES`)? *Default: report only; a
  backlog item for the docs page.*
- **Q4 empty states.** Should a notice-only plot show its Question too? *Default: no change.*
- **Q5 V13 (b) under M12-05a.** This spec assumes M12-05a draws every candidate as a locus view
  sorted by `class_fdr` (q), unestimated last. If M12-05a lists them in a table instead, the
  Panels (b) item becomes "a table of every candidate …". *Default: the M12-05a author adapts the
  V13 entry.*
- **Finding (not in scope):** route E's V13 panel (a) prints raw class tokens
  (`nterm_extension`, `cterm_extension`, `cds_outframe`, `sorf`) on the axis; M12-04 AC 2.1 checks
  report prose, not SVG text. Suggest M12-05a renders them through `words.plain`.

## Lead decisions (2026-11-09)
- Khoa: Q1 and Q2 "Five labels, visible": Question, Panels, How to read, Caveats, Data, always shown, with no extra heading and no `<details>`. Q3 to Q5: the defaults above. M12-05a merges first; this spec's V13 guide follows M12-05a (every candidate, class q first, then best peptide q).
