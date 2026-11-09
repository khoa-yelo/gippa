# M12-05a: the novel ORF plot shows every candidate, lowest q first

Milestone M12 (1.1.0) · depends on: M9-12 (V13), M12-04 (PR 4 empty state, PR 5 locus palette,
PR 11c numeric `class_fdr`, PR 13 linked locus views) · estimate ≈ 160 code lines (implementation)
and ≈ 350 test lines, in 2 PRs (§ 7); regenerated baselines are not counted · roles: spec-writer
(this file), test-author (`tests/unit/viz/test_m12_05a_*.py`,
`tests/unit/report/test_m12_05a_*.py`, updates to the superseded V13 tests listed in § 7),
implementer (`src/gippa/viz/novel_orfs.py`, `src/gippa/viz/locus.py`,
`src/gippa/report/content.py`, `src/gippa/report/build.py`), lead (baseline hand-off to Khoa) ·
not statistical (no estimate changes; ordering and drawing only).

**Plot-only change.** No table, schema, writer or step changes. `tables/novel_orfs.tsv` is read
as it is today. The only files whose content changes are V13's own plot files
(`plots/V13_*`, `plots/data/V13_*`) and the V13 slot of `report/report.html`.

## Source: Khoa's review of the 1.1.0rc1 outputs (2026-11-09)

> "Should visualize the novel orfs regardless of their q value (sorted by low)."

Read as three requirements:

1. **Every candidate is drawn**, whether it is called or not.
2. **Order is by q, lowest first.**
3. **Called and uncalled candidates must look different**, because the plot no longer shows only
   calls.

## Evidence run (read-only)

Route E of the example bundle, 1.1.0rc1:
`local/example/routes/1945698/gippa-example/runs/E/`.

- `tables/novel_orfs.tsv` has 8 candidates and 0 called. In the per-class picked FDR, `class_fdr`
  is 0.4 for the 5 antisense rows and 1 for the `nterm_extension`, `intergenic` and
  `cds_outframe` rows. These are the class q floors (D + 1) / T. Every row has `n_peptides = 1` and
  `rank` NA.
- `plots/` has `V13_novel_orfs.{svg,png,pdf}` (panel (a) only) and
  `plots/data/V13_novel_orfs.tsv`. There is no locus view and no index, so **the report shows no
  candidate at all**, only four bar counts.
- `report/report.html` › Novel ORFs gives the caption "… panel (b) draws the top called candidates
  in their genomic neighbourhood." Panel (b) is absent.

## Current behaviour (confirmed in code at e6f4aa3c)

- `src/gippa/viz/novel_orfs.py::novel_orfs`:
  - Panel (a) (`_panel_a`) counts **every** candidate per class × status (`called`,
    `not_called`, `unestimated`). This spec leaves it unchanged.
  - Panel (b) keeps only `called = true` rows (`called.fillna(False)`). It sorts them by `rank`,
    with NA rank last, then by `orf_id`, and keeps the first `TOP_N = 20`.
  - When no row is called, `novel_orfs` returns panel (a) only (`if top.empty: return out`).
    This is why route E has no views.
  - V13 does not filter on q; it filters on `called`. The q (`class_fdr`) enters only through the
    call.
- `_index` draws a text-only image, "`<rank>. <orf_id>`" per row, titled "V13 novel ORF loci: top
  `<n>` by rank". Its data TSV `V13_novel_orfs_locus.tsv` has `rank`, `orf_id` and `rendered`
  (always true).
- `src/gippa/viz/locus.py::_draw_orf` draws one arrow style: face `NOVEL_ORF` (`#CC79A7`) at
  α 0.35, a solid edge, lw 1.0. `_orf_key` gives one legend entry, `ORF_KEY` = "novel ORF; tick:
  OSTIR start; purple pileup: its peptides".
- `src/gippa/report/content.py`: `CAPTIONS["V13"]` (one caption for panel (a) and every locus
  view) and `INDEX_CAPTIONS["V13"]` = "This index lists the top called candidates drawn as locus
  views, by rank."
- **Report order.** `src/gippa/io/plot_files.py::read_plot_files` returns instances
  `sorted(stems)`, so the report shows V13 locus views in **alphabetical `orf_id` order**, not in
  display order. Route E would show `…_130563_m` before `…_20324_m` whatever order V13 drew them
  in. (Backlog, V03 order review Q4: "show the V03 and V13 locus plots in the report in index-rank
  order"; this spec does it for V13 only.)

## Spec sections touched

- GIPPA_SPEC.md §10.2 V13, "(b) locus views of the top 20 (V3 renderer)". **Unchanged in
  wording.** "Top 20" now means the 20 candidates first in the § 2 order, not the 20 best-ranked
  calls.
- It supersedes M9-12 "Panel (b)" ("The 20 `called = true` rows with the smallest `rank`") and its
  degenerate-case row "no called ORF but unestimated candidates → panel (a) only". It also
  supersedes M7-09 Q5's reading "'top 20' (V13) means top 20 calls". The `rank` column itself
  (§7.6.5, called only) does not change.
- §7.6.5 `rank_key` (G13 + ADR-0012 d4) is reused as a tie-break, not redefined.
- No ADR is needed: §10.2 is not contradicted (lead to confirm, Q6).

## Not in scope (parallel specs; reference only)

- **M12-05b** adds host and flanking gene columns to `novel_orfs.tsv`. This spec reads no new
  column and does not show host genes in titles or the index. If 05b merges first, nothing here
  changes.
- **M12-05c** restructures every plot's caption into labelled parts. This spec changes only the
  **content** of the two V13 caption strings (§ 5). If 05c merges first, the same sentences go
  into its parts for V13, with the panel (b) sentence in its "what is drawn" part and the style
  sentence in its "how to read it" part. Whichever PR merges second resolves the V13 entry.
  Neither spec changes the other's mechanism.
- Panel (a), its title, its data TSV and its baseline (`sim_plots-V13_novel_orfs.png`) are
  unchanged.
- V03's report order is out of scope.

## What already exists (reuse, do not duplicate)

- `gippa.novel.evidence.rank_key(row)`: the §7.6.5 key with NA last and `orf_id` as the final
  tie-break. It is pure and takes a `Mapping` with `n_peptides`, `best_peptide_q`,
  `start_codon_best`, `tir_ostir`, `promoter_score`, `terminator_conf`, `length_aa_best` and
  `orf_id`, which are all `novel_orfs.tsv` columns. `gippa.viz.variants` already imports from
  `gippa.variants`, so a viz → novel import has a precedent.
- `novel_orfs._fdr` (finite float, or `None` for NA / `unestimated`), `_status`, `_span`,
  `_tick`, `_view`, `TOP_N = 20`.
- `locus.render_locus(..., orf=OrfMark)`, `_draw_orf`, `_orf_key`, `NOVEL_ORF`.
- The cap wording of V12: `f" (top {SHOWN_SITES} of {len(alt)}; all in tables/variants.tsv)"`
  (`gippa/viz/variants.py::title`).
- `style.fmt_num`, `style.fmt_q` (M12-04 PR 0).
- `report/build.py::LINKED_PREFIX` (the linked-view mechanism; unchanged).

---

## 1. What is drawn

- **Every row of `novel_orfs.tsv` is a candidate for a locus view**, whatever its `called`,
  `class_fdr`, `rank` or `flags`. `called` NA counts as not called (as today).
- V13 returns these instances:
  - panel (a), unchanged;
  - one locus view (`V13_novel_orfs_locus_<orf_id>`) for each of the first `N = min(TOP_N, M)`
    candidates in the § 2 order, where `M = len(novel_orfs.tsv)`;
  - the index instance (`V13_novel_orfs_locus`) whenever `M ≥ 1`.
- With no row called, the views and the index are still drawn. **The M9-12 "panel (a) only" rule
  is removed.**
- A table with no rows still gives the M12-04 PR 4 empty state, unchanged: one instance, no index,
  no locus files.

## 2. Order (the display order)

Candidates are sorted, stably, by this key:

1. **`class_fdr` ascending.** This is the candidate's q, the value its call is tested on (§7.6.4:
   "a call needs class q ≤ 0.01"). It is read with `_fdr`, so NA, the literal `unestimated` and
   non-finite values are "no q" and sort **after every finite q**. Q1 asks whether Khoa meant
   this q.
2. **`best_peptide_q` ascending, NA last** among equal class q (Khoa, Q1: "Class q, then peptide
   q").
3. **Called before not called** among equal class q and best peptide q (including among "no q").
4. **`rank` ascending, NA last.** Among called rows this is the pipeline's own §7.6.5 order.
   Uncalled rows all have `rank` NA, so they tie here.
5. **`rank_key(row)`** (§7.6.5): peptides descending, `best_peptide_q` ascending, start codon,
   TIR descending, promoter descending, terminator descending, length descending, and finally
   `orf_id` ascending. NA comes last in every component.

The key is total, because `orf_id` is unique. The same table therefore gives the same order on
every render, whatever its row order.

Why these tie-breaks:

- `class_fdr` is pooled per class, so ties are the common case. All 5 antisense rows of route E
  share 0.4, and 13 pseudogene rows of the sim fixture share 0.00332862. Calling is not decided by
  q alone (for example, ≥ 2 peptides for `antisense` and `cds_outframe`). Putting calls first
  within a q therefore puts the stronger evidence first.
- `rank` and then `rank_key` reuse the pipeline's documented evidence order, so V13 introduces no
  new ranking rule. For called rows of equal q, `rank` and `rank_key` agree, because `rank` is
  computed with `rank_key` over the called rows. `rank` is used first so that hand-built tables
  with an explicit `rank` keep the existing review test (`test_v13_review.py`, called row with NA
  rank last).

## 3. How called and uncalled candidates differ

| | called (`called = true`) | not called (`false` or NA) |
|---|---|---|
| ORF arrow face | `NOVEL_ORF` at α 0.35 (unchanged) | none (fully transparent) |
| ORF arrow edge | `NOVEL_ORF`, solid, lw 1.0 (unchanged) | `NOVEL_ORF`, **dashed** `(0, (3, 2))`, lw 1.0 |
| peptide pileup overlay, OSTIR tick | `NOVEL_ORF` (unchanged) | `NOVEL_ORF` (unchanged) |
| key entry (exactly one per view) | patch as the arrow, text `ORF_KEY_CALLED` | patch as the arrow, text `ORF_KEY_UNCALLED` |
| title | "`<orf_id>` (`<class>`, called): …" | "`<orf_id>` (`<class>`, not called): …" |

- `ORF_KEY_CALLED` = "called novel ORF; tick: OSTIR start; purple pileup: its peptides".
- `ORF_KEY_UNCALLED` = "novel ORF candidate, not called (dashed); tick: OSTIR start; purple
  pileup: its peptides".

**Why this style.**

- **One hue.** The arrow, its tick and its pileup overlay stay one colour. That colour is what
  links the ORF to its peptides, and the reader still finds the ORF track by colour.
- **Fill and dash encode the call**, which shows twice and survives greyscale printing and colour
  vision deficiency. A hollow, dashed outline reads as "proposed, not established". That matches
  the M12-04 P-01 idiom: undetected genes are white with a grey outline.
- **No new colour.** Panel (a)'s status colours (`#0072B2` called, `#BBBBBB` not called,
  `#E69F00` unestimated) are not reused in the locus view. There, blue and grey already mean
  pileup and undetected genes, and the cividis gene scale already uses blue/yellow, so they would
  be ambiguous.
- **One key entry per view, describing that view's ORF.** This keeps the key at one entry, as the
  existing review test requires. The title says the status in words as well, so the two styles do
  not have to be compared side by side.
- **Two styles, not three.** An uncalled row with no q (panel (a)'s `unestimated`) is drawn as not
  called. Its title already says "class FDR unestimated" (unchanged `_fmt`).

## 4. Cap for large runs

- `N = TOP_N = 20`. This is the §10.2 V13 number and today's `TOP_N`. The siblings are V12 at 15
  (`SHOWN_SITES`) and V03 at 50 (`report.max_locus_plots`). V13 keeps its own §10.2 figure, so
  the report never draws more V13 views than today's maximum and its size budget (M12-04 PR 13) is
  unchanged.
- The cap applies to the § 2 order, so the views are "the 20 lowest-q candidates".
- **The cap note** goes in the **index image title**, worded as in V12:
  - `M > N`: "Novel ORF loci, lowest class FDR first (top 20 of `<M>`; all in
    tables/novel_orfs.tsv)";
  - `1 < M ≤ N`: "Novel ORF loci, lowest class FDR first (`<M>` candidates)";
  - `M = 1`: "Novel ORF loci, lowest class FDR first (1 candidate)".

  The report caption is static per plot, so the count goes in the figure title (as in V12) and
  not in the caption.
- The index image lists the first N rows of the order. Its **data TSV lists all M** (§ 6).

## 5. Locus views, index image and captions

- **Locus views follow the § 2 order and the § 4 cap.** V13 returns them in that order.
- **The report shows them in that order** (§ 6.3). Span, window, origin unwrapping, the
  `start_longest`-NA fallback and the minus-strand +1 rule do not change, and apply to uncalled
  rows exactly as to called ones.
- **Locus title.** The only change is the status word in the parentheses:
  "`<orf_id>` (`<class>`, called|not called): `<n>` peptides, best q `<fmt_q>`, class FDR
  `<fmt_num>`|unestimated".
- **Index image** (`V13_novel_orfs_locus.{svg,png,pdf}`):
  - double column width (`style.new_figure("double", …)`), because route E's `orf_id`s are 54
    characters long;
  - one line per shown row, in four text columns: "`<order>`. `<orf_id>`", `<class>`,
    "called"|"not called", and "class FDR `<fmt_num>`"|"class FDR unestimated";
  - the title from § 4.
- **Captions** (`report/content.py`, content only; see M12-05c above):
  - `CAPTIONS["V13"]` = "Which new genes or upstream extensions are supported? Panel (a) counts
    the candidates of each class by call status, and panel (b) draws up to 20 candidates, called
    or not, in their genomic neighbourhood, lowest class FDR first. A filled arrow is a called
    ORF, a dashed outline a candidate that was not called."
  - `INDEX_CAPTIONS["V13"]` = "This index lists the candidates drawn as locus views, lowest class
    FDR first, with their call status."
  - Both strings pass the M12-04 rule 2 checks: no ADR, §, raw flag token or path other than
    `tables/…`. `tables/novel_orfs.tsv` appears only in the figure title, as in V12.

## 6. Data TSVs and report order

### 6.1 `plots/data/V13_novel_orfs_locus.tsv` (the index data; changed)

Columns, in this order:

`order`, `orf_id`, `class`, `called`, `class_fdr`, `best_peptide_q`, `n_peptides`, `rank`,
`rendered`

- There is one row per candidate (**all M**, not just N), in the § 2 order. `order` runs from 1
  to M.
- `called` is a bool, with NA → false. `class_fdr` is a float, with NA when `_fdr` is `None`.
- `rank` is an int or NA, copied from the table. `rendered` is true for `order ≤ N`.
- It is written with the house TSV writer (`%.6g`, NA as `NA`).
- This replaces the M9-12 columns `rank, orf_id, rendered`. No schema or reference page lists this
  file, so the change is internal to V13 (no changelog `removal`). The PR body says so.

### 6.2 Other V13 data TSVs

- `V13_novel_orfs.tsv` (panel (a)) is unchanged.
- `V13_novel_orfs_locus_<orf_id>.tsv` keeps the M9-10 layout. Its `novel_orf` row gets the call in
  its `category` cell (`called` or `not_called`; today the cell is NA). This is the one
  machine-readable trace of the style. It needs no new column.

### 6.3 Report order (`src/gippa/report/`)

- In the V13 slot, the instances are shown as:
  1. panel (a) (`V13_novel_orfs`);
  2. the index (`V13_novel_orfs_locus`);
  3. the locus views, in the `order` of `plots/data/V13_novel_orfs_locus.tsv` among its
     `rendered = true` rows.
- The ordering is a pure function in `report/content.py`, for example
  `display_order(stems, index_rows) -> list[str]`.
- The index TSV is read with `gippa.io.plot_files.read_tsv`, where `build.py` already holds the
  plots directory. File I/O stays in `io/`.
- A locus stem that is missing from the index (a stale file) goes after the listed ones, in
  alphabetical order. A missing or unreadable index TSV leaves today's alphabetical order.
- Other plots are untouched.
- Linking (`LINKED_PREFIX`, PR 13) is unchanged and applies to the reordered views.

---

## 7. Examples (input → expected output)

`…` abbreviates `NOVEL_nissle1917_NZ_CP058217_1_1400000_1600000_`.

| # | Input | Expected |
|---|---|---|
| E1 | Route E `novel_orfs.tsv`: 8 rows, 0 called. Antisense q 0.4 ×5, other classes q 1 ×3; best q in the evidence run | 8 locus views, all dashed and unfilled, in this order: `…130563_m` (0.4, best q 0.000875), `…58745_p` (0.4, 0.00247), `…195190_m` (0.4, 0.00339), `…43000_m` (0.4, 0.00509), `…20948_p` (0.4, 0.00721), `…20324_m` (1, 0.000875), `…25949_p` (1, 0.00721), `…149890_p` (1, 0.00779). Index title "Novel ORF loci, lowest class FDR first (8 candidates)". The index TSV has 8 rows, all `rendered = true`. Panel (a) is byte-identical to today's. The report shows the 8 views in this order |
| E2 | Sim fixture `sim_plots` (57 rows: 40 called, 17 not) | Views 1–20 are all called: the 13 pseudogenes (q 0.00332862) in `rank` order 1, 2, 3, 4, 7, 9, 13, 17, 18, 20, 23, 24, 32, then 7 antisense (q 0.00383386), ranks 5, 6, 8, 10, 11, 12, 14. Index title "… (top 20 of 57; all in tables/novel_orfs.tsv)". The index TSV has 57 rows. Row 39 is `NOVEL_chr1_48669_p` (antisense, not called, q 0.00383386, after the 26 called antisense rows of that q). Rows 40–57 are uncalled, and row 57 is `NOVEL_chr1_74837_p`. Views 1 and 2 are `NOVEL_chr1_45143_m` and `NOVEL_chr1_44225_m`, the two kept baselines, as today |
| E3 | All called, mixed q (E2's called rows) | The order is by q first, so `NOVEL_chr1_269094_p` (rank 32) is view 13. `NOVEL_chr1_185480_p` (rank 15), `NOVEL_chr1_247263_p` (16) and `NOVEL_chr1_130084_m` (19) are not drawn (`rendered = false`). The report-structure stems change accordingly |
| E4 | Tie on q, called vs not: 2 antisense rows, both q 0.0065; A `called = true`, `n_peptides = 2`, rank 1; B `called = false`, `n_peptides = 5` | A then B (both have best peptide q 0.001, so step 3 outranks B's peptide count) |
| E5 | Tie on q and call, both uncalled: q 0.4; C best q 0.0025, D best q 0.0009, both 1 peptide | D then C (`rank_key`: best q ascending) |
| E6 | Full tie: 2 uncalled rows equal in every key column, IDs `NOVEL_chr1_900_p` and `NOVEL_chr1_1000_p` | `NOVEL_chr1_1000_p` then `NOVEL_chr1_900_p` (string order of `orf_id`). The same output follows when the input rows are swapped |
| E7 | NaN q: E called, `class_fdr` NA, rank 1; F not called, q 0.9; G not called, `class_fdr` = literal `unestimated` | F, E, G. Finite q comes first; among "no q", called comes first. E's title says "class FDR unestimated" and its arrow is filled. G's arrow is dashed. In the TSV, `class_fdr` is NA for E and G |
| E8 | 25 rows, none called, q 0.01 + 0.01·i | 20 dashed views (i = 0…19). Index title "… (top 20 of 25; all in tables/novel_orfs.tsv)". The index TSV has 25 rows; rows 21–25 have `rendered = false` |
| E9 | 1 row, not called, minus strand, crossing the origin of a circular 10 kb contig (`start_longest = 150`, `end = 9 900`) | Panel (a), 1 dashed view and the index "… (1 candidate)". The span and the unwrapped window equal those of the same row with `called = true` (today's example 3 rule); only the arrow face, the dash, the key text, the title word and the `category` cell differ |
| E10 | Header-only `novel_orfs.tsv` | The M12-04 PR 4 empty state, unchanged: one instance, no `V13_novel_orfs_locus*` files |

## 8. Edge cases

| Case | Behaviour |
|---|---|
| 0 candidates | E10: empty state, no index, no views |
| 1 candidate | E9: 1 view, the singular index title |
| all called | E2/E3: order by q, then `rank`. Every view is filled. The panel (a) title is unchanged |
| none called | E1/E8: every view is dashed. The M9-12 "panel (a) only" rule is gone |
| ties | E4–E6: best peptide q, then called first, then `rank`, then `rank_key`, then `orf_id`. The output does not depend on input row order |
| NaN / `unestimated` q | E7: sorts after every finite q. NaN `best_peptide_q` and other NA `rank_key` components sort last within their step |
| `called` NA | treated as not called, so dashed (as today's panel (b) filter) |
| `flags` has `class_fdr_unestimated` but `class_fdr` is finite | sorts by the finite value (the sort reads `class_fdr` only). After M12-04 PR 11c this cannot occur, since that flag comes with NA |
| `start_longest` NA (`no_start_codon`, route E `…43000_m`, `…149890_p`) | unchanged first-peptide fallback; drawn as any other uncalled row |
| minus strand / origin crossing / linear contig | unchanged span and window rules (M9-12 decisions 6–7). Style depends only on `called` |
| more than 20 candidates | top 20 drawn; the index title has the V12-form note; the index TSV lists all |
| very many candidates (e.g. 5 000) | still 20 views. The index TSV has 5 000 rows, and the index image only 20 lines |
| a stale `V13_novel_orfs_locus_<id>` file not in the index | shown after the listed views, in alphabetical order (§ 6.3) |
| report over the PR 13 size threshold | the views are linked as today, in the § 6.3 order |

## 9. Acceptance criteria

Fixtures are hand-built `novel_orfs.tsv` tables, as in `tests/unit/viz/test_v13_novel_orfs.py`,
plus the route E table copied into `tests/data/` as a small fixture (8 rows, no peptide sequences
needed beyond the table). The test-author may instead rebuild it from the values in E1.

- **AC 1 (every candidate).** For a table with M rows, V13 returns `1 + min(20, M) + 1` instances
  when M ≥ 1: panel (a), the views and the index. With E1 that is 10 instances, with E8 22, and
  with E10 one (the empty state).
- **AC 2 (order).**
  - The `suffix` sequence of the views equals `["locus_" + id for id in expected]` for E1, E2
    (first 20), E4, E5, E6 and E7.
  - A permutation of the input rows gives the identical sequence. Make this a hypothesis property
    over random tables of up to 30 rows with random q (including NA), random `called` and ties.
  - The sequence is non-decreasing in (`_fdr` with None as +∞, not called).
- **AC 3 (style).** In each view, exactly one `Polygon` in the `genes` axes has edge `#CC79A7`.
  - Called: face alpha 0.35 and a solid line style.
  - Not called: face alpha 0 and a dashed line style.
  - Exactly one legend text contains "novel ORF": `ORF_KEY_CALLED` or `ORF_KEY_UNCALLED`, as
    appropriate. Its handle has the arrow's face, edge and line style.
- **AC 4 (title).** The locus title of E9 is "`<orf_id>` (antisense, not called): 1 peptides …"
  (the `n` wording is as today). For a called row it is "(`<class>`, called)". The rest of the
  string is as `test_locus_title` today.
- **AC 5 (cap note).**
  - Index title for E2: "Novel ORF loci, lowest class FDR first (top 20 of 57; all in
    tables/novel_orfs.tsv)".
  - For E1: "Novel ORF loci, lowest class FDR first (8 candidates)".
  - For E9: "Novel ORF loci, lowest class FDR first (1 candidate)".
  - The index image has `min(20, M)` row lines, each with its order number, `orf_id`, class, call
    word and class FDR text.
- **AC 6 (index TSV).** `plots/data/V13_novel_orfs_locus.tsv` has exactly the § 6.1 header.
  - It has M rows and `order` = 1…M in the § 2 sequence.
  - `rendered` is true exactly for `order ≤ 20`.
  - For E7, `class_fdr` is `NA` for E and G.
  - For E8, there are 25 rows and 5 `false` values.
- **AC 7 (locus TSV).** The `novel_orf` row of each view's data TSV has `category` = `called` or
  `not_called`. Every other row is byte-identical to the same view rendered from a table where
  that row has the other `called` value. Only the `novel_orf` row's `category` differs.
- **AC 8 (geometry unchanged).** For E9 and for today's example 3, the arrow vertices, the window
  limits and the tick labels are equal between `called = true` and `called = false`.
- **AC 9 (panel (a) unchanged).** On E1, E2 and E7, panel (a)'s data TSV and title are equal to
  those of `e6f4aa3c`. The baselines `sim_plots-V13_novel_orfs.png` and `sim_empty-V13_novel_orfs.png`
  still pass unchanged.
- **AC 10 (report order).**
  - `content.display_order` is a pure function. Given stems in alphabetical order and E1's index
    rows, it returns panel (a), the index, then the 8 views in E1 order.
  - A stem missing from the index goes last.
  - With no index rows, it returns the input order.
  - In `build.py`, the V13 slot of a report built from an E1-like `OUT/` lists
    `data-stem`/`<figcaption data-stem>` in that order. This also holds when the views are linked
    (PR 13 fixture with padded SVGs).
- **AC 11 (captions).** `content.caption("V13", "V13_novel_orfs")` and
  `caption("V13", "V13_novel_orfs_locus_<id>")` equal the § 5 `CAPTIONS["V13"]`.
  `caption("V13", "V13_novel_orfs_locus")` equals the Question followed by the § 5
  `INDEX_CAPTIONS["V13"]`. Both pass the M12-04 rule 2 prose checks.
- **AC 12 (plot-only).** The PR diffs touch no file under `src/gippa/novel/`, `src/gippa/schemas/`,
  `src/gippa/steps/` or `src/gippa/io/` (only `read_tsv` is called, which already exists), and no
  `tables/*` golden. Tests check this with `git diff --name-only origin/main...HEAD` in review, not
  in CI.
- **AC 13 (determinism and purity).** Two renders of E1 give byte-identical files (as
  `test_files_are_byte_identical_across_renders`). There is no pyplot figure, and every text is at
  least 6 pt (existing V13 checks, extended to the new index image).
- **AC 14 (route E, validation, not CI).** On a copy of route E's `OUT/` (never in place), run
  `gippa report --out <copy> --plots V13` on the implementation branch, inside an allocation or
  sbatch.
  - It writes 8 `V13_novel_orfs_locus_*` views and the index.
  - The report's Novel ORFs section shows them in E1 order.
  - The PR body gives scratchpad PNG paths of panel (a), the index and two views (one antisense,
    one `nterm_extension`) for Khoa.

## 10. Work split

Both PRs are ≤ ~400 lines (baselines and `report_structure.json` not counted).

| PR | Title | Role | Contents | Lines ≈ | Baselines / fixtures |
|---|---|---|---|---|---|
| A | `test(viz): M12-05a novel ORF plot shows every candidate` | test-author | New `tests/unit/viz/test_m12_05a_v13_all_candidates.py` (AC 1–9, 12, 13; the E1–E10 fixtures; the AC 2 hypothesis property) and `tests/unit/report/test_m12_05a_v13_order.py` (AC 10, 11). Updates the superseded tests (below), marked xfail-strict against `main` where the old behaviour is still live, so that PR A is green on its own | 350 | none |
| B | `feat(viz): V13 draws every novel ORF candidate, lowest q first` | implementer | `viz/novel_orfs.py`: the order, cap, index image and TSV. `viz/locus.py`: the `OrfMark.called` field, the uncalled style, the two key texts and `category`. `report/content.py`: captions and `display_order`. `report/build.py`: the V13 slot order. Plus `changelog.d/<PR>.feature.md`. Removes PR A's xfail marks | 160 | regenerates `sim_plots-V13_novel_orfs_locus.png` (index), `sim_plots-V13_novel_orfs_locus_NOVEL_chr1_45143_m.png` and `…_44225_m.png` (key text and title word), and `tests/snapshot/report_structure.json` (the V13 stems of E3 and their order). Any other baseline change is a review failure. Khoa approves the baselines (scratchpad PNG paths in the hand-off, no Artifact pages) |

**Existing tests that PR A updates** (the test-author edits them, M12-04 rule 6):

- `tests/unit/viz/test_v13_novel_orfs.py`:
  - `INDEX_COLUMNS` → the § 6.1 header;
  - `test_no_called_but_unestimated_is_panel_a_only` → views are drawn (E1-like);
  - `test_example2_25_called_gives_20_views_ranks_1_to_20` → still 20 views, ordered by § 2
    (equal q gives rank order);
  - `test_uncalled_rows_never_get_a_view_even_with_smaller_rank` → inverted: uncalled rows get
    views, in q order;
  - `test_fewer_than_20_called_gives_that_many_views` → counts all candidates;
  - `test_locus_title` → the status word.
- `tests/unit/viz/test_v13_review.py`:
  - `test_every_locus_view_has_a_novel_orf_key` → accepts either key text, and checks that it
    matches the row;
  - `test_called_row_with_na_rank_sorts_after_every_ranked_row` → should stay valid unchanged
    (same q, all called, so step 3 applies). The test-author confirms it.
- `tests/unit/report/test_m12_04_pr2_captions.py` / `test_m12_04_pr2_review_r1.py`: only where
  they pin the old V13 caption text.
- `tests/snapshot/test_plots.py`: the rich-fixture assertion `called ≥ LOCUS_KEEP` stays. Views 1–2
  are unchanged IDs (E2).

**Order.** A, then B (B is stacked on A). Run targeted tests on the session allocation (the files
above and `tests/snapshot -k V13`). Run the full suite and the snapshot suite through
`scripts/slurm/test_branch.sbatch`. AC 14 runs after B, by sbatch on a copy of route E.

## 11. Questions for the lead (and Khoa)

Each question has a reversible default, which the PRs use unless the lead says otherwise.

| # | Question | Default |
|---|---|---|
| Q1 | **Which q?** The comment says "q value". A candidate carries two: `class_fdr` (the class-level q the call is tested on, §7.6.4; heavily tied, 0.4 / 1 in route E) and `best_peptide_q` (its best peptide's PSM-level q; distinct per row in route E, but NA everywhere in the de novo sim fixtures). | `class_fdr` first, with `best_peptide_q` inside the tie-break (via `rank_key`). On route E, the alternative (sort by `best_peptide_q`) would put `…20324_m` (q 1, best q 0.000875) second instead of sixth |
| Q2 | Sorting by q over ranks across classes: in the sim fixture, rank 32 (pseudogene) is drawn and ranks 15, 16, 19 (antisense, more peptides) are not (E3). Is that intended, or should called rows always come first in `rank` order, followed by uncalled rows by q? | by q (the literal reading of "sorted by low") |
| Q3 | Cap: keep §10.2's 20, or align with `report.max_locus_plots` (50, the V03 cap; would make V13 configurable)? | 20 |
| Q4 | The V13 index stays an image (M12-04 Q5 kept it). Alternatively, it could become an image-free report table like V03's (P-12), which would show all M candidates in the report rather than 20 lines. | keep the image; the TSV holds all M |
| Q5 | Reordering the report's V13 locus views (§ 6.3) is outside "plot-only" in the narrow sense, since it touches `report/build.py`. Without it, the § 2 order is invisible in the report (alphabetical stems). Same change for V03 (backlog Q4)? | V13 now; V03 stays in the backlog |
| Q6 | No ADR: §10.2's "top 20" wording is kept and only the M9-12/M7-09 readings change. Confirm, or record it in an ADR alongside M12-05b/c. | no ADR; changelog `feature` fragment only |

## Lead decisions (2026-11-09)
- Khoa: Q1 "Class q, then peptide q": sort by `class_fdr` ascending, ties by best peptide q (then the rest of the tie-break above). Q2 to Q6: the defaults above.
- Added from the M12-05c review: V13 panel (a) currently prints raw class names (`nterm_extension`, `cterm_extension`, `cds_outframe`, `sorf`) on its axis. Use the plain words the report already uses for these classes (`src/gippa/viz/words.py`, or the report's class labels), and add an acceptance criterion: no raw class identifier appears in the V13 SVG text.
- Review round 1 (lead, 2026-11-09): § 2 now lists best peptide q as step 2, matching Q1. AC 12 has one exception: `src/gippa/steps/plots.py` passes the plots directory to `build_report`. The report step builds its view from the step inputs, so `build.py` cannot find the V13 index without it, and § 6.3 wrongly assumed it could. Tests added: a fixed example where a smaller best peptide q outranks a call, and the M = 20 and M = 21 index titles.

