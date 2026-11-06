# M12-04: outputs revision for 1.1 (plots, published tables, HTML report)

Milestone M12 (1.1.0) · depends on: M9 (plots and report, M9-01 style and framework), M4/M5/M6/M7
table writers, M11-14h (`scripts/release/schema_diff.py`), M11-20 (initiator M, #490), M12-01
(example bundle) · estimate ≈ 4,300 code lines over 17 PRs (table in § 2), each ≤ ~400 code lines;
regenerated goldens, snapshot baselines and generated reference pages are not counted · roles:
spec-writer (this file), lead (**ADR-0075**, routing, baseline hand-off to Khoa), test-author
(`tests/unit/test_m12_04_*.py`, `tests/integration/test_m12_04_*.py`, regenerated goldens),
implementer (`src/gippa/viz/`, `src/gippa/report/`, `src/gippa/schemas/`, the table writers,
`workflow/modules/local/gippa_step/main.nf`, `src/gippa/runner/`, `scripts/release/`,
`scripts/docs/`), validation (the K re-check run on main, PR 10 and PR 13), statistician (sign-off
on PR 10's confidence rule) · source: the outputs review
`local/handoff/outputs_review_1_1.md` (39 items, 13 proposed PRs) and Khoa's decisions of
2026-11-04 (in session).

## Source: Khoa's decisions, 2026-11-04 (in session; binding)

1. **Scope "Everything".** All 39 items go into 1.1.
2. **Columns "1.1, clean break".**
   - Reorder `genes.tsv` into the groups identity → detection → abundance → function → context →
     features → flags (T-03).
   - Rename (T-04): `ostir_tir`→`tir_ostir`, `mw`→`mw_da`, `peptide_len`→`peptide_len_aa`,
     `promoter_between`→`promoter_score_between`. Align the pairs `class`/`flags` with the
     operon_flags `uniop_class`/`flag`, with one consistent name per pair (this spec proposes the names).
   - `class_error_rate` and `class_fdr` become numeric: NA plus a flag when unestimated (T-05).
     T-06 is in scope.
   - **No aliases.** This is an exception to the spec's major-version rule (§9.2 "breaking column
     changes only in major versions"; ADR-0060 §3). The lead writes the ADR (**ADR-0075**).
3. **OUT layout "Move + README" (O-01).** The top-level intermediates move to `OUT/intermediate/`.
   `OUT/README.txt` explains every folder. Published result paths do not change.
4. **Plot style "Accept proposals".**
   - Undetected genes in locus plots are drawn white, with a grey outline and hatch (P-01).
   - SHAP and partial dependence are hidden when the cross-validated R² is ≤ 0 or its interval
     includes 0 (P-07).
   - Every title uses the TU-block null for Moran's I, and names it (P-04).
   - The V03 index image is dropped. It is kept as a TSV plus a report table (P-12).
   - Regenerated snapshot baselines still go to Khoa for approval in each PR.
5. **F-01.** The example bundle runs `core,bgc` without `function`, so it is not affected. The
   empty or contradictory function section at coverage 0 is still fixed.
6. **V-01** (start-codon artefacts as ALT sites) came from an older K-12 run. This spec says
   "re-check on main first" and states the check as an acceptance criterion (AC 10.0).

## Evidence runs (read-only)

| Alias | Run | What it exercises |
|---|---|---|
| A | Nissle 200 kb example, route A, `core,bgc,function` | DIA-NN peptides + iBAQ, tier `mayu`, Bakta light, 158 genes, 99 detected, 774 peptides |
| B | Nissle 200 kb example, route B, `core,bgc` | proteins + MaxLFQ, tier `protein_list`, no peptides, 106 detected |
| C | Nissle 200 kb example, route C, `core,bgc` | Cascadia de novo, tier `two_peptide`, 21 peptides, 4 detected, no abundance |
| K | K-12 MG1655 full genome, `core,mobile,bgc` | 4,373 genes, 2,792 detected. **Run at commit 77df81f5 (2026-10-28), before M11-20 (#490, 2026-10-29)**. Items from K must be re-checked on main |

The run directories are listed in the review (`local/handoff/outputs_review_1_1.md` › Runs).
**R11 applies to tests.** Nothing row-level from the example bundle or K goes into `tests/data/`.
Tests build A-, B-, C- and K-like fixtures with `gippasim` or by hand: peptides + mayu,
proteins-only + `protein_list`, de novo with ≤ 5 detected genes and no abundance, and a fixture
with more than 500 markers. The values quoted below come from the real runs and are not test
expectations, unless an AC says so.

## Corrections to the review (verified on the runs, 2026-11-04)

These change what the PRs must do. The test-author should use these facts, not the review's.

| Review says | Actually | Consequence |
|---|---|---|
| `detected` is column 50, `abundance_log2` ≈ 62 (T-03) | `detected` 48, `detection_status` 50, `abundance_log2` 61, `cog_category` 68 (all four runs) | the migration table in § 5 uses the real positions |
| T-01: many tables have ≥ 8 decimals (peptide_hits 403, genes 14, …) | Counted as **> 6 significant digits** (the `%.6g` rule), only `tables/detection_summary.tsv` and `tables/regions_stats.tsv` break the rule (A: 61 and 36 values; B: 42 and 34; C: 34 and 18). `0.00012345` has 8 decimals but 5 significant digits, so it passes | PR 9 fixes two writers. The AC counts significant digits, not decimals |
| T-02: `cross_strain_protein` probably comes from `genome_relation: related` | It comes from `gippa.mapping.loci.annotate_protein_hits`: "`protein_id` is not a locus tag of the model". B's proteins are `JNOIEE_*` from a related annotation, and the genes are `EOIENK_*` | the PR 9 rule is stated in terms of the hit class, not the config |
| T-04: rename `tir_ostir`/`mw_da` in genes | `genes.tsv` already has `tir_ostir` and `mw_da`. The renames land in `novel_orfs`, `determinants_features`, `peptide_hits`, `pairs` and `operon_flags`. `determinants_features.promotech_max`/`promotech_tu` hold exactly `genes.promoter_score`/`promoter_score_tu` (`determinants/features.py` maps them) | PR 11b covers those tables; genes keeps its names |
| T-03: "Parquet twin too" | `output.parquet` (`--parquet`) is read by no writer, so there is no published Parquet twin. `genetable/genes.parquet` is an intermediate | nothing to reorder; Q11 |
| V-01: `KOIAIB_01749 aa 1 L→M` has no `start_codon` flag | K's VAR_167 (`KOIAIB_01749`, aa 1) **has** `start_codon`. The unflagged ones are the sites outside any gene: VAR_30, VAR_96, VAR_107, VAR_129, VAR_135, VAR_185, … | PR 10 targets the out-of-gene sites; M11-20 should already remove the in-gene ones |
| R-02: "version reads 1.0.0 … check it comes from `gippa.__version__`" | It does (`report/build.py` passes `gippa.__version__`, which is still `"1.0.0"` on main; release prep bumps it). K's manifest says `gippa_version` `0.0.1` | the run box shows both the builder's version and the run's version when they differ (AC 1.6) |
| D-01: genes has "Path: —" | 18 of the 30 tables in `docs/reference/tables.md` show "Path: —": the datapackage has only 12 resources | PR 3 fixes the generator for every table |
| P-07: hide when "its CI includes 0" | `determinants_model` has `ci_low`/`ci_high` = NA with `ci_unestimated` for `cv_r2_mean` in every run (B, K, the snapshot fixture) | the gating rule uses the mean ± SD the title shows (§ PR 7b; Q6) |

## Spec sections touched

- §9.1 output tree: adds `intermediate/` and `README.txt` (PR 12).
- §9.2 versioning ("breaking column changes only in major versions") and ADR-0060 §3 ("renamed,
  retyped or reinterpreted columns … flag tokens" are major). **ADR-0075** (lead, Khoa accepts)
  grants the 1.1 exception for exactly the changes in § 5. It overrides §9.2 for those changes
  only and names the `schema_diff` exception file (PR 11b).
- §9.3 column dictionaries: new `genes.tsv` group order; renamed columns; `class_error_rate`
  and `class_fdr` become `float`.
- §9.4 flag vocabulary: `detectability_constant` moves from Gene to Run. New flags:
  `class_error_unestimated` (Variants), `putative_start` (Variants, Q7), `foreign_protein_ids`
  (Run).
- §10.1 colours: undetected genes in V3/V13 locus views are white with a `#BBBBBB` outline and
  hatch. The §10.1 table value is kept for every other plot. After M9 approval this is a reviewed
  style change, and Khoa decided it (decision 4).
- §10.2: V3 no longer writes an index image (P-12). V12 (b) shows at most 15 sites, not 30 (P-17).
  V15 hides (b) and (c) when the model has no predictive power (P-07). ADR-0075 records all three.
- §10.3: report sections. The report gains a run box, a TOC and a "Not applicable in this run"
  list. Section 3 is renamed "Mapping". The cards adapt to the inputs. The methods are chosen
  per step and tier.

## What already exists (reuse, do not duplicate)

| Piece | Where | Reused for |
|---|---|---|
| `%.6g` TSV writer | `gippa.io.tsv` | PR 9 routes the two offending writers through it |
| empty state, `Skip`, `title_pct` | `gippa.viz.base` | PR 0 adds `Instance.empty`; PR 4 rewrites messages |
| colour constants | `gippa.viz.style` (`DETECTION`, `STATUS_COLOUR`, `REGION`, …) | PR 5 adds the locus "undetected" style, not a new palette |
| binning helper | `gippa.viz.base.binned(…, max_elements)` | PR 7a calls it with `max_elements=500` for marker tracks |
| run flags, modules off | `RunView.run_flags()`, `RunView.modules_off()` | skip reasons, "not computed" notes, the Not applicable list |
| read fallbacks | `gippa.io.run_view._FALLBACK` | PR 12 keeps 1.0 output trees readable |
| methods snippets | `scripts/dev/sync_report_snippets.py`, `src/gippa/report/snippets/` | PR 3 marks more `docs/methods` blocks |
| flag descriptions | the schema `description`/flag reference (`scripts/docs/gen_reference.py render_flags`) | PR 0's plain-text dictionary starts from them |
| detection_summary overall row `detectability / overall / all / detectability_constant / 1` | `gippa.genetable.detectability` | T-06 needs no new row there, only the manifest flag |
| manifest `flags` | `run_manifest.json › flags` (e.g. B: `anti_sd_default; bakta_db_light; cai_reference_missing; layer_missing:tx_rates; no_detection_limit_model`) | T-06 and T-02 run-level flags |
| determinants minimum | `gippa.determinants.models.MIN_ROWS = 100`, `MIN_TUS = 25` | P-14 states the minimum in the message |
| function coverage | event `function.collections` in `logs/steps/function.events.jsonl` (A: every collection `coverage 0.0, n_terms 0`) | F-01 reads coverage from it |
| schema diff test | `scripts/release/schema_diff.py`, `tests/unit/test_m11_14h_schema_diff.py::test_this_repository_is_schema_compatible_with_its_last_release` (tag `v1.0.0`) | PR 11b adds the ADR-0075 exception file |

---

## 1. Rules for every PR

1. **Numbers in figure text and in the report** go through `gippa.viz.style.fmt_q` and `fmt_num`
   (PR 0). Tables keep `%.6g`.
2. **No internal references in user-facing text.** No `ADR-\d{4}`, no `§\d`, no
   `GIPPA_SPEC.md`, no `docs/methods` path, no raw flag token or skip code in prose. Prose means
   figure titles, captions, cards, notices, notes and the README. Table cells and data
   attributes may still hold flag tokens. Every PR that touches prose adds its strings to the
   PR 0 dictionary rather than inlining them.
3. **No site names or internal paths** in any user-facing file: report, README, changelog
   fragments, docs. A test greps generated files of the fixtures for `farmshare`
   (case-insensitive), `/home/`, `local/` and `$GIPPA_ROOT`.
4. **Baselines.** A PR regenerates only the baselines listed for it in § 2. Any other baseline
   that changes is a review failure. The PR body lists every regenerated baseline, and Khoa
   approves them before merge (§10.1, M9-25). Agents never mark a baseline approved.
5. **Changelog.** Each PR touching `src/` or `workflow/` adds `changelog.d/<PR>.<type>.md`:
   `feature`, `bugfix` or, for PRs 11a–11c and 12, `removal` (§ 5).
6. **Roles.** The test-author writes the failing tests from the ACs first, and the implementer
   follows. Existing tests whose expectations change are listed per PR and updated by the
   test-author.

---

## 2. PR list

| PR | Title | Items | Code ≈ | Depends on | Baselines regenerated | Goldens / fixtures changed |
|---|---|---|---|---|---|---|
| 0 | `feat(viz): shared wording and number helpers` | P-18 (helpers), P-13 (label helper), R-05 (dictionary), P-14 (`Instance.empty`) | 250 | none | none | none |
| 1 | `feat(report): run box, input-aware summary cards, TOC` | R-01, R-02, R-03 | 350 | 0 | `report_structure.json` | none |
| 2 | `feat(report): captions, plain wording, curated previews, conditional sections` | R-04, R-05, R-06, R-08, P-14 (report half), P-04 (report label) | 380 | 0, 1 | `report_structure.json` | none |
| 3 | `feat(report): methods per step and tier; reference page fixes` | R-07, D-01 | 260 | 2 | `report_structure.json` | `docs/reference/tables.md`, `flags.md` (generated) |
| 4 | `fix(viz): degenerate and empty states` | P-02, P-03, P-05, P-06, P-11, P-12, P-15, P-14 (plot half) | 380 | 0, 5 | `sim_empty-*` of V05, V06, V10, V11, V12, V13, V15; `sim_plots-V11`; **removes** `sim_plots-V03_locus.png`, `sim_empty-V03_locus.png` | none |
| 5 | `fix(viz): palette, labels and number style` | P-01, P-13, P-18 (call sites), P-19, P-04 (titles) | 330 | 0, 11b | every `sim_plots-*` whose text has a q, a minus or an abundance axis, and the locus views: V01, V02, V03_locus_ISL_1, V03_locus_prophage_1, V07, V08, V10, V11, V13_novel_orfs_locus_NOVEL_* (2), V14, V15, V18 | none |
| 6 | `feat(viz): V07 layout, V18 small-n, V17 labels` | P-08, P-10, P-16 | 360 | 5 | `sim_plots-`/`sim_empty-` V07, V17, V18 | none |
| 7a | `feat(viz): V01 ring key and marker binning` | P-09 | 260 | 5 | `sim_plots-V01`, `sim_empty-V01` | none |
| 7b | `feat(viz): V15 gating and labels` | P-07 | 200 | 5 | `sim_plots-V15`, `sim_empty-V15` | none |
| 8 | `fix(function): one notice when no term was assigned` | F-01 | 170 | 2, 5 | `sim_empty-V08`, `sim_empty-V14` (when their messages change) | `report_structure.json` |
| 9 | `fix(tables): %.6g everywhere; protein cross-strain flag` | T-01, T-02 | 200 | none | none | `sim50k/m4/detection_summary.tsv`, `sim50k/m4/regions_stats.tsv`, `sim100k_real/m4/detection_summary.tsv`, `sim100k_real/m4/regions_stats.tsv` (values only) |
| 10 | `fix(variants): out-of-gene start artefacts; I/L rendering; V12 top 15` | V-01, P-17 | 250 | AC 10.0 run, 11c, 5 | `sim_plots-V12`, `sim_empty-V12` | `sim50k_m7/variants.tsv`, `variants.vcf` only if the sim has such sites (the test-author checks) |
| 11a | `feat(tables)!: genes.tsv column groups; run-level flags` | T-03, T-06 | 220 | ADR-0075 PR | none | `sim50k/m4/genes.tsv`, `sim100k_real/m4/genes.tsv` (order; `detectability_constant` off rows); every `tests/data/snapshot/*/OUT/tables/genes.tsv` and `CMP/runs/*/tables/genes.tsv` fixture (order only) |
| 11b | `feat(tables)!: one name per quantity` | T-04 | 300 | ADR-0075 PR | `sim_plots-V15` (the PD label, § PR 11b) | `sim50k_m7/novel_orfs.tsv`, `sim50k/mapping/peptide_hits.tsv` (headers); snapshot fixtures `novel_orfs.tsv`, `peptide_hits.tsv`, `pairs.tsv`, `operon_flags.tsv`, `determinants_features.tsv`, `determinants_model.tsv`, `plots/support/determinants_pd.tsv` |
| 11c | `feat(tables)!: numeric class_error_rate and class_fdr` | T-05 | 180 | 11b | none | `sim50k_m7/variants.tsv`, `sim50k_m7/novel_orfs.tsv` (`unestimated` → NA, flag added); snapshot fixtures `variants.tsv`, `novel_orfs.tsv` |
| 12 | `refactor(outputs)!: intermediates under OUT/intermediate/, OUT/README.txt` | O-01, D-01 (outputs.md) | 300 | none | none | none (published paths are unchanged) |
| 13 | `feat(report): detection views; shareable reports; report size` | T-07, R-09 | 220 | 2 | `report_structure.json` | none |

Every item of the review is covered exactly once. The review's PR 7 is split into 7a and 7b, and
its PR 11 into 11a, 11b and 11c. PR 0 is new: it holds the helpers that PRs 1, 2, 4, 5, 6, 7b, 8
and 10 share, so they can run in parallel without duplicating a dictionary.

---

## 3. PRs

Each PR section gives the behaviour, before → after examples from the evidence runs, the edge
cases, and the acceptance criteria (AC *n.m*). "Fixture P/R/D/K" means the test-author's
simulated equivalents of A (peptides, mayu), B (proteins, protein_list), C (de novo, few genes,
no abundance) and K (large: more than 500 markers, many ALT sites).

### PR 0: shared wording and number helpers

**Behaviour.** This PR adds pure helpers only. No call site changes, so no image or HTML changes.

- `gippa.viz.style.fmt_q(q) -> str`:
  - q < 0.001 → `"< 0.001"`;
  - otherwise 2 significant digits, keeping trailing zeros: 0.432 → `"0.43"`, 1 → `"1.0"`, 0.05 →
    `"0.050"`;
  - NA → `"NA"`.
- `fmt_num(x, sig=2) -> str`:
  - significant digits, with a true minus `−` (U+2212): −0.147806 → `"−0.15"`;
  - `sig=3` for log2 abundances (14.2);
  - integers ≥ 10,000 get a thin-space or comma grouping (`4,373`);
  - NA → `"NA"`.
- `fmt_pct(n, d)` replaces `base.title_pct`, which stays as an alias: `"63 %"`, and d = 0 → `"NA"`.
- `abundance_label(type, scale) -> str`. The table below gives the outputs:

  | `abundance_type` | `abundance_scale` | label |
  |---|---|---|
  | `ibaq` | `molar` | `log2 iBAQ (molar)` |
  | `maxlfq` | `relative` | `log2 MaxLFQ (relative)` |
  | `top3` | any | `log2 Top3 (<scale>)` |
  | `intensity` | any | `log2 intensity (<scale>)` |
  | `unknown` | any | `log2 abundance (measure unknown)` |
  | `none` or NA | – | `None` (the caller draws no abundance axis) |

- `gippa.viz.words` (new module, pure):
  - `FLAG_TEXT`: flag token → short sentence. It covers every token of every schema `flags`
    pattern plus the run flags, starting from the flag reference descriptions and rewritten as
    plain words, for example `layer_missing:codon` → "codon usage not computed", `c_assumed` →
    "de novo precision assumed, not measured".
  - `SKIP_TEXT`: skip-reason prefix → sentence, for example `module_off:function` → "the
    function module was not run", and `input_absent:tables/variants.tsv` → "no variant table in
    this run".
  - `FEATURE_LABEL`: feature column → human name, for example `tir_ostir` → "translation
    initiation rate (OSTIR)", `cpb` → "codon pair bias", `at_2_10` → "AT content, codons 2–10",
    `terminator_conf_downstream` → "downstream terminator confidence".
  - `TIER_TEXT`: tier → sentence; see PR 1.
  - `plain(token)` returns the sentence, or the token itself with a logged warning when it is
    unknown. It never raises.
- `gippa.viz.base.Instance` gains `empty: bool = False`, and `empty_state()` sets it to `True`.

**Edge cases.** q exactly 0.001 gives `"0.0010"`. q = 0 gives `"< 0.001"`. A negative zero gives
`"0.0"`, with no minus. A flag token with a parameter (`layer_missing:<layer>`) is looked up by
its prefix and then filled in.

**AC.**
- 0.1: `fmt_q` and `fmt_num` return the strings above. A Hypothesis property: the output never
  holds the ASCII hyphen-minus as a sign, and it parses back, with `−` replaced by `-`, to within
  the significant-digit rounding.
- 0.2: every token allowed by the `flags` pattern of every `src/gippa/schemas/tables/*.schema.json`,
  and every Run flag of §9.4, has an entry in `FLAG_TEXT`. A parametrised test runs over the
  schema patterns, so a new flag without a sentence fails CI.
- 0.3: no `FLAG_TEXT`, `SKIP_TEXT` or `TIER_TEXT` value matches `ADR-\d|§|GIPPA_SPEC|_` (no
  underscores, so no raw tokens).
- 0.4: `empty_state(...).empty is True`, and an ordinary `Instance` has `empty is False`.
- 0.5: `abundance_label` gives the table above, and `("none", NA)` gives `None`.

### PR 1: run box, input-aware summary cards, TOC (R-01, R-02, R-03)

**Behaviour.**

1. **"This run" box**, the first element after the `<h1>`, `id="this-run"`. It has these lines:
   - genome: the stem of the genome file name, with `.gz`, `.fna`, `.fa` or `.fasta` removed;
   - replicons: count and total length;
   - genes;
   - input kind: `peptides (<n>)`, `proteins (<n>)` or `de novo peptides (<n>)`, from the
     manifest inputs and `peptide_hits.source`;
   - abundance: `abundance_label`, or "none";
   - modules;
   - the detection tier in plain words, from `TIER_TEXT` with the run's alpha:

     | Tier | Sentence |
     |---|---|
     | `decoy_picked` | "decoy_picked: a gene counts as detected when its gene-level q-value (target–decoy competition) is ≤ 0.01." |
     | `mayu` | "mayu: a gene counts as detected when it has an accepted peptide at q ≤ 0.01; the gene-level FDR of that list is estimated with Mayu." |
     | `two_peptide` | "two_peptide: a gene counts as detected when ≥ 2 distinct accepted peptides map to it; no FDR is estimated." |
     | `protein_list` | "protein_list: a gene counts as detected when a protein in your list maps to it; no FDR is estimated." |

   Below the lines come **up to 3 key findings**, the finding titles of V01, V07 and V02 taken
   verbatim from their `Instance.title`. A plot that was skipped or is empty contributes no
   sentence.
2. **Version line.** "Report built by GIPPA <`gippa.__version__`>". When the manifest
   `gippa_version` differs, the line adds "from a run made with GIPPA <x>".
3. **Cards** adapt to the inputs:
   - The peptide card appears only when the run had peptides, and the protein card only when it
     had proteins. A "0 of 0" card is never shown.
   - The "also pass the two-peptide rule" line appears only in tiers where it is a separate
     statement: `mayu` and `decoy_picked`.
   - The mayu card says what the FDR is: "Mayu gene-level FDR 0.012, from peptides at
     q ≤ 0.01".
   - Numbers go through `fmt_num` (3 significant digits on cards) and `fmt_q`.
   - The per-stratum FDR string moves off the cards. It stays in V06 and `mapping_fdr.tsv`.
   - Model R² card:
     - label first: "Model R²: −0.017 ± 0.066 (CV)";
     - "no predictive power" is added when the mean is ≤ 0 or mean − SD ≤ 0 (the PR 7b rule);
     - when no model ran: "not fitted: needs ≥ 100 genes with abundance in ≥ 25 TUs (this run:
       99 genes, 34 TUs)", or "not fitted: no abundance in this run".
4. **TOC.** A `<nav id="toc">` comes after the run box. It has one link per section, `#sec-<key>`,
   and a count in brackets: the row count of the section's main table, for example
   "Variants (0)" or "Gene detection (158)". It is sticky through CSS only, with no script.

**Examples (before → after).**

| Run | Before | After |
|---|---|---|
| A | page opens "GIPPA report / GIPPA 1.0.0. Every number is shown…" then cards | run box: "Genome nissle_r200k · 1 replicon, 200,000 bp · 158 genes · peptides (774) · log2 iBAQ (molar) · modules core, bgc, function · mayu: a gene counts as detected when…"; findings: "99 of 158 genes detected (63 %); 0 islands at FWER ≤ 0.05; 0 ALT sites, 0 novel ORFs" (V01), V07's and V02's titles |
| A | card "99 genes detected (Mayu FDR 0.0124093 at q ≤ 0.01) 91 also pass the two-peptide rule; 158 genes" | "99 of 158 genes detected · Mayu gene-level FDR 0.012, from peptides at q ≤ 0.01 · 91 with ≥ 2 distinct peptides" |
| A | card "0 proteins placed of 0 input unestimated" | no protein card |
| B | "106 genes detected unestimated 0 also pass the two-peptide rule; 158 genes" | "106 of 158 genes detected · protein list: no FDR estimated" (no two-peptide line) |
| B | "0 peptides placed of 0 input …" | no peptide card; "106 of 106 proteins placed" |
| B | "Model R² -0.0174037 ± 0.0656274 (CV); null p 0.414585 (uncalibrated)" | "Model R²: −0.017 ± 0.066 (CV), no predictive power; null p 0.41 (uncalibrated)" |
| C | peptide card estimate: "placement q ≤ 0.01 per family (0 / ≥ 1 difference, per decoy type); descriptive stratum FDR (not q-values): shuffle 0 diff 7: 0, shuffle 0 diff 8: 0, …" (≈ 300 characters) | "21 of 21 peptides placed · placement q ≤ 0.01" |
| A, C | "Model R² NA (CV); unestimated" | A: "Model R²: not fitted: needs ≥ 100 genes with abundance in ≥ 25 TUs (this run: 99 genes, 34 TUs)"; C: "Model R²: not fitted: no abundance in this run" |
| K | manifest `gippa_version` 0.0.1, report built later | "Report built by GIPPA 1.1.0 from a run made with GIPPA 0.0.1" |

**Edge cases.**
- Peptides and proteins both given: both cards are shown.
- A de novo run (C) uses "de novo peptides (21)" as the input kind.
- A manifest without `inputs` (an old run) gives "input: not recorded".
- No V01 title (the plot was skipped) gives findings from V07 and V02 only. If none of the three
  is available, the findings list is omitted, not left empty.
- An empty module set (`--modules core`) lists "core".
- A section whose table is absent has no count, rather than "(0)".

**AC.**
- 1.1: on fixture R (proteins only), `cards()` returns no card whose name is `peptides`, and no
  card's text contains "two-peptide". On fixture P (peptides only), there is no `proteins` card.
- 1.2: no card text on fixtures P, R and D matches `\b0 (peptides|proteins) placed of 0\b`.
- 1.3: on fixture P (mayu), the gene card contains "gene-level FDR" and "q ≤ 0.01", and every
  number in it has ≤ 3 significant digits.
- 1.4: no card is longer than 160 characters on any fixture. On fixture D, the card text does not
  contain "stratum".
- 1.5: the R² card starts with "Model R²:". On a fixture with fewer than `MIN_ROWS` genes, it
  contains "not fitted" and both minimums. With `abundance_type = none`, it contains "no
  abundance".
- 1.6: `#this-run` exists on every fixture, and it holds the tier sentence from `TIER_TEXT` for
  the fixture's tier. With manifest `gippa_version` ≠ `gippa.__version__`, it holds both
  versions.
- 1.7: `nav#toc` holds one `a[href="#sec-<key>"]` per `<section>` in document order, and every
  target exists. This is checked by the existing link checker.
- 1.8: `tests/snapshot/report_structure.json` is regenerated, and the diff adds only `nav`,
  `#this-run` and the section ids.

### PR 2: captions, plain wording, curated previews, conditional sections (R-04, R-05, R-06, R-08, P-14 report half, P-04 report label)

**Behaviour.**

1. **Captions (R-04).** Each `<figcaption>` starts with the §10.2 "Question" of the plot, then one
   "how to read" sentence from a `CAPTIONS` mapping in `gippa.report.content`. The mapping has 16
   entries, one per plot ID, and the links (SVG PNG PDF data TSV) follow. For example, V02:
   "Is abundance clustered along the chromosome? Each point is a detected gene's log2 abundance
   at its position; the line is a LOESS trend; shaded bands are abundance islands."
2. **Plain wording (R-05).**
   - ADR and § references become a short caveat plus a link to the docs methods page. The link is
     an absolute `https://khoa-yelo.github.io/gippa/methods/<topic>/` `<a href>` (Q12).
   - Flags and skip reasons in prose go through `plain()`. The report's decoy table headers
     `N_NO`, `N_cv`, `N_sw` become "non-operon pairs", "convergent pairs" and "sandwich pairs".
   - The sentence "generated from the GIPPA 1.0.0 documentation (docs/methods)" is removed.
   - The Methods section gets a glossary `<dl id="glossary">` with: FDR, FWER, q-value, local
     FDR (lfdr), O − E, π₀, block length (ℓ̂) and BB join count.
3. **Curated genes preview (R-06).**
   - Columns: `locus_tag`, `gene`, `product`, `detection_status`, `n_peptides_unique`,
     `abundance_log2`, `region_ids`, `flags`.
   - Sort: by `abundance_log2` descending, NA last, ties by `locus_tag`. When there is no
     abundance, sort by `detected` descending, then `n_peptides_unique` descending, then
     `locus_tag`.
   - The link to the full TSV stays. The table itself does not change.
4. **Conditional sections (R-08).**
   - The heading "Peptide mapping" becomes "Mapping".
   - The Determinants caveats print only when a model was fitted (`determinants_model` has a
     finite `cv_r2_mean`).
   - Sensitivity: when every family has `n_primary = n_variant = 0`, one sentence replaces the
     tables: "No flagged results, so there is nothing to compare across sensitivity reruns."
   - 0-row table previews are hidden, and the section shows "none" with the link.
   - Inapplicable sections (module off, or a needed input absent) move to one
     `<section id="sec-not-applicable">` at the end. It lists "<Section>: <plain reason>", for
     example "Variants: no variant table in this run (variant calling needs peptides)".
5. **Empty states in the report (P-14 report half).** An `Instance` with `empty=True` is shown as
   a `<p class="notice">` holding its title, not as an embedded image. The image files are still
   written to `plots/`.
6. **Island diagnostics label (P-04 report half).** "Î₁ (lag-1 Moran's I)" becomes "lag-1
   Moran's I of the detrended TU scores (used to choose the block length)".

**Examples.**

| Run | Before | After |
|---|---|---|
| A | figcaption "V17_qc SVG PNG PDF data TSV" | "Are the inputs sane? Panels show the preflight summary, contig lengths, peptide lengths, abundances and missing values per feature. SVG PNG PDF data TSV" |
| A | "ADR-00xx" appears 6 times in the report; "§7.7.1" once | 0 and 0 |
| A | genes preview: 20 rows × 73 columns in table order | 20 rows × 8 columns, highest `abundance_log2` first |
| C | genes preview sort (no abundance) | the 4 detected genes first |
| B | heading "Peptide mapping" | "Mapping" |
| C | Determinants: 6 caveats although no model ran | no caveats; the section holds the V15 notice |
| A | "islands.tsv (first 0 rows; full table)" with an empty sortable table | "Islands: none (full table)" |
| B, C | "V12: not run: input variants absent" inline | listed under "Not applicable in this run: Variants: no variant table in this run" |
| A | V15 image "Too few genes for the determinants model (99 rows, 34 TUs)" | `<p class="notice">` with the PR 4 message |
| A | "Î₁ (lag-1 Moran's I) -0.399" | "lag-1 Moran's I of the detrended TU scores (used to choose the block length) −0.40" |

**Edge cases.**
- A run with every optional module off: the not-applicable list holds Function, and Mobile
  when mobile is off.
- The `mobile` module off, as in A, B and C, gives a "prophage and plasmid regions not searched"
  reason in the Genome context notes, not a whole section.
- Empty-module runs (`--modules core`): Function and the BGC-specific content are listed as not
  applicable.
- The no-KEGG case (#543): the Function section stays, with COG and GO only. This interacts with
  PR 8.

**AC.**
- 2.1: in `report.html` of every fixture, the text of `p, li, figcaption, .card, .notice, h2, h3,
  dd` matches none of `ADR-\d{4}`, `§\s?\d`, `GIPPA_SPEC`, `docs/methods`. It also contains no
  token of the flag vocabulary as a whole word. (Lead note, 2026-11-04: "token" means a
  vocabulary key containing `_` or `:`, such as `top3_lt3` or `uncalibrated:operon_distance`;
  plain English keys such as `partial` or `pseudogene` may appear in prose.)
- 2.2: every `figure` has a `figcaption` whose text starts with its plot's §10.2 Question
  (`CAPTIONS` keys = registered plot IDs, tested against `viz.base.registered()`).
- 2.3: the genes preview has exactly the 8 columns in that order, and is sorted as stated, on
  fixtures P and D.
- 2.4: on fixture R, no `h2` text is "Peptide mapping".
- 2.5: on a fixture with no fitted model, the Determinants section has no `.caveat`.
- 2.6: no `table` with 0 body rows is rendered. Each hidden one has the text "none" and a working
  link.
- 2.7: `#sec-not-applicable` exists if and only if at least one section is inapplicable, and it
  lists each with a non-empty plain reason.
- 2.8: an `empty=True` instance is rendered as `.notice`, with no `<svg>` or `<img>` for that plot.
- 2.9: `#glossary` defines the 8 terms.

### PR 3: methods per step and tier; reference page fixes (R-07, D-01)

**Behaviour.**

1. **Methods snippets (R-07).**
   - `docs/methods/` pages gain `<!-- report:begin [key] -->` blocks of ≤ 120 words, in
     `detection-tiers.md` (one keyed block per tier), `placement-fdr.md`, `denovo-precision.md`,
     `abundance.md`, `operons.md`, `variants.md`, `novel-orfs.md`, `function.md` and
     `determinants.md`. The existing `islands.md` block stays.
   - `sync_report_snippets.py` accepts an optional key, which writes `<topic>.<key>.html`.
   - The report includes the snippets of the steps that ran (manifest `steps`) and of the run's
     tier. Each snippet ends with "Full note: <link>".
2. **Reference pages (D-01).**
   - `gen_reference.py` gives a path for every table. The datapackage `path` comes first, then
     `gippa.io.run_view.PATHS`. The compare tables are `<compare outdir>/tables/compare_*.tsv`.
   - Descriptions are rendered with "(GIPPA_SPEC.md §…)" parentheticals removed. The schema files
     are not edited, which avoids conflicts with PRs 11a–11c.
   - `plots-and-report.md` "GIPPA 1.0 builds 16 plots" becomes "GIPPA builds 16 plots" (V1–V3,
     V5–V8 and V10–V18).

**Examples.** A: the Methods section held only "islands". After the change it holds the tier
snippet (mayu), placement FDR, abundance, operons, function and islands. C: the tier snippet
(two_peptide), placement FDR, de novo precision and islands, with no abundance snippet.
`tables.md › genes`: "Path: —" → "Path: `tables/genes.tsv`".

**Edge cases.**
- A step with no snippet contributes nothing.
- Tier `decoy_picked` has its own block.
- An old run whose manifest has no `steps` falls back to the snippets of the tables present.

**AC.**
- 3.1: `sync_report_snippets.py --check` passes. Every keyed block has ≤ 120 words.
- 3.2: on fixture D, the Methods section holds the `detection-tiers.two_peptide` snippet and
  `denovo-precision`, and no `abundance` snippet. On fixture P, it holds `detection-tiers.mayu`,
  and not `detection-tiers.two_peptide`.
- 3.3: `docs/reference/tables.md` has no "Path: —". It contains no `GIPPA_SPEC.md`, and the
  generator `--check` passes in CI.
- 3.4: no page under `docs/user-guide/` contains `GIPPA 1\.0 builds`.

### PR 4: degenerate and empty states (P-02, P-03, P-05, P-06, P-11, P-12, P-15, P-14 plot half)

**Behaviour and examples.**

| Item | Plot | Condition | Before (run) | After |
|---|---|---|---|---|
| P-02 | V06 (a)(b) | every `mapping_fdr` placement row has `targets = 0` (no tolerant targets) | A: "0 FDR-tested placements accepted; strata 0, 1, 2 unestimated (no targets)"; (b) three overlapping "NA (no targets)" labels; (a) log axis with 2 decoys | A: if there are no de novo rows either, the whole plot is an empty state: "All 774 peptides placed by exact match; the tolerant-route FDR is not needed." With de novo rows (C), (a) and (b) collapse to one text line and (c) is drawn. When (b) does draw, NA strata are listed in one stacked note, not as overlapping labels |
| P-03 | V06 (c) | every `denovo_precision` row has flag `c_assumed` | C: 10 bins each at precision 1.00 with CI 1–1, axis "precision (CI)", footnote "c assumed 0.99 (no reference peptides)" | no points; the panel holds the text "De novo precision not measured (no reference peptides); assumed c = 0.99". The axis label "(CI)" appears only when a CI is drawn |
| P-05 | V02 | `abundance_type = none` or no finite `abundance_log2` | C: empty log2 axis −0.04 to 0.04 with a rug | a detection strip (one row of ticks coloured by `detection_status`, no y axis), titled "No abundance in this run: detection along the replicon (4 of 158 detected)" |
| P-06 | V10 | every `morans_i` row of `autocorrelation.tsv` is NA or flagged `not_estimable`/`few_values` | C: "not estimable (ci_unestimated, few_values, not_estimable, tu_block_tu_error)" in (a), and (b) drawn from 4 genes with spikes at lags 10 and 14 | the whole plot is an empty state: "Too few detected genes (4 of 158) to estimate autocorrelation." Flags in partial panels go through `plain()` |
| P-11 | V11 | every panel would be empty | C: three panels "no measured pairs / no lfdr values / no measured same-strand pairs" | one empty state: "No measured operon pairs: partner comparisons need abundance." |
| P-11 | V11 | title vs drawn | A: title "too few measured pairs for densities" while (a) draws densities; (c) titled "partners (OP, NO)" shows only OP; legend "(c) OP" | the title is built from what is drawn; (c)'s title names only the classes drawn; the legend entry is "OP (same operon)"; NO gets `#009E73` (not the detection "unestimated" orange `#E69F00`); OP is drawn on top with alpha 0.6 |
| P-12 | V03 | always | `plots/V03_locus.png` (text-only image "Order alternates between test families…") | no `V03_locus.{svg,png,pdf}`; `plots/data/V03_locus.tsv` is still written; the report shows it as a table "Locus views, in display order" with the columns rank, locus, family, estimate, q |
| P-15 | V05 | exactly one route and one outcome with non-zero flow | A: two identical blue bands at full height ("100 % of 774 input peptides placed in the genome"); B: same for 106 proteins | one horizontal bar at ⅓ of the height, titled "All 774 peptides placed by exact match" (A) / "All 106 proteins placed (miniprot)" (B) |
| P-14 | V12–V15 empty states | any empty state | A V15 "Too few genes for the determinants model (99 rows, 34 TUs)"; A V14 "No GSEA or detection ORA rows in this run" | the pattern "Nothing to show: <reason>. <what would change it>", e.g. "Nothing to show: the determinants model needs ≥ 100 genes with abundance in ≥ 25 TUs (this run: 99 genes, 34 TUs)." V14 at coverage 0 is PR 8's |

**Edge cases.**
- Fixture K-like with tolerant targets: V06 is unchanged apart from the stacked NA note.
- A run with peptides and proteins: V05 has two routes, so it draws the full plot.
- A linear contig in V02 strip mode: positions run 1…length with no wrap. Origin-crossing genes
  in the strip use the M9 convention, a tick at the start.
- Minus-strand genes in the strip: one row only. The strand is not encoded, because the strip
  shows detection.
- In V10, a partly estimable correlogram (some lags estimable) keeps both panels, and the NA lags
  are noted with `plain()`.
- No V13 change here (Q5).

**AC.**
- 4.1: on a fixture with exact-only placements and no de novo rows, V06 returns a single
  `Instance` with `empty=True`, whose title contains "placed by exact match" and the peptide
  count.
- 4.2: on fixture D, V06 panel (c) draws no markers (the data TSV `panel = c` rows keep the
  precision values and add `drawn = false`). No axis label contains "(CI)" unless a CI is drawn.
- 4.3: on fixture D, V02's data TSV has no `abundance_log2` column values, and the title matches
  `No abundance in this run: detection along the replicon \(\d+ of \d+ detected\)`.
- 4.4: on fixture D, V10 is one empty instance whose title holds the detected count. No figure
  text contains `tu_block_tu_error` or any raw flag token.
- 4.5: V11 on fixture D is one empty instance. On fixture P, the legend labels are exactly the
  classes drawn, and "OP" appears as "OP (same operon)". The NO colour is not in
  `style.DETECTION.values()`.
- 4.6: `registered()["V03"]` returns no `Instance` with `suffix is None`, and
  `plots/data/V03_locus.tsv` is still written. The report holds a table built from it. The
  baselines `sim_plots-V03_locus.png` and `sim_empty-V03_locus.png` are deleted, and the snapshot
  test no longer lists them.
- 4.7: V05 with one route and one outcome has figure height ≤ ⅓ of the default V05 height.
- 4.8: every empty-state title of V12–V15 starts with "Nothing to show:" and contains one more
  sentence.

### PR 5: palette, labels and number style (P-01, P-13, P-18 call sites, P-19, P-04 titles)

**Behaviour.**

- **P-01.** In V03 and V13 locus views, a gene with `detection_status ∈ {not_detected,
  partial_excluded}` is drawn with facecolor `white`, edgecolor `#BBBBBB` and hatch `////` in
  `#BBBBBB`. `single_peptide_unestimated` keeps the `#E69F00` outline. Detected genes keep the
  cividis fill. The legend entry "not detected" shows the new style. `style.LOCUS_UNDETECTED`
  holds the three values.
- **P-13.** Every abundance axis (V02, V03 colour bar, V08, V11, V15) uses `abundance_label`.
- **P-18.** Every title and in-figure number uses `fmt_q`/`fmt_num`.
- **P-19.** In V08 the "none" category is `#E0E0E0` with hatch, and "poorly characterised" (S)
  keeps the COG grey. The rotated counts get an axis note "detected / total". The title's n is
  the same count as the label over the column.
- **P-04.**
  - V02 and V10 titles report the lag-1 Moran's I with the **TU-block** null and name it.
  - The gene-label row stays in `autocorrelation.tsv` and in V10's envelopes.

**Examples.**

| Run | Before | After |
|---|---|---|
| A | V03_locus_bgc_1: genes at 13–14 log2 near 40–45 kb drawn mid-grey, the same as undetected `#BBBBBB` | undetected: white, grey outline, hatched; the 13–14 log2 genes keep their cividis colour |
| A | V02 title "lag-1 Moran's I -0.15 (q = 0.73)" (gene-label null) | "lag-1 Moran's I −0.15 (q = 0.92, TU-block null)" |
| A | V10 title "Lag-1 Moran's I -0.15: within the TU-block envelope (q = 0.92)" | "Lag-1 Moran's I −0.15: within the TU-block envelope (q = 0.92, TU-block null)" |
| A / B / K | axis "log2 abundance" | "log2 iBAQ (molar)" / "log2 MaxLFQ (relative)" / "log2 iBAQ (molar)" |
| any | "best q 4.9e-07" (V13 locus); "q = 0.432" (V03) vs "q = 0.43" (V07) | "best q < 0.001"; "q = 0.43" in both |

**Edge cases.**
- No abundance (C): V03 has no colour bar, every detected gene is drawn in `#0072B2`, and
  undetected genes are white and hatched.
- The minus strand: arrows point left, and the fill rule is the same.
- An origin-crossing gene in a locus view keeps its M9 split.
- A TU-block row missing from `autocorrelation.tsv` (old runs): the title falls back to the
  gene-label q and names "gene-label null".

**AC.**
- 5.1: in the V03 and V13 locus data TSVs, the `fill` column of undetected genes is `white`, and
  the drawn patch has `hatch == "////"` and edgecolor `#BBBBBB` (checked on the figure artists).
- 5.2: no figure text of any plot on fixtures P and K contains a q-value written with more than
  2 significant digits, the ASCII `-` before a digit, or `e-0`.
- 5.3: the V02 and V10 titles on fixture P contain "TU-block null", and their q equals
  `fmt_q` of the `tu_block` lag-1 `morans_i` q.
- 5.4: every abundance axis label equals `abundance_label(type, scale)` of the fixture.
- 5.5: in V08, the title's n equals the count printed over its column.

### PR 6: V07 layout, V18 small-n, V17 labels (P-08, P-10, P-16)

**Behaviour.**

- **V07 (P-08).**
  - The figure is double width (7.2 in).
  - The groups get bold headers: Replicon, Feature class, COG category, Region type, Region.
  - The x axis is O − E in percentage points of the group's eligible genes, (O − E) / n × 100,
    with the CI scaled the same way. The data TSV adds `o_minus_e_pp`, `ci_low_pp` and
    `ci_high_pp`.
  - COG rows show "<letter> <short name>", from `cog_categories.tsv`.
  - Rows with n = 0 are dropped.
  - The title counts the tested items by kind.
- **V18 (P-10).**
  - Families with n < 10 p-values are not drawn. A note lists them: "Not shown (fewer than 10
    tests): region O − E (2)".
  - Panel titles use human names, for example "Region tests (O − E)" and "Join count, TU-block
    null".
  - Count ticks are integers.
  - Contig names in titles are shortened with a middle ellipsis to 24 characters.
  - With no calibration curve, the reliability panel collapses to one text line.
  - Panel letters use the "(a)" style.
- **V17 (P-16).**
  - Panel (e) uses `FEATURE_LABEL`, and features that are 100 % missing carry the reason from
    the run flags, for example "not computed: transcription-rate layer off" for
    `layer_missing:tx_rates`, and "not computed: no CAI reference" for `cai_reference_missing`.
  - Plurals are fixed: "1 preflight warning".
  - In panel (b), a contig label longer than the bar goes inside the bar or is shortened.
  - The title states both quantities: "65 % with measured abundance (63 % detected)".

**Examples.**

| Run | Before | After |
|---|---|---|
| A | V07 x label clipped "observed − expected detected ge…"; title "No region departs from detectability at q ≤ 0.05; 2 regions tested" | full label "observed − expected detected (percentage points)"; title "No group departs from detectability at q ≤ 0.05 (2 regions, 3 classes tested)" (counts from the data) |
| C | tick labels "-2.50.02.5" | non-overlapping ticks (double width, `MaxNLocator(5)`) |
| K | replicon row n = 4,367 with CI ±70 squashes every category; "crispr_2 (n = 0)" drawn | all rows on one pp scale; no n = 0 row |
| A | V18 "regions_stats hl_shift, n = 2" histogram with ticks 0.25, 0.75 | not drawn; listed in the note |
| A | V18 last title clipped "nissle1917_NZ_CP058217_1_1965000_2165…" | "nissle1917_NZ…1965000_2165000" (≤ 24 characters) |
| C | V18 large empty top half "detectability constant: no calibration curve" | one line of text, and the panel height shrinks |
| A | V17 (e) "at_2_10", "cpb", "fop", "tx_rate_max_tu" at 100 % missing, no reason; "1 preflight warnings" | human names; "not computed: transcription-rate layer off" next to tx_rate; "1 preflight warning" |

**Edge cases.**
- Fixture K with 4 groups and more than 40 rows: the figure height grows by row count up to a
  maximum, beyond which rows are shown top-down by |pp| and a note says "<n> more rows in the
  data TSV".
- No function (B, C): there is no COG group header.
- In V18 every family can have n < 10 (C). Then the plot is an empty state: "Nothing to show: no
  test family has 10 or more tests in this run."

**AC.**
- 6.1: the V07 figure width is `DOUBLE_IN`. Its data TSV has the three `_pp` columns, and
  `o_minus_e_pp = 100 · o_minus_e / n` (within %.6g).
- 6.2: V07 has no row with `n == 0`, and every group present has a header text artist.
- 6.3: V18 draws no histogram for a family with fewer than 10 p-values, and the note lists it.
  Every y tick of the count axes is an integer.
- 6.4: no V07, V17 or V18 text artist extends beyond the figure bounds. This is checked with
  `get_window_extent` after a draw, the same check M9 uses for clipping.
- 6.5: no V17 text contains a raw column name for which `FEATURE_LABEL` has an entry. "1
  preflight warnings" does not occur.

### PR 7a: V01 ring key and marker binning (P-09)

**Behaviour.**

- **Ring key.** Every ring gets a number at the 0° tick, and a key lists them, for example
  "1 genes (+), 2 genes (−), 3 abundance, 4 peptide density, 5 TUs, 6 promoters and terminators,
  7 GC skew".
- **Colour bar.** The abundance ring gets a colour bar labelled with `abundance_label`.
- **Binning.** Each marker track (promoter, terminator, novel-ORF and ALT markers) with more than
  500 elements is drawn as a density, using `binned(…, max_elements=500)`. The data TSV records
  the bins.
- **Colours.** Terminator ticks are black, `#000000`, not the prophage `#D55E00`.
- **End tick.** The end tick label is dropped when it equals the replicon length rounded to the
  same unit as the 0 tick.

**Examples.**
- K: promoter, terminator, novel-ORF and ALT markers merge into solid rings. After the change,
  each such track is a density ring, since K has 781 novel ORFs.
- A: "200 kb" was drawn on the 0 tick; after the change only "0" is shown.
- The GC skew line is labelled 7 in the key.

**Edge cases.**
- Several contigs: the key is shared, and the numbers are on the first sector only.
- A linear replicon: there is no closing tick and the end label is kept.
- Exactly 500 markers: no binning (strictly more than 500 bins).
- No abundance (C): no abundance ring and no colour bar, and the remaining rings are renumbered.

**AC.**
- 7a.1: the V01 legend or key holds one numbered entry per drawn ring. The count equals the
  number of rings in the data TSV's `track` column.
- 7a.2: on fixture K (any marker track with more than 500 elements), that track's data TSV rows
  are bins, `kind = bin` with `count`, and not elements.
- 7a.3: the terminator colour is `#000000`, and no marker track uses a `style.REGION` colour.
- 7a.4: on a circular single-replicon fixture whose length is a round number of kb, the tick
  labels do not contain the length label twice at angle 0.

### PR 7b: V15 gating and labels (P-07)

**Behaviour.**

- **Gating rule.** Let m = elastic-net `cv_r2_mean` and s = `cv_r2_sd`. The model has **no
  predictive power** when m is NA, m ≤ 0, or m − s ≤ 0 (s NA counts as m ≤ 0 only). Q6 records why
  the rule uses m ± s.
- **No predictive power.** Only panels (a) block ΔR² and (d) observed vs predicted are drawn,
  with a banner "The model does not predict abundance in this run; feature effects are not
  shown."
- **Labels.**
  - Panel titles wrap at the axes width.
  - The partial-dependence labels use `FEATURE_LABEL`, without the column names.
  - "CAI: not available (layer_missing:codon)" becomes "CAI: not computed (codon usage
    unavailable)", placed in the panel, not on a tick.
  - Footers citing "ADR-0031/0032/0033/0034" become plain caveats, for example "Feature blocks
    are correlated; ΔR² is the loss when a block is removed", with no link (figures carry no
    links).
  - The (d) legend moves outside the axes.

**Examples.**

| Run | Before | After |
|---|---|---|
| B | "Elastic net CV R² −0.017 ± 0.066"; full SHAP beeswarm and PD panels drawn | panels (a), (d) and the banner; title "Elastic net CV R² −0.017 ± 0.066: no predictive power" |
| K | 0.133 ± 0.035 (m − s = 0.098 > 0) | all four panels; (d) legend "ribosomal reference (not in R²)" outside the axes |
| snapshot fixture | 0.0359 ± 0.0043 | all four panels (unchanged gating) |

**Edge cases.**
- s NA with m > 0: drawn.
- Fewer than `MIN_ROWS`: the PR 4 empty state applies.
- No abundance (C): V15 is an empty state, and the "not available (no determinants_pd support
  file)" panel texts are gone.

**AC.**
- 7b.1: with `cv_r2_mean = −0.0174037`, `cv_r2_sd = 0.0656274`, V15 draws exactly the panels a
  and d. Its data TSV has no `panel ∈ {b, c}` rows, and the figure holds the banner text.
- 7b.2: with m = 0.05 and s = 0.06 (m − s < 0), the result is the same as 7b.1. With m = 0.0359
  and s = 0.0043, all four panels are drawn.
- 7b.3: no V15 text matches `ADR-\d|layer_missing|promotech_`.
- 7b.4: no V15 text artist extends beyond the figure bounds.

### PR 8: one notice when no function term was assigned (F-01)

**Behaviour.**

- **Coverage source.** Coverage is read from the `function.collections` event. The order of
  sources is: the event, then `category_summary.coverage`, then "not recorded".
- **Every collection at coverage 0** (or `n_terms = 0`):
  - The Function section is replaced by one notice: "No COG, GO or KEGG terms were assigned to
    any gene. With the Bakta light database this is expected: rerun with the full database to
    get function results."
  - The second sentence appears only when `bakta_db_light` is in the run flags. Otherwise it is
    "Check the annotation source."
  - V08 returns `Skip("no_terms:function")`. `SKIP_TEXT` gets the entry "no function terms in
    this run".
  - V14 is an empty state with the same first sentence.
- **Licence notes** (KEGG academic use and retrieval date) print only for collections with
  `n_terms > 0`.
- **Partial coverage:** the section stays, and the notes list coverage per collection, for
  example "COG 72 %, GO 41 %, KEGG not installed".

**Examples (A).**

| Before | After |
|---|---|
| V08 full-width grey violin "none", title "Median abundance 14.2 log2 (99 genes); no COG annotation" | no V08 files; "Function: no function terms in this run" in the run's plot list |
| notes "coverage per collection: not recorded" while the event says `coverage 0.0` | the notice above (coverage known: 0) |
| "KEGG data: academic use only (rest.kegg.jp terms). KEGG retrieval date: 2026-10-15." | not printed (KEGG `n_terms = 0`) |
| empty V14 "No GSEA or detection ORA rows in this run" | V14 empty state: "Nothing to show: no gene has a COG, GO or KEGG term." |

**Edge cases.**
- No KEGG installed (#543): the event has no `kegg` key. KEGG is reported as "not installed"
  and no KEGG note is printed. If COG and GO have coverage 0 too, the notice applies and names
  only COG and GO.
- `function` off (B, C): the section is listed under Not applicable (PR 2), and this PR does
  nothing.
- The event is missing (an old run) and `category_summary` is header-only: coverage is unknown,
  so the section shows "coverage not recorded" and no notice.
- Full DB with real terms: unchanged.

**AC.**
- 8.1: on a fixture whose function event has every coverage 0, the Function section holds
  exactly one `.notice` and no `figure`. `plots/V08_cog_abundance.*` is absent, and V14 is
  `empty=True`.
- 8.2: the notes take coverage from the event when `category_summary` is header-only. The string
  "not recorded" appears only when both sources are absent.
- 8.3: the KEGG licence note is absent when the KEGG `n_terms` is 0 or KEGG is absent from the
  event, and present when it is > 0.
- 8.4: on a no-KEGG fixture with COG coverage > 0, V08 is drawn and the notes say "KEGG not
  installed".

### PR 9: `%.6g` everywhere; protein cross-strain flag (T-01, T-02)

**Behaviour.**
- **T-01.** The writers of `tables/detection_summary.tsv` and `tables/regions_stats.tsv` write
  through `gippa.io.tsv` with `%.6g`, and so does every other published writer found by the AC
  scan.
- **T-02.** `cross_strain_protein` is set on a protein hit only when its `protein_id` is not a
  locus tag of the model **and** its `class ≠ identical`. When at least one protein ID is not a
  locus tag, the run flag `foreign_protein_ids` goes into the manifest `flags` once, whatever
  the class. Its plain text is "proteins come from a related annotation; identical matches are
  not flagged per gene". The gene-level flag in `genes.tsv` follows the hit flag, as today.

**Examples.**

| Table | Before | After |
|---|---|---|
| A `detection_summary.tsv` tier fdr at q 0.001 | `0.0006486894965844569` | `0.000648689` |
| A `detection_summary.tsv` fdr_2pep | `5.870803769695086e-07` | `5.8708e-07` |
| A `regions_stats.tsv` | `42.38390651716451`, `-12.383906517164503` | `42.3839`, `-12.3839` |
| C `detection_summary.tsv` | `4.000000000000001` | `4` |
| B `protein_hits.tsv` (106 rows, all `identity 1`, `query_coverage 1`, `class identical`) | `cross_strain_protein` on all 106; `genes.tsv` flag on all 106 detected genes | no row flag; manifest `flags` gains `foreign_protein_ids` |

**Edge cases.**
- A protein ID that is not a locus tag, with `class = near`: flagged.
- A protein ID that is a locus tag, with `class = divergent`: not flagged (same annotation, so no
  cross-strain claim).
- Integer columns with 7 or more digits (K coordinates such as `1114548`) are not subject to the
  rule.
- `NA`, `inf` and `-inf` are written as `NA`/`inf`/`-inf`, as `gippa.io.tsv` does today.

**AC.**
- 9.1: in an integration run on the sim, no value in a `number`-typed column (per its schema) of
  any published TSV under `tables/` or `mapping/`, or in `plots/data/*.tsv`, has more than 6
  significant digits. Integer and string columns are exempt.
- 9.2: on fixture R, built with foreign IDs and all-identical hits, no protein hit and no gene
  has `cross_strain_protein`, and the manifest `flags` contains `foreign_protein_ids` exactly
  once.
- 9.3: a foreign ID with `class = near` keeps `cross_strain_protein`.
- 9.4: `foreign_protein_ids` is in the §9.4 Run vocabulary, in the manifest schema's flag
  pattern if it has one, and in `FLAG_TEXT`.

### PR 10: out-of-gene start artefacts; I/L rendering; V12 top 15 (V-01, P-17)

**10.0, the re-check on main first (an AC, before any test is written).** Validation runs K on
current main with `sbatch scripts/slurm/gippa_run.sbatch`, using K's inputs and modules
`core,mobile,bgc`, and records the job ID in PROGRESS.md. The first commit of the PR adds a short
table to its PR body. For `tables/variants.tsv`, it gives the counts of rows with `alt_aa = M`
and `ref_codon ∈ {GTG, TTG, CTG, ATT, ATC, ATA}` (the table-11 alternative starts), split by
`confidence` × location (in a gene at `aa_pos = 1` / in a gene, internal / no gene) ×
`start_codon` flag. K at 77df81f5 gives the "before" counts:

| confidence | location | `start_codon` | rows (K, 77df81f5) |
|---|---|---|---|
| high | no gene | no | 9 (e.g. VAR_30 `GTG V→M`, VAR_96, VAR_107 `1114548 (−) TTG L→M`, VAR_129, VAR_135, VAR_185) |
| high | gene, aa 1 | yes | 6 (e.g. VAR_25 `KOIAIB_00334`, VAR_89, VAR_167 `KOIAIB_01749`) |
| high | gene, internal | no | 1 |
| candidate | no gene | no | 117 |
| candidate | gene, aa 1 | yes | 29 |
| candidate | gene, internal | no | 7 |

The expectation on main is that M11-20 credits an initiator M at annotated non-ATG starts. If so,
the "gene, aa 1" rows are 0, and only the out-of-gene rows remain. **The scope of this PR depends
on the result.**

- If the out-of-gene rows are 0 on main, V-01 needs no variants change. The PR keeps only the
  rendering part (P-17, J) and records the result.
- If out-of-gene rows remain, the rule below applies. The statistician signs off, because it
  changes `confidence` values.

The same K run on main gives R-09's report size (PR 13).

**Rule (statistician sign-off with conditions, 2026-11-05; supersedes the first draft).**
AC 10.0 on main: out-of-gene rows remain (high 9, candidate 117; gene aa 1: 0), so the rule applies.
Evidence: K job 1933797 (variants task 1934116; report re-run 1934834; the driver ended FAILED
after the variants step; manifest commit cae1a731, variants code identical to main).

- **Which rows.** An ALT row gets the flag **`putative_start`** (Variants vocabulary) and
  `confidence = candidate`, never `high`, iff all hold:
  1. `alt_type = substitution` and `alt_aa = M`;
  2. `ref_codon` ∈ `annotation.start_codons` minus ATG (the run's configured set, checked against
     `translation_table`, as M11-20 uses it; never a hard-coded list or the NCBI start string);
  3. the codon is not the first codon of an annotated CDS/sORF (those carry `start_codon`);
  4. every alignment at this codon of every ALT peptide of the row has its substitution at peptide
     position 0. "At this codon": same contig and strand, `codons[0] mod L` on circular replicons.
     Positions are in peptide (transcription) order;
  5. the site has no gene, or the reference residue of the preceding codon in transcription order
     is not K/R (trypsin guard). A site with no preceding codon on a linear contig counts as not K/R.
- The rule applies whatever the `call`. `n_ref`, `n_alt`, `call` and the class error rates do not
  change; VCF `CONF` equals the TSV value.
- **id_mixups.** `putative_start` is **not** added to `variants.qc.EXEMPT_FLAGS` (stays
  `{start_codon, stop_codon}`): a same-strain database peptide cannot initiate away from an
  annotated start, so such a row is a mix-up worth counting.
- **Q7** records the rejected options: reusing `start_codon`; excluding these rows; routing them to
  novel-ORF evidence (backlog, after its own chance-mass review).
- **Statistical evidence (C3–C7).**
  - The lead re-scopes `tests/property/test_variant_flags_props.py` and
    `test_flags_adr49_props.py`: M7-05 flags still never change confidence. On the chain, rows
    equal the unflagged chain except `putative_start` rows, where only high→candidate differs.
    The cap never promotes and never changes counts.
  - `class_error_rate` and `class_error_n` are identical per row with the cap on and off.
  - gippasim scenario `alt_start_initiation`, and attack `test_attack_putative_start.py` with
    mutants 1–7. Both are as in the statistician's report; thresholds are in
    `local/handoff/pr10_statistician.md`.
  - The variants calsuite fast tier runs on Slurm before and after. G2, A6a–c and the
    reference-proteome records must be identical. G1 changes are explained by capped rows only:
    pooled |ΔF1| ≤ 0.002, 0 capped truth sites, and G1 stays a strict xfail.
- **Docs.** `docs/methods/variants.md` gets a "Putative starts" note. The `flags.py` and
  `run_variants` docstrings no longer say flags never change confidence. `putative_start` joins
  §9.4 Variants and `FLAG_TEXT`.

**Rendering.**
- V12 (b) writes residue J as "I/L" in figure text ("F→I/L"). The table keeps `J`, as its schema
  says.
- V12 (b) shows the top 15 sites, not 30, and the title adds "(top 15 of <n>; all in
  tables/variants.tsv)".
- "class error unestimated" is stated once in the title when it holds for every shown site, and
  is not repeated per line.

**Examples.**
- K V12 (b) before: 2160 × 3735 px, 30 lines each ending "class error unestimated", "F→J" and
  "M→J". After: 15 lines, "F→I/L", and the title "… ALT sites (top 15 of N; all in tables/variants.tsv)
  (class error not estimated)" (amended at review: clause order kept, 11c wording).
- K VAR_107 (`U00096.3:1114548 (−)`, no gene, `TTG L→M`, high) before. After, if the site
  survives on main: `confidence candidate`, flags `snv;strain_difference;database_derived;putative_start`.

**Edge cases.**
- Minus-strand sites: the "first residue" test is done in peptide order, which is transcription
  order on the minus strand. VAR_107 is minus-strand and is the reference example.
- An origin-crossing codon: the site key uses the unwrapped `codon_start`.
- Genetic codes 4 and 25 (M11-16): the start set comes from the table, not a constant.
- A peptide with M at position 1 but the site at position 2 or later: not affected.
- No variants (A, B, C): V12 is unchanged (empty state or skip).

**AC.**
- 10.0: the PR body holds the K-on-main count table, with the job ID. The PR is not ready for
  review without it.
- 10.1 (conditional on 10.0 rows > 0): on a sim fixture with one out-of-gene GTG site read as
  M by 2 N-terminal peptides on the minus strand, the site has `putative_start` and
  `confidence = candidate`. The same site with the M at peptide position 3 has no
  `putative_start`.
- 10.2: no V12 figure text contains `→J` or `J→`, and "I/L" is used instead.
- 10.3: V12 (b) draws at most 15 site blocks, and the title holds "top 15 of".
- 10.4: when every shown site's class error is NA, the string "class error unestimated" occurs
  exactly once in the figure text.

### PR 11a: `genes.tsv` column groups; run-level flags (T-03, T-06)

**Behaviour.**
- **T-03.** `genes.schema.json` `fields` are reordered into the groups of § 5.2. The writer
  emits columns in schema order, which is already its rule; the test pins it. The column set,
  names and types do not change.
- **T-06.**
  - Columns whose value is the same on every row by construction (`detection_tier`,
    `abundance_type`, `abundance_scale`, `imputation_policy`, `annotation_source`) get
    descriptions starting "Run-level:". They also get the schema property `"gippa:run_level":
    true`, which the reference page renders as "run-level".
  - The gene flag `detectability_constant` is no longer written on rows. It is written once into
    the manifest `flags`. The `detection_summary` row `detectability / overall / all /
    detectability_constant` already exists and stays.
  - §9.4 moves the token from Gene to Run.

**Examples.**
- A: `detected` moves 48 → 15, `detection_status` 50 → 16, `abundance_log2` 61 → 27 and
  `cog_category` 68 → 35.
- C: all 158 rows had `flags = detectability_constant`. After the change they have `flags = NA`,
  and the manifest `flags` contains `detectability_constant`, which it lacked before (C's
  manifest flags: `anti_sd_default, bakta_db_light, c_assumed, cai_reference_missing,
  layer_missing:tx_rates, no_detection_limit_model`).

**Edge cases.**
- K (full genome) has no constant-detectability flag, so nothing moves.
- `gippa compare` reads columns by name, so a 1.0 run and a 1.1 run compare fine (AC 11a.4).
- An all-NA column (11–18 per example run, such as `island_id` and `replichore` on the 200 kb
  fragment) stays, as T-06 says.

**AC.**
- 11a.1: the `genes.tsv` header of an integration run equals the § 5.2 list exactly. The schema
  `fields` order equals it too.
- 11a.2: no row of any `genes.tsv` contains `detectability_constant`. On a constant-detectability
  fixture, the manifest `flags` contains it once.
- 11a.3: the five run-level columns carry `"gippa:run_level": true`, and `tables.md` renders
  "run-level" for them.
- 11a.4: `gippa compare` of a fixture in the 1.0 column order against one in the 1.1 order
  succeeds and gives the same `compare_detection.tsv` as two 1.1 runs.
- 11a.5: the frictionless validation of an integration run passes. The descriptor order is the
  file order.

### PR 11b: one name per quantity (T-04)

**Behaviour.** The renames are in § 5.1. Every writer, reader (viz, report, compare, determinants
support), schema, fixture and golden header follows. Also in this PR:

- `determinants_model.feature_or_block` values `promotech_max`, `promotech_tu` and `mw` become
  `promoter_score`, `promoter_score_tu` and `mw_da`. The `feature` values of the
  `determinants_pd` and `determinants_shap` support files change the same way.
- The V15 PD label "promoter score (promotech_max)" becomes "promoter score". The baseline
  `sim_plots-V15` is regenerated in this PR.
- **Schema-diff exception.**
  - `scripts/release/schema_exceptions.json` maps `{"v1.0.0": {"adr": "ADR-0075", "renamed":
    {"<table>": {"<old>": "<new>", …}}, "retyped": {…}, "reordered": ["genes"]}}`.
  - `schema_diff.check` accepts a removed column without a major bump only if the column is in
    `renamed` for the tag compared against, its new name exists, and `docs/adr/ADR-0075-*.md`
    exists.
  - `retyped` and `reordered` are documentation: the check reads names only.
- `OUTPUT_SCHEMA_VERSION` and `datapackage.json` `output_schema_version` go from `1.2` to `1.3`.
  The major stays 1 under ADR-0075 (Q14).

**Edge cases.**
- `operon_flags` joined to `pairs` on `pair_id`: both tables now have `pair_class` and `flags`.
  A join needs suffixes. The migration note says so (Q2).
- `gippa report` on a 1.0 output tree: the RunView reads with the migration map (Q3; default
  yes). Without it, the M9-01 behaviour 4 `ValueError` names the table and the old column.

**AC.**
- 11b.1: no schema under `src/gippa/schemas/tables/` has a field named `ostir_tir`, `mw`,
  `promotech_max`, `promotech_tu`, `peptide_len`, `class`, `uniop_class`, `flag`,
  `promoter_between` or `terminator_between`. Each new name in § 5.1 exists in its table.
- 11b.2: `test_this_repository_is_schema_compatible_with_its_last_release` passes, unchanged,
  through the exception file.
- 11b.3 (new tests in `test_m12_04_schema_exceptions.py`):
  - an exception naming a missing ADR file fails;
  - an exception listed for another tag does not apply;
  - a removed column not listed still fails.
- 11b.4: every `src/gippa` reference to an old name is gone. This is a grep test over `src/` and
  `workflow/` for the old names as column keys, with an allow-list for the migration map itself.
- 11b.5: if Q3 is accepted, `RunView.from_outdir` on a fixture in the 1.0 layout returns frames
  with the 1.1 names.

### PR 11c: numeric `class_error_rate` and `class_fdr` (T-05)

**Behaviour.**
- `variants.class_error_rate` and `novel_orfs.class_fdr` become `number`. The value is NA
  exactly when it was `unestimated`. Each such row carries the flag: `class_error_unestimated`
  (new, Variants vocabulary) or `class_fdr_unestimated` (exists).
- Readers that parsed the string (the V12 lines, the report cards `_alt_estimate` and
  `_novel_fdr`, the VCF INFO field if it carries it) read the number and use `plain()` for NA.

**Examples.**
- K variants: all 426 rows had `class_error_rate = unestimated`. After the change they are NA,
  with `class_error_unestimated` added. VAR_25 flags become
  `mass_equivalent;snv;strain_difference;database_derived;start_codon;class_error_unestimated`.
- K novel_orfs: 781 rows `unestimated` → NA (every row already has `class_fdr_unestimated`).
- C novel_orfs: `0` → `0`, unchanged.
- The sim golden `sim50k_m7/novel_orfs.tsv` (`0.00646611`, `0.107499`, …): unchanged.

**Edge cases.**
- A class error of exactly 1 (the V12 footnote "not measured if `--reference-peptides` searched
  the reference") stays the number 1, and the footnote logic is unchanged.
- 0 is a number, not NA.
- In the VCF, the INFO value for NA is `.`.

**AC.**
- 11c.1: both columns have schema type `number`, and pandas reads them as float from every
  fixture.
- 11c.2: on every fixture and golden row, the column is NA if and only if the row's `flags`
  contains the matching `_unestimated` flag.
- 11c.3: `class_error_unestimated` is in the variants `flags` pattern, in `FLAG_TEXT` and in the
  §9.4 table.
- 11c.4: the report cards and V12 show "unestimated" in plain words for NA, and no reader
  compares the column with the string `"unestimated"` (grep test).

### PR 12: intermediates under `OUT/intermediate/`, `OUT/README.txt` (O-01)

**Behaviour.**

- **Published top level.** The top level of `OUT/` holds exactly: `README.txt`,
  `run_manifest.json`, `resolved_config.yaml`, `preflight/`, `annotation/`, `mapping/`,
  `tables/`, `plots/`, `report/`, `schemas/`, `logs/`, `intermediate/` and `work/` (Nextflow's
  work directory, as today).
- **Publish rule.**
  - `publishAs` in `workflow/modules/local/gippa_step/main.nf` maps every output whose first path
    component is not one of the published names above to `intermediate/<fn>`.
  - `gippa.runner.local.publish` does the same for `PUBLISHED_FILES` entries outside the
    published folders: `genetable/…`, `context/…`, `novel/…`.
  - The existing `work/` rule (not published) is unchanged.
- **Reads.** `RunView` `PATHS` for those keys point to `intermediate/…`, and `_FALLBACK` adds the
  1.0 locations, so `gippa report` and `gippa compare` still read 1.0 output trees.
  `workflow/tests` `publish_*` tests follow.
- **README.**
  - `OUT/README.txt` is written at the end of `gippa run` and `gippa demo` by `gippa.io`. It is
    plain ASCII text, ≤ 60 lines, with no timestamp and no absolute path.
  - It lists only the entries present, in the fixed order below, each with its fixed sentence.
    It ends with the docs link and the GIPPA version.
- **Docs.** `docs/user-guide/outputs.md` updates the tree and replaces the "Some files sit
  outside this list" paragraph (D-01).

**Example: A's top level.**

Before: 35 files and 15 folders. The files included `annotation_meta.json`, `codon.parquet`,
`features.parquet`, `genes.coords`, `genes.parquet`, `hits_exact.parquet`,
`hits_tolerant.parquet`, `layers.parquet`, `mapping_fdr.parquet`, `novel_candidates.parquet`,
`operon_meta.json`, `ostir_input.csv`, `pairs_raw.parquet`, `placement.parquet`, `rbs.parquet`,
`rrna16s.fa`, `sixframe.npz`, `trna_counts.tsv`, `uniop_cds.gff`, `xrefs.parquet`, … The folders
included `context/`, `genetable/`, `novel/`, `novel_windows/`, `tu_upstream/` and `windows/`.

After:

```text
OUT/
├── README.txt
├── run_manifest.json
├── resolved_config.yaml
├── preflight/  annotation/  mapping/  tables/  plots/  report/  schemas/  logs/
├── intermediate/   codon.parquet, features.parquet, …, genetable/, context/, novel/, windows/, …
└── work/
```

`OUT/README.txt` for A:

```text
GIPPA results

Start here:
  report/report.html    the report: open it in a web browser
  tables/genes.tsv      one row per gene: detection, abundance, function, context, features

Folders:
  tables/        result tables (TSV); every column is described in schemas/datapackage.json
  mapping/       where each peptide or protein was placed in the genome, and what was not placed
  annotation/    the gene model and annotation layers this run used; tools/ holds native tool outputs
  plots/         figures (SVG, PNG, PDF); plots/data/ holds the numbers behind each figure
  report/        report.html, self-contained
  preflight/     the input checks made before the analysis
  schemas/       a machine-readable description of every table (Frictionless Data Package)
  logs/          step logs, events and timings
  intermediate/  files passed between steps; not results; their layout may change between versions
  work/          the workflow engine's working files; needed only to resume this run

Files:
  run_manifest.json     versions, databases, parameters, input checksums and seeds
  resolved_config.yaml  the configuration this run used

Coordinates: TSV, GFF3 and VCF are 1-based closed; BED and proBed are 0-based half-open.
Table reference: https://khoa-yelo.github.io/gippa/reference/tables/
Written by GIPPA 1.1.0.
```

**Edge cases.**
- `--resume` into a 1.0 output directory: the old top-level intermediates are left in place
  (GIPPA never deletes them). One warning event, `outputs.legacy_intermediates`, gives their
  count and says they can be removed.
- The local developer executor: same layout.
- A failed run: `README.txt` is still written, listing what exists. `report/` may be absent.
- `gippa demo` output: same layout.
- `annotation/tools/` is unchanged.

**AC.**
- 12.1: after an integration (stub) run, the set of top-level entries of `OUT/` is a subset of
  the published list above, and every other file of the run is under `intermediate/`.
- 12.2: every path a user-facing doc names under `tables/`, `mapping/`, `annotation/`, `plots/`
  and `report/` exists at the same path as in 1.0. A test lists the published files of a stub run
  on main@1.0 layout vs this branch, excluding `intermediate/` and the moved files.
- 12.3: `README.txt` names every top-level entry that exists, except itself, and no entry that
  does not exist. It has no timestamp, no `/home`, no `farmshare`, and it is ASCII only.
- 12.4: `RunView.from_outdir` on a 1.0-layout fixture (intermediates at the top level) finds
  `gene_table`, `gene_context` and `island_diagnostics` through the fallbacks.
  `gippa report` on it succeeds.
- 12.5: `publishAs` in the nf-test maps `codon.parquet` → `intermediate/codon.parquet`,
  `genetable/genes.parquet` → `intermediate/genetable/genes.parquet`, and `tables/genes.tsv` →
  `tables/genes.tsv`.

### PR 13: detection views; shareable reports; report size (T-07, R-09)

**Behaviour.**

- **Detection views (T-07).** The detection section shows two pivoted tables built from
  `detection_summary.tsv`:
  - **Tier:** rows are `q_threshold`, and the columns are `n_detected`, `n_genes`, `fdr`,
    `fdr_2pep` and `fp`.
  - **Detectability:** O/E by group.

  The long file itself is unchanged, and the raw 20-row preview of it is removed. No new table
  file is added (Q8).
- **Shareable reports (R-09).** `gippa report --redact` rewrites the provenance command line:
  - every absolute path becomes `<input>/<basename>` (`<outdir>` for the output directory);
  - the `--profile` value becomes `<profile>`.

  Projects that publish a report use it, for example docs pages. The default report is unchanged.
- **Report size.**
  - The report size is measured on the K-on-main run of AC 10.0.
  - If `report.html` is ≥ 20 MB (80 % of the 25 MB §10.3 target), V13 locus SVGs are linked
    rather than inlined, through `plots/V13_novel_orfs_locus_*.svg` relative links.

**Examples.**
- A tier view:

  | q | detected | genes | FDR | FDR (2-pep) | FP |
  |---|---|---|---|---|---|
  | 0.001 | 94 | 94 | 0.00065 | 5.9e-07 | 0.061 |
  | 0.005 | 98 | 99 | 0.0061 | 6.1e-05 | 0.60 |
  | 0.01 | … | … | 0.012 | 0.00027 | … |

- A provenance before: `gippa run --profile <site profile>,apptainer --genome /…/genome/nissle_r200k.fna.gz …`.
  After `--redact`: `gippa run --profile <profile> --genome <input>/nissle_r200k.fna.gz …`.
- K report: 19.8 MB at 77df81f5; re-measure on main.

**Edge cases.**
- Tier `protein_list` or `two_peptide` has no FDR rows, so the tier view shows the counts only.
- `detection_summary` absent: the views are omitted.
- `--redact` on a manifest without `command_line`: nothing to rewrite.

**AC.**
- 13.1: on fixture P, the detection section has a table with one row per distinct
  `q_threshold` of the `tier` section. The values equal the long file's values through
  `fmt_num`.
- 13.2: the `report --redact` output contains no `/` followed by a path component of the
  fixture's absolute input paths, and no profile name from the manifest.
- 13.3: with a fixture whose inlined V13 SVGs push the report to ≥ 20 MB (synthetic large
  SVGs), V13 locus views are `<a>`/`<img>` with relative `src` to `plots/`. With small SVGs, they
  stay inlined. The size threshold is a named constant.

---

## 4. Dependency order and parallelism

**The column PRs land first.**
- They change table goldens and snapshot-fixture TSVs, not plot images, apart from the V15 label
  in 11b.
- Several plot and report PRs read the renamed or retyped columns: V11 reads `pairs.class`, V15
  reads `promotech_max`, and V12 and the cards read `class_error_rate`/`class_fdr`.
- If the column PRs landed last, each of those readers would be edited twice and the plot
  fixtures regenerated after their baselines.
- Landing them first means every later PR is written against the final names, and Khoa approves
  each plot baseline set against final data.

```text
Wave 1 (parallel, no dependencies):
  ADR-0075 (lead, docs only)        PR 0 (helpers)        PR 9 (%.6g, cross-strain)        PR 12 (layout)
Wave 2 (after ADR-0075 is merged as Proposed or Accepted):
  PR 11a (genes order, run-level)   PR 11b (renames)      PR 1 (run box, cards, TOC; after PR 0)
Wave 3:
  PR 11c (after 11b)                PR 5 (style; after 0 and 11b)        PR 2 (after 1)
Wave 4 (parallel; each regenerates only its own baselines):
  PR 4 (after 5)   PR 6 (after 5)   PR 7a (after 5)   PR 7b (after 5)
  PR 8 (after 2 and 5)              PR 3 (after 2)    PR 13 (after 2)
  PR 10 (after AC 10.0 run, 11c and 5)
```

- **Ordering rules.**
  - PR 5 goes first in the plot wave: it touches nearly every `sim_plots` baseline (titles, q
    format, abundance axes, locus fill). PRs 4, 6, 7a, 7b, 8 and 10 then regenerate only their
    own plots.
  - PRs 1 → 2 → {3, 8, 13} are serial because they edit the same template and
    `report_structure.json`.
  - PR 11b → 11c is serial (both edit `novel_orfs.schema.json`). PR 11a is independent of 11b.
- **Shared goldens.**
  - PR 9 and PR 11a both change `sim50k/m4` and `sim100k_real/m4` goldens: PR 9 changes the
    `detection_summary`/`regions_stats` values, and PR 11a the `genes.tsv` order. These are
    different files, so there is no conflict.
  - `genes.tsv` is touched only by 11a.
- **The K re-check (AC 10.0)** goes to Slurm in wave 1, so PR 10 and PR 13's size decision are
  not blocked later.
- **Parallel cap.** At most 6 PRs are in flight at once: wave 4 has 8 PRs, so the lead queues 2.
- **Rebase policy.** Branches merge `origin/main` and push plainly (no rebase). A merge that
  conflicts on a baseline PNG is resolved by regenerating that baseline, which goes back to Khoa
  for approval.

---

## 5. Migration note (users): 1.0 → 1.1

This text goes into `changelog.d/<PR>.removal.md` of PRs 11a, 11b, 11c and 12, each with its
part, and into `docs/reference/tables.md` as a hand-written page section. That section is
outside the generated block: `gen_reference.py` gets a `MIGRATION` include, or the section goes
into `docs/user-guide/outputs.md` › "Changes in 1.1" with a link from `tables.md`. The
implementer picks whichever keeps `gen_reference.py --check` passing. User-facing wording:

> **GIPPA 1.1 changes some table columns.** This is a one-time exception to the rule that column
> changes wait for a major version. Scripts that select columns by name keep working, except for
> the renamed columns below. Scripts that select columns by position must be updated for
> `genes.tsv`. Old names are not kept as aliases.

### 5.1 Renamed columns

| Table | 1.0 name | 1.1 name | Why |
|---|---|---|---|
| `tables/novel_orfs.tsv` | `ostir_tir` | `tir_ostir` | same name as in `genes.tsv` and `determinants_features.tsv` |
| `tables/determinants_features.tsv` | `mw` | `mw_da` | unit in the name, same as `genes.tsv` |
| `tables/determinants_features.tsv` | `promotech_max` | `promoter_score` | the same value as `genes.tsv › promoter_score` |
| `tables/determinants_features.tsv` | `promotech_tu` | `promoter_score_tu` | the same value as `genes.tsv › promoter_score_tu` |
| `mapping/peptide_hits.tsv` | `peptide_len` | `peptide_len_aa` | unit in the name |
| `tables/pairs.tsv` | `class` | `pair_class` | one name in `pairs.tsv` and `operon_flags.tsv` (Q1) |
| `tables/operon_flags.tsv` | `uniop_class` | `pair_class` | as above; the class is not from UniOp alone (distance-model fallback, SW decoys) |
| `tables/pairs.tsv` | `promoter_between` | `promoter_score_between` | it is a score; `operon_flags.promoter_between_yes` is the yes/no call |
| `tables/pairs.tsv` | `terminator_between` | `terminator_score_between` | same reason (Q15) |
| `tables/operon_flags.tsv` | `flag` | `flags` | every table's flag column is `flags` (Q2) |

Values that name a feature change the same way: in `tables/determinants_model.tsv ›
feature_or_block`, `promotech_max` → `promoter_score`, `promotech_tu` → `promoter_score_tu` and
`mw` → `mw_da`. The same holds for the `feature` column of `plots/support/determinants_pd.tsv`
and `determinants_shap.tsv`.

Unchanged on purpose: `operon_flags.promoter_between_yes` and `terminator_between_yes` (yes/no
calls) and `pairs.tx_rate_between`.

### 5.2 `tables/genes.tsv` column order

The same 73 columns, in groups: identity (1–14) → detection (15–26) → abundance (27–33) →
function (34–39) → context (40–49) → regulatory and sequence features (50–72) → flags (73).

| Group | New position | Column | 1.0 position |
|---|---|---|---|
| identity | 1 | `locus_tag` | 1 |
| identity | 2 | `source_id` | 2 |
| identity | 3 | `gene` | 3 |
| identity | 4 | `product` | 4 |
| identity | 5 | `feature_type` | 5 |
| identity | 6 | `hypothetical` | 6 |
| identity | 7 | `contig` | 7 |
| identity | 8 | `start` | 8 |
| identity | 9 | `end` | 9 |
| identity | 10 | `strand` | 10 |
| identity | 11 | `length_aa` | 11 |
| identity | 12 | `partial` | 12 |
| identity | 13 | `crosses_origin` | 13 |
| identity | 14 | `replicon_type` | 14 |
| detection | 15 | `detected` | 48 |
| detection | 16 | `detection_status` | 50 |
| detection | 17 | `detection_tier` | 49 |
| detection | 18 | `gene_q_value` | 51 |
| detection | 19 | `detected_2pep` | 52 |
| detection | 20 | `n_peptides_unique` | 53 |
| detection | 21 | `n_peptides_multi_identical` | 54 |
| detection | 22 | `n_peptides_shared` | 55 |
| detection | 23 | `n_peptides_tolerant` | 56 |
| detection | 24 | `coverage_frac` | 57 |
| detection | 25 | `best_peptide_q` | 58 |
| detection | 26 | `p_detect` | 59 |
| abundance | 27 | `abundance_log2` | 61 |
| abundance | 28 | `abundance_raw` | 60 |
| abundance | 29 | `top3_log2` | 66 |
| abundance | 30 | `abundance_imputed` | 62 |
| abundance | 31 | `abundance_type` | 64 |
| abundance | 32 | `abundance_scale` | 65 |
| abundance | 33 | `imputation_policy` | 63 |
| function | 34 | `cog_id` | 67 |
| function | 35 | `cog_category` | 68 |
| function | 36 | `go_terms` | 69 |
| function | 37 | `kegg_ko` | 70 |
| function | 38 | `ec` | 71 |
| function | 39 | `annotation_source` | 72 |
| context | 40 | `tu_id` | 15 |
| context | 41 | `tu_position` | 16 |
| context | 42 | `tu_size` | 17 |
| context | 43 | `uniop_p_prev` | 18 |
| context | 44 | `uniop_p_next` | 19 |
| context | 45 | `region_ids` | 20 |
| context | 46 | `island_id` | 21 |
| context | 47 | `replichore` | 22 |
| context | 48 | `replichore_position` | 23 |
| context | 49 | `gc_local_dev` | 24 |
| features | 50 | `promoter_score` | 25 |
| features | 51 | `promoter_score_tu` | 26 |
| features | 52 | `tx_rate_max_tu` | 27 |
| features | 53 | `tx_rate_sum_tu` | 28 |
| features | 54 | `terminator_conf_downstream` | 29 |
| features | 55 | `tir_ostir` | 30 |
| features | 56 | `sd_dg` | 31 |
| features | 57 | `mrna_dg` | 32 |
| features | 58 | `spacer_bp` | 33 |
| features | 59 | `start_codon` | 34 |
| features | 60 | `cai` | 35 |
| features | 61 | `tai` | 36 |
| features | 62 | `enc` | 37 |
| features | 63 | `fop` | 38 |
| features | 64 | `gc3` | 39 |
| features | 65 | `cpb` | 40 |
| features | 66 | `at_2_10` | 41 |
| features | 67 | `rare_codons_2_10` | 42 |
| features | 68 | `mw_da` | 43 |
| features | 69 | `pi` | 44 |
| features | 70 | `gravy` | 45 |
| features | 71 | `tm_segments` | 46 |
| features | 72 | `n_observable_peptides` | 47 |
| flags | 73 | `flags` | 73 |

### 5.3 Type changes

| Table | Column | 1.0 | 1.1 |
|---|---|---|---|
| `tables/variants.tsv` | `class_error_rate` | number or the word `unestimated` | number; `NA` when unestimated, and the row has the flag `class_error_unestimated` |
| `tables/novel_orfs.tsv` | `class_fdr` | number or the word `unestimated` | number; `NA` when unestimated, and the row has the flag `class_fdr_unestimated` |

### 5.4 Flag changes

| Flag | 1.0 | 1.1 |
|---|---|---|
| `detectability_constant` | on every row of `genes.tsv` | once in `run_manifest.json › flags` (and the `detection_summary.tsv` overall row, as before) |
| `cross_strain_protein` | on every protein hit whose ID is not a locus tag of the genome | only when such a hit is also not identical (`class` near or divergent); the run flag `foreign_protein_ids` records related-annotation IDs once |
| `class_error_unestimated` | – | new (variants) |
| `putative_start` | – | new (variants), if PR 10 needs it |
| `foreign_protein_ids` | – | new (run) |

### 5.5 Output folder

Step intermediates that sat at the top of the output folder (`codon.parquet`, `features.parquet`,
`genetable/`, `context/`, `novel/`, …) are now under `intermediate/`. Published results (`tables/`,
`mapping/`, `annotation/`, `plots/`, `report/`, `schemas/`, `logs/`, `preflight/`,
`run_manifest.json`, `resolved_config.yaml`) keep their paths. A new `README.txt` explains each
folder. `gippa report` and `gippa compare` still read 1.0 output folders.

---

## 6. Questions for the lead

Every question has a reversible default, which the PRs use unless the lead says otherwise.

| # | Question | Default |
|---|---|---|
| Q1 | The pair-class name. The review proposed `uniop_class` for both tables. This spec proposes `pair_class`, because the class also comes from the distance-model fallback (ADR-0061) and includes SW decoys and CV/DV by strand. `class` is also awkward as an attribute name in pandas. | `pair_class` |
| Q2 | `operon_flags.flag` → `flags` (the §9.2 rule: every flag column is `flags`, and its values are §9.4 Operons tokens, one per row) or `operon_flag` (no collision when joined with `pairs.flags` on `pair_id`). | `flags` |
| Q3 | Should `RunView` rename 1.0 columns on read (one migration map), so that `gippa report`/`compare` work on 1.0 outputs? This is not an output alias. Without it, `gippa report` on a 1.0 run fails on `pairs.tsv` with a clear `ValueError`. | yes, read-side map |
| Q4 | The schema-diff exception mechanism (`scripts/release/schema_exceptions.json` keyed by tag + ADR id) has to be named in ADR-0075. Does the lead accept that file and its location? | as specified (PR 11b) |
| Q5 | V13 also writes a text-only index image (`V13_novel_orfs_locus.png`, baseline `sim_plots-V13_novel_orfs_locus.png`). Khoa's decision covered V03 only. Should V13's go the same way? | keep V13's (outside the decision) |
| Q6 | P-07 "interval includes 0": the CV R² interval is never estimated (`ci_unestimated` in B, K and the fixture). The rule uses the mean − SD that the title shows. An alternative is a t-interval over CV repeats, which would be a statistic (viz is in `NON_STATS_PATHS`, so the m ± s display rule needs no stats review). | m ≤ 0 or m − s ≤ 0 |
| Q7 | V-01, if out-of-gene sites remain on main: a new flag `putative_start` with confidence capped at `candidate` (minor), reusing `start_codon` (a reinterpreted token, so breaking), or routing the sites to novel-ORF start evidence (scope growth). Needs statistician sign-off. | `putative_start` + candidate |
| Q8 | T-07: add a wide `tables/detection_overview.tsv` (a new table, minor)? | no new file; report views only |
| Q9 | R-09: the redaction is a new `gippa report --redact` option (a minor CLI change) rather than a step in the docs or example build scripts. No gallery exists yet. | `--redact` option |
| Q10 | `datapackage.json` has resources for only 12 of the ~30 tables, while §9.5 says one per table. Adding the rest widens frictionless validation to every table, and may surface schema mismatches. | generator fallback only in 1.1; backlog item for the resources |
| Q11 | `--parquet` / `output.parquet` is accepted but no writer reads it, so there is no Parquet twin. Implement it, or remove the option (removal is breaking)? | backlog; not in 1.1 |
| Q12 | Methods links in the report: the review said "relative link", but the report is offline. This spec uses absolute docs-site `<a href>` links (not loaded at view time). Confirm that the existing "no external URL" test only forbids loads. | absolute docs links |
| Q13 | Baseline approvals. PR 5 regenerates most `sim_plots` baselines, and PRs 4/6/7a/7b/8/10 then regenerate subsets again. Khoa could instead approve once on a stacked viz branch at the end. | approval per PR (decision 4) |
| Q14 | `output_schema_version` 1.2 → 1.3 with the major kept at 1 under ADR-0075. `gippa compare` then treats 1.0 and 1.1 runs as the same major and does not refuse them. | 1.3 |
| Q15 | `pairs.terminator_between` → `terminator_score_between`, and `determinants_features.promotech_max/promotech_tu` → `promoter_score/promoter_score_tu`, go slightly beyond Khoa's explicit list. They follow the same rule, and the review lists `promotech_*` under T-04. | include |
| Q16 | ADR-0075 also has to record the §10.2 changes (V3 index dropped, V12 top 15, V15 gating) and the §10.1 locus fill. Should one ADR cover both the table exception and the plot changes, or two ADRs? | one ADR-0075 |
