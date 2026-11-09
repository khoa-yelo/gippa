# M12-05b: host and flanking genes of novel ORFs

Status: draft (spec-writer). Roles: spec-writer (this file), test-author (PR A), implementer
(PR B), reviewer, lead (routing, questions to Khoa).

## Source

Khoa's review of the 1.1.0rc1 outputs (2026-11-09): "novel orfs should contain gene that the orfs
land on and/or flanked by." Evidence: route E,
`local/example/routes/1945698/gippa-example/runs/E/tables/novel_orfs.tsv` (8 rows: 1
`nterm_extension`, 5 `antisense`, 1 `cds_outframe`, 1 `intergenic`). `overlapping_gene` holds a
bare locus tag (`MPHIBC_0020`) and is `NA` for the intergenic row, so the row says nothing about
where the ORF sits. No row gives a gene name, and no row gives a neighbour.

## Spec sections touched

- GIPPA_SPEC.md §9.3 `tables/novel_orfs.tsv`: eight columns added after `overlapping_gene`.
  `flags` stays last. Adding columns is a minor change (§9.2 "Versioning"). ADR-0075 already
  allows column changes in 1.1, so this spec needs no new ADR.
- §6.5: internal 0-based half-open coordinates. The new columns hold no coordinates, only locus
  tags, labels and distances, so the writer converts nothing new.
- §9.2: `NA` stays the missing-value token in the TSV (see rule 6).
- Related, referenced only: **M12-05a** (V13 ordering) and **M12-05c** (captions). The parallel
  viz investigation of the literal "NA" gene names shares rule 5.

## What already exists (reuse)

| Need | Existing code |
|---|---|
| host gene | `gippa.novel.candidates._overlapping_gene` gives `Orf.overlapping_gene`. It uses the C11 `gene` for extensions, pseudogenes and sORFs. For `cds_outframe` and `antisense` it takes the smallest same-strand or opposite-strand overlapping locus tag. It is `None` for `intergenic`. |
| ORF interval | `gippa.novel.outputs._span(r, start, n_aa)`, also used to rebuild the GFF3 intervals (origin-aware, `end > length`) |
| overlap query | `GeneModel.overlapping(iv)` (origin-aware, either strand, model order) |
| gene name | `Gene.gene` (`str \| None`) and `genes.tsv` `gene` |
| reading old trees | `gippa.io.migration` and `RunView.read_schema_tsv`. The reader raises when a schema column is missing, so older trees need a fill (rule 8). |
| report preview | `gippa.report.build._table` / `GENES_PREVIEW` (a curated column subset for genes) |

## Rules

1. **`overlapping_gene` is kept, unchanged in meaning and format.** It is the host locus tag, and
   `evidence.py`, `fdr.py` and `patch_features` read it. This spec does not duplicate or rename
   it. It only adds columns next to it.

2. **New columns**, in this order, right after `overlapping_gene`:

   | Column | Type | Meaning |
   |---|---|---|
   | `overlapping_gene_label` | string | label (rule 5) of `overlapping_gene`; `NA` iff `overlapping_gene` is `NA` |
   | `overlap_relation` | enum{`same_strand`, `antisense`, `none`} | strand of the host gene relative to the ORF; `none` when there is no host gene |
   | `other_overlapping_genes` | list (`;`) | locus tags of every other gene that overlaps the ORF span (rule 3), on either strand, in gene-model order; `NA` when none |
   | `upstream_gene` | string | locus tag of the nearest gene on the ORF's 5′ side (rule 4) |
   | `upstream_gene_label` | string | its label (rule 5) |
   | `upstream_distance_nt` | integer ≥ 0 | bases strictly between the ORF span and that gene |
   | `downstream_gene` | string | locus tag of the nearest gene on the ORF's 3′ side |
   | `downstream_gene_label` | string | its label |
   | `downstream_distance_nt` | integer ≥ 0 | as upstream |

   `overlap_relation` for each class: `cds_outframe`, `nterm_extension`, `cterm_extension`,
   `pseudogene` and `sorf` give `same_strand`, and `antisense` gives `antisense`. It is computed
   from the strands, not from the class, so a disagreement between the two shows up and does not
   hide.

3. **ORF span.** The span is the longest interval, from `start_longest` through the stop codon.
   Without a start codon it starts at the first peptide's codon, which is the interval that
   `orf_id` names: `_span(r, start_longest, length_aa_longest)`. Origin-crossing spans keep
   `end > length`.

4. **Flanking genes.** Upstream and downstream are relative to the ORF's strand, in transcription
   direction. On `+`, upstream means lower coordinates and downstream means higher. On `-` they
   swap. The gene's own strand does not matter.
   - A candidate is any gene that does not overlap the span. Genes that overlap the span are
     reported under rule 2 and are never flanks.
   - Walk from the span edge in that direction. The flank is the gene whose near edge comes
     first. The distance is the gap: on `+` upstream, `span.start − gene.end`.
   - **Linear replicon:** the walk stops at the replicon end, and with no candidate before it the
     three columns are `NA`.
   - **Circular replicon:** the walk wraps through the origin, at most once around, and the
     distance includes the wrap. On a circular replicon whose only non-overlapping gene is G, G is
     both the upstream and the downstream gene, with different distances.
   - The walk never crosses into another replicon. Plasmids follow the same rules, and their
     topology comes from `Replicon.topology`.
   - **Ties** (equal distance): the smallest locus tag wins. This is the same rule M3-14 uses for
     `overlapping_gene`.

5. **Label** = `Gene.gene` when it is a usable name, otherwise the locus tag. A name is unusable
   when it is `None`, empty or whitespace, or the literal `NA` (any case). A label is never the
   literal `NA` when a gene exists. Put this in one pure helper,
   `gippa.model.genes.gene_label(gene) -> str`. If the viz investigation lands an equivalent
   helper first, reuse that one and do not add a second.

6. **What replaces `NA` when there is no gene.**
   - *Table:* the TSV keeps `NA` for `overlapping_gene`, its label, `other_overlapping_genes` and
     a missing flank (§9.2, unchanged). The meaning is carried by a non-missing cell: for "no host
     gene" that cell is `overlap_relation = none`. This is a known state, not missing data, so it
     is an enum value and not `NA`. A missing flank reads unambiguously from its `NA` distance: it
     can only mean the walk reached a linear end, or the replicon has no other gene.
   - *Report preview:* rule 9.

7. **Gene set.** The genes are the gene model's `Gene` objects (CDS, sORF, pseudogene), which are
   exactly the rows of `genes.tsv`, so every tag joins to `genes.tsv.locus_tag`. Non-coding
   features (tRNA, rRNA, CRISPR) and other novel ORFs are not genes here (Q1). Partial genes at a
   linear end count, at their stored interval.

8. **Reading older trees.** A 1.0 or 1.1.0rc1 `novel_orfs.tsv` lacks the eight columns.
   `gippa.io.migration` gains `ADDED_COLUMNS = {"novel_orfs": (...)}`, and
   `RunView.read_schema_tsv` fills every absent added column with `NA`, `overlap_relation`
   included, because the relation was not computed. A tree that lacks some other schema column
   still raises. `read_novel_orfs`, which is internal and same-version, keeps its exact-header
   check.

9. **Report preview and V13.**
   - The Novel ORFs section previews a curated subset, `NOVEL_PREVIEW`, in the way `GENES_PREVIEW`
     does: `orf_id`, `class`, `called`, `class_fdr`, `strand`, `length_aa_best`,
     `overlapping_gene_label`, `overlap_relation`, `upstream_gene_label`, `upstream_distance_nt`,
     `downstream_gene_label`, `downstream_distance_nt`, `n_peptides`, `flags`. It uses only the
     columns present and keeps the rows in table order.
   - In the preview, an `NA` label cell shows `none` for `overlapping_gene_label` and `—` for a
     flank label. Every other cell is shown as written.
   - V13 is unchanged by this spec. M12-05a owns V13's ordering, and M12-05c owns the captions. A
     caption may cite these columns once both land.

10. **Writer.** The new cells are strings or integers, so `write_novel_orfs` adds nothing to
    `POSITION_COLUMNS` and the distances are not shifted. The new cells are computed by a pure
    function, `gippa.novel.gene_context.gene_context(table, gm) -> pd.DataFrame`. It returns the
    eight columns, indexed like `table`, and `steps/novel.py` calls it after `orf_table`.
    `orf_table` keeps its signature.

## Worked examples

The examples use a synthetic contig `c1` of length 10 000, circular unless stated. The table
gives genes as internal `[start, end)` intervals, with the TSV `start–end` in brackets:

| locus_tag | gene | strand | internal | TSV |
|---|---|---|---|---|
| T_0001 | dnaA | + | [100, 400) | 101–400 |
| T_0002 | NA (literal) | − | [500, 1100) | 501–1100 |
| T_0003 | lacZ | + | [1300, 2500) | 1301–2500 |
| T_0004 | None | + | [9500, 9900) | 9501–9900 |

In the outputs, columns are written as `ovl / label / relation / other ; up / label / dist ;
down / label / dist`, and `–` means `NA`.

| # | ORF (class, strand, internal span; TSV start_longest, end) | Expected |
|---|---|---|
| E1 | intergenic, +, [2600, 2900); 2601, 2900 | – / – / none / – ; T_0003 / lacZ / 100 ; T_0004 / T_0004 / 6600 |
| E2 | intergenic, −, [2600, 2900); 2900, 2601 (strand swaps the sides) | – / – / none / – ; T_0004 / T_0004 / 6600 ; T_0003 / lacZ / 100 |
| E3 | antisense, −, [1400, 1700) inside T_0003 | T_0003 / lacZ / antisense / – ; T_0004 / T_0004 / 7800 ; T_0002 / T_0002 / 300 |
| E4 | cds_outframe, +, [1000, 1450), overlaps T_0002 (−) and T_0003 (+) | T_0003 / lacZ / same_strand / T_0002 ; T_0001 / dnaA / 600 ; T_0004 / T_0004 / 8050 |
| E5 | nterm_extension of T_0002, −, [500, 1250); 1250, 501 | T_0002 / T_0002 / same_strand / – ; T_0003 / lacZ / 50 ; T_0001 / dnaA / 100 |
| E6 | intergenic, +, [9950, 10060) across the origin; 9951, 60 | – / – / none / – ; T_0004 / T_0004 / 50 ; T_0001 / dnaA / 40 |
| E7 | `c1` **linear**; intergenic, +, [9950, 9990) | – / – / none / – ; T_0004 / T_0004 / 50 ; – / – / – (with `c1` circular, downstream is T_0001 / dnaA / 110) |
| E8 | plasmid `p1` (circular, 3000) with one gene P_0001 `repA` +, [200, 1200); intergenic, +, [1500, 1800) | – / – / none / – ; P_0001 / repA / 300 ; P_0001 / repA / 1400 (linear `p1`: downstream – / – / –) |
| E9 | contig `c2` with no genes; any ORF | – / – / none / – ; – / – / – ; – / – / – |
| E10 | the four genes plus T_0005 +, [5000, 5400) and T_0006 −, [4800, 5400); intergenic, +, [5500, 5700) | upstream tie at gap 100: T_0005 / T_0005 / 100 (smallest tag); downstream T_0004 / T_0004 / 3800 |

Checks for E1: lacZ ends at TSV 2500 and the ORF starts at 2601, so bases 2501–2600 lie between
them, which is 100. E6 shows that a distance never counts the origin itself: the ORF ends at
internal 60 and T_0001 starts at 100, so the distance is 40.

## Edge cases

- **Empty table** (no ORF, `sim_empty`): the header has all 34 columns and there are no rows.
- **Replicon with no genes, or with only the host gene** on a linear replicon: the flanks are
  `NA`. On a circular replicon, the other genes are used as in E8.
- **ORF overlapping several genes:** the host is still `overlapping_gene`, and every other
  overlapping gene goes in `other_overlapping_genes` (E4). A gene that overlaps the span by one
  base counts.
- **No start codon** (`start_longest` NA): the span is taken from the first peptide (rule 3), as
  `orf_id` and the GFF3 already do.
- **Partial ORF open at a linear end:** the span ends at the contig end, and the flank on that
  side is `NA`.
- **Origin-crossing gene or ORF:** the overlap and the gaps are computed on the split
  intervals; E6 is the case to test.
- **Minus strand:** the sides swap (E2, E3, E5). The distances are the same numbers.
- **Ties:** the smallest locus tag, for the flanks (E10) and for the host (existing rule).
- **Literal `NA` or blank gene name:** the label is the locus tag (T_0002, T_0004).

## Acceptance criteria (test-checkable)

- AC1: `gene_context` on the synthetic model returns E1–E10 exactly.
- AC2: `gene_label` returns `dnaA` for a named gene, the locus tag for `None`, `""`, `"  "`,
  `"NA"` and `"na"`, and never returns `"NA"`.
- AC3: `NOVEL_ORFS_COLUMNS` equals the 1.1.0rc1 list with the eight columns inserted after
  `overlapping_gene`. `novel_orfs.schema.json` lists them with the types in rule 2, the enum for
  `overlap_relation` and `minimum: 0` for the distances.
  `test_this_repository_is_schema_compatible_with_its_last_release` passes.
- AC4: `write_novel_orfs` → `read_novel_orfs` round-trips a frame holding every example, and the
  distances are not shifted by 1.
- AC5: on every row, `overlapping_gene_label` is NA if and only if `overlapping_gene` is NA, and
  `overlap_relation == "none"` if and only if `overlapping_gene` is NA (on written runs).
- AC6: every non-NA tag in `overlapping_gene`, `other_overlapping_genes`, `upstream_gene` and
  `downstream_gene` is a `genes.tsv` `locus_tag` of the same run (sim golden run).
- AC7: `RunView` reads `novel_orfs.tsv` with the 1.1.0rc1 header (route E copy, or a fixture
  without the eight columns) and the 1.0 header. The eight columns come back NA. A file missing
  `class` still raises.
- AC8: the report's Novel ORFs preview header equals `NOVEL_PREVIEW`, restricted to the columns
  present. NA labels render as `none` and `—` (rule 9).
- AC9: `scripts/docs/gen_reference.py` regenerates `docs/reference/tables.md` with the eight rows,
  and the reference-sources test passes.

## Changes beyond the table (PR B)

- `src/gippa/schemas/tables/novel_orfs.schema.json`: eight fields, with descriptions that state
  the strand convention. `datapackage.json` `output_schema_version` stays `1.3`, because 1.3 is
  still unreleased and already above v1.0.0, which is what the schema-diff check needs (Q3).
- `src/gippa/io/migration.py`: `ADDED_COLUMNS` and its docstring line. `run_view.py` does the
  fill.
- `docs/user-guide/outputs.md` › Changes in 1.1: one "Added columns" row for `novel_orfs`.
  `docs/reference/tables.md` is regenerated. `docs/methods/novel-orfs.md`: one paragraph with
  rule 4.
- `changelog.d/<PR>.added.md`: "novel_orfs.tsv names the host gene and the nearest genes on each
  side, with distances."
- Goldens and fixtures: `tests/data/golden/sim50k_m7/novel_orfs.tsv` and
  `tests/data/snapshot/*/OUT/tables/novel_orfs.tsv` gain the columns. They are regenerated, not
  hand-edited.

## Work split

| PR | Title | Role | Contents | ~lines |
|---|---|---|---|---|
| A | `test(novel): host and flanking genes (M12-05b)` | test-author | `tests/unit/test_m12_05b_gene_context.py` (AC1, AC2, AC4, AC5 on the synthetic model of the examples); `tests/unit/test_m12_05b_outputs.py` (AC3, AC7, AC8); xfail-strict until PR B | 300 |
| B | `feat(tables): host and flanking genes in novel_orfs (M12-05b)` | implementer | `novel/gene_context.py`, `model/genes.gene_label`, `outputs.py` columns, `steps/novel.py` call, schema, migration and `run_view` fill, `NOVEL_PREVIEW`, docs, changelog; goldens and fixtures regenerated (data diffs are not counted in the line budget); AC6 and AC9 run on Slurm | 350 + data |

B depends on A. Neither touches V13 code, so the two run in parallel with M12-05a and M12-05c. If
the viz label helper merges first, B imports it.

## Questions for the lead (and Khoa), with defaults

| # | Question | Default |
|---|---|---|
| Q1 | Should tRNA, rRNA and other non-coding features count as overlapping or flanking genes? They are not in `genes.tsv`, so their tags would not join. | no: genes only (rule 7) |
| Q2 | Distance unit name: `_distance_nt` matches `terminator_distance_nt` in the same table, while §9.2's example says `distance_bp`. | `_nt` |
| Q3 | Should `output_schema_version` stay at 1.3, or go to 1.4 for an addition after rc1? | stay 1.3 (1.1 unreleased) |
| Q4 | Should `other_overlapping_genes` carry the strand relation per tag (e.g. `T_0002:antisense`)? | no: plain tags, strand in `genes.tsv` |
| Q5 | Should the flanks also give the flanking gene's strand (operon context)? | no (keep light); add later if asked |
| Q6 | Preview tokens `none` / `—` for NA labels: should M12-05c's caption wording decide them instead? | as in rule 9 |

## Lead decisions (2026-11-09)
- Khoa: Q1 "Protein genes only": only genes in `genes.tsv` count as overlapping or flanking genes. Q2 to Q6: the defaults above.
- **Genes preview (Khoa, 2026-11-09, "Locus tag in italics").** In the report's genes preview (`src/gippa/report/build.py` `_genes_preview`), the `gene` cell shows the gene name when there is one. Otherwise it shows the locus tag in italics (`<em>`), with a `title` tooltip of the product and, when Bakta gives one, the UniRef50 ID. This uses the same `gene_label` helper. `genes.tsv` keeps `NA` (§9.2). The preview's sorter (`src/gippa/report/tables.py`, `templates/tables.js`) sorts by the shown label, so these cells are never treated as missing.
  - Why: in the example (Bakta light), 74 of 158 route A CDS have no symbol. Bakta's genome.json holds no gene name for them in any field: 59 are UniRef50 (PSCC) hits with a product only, 17 are hypothetical, 9 have Pfam domains.
  - Acceptance criterion: no preview cell in the `gene` column reads `NA`; a gene with no name shows `<em>{locus_tag}</em>` with the product in its title.
  - This goes in PR A (tests) and PR B (implementation) of this spec.
