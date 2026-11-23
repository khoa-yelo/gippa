# Phase 2 report: PXD066008 (*Haloferax volcanii* H26 on the DS2 genome)

Demo and discovery campaign, Phase 2, group C, dataset 9. The pre-registered answer key is
`ANSWER_KEY.md`. It was committed at `9070f5be` before any spectrum was searched and has not been
edited since. Code: the detached worktree `local/worktrees/rel111-pin` at `bfc74ff6` (GIPPA
1.1.1). Both run manifests record `git_dirty` false. This report holds aggregates only: counts,
shares, q-values, locus tags and 1-based coordinates. No peptide sequence and no DIA-NN row is in
git. Machine-readable numbers are in `metrics.json` (written by `p2c/p2_metrics.py`), `annotation_check.json`
(`p2c/annotation_check.py`) and `provenance.json` (`p2c/collect_provenance.py`).

All shares carry a Wilson 95 % CI in brackets. "unestimated" means neither GIPPA nor this report
gives an error estimate for that number. Input: 5 dia-PASEF runs (2 timsTOF Pro 2, 3 timsTOF
HT). Genome: GCF_000025685.1, 5 replicons, all tagged circular.

## Archaeal annotation check (before any search)

Bakta does not officially support archaea, so the answer key requires two things. First, the
Bakta CDS + sORF count must be within 10 % of the RefSeq non-pseudo CDS count. Second, at least
90 % of RefSeq CDS must share their stop codon (contig, strand, 3' end) with a Bakta CDS.

| | value | rule | result |
|---|---|---|---|
| Bakta CDS + sORF / RefSeq non-pseudo CDS | 4,023 / 3,897 = 1.032 | 0.90-1.10 | pass |
| RefSeq CDS with a Bakta stop | 3,840 / 3,897 = 98.5 % [98.1, 98.9] | ≥ 90 % | pass |

The run went ahead with Bakta. 57 RefSeq CDS have no Bakta stop match (17 of them ≤ 70 aa); see
the exploratory sections. Bakta names almost no ribosomal protein. GIPPA's own ribosomal set
(Bakta names) is empty, and only 1 gene qualifies under the §7.1.2 rule. Row 2 therefore uses the
pre-registered set: 55 RefSeq "(30S|50S) ribosomal protein" CDS, carried to Bakta genes by a
shared stop codon. All 55 map.

## Verdicts

| # | Check | Verdict | Number (error estimate) | Cause / note |
|---|---|---|---|---|
| 1 | Mapping, TE1 (≥ 99.5 % of route-A peptides exact in frame) | **pass** | 31,957 / 31,957 distinct accepted route-A peptides map exactly in frame at their own gene: 100 % [99.988, 100] | 37 peptides with a non-ATG N-terminus, 0 misses (excluding them: 31,920 / 31,920 [99.988, 100]). 2 peptides masked by the pilot's te_proxies rule. 0 q-crossings, 0 exceptions. Placement FDR is not applicable (`--search-db same`): `mapping_fdr.tsv` holds no accepted placement rows |
| 2 | Expression, TE3 (ribosomal ≥ 90 % detected; median abundance in the top quartile) | **pass** | 53 / 55 RefSeq ribosomal proteins detected: 96.4 % [87.7, 99.0]. Median iBAQ log2 20.06 vs Q3 16.75, percentile rank 0.979 (53 measured of 2,988) | Not detected: L40e (FADMKH_00445 / HVO_RS06810) and P1 (FADMKH_02782 / HVO_RS18015). Reported alongside: GIPPA's own set is `not_evaluable` (`no_ribosomal_genes`). The §7.1.2 set is 1 / 1 [20.7, 100]. The pilot's strict product rule on RefSeq names keeps 32 genes: 32 / 32 [89.3, 100], percentile rank 0.978 |
| 3 | Stoichiometry, TE4 | **not testable** (by construction) | rplL/rplJ and atpA/atpD: `gene_missing` for all four | No bacterial L7/L12-L10 or F1 α/β pair exists in archaea, and Bakta names neither, as the answer key stated |
| 4 | Contaminants, TE5 | **pass** | 474 / 474 trypsin (6) and keratin (468) peptides in the contaminant bucket: 100 % [99.2, 100]. 0 placed on the genome | 0 placed tolerant, 0 labelled unmapped, 0 absent. 0 keratin-associated peptides |
| 5 | Gene-level FDR (≤ 0.01, tier `decoy_picked`) | **pass** | 2,978 route-A genes at q ≤ 0.01. Decoy-estimated FDR 0.00974 | 349 paired and 3 unpaired decoys (paired share 0.991). Decoy win weight 1. q floor 3.3e-4. 2,684 genes have ≥ 2 peptides |
| 6 | Operons (RefSeq-named ribosomal and A-type ATP synthase units) | **pass** | 65 genes, 38 adjacent pairs (36 OP, 2 NO). 0 pairs with within-operon q ≤ 0.05 (min q 0.818). 0 `operon_flags` rows, 0 discordant | Median abs Δlog2 0.689 |
| 7 | Genome context (geNomad provirus regions: O < E and q ≤ 0.05) | **not testable** | 0 provirus regions. Regions called: 8 plasmid, 10 BGC, 8 CRISPR. 0 islands | Per-replicon detection is listed below. In route B, one region has q ≤ 0.05 (exploratory) |
| 8 | Negative control: pHV2 (NC_013965.1) | **pass** | 0 / 8 genes detected at q ≤ 0.01 in route A, and 0 / 8 in route B: 0 % [0, 32.4] each | Consistent with H26 being cured of pHV2 |
| 9 | Negative control: pyrE2 (FADMKH_00337, NC_013967.1:301,751-302,281 +) | **pass** | `not_detected` in routes A and B: 0 unique, 0 shared peptides, coverage 0 | Consistent with ΔpyrE2 |
| 10 | Novel ORFs vs RefSeq (share of called ORFs matching a CDS or pseudogene) | **not testable** | 0 ORFs called (142 candidates) | No call to compare, so no known-RefSeq / known-paper / new classing. The exploratory comparison of all candidates is below |
| 11 | µ-proteins Bakta missed (17 RefSeq CDS ≤ 70 aa without a Bakta stop) | **reported** | 3 / 17 have route-B peptide evidence: 17.6 % [6.2, 41.0] (Wilson computed from k, n in `metrics.json`) | All 3 overlap a route-B candidate, and all 3 have exact six-frame peptide hits. Rows are below |
| 12 | Class q floor | **pass** | Route B warns `novel_orfs.class_q_floor_above_alpha`. The four classes have 16-55 candidates (< 100). Min class FDR: 0.048 nterm_extension, 0.083 intergenic and cds_outframe, 0.167 antisense | Alpha 0.01 was not loosened, so there are 0 calls |
| 13 | Route B cost (keeps ≥ 95 % of route-A genes) | **pass** | Route B keeps 2,895 / 2,978 route-A genes: 97.2 % [96.6, 97.7]. 12 genes are detected only in route B | Phase 1: 98.4 %. Route B detects 2,907 genes, FDR 0.0031. Costs are in the route table below |
| 14 | Variants | **not testable** (by construction) | 0 `variants.tsv` rows in route A and in route B | Same-genome database search, no de novo arm |
| 15 | Second round (`makedb --patch`) | **not testable** | not run | Conditional on a route-B call. 0 ORFs were called |

Counts: **9 pass, 0 fail, 5 not testable, 1 reported** (15 pre-registered rows).

## Key numbers

| | Route A (genes) | Route B (+ six-frame ORFs ≥ 20 aa) |
|---|---|---|
| makedb records (incl. 381 `CONT_`) / FASTA bytes | 4,404 / 1,484,825 | 80,841 / 11,936,567 (8.0x bytes) |
| Predicted library: precursors / bytes | 493,854 / 120,777,612 | 3,531,984 / 857,470,627 |
| DIA-NN precursors at 1 % per run (5 runs) | 19,770-29,770 | 18,806-28,320 |
| DIA-NN protein groups, global q ≤ 0.01 | 3,000 | 3,019 |
| GIPPA genes detected (of 4,023; 2 `partial_excluded`) | 2,978 (`decoy_picked`) | 2,907 (`decoy_picked`) |
| Gene FDR at q ≤ 0.01 / detected with ≥ 2 peptides | 0.00974 / 2,684 | 0.00310 / 2,636 |
| Novel-ORF candidates / called | - | 142 / 0 |
| Run warnings | `novel_orfs.decoys_not_makedb`, `novel_orfs.layer_missing` (`rbs_novel`) | `novel_orfs.class_q_floor_above_alpha` |

Abundance agreement between A and B (`gippa compare`): Spearman 0.992 [0.990, 0.994] over 2,887
shared genes (2,205 units). Bland-Altman bias −0.12 log2, limits of agreement [−0.99, 0.74]. All
five stability families are empty in both routes: islands, operon flags, enriched categories,
novel ORFs and alt sites.

## Exploratory: per-replicon detection (no threshold)

| Replicon | Route A | Route B |
|---|---|---|
| NC_013967.1 | 2,326 / 2,981 = 78.0 % [76.5, 79.5] | 2,274 / 2,981 = 76.3 % [74.7, 77.8] |
| NC_013964.1 | 262 / 370 = 70.8 % [66.0, 75.2] | 253 / 370 = 68.4 % [63.5, 72.9] |
| NC_013966.1 | 358 / 590 = 60.7 % [56.7, 64.5] | 348 / 590 = 59.0 % [55.0, 62.9] |
| NC_013968.1 | 32 / 74 = 43.2 % [32.6, 54.6] | 32 / 74 = 43.2 % [32.6, 54.6] |
| NC_013965.1 (pHV2) | 0 / 8 = 0 % [0, 32.4] | 0 / 8 = 0 % [0, 32.4] |

Route B has one region row with q ≤ 0.05: crispr_3 (circular-shift Hodges-Lehmann, effect 4.32,
q 0.049). Route A has none. This is exploratory and borderline.

## Exploratory: novel-ORF candidates vs RefSeq

The same ≥ 50 % same-strand overlap rule was applied to all 142 candidates, called or not:
74 / 142 = 52.1 % [43.9, 60.2] overlap a RefSeq CDS or pseudogene, 4 of them pseudogenes. By
class: nterm_extension 29 / 41, intergenic 21 / 30, cds_outframe 11 / 16, antisense 13 / 55.
Among candidates with ≥ 2 peptides, 23 / 26 match: 88.5 % [71.0, 96.0], none of them a
pseudogene.

As in Phase 1, `phase1_lib._match` counts a candidate with no start codon as a non-match. One
such candidate has ≥ 2 peptides: NOVEL_NC_013967.1_2049171_p (nterm_extension, overlaps
FADMKH_02185, flag `no_start_codon`). The share is therefore conservative.

There are 26 candidates with ≥ 2 peptides: 11 intergenic, 10 nterm_extension, 4 antisense and 1
cds_outframe. Best peptide q runs from 4.0e-7 to 9.4e-4, and class FDR from 0.048 to 0.625.
Rows are in `metrics.json` › `candidates_ge2_peptides`. None is called.

## Exploratory: RefSeq CDS that Bakta missed and route B covers

Of the 57 RefSeq CDS without a Bakta stop match, 16 overlap a route-B candidate (28.1 %
[18.1, 40.8]). 14 have ≥ 1 exact route-B peptide inside the CDS (24.6 % [15.2, 37.1]). Wilson
intervals were computed from k, n in `metrics.json`. All 16 lie on NC_013967.1 (8), NC_013966.1
(7) or NC_013968.1 (1).

µ-proteins (≤ 70 aa, row 11), the 3 of 17 with evidence:

| RefSeq | coordinates (strand) | aa | product | route-B candidate | exact peptides |
|---|---|---|---|---|---|
| HVO_RS20940 | NC_013967.1:1,299,138-1,299,308 (+) | 56 | hypothetical protein | NOVEL_NC_013967.1_1299120_p (intergenic, class FDR 0.091) | 2 |
| HVO_RS21530 | NC_013967.1:1,531,014-1,531,148 (+) | 44 | rubrerythrin-like domain-containing protein | NOVEL_NC_013967.1_1531014_p (cds_outframe, 0.083) | 2 |
| HVO_RS20980 | NC_013967.1:1,660,355-1,660,495 (+) | 46 | hypothetical protein | NOVEL_NC_013967.1_1660307_p | 1 |

The other 13 rows are above 70 aa (71-293 aa). Examples: HVO_RS15475 (71 aa, 4 peptides),
HVO_RS21365 (145 aa, 5), HVO_RS21450 (114 aa, 5), HVO_RS19995 (229 aa, 15) and HVO_RS03855
(universal stress protein, 293 aa, 15). HVO_RS19780 and HVO_RS20470 overlap a candidate but have 0
exact peptides inside the CDS. Full rows: `metrics.json` › `refseq_cds_missed_by_bakta`.

The answer key's "known-paper" class needs the paper's µ-protein list. Because no ORF was called,
no candidate was classed against it.

## Command log (pinned commit `bfc74ff6`, GIPPA 1.1.1)

All commands ran from `local/worktrees/rel111-pin` with `PYTHONPATH=<pin>/src`.
`$D` = `$GIPPA_ROOT/campaign/PXD066008`. Bakta 1.12.2 (DB 6.0 full) ran in job 2001987_15
(from `ANSWER_KEY.md`). The route jobs run steps 3-7 of `p2c/route.sbatch`: makedb, DIA-NN
(`p2c/diann_search.sbatch`, flags from `p2c/diann_campaign.cfg`), convert twice plus preflight
(`p2c/convert_preflight.sbatch`) and `gippa run`.

```bash
# job 2008625 (route A, wheat-06, 16 CPUs): total 16,267 s
gippa makedb --annotation $J --genome $G --outdir $D/db_a                                  # 3 s
#   DIA-NN 2.3.1: step 1 (library) 88 s, step 2 (5 runs) 1,567 s, exit 0
gippa convert --from diann $DN/report.parquet --out $IN/peptides.tsv.gz \
  --abundance-out $IN/abundance.tsv.gz --search-db same          # twice, deterministic (byte-identical), 5 s
gippa preflight ... --strict                                     # exit 0, 6 s
gippa run --profile farmshare_local,apptainer --genome $G --annotation $J \
  --peptides $IN/peptides.tsv.gz --db-dir $GIPPA_DB_DIR --threads 14 --memory 54.GB --outdir $D/out_a   # exit 0, 14,486 s
# job 2009031 (route B, wheat-04, 16 CPUs): total 20,550 s
gippa makedb --annotation $J --genome $G --six-frame-orfs --min-orf-aa 20 --outdir $D/db_b  # 9 s
#   DIA-NN 2.3.1: step 1 542 s, step 2 5,543 s, exit 0
gippa convert ... --search-db same                               # twice, deterministic, 5 s
gippa preflight ... --strict                                     # exit 0, 7 s
gippa run ... --outdir $D/out_b                                  # exit 0, 14,415 s
# p2c/post.sbatch steps, run in 2005379 (iron-04, the session allocation; see Deviations)
gippa compare $D/out_a $D/out_b --labels genes,sixframe --outdir $D/cmp_ab   # exit 0, 10 s
gippa report $D/out_a --redact                                   # exit 0, 47 s, tables/ unchanged (sha256)
gippa report $D/out_b --redact                                   # exit 0, 71 s, tables/ unchanged (sha256)
```

Both runs: all modules (`modules_off` empty), master seed 20261014, no `--abundance`, no
`--strict`, exit 0.

## Provenance

| item | value |
|---|---|
| GIPPA | `bfc74ff6c20acc2cfcf298cdf23cbb30d8fae0dc` (1.1.1), `git_dirty` false |
| genome | sha256 `94be40697a6f799d4f66c8ead03cfaeec62f78e065b60c7e8bddbcb71431fc56` |
| Bakta JSON | sha256 `3fe1e7fb6098e196e083f8a498ba5c047f93d17136399fdad3fadee00f12823f` |
| DIA-NN | 2.3.1 Academia, image sha256 `1946452d56d2fb0f66e05ced473df88e53484e28040f1ecbd83c4a43405a6381`; cfg sha256 `dad24b38e315b68cb6a905dfda877404535310f030ed9f92f8e9370dae0d6a1f` |
| FASTA A / B | `3fafdca7c829944488988523f8b69b67a927fd0ae2afa0a35d2d37c4920b4d06` / `e187d7de88d40e391d837ee065c25c2e78f45a7c891ed450c3b679dd2e802087` |
| DIA-NN report A / B | `8dfa3376e4425e3efa2d787daa78426ca96c8ba07da7eb024b5bb873ac6fc107` / `6cf2230149215d5e22820bbbc138c1c145b55ec0dc1c93d017f119300ccfb6dc` |
| peptides A / B | `46cfeab59deae3d191f7dcda0b027c0b69252031897b624545c7cda6185d2d26` / `0ded6e2aa415990989f7477548d236120f122f8de34b2e2965e1a8db61515649` |
| abundance A / B | `4d91c5a2580180de6dd36412324f3477163f1bcbcb684ae8fe30faae03288f95` / `49bfbd497e17d931e57745bf4548a5f568d80b51d2fcb8b4938172ad4f8768e8` |
| raw (5 `.d.tar`) | sha256 per file in `provenance.json` › `raw_sha256` |
| containers (same in A and B) | antiSMASH 8.0.4 `sha256:d3f390e6…a29d`; geNomad 1.12.0 `sha256:d54652ff…e5c`; OSTIR 1.1.3 and ViennaRNA 2.6.4 `sha256:d9107b9b…9292`; Promotech 56251ad9 `sha256:7366eaf6…53ac`; Promoter Calculator 1.2.4 `sha256:c5bd0662…48df`; TerminatorNet 72180309 `sha256:ec415757…7a`; TransTermHP 2.09 `sha256:fc6fe040…fd03`; UniOp 3ed14dfc (no container digest). Full digests are in `provenance.json` › `gippa_run.tools` |

## Deviations and problems

- **TE3 set (as pre-registered).** Bakta names almost nothing here: GIPPA's ribosomal set is
  empty, and only 1 gene qualifies under the §7.1.2 rule. So TE3 uses the 55 RefSeq ribosomal
  proteins carried to Bakta genes by a shared stop codon (same contig, strand and 3' end), as
  `ANSWER_KEY.md` fixed before the search. The operon set (row 6) uses the same carried names.
  This is not GIPPA's own TE3 rule, and GIPPA's own result (`not_evaluable`) is reported alongside.
- **Compare and report in the session allocation.** `gippa compare` and both `gippa report
  --redact` calls (the `p2c/post.sbatch` steps) ran in job 2005379 on iron-04. That is the
  interactive session allocation, not a separate `gippa-*` batch job. Total wall time was
  about 2 min (10 s + 47 s + 71 s). `tables/` were checked unchanged by sha256 after each report.
- **TE7 rerun not done.** No determinism rerun of route A was made for this dataset. Determinism is
  covered only by convert running twice byte-identical in each route.
- **No second novel-ORF round.** Route B called 0 ORFs, so `makedb --patch` and the round-2
  search (row 15) were not run.
- **Contaminants.** As in Phase 1, `gippa makedb` appended the 381-protein library as `CONT_`
  records, so DIA-NN got one FASTA. The TE5 families were assigned by the pilot rule on
  `contaminants.fasta`, with `CONT_` IDs.
- **Route A warnings.** `novel_orfs.decoys_not_makedb` (route A has no ORF decoys) is expected.
  Route A also warns `novel_orfs.layer_missing` (`rbs_novel`), which E1 in Phase 1 did not. It
  touches no route-A verdict. Route B has no missing layer.
- **Profile.** Both runs used `farmshare_local,apptainer`: all tasks inside one 16-CPU job, with
  the `--threads 14 --memory 54.GB` settings shown in the command log. The route A vs B wall
  times are comparable with each other. They are not comparable with Phase 1's 1.1.0 times
  (distributed `farmshare` profile).
- **Raw staging.** Route B's `raw_s` is 0 because route A's job had already staged the runs
  (101 s).
