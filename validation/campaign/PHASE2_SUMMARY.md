# Phase 2 summary: ten-dataset campaign (cross-dataset)

Campaign Phase 2 (handover plan of 2026-11-17). This file summarises the ten dataset reports. It adds no new runs or analyses and changes no threshold. It holds aggregates only: counts, shares, q-values, locus tags and 1-based coordinates. It contains no peptide sequences.

**Sources.** Paths are relative to `validation/campaign/`. At the time of writing, group A files are in the `campaign-p2a` worktree, group B in `campaign-p2b`, group C in `campaign-p2c`, and Phase 1 in `campaign-ecoli`.
- Group A: `PXD065022/` (T. pallidum, 7 strains), `PXD069004/` (C. crescentus), `PXD055395/` (B. thetaiotaomicron).
- Group B: `PXD072227/` (Synechococcus CC9311), `PXD049476/` (V. natriegens), `PXD060688/` (M. xanthus), `PXD065346/` (C. difficile).
- Group C: `PXD057335/` (M. tuberculosis H37Rv), `PXD066008/` (H. volcanii), `PXD079213/` (S. islandicus).
- Each dataset directory holds `ANSWER_KEY.md`, `REPORT.md` and `metrics.json`.
- Phase 1: `ecoli/REPORT.md`, `ecoli/REPORT_1b.md`.
- Short forms used below: "[Tpa R]" = `PXD065022/REPORT.md`, "[Tpa M]" = `PXD065022/metrics.json`, and likewise Ccre, Bthe, Syn, Vnat, Mxan, Cdif, Mtb, Hvo, Sisl; "[P1]" / "[P1b]" = the Phase 1 reports.

**Shares** carry Wilson 95 % CIs as given in the reports. "computed here" marks a number derived in this summary from counts in a report or metrics.json. "unestimated" means no error estimate exists.

**Code.**
- 1.1.1 release code `bfc74ff6`: C. crescentus, B. theta, all of group B, all of group C.
- main `38e1cdba` (1.1.1 plus the gap fixes #612, #616, #617): T. pallidum. The gapped assemblies could not finish at `bfc74ff6` or `cff0fd68` [Tpa R, Problems]. Its manifests report `gippa_version` 1.1.0, because the version bump is in release PR #609.

## 1. Verdict table

The answer keys number their checks differently by group, so the columns below are check families.

| Family | Group A ID | Group B ID | Group C ID (Mtb / Hvo / Sisl) |
|---|---|---|---|
| TE1 mapping | S1 | S1 | 1 / 1 / 1 |
| TE3 expression | S2 | S2 | 2 / 2 / 2 |
| TE4 stoichiometry | S3 | S3 | 3 / 3 / 3 |
| TE5 contaminants | S4 | S4 | 4 / 4 / 4 |
| Gene FDR | S10 (reported) | S5 (pass/fail) | 5 / 5 / 5 |
| Operons | S5 | S7 | 6 / 6 / 6 |
| Genome context | S6 (reported) | S6 (prophage rule) | 9 / 7 / 7 (prophage rule) |
| Novel ORFs vs RefSeq | S7 | S8 | 10 / 10 / 9 |
| Class q floor | S8 | S9 | 11 / 12 / 10 |
| Route B cost | S9 (reported) | S10 (reported) | 12 / 13 / 11 (pass if ≥ 95 % kept) |
| Second round | S11 | S11 | 14 / 15 / 13 |
| Variants (same genome) | not in key | S12 | 13 / 14 / 12 (not testable by design) |
| TE7 determinism | T7 (Tpa only) | S13 (Syn only) | not in key |

Abbreviations: P pass, F fail, NT not testable, NR not run, R reported, — not in that key.

| Dataset | TE1 | TE3 | TE4 | TE5 | Gene FDR | Operons | Genome context | Novel vs RefSeq | Class floor | Route B | 2nd round | Variants | TE7 | Dataset checks | P / F / NT+NR / R | Commit |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| T. pallidum (7 strains) | P 7/7 | P 7/7 | **F** (6/7; Atp NT) | P 7/7 | R | P 7/7 | R | NT | P 7/7 | R | NR | — (T6 placement FDR P 6/6) | P | T1 P, T2 P, T3 P, **T4 F**, **T5 F**, T6 P | 10 / 3 / 2 / 3 | `38e1cdba` |
| C. crescentus | P | P | P | P | R | P | R | NT | P | R | NR | — | — | C1 P, C2 P, C3 R, C4 P, C5 R | 10 / 0 / 2 / 5 | `bfc74ff6` |
| B. thetaiotaomicron | P | P | P | P | R | P | R | NT | P | R | NR | — | — | B1 P, B2 P, B3 R, B4 P, B5 R | 9 / 0 / 2 / 5 | `bfc74ff6` |
| Synechococcus CC9311 | P | P | P | P | P | P | NT | NT | P | R | NR | P | P | Y1 P, Y2 NR | 10 / 0 / 4 / 1 | `bfc74ff6` |
| V. natriegens | P | P | **F** | P | P | P | P | NT | P | R | NR | P | — | V1 P, V2 P, V3 R, V4 NR | 10 / 1 / 3 / 2 | `bfc74ff6` |
| M. xanthus | P | P | **F** | P | P | P | P | NT | P | R | NR | P | — | M1 R, M2 R, M3 P | 9 / 1 / 2 / 3 | `bfc74ff6` |
| C. difficile | P | P | **F** | P | P | P | **F** | NT | P | R | NR | P | — | C1 R, C2 NT, C3 NT, C4 NT, C5 R, **C6 F** | 7 / 3 / 5 / 3 | `bfc74ff6` |
| M. tuberculosis | P | P | **F** | P | P | P | NT | NT | P | P | NT | NT | — | 7 ESX-1 P, 8 PE/PPE P, 15 topology P | 10 / 1 / 4 / 0 | `bfc74ff6` |
| H. volcanii | P | P | NT | P | P | P | NT | NT | P | P | NT | NT | — | 8 pHV2 P, 9 pyrE2 P, 11 R | 9 / 0 / 5 / 1 | `bfc74ff6` |
| S. islandicus | P | P | NT | P | P | P | NT | NT | P | P | NT | NT | — | **8 neg. controls F** | 7 / 1 / 5 / 0 | `bfc74ff6` |
| **Total** | | | | | | | | | | | | | | | **91 / 10 / 34 / 23** (158 checks) | |
| *Phase 1 E. coli, for reference* | P | P | F | P | R | P | F | NT | P | R | NR | F (E4, TE2) | P | | 10 / 3 / 1 NT, 1 dropped / 1 | `a2a95c76`, rechecked at `bfc74ff6` |

Counts are copied from each report's "Counts" line.

What the 10 Phase 2 fails are:

| Fail | Kind | Where explained |
|---|---|---|
| TE4 RplL/RplJ: Tpa (Madras only), Vnat, Mxan, Cdif, Mtb | known-risk ratio, cause not established | §2.5 |
| Cdif S6 prophage enrichment | test-design limit (whole-replicon region) plus small regions | §3 |
| Cdif C6 PrdB C-terminus | answer-key premise error | §3 |
| Sisl row 8 pyrE/pyrF detected | vector caveat stated in the key; post hoc check fits it | §3 |
| Tpa T4 pseudogene distal parts (0/2) | one half untestable by construction (Bakta pseudogene); the other half has 0 peptides | §3 |
| Tpa T5 Nichols vs SS14 abundance (3/5) | both misses are small differences (−0.18, −0.12 log2), error unestimated | [Tpa R] T5 |

## 2. Cross-dataset patterns

### 2.1 TE1 mapping

- 100 % of accepted route-A peptides map exactly in frame at their own gene, in all 16 genomes. That is 302,396 peptides in total (computed here: 44,971 over the 7 T. pallidum strains plus 257,425 over the other 9). There are 0 misses of any kind.
- Lowest Wilson lower bound: 99.860 % (T. pallidum Grady, 2,746 peptides) [Tpa R S1]. Largest set: M. xanthus, 54,278 / 54,278 [99.993, 100] [Mxan R S1].
- Peptides with a non-ATG N-terminus all mapped: 15 to 37 per T. pallidum strain, and up to 111 in M. tuberculosis [Mtb R row 1].
- No tolerant placement was accepted in any same-genome run. Example tolerant strata: C. crescentus 1 target / 28 decoys [Ccre R S1]; S. islandicus 1 / 16 [Sisl R row 1].

### 2.2 Gene-level FDR tiers and decoy FDR

Tier `decoy_picked` in every run (34 runs: A and B per genome).

| Dataset | Route A genes / FDR at q ≤ 0.01 | Route B genes / FDR | Source |
|---|---|---|---|
| T. pallidum (7) | 668-880 / 0.0011-0.0060 | 663-882 / 0.0011-0.0030 | [Tpa R S10] |
| C. crescentus | 3,005 / 0.0093 | 2,938 / 0.0037 | [Ccre R S10] |
| B. theta | 2,244 / 0.0098 | 2,235 / 0.0094 | [Bthe R S10] |
| Synechococcus | 1,897 / 0.0084 | 1,857 / 0.0054 | [Syn R S5] |
| V. natriegens | 2,515 / 0.0099 | 2,468 / 0.0093 | [Vnat R S5] |
| M. xanthus | 4,864 / 0.0097 | 4,941 / 0.0091 | [Mxan R S5] |
| C. difficile | 1,928 / 0.0099 | 1,968 / 0.0097 | [Cdif R S5] |
| M. tuberculosis | 2,942 / 0.0099 | 2,889 / 0.0066 | [Mtb R row 5] |
| H. volcanii | 2,978 / 0.00974 | 2,907 / 0.00310 | [Hvo R row 5] |
| S. islandicus | 2,109 / 0.0085 | 2,076 / 0.0053 | [Sisl R row 5] |

- Every route-A estimate is ≤ 0.0099. Route B's estimate is equal or lower in every dataset.
- The decoy paired share is 0.990 to 0.993 where reported (Mtb, Hvo, Sisl), and the decoy win weight is 1 (ADR-0078) [Mtb R row 5; Hvo R row 5; Sisl R row 5].
- T. pallidum route-A FDRs are low (0.0011-0.0060) in a search dominated by rabbit proteins (§2.3). These reports do not analyse why.
- Placement decoy FDR (pooled TE8 estimator) in the six T. pallidum cross-strain arms: 0.059 % to 0.101 %, all ≤ 1 %. Every accepted placement was exact [Tpa R T6; Tpa M `T6`].

### 2.3 Detection depth

| Dataset | Runs, instrument | Precursors at 1 % per run (A) | Route-A detected / coding genes | Source |
|---|---|---|---|---|
| T. pallidum (7) | 3 per strain, timsTOF Pro | 3,080-33,704 (rabbit-dominated) | 68.5 % (Grady) to 90.3 % (Haiti B) | [Tpa R T2] |
| C. crescentus | 3, Exploris 480 | 37,937-39,493 | 3,005 / 3,723 = 80.7 % [79.4, 81.9] | [Ccre R C3] |
| B. theta | 16, timsTOF Pro, 11.5 min | 14,722-21,355 | 2,244 / 4,931 = 45.5 % [44.1, 46.9] | [Bthe R B1, B4] |
| Synechococcus | 5, Q Exactive | 15,914-17,535 | 1,897 / 2,859 = 66.4 % [64.6, 68.1] (computed here) | [Syn R key numbers] |
| V. natriegens | 4, Exploris 480 | 23,231-27,522 | 2,515 / 4,551 = 55.3 % [53.8, 56.7] (computed here) | [Vnat R key numbers] |
| M. xanthus | 4, Exploris 480 | 62,329-65,562 | 4,864 / 7,318 = 66.5 % [65.4, 67.5] (computed here) | [Mxan R key numbers] |
| C. difficile | 3, DIA | 16,951-22,839 | 1,928 / 3,824 = 50.4 % [48.8, 52.0] (computed here) | [Cdif R key numbers] |
| M. tuberculosis | 3, Exploris 480 | 24,264-31,130 | 2,942 / 4,075 = 72.2 % [70.8, 73.6] | [Mtb R key numbers] |
| H. volcanii | 5, timsTOF Pro 2 / HT | 19,770-29,770 | 2,978 / 4,023 = 74.0 % [72.6, 75.4] (computed here; 2 genes `partial_excluded`) | [Hvo R key numbers] |
| S. islandicus | 3, Astral, 8 min | 28,553-28,926 | 2,109 / 2,807 = 75.1 % [73.5, 76.7] | [Sisl R key numbers] |
| *E. coli (Phase 1)* | 4, Eclipse | 35,342-40,070 | 2,691 / 4,368 = 61.6 % [60.2, 63.0] (computed here) | [P1 key numbers] |

- T. pallidum is host-dominated: 65.9 % to 77.3 % of the distinct accepted route-A target peptides carry only `CONT_` proteins (rabbit plus the contaminant library) [Tpa R T1].
- Ribosomal proteins (TE3) passed everywhere. Detected share: 90.2 % (Grady 46/51; M. xanthus 55/61, where the CI's lower bound is 80.2 %) to 100 %. Median iBAQ percentile rank: 0.78 to 0.99 [each report, TE3 row].
- In the archaea, Bakta names almost no ribosomal protein. TE3 used RefSeq names carried to Bakta genes by a shared stop, and GIPPA's own set was `not_evaluable` [Hvo R row 2; Sisl R row 2].

### 2.4 Route B: gene retention and cost

| Dataset | Kept share of route-A genes | B-only genes | Records B/A | FASTA bytes B/A | Library precursors A → B (ratio) | DIA-NN total A → B (ratio) | `gippa run` A / B |
|---|---|---|---|---|---|---|---|
| T. pallidum (7) | 97.3 % (Grady) to 99.9 % (DAL-1) | 3-13 | 23.4x | 6.5x | 5.30 M → 6.20 M (1.17x, rabbit included) | step 2: 2 h 28-3 h 40 → 2 h 51-5 h 00 (about 1.2x) | 35-55 / 35-52 min |
| C. crescentus | 97.1 % [96.4, 97.6] | 21 | 21.8x | 7.9x | 0.69 M → 4.12 M (6.0x) | 13 min 42 s → 48 min 09 s (3.5x) | 3 h 33 / 3 h 21 |
| B. theta | 97.2 % [96.4, 97.8] | 54 | 33.5x | 8.1x | 1.13 M → 5.25 M (4.6x) | 1 h 05 min 53 s → 3 h 48 min 35 s (3.5x) | 2 h 40 / 2 h 47 |
| Synechococcus | 96.8 % [96.0, 97.5] | 20 | 24x | 7.9x | 0.42 M → 2.37 M (5.6x) | 19 min 54 s → 1 h 15 min 18 s (3.8x, computed here) | 1 h 28 / 1 h 25 |
| V. natriegens | 96.9 % [96.1, 97.5] | 32 | 31x | 8.3x | 0.78 M → 4.08 M (5.2x) | 19 min 35 s → 59 min 17 s (3.0x, computed here) | 2 h 10 / 2 h 11 |
| M. xanthus | 98.5 % [98.1, 98.8] | 150 | 22.7x | 8.2x | 1.48 M → 9.80 M (6.6x) | 47 min 16 s → 4 h 05 min 05 s (5.2x, computed here) | 8 h 10 / 7 h 54 |
| C. difficile | 98.8 % [98.2, 99.2] | 63 | 24x | 6.6x | 0.75 M → 2.33 M (3.1x) | 36 min 57 s → 1 h 30 min 11 s (2.4x, computed here) | 3 h 05 / 2 h 41 |
| M. tuberculosis | 97.7 % [97.1, 98.2] | 14 | 23.1x | 8.1x | 0.63 M → 4.30 M (6.8x) | 23 min 25 s → 1 h 06 min 43 s (2.8x) | 5 h 17 / 3 h 12 (see caveats) |
| H. volcanii | 97.2 % [96.6, 97.7] | 12 | 18.4x (computed here) | 8.0x | 0.49 M → 3.53 M (7.2x, computed here) | 27 min 35 s → 1 h 41 min 25 s (3.7x, computed here) | 4 h 01 / 4 h 00 |
| S. islandicus | 98.3 % [97.6, 98.8] | 3 | 19.6x | 6.0x | 0.46 M → 1.42 M (3.1x) | 8 min 16 s → 12 min 08 s (1.5x, computed here) | 1 h 10 / 1 h 01 |
| *E. coli (Phase 1)* | 98.4 % [97.9, 98.8] | 21 | 28x | 8.0x | 0.74 M → 4.09 M (5.5x) | 2.8x | 2 h 32 / 2 h 16 |

Sources: the "Route A vs route B" table and the route-B cost row of each report. DIA-NN totals are step 1 plus step 2 on 16 CPUs.

- Route B kept at least 96.8 % of route-A genes everywhere. The ≥ 95 % prior or threshold held in every dataset where one existed (group B prior; group C rows 12/13/11).
- In M. xanthus, C. difficile and three T. pallidum strains, route B detects more genes than route A [Mxan R S10; Cdif R key numbers; Tpa R S10].
- Abundance agreement between the routes: Spearman 0.971 (T. pallidum Grady) to 0.996 (S. islandicus). Bland-Altman bias −0.25 (SS14) to +0.08 log2 [S9/S10/row 11-13 of each report].
- Cost is paid in the DIA-NN search, not in GIPPA: `gippa run` times for A and B are similar everywhere except M. tuberculosis (caveats).
- The M. xanthus lane peaked at 60.9 GB MaxRSS within a 64,000 MB allocation [Mxan R M1].

### 2.5 RplL/RplJ stoichiometry (TE4)

Band [1, 3] on log2(RplL/RplJ), iBAQ (molar) on route A. Peptide counts are `n_peptides_unique` / `n_observable_peptides` from each metrics.json (TE4 `secondary` block or `te4_name_based`).

| Genome | log2 iBAQ (top3) | Verdict | rplL uniq / obs | rplJ uniq / obs | AtpA/AtpD log2 |
|---|---|---|---|---|---|
| T. pallidum DAL-1 | 2.72 (2.05) | pass | 11 / 6 | 7 / 10 | NT (V-type) |
| T. pallidum Grady | 2.87 (1.55) | pass | 3 / 6 | 2 / 10 | NT |
| T. pallidum Haiti B | 1.41 (0.98) | pass | 8 / 6 | 13 / 10 | NT |
| T. pallidum Madras | **3.43** (2.63) | **fail** | 7 / 6 | 3 / 10 | NT |
| T. pallidum Mexico A | 2.89 (1.97) | pass | 8 / 6 | 6 / 10 | NT |
| T. pallidum Phil 1 | 2.96 (2.41) | pass | 10 / 6 | 9 / 10 | NT |
| T. pallidum SS14 | 2.76 (1.77) | pass | 4 / 6 | 4 / 10 | NT |
| C. crescentus | 2.68 (1.68) | pass | 8 / 4 | 14 / 9 | 0.89 |
| B. theta | 1.44 (1.38) | pass | 14 / 8 | 10 / 8 | −0.54 |
| Synechococcus | 2.53 (1.54) | pass | 10 / 7 | 11 / 12 | 0.14 |
| V. natriegens | **0.96** (0.91) | **fail** | 6 / 5 | 9 / 5 | 0.19 |
| M. xanthus | **0.16** (−0.25) | **fail** | 10 / 6 | 12 / 10 | 0.03 |
| C. difficile | **−1.22** (−0.98) | **fail** | 16 / 8 | 22 / 10 | 0.12 |
| M. tuberculosis | **0.61** (0.58) | **fail** | 12 / 6 | 18 / 9 | −0.02 |
| H. volcanii, S. islandicus | NT | NT | — | — | NT (archaeal stalk, A-type ATPase) |
| *E. coli (Phase 1)* | −0.33 (−0.12) | fail | – / 7 | – / 8 | −0.30 |

Sources: [Tpa M] `strains.*.metrics.S3_te4`; [Ccre M], [Bthe M] `S3_te4`; [Syn M], [Vnat M], [Mxan M], [Cdif M] `S3_te4`; [Mtb M] `row_te4`; [P1] row 3.

What the data show:
- **AtpA/AtpD passed in every genome where it could be tested** (8 of 8 Phase 2 genomes plus E. coli). The problem is specific to RplL/RplJ.
- **Pass or fail does not follow peptide counts.** C. difficile has the most rplL unique peptides (16) and the lowest ratio (−1.22). T. pallidum Grady has 3 and passes.
- **Grouping by genome.** The pair passes in C. crescentus, B. theta, Synechococcus and 6 of 7 T. pallidum strains. It fails in E. coli, V. natriegens, M. xanthus, C. difficile and M. tuberculosis. Nine genomes are too few to separate a lineage effect from chance, and these data test no mechanism.
- **The verdict depends on the abundance measure near the band edges.** Madras fails on iBAQ (3.43) but its top3 value (2.63) is inside the band. Haiti B passes on iBAQ (1.41) but its top3 value (0.98) is outside. V. natriegens misses by 0.04 log2.

What might explain it (not shown by these data):
- L7 is the N-terminally acetylated form of L12, and N-terminal acetylation is not a search modification here, so part of the L7/L12 signal goes unsearched [P1 row 3; Mxan R S3; Mtb R row 3].
- No Phase 2 run measured the acetylated share, so this remains a candidate explanation only.
- The few-peptides explanation offered in Phase 1 does not fit C. difficile (16 unique rplL peptides) [Cdif R S3].

### 2.6 Novel ORFs: 0 called everywhere

| Dataset | Candidates | With ≥ 2 peptides | Min class FDR by class | Candidates matching RefSeq (≥ 50 % same strand) |
|---|---|---|---|---|
| T. pallidum (7) | 200 (17-41 per strain) | 36 | 0.071 (DAL-1 nterm_extension) to 0.33 | 103 / 200 (51.5 %); ≥ 2 pep 30 / 36 |
| C. crescentus | 164 | 26 | nterm 0.026 (49 cand.), intergenic 0.077, cds_outframe 0.17, antisense 0.76 | 65 / 164, 39.6 % [32.5, 47.3]; ≥ 2 pep 22 / 26 |
| B. theta | 99 | 6 | nterm 0.11, others 0.36 | 22 / 99, 22.2 % [15.2, 31.4]; ≥ 2 pep 5 / 6 |
| Synechococcus | 117 | 26 | nterm 0.029 (43), antisense 0.14, intergenic 0.14, cds_outframe 0.17 | ≥ 2 pep 20 / 26 |
| V. natriegens | 78 | 7 | nterm 0.17, intergenic 0.50, antisense 0.67, cds_outframe 1.0 | ≥ 2 pep 7 / 7 (all inside RefSeq genes) |
| M. xanthus | 265 | 45 | **nterm 0.0127 (115 cand.)**, intergenic 0.14, antisense 0.17, cds_outframe 0.17 | ≥ 2 pep 31 / 45 |
| C. difficile | 53 | 4 | nterm 0.14, antisense 0.50, intergenic 0.50, cds_outframe 1.0 | ≥ 2 pep 4 / 4 |
| M. tuberculosis | 247 | 47 | **nterm 0.0118 (109 cand.)**, others 0.091 | 127 / 247, 51.4 % [45.2, 57.6]; ≥ 2 pep 39 / 47 |
| H. volcanii | 142 | 26 | nterm 0.048, intergenic 0.083, cds_outframe 0.083, antisense 0.167 | 50 / 142, 35.2 % [27.8, 43.4] (corrected 2026-11-24: same-replicon rule; was 74 / 142); ≥ 2 pep 23 / 26 (not rechecked) |
| S. islandicus | 92 | 16 | nterm 0.05, intergenic 0.143, cds_outframe 0.143, antisense 0.857 | 25 / 92, 27.2 % [19.1, 37.0]; ≥ 2 pep 5 / 16 |
| **Total Phase 2** | **1,457** | **239** | none ≤ 0.01 | |
| *E. coli (Phase 1)* | 124 | 6 | 0.14 to 0.67 | 23 / 124, 18.5 % |

Sources: the S7/S8, S8/S9 and rows 10-11 of each report, and their candidate sections. Totals computed here.

- **Why 0 calls.**
  - Every route-B run carries `novel_orfs.class_q_floor_above_alpha`, and alpha 0.01 was not loosened (class-floor check passed in 10/10).
  - Nearly every class has fewer than 100 candidates. With so few, the class FDR cannot reach 0.01.
  - The only classes with ≥ 100 candidates (M. xanthus and M. tuberculosis nterm_extension) came closest, at 0.0127 and 0.0118, still above alpha. Those two classes fall outside the key's "< ~100" wording [Mxan R S9; Mtb R row 11].
- **The second round never ran**, because nothing was called (10/10).
- **Classes behave consistently across genomes.** nterm_extension has the lowest class FDR in every dataset and mostly matches RefSeq starts (T. pallidum 50/56, C. crescentus 39/49, M. tuberculosis 93/109, H. volcanii 29/41). Antisense has the highest class FDRs (0.14-1.0) and rarely matches RefSeq (0/38 in B. theta, 1/65 in C. crescentus, 2/76 in M. tuberculosis, 0/26 in S. islandicus, 7/73 in T. pallidum).
- **The RefSeq match shares are low by construction.** A candidate without a start codon counts as a non-match [Mtb R; Hvo R; Sisl R].

### 2.7 Genome-context enrichment

**Prophage rule** (group B and C keys; also Phase 1): every geNomad prophage region has O < E and q ≤ 0.05.

| Dataset | Prophage regions | O < E | q ≤ 0.05 | Verdict |
|---|---|---|---|---|
| V. natriegens | 2 | 2 / 2 (O/E 0.19, 0.25) | 2 / 2 (q 0.015, 0.015) | pass [Vnat R S6] |
| M. xanthus | 2 | 2 / 2 (O/E 0.15, 0.15) | 2 / 2 (q 0.008, 0.004) | pass [Mxan R S6] |
| C. difficile | 4 (3 chromosome, 1 whole plasmid) | 4 / 4 (O/E 0.00-0.39) | 0 / 4 (q 0.079 ×3, 1.0) | fail [Cdif R S6] |
| Synechococcus, M. tuberculosis, H. volcanii, S. islandicus | 0 called (phiRv1/phiRv2 not called in Mtb) | — | — | NT |
| B. theta (group A, reported only) | 16 | not reported | 0 / 16 | R [Bthe R S6] |
| *E. coli (Phase 1)* | 5 | 5 / 5 | 2 / 5 | fail [P1 row 5] |

- Across every tested region (E. coli, V. natriegens, M. xanthus, C. difficile), the direction holds in 13 / 13 and significance in 6 / 13 (computed here).
- Regions that fail significance are small or a whole replicon (§3).

Other context, reported or exploratory:
- **V. natriegens.** Detection on chromosome 1 vs chromosome 2: 65.3 % vs 38.0 %, difference +27.3 points [24.4, 30.2]. Median abundance difference 1.98 log2 [1.68, 2.31] [Vnat R V1, V2].
- **M. xanthus.** antiSMASH regions are slightly *better* detected: 71.4 % inside vs 65.9 % outside, +5.5 points [+2.1, +8.7]. The pre-registered prior was the opposite. 0 / 27 regions have q ≤ 0.05 [Mxan R M2].
- **C. crescentus.** `bgc_3` is under-detected: O/E 0.33, q 0.032. Two islands (ISL_1 high, ISL_2 low; island q unestimated) [Ccre R S6].
- **B. theta.** Plasmid detection 39.0 % [25.7, 54.3] vs chromosome 45.6 % [44.2, 47.0]; the CIs overlap [Bthe R B3].
- **H. volcanii.** Detection by replicon runs from 78.0 % (main chromosome) to 43.2 % (NC_013968.1). pHV2 0/8 [Hvo R].
- **Exploratory, not pre-registered.** M. tuberculosis `crispr_2` q 0.022; H. volcanii `crispr_3` route B q 0.049; S. islandicus ISL_1 low island (no q) [Mtb R; Hvo R; Sisl R row 7].

### 2.8 Operons

Ribosomal-protein and ATP-synthase adjacent pairs passed in all 10 datasets (16 genomes).
- 27 to 41 pairs per genome, 0 pairs with within-operon q ≤ 0.05.
- Minimum q ranges from 0.65 (C. crescentus) to 0.94 (B. theta).
- No discordant `operon_flags` row among these pairs [S5/S7/row 6 of each report].

Genome-wide `conflict_terminator_in_operon` flags: 1 per T. pallidum strain, B. theta 14, V. natriegens 16, M. xanthus 3, C. difficile 36. C. crescentus, Synechococcus, M. tuberculosis, H. volcanii and S. islandicus have 0 [same rows].

This check only shows the absence of discordance. It does not test the power to detect it.

### 2.9 Placement and variants

- **Same-genome runs.** 0 `variants.tsv` rows in every run that reports it: group B S12 pass in 4/4; group C 0 rows, not testable by design [Syn/Vnat/Mxan/Cdif R S12; Mtb/Hvo/Sisl R].
- **Cross-strain runs** (T. pallidum T6; six strains on the DAL-1 genome, `--search-db other --genome-relation related`):
  - Recall is 0 on covered truth sites in every arm: 0/18, 0/21, 0/2, 0/48, 0/42, 0/28. That is 0 of 159 covered sites (sum computed here).
  - `variants.tsv` is empty, so precision is undefined.
  - Exact placements accepted: 2,722 to 9,254 per arm. Unplaced share: 0.10 % to 0.90 %.
  - Ortholog Spearman 0.9969 to 0.9997 [Tpa R T6].
- **Why recall is 0: the ADR-0079 floor.** With π̂ = 1 and c = 1, no q in a tolerant stratum can fall below 1/T. The figures below are the reverse-decoy strata from [Tpa M] `strains.*.T6_score.n2_mapping_fdr_all`:

| Arm | m = 1: targets / decoys / FDR floor | m = 2: targets / decoys / FDR floor | accepted (m ≥ 1) |
|---|---|---|---|
| Grady | 26 / 0 / 0.038 | 8 / 0 / 0.125 | 0 |
| Haiti B | 29 / 2 / 0.103 | 20 / 0 / 0.05 | 0 |
| Madras | 7 / 0 / 0.143 | 9 / 0 / 0.111 | 0 |
| Mexico A | 36 / 1 / 0.056 | 11 / 0 / 0.091 | 0 |
| Phil 1 | 51 / 1 / 0.039 | 18 / 0 / 0.056 | 0 |
| SS14 | 33 / 0 / 0.030 | 10 / 0 / 0.1 | 0 |

- No tolerant stratum reached 100 targets. Even strata with 0 decoys (Grady, Madras, SS14 at m = 1) accept nothing.
- This repeats Phase 1b N2: 0/11, with the m = 1 stratum at 13 targets [P1b "Why recall is zero"].
- The placement FDR check (T6) passes only on exact placements.

### 2.10 Determinism

- **TE7.** Two route-A reruns with the same inputs and seed: T. pallidum DAL-1 (22/22 `tables/` files byte-identical, job 2014083) and Synechococcus (22/22, job 2005461). Both pass [Tpa R T7; Syn R S13]. Group C made no TE7 rerun, and its key does not list one [Mtb R; Hvo R; Sisl R].
- **Convert.** `gippa convert` ran twice per route and gave byte-identical output in every dataset that records it (all ten reports' command logs).
- **Redaction.** `gippa report --redact` left `tables/` unchanged (sha256) wherever it was checked.

## 3. Answer-key premise errors and test-design limits (not GIPPA defects)

| Dataset, check | What happened | Source |
|---|---|---|
| C. difficile C6, PrdB C-terminus | The key assumed the 90-aa tail after the Sec codon is absent from route A. Bakta annotates most of it as FLEEGO_03445 (81 aa), which route A detects with 6 peptides (q 5.8e-4). Of route B's 8 peptides in the window, 6 are `cds_inframe` on FLEEGO_03445. 2 fall in an `nterm_extension` candidate lying 100 % inside RefSeq prdB. The rule fails 2/8, but GIPPA placed the evidence correctly | [Cdif R C6] |
| C. difficile S6, whole-replicon prophage | geNomad calls the whole 7,881-bp plasmid a prophage. A circular-shift null cannot vary when the region covers the entire replicon, so p = q = 1.0 even with 0/9 genes detected. The three chromosomal regions have p 0.026-0.059 and q 0.079 after BH over 4 | [Cdif R S6 detail] |
| V-type ATPase vs AtpA/AtpD | T. pallidum has only a V-type ATPase (atpA1/A2, atpD1/D2), so AtpA/AtpD is not testable. In C. difficile (V-type copies) and V. natriegens (a second atpD on chromosome 2), the name-based `te4` is ambiguous. The group B keys fixed F0F1 locus tags, so the locus-based verdict worked there | [Tpa R S3; Cdif R S3; Vnat R S3] |
| Archaea, TE3 and TE4 | Bakta names almost no ribosomal protein, so GIPPA's set is `not_evaluable`. The keys pre-registered RefSeq names carried by stop. TE4 is not testable by construction | [Hvo R row 2-3; Sisl R row 2-3] |
| S. islandicus row 8, pyrE/pyrF | pyrE (9 peptides) and pyrF (2) were detected. The key stated the empty-vector caveat in advance. A post hoc check found all 11 peptides also in the S. solfataricus PyrE/PyrF orthologues (vector marker). lacS has 0 peptides. The verdict stays fail against the key | [Sisl R row 8] |
| T. pallidum T4, TP0924a | Bakta types BIIDIP_00910 as a partial pseudogene, so it cannot be in the route-A database. Route B shows no peptide in the region. TP0865a has 0 peptides in either route | [Tpa R T4] |
| T. pallidum T3 extended list | The paper's list of 100 letter-suffix genes was not materialised before the search. Only T3's 12 genes (11 recovered) were compared | [Tpa R] |
| M. tuberculosis row 8, PE/PPE | Measured on 150 Bakta genes (RefSeq product by shared stop) rather than the key's 155 RefSeq CDS. The effect (12.7 % vs 74.5 %, Fisher p 1.1e-54) does not depend on this | [Mtb R row 8] |
| M. tuberculosis topology | NCBI labels NC_000962.3 linear. It was tagged circular and Bakta was rerun (lead decision). Row 15 confirms circular topology | [Mtb R] |
| Group A genome context | S6 was "reported" in group A and pass/fail in groups B/C. B. theta's 16 prophage regions (0 at q ≤ 0.05) were therefore not judged | [Bthe R S6] |
| Not runnable from the deposits | Synechococcus Y2 (WH8102 absent from the deposit), V. natriegens V4 (parental runs only), C. difficile C3 (0 spike-in peptides accepted at 0.5 pg/µg), C. difficile C4 (grdB at q 0.0145) | [Syn R; Vnat R; Cdif R] |
| Priors that did not hold (not errors) | M. xanthus M2 (BGC detection higher inside). C. crescentus C2 clause 1 (12/16 GlcNAc-cluster genes detected; passed on the low-abundance clause, median percentile 0.19) | [Mxan R M2; Ccre R C2] |
| T. pallidum preflight | `--strict` exited 3 on `genome.gaps` alone (N runs in the RefSeq assemblies). The waiver was recorded, and non-strict preflight exited 0 in all 20 arms | [Tpa R Problems] |

## 4. GIPPA defects found and fixed around Phase 2

| Issue → fix | Defect | Effect on campaign runs |
|---|---|---|
| #606 → PR #607 (ADR-0078, #605) | Picked gene FDR was anti-conservative on the peptide-vote route (`--search-db other`, or non-locus-tag IDs without protein hits). Decoys whose protein had no accepted target were left unpaired and dropped from D | Found in Phase 1b N2: 2,399 genes at max q 0.0051, about 0.026-0.028 realised FDR. Fixed: 2,305 genes at max q 0.0100, detection Jaccard with N1 1.000. E1/E2 keep the same genes; only q shifts slightly. Every Phase 2 run used the fixed code. On the vote route in Phase 2 are the T. pallidum T6 arms, e.g. Grady N2 662 vs N1 668 genes [P1b Rerun; P1 Recheck; Tpa M `T6_score`] |
| #611 → PR #612 | Promotech crashed on promoter windows containing N (gapped assemblies) | DAL-1 route B exited 4 at `bfc74ff6`. After the fix, 40-162 windows were skipped per genome [Tpa R gap-fix table] |
| #614 → PR #616 | Promoter Calculator crashed on non-ACGT bases in TU upstream sequences | DAL-1 A/B and Grady A exited 4 at `cff0fd68`. After the fix, 1 TU upstream region was skipped per genome [Tpa R] |
| #615 → PR #617 | TerminatorNet silently shifted coordinates after non-ACGT bases | Calls are now mapped back across 65-324 removed bases per genome. 13 TerminatorNet calls were dropped in DAL-1, 0 elsewhere [Tpa R gap-fix table] |
| ADR-0079 (#608) → PR #618 | No warning when a placement stratum's q floor c/T exceeds α, so zero tolerant placements looked like biology | A new run warning `placement.stratum_q_floor_above_alpha` (acceptance and q unchanged), accepted 2026-11-24 and merged at `77f942b9`. No Phase 2 run includes it (all ran earlier). By the ADR's rule it would fire on every T6 arm (§2.9 floors 0.030-0.143) and on same-genome runs with a small tolerant stratum (e.g. C. crescentus m = 1 with 1 target) |

The three gap fixes affect only genomes with non-ACGT bases. Among the Phase 2 genomes, that means only the T. pallidum assemblies (C. crescentus and B. theta have none; [Ccre R, Bthe R Inputs]).

## 5. Exploratory leads (uncalled; not discoveries)

All of these are below the class q floor. Most classes have an FDR of 0.07-1.0, so each list is expected to contain false hits. They are listed for Khoa's review only.

**RefSeq genes that Bakta did not call, with route-B peptides**

| Gene | Evidence | Source |
|---|---|---|
| rpmB (L28), T. pallidum | candidate in all 7 strains, 2 peptides, class FDR 0.5-1.0 | [Tpa R] |
| L36 Rv3461c, M. tuberculosis | 2 peptides | [Mtb R] |
| L34 Rv3924c, M. tuberculosis | 1 peptide | [Mtb R] |
| Rv3923c RNase P, M. tuberculosis | 4-peptide candidate | [Mtb R] |
| MazF7 and VapB43, M. tuberculosis | 7 and 9 peptides, `cds_outframe` | [Mtb R] |
| L37ae, S. islandicus | 3 peptides, `cds_outframe`, class FDR 0.143 | [Sisl R] |
| L40e, S. islandicus | 2 peptides, same coordinates as RefSeq, class FDR 0.167 | [Sisl R] |
| C. crescentus cdzC, cdzD, cdzI (bacteriocin) | 3, 3 and 2 peptides | [Ccre R C5] |
| H. volcanii µ-proteins ≤ 70 aa | 3 / 17 with peptides; 14 / 57 missed RefSeq CDS have ≥ 1 exact peptide | [Hvo R] |
| Synechococcus intergenic RefSeq genes | 2 genes, incl. a Nif11-like precursor | [Syn R] |

**Recurring candidates with no RefSeq match**
- T. pallidum, antisense near 70.2 kb: DAL-1, Madras, SS14; 2-3 peptides; class FDR 0.125-1.0 [Tpa R].
- T. pallidum, antisense near 443.3 kb: Haiti B, Mexico A; 4 peptides; best q 2.9e-6; class FDR 0.20-0.25 [Tpa R].

**Single-dataset candidates with no RefSeq match and ≥ 2 peptides**
- Synechococcus: 6, of which 4 are intergenic with RBS and promoter support; class FDR 0.14-0.22 [Syn R].
- C. crescentus: 2 (antisense at 1.94 Mb, class FDR 1.0; intergenic at 3.12 Mb, 0.17) [Ccre R].
- M. xanthus: 3 [Mxan R].
- B. theta, V. natriegens and C. difficile: none [Bthe R; Vnat R; Cdif R].

**N-terminal extensions and start differences**
- Most nterm_extension candidates match a RefSeq start upstream of Bakta's (§2.6).
- Recurring biology:
  - the upstream part of the prfB programmed frameshift, in B. theta (4 peptides), V. natriegens and E. coli;
  - Bakta's late infC start, in C. difficile (4 peptides) and E. coli;
  - the PrdB read-through past the Sec codon, in C. difficile.
  [Bthe R; Vnat R; Cdif R; P1 row 8]
- The largest: M. xanthus has 34 nterm_extension candidates with ≥ 2 peptides at class FDR 0.0127, the strongest with 15 peptides (BFAOEK_02911). M. tuberculosis has NOVEL_NC_000962.3_9899_p (11 peptides, class FDR 0.0118) [Mxan R; Mtb R].

**T. pallidum TP0134a** was recovered by one route-B peptide in frame, on an uncalled `cds_outframe` candidate (class FDR 0.20) [Tpa R T3].

## 6. Open items and caveats

- **Novel ORFs.** At α 0.01 and these genome sizes, route B called nothing in 11 genomes in a row, counting Phase 1. Any change to how classes are formed or pooled changes the statistics and needs an ADR and a Khoa decision under R3. Nothing was loosened here.
- **Cross-strain variants.** Recall at defaults is 0 (Phase 1b and T6). ADR-0079 makes this visible but does not fix it. Synechococcus Y2 and the V. natriegens fused-chromosome arm were not run.
- **RplL/RplJ cause is open.** A search that adds N-terminal acetylation would test the L7 explanation. That is a new run, not done here.
- **Unexplained warning.** Route A warns `novel_orfs.layer_missing` (`rbs_novel`) in the Phase 2 runs. Phase 1 E1 at the same commit did not. The cause was not established; no route-A check uses that layer [Hvo R; Sisl R Deviations; P1 Recheck].
- **M. tuberculosis wall time.** Route A's `gippa run` (5 h 17 min) was slower than route B's (3 h 12 min), on different nodes with overlapping runs. The cause was not measured [Mtb R].
- **Wall times are not comparable with Phase 1.** All Phase 2 runs used `farmshare_local` in one 16-CPU job; Phase 1's 1.1.0 runs used the distributed profile.
- **Group C procedure.** `gippa compare` and `report --redact` ran in the session allocation (job 2005379), not as separate `gippa-*` jobs, about 2-3 min per dataset [Mtb/Hvo/Sisl R].
- **Version string.** T. pallidum manifests read `gippa_version` 1.1.0 at `38e1cdba`; makedb at `bfc74ff6` reads 1.1.1 [Tpa R Inputs].
- **Sample identity is assumed** from file names in some deposits: V. natriegens parental in exponential phase, M. xanthus untreated wild type, Synechococcus N or P stress not labelled [Vnat R; Mxan R; Syn R].
- **Inconsistency between reports.** The Synechococcus report calls its RplL/RplJ pass the first in the campaign. C. crescentus and B. theta, run in parallel in group A, also passed [Syn R S3; Ccre R S3; Bthe R S3].
- **Error estimates.** TE4 ratios, T5 differences and island membership are unestimated. Wilson CIs treat peptides and sites as independent, so they are optimistic [P1b].
- **Memory headroom.** M. xanthus used 60.9 / 64 GB. A larger genome would need a bigger allocation [Mxan R M1].
- **Plan deliverable.** The handover plan's `SUMMARY.md` (one table across datasets plus the issues list) is not written yet. §1, §2.3, §2.4 and §4 here supply its contents.
