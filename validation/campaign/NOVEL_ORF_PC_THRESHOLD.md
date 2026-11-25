# Novel ORFs: a positive-control q threshold as an extra reporting tier

Khoa's request (2026-11-24): "Use the prfB gene and other positive control q value as the thresh to report on other runs for novel orfs."

**This is an exploratory reporting tier, not a call.**
- No GIPPA code or α changed (R3). GIPPA's calls at α 0.01 stand: 0 novel ORFs called in all 11 genomes.
- A candidate that passes τ here is "at or below the positive-control q". It is not "called".
- The file holds aggregates only: counts, q-values, ORF ids and 1-based coordinates. It has no peptide sequences and no DIA-NN rows.

**Code and data.**
- Script `novel_pc/novel_pc_threshold.py`; outputs `novel_pc/novel_pc_threshold.json` and `novel_pc/passers.tsv` (every passer of every tier, all columns below). Commit a7792873 on `validation/campaign-p2a`.
- Inputs: route-B `tables/novel_orfs.tsv` and `mapping/peptide_hits.tsv` under `$GIPPA_ROOT/campaign`, plus each genome's RefSeq GFF.
- Runs: E. coli `ecoli_k12/runs/out_e2_fix` (byte-identical novel_orfs.tsv to 1.1.0 E2); T. pallidum `<strain>/out_b_fix2`; the other nine datasets `out_b`.
- 1,581 candidates in total: 1,457 in Phase 2 and 124 in E. coli. Per-run counts equal PHASE2_SUMMARY §2.6.

## 1. Which q

- **Column `class_fdr`.** It is the only candidate-level q in `novel_orfs.tsv`.
- GIPPA's call rule (`gippa.novel.evidence._call`, rule 6) calls a candidate with a numeric `class_fdr` when `class_fdr <= alpha` **and** the class rule holds.
  - Class rule for nterm_extension: start codon, start upstream of Bakta's, RBS supported.
  - Class rule for intergenic: start codon, ≥ 10 codons, RBS supported.
  - Class rule for cds_outframe and antisense: ≥ 2 peptides.
  - A candidate with an NA `class_fdr` instead needs 2 peptides.
- `class_fdr` is a q-value from the class-wise target-decoy competition. It is monotone in the candidate's score within its class, so a q cut-off τ admits the class's best candidates up to that point.
- One candidate has NA `class_fdr`: T. pallidum DAL-1 `NOVEL_NC_016844.1_662748_p`, a pseudogene with 1 peptide, flagged `class_fdr_unestimated`. It passes no tier.
- **The tier uses q only.** The column "class rule" in the lists shows whether the visible part of GIPPA's class rule holds. `start_upstream_nt` is not in the table, so for nterm_extension this column checks only start codon + RBS.

## 2. Positive controls (defined before any passer was looked at)

**Rules** (the `CONTROLS` table in the script). A control is a route-B candidate that one of these known true genes independently confirms:

| Control | Rule | Reason |
|---|---|---|
| **prfB N-terminal extension** (primary; searched in every run) | candidate on prfB's strand whose ORF span overlaps RefSeq prfB, class nterm_extension | RF2's N-terminus is made by a programmed +1 frameshift that Bakta does not model. Peptides upstream of Bakta's start are true by biology |
| prfB leader (secondary) | the same overlap rule, any other class | the 25-codon leader before the frameshift |
| infC late start (E. coli, C. difficile) | nterm_extension overlapping RefSeq infC | IF3 starts at an AUU codon upstream of Bakta's start |
| PrdB read-through (C. difficile) | candidate overlapping RefSeq prdB | selenoprotein; translation continues through the UGA (Sec) codon |
| Bakta-missed RefSeq genes: T. pallidum rpmB (L28); M. tuberculosis L36 Rv3461c and L34 Rv3924c; S. islandicus L37ae and L40e; C. crescentus cdzC/D/I (CCNA_03933/03932/03931) | candidate span covers ≥ 50 % of the RefSeq CDS, same strand | RefSeq proteins (ribosomal proteins; bacteriocins with LC-MS / Ribo-seq evidence in RefSeq) that Bakta did not call |

**Notes on the matches.**
- The prdB rule matched two C. difficile candidates. `NOVEL_NZ_CP016318.1_3749916_m` ends at the prdB stop, so it is the read-through. `NOVEL_NZ_CP016318.1_3750189_m` is the N-terminal part of prdB, ending at the Sec codon. It is kept and labelled. Both have q 0.143, so τ is unaffected.
- The prfB rule also matched a Synechococcus prfB-leader candidate. No earlier report mentioned it.

**Controls found (23):**

| Dataset (run) | Control | Role | ORF id | Class | q (`class_fdr`) | Peptides | Start codon |
|---|---|---|---|---|---|---|---|
| E. coli | prfB | **primary** | NOVEL_U00096.3_3035184_m | nterm_extension | 0.143 | 4 | no |
| E. coli | prfB leader | secondary | NOVEL_U00096.3_3036205_m | intergenic | 0.143 | 1 | yes |
| E. coli | infC late start | secondary | NOVEL_U00096.3_1800096_m | nterm_extension | 0.143 | 2 | no |
| B. theta | prfB | **primary** | NOVEL_NC_004663.1_4594807_p | nterm_extension | 0.111 | 4 | no |
| V. natriegens | prfB | **primary** | NOVEL_NZ_CP009977.1_2396209_m | nterm_extension | 0.167 | 2 | no |
| Synechococcus | prfB leader | secondary | NOVEL_NC_008319.1_2330999_p | intergenic | 0.429 | 1 | yes |
| C. difficile | infC late start | secondary | NOVEL_NZ_CP016318.1_836517_p | nterm_extension | 0.143 | 4 | no |
| C. difficile | PrdB read-through | secondary | NOVEL_NZ_CP016318.1_3749916_m | nterm_extension | 0.143 | 2 | no |
| C. difficile | prdB N-terminal part (to the Sec) | secondary | NOVEL_NZ_CP016318.1_3750189_m | nterm_extension | 0.143 | 1 | no |
| C. crescentus | cdzC | secondary | NOVEL_NC_011916.1_774018_m | intergenic | 0.077 | 3 | yes |
| C. crescentus | cdzD | secondary | NOVEL_NC_011916.1_773733_m | intergenic | 0.077 | 3 | yes |
| C. crescentus | cdzI | secondary | NOVEL_NC_011916.1_773437_m | intergenic | 0.077 | 2 | yes |
| M. tuberculosis | L36 Rv3461c | secondary | NOVEL_NC_000962.3_3880286_m | intergenic | 0.100 | 2 | yes |
| M. tuberculosis | L34 Rv3924c | secondary | NOVEL_NC_000962.3_4410786_m | intergenic | 0.100 | 1 | yes |
| S. islandicus | L37ae | secondary | NOVEL_NC_017276.1_1215184_p | cds_outframe | 0.143 | 3 | yes |
| S. islandicus | L40e | secondary | NOVEL_NC_017276.1_1702242_p | intergenic | 0.167 | 2 | yes |
| T. pallidum DAL-1 | rpmB | secondary | NOVEL_NC_016844.1_385062_p | intergenic | 1.0 | 2 | yes |
| T. pallidum Grady | rpmB | secondary | NOVEL_NZ_CP035104.1_385043_p | intergenic | 1.0 | 2 | yes |
| T. pallidum Haiti B | rpmB | secondary | NOVEL_NZ_CP032623.1_385058_p | intergenic | 1.0 | 2 | yes |
| T. pallidum Madras | rpmB | secondary | NOVEL_NZ_CP078121.1_385073_p | intergenic | 1.0 | 2 | yes |
| T. pallidum Mexico A | rpmB | secondary | NOVEL_NC_018722.1_385032_p | intergenic | 1.0 | 2 | yes |
| T. pallidum Phil 1 | rpmB | secondary | NOVEL_NZ_CP035193.1_385044_p | intergenic | 1.0 | 2 | yes |
| T. pallidum SS14 | rpmB | secondary | NOVEL_NC_021508.1_385044_p | intergenic | 0.5 | 2 | yes |

**What the controls show.**
- prfB is not a route-B candidate in M. xanthus, C. difficile, T. pallidum, M. tuberculosis, H. volcanii, S. islandicus or C. crescentus. In Synechococcus only its leader is a candidate.
- All 3 primary controls, and both infC and both prdB controls, lack a start codon (`no_start_codon`). GIPPA's nterm_extension class rule therefore fails for them at any q. They show only where the q of true signal sits.

## 3. Thresholds

| Threshold | Value | Set by | Controls at or below it |
|---|---|---|---|
| **τ_prfB** = largest q of the primary prfB candidates | **0.1667** (1/6) | V. natriegens prfB | 15 / 23, 65.2 % [44.9, 81.2] |
| τ_all = largest q over all 23 controls | **1.0** | T. pallidum rpmB (6 strains) | 23 / 23 |
| τ per class: nterm_extension | 0.1667 | V. natriegens prfB | |
| τ per class: intergenic | 1.0 | T. pallidum rpmB | |
| τ per class: cds_outframe | 0.1429 | S. islandicus L37ae | |
| τ per class: antisense | none (no control) | | |
| *Sensitivity, post hoc:* τ_all without the T. pallidum rpmB controls | 0.4286 | Synechococcus prfB leader | 16 / 16 |

**τ_all is degenerate.**
- It admits every candidate with a numeric q: 1,580 / 1,581.
- In T. pallidum the classes hold 17-41 candidates, so q comes in coarse steps and rpmB lands at 0.5-1.0.
- **τ_prfB is the tier worth using.** The sensitivity row was chosen after seeing the q values and is labelled as such.

## 4. Results

**Definitions used in every table.**
- *Ctl*: passing controls. *Non-ctl*: passers that are not controls.
- *RefSeq match*: ≥ 50 % of the candidate (start_longest..end) lies on one RefSeq CDS or pseudogene on the same strand **and the same replicon**. As in the earlier reports, a candidate without a start codon counts as no match.
- *Span match*: the same rule, but a no-start candidate is scored by its ORF span (stop and ORF length).
- Shares carry Wilson 95 % CIs.
- *GIPPA class FDR at τ*: the expected number of false passers, Σ over classes of (passers in the class × the largest `class_fdr` among them), and that number as a share of the passers. Because `class_fdr` is a q-value, each class's passing set has an estimated FDR ≤ that class's largest q ≤ τ.

### 4.1 τ_prfB = 0.1667 (main tier)

| Dataset | Candidates | Pass (nterm / intergenic / cds_outframe / antisense) | Ctl | Non-ctl | Non-ctl ≥ 2 pep | Non-ctl RefSeq match | Span match | GIPPA class FDR at τ |
|---|---|---|---|---|---|---|---|---|
| E. coli (Phase 1 E2) | 124 | 14 (7 / 7 / 0 / 0) | 3 | 11 | 1 | 11/11, 100 % [74.1, 100] | 11/11 | 2.0 (14.3 %) |
| T. pallidum DAL-1 | 41 | 28 (14 / 4 / 2 / 8) | 0 | 28 | 5 | 17/28, 60.7 % [42.4, 76.4] | 19/28 | 3.3 (11.9 %) |
| T. pallidum Grady | 18 | 0 | 0 | 0 | 0 | — | — | 0 |
| T. pallidum Haiti B | 34 | 10 (10 / 0 / 0 / 0) | 0 | 10 | 5 | 9/10, 90.0 % [59.6, 98.2] | 10/10 | 1.0 (10 %) |
| T. pallidum Madras | 17 | 7 (0 / 2 / 3 / 2) | 0 | 7 | 2 | 6/7, 85.7 % [48.7, 97.4] | 6/7 | 1.0 (14.3 %) |
| T. pallidum Mexico A | 30 | 6 (6 / 0 / 0 / 0) | 0 | 6 | 0 | 6/6, 100 % [61.0, 100] | 6/6 | 1.0 (16.7 %) |
| T. pallidum Phil 1 | 35 | 6 (0 / 2 / 3 / 1) | 0 | 6 | 3 | 6/6, 100 % [61.0, 100] | 6/6 | 1.0 (16.7 %) |
| T. pallidum SS14 | 25 | 14 (6 / 2 / 4 / 2) | 0 | 14 | 3 | 11/14, 78.6 % [52.4, 92.4] | 13/14 | 2.0 (14.3 %) |
| *T. pallidum, 7 strains* | 200 | 71 | 0 | 71 | 18 | 55/71, 77.5 % [66.5, 85.6] | 60/71, 84.5 % [74.3, 91.1] | 9.3 (13.1 %) |
| C. crescentus | 164 | 67 (48 / 17 / 2 / 0) | 3 | 64 | 22 | 51/64, 79.7 % [68.3, 87.7] | 58/64 | 7.2 (10.7 %) |
| B. theta | 99 | 9 (9 / 0 / 0 / 0) | 1 | 8 | 3 | 6/8, 75.0 % [40.9, 92.9] | 8/8 | 1.0 (11.1 %) |
| Synechococcus | 117 | 62 (43 / 10 / 3 / 6) | 0 | 62 | 24 | 49/62, 79.0 % [67.4, 87.3] | 54/62 | 7.0 (11.3 %) |
| V. natriegens | 78 | 12 (12 / 0 / 0 / 0) | 1 | 11 | 6 | 10/11, 90.9 % [62.3, 98.4] | 11/11 | 2.0 (16.7 %) |
| M. xanthus | 265 | 116 (103 / 10 / 1 / 2) | 0 | 116 | 39 | 86/116, 74.1 % [65.5, 81.2] | 97/116 | 17.2 (14.8 %) |
| C. difficile | 53 | 7 (7 / 0 / 0 / 0) | 3 | 4 | 2 | 0/4, 0 % [0, 49.0] (all 4 lack a start codon) | 4/4 | 1.0 (14.3 %) |
| M. tuberculosis | 247 | 131 (109 / 15 / 6 / 1) | 2 | 129 | 42 | 105/129, 81.4 % [73.8, 87.2] | 119/129 | 15.6 (11.9 %) |
| H. volcanii | 142 | 51 (22 / 20 / 5 / 4) | 0 | 51 | 25 | 33/51, 64.7 % [51.0, 76.4] | 36/51 | 6.4 (12.6 %) |
| S. islandicus | 92 | 41 (28 / 9 / 4 / 0) | 2 | 39 | 13 | 17/39, 43.6 % [29.3, 59.0] (16 lack a start codon) | 32/39 | 6.1 (14.8 %) |
| **Total** | **1,581** | **581 (424 / 98 / 33 / 26)** | **15** | **566** | **195** | **423/566, 74.7 % [71.0, 78.1]** | **490/566, 86.6 % [83.5, 89.1]** | **74.8 (12.9 %)** |

**Non-control passers by class (all runs):**

| Class | Passers | RefSeq match | Span match | ≥ 2 peptides | Class rule met |
|---|---|---|---|---|---|
| nterm_extension | 417 | 337 | 399 | 127 | 325 |
| intergenic | 91 | 52 | 52 | 46 | 81 |
| cds_outframe | 32 | 23 | 28 | 11 | 11 |
| antisense | 26 | 11 | 11 | 11 | 11 |

**Subsets of the non-control passers:**
- With ≥ 2 peptides: 145 / 195 match RefSeq, 74.4 % [67.8, 80.0]; span match 166 / 195, 85.1 % [79.5, 89.4].
- Meeting the visible class rule: 428 passers, of which 373 match, 87.1 % [83.6, 90.0].
- cds_outframe and antisense passers with 1 peptide (36) need ≥ 2 peptides under GIPPA's class rule, so GIPPA would not call them at any q.

### 4.2 τ per class (nterm 0.1667, intergenic 1.0, cds_outframe 0.1429, antisense none)

| Dataset | Pass (nterm / intergenic / cds_outframe) | Ctl | Non-ctl | Non-ctl RefSeq match | GIPPA class FDR at τ |
|---|---|---|---|---|---|
| E. coli | 28 (7 / 21 / 0) | 3 | 25 | 12/25, 48.0 % [30.0, 66.5] | 22.0 |
| T. pallidum DAL-1 / Grady / Haiti B / Madras / Mexico A / Phil 1 / SS14 | 21 / 2 / 14 / 7 / 9 / 3 / 16 | 1 each | 20 / 1 / 13 / 6 / 8 / 2 / 15 | 14/20, 1/1, 10/13, 5/6, 8/8, 2/2, 11/15 | 8.0 / 2.0 / 5.0 / 4.4 / 4.0 / 3.0 / 4.5 |
| C. crescentus | 77 (48 / 29 / 0) | 3 | 74 | 52/74, 70.3 % [59.1, 79.5] | 33.0 |
| B. theta | 21 (9 / 12 / 0) | 1 | 20 | 8/20, 40.0 % [21.9, 61.3] | 13.0 |
| Synechococcus | 61 (43 / 18 / 0) | 1 | 60 | 45/60, 75.0 % [62.8, 84.2] | 18.1 |
| V. natriegens | 24 (12 / 12 / 0) | 1 | 23 | 12/23, 52.2 % [33.0, 70.8] | 14.0 |
| M. xanthus | 129 (103 / 26 / 0) | 0 | 129 | 86/129, 66.7 % [58.2, 74.2] | 41.0 |
| C. difficile | 20 (7 / 13 / 0) | 3 | 17 | 1/17, 5.9 % [1.0, 27.0] | 14.0 |
| M. tuberculosis | 140 (109 / 26 / 5) | 2 | 138 | 107/138, 77.5 % [69.9, 83.7] | 39.5 |
| H. volcanii | 57 (22 / 30 / 5) | 0 | 57 | 33/57, 57.9 % [45.0, 69.8] | 31.4 |
| S. islandicus | 50 (28 / 18 / 4) | 2 | 48 | 17/48, 35.4 % [23.4, 49.6] | 22.6 |
| **Total** | **679 (424 / 234 / 21)** | **23** | **656** | **424/656, 64.6 % [60.9, 68.2]** | **279.5 (41.2 %)** |

The intergenic τ of 1.0 admits every intergenic candidate. That is why the precision proxy drops and GIPPA's class FDR rises to 41 %.

### 4.3 τ_all = 1.0

- 1,580 / 1,581 candidates pass, i.e. every candidate with a numeric q. All 23 controls pass.
- Non-control RefSeq match: 597 / 1,557, 38.3 % [36.0, 40.8]; span match 815 / 1,557, 52.3 % [49.9, 54.8].
- GIPPA class FDR at τ: 1,196.6 expected false (75.7 %).
- This tier carries no information beyond "all candidates". Per-run rows are in `novel_pc_threshold.json`.

### 4.4 Sensitivity (post hoc): τ = 0.4286 (τ_all without T. pallidum rpmB)

- 789 pass, 16 of them controls.
- Non-control RefSeq match: 513 / 773, 66.4 % [63.0, 69.6]; span match 633 / 773, 81.9 % [79.0, 84.4].
- GIPPA class FDR at τ: 206.1 (26.1 %). Per-run rows are in the JSON.

### 4.5 Non-control τ_prfB passers with no RefSeq match even by span (76 of 566)

- The other 490 non-control passers overlap a RefSeq CDS or pseudogene. Most are nterm_extension at RefSeq starts upstream of Bakta's.
- Every passer, with its RefSeq gene and overlap fraction, is in `novel_pc/passers.tsv` (`threshold == tau_prfB`).
- Coordinates are 1-based closed. "~no start" means the coordinates are the ORF span estimated from the stop and the ORF length.
- "Needs ≥ 2 pep" = the class is cds_outframe or antisense. "Class rule" = the visible part of GIPPA's rule 4.
- No passers of this kind in E. coli, T. pallidum Grady / Haiti B / Mexico A / Phil 1, B. theta, V. natriegens or C. difficile.

| Run | ORF id | Class | Coordinates | aa | Peptides | q | RefSeq match | Needs ≥ 2 pep | Class rule |
|---|---|---|---|---|---|---|---|---|---|
| tpa_dal1 | NOVEL_NC_016844.1_150595_m | intergenic | 150,595-151,020 (-) | 141 | 1 | 0.125 | no | no | yes |
| tpa_dal1 | NOVEL_NC_016844.1_309659_m | antisense | 309,659-309,784 (-) | 41 | 1 | 0.125 | no | yes | no |
| tpa_dal1 | NOVEL_NC_016844.1_310151_m | antisense | 310,151-310,300 (-) ~no start | 49 | 1 | 0.125 | no | yes | no |
| tpa_dal1 | NOVEL_NC_016844.1_368101_p | antisense | 368,101-368,196 (+) | 31 | 1 | 0.125 | no | yes | no |
| tpa_dal1 | NOVEL_NC_016844.1_430089_p | antisense | 430,089-430,259 (+) | 56 | 1 | 0.125 | no | yes | no |
| tpa_dal1 | NOVEL_NC_016844.1_448858_m | intergenic | 448,858-449,037 (-) ~no start | 59 | 1 | 0.125 | no | no | no |
| tpa_dal1 | NOVEL_NC_016844.1_70185_p | antisense | 70,185-70,397 (+) | 70 | 2 | 0.125 | no | yes | yes |
| tpa_dal1 | NOVEL_NC_016844.1_915147_p | antisense | 915,147-915,314 (+) | 55 | 1 | 0.125 | no | yes | no |
| tpa_dal1 | NOVEL_NC_016844.1_443419_p | antisense | 443,419-443,937 (+) | 172 | 1 | 0.167 | no | yes | no |
| tpa_madras | NOVEL_NZ_CP078121.1_443431_p | antisense | 443,431-443,949 (+) | 172 | 1 | 0.143 | no | yes | no |
| tpa_ss14 | NOVEL_NC_021508.1_443398_p | antisense | 443,398-443,916 (+) | 172 | 1 | 0.125 | no | yes | no |
| ccre | NOVEL_NC_011916.1_1844245_m | nterm_extension | 1,844,245-1,845,300 (-) | 351 | 1 | 0.0256 | no | no | yes |
| ccre | NOVEL_NC_011916.1_519729_m | nterm_extension | 519,729-520,850 (-) | 373 | 1 | 0.0256 | no | no | yes |
| ccre | NOVEL_NC_011916.1_520158_p | intergenic | 520,158-520,391 (+) | 77 | 1 | 0.0769 | no | no | yes |
| ccre | NOVEL_NC_011916.1_977024_m | nterm_extension | 977,024-980,482 (-) | 1152 | 1 | 0.0833 | no | no | yes |
| ccre | NOVEL_NC_011916.1_2706705_m | intergenic | 2,706,705-2,706,860 (-) | 51 | 1 | 0.167 | no | no | yes |
| ccre | NOVEL_NC_011916.1_3118346_p | intergenic | 3,118,346-3,119,035 (+) | 229 | 3 | 0.167 | no | no | yes |
| syn | NOVEL_NC_008319.1_1096723_m | antisense | 1,096,723-1,096,791 (-) ~no start | 22 | 1 | 0.143 | no | yes | no |
| syn | NOVEL_NC_008319.1_1363427_p | antisense | 1,363,427-1,363,600 (+) | 57 | 2 | 0.143 | no | yes | yes |
| syn | NOVEL_NC_008319.1_1135321_m | intergenic | 1,135,321-1,135,530 (-) | 69 | 2 | 0.167 | no | no | yes |
| syn | NOVEL_NC_008319.1_1601494_p | intergenic | 1,601,494-1,601,721 (+) | 75 | 1 | 0.167 | no | no | yes |
| syn | NOVEL_NC_008319.1_1812320_m | intergenic | 1,812,320-1,812,400 (-) | 26 | 1 | 0.167 | no | no | yes |
| syn | NOVEL_NC_008319.1_2131444_p | intergenic | 2,131,444-2,131,620 (+) | 58 | 2 | 0.167 | no | no | yes |
| syn | NOVEL_NC_008319.1_2552463_p | intergenic | 2,552,463-2,552,597 (+) | 44 | 2 | 0.167 | no | no | yes |
| syn | NOVEL_NC_008319.1_973886_p | cds_outframe | 973,886-974,275 (+) | 129 | 2 | 0.167 | no | yes | yes |
| mxan | NOVEL_NC_008095.1_1000928_m | nterm_extension | 1,000,928-1,002,991 (-) | 687 | 1 | 0.0127 | no | no | yes |
| mxan | NOVEL_NC_008095.1_1318343_m | nterm_extension | 1,318,343-1,320,055 (-) | 570 | 9 | 0.0127 | no | no | yes |
| mxan | NOVEL_NC_008095.1_1753397_p | nterm_extension | 1,753,397-1,754,083 (+) | 228 | 3 | 0.0127 | no | no | yes |
| mxan | NOVEL_NC_008095.1_2493874_m | nterm_extension | 2,493,874-2,494,884 (-) | 336 | 2 | 0.0127 | no | no | yes |
| mxan | NOVEL_NC_008095.1_2758541_m | nterm_extension | 2,758,541-2,760,583 (-) | 680 | 2 | 0.0127 | no | no | yes |
| mxan | NOVEL_NC_008095.1_3197964_m | nterm_extension | 3,197,964-3,200,543 (-) | 859 | 12 | 0.0127 | no | no | yes |
| mxan | NOVEL_NC_008095.1_3283400_m | nterm_extension | 3,283,400-3,286,846 (-) | 1148 | 1 | 0.0127 | no | no | yes |
| mxan | NOVEL_NC_008095.1_3426703_p | nterm_extension | 3,426,703-3,433,113 (+) | 2136 | 15 | 0.0127 | no | no | yes |
| mxan | NOVEL_NC_008095.1_5859816_p | nterm_extension | 5,859,816-5,862,137 (+) | 773 | 2 | 0.0127 | no | no | yes |
| mxan | NOVEL_NC_008095.1_7192637_p | nterm_extension | 7,192,637-7,196,518 (+) | 1293 | 2 | 0.0127 | no | no | no |
| mxan | NOVEL_NC_008095.1_8627760_m | nterm_extension | 8,627,760-8,628,566 (-) | 268 | 2 | 0.0127 | no | no | yes |
| mxan | NOVEL_NC_008095.1_3204813_p | intergenic | 3,204,813-3,204,992 (+) | 59 | 1 | 0.143 | no | no | yes |
| mxan | NOVEL_NC_008095.1_4120117_p | intergenic | 4,120,117-4,120,224 (+) ~no start | 35 | 1 | 0.143 | no | no | no |
| mxan | NOVEL_NC_008095.1_6155695_m | intergenic | 6,155,695-6,155,763 (-) ~no start | 22 | 1 | 0.143 | no | no | no |
| mxan | NOVEL_NC_008095.1_7591920_m | intergenic | 7,591,920-7,592,099 (-) | 59 | 1 | 0.143 | no | no | no |
| mxan | NOVEL_NC_008095.1_7770733_m | intergenic | 7,770,733-7,771,260 (-) | 175 | 2 | 0.143 | no | no | yes |
| mxan | NOVEL_NC_008095.1_6895378_m | nterm_extension | 6,895,378-6,896,796 (-) | 472 | 1 | 0.146 | no | no | yes |
| mxan | NOVEL_NC_008095.1_410191_m | antisense | 410,191-412,728 (-) | 845 | 3 | 0.167 | no | yes | yes |
| mxan | NOVEL_NC_008095.1_8206277_m | intergenic | 8,206,277-8,206,765 (-) | 162 | 7 | 0.167 | no | no | yes |
| mtb | NOVEL_NC_000962.3_3429969_p | nterm_extension | 3,429,969-3,431,315 (+) | 448 | 1 | 0.0118 | no | no | yes |
| mtb | NOVEL_NC_000962.3_2181262_m | intergenic | 2,181,262-2,182,560 (-) | 432 | 5 | 0.0909 | no | no | yes |
| mtb | NOVEL_NC_000962.3_2558671_m | antisense | 2,558,671-2,558,865 (-) | 64 | 1 | 0.0909 | no | yes | no |
| mtb | NOVEL_NC_000962.3_3182268_p | cds_outframe | 3,182,268-3,183,395 (+) | 375 | 9 | 0.0909 | no | yes | yes |
| mtb | NOVEL_NC_000962.3_3414936_m | intergenic | 3,414,936-3,415,280 (-) | 114 | 1 | 0.0909 | no | no | yes |
| mtb | NOVEL_NC_000962.3_3580526_p | intergenic | 3,580,526-3,580,684 (+) | 52 | 3 | 0.0909 | no | no | yes |
| mtb | NOVEL_NC_000962.3_1697604_m | intergenic | 1,697,604-1,697,684 (-) | 26 | 1 | 0.1 | no | no | yes |
| mtb | NOVEL_NC_000962.3_2124830_m | intergenic | 2,124,830-2,125,048 (-) | 72 | 2 | 0.1 | no | no | yes |
| mtb | NOVEL_NC_000962.3_2684357_p | intergenic | 2,684,357-2,684,596 (+) | 79 | 1 | 0.1 | no | no | yes |
| mtb | NOVEL_NC_000962.3_2013748_p | nterm_extension | 2,013,748-2,015,262 (+) | 504 | 1 | 0.119 | no | no | yes |
| hvo | NOVEL_NC_013966.1_326705_p | cds_outframe | 326,705-327,742 (+) | 345 | 1 | 0.0833 | no | yes | no |
| hvo | NOVEL_NC_013966.1_36821_m | intergenic | 36,821-37,267 (-) | 148 | 6 | 0.0833 | no | no | yes |
| hvo | NOVEL_NC_013967.1_1911333_m | intergenic | 1,911,333-1,911,362 (-) ~no start | 9 | 1 | 0.0833 | no | no | no |
| hvo | NOVEL_NC_013967.1_339300_m | intergenic | 339,300-339,581 (-) | 93 | 3 | 0.0833 | no | no | yes |
| hvo | NOVEL_NC_013964.1_353800_p | intergenic | 353,800-354,093 (+) ~no start | 97 | 1 | 0.0909 | no | no | no |
| hvo | NOVEL_NC_013966.1_129718_p | intergenic | 129,718-130,038 (+) | 106 | 1 | 0.0909 | no | no | yes |
| hvo | NOVEL_NC_013966.1_245790_p | intergenic | 245,790-246,359 (+) | 189 | 6 | 0.0909 | no | no | yes |
| hvo | NOVEL_NC_013966.1_361252_m | intergenic | 361,252-361,485 (-) | 77 | 1 | 0.0909 | no | no | yes |
| hvo | NOVEL_NC_013966.1_438089_p | intergenic | 438,089-438,697 (+) | 202 | 1 | 0.0909 | no | no | yes |
| hvo | NOVEL_NC_013966.1_634653_m | nterm_extension | 634,653-635,249 (-) | 198 | 1 | 0.0909 | no | no | yes |
| hvo | NOVEL_NC_013967.1_1693214_m | intergenic | 1,693,214-1,693,417 (-) | 67 | 1 | 0.0909 | no | no | yes |
| hvo | NOVEL_NC_013967.1_1693430_m | intergenic | 1,693,430-1,693,489 (-) ~no start | 19 | 1 | 0.0909 | no | no | no |
| hvo | NOVEL_NC_013966.1_243685_p | antisense | 243,685-244,197 (+) | 170 | 2 | 0.167 | no | yes | yes |
| hvo | NOVEL_NC_013967.1_2099838_p | intergenic | 2,099,838-2,101,568 (+) | 576 | 4 | 0.167 | no | no | yes |
| hvo | NOVEL_NC_013967.1_2527630_m | antisense | 2,527,630-2,528,166 (-) | 178 | 1 | 0.167 | no | yes | no |
| sisl | NOVEL_NC_017276.1_1148017_p | cds_outframe | 1,148,017-1,148,364 (+) | 115 | 2 | 0.143 | no | yes | yes |
| sisl | NOVEL_NC_017276.1_1627723_m | intergenic | 1,627,723-1,627,941 (-) | 72 | 1 | 0.143 | no | no | yes |
| sisl | NOVEL_NC_017276.1_1243417_m | intergenic | 1,243,417-1,243,512 (-) | 31 | 1 | 0.167 | no | no | yes |
| sisl | NOVEL_NC_017276.1_1261419_m | intergenic | 1,261,419-1,261,487 (-) | 22 | 1 | 0.167 | no | no | yes |
| sisl | NOVEL_NC_017276.1_1655651_m | intergenic | 1,655,651-1,655,725 (-) | 24 | 1 | 0.167 | no | no | yes |
| sisl | NOVEL_NC_017276.1_1786001_p | intergenic | 1,786,001-1,786,051 (+) ~no start | 16 | 1 | 0.167 | no | no | no |
| sisl | NOVEL_NC_017276.1_2193084_p | intergenic | 2,193,084-2,193,170 (+) | 28 | 1 | 0.167 | no | no | yes |

**What recurs or stands out.**
- The T. pallidum antisense candidate near 443.4 kb (172 aa, 1 peptide here) passes in DAL-1, Madras and SS14.
- M. xanthus nterm_extension candidates with 9, 12 and 15 peptides have no RefSeq match. They sit in the class with q 0.0127.
- Among non-RefSeq passers with ≥ 5 peptides:
  - M. tuberculosis `3182268_p` (cds_outframe, 9 peptides) and `2181262_m` (intergenic, 5 peptides).
  - H. volcanii `36821_m` and `245790_p` (intergenic, 6 peptides each).
  - M. xanthus `8206277_m` (intergenic, 7 peptides).

## 5. Caveats

- **The threshold is chosen from the same data, after the fact.** The controls are candidates in these same runs, and τ is their largest q. There is no independent calibration set.
- **It is not an FDR guarantee.** τ_prfB = 0.1667 is about 17 × α.
- **GIPPA's own estimate of the passing set's class FDR.** At τ_prfB, GIPPA's class FDRs imply about 74.8 false of 581 passers (12.9 %), and per dataset 10.7-16.7 % (§4.1). That estimate assumes GIPPA's class-wise target-decoy model holds in these small classes. It is unchecked here, and the q-values are coarse (steps of 1/6, 1/7, … in classes of 10-60 candidates).
- **RefSeq match is a weak precision proxy, not ground truth.**
  - A RefSeq match confirms that a protein exists there, not that the candidate's start or extension is right.
  - A non-match is not proof of a false hit.
  - The no-start rule deflates the strict share (C. difficile 0/4 strict, 4/4 by span).
  - Peptides are treated as independent in the Wilson CIs, so the CIs are optimistic.
- **The controls are biased toward easy, abundant cases.** prfB, infC and ribosomal proteins are highly expressed, and 15 of 23 controls are nterm_extension or intergenic. True novel ORFs that are short, rare or antisense would sit at higher q. antisense has no control at all.
- **τ_all is driven by T. pallidum's tiny classes.** It equals 1.0 and admits everything, so it is not a usable tier. The per-class intergenic τ (1.0) has the same problem.
- **The primary controls would not be called at any q.** All 3 prfB nterm_extension controls lack a start codon, so GIPPA's nterm_extension class rule rejects them.
- **None of this changes a GIPPA call.** GIPPA called 0 novel ORFs at α 0.01 in all 11 genomes, and α was not loosened (R3). Any change to classes, pooling or α needs an ADR and Khoa's decision.

## 6. Side finding: H. volcanii RefSeq-match share in the existing report

- `PXD066008/REPORT.md` reports 74 / 142 candidates matching RefSeq. Its matching code does not check the replicon, and H. volcanii has 5 replicons.
- Requiring the same replicon gives 50 / 142: nterm_extension 28 / 41, intergenic 11 / 30, cds_outframe 9 / 16, antisense 2 / 55 (the report had 13 / 55).
- PHASE2_SUMMARY §2.6 copies the 74 / 142.
- Single-replicon datasets are unaffected: E. coli 23, T. pallidum 103, C. crescentus 65, B. theta 22, M. tuberculosis 127 and S. islandicus 25 all reproduce exactly. The V. natriegens report gave only its ≥ 2-peptide share, which was not rechecked.
- This needs a correction in the H. volcanii report and the summary. I did not edit those files.
