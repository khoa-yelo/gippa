# PXD065346 report: Clostridioides difficile 630Δerm (campaign Phase 2, group B)

Campaign Phase 2, group B (Khoa, 2026-11-21: "run phase 2 with 1.1.1"). The answer key `ANSWER_KEY.md` was committed and pushed before any spectrum was searched, and it has not been edited. This report holds aggregates only. Numbers: `metrics.json`; provenance: `provenance.json`. Shares carry Wilson 95 % CIs.

The six H. volcanii spike-in proteins were searched as a separate `spike.faa` (`SPIKE_` records) next to the makedb database in both routes. They are excluded from every C. difficile claim.

## Verdicts

| # | Check | Verdict | Number (error estimate) | Cause / note |
|---|---|---|---|---|
| S1 | Mapping, TE1 | **pass** | 20,202 / 20,202 accepted route-A peptides exact in frame, 100 % [99.981, 100] | 21 non-ATG N-terminal peptides, all mapped; 3 masked. Placement FDR not applicable |
| S2 | Expression, TE3 | **pass** | 53 / 55 ribosomal proteins detected, 96.4 % [87.7, 99.0]; median iBAQ log2 30.75 vs Q3 24.64 (percentile rank 0.985) | |
| S3 | Stoichiometry, TE4 | **fail** (1 of 2 pairs) | log2(RplL/RplJ) = −1.22 (FLEEGO_00126 / FLEEGO_00125; top3 −0.98), outside [1, 3]. log2(AtpA/AtpD) = 0.12 (F0F1 FLEEGO_03723 / FLEEGO_03721), inside [−1, 1]. Error unestimated | Known-risk pair, failing as in Nissle, K-12, V. natriegens and M. xanthus. RplL has 16 unique peptides here, so this is not a few-peptide artefact. The name-based `te4` cannot evaluate Atp (V-type copies) |
| S4 | Contaminants, TE5 | **pass** | 248 / 248 trypsin (6) and keratin (242) peptides in the bucket, 100 % [98.5, 100]; 0 on the genome | |
| S5 | Gene-level FDR | **pass** | 0.0099 (A), 0.0097 (B), `decoy_picked` | |
| S6 | Genome enrichment | **fail** | 4 geNomad prophage regions, all observed < expected, none with q ≤ 0.05 | Direction holds in 4 / 4, significance in 0 / 4 (details below) |
| S7 | Operons | **pass** | 38 adjacent pairs among 65 ribosomal-protein/atp genes (37 OP, 1 NO); 0 with q ≤ 0.05 (min 0.91). One flag, atpE-atpF `conflict_terminator_in_operon`, has a concordant pattern (q 0.99) | Median abs. Δlog2 0.52; 36 terminator-conflict flags genome-wide |
| S8 | Novel ORFs vs RefSeq | **not testable** | 0 called (53 candidates) | |
| S9 | Class q floor | **pass** | Warning present. All classes < 100 with FDR > 0.01: antisense 17 (0.50), cds_outframe 11 (1.0), intergenic 13 (0.50), nterm_extension 12 (0.14) | |
| S10 | Route B cost | **reported** | B keeps 1,905 / 1,928 A genes, 98.8 % [98.2, 99.2]; 63 only in B. Spearman 0.991 [0.988, 0.993] over 1,883 genes (985 TUs), Bland-Altman bias −0.02 log2 | Prior met |
| S11 | Second round | **not run** | 0 called | |
| S12 | Variants | **pass** | 0 rows in A and B | |
| C1 | Small proteins (≤ 100 aa) | **reported** | A 94 / 419 = 22.4 % [18.7, 26.7]; B 97 / 419 = 23.2 % [19.4, 27.4] | Prior (below the paper's pooled 61.4 %) holds |
| C2 | Novel ORFs vs the paper | **not testable** | 0 called | |
| C3 | Spike-in proteins | **not testable** | 0 accepted `SPIKE_` peptides in A and B | `spike.faa` was loaded in both DIA-NN searches (step 1 and step 2 logs), but no spike-in precursor appears in either report (0 rows in `report.parquet`, `pr_matrix`, `pg_matrix`). This matches the prior that few or none are found at 0.5 pg/µg |
| C4 | Selenoprotein GrdB | **not testable** | grdB (FLEEGO_02514) not detected: 1 unique peptide, gene q 0.0145 | That one peptide lies on the N-terminal side of the Sec TGA; 0 across it |
| C5 | Selenoprotein GrdA | **reported** | 0 peptides in the N-terminal fragment (FLEEGO_02516), the 5-codon gap or the C-terminal fragment (FLEEGO_02515), in A and in B; no route-B candidate there | Neither fragment detected |
| C6 | PrdB C-terminus | **fail** | Route B places 8 peptides in 3,749,919-3,750,188 (−); only 2 belong to a candidate | The answer key's premise was wrong; not a GIPPA defect (below) |

Counts: **7 pass, 3 fail (S3, S6, C6), 4 not testable (S8, C2, C3, C4), 1 not run (S11), 3 reported (S10, C1, C5).** S13 (TE7) is run on PXD072227 for this group.

### S6 detail

The rule, fixed before reading: every prophage region has O < E and q ≤ 0.05.

| region | coordinates | genes | O | E | O/E | q |
|---|---|---|---|---|---|---|
| prophage_1 | NZ_CP016318.1 1,093,454-1,149,677 | 87 | 14 | 36.0 | 0.39 | 0.079 |
| prophage_2 | NZ_CP016318.1 1,579,890-1,601,066 | 24 | 1 | 11.3 | 0.09 | 0.079 |
| prophage_3 | NZ_CP016318.1 3,336,043-3,386,079 | 77 | 12 | 32.5 | 0.37 | 0.079 |
| prophage_4 | NZ_CP016319.1 1-7,881 (the whole plasmid) | 9 | 0 | 3.6 | 0.00 | 1.0 |

- **Chromosome regions:** the three have p 0.026-0.059 and q 0.079 after BH over 4 regions.
- **Plasmid:** geNomad calls the whole 7,881-bp plasmid a prophage. A circular-shift null on a replicon the region covers entirely cannot vary, so its p (and q) is 1.0 even with 0 of 9 genes detected. That entry cannot reach q ≤ 0.05, so the rule could not pass here. This is a limit of the test design for whole-replicon regions, not evidence against under-detection.
- No island; 5 antiSMASH regions (q ≥ 0.61).

### C6 detail (answer-key premise)

The key assumed the 90-aa PrdB tail after the Sec codon (3,749,919-3,750,188, −) is absent from route A. It is not:
- Bakta annotates most of it as its own CDS, FLEEGO_03445 (3,749,916-3,750,161, −, 81 aa, "D-proline reductase subunit gamma").
- Route A detects FLEEGO_03445 with 6 peptides (gene q 5.8e-4, coverage 0.91).

Route B's 8 peptides in the window split as follows:
- 6 are `cds_inframe` on FLEEGO_03445.
- 2 are `nterm_extension` peptides at 3,750,093-3,750,173. They extend FLEEGO_03445 upstream towards the Sec codon, in candidate NOVEL_NZ_CP016318.1_3749916_m (85 aa, `no_start_codon`, class FDR 0.14, not called).

That candidate lies 100 % inside RefSeq prdB (CDIF630erm_RS17620), which RefSeq annotates through the Sec to the downstream stop. So GIPPA placed the read-through evidence where it belongs. The pre-registered rule ("every route-B peptide in the window belongs to a candidate") still fails, 2 / 8, and the verdict is kept as fail.

## Key numbers

| | Route A | Route B |
|---|---|---|
| DIA-NN precursors at 1 % (per run, 3 runs) | 16,951-22,839 | 16,421-22,094 |
| DIA-NN proteins identified (per run) | 831-911 | 832-869 |
| GIPPA genes detected (CDS + sORF, 3,824) | 1,928 | 1,968 |
| Gene-level FDR at q ≤ 0.01 | 0.0099 | 0.0097 |
| Detected with ≥ 2 peptides | 1,727 | 1,716 |
| Novel-ORF candidates / called | - | 53 / 0 |
| `gippa compare` enriched categories | 229 vs 208, Jaccard 0.884 | |

## Candidate ORFs (exploratory, none called)

4 route-B candidates have ≥ 2 peptides. All are `nterm_extension` (class FDR 0.14) and all lie inside RefSeq genes:

| ORF | coordinates (strand) | aa | peptides | RefSeq gene | reading |
|---|---|---|---|---|---|
| NOVEL_NZ_CP016318.1_836517_p | 836,517-837,017 (+) | 166 | 4 | infC (IF-3), `no_start_codon` | Bakta starts IF-3 late, as in K-12 (infC uses the AUU start) |
| NOVEL_NZ_CP016318.1_2000110_p | 2,000,110-2,001,279 (+) | 389 | 2 | hydantoinase/oxoprolinase WP_262449342.1 | start upstream of Bakta's |
| NOVEL_NZ_CP016318.1_2411505_p | 2,411,505-2,412,539 (+) | 344 | 2 | NAD(P)H-dependent flavin oxidoreductase | start upstream of Bakta's |
| NOVEL_NZ_CP016318.1_3749916_m | 3,749,916-3,750,173 (−) | 85 | 2 | prdB, past the Sec | PrdB C-terminus read-through (C6) |

## Route A vs route B

| | Route A | Route B | B / A |
|---|---|---|---|
| makedb records | 4,205 (incl. 381 CONT_) | 100,852 (+ 96,647 six-frame ORFs) | 24x |
| FASTA bytes | 1,560,612 | 10,273,231 | 6.6x |
| spike-in FASTA (both routes) | 6 `SPIKE_` records | same | |
| predicted library | 754,760 precursors, 184,384,295 bytes | 2,333,967 precursors, 585,116,097 bytes | 3.1x |
| DIA-NN step 1 | 3 min 33 s | 7 min 52 s | 2.2x |
| DIA-NN step 2, 3 runs | 33 min 24 s | 1 h 22 min 19 s | 2.5x |
| `gippa run` (farmshare_local, 16 CPUs) | 3 h 04 min 37 s | 2 h 40 min 54 s | |
| genes detected | 1,928 | 1,968 | B keeps 98.8 % |

## Command log (pinned commit `bfc74ff6`, GIPPA 1.1.1)

The whole recipe ran in one lane job: `sbatch -J gippa-camp-p2b-cdif validation/campaign/p2b/lane.sbatch PXD065346` (makedb … post). That is job **2011418** on wheat-01, 7 h 57 min 54 s, MaxRSS 43.7 GB, every step exit 0. DIA-NN got `db_<r>/strain_db.faa` plus `spike/spike.faa` (`spike_fasta.py`) in both routes.

| step | wall |
|---|---|
| makedb A and B | 27 s |
| lib A / B | 225 s / 492 s |
| DIA-NN A / B | 33 min 32 s / 1 h 22 min 38 s |
| convert + preflight `--strict` (exit 0) A / B | 24 s / 17 s |
| `gippa run` A / B | 3 h 04 min 37 s / 2 h 40 min 54 s |
| post: compare (exit 0), `report --redact` A and B (exit 0; A's `tables/` unchanged) | 177 s |

Scoring was done on 2026-11-23: `p2b_metrics.py --acc PXD065346` (`PYTHONPATH=<pin>:<pin>/src`), then `p2b_provenance.py` (commit c5aa6a60).

## Deviations and notes

- **`p2b_metrics.py` was corrected for the C. difficile rows before this dataset was scored (commit c5aa6a60).**
  - `novel_overlapping` now takes a candidate's span from `_orf_span`. Candidates without a start codon (`start_longest` NA) used to collapse to their stop coordinate and drop out of window queries.
  - C6 is now scored as ANSWER_KEY.md words it: every route-B peptide in the window must belong to a candidate's peptide set, matched I/L-blind, because `novel_orfs.tsv` and `peptide_hits.tsv` may spell I/L differently. The earlier code only checked that no peptide had locus class `cds`/`annotated`. That would have passed C6, since the class here is `cds_inframe`.
  - No threshold changed. The other three datasets do not use these functions.
- **Route A warnings.** `novel_orfs.decoys_not_makedb` and `layer_missing:rbs_novel` (no candidates in route A).
- **Raw files.** The runs are the lowest spike-in level (0.5 pg/µg), technical triplicate (Khoa, question 20), SHA-1 verified (`PREP.md`).

## Provenance

| item | value |
|---|---|
| GIPPA | `bfc74ff6c20acc2cfcf298cdf23cbb30d8fae0dc` (1.1.1), clean |
| DIA-NN | 2.3.1 Academia, image sha256 `1946452d56d2…`; cfg sha256 `f89dfab40ef5…` |
| genome | `cdif_630derm.fna.gz` sha256 `49065d7c83c1…` |
| Bakta JSON | sha256 `7ae26fcc437b…` |
| spike FASTA | sha256 `793fcdf08f7a…` |
| DIA-NN reports A / B | `392ec4b77145…` / `59350567b806…` |
| peptides A | `1940ffa3c7d5…` |
| raw | 3 runs, SHA-1 verified; sha256 in `provenance.json` |
