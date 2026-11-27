# PXD049476 report: Vibrio natriegens ATCC 14048 (campaign Phase 2, group B)

Campaign Phase 2, group B (Khoa, 2026-11-21: "run phase 2 with 1.1.1"). The answer key `ANSWER_KEY.md` was committed and pushed before any spectrum was searched, and it has not been edited. This report holds aggregates only: counts, shares, q-values, locus tags and 1-based coordinates. No peptide sequence and no DIA-NN row is in git.

- Machine-readable numbers: `metrics.json` (`validation/campaign/p2b/p2b_metrics.py`).
- Provenance: `provenance.json`.
- Shares carry a Wilson 95 % CI; differences of shares a Newcombe 95 % CI.

## Verdicts

| # | Check | Verdict | Number (error estimate) | Cause / note |
|---|---|---|---|---|
| S1 | Mapping, TE1 | **pass** | 22,438 / 22,438 accepted route-A peptides exact in frame, 100 % [99.983, 100] | 29 non-ATG N-terminal peptides, all mapped; 5 masked. Placement FDR not applicable (`same_db`) |
| S2 | Expression, TE3 | **pass** | 56 / 56 ribosomal proteins detected, 100 % [93.6, 100]; median iBAQ log2 32.65 vs Q3 26.12 (percentile rank 0.989) | |
| S3 | Stoichiometry, TE4 | **fail** (1 of 2 pairs) | log2(RplL/RplJ) = 0.96 (KGNFJJ_00023 / KGNFJJ_00022; top3 0.91), just outside [1, 3]. log2(AtpA/AtpD) = 0.19 (F0F1 copies KGNFJJ_02722 / KGNFJJ_02724), inside [−1, 1]. Error unestimated | Same pair that failed on Nissle and K-12 (a known risk in the key; waived for 1.1.0, ADR-0076). Here it misses by 0.04 log2. The name-based `te4` cannot evaluate AtpA/AtpD (`gene_ambiguous:atpD`; second copy KGNFJJ_04534 on chromosome 2), which is why the key fixed locus tags |
| S4 | Contaminants, TE5 | **pass** | 97 / 97 trypsin (4) and keratin (93) peptides in the contaminant bucket, 100 % [96.2, 100]; 0 on the genome | |
| S5 | Gene-level FDR | **pass** | 0.0099 (A), 0.0093 (B), `decoy_picked` | |
| S6 | Genome enrichment | **pass** | 2 geNomad prophage regions on chromosome 1, both under-detected with q ≤ 0.05: 936,639-972,714 (51 genes, O 5 vs E 26.1, O/E 0.19, q 0.015) and 1,495,937-1,532,177 (44 genes, O 6 vs E 24.3, O/E 0.25, q 0.015) | 8 antiSMASH regions, none with q ≤ 0.05 (q 0.81-0.99); no island |
| S7 | Operons | **pass** | 41 adjacent pairs among 65 ribosomal-protein/atp genes (36 OP, 5 NO); 0 with q ≤ 0.05 (min 0.80); no flag on these pairs | 16 `conflict_terminator_in_operon` flags elsewhere in the genome |
| S8 | Novel ORFs vs RefSeq | **not testable** | 0 called (78 candidates) | Exploratory below |
| S9 | Class q floor | **pass** | Warning present. All classes < 100 candidates and FDR > 0.01: antisense 34 (0.67), cds_outframe 17 (1.0), nterm_extension 15 (0.17), intergenic 12 (0.50) | |
| S10 | Route B cost | **reported** | B keeps 2,436 / 2,515 A genes, 96.9 % [96.1, 97.5]; 32 only in B. Spearman 0.994 [0.992, 0.996] over 2,436 genes (1,747 TUs), Bland-Altman bias 0.05 log2 | Prior met |
| S11 | Second round | **not run** | 0 called | |
| S12 | Variants | **pass** | 0 rows in A and B | |
| V1 | Detection per chromosome | **pass** | chr1 1,877 / 2,873 = 65.3 % [63.6, 67.1] vs chr2 638 / 1,678 = 38.0 % [35.7, 40.4]; difference +27.3 points [24.4, 30.2] | |
| V2 | Abundance per chromosome | **pass** | median abundance log2 of detected genes 24.48 (chr1, n 1,877) vs 22.50 (chr2, n 638); difference 1.98 [1.68, 2.31] (2,000 gene bootstraps, seed 20261014) | |
| V3 | Route B per chromosome | **reported** | kept chr1 1,829 / 1,877 = 97.4 % [96.6, 98.1]; chr2 607 / 638 = 95.1 % [93.2, 96.6] | |
| V4 | Fused-chromosome strain | **not run** | | As pre-registered: parental runs only |

Counts: **10 pass, 1 fail (S3), 1 not testable (S8), 2 not run (S11, V4), 2 reported (S10, V3).** S13 (TE7) is run on PXD072227 for this group.

## Key numbers

| | Route A | Route B |
|---|---|---|
| DIA-NN precursors at 1 % (per run, 4 runs) | 23,231-27,522 | 22,159-26,748 |
| DIA-NN proteins identified (per run) | 1,173-1,256 | 1,146-1,219 |
| GIPPA genes detected (CDS + sORF, 4,551; 1 partial excluded) | 2,515 | 2,468 |
| Gene-level FDR at q ≤ 0.01 | 0.0099 | 0.0093 |
| Detected with ≥ 2 peptides | 2,199 | 2,129 |
| Novel-ORF candidates / called | - | 78 / 0 |
| `gippa compare` enriched categories | 364 vs 374, Jaccard 0.907 | |

## Candidate ORFs (exploratory, none called)

7 route-B candidates have ≥ 2 peptides. All are `nterm_extension` and all lie ≥ 96 % inside a RefSeq gene:
- prfB (2,396,209-2,397,213 −, `no_start_codon`). This is the expected programmed-frameshift N-terminus, as in K-12.
- iscR, dusA, dusB, ilvY, rfbA, and a 4'-phosphopantetheinyl transferase (KGNFJJ_01839). These suggest a start upstream of Bakta's.

No candidate without a RefSeq match has ≥ 2 peptides.

## Route A vs route B

| | Route A | Route B | B / A |
|---|---|---|---|
| makedb records | 4,932 (incl. 381 CONT_) | 154,248 (+ 149,316 six-frame ORFs) | 31x |
| FASTA bytes | 1,890,882 | 15,657,076 | 8.3x |
| predicted library | 782,456 precursors, 191,949,806 bytes | 4,075,479 precursors, 1,018,493,240 bytes | 5.2x |
| DIA-NN step 1 | 2 min 02 s | 9 min 02 s | 4.4x |
| DIA-NN step 2, 4 runs | 17 min 33 s | 50 min 15 s | 2.9x |
| `gippa run` (farmshare_local, 16 CPUs) | 2 h 09 min 47 s | 2 h 10 min 46 s | |
| genes detected | 2,515 | 2,468 | B keeps 96.9 % |

## Command log (pinned commit `bfc74ff6`, GIPPA 1.1.1)

The whole recipe ran in one lane job: `sbatch -J gippa-camp-p2b-vnat validation/campaign/p2b/lane.sbatch PXD049476` (steps makedb … post). That is job **2008649** on wheat-05, 5 h 43 min 18 s, MaxRSS 51.7 GB, every step exit 0.

| step | wall |
|---|---|
| makedb A and B | 13 s |
| lib A / B | 129 s / 552 s |
| DIA-NN A / B | 17 min 39 s / 50 min 25 s |
| convert + preflight `--strict` (exit 0) A / B | 17 s / 13 s |
| `gippa run` A / B (`farmshare_local,apptainer`, all modules, seed 20261014, `--threads 14 --memory 54.GB`) | 2 h 09 min 47 s / 2 h 10 min 46 s |
| post: `gippa compare` (exit 0), `gippa report --redact` A and B (exit 0; A's `tables/` unchanged) | 149 s |

The commands are as in the PXD072227 report. Scoring: `p2b_metrics.py --acc PXD049476`, then `p2b_provenance.py`.

## Deviations and notes

- **Route A warnings.** `novel_orfs.decoys_not_makedb` and `layer_missing:rbs_novel`: route A has no candidates. No check depends on them.
- **Profile.** `farmshare_local` profile (one 16-CPU allocation).
- **Strain and phase.** The runs are the presumed parental strain in exponential phase (Khoa, question 18). The deposit has no sample key beyond the file names.
- **Raw files.** SHA-1 verified against the deposit (`PREP.md`).

## Provenance

| item | value |
|---|---|
| GIPPA | `bfc74ff6c20acc2cfcf298cdf23cbb30d8fae0dc` (1.1.1), clean |
| DIA-NN | 2.3.1 Academia, image sha256 `1946452d56d2…` |
| genome | `vnat_atcc14048.fna.gz` sha256 `abffc6c36484…` |
| Bakta JSON | sha256 `e96c6a7eacd3…` |
| DIA-NN reports A / B | `2ad6225bb9be…` / `e1508dd6a268…` |
| raw | 4 runs, SHA-1 verified; sha256 in `provenance.json` |
