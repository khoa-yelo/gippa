# PXD072227 report: Synechococcus sp. CC9311 (campaign Phase 2, group B)

Campaign Phase 2, group B (Khoa, 2026-11-21: "run phase 2 with 1.1.1"). The answer key `ANSWER_KEY.md` was committed and pushed before any spectrum was searched, and it has not been edited. This report holds aggregates only: counts, shares, q-values, locus tags and 1-based coordinates. No peptide sequence and no DIA-NN row is in git.

- Machine-readable numbers: `metrics.json` (written by `validation/campaign/p2b/p2b_metrics.py`).
- Provenance: `provenance.json` (`p2b_provenance.py`).
- Every share carries a Wilson 95 % CI in brackets. "unestimated" means no error estimate is given.

## Verdicts

| # | Check | Verdict | Number (error estimate) | Cause / note |
|---|---|---|---|---|
| S1 | Mapping, TE1 | **pass** | 15,437 / 15,437 accepted route-A peptides exact in frame at their own gene: 100 % [99.975, 100] | 15 peptides with a non-ATG N-terminus, all mapped; 5 masked by the pilot rule. Placement FDR not applicable (`--search-db same`) |
| S2 | Expression, TE3 | **pass** | 46 / 47 ribosomal proteins detected, 97.9 % [88.9, 99.6]; median iBAQ log2 28.42 vs Q3 25.46 (percentile rank 0.95) | |
| S3 | Stoichiometry, TE4 | **pass** | log2(RplL/RplJ) = 2.53 (AAHINN_02716 / AAHINN_02717; top3 1.54), inside [1, 3]; log2(AtpA/AtpD) = 0.14 (AAHINN_02307 / AAHINN_02272; top3 0.58), inside [−1, 1]. Error unestimated | First pass of the RplL/RplJ pair in the campaign (it failed on Nissle, K-12 and the other group B datasets). The name-based `te4` agrees |
| S4 | Contaminants, TE5 | **pass** | 183 / 183 trypsin (6) and keratin (177) peptides in the contaminant bucket, 100 % [97.9, 100]; 0 placed on the genome | |
| S5 | Gene-level FDR | **pass** | decoy-estimated FDR of the detected set 0.0084 (A), 0.0054 (B); tier `decoy_picked` in both | |
| S6 | Genome enrichment | **not testable** | geNomad called no prophage region; `islands.tsv` empty | Reported: 3 antiSMASH regions, O/E 1.09 (q 0.54), 0.51 (q 0.060), 1.11 (q 0.54) |
| S7 | Operons | **pass** | 27 adjacent pairs among 53 ribosomal-protein/atp genes (14 OP, 12 NO, 1 SW); 0 with within-operon q ≤ 0.05 (min q 0.86); no `operon_flags` row on these pairs | Median abs. Δlog2 1.06. No operon flag anywhere in the genome |
| S8 | Novel ORFs vs RefSeq | **not testable** | 0 ORFs called (117 candidates) | Exploratory list below |
| S9 | Class q floor | **pass** | B warns `novel_orfs.class_q_floor_above_alpha`. Each class has fewer than 100 candidates with min class FDR > 0.01: nterm_extension 43 (0.029), antisense 37 (0.14), cds_outframe 19 (0.17), intergenic 18 (0.14). All are covered by the warning | Alpha 0.01 not loosened; hence 0 calls |
| S10 | Route B cost | **reported** | B keeps 1,837 / 1,897 A genes, 96.8 % [96.0, 97.5]; 20 genes only in B. Abundance Spearman 0.983 [0.976, 0.988] over 1,827 shared genes (1,383 TUs), Bland-Altman bias 0.08 log2 | Prior (kept ≥ 0.95) met. Cost table below |
| S11 | Second round | **not run** | 0 ORFs called | |
| S12 | Variants | **pass** | 0 rows in `variants.tsv` of A and B | |
| S13 | Determinism, TE7 | **pass** | Route-A rerun in the same lane job (2005461, seed 20261014): `tables/` byte-identical, 22 / 22 files (`diff -r`, empty diff) | `gippa report --redact` did not change A's `tables/` (sha256 checked) |
| Y1 | Photosynthesis core | **pass** | 7 / 7 markers detected (psaA, psaB, psbB, cbbL, cbbS, mpeA, mpeB; gene q 6.7e-4 each); median abundance log2 28.77 vs Q75 25.42 (percentile 0.966) | psaA, psaB and mpeB also carry shared peptides (`shared_evidence`) |
| Y2 | Cross-strain placement | **not run** | | As pre-registered: only CC9311 runs, WH8102 absent from the deposit |

Counts: **10 pass, 0 fail, 2 not testable (S6, S8), 2 not run (S11, Y2), 1 reported (S10).**

## Key numbers

| | Route A (genes) | Route B (six-frame) |
|---|---|---|
| DIA-NN precursors at 1 % (per run, 5 runs) | 15,914-17,535 | 14,830-16,450 |
| DIA-NN proteins identified (per run) | 906-966 | 872-958 |
| GIPPA genes detected (CDS + sORF, 2,859) | 1,897 (`decoy_picked`) | 1,857 (`decoy_picked`) |
| Gene-level FDR at q ≤ 0.01 | 0.0084 | 0.0054 |
| Detected with ≥ 2 peptides | 1,650 | 1,607 |
| Novel-ORF candidates / called | - | 117 / 0 |
| `gippa compare` enriched categories | 180 vs 177, Jaccard 0.785 | |

## Candidate ORFs for Khoa's review (exploratory, none called)

26 route-B candidates have ≥ 2 peptides: 11 nterm_extension, 9 intergenic, 4 antisense and 2 cds_outframe. 20 of them overlap a RefSeq CDS or pseudogene on the same strand by ≥ 50 %.

The nterm_extensions sit on known genes. Examples: infB (16 peptides), hisS, bchI, dnaA, radA, thrC, moaB and a BMC-domain protein (7 peptides). These suggest Bakta chose a later start than the protein uses.

Two intergenic candidates match RefSeq genes that Bakta did not call:
- 777,813-778,355 (+), 5 peptides, WP_041426407.1;
- 494,709-495,101 (+), a Nif11-like leader-peptide precursor.

The six candidates with no RefSeq match:

| ORF | class | coordinates (strand) | aa | peptides | best peptide q | class FDR | note |
|---|---|---|---|---|---|---|---|
| NOVEL_NC_008319.1_859779_p | intergenic | 859,779-859,886 (+) | 35 | 2 | 2.0e-3 | 0.22 | RBS, promoter supported |
| NOVEL_NC_008319.1_1135321_m | intergenic | 1,135,321-1,135,530 (−) | 69 | 2 | 1.2e-6 | 0.17 | RBS, promoter supported |
| NOVEL_NC_008319.1_2131444_p | intergenic | 2,131,444-2,131,620 (+) | 58 | 2 | 4.1e-4 | 0.17 | RBS, promoter supported |
| NOVEL_NC_008319.1_2552463_p | intergenic | 2,552,463-2,552,597 (+) | 44 | 2 | 1.5e-6 | 0.17 | RBS supported |
| NOVEL_NC_008319.1_1363427_p | antisense | 1,363,427-1,363,600 (+) | 57 | 2 | 1.0e-6 | 0.14 | antisense to AAHINN_01504 |
| NOVEL_NC_008319.1_973886_p | cds_outframe | 973,886-974,275 (+) | 129 | 2 | 1.0e-6 | 0.17 | 41 % overlap with RefSeq WP_167897175.1 |

The class FDRs are 0.14-0.22, so each class list is expected to contain false hits. These are leads, not discoveries.

## Route A vs route B

| | Route A | Route B | B / A |
|---|---|---|---|
| makedb records | 3,240 (incl. 381 CONT_) | 78,782 (+ 75,542 six-frame ORFs ≥ 20 aa) | 24x |
| FASTA bytes | 1,077,512 | 8,517,883 | 7.9x |
| predicted library | 424,472 precursors, 104,123,478 bytes | 2,372,086 precursors, 590,605,341 bytes | 5.6x |
| DIA-NN step 1 (16 CPUs) | 1 min 28 s | 7 min 09 s | 4.9x |
| DIA-NN step 2, 5 runs | 18 min 26 s | 1 h 08 min 09 s | 3.7x |
| `gippa run` (farmshare_local, 16 CPUs) | 1 h 28 min 14 s | 1 h 24 min 40 s | |
| genes detected | 1,897 | 1,857 | B keeps 96.8 % of A's |

## Command log (pinned commit `bfc74ff6`, GIPPA 1.1.1)

The whole recipe ran in one 16-CPU lane job: `sbatch -J gippa-camp-p2b-syn validation/campaign/p2b/lane.sbatch PXD072227 makedb lib diann_a diann_b conv_a conv_b run_a run_b post te7`. That is job **2005461** on wheat-07, 5 h 46 min 34 s, MaxRSS 49.1 GB, every step exit 0.

The code was the detached worktree `local/worktrees/rel111-pin`, used with `PYTHONPATH=<pin>/src` and the venv `bin/` on PATH. Every manifest records `git_dirty` false and `gippa_version` 1.1.1.

| step | command | wall |
|---|---|---|
| makedb | `gippa makedb --annotation J --genome G --outdir db_a`; the same with `--six-frame-orfs --min-orf-aa 20 --outdir db_b` | 9 s |
| lib | DIA-NN 2.3.1 step 1 (in-silico library) for A and B, `diann_campaign.cfg` | 95 s / 442 s |
| diann | DIA-NN step 2 on the 5 runs, A and B | 18 min 32 s / 1 h 08 min 18 s |
| conv | `gippa convert --from diann … --search-db same` (twice, byte-identical) + `gippa preflight --strict` (exit 0) | 15 s / 13 s |
| run | `gippa run --profile farmshare_local,apptainer … --threads 14 --memory 54.GB` (all modules, seed 20261014) | 1 h 28 min / 1 h 25 min |
| post | `gippa compare out_a out_b --labels genes,sixframe` (exit 0); `gippa report --redact` on A and B (exit 0) | 115 s |
| te7 | route-A rerun into `out_a_te7`, `diff -r` of `tables/` | 1 h 15 min 13 s |

Scoring was done afterwards:
- `p2b_metrics.py --acc PXD072227`, with `PYTHONPATH=<pin>:<pin>/src`;
- `p2b_provenance.py`.

## Deviations and notes

- **Route A warnings.** Route A warns `novel_orfs.decoys_not_makedb` and `novel_orfs.layer_missing` (`rbs_novel`). Route A has no six-frame entries, so `novel_orfs.tsv` is empty and the novel RBS layer gets no input. No check reads route A's novel ORFs.
- **Profile.** All Nextflow tasks ran inside the one 16-CPU allocation (`farmshare_local`), so the GIPPA wall times are not comparable with Phase 1's distributed runs.
- **Raw files.** The deposit has no checksum file. Byte sizes were checked against PRIDE and sha256 recorded (`PREP.md`). The N/P treatment of runs 23-27 is not labelled in the deposit.

## Provenance

| item | value |
|---|---|
| GIPPA | `bfc74ff6c20acc2cfcf298cdf23cbb30d8fae0dc` (1.1.1), clean |
| DIA-NN | 2.3.1 Academia, image sha256 `1946452d56d2…`; cfg sha256 `f89dfab40ef5…` |
| genome | `syn_cc9311.fna.gz` sha256 `675d870a5605…` |
| Bakta JSON | sha256 `700356c0f46b…` |
| DIA-NN reports A / B | `e981ae11df3e…` / `98e63e0848c9…` |
| raw | 5 runs, sha256 in `provenance.json` › `raw_sha256` |
| container digests, DB sizes, per-run DIA-NN counts | `provenance.json` |
