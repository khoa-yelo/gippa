# PXD060688 report: Myxococcus xanthus DK1622 (campaign Phase 2, group B)

Campaign Phase 2, group B (Khoa, 2026-11-21: "run phase 2 with 1.1.1"). The answer key `ANSWER_KEY.md` was committed and pushed before any spectrum was searched, and it has not been edited. This report holds aggregates only. Numbers: `metrics.json`; provenance: `provenance.json`. Shares carry Wilson 95 % CIs; differences a Newcombe CI.

## Verdicts

| # | Check | Verdict | Number (error estimate) | Cause / note |
|---|---|---|---|---|
| S1 | Mapping, TE1 | **pass** | 54,278 / 54,278 accepted route-A peptides exact in frame, 100 % [99.993, 100] | 56 non-ATG N-terminal peptides, all mapped; 22 masked. Placement FDR not applicable |
| S2 | Expression, TE3 | **pass** | 55 / 61 ribosomal proteins detected, 90.2 % [80.2, 95.4]; median iBAQ log2 30.92 vs Q3 26.25 (percentile rank 0.986) | The point estimate clears the 90 % threshold, but the CI's lower bound does not. The rule judges the point estimate |
| S3 | Stoichiometry, TE4 | **fail** (1 of 2 pairs) | log2(RplL/RplJ) = 0.16 (BFAOEK_03063 / BFAOEK_03062; top3 −0.25), outside [1, 3]. log2(AtpA/AtpD) = 0.03 (BFAOEK_06952 / BFAOEK_06849), inside [−1, 1]. Error unestimated | Same known-risk pair as Nissle, K-12, V. natriegens and C. difficile. Cause not established (L7 N-acetylation not searched; few L12 peptides) |
| S4 | Contaminants, TE5 | **pass** | 277 / 277 trypsin (5) and keratin (272) peptides in the bucket, 100 % [98.6, 100]; 0 on the genome | |
| S5 | Gene-level FDR | **pass** | 0.0097 (A), 0.0091 (B), `decoy_picked` | |
| S6 | Genome enrichment | **pass** | 2 geNomad prophage regions, both under-detected with q ≤ 0.05: 1,407,056-1,441,584 (32 genes, O 3 vs E 20.3, O/E 0.15, q 0.008) and 2,151,086-2,218,747 (80 genes, O 7 vs E 47.7, O/E 0.15, q 0.004) | 27 antiSMASH and 5 CRISPR regions, none with q ≤ 0.05; no island |
| S7 | Operons | **pass** | 33 adjacent pairs among 66 ribosomal-protein/atp genes (14 OP, 19 NO); 0 with q ≤ 0.05 (min 0.89); no flag on these pairs | 3 `conflict_terminator_in_operon` flags elsewhere |
| S8 | Novel ORFs vs RefSeq | **not testable** | 0 called (265 candidates) | Exploratory below |
| S9 | Class q floor | **pass** | Warning present. Classes < 100 candidates with FDR > 0.01: antisense 96 (0.17), cds_outframe 28 (0.17), intergenic 26 (0.14), all covered | nterm_extension has 115 candidates (≥ 100, outside the rule) and min class FDR 0.0127. That is the closest any group B class came to alpha 0.01, so it has no call |
| S10 | Route B cost | **reported** | B keeps 4,791 / 4,864 A genes, 98.5 % [98.1, 98.8]; 150 only in B. Spearman 0.9925 [0.9908, 0.9939] over 4,787 genes (3,372 TUs), Bland-Altman bias −0.08 log2 | Prior met. B detects more genes than A (4,941 vs 4,864) |
| S11 | Second round | **not run** | 0 called | |
| S12 | Variants | **pass** | 0 rows in A and B | |
| M1 | Library size first | **reported** | Libraries predicted before any search: A 1,479,133 precursors / 359,929,707 bytes (step 1 4 min 34 s); B 9,802,626 precursors / 2,396,986,857 bytes (step 1 27 min 29 s), 6.6x. Route B is the largest six-frame database of the campaign (174,908 records, 27.1 MB). Step 1 succeeded in both | Per-step peak memory was not recorded. The whole lane peaked at MaxRSS 60.9 GB (sacct) within the 64,000 MB allocation, so there is little headroom for a larger genome |
| M2 | Detection inside antiSMASH regions | **reported** | pooled detected share inside 27 regions 573 / 803 = 71.4 % [68.1, 74.4] vs outside 4,291 / 6,515 = 65.9 % [64.7, 67.0]; difference +5.5 points [+2.1, +8.7]. 27 / 27 regions have detected genes; 0 regions with q ≤ 0.05 (min q 0.63) | The prior (lower inside) did not hold: biosynthetic clusters are, if anything, slightly better detected in these runs |
| M3 | antiSMASH ran | **pass** | `layers_missing` lists no antiSMASH/bgc layer in A or B; 27 bgc regions | |

Counts: **9 pass, 1 fail (S3), 1 not testable (S8), 1 not run (S11), 3 reported (S10, M1, M2).** S13 (TE7) is run on PXD072227 for this group.

## Key numbers

| | Route A | Route B |
|---|---|---|
| DIA-NN precursors at 1 % (per run, 4 runs) | 62,329-65,562 | 63,070-67,234 |
| DIA-NN proteins identified (per run) | 1,661-1,672 | 1,636-1,654 |
| GIPPA genes detected (CDS + sORF, 7,318) | 4,864 | 4,941 |
| Gene-level FDR at q ≤ 0.01 | 0.0097 | 0.0091 |
| Detected with ≥ 2 peptides | 4,471 | 4,480 |
| Novel-ORF candidates / called | - | 265 / 0 |
| `gippa compare` enriched categories | 312 vs 307, Jaccard 0.870 | |

## Candidate ORFs (exploratory, none called)

45 route-B candidates have ≥ 2 peptides: 34 nterm_extension, 5 antisense, 5 intergenic and 1 cds_outframe. 31 of them overlap a RefSeq CDS or pseudogene on the same strand by ≥ 50 %.

The nterm_extension candidates all have class FDR 0.0127. The strongest:

| ORF | coordinates (strand) | aa (longest) | peptides | gene | RefSeq overlap |
|---|---|---|---|---|---|
| NOVEL_NC_008095.1_3426703_p | 3,426,703-3,433,113 (+) | 2,136 | 15 | BFAOEK_02911 | 19 % (hypothetical WP_020480776.1) |
| NOVEL_NC_008095.1_3197964_m | 3,197,964-3,200,543 (−) | 859 | 12 | BFAOEK_02731 | 34 % (J-domain protein) |
| NOVEL_NC_008095.1_7170347_p | 7,170,347-7,171,534 (+) | 395 | 12 | BFAOEK_05744 | 76 % |
| NOVEL_NC_008095.1_1318343_m | 1,318,343-1,320,055 (−) | 570 | 9 | BFAOEK_01107 | 43 % (RDD family) |
| NOVEL_NC_008095.1_7205941_m | 7,205,941-7,208,304 (−) | 787 | 7 | BFAOEK_05766 | 100 % |

For an nterm_extension, the overlap is measured over the longest ORF, from its most upstream start. A low share therefore mostly reflects a far-upstream start codon, not peptides outside the RefSeq gene.

Candidates with no RefSeq match and ≥ 2 peptides:
- intergenic 7,770,733-7,771,260 (−), 175 aa, 2 peptides, class FDR 0.14;
- antisense 1,273,165-1,273,824 (+), `no_start_codon`, class FDR 0.48;
- nterm_extension of BFAOEK_06976 at 8,627,760-8,628,566 (−).

The 7-peptide intergenic candidate at 8,206,277-8,206,765 (−) overlaps an RHS-repeat protein by 31 %.

## Route A vs route B

| | Route A | Route B | B / A |
|---|---|---|---|
| makedb records | 7,699 (incl. 381 CONT_) | 174,908 (+ 167,209 six-frame ORFs) | 22.7x |
| FASTA bytes | 3,312,429 | 27,084,255 | 8.2x |
| predicted library | 1,479,133 precursors, 359,929,707 bytes | 9,802,626 precursors, 2,396,986,857 bytes | 6.6x |
| DIA-NN step 1 | 4 min 34 s | 27 min 29 s | 6.0x |
| DIA-NN step 2, 4 runs | 42 min 42 s | 3 h 37 min 36 s | 5.1x |
| `gippa run` (farmshare_local, 16 CPUs) | 8 h 09 min 52 s | 7 h 54 min 15 s | |
| genes detected | 4,864 | 4,941 | B keeps 98.5 % |

## Command log (pinned commit `bfc74ff6`, GIPPA 1.1.1)

The whole recipe ran in one lane job: `sbatch -J gippa-camp-p2b-mxan validation/campaign/p2b/lane.sbatch PXD060688` (makedb … post). That is job **2005509** on wheat-08, 21 h 04 min 53 s, MaxRSS 60.9 GB, every step exit 0.

| step | wall |
|---|---|
| makedb A and B | 20 s |
| lib A / B | 286 s / 1,704 s |
| DIA-NN A / B | 42 min 49 s / 3 h 37 min 53 s |
| convert + preflight `--strict` (exit 0) A / B | 34 s / 27 s |
| `gippa run` A / B | 8 h 09 min 52 s / 7 h 54 min 15 s |
| post: compare (exit 0), `report --redact` A and B (exit 0; A's `tables/` unchanged) | 322 s |

Scoring: `p2b_metrics.py --acc PXD060688`, then `p2b_provenance.py`.

## Deviations and notes

- **Report size.** Both runs warn `report.size_exceeded`. The HTML report of a 7,318-gene genome is over GIPPA's size budget; the tables are unaffected.
- **Route A warnings.** `novel_orfs.decoys_not_makedb` and `layer_missing:rbs_novel` (no candidates in route A).
- **Run times.** The longest lane of the group, since both `gippa run`s take about 8 h inside 16 CPUs (`farmshare_local`).
- **Raw files.** No file-to-sample key in the deposit: the runs are presumed wild type, untreated, exponential. SHA-1 verified (`PREP.md`).

## Provenance

| item | value |
|---|---|
| GIPPA | `bfc74ff6c20acc2cfcf298cdf23cbb30d8fae0dc` (1.1.1), clean |
| DIA-NN | 2.3.1 Academia, image sha256 `1946452d56d2…` |
| genome | `mxan_dk1622.fna.gz` sha256 `a36130a65c51…` |
| Bakta JSON | sha256 `af6d14c5415a…` |
| DIA-NN reports A / B | `0fc701fae11a…` / `cbdfba879598…` |
| raw | 4 runs, SHA-1 verified; sha256 in `provenance.json` |
