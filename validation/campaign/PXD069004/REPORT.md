# PXD069004 report: Caulobacter crescentus NA1000 (Phase 2, group A, dataset 3)

Checks pre-registered in `ANSWER_KEY.md` (commit `854c143c`, before any search). Aggregates only:
counts, shares, q-values, locus tags and 1-based coordinates. No peptide sequence and no DIA-NN row
is in git. Machine-readable numbers: `metrics.json` (`p2a/p2_metrics.py`). Per-route provenance:
`provenance/{a,b}.json` (`p2a/provenance.py`). Shares carry a Wilson 95 % CI in brackets;
"unestimated" means no error estimate exists. `$D` = `$GIPPA_ROOT/campaign/PXD069004`.

## Inputs

| Item | Value |
|---|---|
| Accession, runs | PXD069004 (eLife 2026, PMC13046382), `AmiX-5-4/5/6.raw`: wild type in M2G without GlcNAc, Orbitrap Exploris 480 DIA, SHA-1 verified (`PREP.md`) |
| Genome | GCF_000022005.1, NC_011916.1, 4,042,929 bp, circular, GC 67.2 %, no N or IUPAC bases |
| Annotation | Bakta 1.12.2, DB 6.0 full (job 2001987): 3,722 CDS + 2 sORF (3,723 coding rows in `genes.tsv`) vs RefSeq 3,886 CDS (ratio 0.958) |
| GIPPA | `bfc74ff6c20acc2cfcf298cdf23cbb30d8fae0dc` (1.1.1 release code, `git_dirty` false; provenance says 1.1.1, not yet published). Not rerun at `cff0fd68` (#612 Promotech gap fix): the genome has no non-ACGT base, so the fix cannot change the outputs |
| DIA-NN | 2.3.1 Academia, image sha256 `1946452d56d2fb0f66e05ced473df88e53484e28040f1ecbd83c4a43405a6381`, flags of `validation/pilot/diann_A2.cfg` (`p2a/diann_campaign.cfg`), `.raw` read directly |
| Containers | digests per tool in `provenance/{a,b}.json` › `gippa_run.tools` (antiSMASH 8.0.4, geNomad 1.12.0, TerminatorNet sharded, TransTermHP 2.09, Promotech, Promoter Calculator 1.2.4, OSTIR 1.1.3, UniOp) |

## Verdicts

| # | Check | Verdict | Number (error estimate) | Note |
|---|---|---|---|---|
| S1 | Mapping, TE1 (≥ 99.5 % exact in frame) | **pass** | 32,356 / 32,356 route-A accepted peptides exact in frame at their own gene: 100 % [99.988, 100] | 34 with a non-ATG N-terminus, all mapped. Placement FDR not applicable (`same_db`); `mapping_fdr.tsv` accepts no tolerant placement (1-difference stratum: 1 target / 28 decoys) |
| S2 | Expression, TE3 | **pass** | 54 / 54 ribosomal proteins detected, 100 % [93.4, 100]; median iBAQ log2 30.09 vs Q3 25.16 (percentile rank 0.985) | |
| S3 | Stoichiometry, TE4 | **pass** | log2(RplL/RplJ) = 2.68 (top3 1.68), in [1, 3]; log2(AtpA/AtpD) = 0.89 (top3 0.26), in [−1, 1]. Error unestimated | First campaign dataset where RplL/RplJ passes (it failed in the pilot and in Phase 1) |
| S4 | Contaminants, TE5 | **pass** | 106 trypsin (4) and keratin (102) peptides: 105 / 105 non-genomic ones in the contaminant bucket, 100 % [96.5, 100]; 1 matches the genome exactly (`genome_exact`, not a violation under the pilot rule) | Unmapped: 126 peptides, all labelled contaminant (84 `no_hit`, 41 `too_short`, 1 `rejected_placement_fdr`) |
| S5 | Operons | **pass** | 35 adjacent ribosomal/atp pairs (29 OP, 6 NO), all both detected; 0 with within-operon q ≤ 0.05 (min q 0.65); 0 `operon_flags` rows | Median abs Δlog2 0.71 |
| S6 | Genome context | **reported** | 12 region tests (10 antiSMASH BGC, 2 AMR). 1 at q ≤ 0.05: `bgc_3`, 2,943,501-2,965,908, 20 genes, O/E 0.33 (5 vs 15.0 expected), q 0.032 (route B q 0.047). Islands (route A): ISL_1 high, 1,430,061-1,445,888 (28 genes, all detected); ISL_2 low, 2,010,295-2,039,776 (23 genes, 26 % detected). Island q unestimated (`q_value` empty). 0 operon flags | Route B keeps ISL_2 only (`gippa compare` islands Jaccard 0.45) |
| S7 | Novel ORFs vs RefSeq | **not testable** | 0 ORFs called (164 candidates) | Exploratory, all candidates: 65 / 164 match a RefSeq CDS (39.6 % [32.5, 47.3]). By class: nterm_extension 39/49, intergenic 16/29, cds_outframe 9/21, antisense 1/65. Candidates with ≥ 2 peptides: 22 / 26 match (84.6 % [66.5, 93.8]); table below |
| S8 | Class q floor | **pass** | warning `novel_orfs.class_q_floor_above_alpha`; classes of 21-65 candidates; minimum class FDR 0.026 (nterm_extension), 0.077 (intergenic), 0.17 (cds_outframe), 0.76 (antisense); none ≤ 0.01, none called | Alpha 0.01 not loosened. nterm_extension (49 candidates) came closest |
| S9 | Route B cost | **reported** | route B keeps 2,917 / 3,005 route-A genes, 97.1 % [96.4, 97.6]; 21 genes detected only in route B. Abundance Spearman 0.991 [0.987, 0.994] over 2,900 shared genes (2,083 TUs, TU bootstrap); Bland-Altman bias 0.026 log2 | Cost table below |
| S10 | Gene-level FDR | **reported** | route A: 3,005 detected (2,722 with ≥ 2 peptides), tier `decoy_picked`, decoy FDR at q ≤ 0.01 0.0093. Route B: 2,938 (2,687), FDR 0.0037 | |
| S11 | Second round | **not run** | 0 ORFs called | |
| C1 | PG recycling enzymes in wild type (≥ 6 of 10) | **pass** | 10 / 10 detected at gene q ≤ 0.01, 100 % [72.2, 100]: GlmM, AmpG, AnmK, NagZ, GlmU, MupP, AmpD, AmgK, MurU, NagK. Abundance percentiles 0.25 (AmpG) to 0.82 (GlmM) | |
| C2 | GlcNAc clusters uninduced (≤ 4 of 16 detected, or median percentile of the detected ≤ 0.50) | **pass** (second clause) | 12 / 16 detected (75 % [50.5, 89.8]), so clause 1 fails. The median abundance percentile of the 12 is 0.19 (cluster 1: 0.29; cluster 2: 0.09). Both nagE PTS genes (CCNA_00457, CCNA_00571) carry only peptides shared with each other (23 each) and are not detected | Low but present without GlcNAc. CCNA_00458 (ptsG2) is at percentile 0.03 |
| C3 | Depth | **reported** | route A 3,005 / 3,723 coding genes, 80.7 % [79.4, 81.9]; prior ≥ 0.50 met | |
| C4 | High-GC genome | **pass** | both runs exit 0; TerminatorNet ran sharded (`ANNOTATE:TERMINATORNET_SHARDED`). `layers_missing`: route A `rbs_novel` only (no novel ORFs on route A), route B none | |
| C5 | Annotation vs RefSeq | **reported** | 16 of the 29 `intergenic` candidates lie on RefSeq CDSs that Bakta did not call, e.g. the bacteriocin genes cdzC (CCNA_03933), cdzD (CCNA_03932) and cdzI (CCNA_03931), with 3, 3 and 2 peptides | Known in RefSeq; not new |

Counts: **10 pass, 0 fail, 2 not testable or not run (S7, S11), 5 reported** (S6, S9, S10, C3, C5).

## Candidates with ≥ 2 peptides (route B; none called)

The class FDR is the minimum of the candidate's class. Coordinates are 1-based, from the longest-ORF
start to the stop. "RefSeq" is the same-strand RefSeq CDS that overlaps most of the ORF
(answer-key annotation).

| ORF | class | coordinates (strand) | aa | peptides | best q | class FDR | RefSeq (share of ORF) | Reading |
|---|---|---|---|---|---|---|---|---|
| NOVEL_NC_011916.1_271910_m | nterm_extension | 271,910-272,620 (−) | 236 | 2 | 6.4e-4 | 0.026 | CCNA_00259 (1.0) | known start difference |
| NOVEL_NC_011916.1_344130_m | nterm_extension | 344,130-345,041 (−) | 303 | 2 | 1.1e-4 | 0.026 | CCNA_00333 (0.91) | known start difference |
| NOVEL_NC_011916.1_435633_m | intergenic | 435,633-436,094 (−) | 153 | 4 | 4.4e-7 | 0.077 | CCNA_00425 (1.0) | known: RefSeq CDS not called by Bakta |
| NOVEL_NC_011916.1_494406_m | nterm_extension | 494,406-495,143 (−) | 245 | 3 | 4.4e-7 | 0.026 | CCNA_00478 (1.0) | known start difference |
| NOVEL_NC_011916.1_651319_p | nterm_extension | 651,319-652,815 (+) | 498 | 11 | 4.4e-7 | 0.026 | CCNA_00611 (0.98) | known start difference |
| NOVEL_NC_011916.1_712528_p | nterm_extension | 712,528-713,880 (+) | 450 | 3 | 8.5e-6 | 0.026 | CCNA_00657 (1.0) | known start difference |
| NOVEL_NC_011916.1_773437_m | intergenic | 773,437-773,736 (−) | 99 | 2 | 8.1e-4 | 0.077 | cdzI CCNA_03931 (0.98) | known (RefSeq) |
| NOVEL_NC_011916.1_773733_m | intergenic | 773,733-773,987 (−) | 84 | 3 | 4.4e-7 | 0.077 | cdzD CCNA_03932 (1.0) | known (RefSeq) |
| NOVEL_NC_011916.1_774018_m | intergenic | 774,018-774,311 (−) | 97 | 3 | 4.4e-7 | 0.077 | cdzC CCNA_03933 (0.89) | known (RefSeq) |
| NOVEL_NC_011916.1_829105_p | intergenic | 829,105-829,290 (+) | 61 | 4 | 3.6e-4 | 0.077 | CCNA_03937 (1.0) | known (RefSeq small CDS) |
| NOVEL_NC_011916.1_1111617_m | intergenic | 1,111,617-1,112,111 (−) | 164 | 8 | 4.4e-7 | 0.077 | CCNA_01024 (1.0) | known (RefSeq) |
| NOVEL_NC_011916.1_1246920_p | intergenic | 1,246,920-1,247,396 (+) | 158 | 2 | 1.9e-3 | 0.077 | CCNA_01142 (1.0) | known (RefSeq) |
| NOVEL_NC_011916.1_1422154_m | nterm_extension | 1,422,154-1,422,990 (−) | 278 | 3 | 4.4e-7 | 0.026 | CCNA_01297 (1.0) | known start difference |
| NOVEL_NC_011916.1_1634200_p | nterm_extension | 1,634,200-1,634,691 (+) | 163 | 5 | 4.4e-7 | 0.026 | flbT CCNA_01525 (0.87) | known start difference |
| NOVEL_NC_011916.1_1638743_p | nterm_extension | 1,638,743-1,639,297 (+) | 184 | 3 | 4.4e-7 | 0.026 | CCNA_01531 (1.0) | known start difference |
| NOVEL_NC_011916.1_1751025_m | intergenic | 1,751,025-1,751,486 (−) | 153 | 5 | 7.9e-7 | 0.077 | CCNA_01632 (1.0) | known (RefSeq) |
| NOVEL_NC_011916.1_1849629_p | nterm_extension | ..-1,850,426 (+), no start codon | 265 | 2 | 4.1e-5 | 0.026 | CCNA_01723 (1.0) | upstream of a RefSeq gene, no start codon: exploratory |
| NOVEL_NC_011916.1_1912227_m | intergenic | 1,912,227-1,912,619 (−) | 130 | 5 | 1.1e-4 | 0.17 | CCNA_01784 (1.0) | known (RefSeq) |
| NOVEL_NC_011916.1_1937614_m | antisense | 1,937,614-1,938,561 (−) | 315 | 2 | 3.5e-3 | 1.0 | none (antisense to CCNA_01813) | **new candidate**; class FDR 1.0, likely false |
| NOVEL_NC_011916.1_2130432_m | nterm_extension | 2,130,432-2,132,201 (−) | 589 | 5 | 2.1e-6 | 0.026 | CCNA_01985 (1.0) | known start difference |
| NOVEL_NC_011916.1_2148110_m | nterm_extension | 2,148,110-.. (−), no start codon | 215 | 2 | 5.2e-7 | 0.026 | CCNA_02002 (1.0) | upstream of a RefSeq gene, no start codon: exploratory |
| NOVEL_NC_011916.1_2933516_p | intergenic | 2,933,516-2,934,097 (+) | 193 | 5 | 4.4e-7 | 0.17 | CCNA_02774 (0.79) | known (RefSeq) |
| NOVEL_NC_011916.1_2939350_p | intergenic | 2,939,350-2,939,946 (+) | 198 | 2 | 4.4e-7 | 0.077 | CCNA_02781 (1.0) | known (RefSeq) |
| NOVEL_NC_011916.1_3101435_p | nterm_extension | 3,101,435-3,102,016 (+) | 193 | 2 | 4.4e-7 | 0.026 | greA CCNA_02941 (0.81) | known start difference |
| NOVEL_NC_011916.1_3118346_p | intergenic | 3,118,346-3,119,035 (+) | 229 | 3 | 2.4e-5 | 0.17 | CCNA_03983 (0.38) | overlap < 50 %: **new candidate**, next to CCNA_02961 |
| NOVEL_NC_011916.1_3658588_m | nterm_extension | 3,658,588-3,659,565 (−) | 325 | 2 | 4.4e-7 | 0.026 | CCNA_03501 (1.0) | known start difference |

The paper does no proteogenomics, so no candidate is "known in the paper". RefSeq does not explain
two candidates (the antisense one at 1.94 Mb and the intergenic one at 3.12 Mb). Both are listed for
Khoa's review, not as discoveries.

## Route A vs route B

| | Route A | Route B | B / A |
|---|---|---|---|
| makedb records | 4,104 (CDS + sORF + 381 `CONT_`) | 89,471 (+ six-frame ORFs ≥ 20 aa) | 21.8x |
| FASTA bytes | 1,568,123 | 12,455,297 | 7.9x |
| predicted library | 688,052 precursors, 167,920,698 bytes | 4,118,665 precursors, 1,003,062,992 bytes | 6.0x |
| DIA-NN step 1 / step 2 (16 CPUs, 3 runs) | 1 min 51 s / 11 min 51 s | 9 min 23 s / 38 min 46 s | 3.5x total |
| precursors at 1 % per run | 39,493 / 38,287 / 37,937 | 38,801 / 37,418 / 37,409 | |
| protein groups, global q ≤ 0.01 | 2,924 | 2,917 | |
| GIPPA detected genes | 3,005 | 2,938 | keeps 97.1 % [96.4, 97.6] |
| `gippa run` wall time (14 threads, 54 GB, `farmshare_local` in a 16-CPU job) | 3 h 33 min | 3 h 21 min | |

## Command log (pinned commit `bfc74ff6`, `PYTHONPATH=<pin>/src`, venv Python)

```bash
gippa makedb --annotation $D/bakta/ccre_na1000/out/ccre_na1000.json --outdir $D/db/ccre_na1000/db_a      # job 2005458: exit 0, 2 s
gippa makedb --annotation $D/bakta/ccre_na1000/out/ccre_na1000.json --six-frame-orfs --min-orf-aa 20 \
  --outdir $D/db/ccre_na1000/db_b                                                                       # job 2005458: exit 0, 7 s
# DIA-NN 2.3.1, two calls per route (p2a/diann_search.sh), 3 runs                                       # job 2006169
gippa convert --from diann $D/ccre_na1000/diann_a/report.parquet --out $D/ccre_na1000/in_a/peptides.tsv.gz \
  --abundance-out $D/ccre_na1000/in_a/abundance.tsv.gz --search-db same          # job 2006169, twice, byte-identical
gippa preflight --genome $D/genome/ccre_na1000/ccre_na1000.fna.gz --annotation .../ccre_na1000.json \
  --peptides $D/ccre_na1000/in_a/peptides.tsv.gz --outdir $D/ccre_na1000/in_a/preflight --strict   # exit 0
gippa convert ... diann_b/report.parquet ... --search-db same ; gippa preflight ... in_b ... --strict  # job 2006169: exit 0
gippa run --profile farmshare_local,apptainer --genome ... --annotation ... --peptides $D/ccre_na1000/in_a/peptides.tsv.gz \
  --db-dir $GIPPA_DB_DIR --threads 14 --memory 54.GB --outdir $D/ccre_na1000/out_a   # job 2006169: exit 0, 3 h 33 min
gippa run ... --peptides $D/ccre_na1000/in_b/peptides.tsv.gz ... --outdir $D/ccre_na1000/out_b   # job 2006169: exit 0, 3 h 21 min
gippa compare $D/ccre_na1000/out_a $D/ccre_na1000/out_b --labels genes,sixframe --outdir $D/ccre_na1000/cmp_ab   # allocation: exit 0, 8 s
gippa report $D/ccre_na1000/out_b --redact                                       # allocation: exit 0, 77 s; tables/ unchanged (sha256)
```

All modules on (`modules_off` empty), default seed, no `--abundance`. Each run had 110-116 Nextflow
tasks, 0 failed.

## Problems and deviations

- `gippa run` used the `farmshare_local` profile inside a 16-CPU lane job (`p2a/lane.sbatch`, as
  Phase 1b `gippa_arm.sbatch`) rather than `scripts/slurm/gippa_run.sbatch` (profile `farmshare`).
  This kept the agent within its 32-CPU share; the command is otherwise the same.
- The lane log records `exit=2` for route B because `arm.sh` was edited while the lane was running,
  and bash read the changed file at its last line. Every step had already finished with exit 0
  (`run_manifest.json` exit 0, `runtime_out_b.txt` run_exit=0). Scripts are now never edited while
  a lane runs.
- Route A warns `novel_orfs.decoys_not_makedb` and `novel_orfs.layer_missing`; route A has no ORF
  decoys, so both are expected.
- TE7 determinism for group A is pre-registered on PXD065022 DAL-1 (see that report).
