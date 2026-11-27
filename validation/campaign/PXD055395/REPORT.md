# PXD055395 report: Bacteroides thetaiotaomicron VPI-5482 (Phase 2, group A, dataset 2)

Checks pre-registered in `ANSWER_KEY.md` (commit `854c143c`, before any search). Aggregates only.
Machine-readable numbers: `metrics.json`; per-route provenance: `provenance/{a,b}.json`. Shares
carry a Wilson 95 % CI in brackets. `$D` = `$GIPPA_ROOT/campaign/PXD055395`.

## Inputs

| Item | Value |
|---|---|
| Accession, runs | PXD055395 (PMC13253324), the 16 QC-passed B. thetaiotaomicron monoculture runs in mGAM, timsTOF Pro dia-PASEF, Evosep 100 SPD (11.5-min gradient); SHA-1 verified against the deposit (`PREP.md`; one archive refetched in prep job 2003431) |
| Genome | GCF_000011065.1: NC_004663.1 (6,260,361 bp) + plasmid p5482 NC_004703.1 (33,038 bp), both circular, no N or IUPAC bases |
| Annotation | Bakta 1.12.2, DB 6.0 full (job 2001987): 4,939 CDS + 2 sORF (4,931 coding rows in `genes.tsv`) vs RefSeq 4,687 CDS + 114 pseudogenes (ratio 1.054) |
| GIPPA | `bfc74ff6c20acc2cfcf298cdf23cbb30d8fae0dc` (1.1.1 release code, clean). Not rerun at `cff0fd68` (#612): the genome has no non-ACGT base, so the outputs would be identical |
| DIA-NN | 2.3.1 Academia, image sha256 `1946452d…a6381`, flags of `validation/pilot/diann_A2.cfg`; Bruker `.d` folders unpacked from the deposit zips (CRC-tested, prep job 2005458) |
| Containers | digests in `provenance/{a,b}.json` › `gippa_run.tools` |

## Verdicts

| # | Check | Verdict | Number (error estimate) | Note |
|---|---|---|---|---|
| S1 | Mapping, TE1 | **pass** | 24,854 / 24,854 route-A accepted peptides exact in frame at their own gene: 100 % [99.985, 100] | No miss of any kind. Placement FDR not applicable (`same_db`); no tolerant placement accepted |
| S2 | Expression, TE3 | **pass** | 53 / 55 ribosomal proteins detected, 96.4 % [87.7, 99.0]; median iBAQ log2 17.93 vs Q3 14.33 (percentile rank 0.969) | |
| S3 | Stoichiometry, TE4 | **pass** | log2(RplL/RplJ) = 1.44 (top3 1.38), in [1, 3]; log2(AtpA/AtpD) = −0.54 (top3 −0.53), in [−1, 1]. Error unestimated | |
| S4 | Contaminants, TE5 | **pass** | 18 / 18 trypsin (5) and keratin (13) peptides in the contaminant bucket, 100 % [82.4, 100]; 0 on the genome | Unmapped: 133 peptides, all labelled contaminant (107 `no_hit`, 26 `too_short`) |
| S5 | Operons | **pass** | 39 adjacent ribosomal/atp pairs (35 OP, 4 NO), all both detected; 0 with within-operon q ≤ 0.05 (min q 0.94); no discordant `operon_flags` row among them | Run-wide: 14 `conflict_terminator_in_operon` flags (route B: 15) |
| S6 | Genome context | **reported** | 28 region tests (16 geNomad prophage, 2 plasmid, 8 BGC, 2 AMR); 0 at q ≤ 0.05; no island in either route | |
| S7 | Novel ORFs vs RefSeq | **not testable** | 0 ORFs called (99 candidates) | Exploratory: 22 / 99 match a RefSeq CDS or pseudogene (22.2 % [15.2, 31.4]). By class: nterm_extension 9/14, cds_outframe 11/35, intergenic 2/12, antisense 0/38. 2 matches are RefSeq pseudogenes. Candidates with ≥ 2 peptides: 5 / 6 match (83.3 % [43.6, 97.0]) |
| S8 | Class q floor | **pass** | warning `novel_orfs.class_q_floor_above_alpha`; classes of 12-38 candidates; minimum class FDR 0.11 (nterm_extension), 0.36 (the other three); none called | Alpha not loosened |
| S9 | Route B cost | **reported** | route B keeps 2,181 / 2,244 route-A genes, 97.2 % [96.4, 97.8]; 54 genes detected only in route B. Spearman 0.989 [0.984, 0.992] over 2,174 shared genes (1,273 TUs); Bland-Altman bias −0.062 log2 | Cost table below |
| S10 | Gene-level FDR | **reported** | route A: 2,244 detected (2,011 with ≥ 2 peptides), `decoy_picked`, decoy FDR at q ≤ 0.01 0.0098. Route B: 2,235 (1,949), FDR 0.0094 | |
| S11 | Second round | **not run** | 0 ORFs called | |
| B1 | Depth with short gradients ([0.25, 0.70] of Bakta CDS+sORF) | **pass** | 2,244 / 4,931 = 45.5 % [44.1, 46.9] | The paper quantified 37 % on average per species (≥ 2 peptides in ≥ 2 replicates, DIA-NN 1.8.2 with MBR) |
| B2 | Sus operon expressed in mGAM (≥ 4 of 5) | **pass** | 5 / 5 detected at gene q ≤ 0.01 (100 % [56.6, 100]): susC (44 unique peptides), susD (26), susE (10), susF (12), susG (20); abundance percentiles 0.55 (susG) to 0.90 | |
| B3 | Plasmid p5482 | **reported** | plasmid 16 / 41 genes detected, 39.0 % [25.7, 54.3]; chromosome 2,228 / 4,890, 45.6 % [44.2, 47.0] | Prior (plasmid lower) holds in direction; the CIs overlap |
| B4 | Replicate depth (every run ≥ 4,000 precursors) | **pass** | route A precursors at 1 % per run: 14,722-21,355 (16 / 16 runs ≥ 4,000) | Route B: 14,772-21,110 |
| B5 | Annotation vs RefSeq | **reported** | Bakta 4,941 vs RefSeq 4,687 CDS; 2 route-B candidates overlap RefSeq pseudogenes | |

Counts: **9 pass, 0 fail, 2 not testable or not run (S7, S11), 5 reported** (S6, S9, S10, B3, B5).

## Candidates with ≥ 2 peptides (route B; none called)

| ORF | class | coordinates (strand) | aa | peptides | best q | class FDR | RefSeq (share of ORF) | Reading |
|---|---|---|---|---|---|---|---|---|
| NOVEL_NC_004663.1_3895716_m | nterm_extension | 3,895,716-3,896,999 (−) | 427 | 2 | 7.5e-7 | 0.11 | BT_RS15575 (1.0) | known: RefSeq start upstream of Bakta's |
| NOVEL_NC_004663.1_4240031_m | nterm_extension | 4,240,031-4,241,974 (−) | 647 | 2 | 9.9e-7 | 0.11 | BT_RS16730 (1.0) | known start difference |
| NOVEL_NC_004663.1_4469216_m | nterm_extension | 4,469,216-4,470,634 (−) | 472 | 2 | 4.3e-6 | 0.11 | BT_RS17555 (1.0) | known start difference |
| NOVEL_NC_004663.1_4594807_p | nterm_extension | ..-4,595,844 (+), no start codon | 345 | 4 | 6.7e-7 | 0.11 | prfB BT_RS17935 (1.0) | known: prfB programmed +1 frameshift upstream part, as in E. coli (Phase 1 row 8) |
| NOVEL_NC_004663.1_5962864_p | nterm_extension | 5,962,864-5,963,418 (+) | 184 | 2 | 3.3e-3 | 0.43 | BT_RS22900 (0.96) | known start difference |
| NOVEL_NC_004663.1_5983989_m | cds_outframe | 5,983,989-5,984,540 (−) | 183 | 3 | 8.3e-7 | 0.36 | BT_RS23000 (1.0) | known: a RefSeq CDS that Bakta calls in another frame |

No new candidate with ≥ 2 peptides. The paper does no proteogenomics.

## Route A vs route B

| | Route A | Route B | B / A |
|---|---|---|---|
| makedb records | 5,312 | 178,059 | 33.5x |
| FASTA bytes | 2,315,812 | 18,723,354 | 8.1x |
| predicted library | 1,128,893 precursors, 276,327,099 bytes | 5,246,413 precursors, 1,307,699,714 bytes | 4.6x |
| DIA-NN step 1 / step 2 (16 CPUs, 16 runs) | 3 min 02 s / 62 min 51 s | 15 min 44 s / 3 h 32 min 51 s | 3.5x total |
| protein groups, global q ≤ 0.01 | 2,241 | 2,185 | |
| GIPPA detected genes | 2,244 | 2,235 | keeps 97.2 % [96.4, 97.8] |
| `gippa run` wall time | 2 h 40 min | 2 h 47 min | |

## Command log (pinned commit `bfc74ff6`)

```bash
gippa makedb --annotation $D/bakta/bthe_vpi5482/out/bthe_vpi5482.json --outdir $D/db/bthe_vpi5482/db_a   # job 2005458: exit 0, 2 s
gippa makedb --annotation ... --six-frame-orfs --min-orf-aa 20 --outdir $D/db/bthe_vpi5482/db_b          # job 2005458: exit 0, 10 s
# DIA-NN 2.3.1 (p2a/diann_search.sh), 16 runs per route                                                   # job 2006169
gippa convert --from diann $D/bthe_vpi5482/diann_a/report.parquet --out .../in_a/peptides.tsv.gz \
  --abundance-out .../in_a/abundance.tsv.gz --search-db same     # job 2006169: twice, byte-identical
gippa preflight ... --peptides .../in_a/peptides.tsv.gz --strict   # exit 0; route B likewise: exit 0
gippa run --profile farmshare_local,apptainer ... --threads 14 --memory 54.GB --outdir $D/bthe_vpi5482/out_a   # job 2006169: exit 0, 2 h 40 min
gippa run ... --outdir $D/bthe_vpi5482/out_b                                                                  # job 2006169: exit 0, 2 h 47 min
gippa compare $D/bthe_vpi5482/out_a $D/bthe_vpi5482/out_b --labels genes,sixframe --outdir $D/bthe_vpi5482/cmp_ab   # allocation: exit 0, 8 s
gippa report $D/bthe_vpi5482/out_b --redact      # allocation: exit 0, 75 s; tables/ unchanged (sha256)
```

## Problems and deviations

- Same `farmshare_local` lane setup as PXD069004.
- DIA-NN advises MBR for a large library and fixed timsTOF mass accuracies. Both were left as in the
  pilot cfg (no flag changes).
- Route A warns `novel_orfs.decoys_not_makedb` / `layer_missing` (expected).
