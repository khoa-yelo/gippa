# PXD065022 report: Treponema pallidum, 7 strains (Phase 2, group A, dataset 1)

Checks pre-registered in `ANSWER_KEY.md` (commit `854c143c`, before any search). Aggregates only:
counts, shares, q-values, locus tags and 1-based coordinates. No peptide sequence and no DIA-NN row
is in git. Machine-readable numbers: `metrics.json` (`p2a/tpa_score.sbatch` → `p2a/p2_metrics.py --host`,
`p2a/score.py`, `p2a/placement_fdr.py`, `p2a/provenance.py`; assembled by `p2a/tpa_metrics.py`).
Shares carry a Wilson 95 % CI in brackets; "unestimated" means no error estimate exists.
`$D` = `$GIPPA_ROOT/campaign/PXD065022`. Strains: DAL-1, Grady, Haiti B, Madras, Mexico A,
Philadelphia 1 (Phil 1), SS14.

## Inputs

| Item | Value |
|---|---|
| Accession, runs | PXD065022 (J Proteome Res 2025, PMC12687308), timsTOF Pro DIA, 3 replicates per strain (21 `.d` runs); sha256 in `files.tsv` / `PREP.md` (`PHIL_A` refetched after a CRC failure) |
| Genomes | each strain on its own RefSeq assembly (`PREP.md`), one circular chromosome of 1,139,058-1,140,038 bp; every assembly has N runs (64-324 N; DAL-1 also one `R`) |
| Annotation | Bakta 1.12.2, DB 6.0 full (job 2001987): 975-979 CDS, 0 sORF per strain; 972-976 coding rows in `genes.tsv` |
| Host | UniProt UP000001811 (rabbit, release 2026_04), 21,205 sequences, MD5 `97516cae5c749aeb11cf5887864ab31a` = `RELEASE.metalink`; headers `CONT_RABIT_<acc>` |
| GIPPA | main `38e1cdba100662f746729c3e2e2eb0eb9a99828c` (1.1.1 code + #611/#612 Promotech, #614/#616 Promoter Calculator, #615 TerminatorNet gap fixes), detached pin `$GIPPA_ROOT/worktrees/rel111c-pin`, `git_dirty` false. The run manifests report `gippa_version` 1.1.0 because the version bump to 1.1.1 lives in release PR #609, not yet on main; makedb (run at `bfc74ff6`, the release-prep pin) reports 1.1.1 |
| DIA-NN | 2.3.1 Academia, image sha256 `1946452d…a6381`, flags of `validation/pilot/diann_A2.cfg` (`p2a/diann_campaign.cfg`); FASTA = strain `makedb` DB + rabbit FASTA |
| Containers | digests per tool in `metrics.json` › `strains.<tag>.prov_{a,b,x,te7}.gippa_run.tools` |

## Verdicts

| # | Check | Verdict | Number (error estimate) | Note |
|---|---|---|---|---|
| S1 | Mapping, TE1 (≥ 99.5 % exact in frame) | **pass** (7/7) | 100 % in every strain: DAL-1 8,350/8,350 [99.954, 100], Grady 2,746/2,746 [99.860, 100], Haiti B 9,284/9,284, Madras 4,920/4,920, Mexico A 5,567/5,567, Phil 1 9,010/9,010, SS14 5,094/5,094 (lowest CI bound 99.86 %) | 0 misses of any kind; 15-37 non-ATG N-terminal peptides per strain, all mapped. No tolerant placement accepted (1-difference stratum 5-14 targets, 0-2 decoys) |
| S2 | Expression, TE3 | **pass** (7/7) | ribosomal proteins detected: 51/51 (DAL-1, Madras), 50/51 (Haiti B, Mexico A, Phil 1, SS14), 46/51 = 90.2 % [79.0, 95.7] (Grady); median iBAQ log2 above Q3 in all (percentile rank 0.78-0.86) | Grady at the 90 % threshold by point estimate |
| S3 | Stoichiometry, TE4 | **fail** (6 pass, 1 fail) | log2(RplL/RplJ): DAL-1 2.72, Grady 2.87, Haiti B 1.41, **Madras 3.43**, Mexico A 2.89, Phil 1 2.96, SS14 2.76 (top3: 0.98-2.63). AtpA/AtpD **not testable** in all strains. Error unestimated | T. pallidum has a V-type ATPase (Bakta: atpA1/atpA2, atpD1/atpD2, ntp*), so no `atpA`/`atpD` pair is named. Madras top3 ratio 2.63 would be inside the band |
| S4 | Contaminants, TE5 | **pass** (7/7) | trypsin + keratin peptides in the contaminant bucket: 330-675 per strain (5-8 trypsin, 324-667 keratin), 100 % in each (lowest CI bound 98.8 %); 0 on the genome | |
| S5 | Operons | **pass** (7/7) | 33 adjacent ribosomal/ATPase pairs (29 OP, 4 NO), 27-33 both detected; 0 with within-operon q ≤ 0.05 (min q 0.71-0.91); no discordant `operon_flags` row among them | Each run has 1 operon flag (`conflict_terminator_in_operon`), not among these pairs |
| S6 | Genome context | **reported** | 2 region tests per strain (antiSMASH BGC), 0 at q ≤ 0.05, in both routes; no island in any run | 1.14 Mb genome: no geNomad/plasmid/AMR regions |
| S7 | Novel ORFs vs RefSeq | **not testable** | 0 ORFs called in all 7 strains (200 candidates: DAL-1 41, Grady 18, Haiti B 34, Madras 17, Mexico A 30, Phil 1 35, SS14 25) | Exploratory: 103/200 candidates match a RefSeq CDS/pseudogene (51.5 %); with ≥ 2 peptides 30/36. Table below |
| S8 | Class q floor | **pass** (7/7) | warning `novel_orfs.class_q_floor_above_alpha` in every route-B run; minimum class FDR 0.071 (DAL-1 nterm_extension) to 0.33; none ≤ 0.01, none called | Alpha 0.01 not loosened |
| S9 | Route B cost | **reported** | route B keeps route-A genes: DAL-1 868/869 99.9 % [99.4, 100.0], Grady 650/668 97.3 % [95.8, 98.3], Haiti B 876/879 99.7 % [99.0, 99.9], Madras 783/789 99.2 % [98.4, 99.7], Mexico A 798/815 97.9 % [96.7, 98.7], Phil 1 876/880 99.5 % [98.8, 99.8], SS14 792/805 98.4 % [97.3, 99.1]. B-only genes 3-13. Spearman (TU bootstrap) 0.971 (Grady) to 0.994; Bland-Altman bias −0.091 to +0.061 log2, SS14 −0.246 | Cost table below |
| S10 | Gene-level FDR | **reported** | route A detected / decoy FDR at q ≤ 0.01: DAL-1 869 / 0.0012, Grady 668 / 0.0060, Haiti B 879 / 0.0011, Madras 789 / 0.0038, Mexico A 815 / 0.0025, Phil 1 880 / 0.0011, SS14 805 / 0.0037. Route B: 871 / 0.0011, 663 / 0.0030, 882 / 0.0023, 792 / 0.0025, 808 / 0.0025, 879 / 0.0011, 798 / 0.0025 | Tier `decoy_picked` in all 14 runs |
| S11 | Second round | **not run** | 0 ORFs called | |
| T1 | Host share ≥ 50 % in every strain | **pass** (7/7) | distinct accepted route-A target peptides with only `CONT_` proteins: DAL-1 16,160/24,510 = 65.9 % [65.3, 66.5]; Grady 9,371/12,117 = 77.3 % [76.6, 78.1]; Haiti B 23,018/32,302 = 71.3 % [70.8, 71.7]; Madras 11,612/16,532 = 70.2 % [69.5, 70.9]; Mexico A 11,586/17,153 = 67.5 % [66.8, 68.2]; Phil 1 20,762/29,772 = 69.7 % [69.2, 70.3]; SS14 9,917/15,011 = 66.1 % [65.3, 66.8] | Rabbit-only 53.7-64.3 %; contaminant-library-only 6.9-9.4 % |
| T2 | Detection per strain | **pass** | (a) route-A detected / coding genes: DAL-1 869/973 = 89.3 % [87.2, 91.1], Grady 668/975 = 68.5 % [65.5, 71.4], Haiti B 879/973 = 90.3 % [88.3, 92.0], Madras 789/972 = 81.2 % [78.6, 83.5], Mexico A 815/973 = 83.8 % [81.3, 85.9], Phil 1 880/976 = 90.2 % [88.1, 91.9], SS14 805/976 = 82.5 % [80.0, 84.7], all in [0.50, 0.95]; (b) Grady fewest | Paper: 62.4-88.9 % of annotated proteins, Grady lowest |
| T3 | 12 newly annotated genes (≥ 9) | **pass** | 11/12 recovered, 91.7 % [64.6, 98.5] | Table below; miss: TP0126d (only shared peptides) |
| T4 | Translated pseudogene distal parts | **fail** (0/2) | TP0865a (Haiti B FOKJEF_00855): not detected, 0 peptides in 944,868-945,575 (−) in either route. TP0924a (Madras BIIDIP_00910): `partial_excluded`; Bakta types it as a partial pseudogene, so it is not in the route-A database; 0 peptides in 1,006,278-1,006,985 (+) in route A or B | Madras has 4 route-A peptides in the proximal part (BIIDIP_00909, 1,005,069-1,005,674), none in the distal part |
| T5 | Nichols-like vs SS14-like (≥ 4 of 5) | **fail** (3/5) | all 5 testable (7/7 strains detected); 653 genes detected in all seven used for centring. Agree: TP0487 (+0.60), TP0594 (−2.57), TP0732 folD (−2.42). Disagree: TP0477 pgl (−0.18), TP0597 (−0.12) | Sign of mean(Nichols) − mean(SS14) of centred route-A log2 abundance; error unestimated. Both misses are small differences |
| T6 | Cross-strain variants on DAL-1 | **pass** (6/6, placement FDR) | pooled TE8 placement FDR: Grady 0.073 %, Haiti B 0.065 %, Madras 0.061 %, Mexico A 0.091 %, Phil 1 0.101 %, SS14 0.059 % (all ≤ 1 %). Recall on covered sites 0 in every strain; `variants.tsv` empty; precision of `high` ALT undefined (0 calls) | Table below. Prior met: no tolerant placement accepted, so no variant can be called |
| T7 | Determinism, TE7 | **pass** | DAL-1 route A rerun (`out_a_fix2_te7`): `diff -r tables/` 0 lines, 22 files byte-identical (sha256 lists in `$D/score_fix2/tpa_dal1/`) | Same inputs and seed, same node type, run straight after route A in job 2014083 |

Counts: **10 pass** (S1, S2, S4, S5, S8, T1, T2, T3, T6, T7), **3 fail** (S3, T4, T5), **2 not
testable or not run** (S7, S11), **3 reported** (S6, S9, S10).

## T3: the 12 newly annotated genes

| Gene | Strain | Bakta locus | Route A | Route B | Recovered |
|---|---|---|---|---|---|
| TP0005a | DAL-1 | CLKMHN_00006 | detected, q 0.0012, 1 unique peptide | detected | yes |
| TP0134a | DAL-1 | (no same-stop CDS) | — | 1 accepted exact peptide inside 156,669-157,301 (+), in frame (157,113-157,148; q 0.0023), on candidate NOVEL_NC_016844.1_156594_p (cds_outframe, class FDR 0.20, not called) | yes |
| TP0126d | DAL-1 | CLKMHN_00129 | not detected (0 unique, 3 shared) | not detected | **no** |
| TP0127a | DAL-1 | CLKMHN_00130 | detected, 2 unique | detected | yes |
| TP0206a (rpmD) | DAL-1 | CLKMHN_00204 | detected, 2 unique | detected | yes |
| TP0409a (yajC) | DAL-1 | CLKMHN_00405 | detected, 4 unique | detected | yes |
| TP0453a | DAL-1 | CLKMHN_00447 | detected, 3 unique | detected, 4 unique | yes |
| TP0913a | DAL-1 | CLKMHN_00900 | detected, 2 unique | detected | yes |
| TP0921a | DAL-1 | CLKMHN_00908 (glmU) | detected, 3 unique | detected | yes |
| TP0950a (rnpA) | DAL-1 | CLKMHN_00937 | detected, 2 unique | detected | yes |
| TP0954a | DAL-1 | CLKMHN_00942 | detected, 5 unique | detected | yes |
| TP0250b (rpsT) | Grady | GKICKP_00254 | detected, q 0.0016, 3 unique | detected, 2 unique | yes |

Gene q-values of all detected rows: 0.0011-0.0016.

## T6: cross-strain runs on DAL-1 (`--search-db other`, `--genome-relation related`)

| Strain | truth sites | covered / eligible | recall, covered (any call) | unplaced N1 peptides | placed exact | 1-diff stratum (targets / decoys / accepted) | pooled placement FDR | N1 detected orthologs kept | Spearman (orthologs) |
|---|---|---|---|---|---|---|---|---|---|
| Grady | 218 | 18 / 14 | 0/18 [0, 17.6] | 24/2,753 = 0.87 % [0.59, 1.29] | 2,722 | 26 / 0 / 0 | 0.073 % | 658/663 | 0.9997 [0.9992, 1.0] |
| Haiti B | 82 | 21 / 13 | 0/21 [0, 15.5] | 30/9,296 = 0.32 % [0.23, 0.46] | 9,254 | 29 / 2 / 0 | 0.065 % | 871/872 | 0.9997 [0.9989, 1.0] |
| Madras | 26 | 2 / 1 | 0/2 [0, 65.8] | 5/4,925 = 0.10 % [0.04, 0.24] | 4,915 | 7 / 0 / 0 | 0.061 % | 786/786 | 0.9969 [0.9900, 1.0] |
| Mexico A | 247 | 48 / 36 | 0/48 [0, 7.4] | 50/5,571 = 0.90 % [0.68, 1.18] | 5,517 | 36 / 1 / 0 | 0.091 % | 808/808 | 0.9980 [0.9944, 0.9998] |
| Phil 1 | 221 | 42 / 32 | 0/42 [0, 8.4] | 56/9,022 = 0.62 % [0.48, 0.81] | 8,954 | 51 / 1 / 0 | 0.101 % | 872/873 | 0.9992 [0.9979, 0.9998] |
| SS14 | 224 | 28 / 22 | 0/28 [0, 12.1] | 38/5,098 = 0.75 % [0.54, 1.02] | 5,056 | 33 / 0 / 0 | 0.059 % | 797/799 | 0.9971 [0.9915, 0.9998] |

Truth sets: blastp RBH ≥ 90 % identity, ≥ 80 % coverage (job 2005467, built before any
cross-strain run). Every run's `mapping_fdr.tsv` sets `pi_wrong` = 1 in the 1- and 2-difference
strata (0-2 decoys), so tolerant placements are rejected (`rejected_placement_fdr`: 2-38 peptides
per strain) and no variant row is written. No non-truth `high` ALT exists to classify. Same
behaviour as Phase 1b.

## Candidates with ≥ 2 peptides (route B; none called)

Most candidates recur in several strains at the same position. Coordinates are 1-based, from the
longest-ORF start to the stop, given for the first strain listed. "RefSeq" is the same-strand
RefSeq CDS that overlaps most of the ORF (share of the ORF).

| Locus (strains) | class | coordinates (strand) | aa | peptides | class FDR | RefSeq | Reading |
|---|---|---|---|---|---|---|---|
| ~70.2 kb (DAL-1, Madras, SS14) | antisense | 70,185-70,397 (+) | 70 | 2-3 | 0.125-1.0 | none (antisense to TPADAL_RS00300) | **new candidate**, not called |
| ~195.5 kb (Haiti B, Phil 1) | nterm_extension | 196,915-195,503 (−) | 470 | 2 | 0.10-0.20 | WP_014342314.1 (1.0) | known start difference |
| ~212.2 kb (all but Grady) | intergenic | 212,241-212,426 (+) | 61 | 2 | 0.125-0.25 | WP_010881649.1 (1.0) | known: RefSeq CDS not called by Bakta |
| ~228.1 kb (all but Grady) | cds_outframe | 228,143-228,550 (+) | 135 | 2-5 | 0.125-0.25 | WP_010881670.1 (0.88-0.92) | known: RefSeq CDS in another frame |
| ~385.0 kb (all 7) | intergenic | 385,062-385,298 (+) | 78 | 2 | 0.5-1.0 | rpmB (0.79-1.0) | known: rpmB missed by Bakta in every strain |
| ~430.6 kb (Haiti B, Phil 1) | antisense (to Bakta) | 431,120-430,596 (−) | 174 | 2 | 0.17-0.25 | WP_237249931.1 (0.83) | known (RefSeq same strand) |
| ~443.3 kb (Haiti B, Mexico A) | antisense | 443,915-443,334 (−) | 193 | 4 | 0.20-0.25 | ffh (0.03); antisense to CRX41_RS02135 | **new candidate**, best q 2.9e-6, not called |
| ~467.2 kb (DAL-1, Haiti B, Madras, SS14) | nterm_extension | 467,277-468,518 (+) | 413 | 3-4 | 0.07-0.20 | WP_010881888.1 (0.93) | known start difference |
| 1,070,658 (DAL-1) | nterm_extension | 1,070,658-1,072,466 (+) | 602 | 3 | 0.071 | aspS (1.0) | known start difference |
| Haiti B only | nterm_extension | 496,572-495,403 (−); 729,561-729,268 (−); ..-1,080,937 (+, no start) | 389; 97; 289 | 2 each | 0.10 | WP_010881915.1, WP_010882111.1, WP_235214132.1 (1.0) | known start differences |

Exploratory, all candidates against RefSeq (≥ 50 % same-strand overlap): DAL-1 20/41, Grady 10/18,
Haiti B 15/34, Madras 13/17, Mexico A 16/30, Phil 1 16/35, SS14 13/25. By class, nterm_extension
almost always matches (50/56) and antisense rarely does (7/73). Against the RefSeq genes whose old
locus tag has a letter suffix (13 in DAL-1, 7 in Grady; the paper's full list of 100 is not
materialised as a file), only two DAL-1 candidates overlap one: the TP0134a candidate above
(TPADAL_0134c, 0.40 of the ORF) and NOVEL_NC_016844.1_214232_p (TPADAL_0206a, 0.26; 1 peptide).

## Route A vs route B

Ranges over the 7 strains (per-strain numbers in `metrics.json` › `strains.<tag>.prov_{a,b}`).

| | Route A | Route B | B / A |
|---|---|---|---|
| makedb records | 1,353-1,357 | 31,751-31,759 | 23.4x |
| FASTA bytes (strain DB) | 572,893-573,597 | 3,679,870-3,740,830 | 6.5x |
| predicted library (with rabbit) | 5.295-5.296 M precursors, 1.300 GB | 6.198-6.200 M precursors, 1.525 GB | 1.17x |
| DIA-NN step 1 / step 2 (16 CPUs, 3 runs) | 14-18 min / 2 h 28 min-3 h 40 min | 15-22 min / 2 h 51 min-5 h 00 min | about 1.2x |
| protein groups, global q ≤ 0.01 (rabbit + treponemal) | DAL-1 4,322, Grady 3,341, Haiti B 4,702, Madras 3,915, Mexico A 3,726, Phil 1 4,826, SS14 3,250 | 4,257, 3,286, 4,656, 3,862, 3,583, 4,736, 3,256 | |
| precursors at 1 % per run | 3,080-33,704 | 3,380-33,624 | |
| GIPPA detected genes | 668-880 | 663-882 | keeps 97.3-99.9 % |
| `gippa run` wall time (14 threads, 54 GB) | 35-55 min | 35-52 min | |

Most of the search space is the rabbit FASTA, so route B adds only 17 % precursors. Its cost is
much smaller than for the other two datasets (4.6-6.0x).

## Gap-fix warnings (non-ACGT bases; #611, #614, #615)

Read from the Nextflow task logs (`.command.err`) of every run; identical in route A, route B
and TE7 of the same genome. The cross-strain runs use the DAL-1 genome and show the DAL-1 values.

| Genome | non-ACGT bases removed by TerminatorNet | TerminatorNet calls dropped | TU upstream regions skipped (Promoter Calculator) | Promotech windows skipped |
|---|---|---|---|---|
| DAL-1 | 65 | 13 | 1 | 40 |
| Grady | 247 | 0 | 1 | 127 |
| Haiti B | 255 | 0 | 1 | 139 |
| Madras | 262 | 0 | 1 | 140 |
| Mexico A | 324 | 0 | 1 | 162 |
| Phil 1 | 262 | 0 | 1 | 127 |
| SS14 | 262 | 0 | 1 | 127 |

Every run also carries the `preflight:genome.gaps` warning in `run_manifest.json`.

## Command log (main `38e1cdba`, pin `rel111c-pin`, `PYTHONPATH=<pin>/src`, venv Python)

```bash
# makedb (at bfc74ff6; code unchanged since) and DIA-NN: jobs 2005458, 2006168, 2009201, 2011438, 2013337
# per strain and route, p2a/arm_v3.sh with SUF=_fix2, inside 16-CPU lane jobs (p2a/lane_v2.sbatch):
gippa convert --from diann $D/<tag>/diann_<r>/report.parquet --out $D/<tag>/in_<r>_fix2/peptides.tsv.gz \
  --abundance-out .../abundance.tsv.gz --search-db same          # twice, byte-identical, exit 0
gippa preflight ... --strict                                     # exit 3 (genome.gaps only) -> non-strict exit 0
gippa run --profile farmshare_local,apptainer --genome ... --annotation ... --peptides .../in_<r>_fix2/peptides.tsv.gz \
  --db-dir $GIPPA_DB_DIR --threads 14 --memory 54.GB --outdir $D/<tag>/out_<r>_fix2   # exit 0 (all 14)
# cross-strain (T6): convert --search-db other on DAL-1, gippa run ... --genome-relation related -> out_x_tpa_dal1_fix2
# TE7: gippa run (same inputs, same seed) -> tpa_dal1/out_a_fix2_te7
#   lanes 2014083 (DAL-1 A, TE7, Madras A), 2014084 (DAL-1 B, Grady B, Mexico A A), 2014085 (Grady A, Haiti B A),
#   2015973-2015976 (Phil 1 A/B, Haiti B B, Madras B, Mexico A B, SS14 A, 6 cross-strain arms), 2020015 (SS14 B): every arm exit 0
# scoring (p2a/tpa_score.sbatch, 8 CPUs): jobs 2020034, 2020035, 2020161, exit 0
gippa compare $D/<tag>/out_a_fix2 $D/<tag>/out_b_fix2 --labels genes,sixframe --outdir $D/<tag>/cmp_ab_fix2   # exit 0, 3-16 s
p2_metrics.py --a out_a_fix2 --b out_b_fix2 --cmp cmp_ab_fix2 --host --loci answer_loci.tsv ...
score.py $D/truth/<tag> out_a_fix2 out_x_tpa_dal1_fix2 ...      # T6, seed 20261014
placement_fdr.py in_x_tpa_dal1_fix2/peptides.tsv.gz <DAL-1 genome> <DAL-1 Bakta JSON> ...
diff -r tpa_dal1/out_a_fix2/tables tpa_dal1/out_a_fix2_te7/tables   # 0 lines
tpa_metrics.py metrics.json
```

All modules on (`modules_off` empty), default seed (master 20261014), no `--abundance`.
`layers_missing` is empty in route B; route A warns `novel_orfs.decoys_not_makedb` and
`novel_orfs.layer_missing` (expected: route A has no ORF decoys).

## Problems and deviations

- **Preflight strict waiver.** `--strict` turns any warning into exit 3. The only warning is
  `genome.gaps`: N runs in the RefSeq assemblies themselves, which no permitted input change can
  remove. `p2a/convert_preflight.sh` then reruns preflight without `--strict`, requires exit 0 and
  records `preflight_strict_exit` 3, `preflight_exit` 0 and `strict_waived` `genome.gaps` in
  `convert_provenance.json` (all 20 arms). Any other warning or error would still have stopped the arm.
- **Host as contaminant.** Rabbit proteins enter DIA-NN as `CONT_RABIT_*` records, so GIPPA's default
  contaminant prefixes bucket them; no GIPPA setting was changed (pre-registered).
- **GIPPA commit differs from the answer key.** The key fixed `bfc74ff6` (1.1.1 release code). That
  code could not finish a run on these genomes:
  - **`bfc74ff6` (`out_b`):** DAL-1 route B ended with exit 4 in `ANNOTATE:PROMOTECH`, on 40 windows
    with N (#611, fixed in #612).
  - **`cff0fd68` (`rel111b-pin`, `out_*_fix`):** DAL-1 A, DAL-1 B and Grady A ended with exit 4 in
    `ANNOTATE:PROMOTER_CALCULATOR` on a non-ACGT base (#614, fixed in #616).
  - **What counts:** all results above come from main `38e1cdba` (`_fix2` outdirs). Those failed
    outputs are history and are not scored. makedb, DIA-NN and convert code are unchanged across
    these commits, so the DIA-NN reports were reused.
- **S3 AtpA/AtpD.** Not testable by the pre-registered rule: T. pallidum has no F-type `atpA`/`atpD`
  pair (V-type ATPase). S3 is judged on RplL/RplJ alone and fails on Madras (3.43 > 3).
- **T4 Madras.** BIIDIP_00910 is typed as a pseudogene by Bakta, so the route-A database cannot
  contain it. Route B's six-frame ORFs cover the region but show no peptide there.
- **T3 not-gated list.** The paper's list of 100 letter-suffix genes was not materialised as a file
  before the search. Only the 12 genes of T3 and the RefSeq `old_locus_tag` letter-suffix genes are
  compared.
- **Slurm.** SS14 route B (job 2020015) was queued with a 12 h limit and waited on priority. The
  limit was lowered to 3 h (`scontrol update`); the job then started and used 37 min.
- **Bugs:** #611 (Promotech) and #614 (Promoter Calculator), both opened by this agent, are fixed
  at the commit used.
- **DIA-NN cost.** The search is dominated by the rabbit library (5.3 M precursors): about 2.5-5 h
  per strain and route on 16 CPUs, against 12-63 min for the other two datasets.
