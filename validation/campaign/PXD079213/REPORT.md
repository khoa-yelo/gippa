# Phase 2 report: PXD079213 (*Saccharolobus islandicus* E233S on REY15A)

Campaign Phase 2, group C, dataset 10. The data are three Orbitrap Astral DIA runs (`cont_1..3`, 8 min gradients) of the E233S control strain from iProX. They were searched against the complete REY15A genome GCF_000189555.1 (NC_017276.1, circular). The answer key is `ANSWER_KEY.md`. It was pre-registered in commit `9070f5be` before any spectrum was searched, and it has not been edited since. Code: the detached pin `local/worktrees/rel111-pin` at `bfc74ff6` (GIPPA 1.1.1). Every manifest records `git_dirty` false.

This report holds aggregates only: counts, shares, q-values, locus tags and 1-based coordinates. It contains no peptide sequence and no DIA-NN row. The numbers come from `metrics.json` (written by `p2c/p2_metrics.py`), `provenance.json` (written by `p2c/collect_provenance.py`) and `annotation_check.json` (written by `p2c/annotation_check.py`). Every share carries a Wilson 95 % CI in brackets. "unestimated" means neither GIPPA nor this report gives an error estimate for that number.

## Archaeal annotation check (run before any search)

| | value | rule | result |
|---|---|---|---|
| Bakta CDS + sORF / RefSeq non-pseudo CDS | 2,808 / 2,578 = 1.089 | 0.90-1.10 | pass |
| RefSeq CDS sharing their stop with a Bakta CDS | 2,536 / 2,578 = 98.4 % [97.8, 98.8] | ≥ 90 % | pass |

**Pass: the run used Bakta.** Bakta calls 2,664 of its 2,808 CDS hypothetical. Name-based checks therefore use RefSeq names carried over to the Bakta CDS that has the same stop (same contig, strand and 3' end). That covers 2,536 genes and 62 of the 64 RefSeq ribosomal proteins. The two that do not map are L37ae and L40e (see the exploratory section).

## Verdicts

| # | Check | Verdict | Number (error estimate) | Cause / note |
|---|---|---|---|---|
| 1 | Mapping, TE1 (≥ 99.5 % of route-A peptides exact in frame) | **pass** | 23,052 / 23,052 distinct accepted peptides exact in frame at their own gene: 100 % [99.983, 100] | 57 peptides have a non-ATG N-terminus, all mapped. 0 masked, 0 misses. Placement FDR does not apply (`--search-db same`): 0 placements accepted (n_diff 1: 1 target / 16 decoys; n_diff 2: 0 / 1) |
| 2 | Expression, TE3 (ribosomal ≥ 90 % detected; median abundance in the top quartile) | **pass** | Answer-key set (RefSeq "30S/50S ribosomal protein", carried by stop): 60 / 62 detected, 96.8 % [89.0, 99.1]. Median iBAQ log2 28.43 vs Q3 25.46 (percentile rank 0.960, 60 measured of 2,107) | Not detected: S30e (CGDCIA_01228 / SIRE_RS05825) and L39e (CGDCIA_01755 / SIRE_RS08375). The pilot's strict product rule on the same names keeps 32 proteins: 32 / 32, 100 % [89.3, 100], median 28.66 (rank 0.967). GIPPA's own Bakta-named set is empty (`not_evaluable`, `no_ribosomal_genes`). Its §7.1.2 set holds 1 gene: 1 / 1 [20.7, 100] |
| 3 | Stoichiometry, TE4 | **not testable** (by construction) | rplL/rplJ and atpA/atpD: `gene_missing` for all four | Archaeal stalk proteins and an A-type ATP synthase, with no matching Bakta names. Stated in the key before the search |
| 4 | Contaminants, TE5 | **pass** | 265 / 265 trypsin (4) and keratin (261) peptides in the contaminant bucket: 100 % [98.6, 100]. 0 placed on the genome | 0 placed by the tolerant route, 0 absent |
| 5 | Gene-level FDR (route A ≤ 0.01) | **pass** | 2,109 genes at q ≤ 0.01 (tier `decoy_picked`), decoy-estimated FDR 0.0085. Route B: 2,076 genes, 0.0053 | Route A: 275 gene-paired and 2 unpaired decoys, paired share 0.993, q floor 4.7e-4. 1,928 genes have ≥ 2 peptides |
| 6 | Operons (ribosomal-protein and ATP-synthase units, RefSeq names) | **pass** | 37 adjacent pairs among 70 set genes (30 OP, 3 NO, 2 DV, 2 SW). 0 pairs with within-operon q ≤ 0.05 (min q 0.82). 0 `operon_flags` rows, 0 discordant | Median abs(Δlog2) 1.13 |
| 7 | Genome context (geNomad provirus regions) | **not testable** | 0 provirus regions. Region layer: 6 BGC, 4 CRISPR, 2 AMR | One island, ISL_1: low, NC_017276.1 1,653,011-1,665,971, 21 genes, detected fraction 0.286 (B 0.238), effect −1.18. The island carries no q-value |
| 8 | Negative controls (pyrE, pyrF, lacS not detected in A and B) | **fail** (vector caveat) | lacS (CGDCIA_02367) not detected in A or B (0 peptides). pyrE (CGDCIA_01488) detected: 9 peptides in A, 11 in B, gene q 4.9e-4, coverage 0.37 / 0.47. pyrF (CGDCIA_01487) detected: 2 peptides in A and B, q 4.9e-4, flag `top3_lt3` | Fits the empty-vector caveat stated in the key before the search. Exploratory, post hoc: see below. Under the key, this is not called a GIPPA error, and the verdict stays fail against the key |
| 9 | Novel ORFs vs RefSeq | **not testable** | 0 ORFs called (92 candidates) | No call to compare. Exploratory figures for all candidates are below |
| 10 | Class q floor | **pass** | Route B warns `novel_orfs.class_q_floor_above_alpha`. The classes hold 18-30 candidates (< 100). Minimum class FDR: 0.05 nterm_extension, 0.143 intergenic and cds_outframe, 0.857 antisense | Alpha 0.01 was not loosened, so there are 0 calls |
| 11 | Route B cost (keeps ≥ 95 % of route-A genes) | **pass** | B keeps 2,073 / 2,109 route-A genes: 98.3 % [97.6, 98.8]. 3 genes are detected only in B | Database, library and wall times are below. Abundance agreement over 2,071 shared genes (1,400 units): Spearman 0.996 [0.994, 0.997]. Bland-Altman bias −0.13 log2, limits of agreement [−1.26, 1.01] |
| 12 | Variants | **not testable** (same-genome search) | 0 rows in `variants.tsv` (A and B) | |
| 13 | Second round | **not testable** | not run | Route B called no ORF |

Counts: **7 pass, 1 fail (row 8, under the vector caveat), 5 not testable**, out of 13 checks.

## Key numbers

| | Route A | Route B |
|---|---|---|
| DIA-NN precursors at 1 % (per run) | 28,553-28,926 | 27,555-27,683 |
| DIA-NN protein groups, global q ≤ 0.01 | 2,129 | 2,146 |
| GIPPA genes detected (of 2,807 coding genes) | 2,109: 75.1 % [73.5, 76.7] | 2,076: 74.0 % [72.3, 75.5] |
| Detected with ≥ 2 peptides | 1,928 | 1,909 |
| Gene-level FDR at q ≤ 0.01 | 0.0085 | 0.0053 |
| Novel-ORF candidates / called | - | 92 / 0 |

These are 8-min gradients, so the answer key expected fewer detected genes than in the other datasets.

## Negative controls: exploratory check (post hoc)

This check was not pre-registered and was done after the pyrE and pyrF rows had been read. It is exploratory only. On 2026-11-22, NCBI E-utilities served the *S. solfataricus* orthologues PyrE (WP_009991146.1, orotate phosphoribosyltransferase) and PyrF (WP_009991148.1, orotidine-5'-phosphate decarboxylase). They are kept only under `$GIPPA_ROOT/campaign/PXD079213/sso_markers/` (sha256 `effe58af7137…` in `SHA256SUMS`), not in git. All 11 route-A peptides (pyrE 9, pyrF 2) also occur in these orthologues. So the evidence fits the expected source: the *S. solfataricus* pyrEF marker on an empty expression vector in the control strain. It does not show that the chromosomal REY15A genes are expressed. lacS (also a usual vector marker) has 0 peptides. Route B's two additional pyrE peptides are not part of the count of 11. None of this changes the verdict.

## Exploratory: what route B recovers

There are 42 RefSeq non-pseudo CDS without a Bakta stop, 31 of them ≤ 70 aa. Route B has exact peptides in 6 of them (3 of the ≤ 70 aa ones), and 4 of the 42 overlap a route-B candidate. None is called (class q floor, row 10).

| RefSeq | product | coordinates (strand), aa | route B | Reading |
|---|---|---|---|---|
| SIRE_RS14015 | 50S ribosomal protein L37ae | 1,215,193-1,215,405 (+), 70 | 3 peptides. Candidate NOVEL_NC_017276.1_1215184_p: cds_outframe, 73 aa, best q 3.2e-6, class FDR 0.143, RBS and promoter supported | A ribosomal protein that Bakta misses. Bakta calls CGDCIA_01368 over it in another frame |
| SIRE_RS14085 | 50S ribosomal protein L40e | 1,702,242-1,702,412 (+), 56 | 2 peptides. Candidate NOVEL_NC_017276.1_1702242_p: intergenic, 56 aa, best q 4.2e-6, class FDR 0.167, RBS and promoter supported | A ribosomal protein that Bakta misses. The candidate has the same coordinates as RefSeq |
| SIRE_RS15525 | TRASH domain-containing protein | 1,453,963-1,454,172 (−), 69 | 3 peptides. NOVEL_NC_017276.1_1453963_m: intergenic, 72 aa, class FDR 0.143 | Known RefSeq small CDS |
| SIRE_RS09225 | Mut7-C RNAse domain-containing protein | 1,678,259-1,678,789 (+), 176 | 2 peptides. NOVEL_NC_017276.1_1678217_p: intergenic, 190 aa, class FDR 0.143 | Known RefSeq CDS. The candidate starts 42 nt upstream |
| SIRE_RS14850 | FtsX-like permease family protein | 617,959-618,246 (−), 95 | 1 peptide, no candidate | |
| SIRE_RS15205 | hypothetical protein | 994,255-996,168 (−), 637 | 1 peptide, no candidate | |

L37ae and L40e are the two RefSeq ribosomal proteins left out of the TE3 set (62 of 64 map by stop). Both carry route-B peptides.

All 92 candidates compared with RefSeq (≥ 50 % same-strand overlap with a RefSeq CDS or pseudogene):
- Overall: 25 / 92, 27.2 % [19.1, 37.0], including 1 pseudogene match.
- By class: nterm_extension 14 / 30, cds_outframe 8 / 18, intergenic 3 / 18, antisense 0 / 26.
- Candidates with ≥ 2 peptides: 5 / 16, 31.3 % [14.2, 55.6]. These are 10 nterm_extension (9 flagged `no_start_codon`, class FDR 0.05), 3 intergenic (the three above), 2 cds_outframe (L37ae and NOVEL_NC_017276.1_1148017_p inside CGDCIA_01276) and 1 antisense (29 aa, class FDR 0.857).

As in Phase 1, a candidate without a start codon counts as a non-match, which pulls these shares down.

## Route A vs route B

| | Route A | Route B | B / A |
|---|---|---|---|
| makedb records | 3,188 (2,807 genes + 381 `CONT_`) | 62,578 (+ 59,390 six-frame ORFs ≥ 20 aa) | 19.6x |
| FASTA size | 1,025,394 bytes | 6,194,840 bytes | 6.0x |
| makedb wall time | 3 s | 6 s | |
| predicted library | 457,425 precursors, 111,157,989 bytes | 1,418,541 precursors, 354,197,030 bytes | 3.1x |
| DIA-NN step 1 / step 2 (3 runs, 16 CPUs) | 81 s / 415 s | 194 s / 534 s | 2.4x / 1.3x |
| `gippa run` (14 threads, 54 GB, `farmshare_local`) | 4,185 s (1 h 09 min 45 s) | 3,648 s (1 h 00 min 48 s) | |
| whole route job | 4,778 s, wheat-06 | 4,408 s, wheat-05 | |
| genes detected | 2,109 | 2,076 | 98.3 % kept |

`gippa compare` stability between the routes:
- islands: Jaccard 1.0 (21 / 21 genes shared);
- operon flags, enriched categories, novel ORFs and ALT sites: empty in both routes.

## Command log

All commands ran from the pin at `bfc74ff6` with `PYTHONPATH=<pin>/src`. `$D` = `$GIPPA_ROOT/campaign/PXD079213`, and the scripts are in `validation/campaign/p2c/`.

```bash
# fetch (iProX, download.iprox.cn): job 2001986 task 13 refused by the host allowlist (fixed); retry 2003431_13
# Bakta 1.12.2, DB 6.0 full: job 2001987_16 (wheat-07, 1,033 s)
python annotation_check.py $D/bakta/sisl_rey15a/out/sisl_rey15a.json GCF_000189555.1_ASM18955v1_genomic.gff.gz
sbatch -J gippa-camp-p2c-sisl-a route.sbatch PXD079213 sisl_rey15a a      # 2005460: exit 0, 1 h 19 min 45 s
sbatch -J gippa-camp-p2c-sisl-b route.sbatch PXD079213 sisl_rey15a b      # 2006191: exit 0, 1 h 13 min 36 s
#  each route job: sha256sum -c raw/SHA256SUMS; gippa makedb [--six-frame-orfs --min-orf-aa 20];
#  DIA-NN 2.3.1 (diann_search.sbatch, two calls); gippa convert --search-db same (twice, byte-identical);
#  gippa preflight --strict (exit 0); gippa run --profile farmshare_local,apptainer --threads 14 --memory 54.GB
bash post.sbatch PXD079213   # in the session allocation 2005379 (iron-04):
#  gippa compare out_a out_b --labels genes,sixframe   exit 0, 8 s
#  gippa report out_a --redact   exit 0, 34 s, tables unchanged (sha256)
#  gippa report out_b --redact   exit 0, 56 s, tables unchanged (sha256)
python p2_metrics.py --dataset sisl ... ; python collect_provenance.py $D provenance.json
```

Both runs used all modules (`modules_off` empty), master seed 20261014, and exited 0.
- Route A warnings: `novel_orfs.decoys_not_makedb` and `novel_orfs.layer_missing` (`rbs_novel`).
- Route B warning: `novel_orfs.class_q_floor_above_alpha`, with no layer missing.

## Provenance

| item | value |
|---|---|
| GIPPA | `bfc74ff6c20acc2cfcf298cdf23cbb30d8fae0dc` (1.1.1), clean |
| DIA-NN | 2.3.1 Academia, image sha256 `1946452d56d2fb0f66e05ced473df88e53484e28040f1ecbd83c4a43405a6381`. cfg sha256 `dad24b38e315…` |
| Bakta | 1.12.2, DB 6.0 full, image `sha256:1e24ae2f435c46f876a322298ab96d994d15aadc2824e1f16a77fa9a1bc07f18` |
| genome | `sisl_rey15a.fna.gz` sha256 `5f5698b621c3e08a90c7c81c38e6e1dff6184ec93e8f1dc9c40defc5dd957369`. RefSeq FASTA and GFF MD5 checked against `md5checksums.txt` |
| Bakta JSON | `180c6a955535a52fa0017473d1892b14b4f685fcee777ffb24212e3a18eb1522` |
| raw | cont_1/2/3.raw: 1,580,759,346 / 1,560,198,727 / 1,564,899,673 bytes. iProX publishes no checksum, so the files were verified by byte size against the iProX file list. sha256 `d01cbd6c227e…` / `4feb3514c01e…` / `2a52fed961d9…` (full values in `provenance.json`), rechecked by `sha256sum -c` in job 2005460 |
| DB FASTA A / B | `7d248cb37202…` / `e2c6d572d70a…` |
| peptides A / B | `f5190034e8bb…` / `92dc595f14cf…` |
| containers | antiSMASH 8.0.4 `d3f390e687cf…`; geNomad 1.12.0 `d54652ff868a…`; OSTIR 1.1.3 and ViennaRNA 2.6.4 `d9107b9b87fd…`; Promotech `7366eaf67d9b…`; Promoter Calculator 1.2.4 `c5bd0662ceed…`; TerminatorNet `ec415757b6b7…`; TransTermHP 2.09 `fc6fe0404be9…`; UniOp: no digest recorded (full digests in `provenance.json`) |

## Deviations and problems

- **TE3 set.** Bakta names no ribosomal protein, so GIPPA's own set is empty (`not_evaluable`). As the key fixed in advance, the verdict uses RefSeq names carried over to Bakta CDS by shared stop codon (62 proteins). The pilot's stricter product rule keeps only 32, because it drops archaeal names such as "S19e". Both sets pass.
- **TE7 rerun not done.** No determinism rerun was made for this dataset.
- **No second novel-ORF round.** Route B called no ORF, so there was no `makedb --patch` and no round 2.
- **Session allocation.** `post.sbatch` (compare and both `report --redact`, 98 s in total) ran as a script inside the session allocation 2005379 (iron-04), not as its own job. The metrics and provenance scripts (`annotation_check.py`, `p2_metrics.py`, `collect_provenance.py`) are short Python steps; where they ran is not recorded in the JSON files.
- **Route A warning `novel_orfs.layer_missing` (`rbs_novel`).** It is not in route B, and it did not appear in Phase 1 route A. The cause was not established, and it does not affect any route-A check.
- **Genome.** The search used the complete REY15A genome (one circular chromosome), not the 141-contig assembly the paper cites. E233S derives from REY15A and lacks pyrEF and lacS.
- **Fetch.** iProX first failed on the fetch script's host allowlist (job 2001986 task 13, 11 s). After the fix, job 2003431_13 fetched the files. There is no deposit checksum (see Provenance).
- **Contaminants.** As in Phase 1, `makedb` appends 381 `CONT_` records, so DIA-NN got one FASTA. Otherwise the flags are the pilot's `diann_A2.cfg`.
