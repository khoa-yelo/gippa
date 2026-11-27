# GIPPA preprint: supplementary tables

<!-- Tables moved from the main text in round 2 (review VIZ: redundancy). The src comment after
each table names its evidence; abbreviations as in paper.md (CAL111, RD, P2S, NPC). -->

**Table S1.** Commands and options of `gippa` 1.1.1 (from `gippa <command> --help`). Every command
also takes `--config FILE`, `--set KEY=VALUE` (any configuration key, repeatable), `--log-level`
and `--version`, except `demo` and `example`.

| Command | Purpose | Options |
|---|---|---|
| `setup` | download databases and models; writes `db_manifest.json` | `--db-dir`; `--bakta light\|full\|skip`; `--genomad/--no-genomad` (non-commercial use); `--antismash/--no-antismash`; `--promotech-models/--no-promotech-models`; `--go/--no-go`; `--kegg/--no-kegg` (academic use); `--cog/--no-cog`; `--contaminants/--no-contaminants`; `--minimal` (no geNOMAD, no antiSMASH); `--check`; `--images [all]`, `--no-images`; `--profile auto\|apptainer\|docker`; `--cache-dir`; `--install-uniop-from-git`, `--uniop`, `--no-uniop`; `--nextflow/--no-nextflow`; `--accept-download-plan` (accepts no licence terms) |
| `preflight` | validate inputs and installed tools | `--genome` (required); `--proteins`; `--peptides`; `--abundance`; `--annotation`; `--outdir`; `--strict` (warnings exit 3); `--profile` (`none`: inputs only) |
| `run` | run the analysis | `--genome` (required); `--outdir` (required); `--peptides`; `--proteins`; `--abundance`; `--annotation`; `--reference-peptides`; `--replicon-depth`; `--db-dir`; `--threads`; `--memory` (default 64.GB); `--modules core,mobile,bgc,function`; `--genome-relation related\|same`; `--abundance-missing median\|exclude\|min\|zero\|mean`; `--abundance-type ibaq\|top3\|maxlfq\|intensity\|unknown`; `--abundance-scale auto\|log2\|linear`; `--promoters targeted\|genome`; `--terminators terminatornet\|transtermhp`; `--add-bactermfinder`; `--publication` (10,000 permutations and bootstrap draws); `--seed`; `--parquet`; `--sensitivity basic\|full\|off`; `--strict`; `--profile` (standard, farmshare, farmshare_local, slurm, apptainer, conda, docker; default auto); `--resume`; `--executor nextflow\|local` |
| `makedb` | build a strain search database (no FDR is built on it) | `--annotation`; `--genome`; `--outdir` (required); `--db-dir`; `--six-frame-orfs`; `--min-orf-aa` (default 20); `--patch`; `--no-contaminants` |
| `convert` | search-engine output to the peptide table | `INPUTS…`; `--from diann\|casanovo\|instanovo\|cascadia\|generic` (required; casanovo and instanovo exit "not implemented"); `--out` (required); `--abundance-out`; `--source database\|denovo`; `--search-db same\|other` |
| `compare` | compare runs on the same genome | `RUNS…`; `--outdir` (required); `--labels`; `--allow-cross-major` |
| `report` | rebuild plots and report from tables | `OUT`; `--no-interactive`; `--plots`; `--strict` (size budget); `--redact` (shareable report) |
| `example fetch` | download, verify and unpack the example bundle | `--outdir`; `--version X.Y` (bundle version); `--from-file` (offline) |
| `demo` | bundled 100-kb simulated demo, no Nextflow, containers or network | `--outdir`; `--full` (through `gippa run` with the real tools); `--open` |

<!-- src: $GIPPA_ROOT/.. scratch capture of `gippa <cmd> --help` for run, setup, preflight, makedb,
convert, compare, report, example fetch, demo (gippa 1.1.1, this branch); methods.md §2
(Converters: Casanovo and InstaNovo exit 5 "not implemented") -->

**Table S2.** Observed error against nominal level, representative streams of the 1.1.1 full tier
(α or q = 0.05 unless stated). Bounds are nominal + 2 Monte Carlo SE (FDR rows: empirical SD of the
FDP). R replicates; m pooled tests per replicate.

| Method | Stream | Criterion | Observed | Bound | R (m) |
|---|---|---|---|---|---|
| placement | read-through, shuffled decoys | FDR | 0.0355 | 0.0510 | 300 |
| placement | contaminant attack, q = 0.01 (known failure; placement gate not met) | FDR | 0.0150 | 0.0102 | 300 |
| gene FDR | picked, π = 0.05, q = 0.01 | FDR | 0.00941 | 0.01031 | 300 |
| gene FDR | MAYU (PEP), π = 0.05, q = 0.01 | FDR | 0.00962 | 0.01022 | 300 |
| region tests | HL shift, chromosome pooled | level | 0.0505 | 0.0569 | 1,000 (4) |
| region tests | prophage HL shift (within tolerance, not at level: 26 % above q) | FDR | 0.063 | 0.0654 | 1,000 |
| autocorrelation | Moran's I lag 1, exact TUs | level | 0.055 | 0.0638 | 1,000 |
| autocorrelation | same, UniOP-like TU errors (recorded only; rows flagged) | level | 0.071 | 0.0638 | 1,000 |
| autocorrelation | same, UniOP's K-12 TU error rates (recorded only; rows flagged) | level | 0.432 / 0.659 / 0.758 at α 0.01 / 0.05 / 0.10 | 0.0163 / 0.0638 / 0.119 | 1,000 |
| autocorrelation | same, distance-model operon fallback (known failure; rows flagged) | level | 0.149 / 0.354 / 0.476 at α 0.01 / 0.05 / 0.10 | 0.0163 / 0.0638 / 0.119 | 1,000 |
| islands | AR(1) ρ = 0 / 0.5 / 0.9 | FWER | 0.022 / 0.020 / 0.024 | 0.0695 | 500 |
| islands | ρ = 0.9, short block length, 1,000 / 1,800 genes (outside the claim; recorded) | FWER | 0.096 / 0.130 | 0.0695 | 500 |
| operon discordance | gate null, 1,800 genes | level / FDR | 0.0499 / 0.034 | 0.0510 / 0.0615 | 1,000 (200) |
| operon discordance | step-shaped noise attack (known failure; claim narrowed to smooth noise) | level | 0.0517 / 0.1031 at α 0.05 / 0.10 | 0.0510 / 0.1013 | 1,000 (200) |
| missed operons | shipped decoys | FDR | no calls (vacuous) | – | 300 |
| GSEA | pooled null, R engine | level | 0.0415 | 0.0525 | 300 (100) |
| GSEA | unflagged rows under detection confounding, arm `null_b2` (Python port engine; gated by a 1.5α screen set in the test file, not in the tolerance file or an ADR) | FWER | 0.0833 at α 0.05; 0.1367 at α 0.10 | α + 2 SE: 0.0752 / 0.1346 (both exceeded); 1.5α screen: 0.1054 / 0.1912 (passed) | 300 |
| novel ORFs | class FDR, database regime | FDR | 0.0453 | 0.0511 | 300 |
| novel ORFs | de novo regime: multi-peptide false ORFs (known failure) | FDR | 0.0183 / 0.0788 / 0.1139 at q 0.01 / 0.05 / 0.10 | 0.0109 / 0.0511 / 0.1012 | 300 |
| novel ORFs | de novo regime: exact systematic errors, G4-k3, out of frame (known failure) | FDR | 0.0749 at q 0.01 and 0.05 | 0.0102 / 0.0502 | 300 |
| novel ORFs | de novo regime: G4-k3, antisense (known failure) | FDR | 0.0719 at q 0.01 and 0.05 | 0.0102 / 0.0502 | 300 |
| novel ORFs | de novo regime: G4-k3, intergenic (known failure) | FDR | 0.0153 at q 0.01 | 0.0114 | 300 |
| novel ORFs | de novo regime: exact GC-poor intergenic reads (known failure) | FDR | 0.0126 at q 0.01 | 0.0112 | 300 |
| variants | F1 of `high` ALT calls, scope before ADR-0049 (≥ 0.90; known failure; variant gate not met) | F1 | 0.874 (SE 0.0018) | 0.90 | 100 |
| variants | random / systematic isobaric-swap attacks (known failure) | FWER | 0.97 / 1.0 | 0.0936 | 100 |
| compare | ρ CI non-coverage, TU block | level | 0.048 | 0.0638 | 1,000 |

<!-- src: results.md Table R2 and Table R3b (CAL110 line references), re-read in CAL111 at the
same lines (identical values): lanes/placement-group1/placement.jsonl:2; lanes/gene_fdr-group1/
gene_fdr.jsonl:44,47; lanes/regions-group2/regions.jsonl:21; lanes/regions-group0/regions.jsonl:8;
lanes/autocorrelation-group0/autocorrelation.jsonl:77,44; lanes/islands-group0/islands.jsonl:32,35,38;
lanes/islands-group0/islands_recorded.jsonl:2; lanes/discordance-group1/discordance.jsonl:17,21;
lanes/missed_operons-group0/missed_operons.jsonl:35; lanes/gsea-group0/gsea.jsonl:45;
novel_orfs.merged.jsonl:38, 23 and the strict-xfail records of novel_orfs.merged.jsonl (streams
novel_orfs.attack.denovo.multi_peptide_orfs.full, novel_orfs.gc, novel_orfs.g4k3.intergenic,
novel_orfs.g4k3.cds_outframe, novel_orfs.g4k3.antisense; block 2; tests test_attack_novel_class_fdr.py,
test_attack_m7_gate.py, test_novel_orfs_gate.py; ADR-0050); variants.merged.jsonl:3;
lanes/compare-group0/compare.jsonl:41; GSEA null_b2: CAL111 lanes/gsea-group1/gsea.jsonl:20-21
(screen levels 0.075 and 0.15; observed 0.08333, 0.13667; R 300; α + 2 SE bounds computed:
0.07517, 0.13464); screen constant tests/calibration/attacks/test_attack_gsea_flags.py:367
(SCREEN = 1.5); these records are at the screen levels 0.075 and 0.15, so they are not among the
1,244 records of Fig. 2b; prophage 26 %: results.md §1.2 line 199; contaminant attack:
placement.merged.jsonl:43; UniOP K-12 rates: lanes/autocorrelation-group0/autocorrelation.jsonl:47,49,51;
distance model: :23,25,27 (strict xfail); islands_recorded.jsonl:1-2; discordance step:
lanes/discordance-group*/discordance.jsonl:8,10 (strict xfail); variants.merged.jsonl:1-2 -->

**Table S3.** Campaign datasets and headline results (route A unless stated; shares with Wilson
95 % CI, optimistic because genes and peptides co-vary).

| Organism (genome) | Data | Runs, instrument | Genes detected / coding (share) | Estimated gene FDR at q ≤ 0.01 (route A) | Route B keeps | Novel-ORF candidates (called) |
|---|---|---|---|---|---|---|
| *E. coli* K-12 MG1655 | PXD058422 | 4, Orbitrap Eclipse | 2,691 / 4,368 (0.616 [0.602, 0.630]) | 0.0093 | 0.984 [0.979, 0.988] | 124 (0) |
| *T. pallidum*, 7 strains | PXD065022 | 3 per strain, timsTOF Pro | 668–880 / 972–976 (0.685–0.903) | 0.0011–0.0060 | 0.973–0.999 | 200 (0) |
| *Caulobacter crescentus* NA1000 | PXD069004 | 3, Exploris 480 | 3,005 / 3,723 (0.807 [0.794, 0.819]) | 0.0093 | 0.971 [0.964, 0.976] | 164 (0) |
| *Bacteroides thetaiotaomicron* | PXD055395 | 16, timsTOF Pro | 2,244 / 4,931 (0.455 [0.441, 0.469]) | 0.0098 | 0.972 [0.964, 0.978] | 99 (0) |
| *Synechococcus* sp. CC9311 | PXD072227 | 5, Q Exactive | 1,897 / 2,859 (0.664 [0.646, 0.681]) | 0.0084 | 0.968 [0.960, 0.975] | 117 (0) |
| *Vibrio natriegens* | PXD049476 | 4, Exploris 480 | 2,515 / 4,551 (0.553 [0.538, 0.567]) | 0.0099 | 0.969 [0.961, 0.975] | 78 (0) |
| *Myxococcus xanthus* DK1622 | PXD060688 | 4, Exploris 480 | 4,864 / 7,318 (0.665 [0.654, 0.675]) | 0.0097 | 0.985 [0.981, 0.988] | 265 (0) |
| *Clostridioides difficile* 630Δerm | PXD065346 | 3, DIA | 1,928 / 3,824 (0.504 [0.488, 0.520]) | 0.0099 | 0.988 [0.982, 0.992] | 53 (0) |
| *Mycobacterium tuberculosis* H37Rv | PXD057335 | 3, Exploris 480 | 2,942 / 4,075 (0.722 [0.708, 0.736]) | 0.0099 | 0.977 [0.971, 0.982] | 247 (0) |
| *Haloferax volcanii* H26 (DS2 genome) | PXD066008 | 5, timsTOF Pro 2 / HT | 2,978 / 4,023 (0.740 [0.726, 0.754]) | 0.0097 | 0.972 [0.966, 0.977] | 142 (0) |
| *Saccharolobus islandicus* REY15A | PXD079213 | 3, Orbitrap Astral | 2,109 / 2,807 (0.751 [0.735, 0.767]) | 0.0085 | 0.983 [0.976, 0.988] | 92 (0) |

<!-- src: P2S §2.2 (gene FDR), §2.3 (runs, instrument, detected/coding, Wilson CIs), §2.4 (route B
kept), §2.6 (candidates, 0 called); CAMP/ecoli/REPORT.md Key numbers (E1 2,691, FDR 0.0093 1.1.1
recheck 0.00929), row 14 (98.4 % [97.9, 98.8]), row 9/15 (124 candidates, 0 called); E. coli
share 2,691/4,368 Wilson computed by fig4_campaign.py (0.616 [0.602, 0.630]); T. pallidum per
strain: PXD065022/metrics.json › T2.per_strain and strains.*.metrics.S9_route_b (fig4_campaign.tsv);
accessions: handover plan Phase 2 table -->

**Table S4.** Release readiness of GIPPA 1.1.1 (acceptance rows).

| Row | Criterion | Status | Value |
|---|---|---|---|
| R1 | ≥ 99.5 % of A2 peptides map exactly | pass | 100.00 % |
| R2 | decoy-estimated mapping FDR ≤ 1 % | **waived** (ADR-0071, ADR-0080) | not measured (same-database search) |
| R3 | every statistical method passes its calibration suite | pass (no waiver; expected failures listed) | full tier 558 tests: 462 passed, 90 skipped, 6 xfailed, 0 failed; fast tier 29,049 passed, 149 skipped, 17 xfailed, 3 xpassed; coverage 97 % |
| R4 | every attack null stays at level | pass (expected failures listed) | 284 attack tests: 276 passed, 3 skipped, 5 xfailed. R4 passes with the GSEA detection-confounding arm `null_b2` gated at 1.5α, a screen not recorded in the tolerance file or an ADR; at the α + 2 SE bound that arm fails (Table S2) |
| R5 | every true-example assertion the data allow passes | **waived** for TE4 only (ADR-0070, ADR-0080) | TE1–TE3, TE5, TE7, TE8 pass; TE4 fails |
| R6 | same seed gives byte-identical tables | pass | TE7; integration 967 passed; nf-test 144 tests |
| R7 | standard run on 16 threads and 64 GB ≈ 35 min (16–55) | **waived** (ADR-0072, ADR-0080) | 3 h 41 min |
| R8 | core modules ≈ 25 min | pass (threshold "reported only") | 3 h 01 min |
| R9 | the pilot passes | pass | 5 criteria; TE4 reported, cause not established; the five criteria do not include the pre-registered hypotheses, and H1 was not supported |
| R10 | co-abundance against K-12 operons | **deferred** to full validation | – |

<!-- src: RD rows R1-R10 (status, value, threshold) -->
