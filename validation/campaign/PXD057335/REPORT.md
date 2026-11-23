# Phase 2 report: PXD057335 (*M. tuberculosis* H37Rv, DMSO vehicle, 48 h)

Campaign Phase 2, group C, dataset 7. The answer key `ANSWER_KEY.md` is pre-registered: it was committed
in `9070f5be` before any spectrum was searched and has not been edited since. Code: the detached worktree
`local/worktrees/rel111-pin` at `bfc74ff6` (GIPPA 1.1.1). Both run manifests record `git_dirty` false.

This report holds aggregates only: counts, shares, q-values, locus tags and 1-based coordinates. Git
holds no peptide sequence and no DIA-NN row. Every number below comes from `metrics.json` (written by
`p2c/p2_metrics.py`) or `provenance.json` (written by `p2c/collect_provenance.py`). Shares carry a
Wilson 95 % CI in brackets. "unestimated" means that neither GIPPA nor this report gives an error
estimate for that number.

## Verdicts

| # | Check | Verdict | Number (error estimate) | Cause / note |
|---|---|---|---|---|
| 1 | Mapping, TE1 (≥ 99.5 % of route-A peptides exact in frame) | **pass** | 32,851 / 32,851 accepted route-A peptides exact in frame at their own gene: 100 % [99.988, 100] | 0 misses. 111 peptides have a non-ATG N-terminus, all mapped; excluding them: 32,740 / 32,740 [99.988, 100]. 14 peptides masked by the pilot's te_proxies rule. Placement FDR not applicable (`--search-db same`): `mapping_fdr.tsv` has 0 accepted placement rows |
| 2 | Expression, TE3 (ribosomal ≥ 90 % detected; median abundance in top quartile) | **pass** | 54 / 56 ribosomal proteins (GIPPA set, Bakta names) detected: 96.4 % [87.9, 99.0]. Median iBAQ log2 28.10 vs Q3 23.76 (percentile rank 0.981; 54 of 2,962 measured) | §7.1.2 set: 54 / 58, 93.1 % [83.6, 97.3]. RefSeq-product set: the same 54 / 56. Not detected: L28 (LNBHIF_00117, Rv0105c) and S14 (LNBHIF_02185, Rv2056c). Bakta did not call RefSeq L36 (Rv3461c) or L34 (Rv3924c), so neither is in the set (see "RefSeq CDS without a Bakta call") |
| 3 | Stoichiometry, TE4 | **fail** (1 of 2 pairs) | log2(RplL/RplJ) = 0.61 (iBAQ; top3 0.58), outside [1, 3]. log2(AtpA/AtpD) = −0.02 (top3 −0.03), inside [−1, 1]. Error unestimated | The answer key named this a known risk. It is the same failure as the Nissle pilot and Phase 1 E1, and the cause is not established. Candidate causes: L7 N-acetylation is not searched, and iBAQ rests on few peptides (rplL has 6 observable peptides and 12 unique peptides seen; rplJ has 9 and 18) |
| 4 | Contaminants, TE5 | **pass** | 168 / 168 trypsin (3) and keratin (165) peptides in the contaminant bucket: 100 % [97.8, 100]. 0 placed on the genome | Also the 2 keratin-associated peptides: both in the bucket |
| 5 | Gene-level FDR (≤ 0.01) | **pass** | Decoy-estimated FDR 0.0099 over 2,942 route-A genes at q ≤ 0.01 (tier `decoy_picked`; q floor 3.4e-4) | Decoys: 394 gene-paired, 4 unpaired (paired share 0.990), win weight 1. 2,738 of the genes are detected with ≥ 2 peptides |
| 6 | Operons: ribosomal and ATP synthase partners concordant | **pass** | 38 adjacent pairs among 66 set genes (28 OP, 10 NO). 0 pairs with within-operon q ≤ 0.05 (min q 0.81). 0 `operon_flags` rows, 0 discordant | Median abs Δlog2 0.92. The set was built from Bakta gene names (`rp[lms]*`, `atp*`) plus product rules, as `p2_metrics.py` was fixed before reading |
| 7 | ESX-1: EsxA and EsxB at gene q ≤ 0.01 | **pass** | esxA (LNBHIF_04101, 4,352,609-4,352,896 +): q 3.7e-4, 2 unique peptides, coverage 0.375, flag `top3_lt3`. esxB (LNBHIF_04100, 4,352,274-4,352,576 +): q 3.7e-4, 5 unique peptides, coverage 0.66 | Both Bakta loci share their stop codon with RefSeq Rv3875 / Rv3874 |
| 8 | PE/PPE detected less than other CDS (Fisher p ≤ 0.05) | **pass** | PE/PPE: 19 / 150, 12.7 % [8.3, 18.9]. Other: 2,923 / 3,925, 74.5 % [73.1, 75.8]. Fisher p 1.1e-54, odds ratio 0.050 | **Measured differently from the pre-registration.** The key counted 155 RefSeq CDS. The script tests the 4,075 Bakta coding genes, matching the PE/PPE pattern on the RefSeq product carried by a shared stop. Only 3,789 genes carry one, so 150 genes are PE/PPE. The 286 genes without a RefSeq stop fall into "other". The effect is far beyond the threshold, so the verdict does not depend on this |
| 9 | Genome enrichment: prophage regions | **not testable** | geNomad called 0 prophage regions (phiRv1 and phiRv2 not called). 0 islands | Region types called: 34 BGC, 4 CRISPR, 2 AMR (exploratory, below) |
| 10 | Novel ORFs vs RefSeq | **not testable** | 0 route-B ORFs called (247 candidates), so there is no share to report | Exploratory, all candidates: see below |
| 11 | Class q floor | **pass** | Route B warns `novel_orfs.class_q_floor_above_alpha`. Candidates per class: antisense 76, cds_outframe 36, intergenic 26, nterm_extension 109. Min class FDR: 0.091 (antisense, cds_outframe, intergenic), 0.0118 (nterm_extension) | Alpha 0.01 was not loosened. nterm_extension has 109 candidates (slightly above the key's ~100) and its floor of 0.0118 is just above alpha, so it also has no call |
| 12 | Route B cost (keeps ≥ 95 % of A's genes) | **pass** | Route B keeps 2,875 / 2,942 route-A genes: 97.7 % [97.1, 98.2]. 14 genes are detected only in B | DB, library and wall times are in "Route A vs route B" |
| 13 | Variants | **not testable** (pre-registered) | `variants.tsv`: 0 rows in A, 0 in B | Same-genome search, no de novo arm |
| 14 | Second round | **not testable** | No ORF called in row 10, so no `makedb --patch` round | See Deviations |
| 15 | Circular topology | **pass** | Run `replicons.parquet`: NC_000962.3, 4,411,532 bp, `circular`, chromosome. 0 genes with `crosses_origin` | The deviation took effect. dnaA is LNBHIF_00001 at 1-1,524 (+), detected, q 3.7e-4, 34 unique peptides |

Counts: **10 pass, 1 fail** (row 3), **4 not testable** (rows 9, 10, 13 and 14; row 13 was pre-registered as not testable).

## Key numbers

| | Route A (called genes) | Route B (+ six-frame ORFs ≥ 20 aa) |
|---|---|---|
| DIA-NN precursors at 1 %, per run (3 runs) | 31,130 / 24,264 / 27,358 | 29,800 / 23,101 / 25,594 |
| DIA-NN protein groups, global q ≤ 0.01 | 2,926 | 2,983 |
| GIPPA coding genes / detected (`decoy_picked`) | 4,075 / 2,942: 72.2 % [70.8, 73.6] | 4,075 / 2,889: 70.9 % [69.5, 72.3] |
| Gene FDR at q ≤ 0.01 (decoy estimate) | 0.0099 | 0.0066 |
| Detected with ≥ 2 peptides | 2,738 | 2,677 |
| Novel-ORF candidates / called | - | 247 / 0 |
| Warnings | `novel_orfs.decoys_not_makedb`, `novel_orfs.layer_missing` (`rbs_novel`) | `novel_orfs.class_q_floor_above_alpha` |

`gippa compare` A vs B:
- Abundance: Spearman 0.989 [0.986, 0.992] over 2,852 shared genes (1,954 TUs), 1,000 bootstraps.
  Bland-Altman bias 0.076 log2, limits of agreement −1.12 to 1.27.
- Enriched categories: 141 in A vs 116 in B, 107 shared (Jaccard 0.71).
- Islands, operon flags, novel ORFs and alt sites: both empty.

## Exploratory (not pre-registered)

**Regions.** One region shows a significant observed-minus-expected deficit: CRISPR region
`crispr_2` (circular-shift test; effect −1.59, q 0.022 in A; −1.56, q 0.027 in B). No row
pre-registered it.

**Novel-ORF candidates vs RefSeq (all 247, uncalled).** The rule is the same as row 10: ≥ 50 %
same-strand overlap with a RefSeq CDS or pseudogene.
- All candidates: 127 / 247 match, 51.4 % [45.2, 57.6].
- By class:
  - nterm_extension: 93 / 109, 85.3 % [77.5, 90.8]
  - cds_outframe: 20 / 36, 55.6 % [39.6, 70.5]
  - intergenic: 12 / 26, 46.2 % [28.8, 64.5]
  - antisense: 2 / 76, 2.6 % [0.7, 9.1]
- Candidates with ≥ 2 peptides: 39 / 47 match, 83.0 % [69.9, 91.1]. 0 pseudogene matches.
- The rule counts every candidate without a start codon as a non-match, so these shares are low by construction.

**RefSeq CDS without a Bakta call** (no Bakta CDS shares the stop):
- 117 such CDS, 17 of them ≤ 70 aa.
- 22 have ≥ 1 exact route-B peptide, 4 of them ≤ 70 aa: L36 Rv3461c (37 aa, 2 peptides),
  L34 Rv3924c (47 aa, 1), VapB16 Rv2231B (58 aa, 1) and Rv2256a (64 aa, 1).
- 21 overlap a route-B candidate.

**Candidates with ≥ 2 peptides that coincide with a RefSeq CDS Bakta missed** (none called; for Khoa's review, not confirmed):

| ORF | class | coordinates (strand), aa | peptides | best peptide q | class FDR | RefSeq |
|---|---|---|---|---|---|---|
| NOVEL_NC_000962.3_3880286_m | intergenic | 3,880,286-3,880,399 (−), 37 | 2 | 5.3e-7 | 0.10 | Rv3461c, 50S ribosomal protein L36 |
| NOVEL_NC_000962.3_4410412_m | intergenic | 4,410,412-4,410,789 (−), 125 | 4 | 4.8e-7 | 0.10 | Rv3923c RNase P protein (also spans Rv3922c and Rv3924c L34) |
| NOVEL_NC_000962.3_2321057_p | cds_outframe | 2,321,057-2,321,467 (+), 136 | 7 | 4.8e-7 | 0.091 | Rv2063A, MazF7 (inside Bakta LNBHIF_02193, other frame) |
| NOVEL_NC_000962.3_3182268_p | cds_outframe | 3,182,268-3,183,395 (+), 375 | 9 | 4.8e-7 | 0.091 | Rv2871, VapB43 (85 aa) |
| NOVEL_NC_000962.3_2181262_m | intergenic | 2,181,262-2,182,560 (−), 432 | 5 | 5.3e-7 | 0.091 | Rv1929c (214 aa) |
| NOVEL_NC_000962.3_4116002_p | cds_outframe | 4,116,002-4,116,379 (+), 125 | 4 | 4.1e-4 | 0.091 | Rv3675, membrane protein |
| NOVEL_NC_000962.3_1484982_p | intergenic | 1,484,982-1,485,278 (+), 98 | 3 | 1.9e-4 | 0.10 | Rv1322 |
| NOVEL_NC_000962.3_2808083_p | intergenic | 2,808,083-2,808,304 (+), 73 | 2 | 5.9e-7 | 0.10 | Rv2493, VapB38 |

All 47 candidates with ≥ 2 peptides are listed in `metrics.json` › `candidates_ge2_peptides`. The largest
nterm_extension candidates are NOVEL_NC_000962.3_9899_p (11 peptides, upstream of LNBHIF_00008) and
NOVEL_NC_000962.3_4353010_p (9 peptides, LNBHIF_04103). Both have class FDR 0.0118.

## Route A vs route B

| | Route A | Route B | B / A |
|---|---|---|---|
| makedb records (incl. 381 `CONT_`) | 4,456 | 102,820 (+ 98,364 six-frame ORFs) | 23.1x |
| FASTA size | 1,719,817 bytes | 13,850,157 bytes | 8.1x |
| makedb wall time | 3 s | 9 s | |
| predicted library | 634,544 precursors, 155,204,874 bytes | 4,300,733 precursors, 1,056,212,594 bytes | 6.8x |
| DIA-NN step 1 / step 2 (3 runs, 16 CPUs) | 145 s / 1,251 s | 578 s / 3,414 s | 4.0x / 2.7x |
| DIA-NN total | 1,405 s | 4,003 s | 2.8x |
| convert + preflight | 21 s | 15 s | |
| `gippa run` (`farmshare_local`, 14 threads, 54 GB) | 19,045 s (5 h 17 min), wheat-01 | 11,515 s (3 h 12 min), wheat-05 | 0.60x |
| whole job | 20,699 s (job 2005884) | 15,542 s (job 2007057) | |

Route A's `gippa run` took longer than route B's. The two ran on different nodes and overlapped in
time: A ran 07:11-12:29 UTC, B ran 11:07-14:19 UTC. The cause was not measured, so these wall times do
not establish a cost of route B in GIPPA.

## Command log

Everything ran from the pin (`bfc74ff6`) with `PYTHONPATH=<pin>/src`. `$D` = `$GIPPA_ROOT/campaign/PXD057335`, label `mtb_h37rv_circ`.

```bash
sbatch --array=1 p2c/genome_bakta.sbatch            # 2005457: Bakta 1.12.2 rerun on the circular-tagged genome
sbatch -J gippa-camp-p2c-<ds>-a route.sbatch PXD057335 mtb_h37rv_circ a   # 2005884 (wheat-01)
sbatch -J gippa-camp-p2c-<ds>-b route.sbatch PXD057335 mtb_h37rv_circ b   # 2007057 (wheat-05)
# inside route.sbatch, per route:
gippa makedb --annotation $J --genome $G --outdir $D/db_a                  # route b: + --six-frame-orfs --min-orf-aa 20
bash p2c/diann_search.sbatch $D/diann_<r> $D/db_<r>/strain_db.faa $D/runs.txt   # flags p2c/diann_campaign.cfg; exit 0
gippa convert --from diann ... --search-db same    # twice, byte-identical (deterministic true); exit 0 both routes
gippa preflight ... --strict                        # exit 0 both routes
gippa run --profile farmshare_local,apptainer --genome $G --annotation $J \
  --peptides $D/in_<r>/peptides.tsv.gz --db-dir $GIPPA_DB_DIR --threads 14 --memory 54.GB --outdir $D/out_<r>  # exit 0
# post (session allocation, job 2005379, iron-04; p2c/post.sbatch steps):
gippa compare $D/out_a $D/out_b --labels genes,sixframe --outdir $D/cmp_ab  # exit 0, 9 s
gippa report $D/out_a --redact                                             # exit 0, 75 s, tables/ unchanged (sha256)
gippa report $D/out_b --redact                                             # exit 0, 100 s, tables/ unchanged (sha256)
```

Both runs: all modules (`modules_off` empty), master seed 20261014, `status` completed, exit 0.

| item | value |
|---|---|
| GIPPA | `bfc74ff6c20acc2cfcf298cdf23cbb30d8fae0dc` (1.1.1), clean |
| genome (circular-tagged) | sha256 `03b3edadc33561efcb071f0b69d3836c0665509912138e6ffea7d14aef35ae04` |
| Bakta JSON (rerun 2005457) | sha256 `a17257269010230bf5f42e1ed1e2c26c1bac1067f24bed90ed7f9a1b4020b76c` |
| peptides A / B | `443d101b502eb0c70c907671ce8cf6d554d4f5a86bbd64bad3a96d9ca432e21a` / `6eff66c6578233064481d303bf772925bfc15dec9fefab04bb98789eaf7628bc` |
| abundance A / B | `45066481cf54fb3a3a898c0f2ae90eeab10beb4d8f942d63b5cb36d1eb65058b` / `8afc6b128d4589054510aa954715123eea238eb1ae23383788f636f2c67614e0` |
| DIA-NN report A / B | `ec07e697ad9c13796e5eec524e627d45f51cdff393b4f0c6d080dcf204625560` / `53f59d7b87f703c4b7844974e607568080cd8dffca7aa38b8bb968ba3fcc0f71` |
| strain_db.faa A / B | `3184cb41b9c4cec525450df332f528f7af60742756f4929124de3a2b3853e42d` / `4bd46d488643eead7b5ccfa3a7107359c352976e3d938be68732aadeca8f444c` |
| DIA-NN | 2.3.1 Academia, image sha256 `1946452d56d2fb0f66e05ced473df88e53484e28040f1ecbd83c4a43405a6381`; cfg sha256 `dad24b38e315b68cb6a905dfda877404535310f030ed9f92f8e9370dae0d6a1f` |
| raw `_DMSO_01` / `_02` / `_03` | `9324e626ab3c…` / `0b13121d2b75…` / `4880cd9c1e3a…` (full sha256 in `provenance.json`; size-verified, no deposit checksum) |
| antiSMASH 8.0.4 | `sha256:d3f390e687cfe02f1387be40256dddbf678e91bc56b96c36da07b56ec8bda29d` |
| geNomad 1.12.0 | `sha256:d54652ff868af72604ea0e54c541c5c92168c7ac821a422de05dcbad65509e5c` |
| OSTIR 1.1.3 / ViennaRNA 2.6.4 | `sha256:d9107b9b87fd6f3d3cd329ba1d2ac66b9af9e38512f11d2b8a0210110c6c9292` |
| Promotech 56251ad9 | `sha256:7366eaf67d9b26c6036ed5ececa740fe2728fdcf0cd4b2a12a0d9597b87a53ac` |
| Promoter Calculator 1.2.4 | `sha256:c5bd0662ceed623c00c4bf40b085c6541b7676e1ccc37a390e0efd73e02b48df` |
| TerminatorNet 72180309 | `sha256:ec415757b6b70cce8240f425b35bf3e74dc3abe8ff413049c04d808034592f7a` |
| TransTermHP 2.09 | `sha256:fc6fe0404be9340a883ec7a7b90415dfc3826c2bca7333a0314335c7bfee2d03` |
| UniOP 3ed14dfc | no container digest recorded |

## Deviations and problems

- **Circular tag and Bakta rerun (lead decision, 2026-11-21).**
  - NCBI marks NC_000962.3 linear, so the recipe would have left it untagged. It was tagged
    `[topology=circular]` instead.
  - GIPPA reads topology from the precomputed Bakta JSON, and the prep JSON (job 2001987_13) says
    linear. Bakta was therefore rerun on the tagged genome with the same image and argv (job 2005457).
  - A comparison by the group C agent found one difference from the prep annotation: dnaA starts at 1
    instead of 34. The rerun's dnaA (1-1,524 +) is confirmed in `metrics.json`; the "34" comes from
    that comparison, not from `metrics.json`.
  - Row 15 confirms that the run used circular topology.
- **`gippa compare` and `report --redact` ran in the session allocation** (job 2005379, iron-04), not
  as a separate sbatch.
  - Wall times: compare 9 s, reports 75 s and 100 s (about 3 min in total).
  - Both runs' `tables/` were unchanged by `report --redact` (sha256 checked before and after).
- **TE7 determinism rerun not done** for this dataset. The answer key does not list TE7, so no verdict
  is missing.
- **No second novel-ORF round**: route B called nothing (row 14).
- **Measured differently from the pre-registration:**
  - Row 8: 150 PE/PPE Bakta genes by shared stop, not the key's 155 RefSeq CDS. Details in the row.
  - Row 11: nterm_extension has 109 candidates, slightly above the key's "< ~100".
- **Contaminants** were appended by makedb as `CONT_` records (381), as the key fixed. DIA-NN got one FASTA.
- **Route A warning** `novel_orfs.layer_missing` (`rbs_novel`) and `decoys_not_makedb`: route A has no
  six-frame ORFs, as expected.
