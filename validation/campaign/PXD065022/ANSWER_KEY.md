# PXD065022 answer key: Treponema pallidum, 7 strains (pre-registered before any search)

Campaign Phase 2, group A (handover plan of 2026-11-17, dataset 1; Khoa 2026-11-21 "run phase 2
with 1.1.1"). Committed before any spectrum of this dataset is searched. A result found only after
looking is reported as exploratory. Aggregates only in git.

## Fixed inputs

| Item | Value |
|---|---|
| GIPPA code | detached worktree `local/worktrees/rel111-pin` at `bfc74ff6` (1.1.1 release code, not yet published; provenance says 1.1.1), `PYTHONPATH=<pin>/src`, venv Python, never `uv run` |
| Spectra | PXD065022 (J Proteome Res 2025, PMC12687308), timsTOF Pro DIA, 3 replicates per strain (A, B, C) for DAL-1, Grady, Haiti B, Madras, Mexico A, Philadelphia 1, SS14 (21 `.d` runs; isolates 171 and 176 skipped, Khoa question 22). sha256 in `files.tsv` / `PREP.md` |
| Genomes | each strain on its own RefSeq assembly (`genomes.tsv`): DAL-1 GCF_000246815.1, Grady GCF_005408405.1, Haiti B GCF_003627115.1, Madras GCF_024971635.1, Mexico A GCF_000304295.1, Philadelphia 1 GCF_005408385.1, SS14 GCF_000410555.1; one circular chromosome each, tagged `[topology=circular]` |
| Annotation | Bakta 1.12.2, DB 6.0 full, `--keep-contig-headers` (job 2001987; counts in `PREP.md`) |
| Host | rabbit reference proteome UniProt UP000001811 (taxon 9986, current release), MD5 checked against the release `RELEASE.metalink`; headers rewritten to `CONT_RABIT_<accession>` so host peptides fall in GIPPA's contaminant bucket through the default `input.contaminant_prefixes` (`CONT_`). No GIPPA setting is changed. The paper searched 21,176 rabbit sequences |
| Search | DIA-NN 2.3.1 Academia (pinned image `1946452d…`), the two calls and flags of `validation/pilot/diann_A2.cfg` (`validation/campaign/p2a/diann_campaign.cfg`, flags identical); `--report-decoys` kept; no `--met-excision`, `--reanalyse`, `--rt-profiling`. FASTA = the strain's `gippa makedb` database (contaminants appended as `CONT_`) + the rabbit FASTA. One search per strain and route (3 runs each) |
| Runs | per strain: route A and route B, `convert --search-db same`, `preflight --strict` exit 0, `gippa run` all modules, default seed, `--threads 14 --memory 54.GB`, profile `farmshare_local,apptainer` inside a 16-CPU job; `gippa compare A B` per strain |
| Answer-key loci | `answer_loci.tsv` (built by `p2a/answer_loci.py` from the RefSeq GFFs, the GenBank feature tables CP003115.1 and CP035104.1 and the Bakta JSONs, before any search) |

## Standard checks (every Phase 2 group A dataset)

Route A = `gippa makedb` called genes; route B = called genes + six-frame ORFs ≥ 20 aa. Both runs
use the same DIA-NN settings and the same spectra. Thresholds come from GIPPA_SPEC.md §11.6 and the
handover plan and are not loosened. Verdicts: pass, fail, not testable, or `reported` (no plan
threshold). Shares carry a Wilson 95 % CI; anything without an error estimate says `unestimated`.

| # | Check | Expected (verdict rule) | Read from |
|---|---|---|---|
| S1 | Mapping (TE1) | ≥ 99.5 % of route-A accepted target peptides map exactly in frame at their own gene (pilot rule, `validation/pilot/metrics.py::te1_real`) | `mapping/peptide_hits.tsv`, `unmapped.tsv` |
| S2 | Expression (TE3) | ribosomal proteins ≥ 90 % detected and their median abundance in the top quartile (`te_real.te3`) | `genes.tsv` |
| S3 | Stoichiometry (TE4) | log2(RplL/RplJ) in [1, 3]; log2(AtpA/AtpD) in [−1, 1] (`te_real.te4`); a pair whose genes are not both named is not testable. RplL/RplJ failed in the pilot and in Phase 1: a known risk | `genes.tsv` |
| S4 | Contaminants (TE5) | trypsin and keratin peptides land in the contaminant bucket, none placed on the genome (`te_real.te5_arm`) | `unmapped.tsv` |
| S5 | Operons | no adjacent ribosomal-protein or ATP-synthase gene pair with within-operon q ≤ 0.05 in `pairs.tsv`, and no `operon_flags.tsv` row among them with a non-concordant abundance pattern | `pairs.tsv`, `operon_flags.tsv` |
| S6 | Genome context | significant regions (q ≤ 0.05), islands and operon flags listed with q-values | `regions_stats.tsv`, `islands.tsv`, `operon_flags.tsv`; `reported` |
| S7 | Novel ORFs vs RefSeq | every called route-B ORF compared with the RefSeq GFF (≥ 50 % same-strand overlap with a CDS or pseudogene = known); share known with Wilson CI; each call also classed as known in the paper or new. 0 calls → not testable | `novel_orfs.tsv` vs RefSeq GFF (answer key, never an input) |
| S8 | Class q floor | a class with too few candidates for q 0.01 carries `novel_orfs.class_q_floor_above_alpha`; alpha stays 0.01 → pass when the warning is present wherever the class minimum FDR > 0.01 and no class is called below its floor | `run_manifest.json` warnings, `novel_orfs.tsv` |
| S9 | Route B cost | share of route A's detected genes that route B keeps (Wilson CI); database records and bytes; library precursors; DIA-NN and `gippa run` wall times | `gippa compare` output, provenance; `reported` |
| S10 | Gene-level FDR | detected genes by tier and the run's decoy-estimated gene-level FDR at q ≤ 0.01 | `genes.tsv`, `detection_summary.tsv`; `reported` |
| S11 | Second round | only if route B calls ≥ 1 ORF: `makedb --patch`, search again, every patched protein detected at gene-level q ≤ 0.01 → pass; no call → not run | round-2 `genes.tsv` |

S1-S11 are judged per strain; the dataset verdict of a pass/fail row is pass only if every strain
passes (the per-strain numbers are listed).

## Dataset checks

| # | Check | Expected (verdict rule) | Read from |
|---|---|---|---|
| T1 | Host share | the rabbit host dominates: in every strain, ≥ 50 % of the distinct accepted (q ≤ 0.01) target peptides of route A carry only `CONT_` proteins (rabbit + contaminant library). The paper removed on average 9,592 rabbit and 557 cRAP peptides per sample against 928-6,867 treponemal ones | route-A `peptides.tsv.gz` (aggregate only) |
| T2 | Detection per strain | (a) in every strain, route-A detected genes / Bakta CDS+sORF in [0.50, 0.95] (paper: 62.4-88.9 % of annotated proteins with all annotations combined; 605-863 proteins per strain with RAST); (b) Grady has the fewest detected genes of the seven (paper: Grady lowest) | `genes.tsv` |
| T3 | 12 newly annotated genes with peptides in the paper | TP0005a, TP0134a, TP0126d, TP0127a, TP0206a (rpmD), TP0409a (yajC), TP0453a, TP0913a, TP0921a, TP0950a (rnpA), TP0954a on DAL-1 and TP0250b (rpsT) on Grady (`answer_loci.tsv`). A gene is recovered when its same-stop Bakta CDS is detected (gene q ≤ 0.01) in route A or B of that strain; for TP0134a, which has no same-stop Bakta CDS, when ≥ 1 accepted route-B peptide lies inside 156,669-157,301 (+) in that frame. Pass if ≥ 9 of 12 are recovered (Wilson CI reported) | DAL-1 and Grady `genes.tsv`, `peptide_hits.tsv`, `novel_orfs.tsv` |
| T4 | Translated pseudogene distal parts | TP0865a (Haiti B, Bakta FOKJEF_00855, 944,868-945,575 −) detected in Haiti B, and TP0924a (Madras, Bakta BIIDIP_00910, 1,006,278-1,006,985 +) detected in Madras, each at gene q ≤ 0.01 in route A → pass; one of two → fail | Haiti B / Madras `genes.tsv` |
| T5 | Nichols-like vs SS14-like abundance | for TP0477 (pgl), TP0487, TP0597 (higher in Nichols-like DAL-1, Haiti B, Madras) and TP0594, TP0732 (folD) (higher in SS14-like SS14, Mexico A, Grady, Philadelphia 1): the sign of mean(Nichols) − mean(SS14) of route-A `abundance_log2`, each strain median-centred on the genes detected in all seven strains (orthologs from the T6 reciprocal best hits), agrees with the paper for ≥ 4 of the 5 proteins. A protein detected in fewer than 2 strains of either group is not testable and leaves the denominator; < 3 testable → not testable | route-A `genes.tsv`, ortholog table |
| T6 | Cross-strain variants on DAL-1 (Phase 1b design) | each of the six other strains' route-A DIA-NN report converted again with `--search-db other` and run on DAL-1 with `--genome-relation related`. Truth set: Bakta CDS/sORF proteins of the strain vs DAL-1, blastp (pinned GIPPA BLAST image) both ways, `-seg no -evalue 1e-10 -max_hsps 1`, reciprocal best hits with identity ≥ 90 % and coverage ≥ 80 % of both; every differing non-gap, non-I/L column is a site (`p2a/truth_set.py`, as Phase 1b). Scored as Phase 1b (`ANSWER_KEY_1b.md` rules): recall on covered sites (`reported`; *prior* low, because Phase 1b showed a 1-difference stratum with < ~100 targets accepts nothing), precision of `high` ALT (`reported`), every non-truth `high` ALT listed and classed (`mass_equivalent`, start/stop, unexplained), unplaced share (`reported`), placement FDR ≤ 1 % by the pooled TE8 estimator → **pass/fail** per strain | cross-strain `variants.tsv`, `peptide_hits.tsv`, `mapping_fdr.tsv`, truth set |
| T7 | Determinism (TE7, one per phase for group A) | DAL-1 route A rerun with the same inputs and seed: `tables/` byte-identical | `diff -r` |

Also reported, not gated: novel-ORF candidates in all strains against the paper's list of 100
letter-suffix genes in the CP003115.1 and CP035104.1 feature tables (the paper says 102), the
`novel_orfs.class_q_floor_above_alpha` warning expected on a 1.14 Mb genome, and the route-B cost
with a rabbit FASTA in the search space.
