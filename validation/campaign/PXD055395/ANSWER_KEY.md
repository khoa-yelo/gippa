# PXD055395 answer key: Bacteroides thetaiotaomicron VPI-5482 (pre-registered before any search)

Campaign Phase 2, group A (handover plan of 2026-11-17, dataset 2; Khoa 2026-11-21 "run phase 2
with 1.1.1"). Committed before any spectrum of this dataset is searched. A result found only after
looking is reported as exploratory. Aggregates only in git.

## Fixed inputs

| Item | Value |
|---|---|
| Spectra | PXD055395 (PMC13253324), timsTOF Pro dia-PASEF, Evosep 100 SPD (11.5-min gradient), the 16 QC-passed B. thetaiotaomicron monoculture runs in mGAM (`run_selection.tsv`), SHA-1 verified against the deposit |
| Genome | GCF_000011065.1 (VPI-5482; the deposit's DIA-NN FASTA is UP000001414 = ATCC 29148 / VPI-5482): chromosome NC_004663.1 (6,260,361 bp) and plasmid p5482 NC_004703.1 (33,038 bp), both circular |
| GIPPA code | detached worktree `local/worktrees/rel111-pin` at `bfc74ff6` (1.1.1 release code, not yet published; provenance says 1.1.1), `PYTHONPATH=<pin>/src`, venv Python, never `uv run` |
| Annotation | Bakta 1.12.2, DB 6.0 full, `--keep-contig-headers` (job 2001987; counts in `PREP.md`) |
| Search | DIA-NN 2.3.1 Academia (pinned image `1946452d…`), the two calls and flags of `validation/pilot/diann_A2.cfg` (`validation/campaign/p2a/diann_campaign.cfg`, flags identical); `--report-decoys` kept; no `--met-excision`, `--reanalyse`, `--rt-profiling`. FASTA = the `gippa makedb` database (contaminants appended as `CONT_`); one search per route over all selected runs |
| Runs | route A and route B, `convert --search-db same`, `preflight --strict` exit 0, `gippa run` all modules, default seed, `--threads 14 --memory 54.GB`, profile `farmshare_local,apptainer` inside a 16-CPU job; `gippa compare A B` |
| Answer-key loci | `answer_loci.tsv` (built by `p2a/answer_loci.py` from the RefSeq GFF and the Bakta JSON, before any search) |

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

## Dataset checks

| # | Check | Expected (verdict rule) | Read from |
|---|---|---|---|
| B1 | Depth with short gradients | route-A detected genes / Bakta CDS+sORF (4,941) in [0.25, 0.70]. The paper quantified on average > 1,300 proteins per species (37 % of the genome; ≥ 2 peptides in ≥ 2 replicates) with DIA-NN 1.8.2 | `genes.tsv` |
| B2 | Starch utilisation (Sus) in mGAM | mGAM holds 12 g/l starch and the paper's ΔsusG mutant changed the phenotype in mGAM, so the Sus operon is expressed: ≥ 4 of susC, susD, susE, susF, susG (BT3702-BT3698, Bakta BEBDPP_03874-03870) detected at gene q ≤ 0.01 in route A → pass | `genes.tsv` |
| B3 | Plasmid p5482 | detected share of plasmid genes vs chromosome genes, with Wilson CIs. *Prior*: plasmid share lower than the chromosome's | `genes.tsv`; `reported` |
| B4 | Replicate depth | per-run precursor counts at 1 % (DIA-NN log) and the paper's QC floor of 4,000 precursors per run: every run ≥ 4,000 → pass | DIA-NN step-2 log (counts only) |
| B5 | Annotation vs RefSeq | Bakta CDS+sORF 4,941 vs RefSeq 4,687 CDS (ratio 1.054, `PREP.md`); novel-ORF candidates that overlap RefSeq pseudogenes (114) listed | `PREP.md`, `novel_orfs.tsv`; `reported` |
