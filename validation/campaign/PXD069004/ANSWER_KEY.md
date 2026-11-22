# PXD069004 answer key: Caulobacter crescentus NA1000 (pre-registered before any search)

Campaign Phase 2, group A (handover plan of 2026-11-17, dataset 3; Khoa 2026-11-21 "run phase 2
with 1.1.1"). Committed before any spectrum of this dataset is searched. A result found only after
looking is reported as exploratory. Aggregates only in git.

## Fixed inputs

| Item | Value |
|---|---|
| Spectra | PXD069004 (eLife 2026, PMC13046382, Richter et al., peptidoglycan recycling), Orbitrap Exploris 480 DIA, `AmiX-5-4/5/6.raw` = wild type in M2G without GlcNAc (deposit `AmiX.xlsx`), SHA-1 verified. DIA-NN reads the Thermo `.raw` files directly (Phase 0 raw check, ADR-0077) |
| Genome | GCF_000022005.1 (NA1000): chromosome NC_011916.1 (4,042,929 bp, circular, GC 67.2 %) |
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
| C1 | PG recycling enzymes in wild type | the paper shows a functional PG recycling pathway in wild-type cells: ≥ 6 of its 10 Table 1 enzymes (GlmM CCNA_00116, AmpG 00136, AnmK 01945, NagZ 02085, GlmU 02389, MupP 02390, AmpD 02650, AmgK 03649, MurU 03650, NagK 03849; Bakta loci in `answer_loci.tsv`) detected at gene q ≤ 0.01 in route A → pass | `genes.tsv` |
| C2 | GlcNAc clusters uninduced | the two clusters CCNA_00452-00458 and CCNA_00566-00574 (16 genes) are strongly induced by GlcNAc, which these runs lack: pass if ≤ 4 of the 16 are detected, or the median abundance percentile (among route-A detected genes, `abundance_log2`) of the detected cluster genes is ≤ 0.50 | `genes.tsv` |
| C3 | Depth | route-A detected genes / Bakta CDS+sORF (3,724) with Wilson CI. *Prior*: ≥ 0.50 (3 Exploris DIA runs; Phase 1 E. coli reached 0.62) | `genes.tsv`; `reported` |
| C4 | High-GC genome | the run completes with the default sharded TerminatorNet path (M2-15b) and `layers_missing` lists no TerminatorNet layer → pass | `run_manifest.json` |
| C5 | Annotation vs RefSeq | Bakta 3,724 vs RefSeq 3,886 CDS (ratio 0.958, `PREP.md`): route-B candidates that match RefSeq CDSs Bakta did not call are listed as known | `novel_orfs.tsv` vs RefSeq; `reported` |
