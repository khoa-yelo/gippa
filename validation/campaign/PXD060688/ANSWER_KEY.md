# PXD060688 answer key: Myxococcus xanthus DK1622 (pre-registered before any search)

Campaign Phase 2, group B (Khoa, 2026-11-21: "run phase 2 with 1.1.1"). Committed and pushed
before any spectrum of this dataset is searched; only `gippa makedb` and the in-silico DIA-NN
library prediction (no spectra) may precede it. A result found only after looking is reported as
exploratory. Thresholds come from the handover plan, GIPPA_SPEC.md §11.6 and the Phase 1 answer
key, and are not loosened. Expectations marked *prior* are the agent's own, written down so they
cannot move; they are not gates. Aggregates only in git: no peptide sequences, no DIA-NN rows.
Coordinates are 1-based inclusive.

## Fixed inputs

| Item | Value |
|---|---|
| Spectra | PXD060688 (PMC12288472), Orbitrap Exploris 480 DIA, `DdiA-4-01..04.raw` (DK1622 wild type, untreated, presumed exponential; no file-to-sample key in the deposit), SHA-1 verified (`PREP.md`) |
| Genome | GCF_000012685.1 (ASM1268v1): NC_008095.1, 9,139,763 bp, circular-tagged, GC 68.9 %; `mxan_dk1622.fna.gz` sha256 `a36130a6…` |
| Annotation | Bakta 1.12.2, DB 6.0 full (job 2001987_14): 7,326 CDS, 0 sORF; GBFF kept for antiSMASH |
| Answer-key annotation (never an input) | RefSeq `GCF_000012685.1_ASM1268v1_genomic.gff.gz` (7,234 CDS, 72 pseudogenes) |
| GIPPA code | GIPPA 1.1.1 release, detached worktree `local/worktrees/rel111-pin` at `bfc74ff6`, `PYTHONPATH=<pin>/src`, venv Python (never `uv run`) |
| Search | DIA-NN 2.3.1 Academia, pinned image (`validation/pilot/diann_image.sha256`), the two calls and flags of `validation/pilot/diann_A2.cfg` (`validation/campaign/p2b/diann_campaign.cfg`, identical flags; `--report-decoys` kept; no `--met-excision`, `--reanalyse`, `--rt-profiling`); Thermo `.raw` read directly (ADR-0077 Amendment 1); FASTA = the `gippa makedb` database, which appends the contaminant library as `CONT_` records |
| Run | `gippa run` all modules, default seed 20261014, `--threads 14 --memory 54.GB`, profile `farmshare_local,apptainer` inside one 16-CPU Slurm allocation |

## Loci (Bakta 1.12.2 locus tags)

| Gene | Locus | Coordinates (strand) | Note |
|---|---|---|---|
| rplJ / rplL | BFAOEK_03062 / BFAOEK_03063 | 3,597,887-3,598,408 / 3,598,512-3,598,880 (+) | S3 |
| atpA / atpD | BFAOEK_06952 / BFAOEK_06849 | 8,601,722-8,603,269 (+) / 8,474,315-8,475,760 (−) | S3 |

## Standard checks (every Phase 2 dataset, same rules as Phase 1)

Route A = called genes (`gippa makedb`), route B = called genes + six-frame ORFs ≥ 20 aa
(`--six-frame-orfs --min-orf-aa 20`). Both searched with DIA-NN 2.3.1 on the same runs, converted
with `--search-db same`, `preflight --strict` exit 0, `gippa run` with all modules and the default
seed, then `gippa compare A B --labels genes,sixframe`. Verdict rules are fixed here and computed
by `validation/campaign/p2b/p2b_metrics.py`; the TE functions are the pilot's
(`validation/pilot/metrics.py`, `te_real.py`) at the pinned commit. Every share carries a Wilson
95 % CI; differences of shares a Newcombe (hybrid Wilson) 95 % CI.

| # | Check | Expected (threshold) | Verdict rule | Read from |
|---|---|---|---|---|
| S1 | Mapping, TE1 | ≥ 99.5 % of route-A accepted peptides map exactly in frame | pass/fail (`te1_real`) | A `mapping/`, `mapping_fdr.tsv` |
| S2 | Expression, TE3 | ribosomal proteins ≥ 90 % detected and their median abundance ≥ the 75th percentile of measured genes | pass/fail (`te3`) | A `genes.tsv` |
| S3 | Stoichiometry, TE4 | log2(RplL/RplJ) in [1, 3]; log2(AtpA/AtpD) in [−1, 1] (molar abundance, route A), on the locus tags named in "Loci" below (the F0F1 copies where a name is duplicated). The pair failed on the Nissle pilot and in Phase 1 (RplL/RplJ, waived for 1.1.0, ADR-0076): a known risk | pass if both pairs inside; fail otherwise; `te4` (name-based) reported next to it | A `genes.tsv` |
| S4 | Contaminants, TE5 | trypsin and keratin peptides land in the contaminant bucket, none placed on the genome | pass/fail (`te5_arm`) | A `unmapped.tsv`, `peptide_hits.tsv` |
| S5 | Gene-level FDR | decoy-estimated FDR of the detected set ≤ 0.01 in A and in B (tier `decoy_picked`) | pass/fail | `detection_summary.tsv` (tier, q 0.01, fdr) |
| S6 | Genome enrichment | every geNomad prophage region of route A has observed < expected detection and q ≤ 0.05 (`circular_shift_o_minus_e`), the Phase 1 rule; islands listed with q | pass if all; fail if any O ≥ E or any q > 0.05; not testable if no prophage region | A `regions_stats.tsv`, `islands.tsv` |
| S7 | Operons | no adjacent ribosomal-protein or ATP-synthase gene pair with within-operon q ≤ 0.05 in `pairs.tsv` or a non-concordant `operon_flags.tsv` pattern | pass/fail | A `pairs.tsv`, `operon_flags.tsv` |
| S8 | Novel ORFs vs RefSeq | each called route-B ORF compared with the RefSeq GFF (answer key, never an input): matches when a RefSeq CDS or pseudogene on the same strand overlaps ≥ 50 % of its length; each called ORF classed as known in RefSeq, known in the dataset's paper, or new (candidate for Khoa) | `reported` with Wilson CI; not testable if 0 called | B `novel_orfs.tsv`, `novel_orfs.gff3` |
| S9 | Class q floor | a novel-ORF class with too few candidates reports `novel_orfs.class_q_floor_above_alpha`; alpha 0.01 is not loosened | pass if every class with min class FDR > 0.01 and < 100 candidates is covered by the warning (or no such class) | B `run_manifest.json` warnings, `novel_orfs.tsv` |
| S10 | Route B cost | share of route A's detected genes that route B keeps; DB records and bytes; library precursors and bytes; DIA-NN and GIPPA wall times. *Prior* (not a gate): kept share ≥ 0.95 | `reported` | `gippa compare`, provenance |
| S11 | Second round | only if route B calls ≥ 1 ORF: `makedb --patch out_b/annotation_patch.gff3`, search again; every patched protein detected at gene-level q ≤ 0.01 | pass/fail; not run if 0 calls | round-2 `genes.tsv` |
| S12 | Variants | a same-genome database search cannot produce strain variants: expected 0 `high` ALT rows in `variants.tsv` of A and B (any row is listed with its flags, `mass_equivalent` separately) | pass if 0 `high` ALT outside start/stop codons; `variants` module skip recorded | A/B `variants.tsv`, manifest |

## Dataset checks (handover plan, Phase 2 item 8)

| # | Check | Expected | Verdict rule | Read from |
|---|---|---|---|---|
| M1 | Library size first | the in-silico DIA-NN libraries of routes A and B are predicted and measured (precursors, bytes, wall time, peak memory) before the spectra are searched; route B is the largest six-frame database of the campaign. The search proceeds unless step 1 fails or exceeds the 64 GB allocation | `reported` | `diann_a/`, `diann_b/` step-1 logs and provenance |
| M2 | Detection inside antiSMASH regions (`bgc` module) | per region: genes, detected, O/E and q (`regions_stats.tsv`, `region_type` bgc); pooled detected share inside vs outside regions with a Newcombe CI. *Prior*: the pooled share inside is lower than outside (many clusters are silent in lab culture), and at least one region (e.g. a known vegetative-phase cluster) has detected genes | `reported` (no plan threshold); the prior is stated for honesty | A `regions_stats.tsv`, `genes.tsv` (`region_ids`) |
| M3 | antiSMASH ran | `layers_missing` does not list antiSMASH/bgc in A or B | pass/fail | `run_manifest.json` |

S13 (determinism, TE7) is run on PXD072227 in this group, not here.
