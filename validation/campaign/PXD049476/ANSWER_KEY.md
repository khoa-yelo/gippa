# PXD049476 answer key: Vibrio natriegens ATCC 14048 (pre-registered before any search)

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
| Spectra | PXD049476 (PMC11066055), Orbitrap Exploris 480 DIA, `Vnat-1-1..4.raw` (presumed parental strain, exponential; Khoa, question 18), SHA-1 verified against the deposit (`PREP.md`) |
| Genome | GCF_001456255.1 (ASM145625v1): chromosome 1 NZ_CP009977.1 (3,248,023 bp) and chromosome 2 NZ_CP009978.1 (1,927,130 bp), both circular-tagged; `vnat_atcc14048.fna.gz` sha256 `abffc6c3…` |
| Annotation | Bakta 1.12.2, DB 6.0 full (job 2001987_10): 4,553 CDS, 0 sORF; JSON and GBFF together |
| Answer-key annotation (never an input) | RefSeq `GCF_001456255.1_ASM145625v1_genomic.gff.gz` (4,504 CDS, 19 pseudogenes) |
| GIPPA code | GIPPA 1.1.1 release, detached worktree `local/worktrees/rel111-pin` at `bfc74ff6`, `PYTHONPATH=<pin>/src`, venv Python (never `uv run`) |
| Search | DIA-NN 2.3.1 Academia, pinned image (`validation/pilot/diann_image.sha256`), the two calls and flags of `validation/pilot/diann_A2.cfg` (`validation/campaign/p2b/diann_campaign.cfg`, identical flags; `--report-decoys` kept; no `--met-excision`, `--reanalyse`, `--rt-profiling`); Thermo `.raw` read directly (ADR-0077 Amendment 1); FASTA = the `gippa makedb` database, which appends the contaminant library as `CONT_` records |
| Run | `gippa run` all modules, default seed 20261014, `--threads 14 --memory 54.GB`, profile `farmshare_local,apptainer` inside one 16-CPU Slurm allocation |

## Loci (Bakta 1.12.2 locus tags)

| Gene | Locus | Coordinates | Note |
|---|---|---|---|
| rplJ / rplL | KGNFJJ_00022 / KGNFJJ_00023 | NZ_CP009977.1 15,607-16,101 / 16,159-16,527 (+) | S3 |
| atpA / atpD | KGNFJJ_02722 / KGNFJJ_02724 | NZ_CP009977.1 2,927,838-2,929,379 / 2,930,326-2,931,729 (+) | S3: the F0F1 operon copy; a second `atpD` (KGNFJJ_04534, chromosome 2) makes the name-based `te4` ambiguous |

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

## Dataset checks (handover plan, Phase 2 item 4)

| # | Check | Expected | Verdict rule | Read from |
|---|---|---|---|---|
| V1 | Detection per chromosome | detected share of coding genes (CDS + sORF) on chromosome 1 > chromosome 2 (chromosome 2 carries more condition-specific and accessory genes) | pass if the difference chr1 − chr2 is > 0 and its Newcombe 95 % CI excludes 0; fail otherwise | A `genes.tsv` (`contig`) |
| V2 | Abundance per chromosome | median `abundance_log2` of measured genes on chromosome 1 ≥ chromosome 2; difference with a 95 % CI from 2,000 gene bootstrap draws (seed 20261014) | pass if the median difference > 0 and its CI excludes 0 | A `genes.tsv` |
| V3 | Route B per chromosome | share of route A's detected genes kept by route B, per chromosome | `reported` | A and B `genes.tsv` |
| V4 | Fused-chromosome strain (optional in the plan) | not in this pass: the selected runs are the parental block only | not run | - |

S13 (determinism, TE7) is run on PXD072227 in this group, not here.
