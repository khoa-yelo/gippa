# PXD065346 answer key: Clostridioides difficile 630Δerm (pre-registered before any search)

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
| Spectra | PXD065346 (PMC12850539), Q Exactive DIA, total lysate with the lowest spike-in (0.5 pg/µg), technical triplicate `221217_…_Ctrl_A`, `221219_…_Ctrl_B`, `230106_…_Ctrl_C` (Khoa, question 20), SHA-1 verified (`PREP.md`) |
| Spike-in | the paper spikes six purified *Haloferax volcanii* small proteins (HVO_0758, HVO_2212, HVO_2753, HVO_2922, HVO_2983, HVO_A0101), not C. difficile proteins. They are added to both search FASTAs as `SPIKE_` records (translations of the RefSeq DS2 CDSs, GCF_000025685.1, the PXD066008 prep genome) so their spectra compete in the search; they are excluded from every C. difficile claim |
| Genome | GCF_002080065.1 (ASM208006v1): chromosome NZ_CP016318.1 (4,293,048 bp) and plasmid NZ_CP016319.1 (7,881 bp), both circular-tagged; `cdif_630derm.fna.gz` sha256 `49065d7c…` |
| Annotation | Bakta 1.12.2, DB 6.0 full (job 2001987_11): 3,820 CDS, 6 sORF |
| Answer-key annotation (never an input) | RefSeq `GCF_002080065.1_ASM208006v1_genomic.gff.gz` (3,801 CDS, 49 pseudogenes) |
| GIPPA code | GIPPA 1.1.1 release, detached worktree `local/worktrees/rel111-pin` at `bfc74ff6`, `PYTHONPATH=<pin>/src`, venv Python (never `uv run`) |
| Search | DIA-NN 2.3.1 Academia, pinned image (`validation/pilot/diann_image.sha256`), the two calls and flags of `validation/pilot/diann_A2.cfg` (`validation/campaign/p2b/diann_campaign.cfg`, identical flags; `--report-decoys` kept; no `--met-excision`, `--reanalyse`, `--rt-profiling`); Thermo `.raw` read directly (ADR-0077 Amendment 1); FASTA = the `gippa makedb` database, which appends the contaminant library as `CONT_` records |
| Run | `gippa run` all modules, default seed 20261014, `--threads 14 --memory 54.GB`, profile `farmshare_local,apptainer` inside one 16-CPU Slurm allocation |

## Loci (Bakta 1.12.2 locus tags)

| Gene | Locus | Coordinates (strand) | Note |
|---|---|---|---|
| rplJ / rplL | FLEEGO_00125 / FLEEGO_00126 | 92,876-93,382 / 93,439-93,804 (+) | S3 |
| atpA / atpD | FLEEGO_03723 / FLEEGO_03721 | 4,064,016-4,065,518 / 4,061,703-4,063,097 (−) | S3: the F0F1 copies; V-type `atpA`/`atpD` (FLEEGO_03149/03147) make the name-based `te4` ambiguous |
| grdB | FLEEGO_02514 | 2,677,566-2,678,876 (−); Sec TGA 2,677,827-2,677,829 (U350) | Bakta translates the full protein (U) |
| grdA | FLEEGO_02516 (N-terminal 44 aa) + FLEEGO_02515 (107 aa) | 2,679,247-2,679,381 and 2,678,908-2,679,231 (−); Sec TGA 2,679,247-2,679,249; in-frame gap 2,679,232-2,679,246 (5 codons) | Bakta splits grdA at the Sec codon and restarts at the next ATG |
| prdB | FLEEGO_03446 | 3,750,189-3,750,641 (−), ends at the Sec TGA 3,750,189-3,750,191 | the 90-aa C-terminal part after the Sec (3,749,919-3,750,188, −; stop 3,749,916-3,749,918) is absent from route A |

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

## Dataset checks (handover plan, Phase 2 item 5)

| # | Check | Expected | Verdict rule | Read from |
|---|---|---|---|---|
| C1 | Small proteins | detected share of annotated small proteins (CDS + sORF ≤ 100 aa, spike-ins excluded) in A and in B. *Prior*: below the paper's 61.4 % of predicted small proteins, which pooled every workflow and enrichment | `reported` (Wilson CI) | A/B `genes.tsv` |
| C2 | Novel ORFs vs the paper | the paper searched only the curated annotation (Dannheim et al. 2017, 3,781 entries) and reports no novel small protein, so a called ORF that matches no RefSeq CDS/pseudogene is classed `new`; called ORFs ≤ 100 aa listed separately | `reported`; not testable if 0 called | B `novel_orfs.tsv` |
| C3 | Spike-in proteins | every accepted spike-in (`SPIKE_`) peptide stays off the C. difficile genome (in `unmapped.tsv` or absent from `peptide_hits.tsv`); spike-in peptides counted per protein. *Prior*: at 0.5 pg/µg few or none are identified | pass if 0 spike-in peptides placed on the genome; not testable if none accepted | A/B peptides, `peptide_hits.tsv`, `unmapped.tsv` |
| C4 | Selenoprotein GrdB | if detected, peptides on both sides of the Sec TGA; none across it (makedb writes U) | pass/fail; not testable if not detected | A `peptide_hits.tsv` |
| C5 | Selenoprotein GrdA | if either fragment is detected, report peptides per fragment and in the 5-codon gap; peptides in the gap in route B land in an `nterm_extension` candidate of FLEEGO_02515 | `reported` | A/B `peptide_hits.tsv`, B `novel_orfs.tsv` |
| C6 | Selenoprotein PrdB C-terminus | if route B places peptides in 3,749,919-3,750,188 (−), they lie in a six-frame candidate (any class); route A has none there by construction | pass if every route-B peptide in the window belongs to a candidate; not testable if none | B `peptide_hits.tsv`, `novel_orfs.tsv` |

S13 (determinism, TE7) is run on PXD072227 in this group, not here.
