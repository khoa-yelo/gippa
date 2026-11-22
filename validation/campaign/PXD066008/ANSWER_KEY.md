# PXD066008 answer key (pre-registered before any search)

Campaign Phase 2, dataset 9: *Haloferax volcanii* H26 (wild type of the zinc-finger µ-protein
study, PMC12949519), 5 dia-PASEF runs (2 timsTOF Pro 2, 3 timsTOF HT; `run_selection.tsv`), on
the DS2 genome GCF_000025685.1 (5 replicons, all tagged circular). Committed before any
spectrum is searched; a result found only after looking is reported as exploratory. Aggregates
only.

## Fixed inputs

| Item | Value |
|---|---|
| GIPPA code | detached worktree `local/worktrees/rel111-pin` at `bfc74ff6` (release 1.1.1), `PYTHONPATH=<pin>/src` |
| Genome | `hvol_ds2.fna.gz` (GCF_000025685.1, all 5 replicons `[topology=circular]`), sha256 `94be4069…` (`PREP.md`) |
| Annotation | Bakta 1.12.2, DB 6.0 full, job 2001987_15; JSON sha256 `3fe1e7fb6098…` |
| Answer-key annotation (never an input) | RefSeq GFF `GCF_000025685.1_ASM2568v1_genomic.gff.gz`, MD5 `476dedf9…` |
| Search | DIA-NN 2.3.1 (pinned image), flags of `validation/pilot/diann_A2.cfg` (`p2c/diann_campaign.cfg`); contaminants as makedb `CONT_` records |
| Routes | A: `makedb` called genes; B: `+ --six-frame-orfs --min-orf-aa 20`; both `--search-db same`; `gippa compare` A vs B |

## Archaeal annotation check (handover plan; run before any search)

Bakta does not officially support archaea. Rule (M. hyorhinis checks): Bakta CDS + sORF count
within 10 % of the RefSeq non-pseudo CDS count, and ≥ 90 % of RefSeq CDS share their stop codon
(contig, strand, 3' end) with a Bakta CDS. Result (`annotation_check.json`, script
`p2c/annotation_check.py`), computed and recorded before this key was committed:

| | value | rule | result |
|---|---|---|---|
| Bakta CDS + sORF / RefSeq CDS | 4,023 / 3,897 = 1.032 | 0.90-1.10 | pass |
| RefSeq CDS with a Bakta stop | 3,840 / 3,897 = 98.5 % [98.1, 98.9] | ≥ 90 % | pass |

**Pass: the run proceeds with Bakta.** Bakta names almost nothing here (3,488 of 4,023 CDS
hypothetical; 1 of 55 RefSeq ribosomal proteins carries a ribosomal product name), so checks that
need gene names use RefSeq names carried to the Bakta CDS with the same stop codon (rule fixed
here: same contig, strand and 3' end; 55 / 55 RefSeq ribosomal proteins map).

## Loci (1-based, RefSeq GFF; Bakta locus by shared stop)

| Locus | RefSeq | Bakta | Use |
|---|---|---|---|
| pyrE2 (HVO_0333, HVO_RS06305) | NC_013967.1:301,751-302,281 (+) | FADMKH_00337 (same coordinates, "Orotate phosphoribosyltransferase") | H26 is ΔpyrE2 (paper Methods): not detected |
| pHV2 | NC_013965.1, 6,359 bp | 8 CDS | H26 is cured of pHV2: no gene detected |

## Acceptance checks

Verdicts: pass, fail, not testable (or reported, where no threshold applies). Shares carry a
Wilson 95 % CI; other numbers their GIPPA error estimate or `unestimated`. Thresholds from
GIPPA_SPEC.md §11.6 where they exist; none is loosened.

| # | Question | Expected (threshold) | Read from |
|---|---|---|---|
| 1 | Mapping (TE1) | ≥ 99.5 % of route-A accepted peptides map exactly in frame (pilot `te1_real`) | `mapping/`, `mapping_fdr.tsv` |
| 2 | Expression (TE3) | ribosomal proteins ≥ 90 % detected; median abundance in the top quartile. GIPPA's own ribosomal set (Bakta names) has < 20 members here, so the pre-registered set is the 55 RefSeq ribosomal proteins mapped by stop; GIPPA's set reported alongside | `genes.tsv` |
| 3 | Stoichiometry (TE4) | **not testable** by construction: no bacterial L7/L12-L10 or F1 α/β pair in archaea (A-type ATP synthase, P1/P0 stalk), and Bakta names neither | - |
| 4 | Contaminants (TE5) | every accepted trypsin/keratin peptide in the contaminant bucket, 0 placed on the genome | `unmapped.tsv` |
| 5 | Gene-level FDR | decoy-estimated gene FDR of route-A detected genes ≤ 0.01 (tier `decoy_picked`) | `genes.tsv`, `detection_summary.tsv` |
| 6 | Operons | adjacent pairs within RefSeq-named ribosomal-protein and ATP-synthase (A-type, product "ATP synthase subunit …" / "V-type ATP synthase …") units: no discordant flag and no within-operon q ≤ 0.05 | `pairs.tsv`, `operon_flags.tsv` |
| 7 | Genome context | geNomad provirus regions, if any: observed < expected detection and q ≤ 0.05 in every region (Phase 1 rule); not testable if none. Islands and per-replicon detection reported | `regions_stats.tsv`, `islands.tsv` |
| 8 | Negative control: pHV2 | 0 genes on NC_013965.1 detected at gene q ≤ 0.01 (routes A and B) | `genes.tsv` |
| 9 | Negative control: pyrE2 | FADMKH_00337 not detected (routes A and B) | `genes.tsv` |
| 10 | Novel ORFs vs RefSeq | every called route-B ORF compared with RefSeq CDS/pseudogenes (≥ 50 % same-strand overlap = match); share with Wilson CI (reported). Each call classed known-RefSeq, known-paper (µ-protein list) or new | `novel_orfs.tsv` vs RefSeq |
| 11 | µ-proteins Bakta missed | of the 17 RefSeq CDS ≤ 70 aa without a Bakta stop match, the number with route-B peptide evidence (a candidate, or six-frame ORF hits) (reported, exploratory) | `novel_orfs.tsv`, `peptide_hits` |
| 12 | Class q floor | a class with < ~100 candidates reports `novel_orfs.class_q_floor_above_alpha`; alpha 0.01 not loosened (pass when GIPPA warns for every such class) | run warnings, `novel_orfs.tsv` |
| 13 | Route B cost | route B keeps ≥ 95 % of route A's detected genes (Phase 1: 98.4 %); DB, library and wall times reported | `gippa compare`, provenance |
| 14 | Variants | **not testable**: same-genome database search, no de novo arm; any `variants.tsv` row is reported | `variants.tsv` |
| 15 | Second round | only if route B calls ORFs: `makedb --patch`, search again, patched proteins detected at gene q ≤ 0.01 | round-2 `genes.tsv` |

Expression context, not a check: the paper reports 3,072 proteins over all 20 strains (DIA-NN 1.9,
library-free).
