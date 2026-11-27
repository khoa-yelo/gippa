# PXD079213 answer key (pre-registered before any search)

Campaign Phase 2, dataset 10: *Saccharolobus islandicus* E233S control strain (PMC13736083,
iProX), 3 Orbitrap Astral DIA runs `cont_1..3` (8 min gradients; `run_selection.tsv`), on the
complete REY15A genome GCF_000189555.1 (1 circular chromosome; not the 141-contig assembly the
paper cites). Committed before any spectrum is searched; a result found only after looking is
reported as exploratory. Aggregates only.

## Fixed inputs

| Item | Value |
|---|---|
| GIPPA code | detached worktree `local/worktrees/rel111-pin` at `bfc74ff6` (release 1.1.1), `PYTHONPATH=<pin>/src` |
| Genome | `sisl_rey15a.fna.gz` (NC_017276.1 `[topology=circular]`), sha256 `5f5698b6…` (`PREP.md`) |
| Annotation | Bakta 1.12.2, DB 6.0 full, job 2001987_16; JSON sha256 `180c6a955535…` |
| Answer-key annotation (never an input) | RefSeq GFF `GCF_000189555.1_ASM18955v1_genomic.gff.gz`, MD5 `50d0f202…` |
| Search | DIA-NN 2.3.1 (pinned image), flags of `validation/pilot/diann_A2.cfg` (`p2c/diann_campaign.cfg`); Thermo `.raw` read directly (`raw_check.md`); contaminants as makedb `CONT_` records |
| Routes | A: called genes; B: `+ --six-frame-orfs --min-orf-aa 20`; both `--search-db same`; `gippa compare` A vs B |

## Archaeal annotation check (handover plan; run before any search)

Rule: Bakta CDS + sORF count within 10 % of RefSeq non-pseudo CDS, and ≥ 90 % of RefSeq CDS share
their stop codon with a Bakta CDS. Result (`annotation_check.json`, `p2c/annotation_check.py`):

| | value | rule | result |
|---|---|---|---|
| Bakta CDS + sORF / RefSeq CDS | 2,808 / 2,578 = 1.089 | 0.90-1.10 | pass |
| RefSeq CDS with a Bakta stop | 2,536 / 2,578 = 98.4 % [97.8, 98.8] | ≥ 90 % | pass |

**Pass: the run proceeds with Bakta.** Bakta names almost nothing (2,664 of 2,808 CDS
hypothetical), so name-based checks use RefSeq names carried to the Bakta CDS with the same stop
(same contig, strand, 3' end; 62 / 64 RefSeq ribosomal proteins map). The RefSeq annotation has
191 pseudogenes, many of them IS elements; novel ORFs matching them are classed known-RefSeq.

## Loci (1-based, RefSeq GFF; Bakta locus by shared stop)

| Locus | RefSeq | Bakta | Use |
|---|---|---|---|
| pyrF (orotidine-5'-phosphate decarboxylase) | SIRE_RS07070, 1,307,768-1,308,436 (−) | CGDCIA_01487 | E233S lacks pyrEF: not detected |
| pyrE | SIRE_RS07075, 1,308,423-1,309,010 (−) | CGDCIA_01488 | not detected |
| lacS (bgaS, β-galactosidase) | SIRE_RS11295, 2,067,746-2,069,215 (−) | CGDCIA_02367 | E233S lacks lacS: not detected |

Known risk, stated before looking: the control strain of an overexpression study may carry the
empty expression vector, and the usual S. islandicus vectors carry *S. solfataricus* pyrEF and
lacS as markers. Peptides identical between the orthologues would then make these genes
"detected" without contradicting the deletions. If detected, the report says so, gives the
peptide counts, and does not call it a GIPPA error; the verdict stays **fail** against this key.

## Acceptance checks

Verdicts: pass, fail, not testable (or reported where no threshold applies). Shares carry a Wilson
95 % CI; other numbers their GIPPA error estimate or `unestimated`. No threshold is loosened.

| # | Question | Expected (threshold) | Read from |
|---|---|---|---|
| 1 | Mapping (TE1) | ≥ 99.5 % of route-A accepted peptides map exactly in frame | `mapping/`, `mapping_fdr.tsv` |
| 2 | Expression (TE3) | ribosomal proteins ≥ 90 % detected; median abundance in the top quartile; set = RefSeq ribosomal proteins mapped by stop (GIPPA's Bakta-named set reported alongside) | `genes.tsv` |
| 3 | Stoichiometry (TE4) | **not testable** by construction (archaeal stalk and A-type ATP synthase; no names) | - |
| 4 | Contaminants (TE5) | every accepted trypsin/keratin peptide in the contaminant bucket, 0 placed on the genome | `unmapped.tsv` |
| 5 | Gene-level FDR | decoy-estimated gene FDR of route-A detected genes ≤ 0.01 | `genes.tsv`, `detection_summary.tsv` |
| 6 | Operons | adjacent pairs within RefSeq-named ribosomal-protein and ATP-synthase units: no discordant flag, no within-operon q ≤ 0.05 | `pairs.tsv`, `operon_flags.tsv` |
| 7 | Genome context | geNomad provirus regions, if any: observed < expected and q ≤ 0.05 in every region; not testable if none; islands reported | `regions_stats.tsv`, `islands.tsv` |
| 8 | Negative controls | pyrE, pyrF and lacS (CGDCIA_01488, _01487, _02367) not detected in routes A and B (vector caveat above) | `genes.tsv` |
| 9 | Novel ORFs vs RefSeq | every called route-B ORF compared with RefSeq CDS/pseudogenes (≥ 50 % same-strand overlap); share with Wilson CI (reported); each classed known-RefSeq, known-paper or new | `novel_orfs.tsv` |
| 10 | Class q floor | classes with < ~100 candidates report `novel_orfs.class_q_floor_above_alpha`; alpha not loosened | run warnings |
| 11 | Route B cost | route B keeps ≥ 95 % of route A's detected genes; DB, library, wall times reported | `gippa compare`, provenance |
| 12 | Variants | **not testable** (same-genome search, no de novo arm); `variants.tsv` rows reported | `variants.tsv` |
| 13 | Second round | only if route B calls ORFs: patched proteins detected at gene q ≤ 0.01 | round-2 `genes.tsv` |

Depth context, not a check: 8-min Astral gradients, so fewer detected genes than the other
datasets are expected.
