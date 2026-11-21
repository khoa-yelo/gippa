# Phase 1 answer key (pre-registered before any search)

Pre-registration for the E. coli K-12 MG1655 demonstration of the demo and discovery campaign
(handover plan of 2026-11-17). Committed before any spectrum is searched. A result found only
after looking is reported as exploratory.

## Fixed inputs

| Item | Value |
|---|---|
| GIPPA code | detached worktree at `a2a95c76` (tag v1.1.0; code-identical to release `e43683c7`) |
| Genome | GCF_000005845.2 (NC_000913.3, 4,641,652 bp, circular). `genomic.fna.gz` MD5 `c13d459b5caa702ff7e1f26fe44b8ad7`, matching NCBI `md5checksums.txt`. Its sequence is identical to U00096.3 (`tests/data/genomes/k12_mg1655.fna.gz`): the MD5 of the uppercased sequence is `482a2b04485ec8c4b5f4eaba2c2002da` for both |
| Annotation | the existing Bakta 1.12.2 full DB 6.0 recording on U00096.3 (image `sha256:1e24ae2f…`, job 1849096, `$GIPPA_ROOT/recordings/bakta_full_circ/k12_mg1655/`). Reused, not re-run: same tool, image, database and sequence |
| Answer-key annotation (never an input) | RefSeq GFF `GCF_000005845.2_ASM584v2_genomic.gff.gz`, MD5 `0f52ffc94af5ddf544ff89cc6f546b0c` (matches NCBI) |
| Data | PXD058422 and PXD066241 (Mol Cell Proteomics 2025, PMC12529367); runs chosen in the Phase 0 inventory (ADR-0077 pending) |

Coordinates below are 1-based inclusive, as in the GFF.

## Loci the checks use (from the RefSeq GFF and the Bakta JSON)

| Locus | RefSeq | Bakta 1.12.2 call | Consequence for the check |
|---|---|---|---|
| crl (b0240) | pseudogene 257,829–259,006 (+), split by IS1I (mobile element 257,908–258,675, −) into fragment 1 (257,829–257,899) and fragment 2 (258,676–259,006) | CDS 258,734–259,006 (+), named crl: fragment 2 from an internal start only | The strain used is crl+ (no IS1). Expected: C-terminal peptides map to the Bakta crl CDS; peptides in 258,676–258,733 are an `nterm_extension` candidate of crl; peptides from fragment 1 lie on the + strand opposite the IS1 genes (an `antisense` or `intergenic` candidate, depending on the six-frame ORF); the peptide spanning the fragment 1/2 junction does not map (unmapped) |
| prfB (b2891) | CDS in two parts with ribosomal slippage: 3,036,208–3,036,282 and 3,035,184–3,036,206 (−) | CDS 3,035,184–3,036,065 (−): the +1 frame part from an internal start | GIPPA has no frameshift class (spec §16 puts prfB +1 recoding out of scope). Expected: peptides in 3,036,066–3,036,206 are an `nterm_extension` candidate; peptides from the 25-codon leader (other frame) are a `cds_outframe` candidate or no call. Reported as where the peptides land, not as a frameshift call |
| fdnG (b1474) | CDS 1,547,401–1,550,448 (+), Sec at 1,547,986–1,547,988 | full-length CDS with a SECIS region | if detected: peptides on both sides of the in-frame TGA; any peptide across it is a read-through, not an ALT call |
| fdoG (b3894) | CDS 4,082,772–4,085,822 (−), Sec at 4,085,235–4,085,237 | full-length CDS (Bakta names it fdnG) | as fdnG |
| fdhF (b4079) | CDS 4,297,219–4,299,366 (−), Sec at 4,298,947–4,298,949 | full-length CDS with a SECIS region | as fdnG |

## Acceptance checks

Verdicts are pass, fail or not testable, each with its number and error estimate (or
`unestimated`). Thresholds come from GIPPA_SPEC.md §11.6 and are not loosened.

| # | Question | Expected | Read from |
|---|---|---|---|
| 1 | Mapping (TE1) | ≥ 99.5 % of E1 peptides map exactly in frame | `mapping/`, `mapping_fdr.tsv` |
| 2 | Expression (TE3) | ribosomal proteins ≥ 90 % detected; median abundance in the top quartile | `genes.tsv`, `detection_summary.tsv` |
| 3 | Stoichiometry (TE4) | log2(RplL/RplJ) in [1, 3]; log2(AtpA/AtpD) in [−1, 1]. The RplL/RplJ pair failed on the Nissle pilot (waived for 1.1.0, ADR-0076), so a second failure is a known risk | `genes.tsv` |
| 4 | Contaminants (TE5) | trypsin and keratin peptides land in the contaminant bucket | `unmapped.tsv`, V05 |
| 5 | Genome enrichment | cryptic prophages (DLP12, e14, Rac, Qin, CP4-6, CP4-57 and others) detected less than detectability predicts | `regions_stats.tsv`, `islands.tsv`, V01, V07 |
| 6 | Operons | ribosomal and ATP synthase operon partners concordant | `operons.tsv`, `operon_flags.tsv`, V11 |
| 7 | Crl | as in the loci table: C-terminal peptides on the Bakta crl CDS, an N-terminal-extension and/or fragment-1 candidate, junction peptide unmapped | `genes.tsv`, `novel_orfs.tsv`, `unmapped.tsv`, V13 |
| 8 | PrfB | as in the loci table: where peptides upstream of the Bakta start land (`nterm_extension`, `cds_outframe` or no call) | `novel_orfs.tsv`, `unmapped.tsv` |
| 9 | Novel ORFs vs known genes | every called ORF compared with the RefSeq GFF; share matching an annotated gene or pseudogene, with a Wilson 95 % CI | `novel_orfs.gff3` vs RefSeq |
| 10 | Selenoproteins | if FdoG, FdnG or FdhF is detected, peptides on both sides of the in-frame TGA | `genes.tsv`, stop-codon flags |
| 11 | Negative controls (E3) | MntS, Lon and ClpP not detected in their own knockout runs; resistance-cassette peptides unmapped | E3 `genes.tsv`, `unmapped.tsv` |
| 12 | Variants, same strain (E4, TE2) | no `high` ALT call outside `start_codon` and `stop_codon` sites | `variants.tsv`, V12 |
| 13 | Modification-like calls (E4) | ALT calls flagged `mass_equivalent` listed separately, not counted as strain variants | `variants.tsv` |
| 14 | Route B cost | share of E1's detected genes that E2 keeps; six-frame database size; DIA-NN and GIPPA wall time | `gippa compare`, `provenance.json` |
| 15 | Class q floor | a novel-ORF class with too few candidates reports `novel_orfs.class_q_floor_above_alpha`; alpha is not loosened | `novel_orfs.tsv`, run log |
| 16 | Determinism (TE7) | E1 rerun with the same seed gives byte-identical `tables/` | `cmp` of the two `tables/` trees |
