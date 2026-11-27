# Phase 1b answer key: NCM3722 strain-variant arm (pre-registered before any search)

Pre-registration for Phase 1b of the demo and discovery campaign (handover plan of 2026-11-17,
"Phase 1b: NCM3722 strain-variant arm"). Committed on 2026-11-19 before any genome download,
spectrum search or GIPPA run of this arm, and before the truth set is built. A result found only
after looking is reported as exploratory. Thresholds come from the plan and GIPPA_SPEC.md §11.6
(TE8) and are not loosened; where the plan sets no threshold, the number is reported with its
error estimate and the verdict is `reported`. Expectations marked *prior* are the agent's own,
written down here so they cannot move later; they are not plan gates.

## Fixed inputs

| Item | Value |
|---|---|
| GIPPA code | detached worktree `local/worktrees/campaign-pin` at `a2a95c76` (v1.1.0), `PYTHONPATH=<pin>/src`, venv Python; never `uv run` |
| N1 genome | GCF_001043215.1 (ASM104321v1; NCM3722, isolate K-12): chromosome CP011495.1 (4,678,046 bp) and plasmid F CP011496.1 (67,545 bp), 4,745,591 bp. NCBI `md5checksums.txt` lists `genomic.fna.gz` MD5 `b3a7430b1beda7a8541a6953667cc28e`; NCBI nuccore esummary gives `topology=circular`, `completeness=complete` for both (NZ_ and GenBank accessions), so both headers get ` [topology=circular]` |
| N1 annotation | Bakta 1.12.2, full DB 6.0, `--keep-contig-headers`, the pinned image of `workflow/conf/containers.config` (the command of `gippa.annotate.bakta_cmd.bakta_command`); JSON and GBFF kept together |
| N2 genome | MG1655 U00096.3 (`tests/data/genomes/k12_mg1655.fna.gz`, circular-tagged), sequence-identical to GCF_000005845.2 / NC_000913.3 (uppercased-sequence MD5 `482a2b04485ec8c4b5f4eaba2c2002da` for both, Phase 1 answer key) |
| N2 annotation | the existing Bakta 1.12.2 full DB 6.0 recording on U00096.3 (job 1849096, `$GIPPA_ROOT/recordings/bakta_full_circ/k12_mg1655/bakta/k12_mg1655.{json,gbff}`), reused |
| Spectra | PXD060434 (PNAS 2026, PMC13684143), Orbitrap Exploris 480 DIA, batch B: `BBM_796_P110_C15_01_FRH_{049,056,061}.raw`, NCM3722 wild type, MOPS + 20 mM glucose, steady-state exponential (`run_selection.tsv`). No deposit checksum: byte size checked against the PRIDE file list, sha256 recorded |
| Search | DIA-NN 2.3.1 Academia (pinned image, `validation/pilot/diann_image.sha256`), the two calls and flags of `validation/pilot/diann_A2.cfg` (`--report-decoys` kept; no `--met-excision`, `--reanalyse`, `--rt-profiling`); database = `gippa makedb` route A on the NCM3722 Bakta JSON, which adds the contaminant library as `CONT_` records |

## Arms

| Arm | Peptide table | Genome and annotation | Relation |
|---|---|---|---|
| N1 | `gippa convert --from diann ... --search-db same` of the DIA-NN report | NCM3722 | default (`related`), as the standard recipe |
| N2 | the **same** DIA-NN report converted again with `--search-db other` | MG1655 | `--genome-relation related` |

Both: `gippa preflight --strict` exits 0 first; `gippa run` with all modules, `--threads 14
--memory 54.GB`, default seed.

## Truth set (method fixed here; built before N2 outputs are opened)

1. Proteins: the CDS and sORF translations of the two Bakta 1.12.2 full-DB annotations (same
   tool, version and database on both genomes, so gene calls are comparable).
2. Orthologs: BLAST+ `blastp` from the pinned GIPPA BLAST image
   (`quay.io/biocontainers/blast@sha256:9dfc69f9…`), NCM3722 vs MG1655 and MG1655 vs NCM3722,
   `-seg no -evalue 1e-10 -max_hsps 1`. A pair is kept when it is the reciprocal best hit by
   bitscore, with identity ≥ 90 % and the HSP covering ≥ 80 % of both proteins.
3. Truth sites: every aligned HSP column of a kept pair where the residues differ and neither is a
   gap. A site is (MG1655 locus tag, MG1655 residue position, REF = MG1655 residue,
   ALT = NCM3722 residue), mapped to the MG1655 codon (contig, codon start, strand) from the Bakta
   CDS coordinates, and to the NCM3722 codon the same way.
4. Not scored as substitution sites, counted only: I↔L differences (GIPPA folds I/L, so they are
   undetectable by design); gap columns (indels); residues outside the HSP (start or end
   differences); NCM3722 proteins without a kept pair (plasmid F genes, NCM3722-specific genes);
   pairs below the identity or coverage limits.
5. What is committed: the method, the BLAST version and image digest, and counts (pairs, sites,
   sites by replicon and class). No sequences.

## Scoring rules

- **Covered truth site**: a truth site whose NCM3722 codon lies inside the placement of at least
  one N1 target peptide accepted as evidence (N1 `mapping/peptide_hits.tsv`, exact route, in the
  ortholog's frame). All scores below are restricted to covered sites.
- **Eligible covered site** (secondary stratum): covered by at least one such peptide of length
  ≥ 10 whose number of truth differences on MG1655 is within the tolerant mapper's limit (1 at
  10–11 aa, 2 at 12–40 aa; `mapping.tolerant.max_diff`). Sites outside it cannot be called by
  design; reporting both separates design limits from misses.
- **N2 call matches a truth site**: an N2 `tables/variants.tsv` row (call `ALT` or `mixed`) at the
  same MG1655 contig, strand and codon whose `alt_aa` equals the truth ALT (I/L folded, `J`
  matching I or L). A two-codon `swap` or `block` row matches when its ALT residues equal the
  NCM3722 residues at the codons it spans and at least one of those codons is a truth site.
  Frameshift rows match nothing (no truth frameshift site in this method).
- Uncertainty: Wilson 95 % intervals for every share (peptides and sites are treated as
  independent, which is optimistic for sites in one gene).

## Checks

| # | Check | Expected | Verdict rule | Read from |
|---|---|---|---|---|
| 1b-1 | Recall | share of covered truth sites that N2 calls ALT or mixed (any confidence) with the matching ALT; also `high`-only, and both within the eligible stratum. *Prior*: eligible-stratum recall (any confidence) ≥ 0.8 | `reported` (no plan threshold); the prior is stated for honesty only | N2 `variants.tsv`, truth set, N1 `peptide_hits.tsv` |
| 1b-2 | Precision | share of N2 `high` rows with call `ALT` that match a truth site (`mixed` rows reported separately). *Prior*: ≥ 0.9 | `reported` | N2 `variants.tsv` |
| 1b-3 | Non-truth `high` ALT calls | each listed with contig, codon, strand, gene, REF/ALT, `n_alt`, flags, and classed as `mass_equivalent`, start/stop codon (`start_codon`, `stop_codon`, `putative_start`), or unexplained; candidates for Khoa, not confirmed variants | listed | N2 `variants.tsv` |
| 1b-4 | Unplaced share by origin | share of N1 evidence peptides (distinct, I/L-folded, target, q ≤ 0.01, not contaminant) with no accepted placement on MG1655 (absent from N2 `peptide_hits.tsv`), split by their N1 replicon: chromosome CP011495.1, plasmid CP011496.1, both. *Prior*: plasmid share ≥ 0.5 (F is absent from MG1655 apart from shared IS/Tn elements); chromosome share ≤ 0.05 | `reported` | N1 and N2 `peptide_hits.tsv`, N2 `unmapped.tsv` |
| 1b-5 | Placement FDR (TE8 analogue) | decoy-estimated placement FDR of N2's accepted placements at defaults ≤ 1 % | **pass/fail** on the pooled estimator of TE8 (`validation/proxies/te_proxies.py::pooled_decoy_fdr`, pinned code) applied to N2's peptides on MG1655; the run's own `mapping/mapping_fdr.tsv` strata are reported next to it | TE8-style chain, N2 `mapping_fdr.tsv` |
| 1b-6 | N1 vs N2 comparison | `gippa compare N1 N2`: Spearman ρ of abundance with its CI; Jaccard of islands and operon flags. `gippa compare` requires one genome and one locus-tag set, so it is **expected to exit 1** here (documented check); the result is recorded. Then the same quantities are computed on ortholog pairs from step 2 of the truth set: Spearman ρ of N1 vs N2 gene abundance on genes detected with measured abundance in both, 95 % CI from 1,000 gene bootstrap draws (seed 20261014); Jaccard of island membership and of flagged operon pairs mapped through orthologs. *Prior*: ρ ≥ 0.9 | `reported` | N1/N2 `genes.tsv`, `islands.tsv`, `operon_flags.tsv` |

Also reported, not gated: N1 detected genes by tier and gene-level FDR, N1 TE1-style exact in-frame
share, and wall times per step. Aggregates only in git: no peptide sequences, no DIA-NN rows.
