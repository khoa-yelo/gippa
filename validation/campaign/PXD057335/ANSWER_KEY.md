# PXD057335 answer key (pre-registered before any search)

Campaign Phase 2, dataset 7: *Mycobacterium tuberculosis* H37Rv, DMSO vehicle control, 48 h,
37 °C (PMC12913811), 3 Orbitrap Exploris 480 DIA runs `20230614_Exploris_IB_EXP57_DMSO_01..03`
(`run_selection.tsv`), on GCF_000195955.2 (NC_000962.3, 4,411,532 bp). Committed before any
spectrum is searched; a result found only after looking is reported as exploratory. Aggregates
only.

## Fixed inputs

| Item | Value |
|---|---|
| GIPPA code | detached worktree `local/worktrees/rel111-pin` at `bfc74ff6` (release 1.1.1), `PYTHONPATH=<pin>/src` |
| Genome | `mtb_h37rv_circ.fna.gz`: GCF_000195955.2 with NC_000962.3 tagged `[topology=circular]` (deviation below) |
| Annotation | Bakta 1.12.2, DB 6.0 full, same image and `bakta_command` argv as prep job 2001987_13, rerun on the tagged genome (job 2005457, `p2c/genome_bakta.sbatch`) |
| Answer-key annotation (never an input) | RefSeq GFF `GCF_000195955.2_ASM19595v2_genomic.gff.gz`, MD5 `4d53755c…` |
| Search | DIA-NN 2.3.1 (pinned image), flags of `validation/pilot/diann_A2.cfg` (`p2c/diann_campaign.cfg`); Thermo `.raw` read directly; contaminants as makedb `CONT_` records |
| Routes | A: called genes; B: `+ --six-frame-orfs --min-orf-aa 20`; both `--search-db same`; `gippa compare` A vs B |
| TerminatorNet | sharded path (M2-15b) is the 1.1.1 default (`params.terminatornet_sharded = true`); nothing to set |

## Deviation from the recipe (lead decision, 2026-11-21)

NCBI nuccore marks NC_000962.3 linear, so the recipe would leave it untagged; the chromosome is
circular. It is tagged `[topology=circular]` for this run. GIPPA reads replicon topology from the
Bakta JSON when the annotation is precomputed (`gippa.annotate.bakta`, `topology` of the sequence
record), and the prep JSON says linear, so tagging the FASTA alone would leave GIPPA's run linear
while preflight saw circular. Bakta was therefore rerun on the tagged genome with the same image
and argv (Bakta replicon table topology `c`). The prep annotation is kept on disk, unused.

## Loci (1-based, RefSeq GFF; Bakta locus by shared stop)

| Locus | RefSeq | Use |
|---|---|---|
| esxA (Rv3875, ESAT-6) | 4,352,609-4,352,896 (+) | detected (row 7) |
| esxB (Rv3874, CFP-10) | 4,352,274-4,352,576 (+) | detected (row 7) |
| rplJ / rplL (Rv0651 / Rv0652) | 748,276-748,812 / 748,849-749,241 (+) | TE4 |
| atpA / atpD (Rv1308 / Rv1310) | 1,463,228-1,464,877 / 1,465,841-1,467,301 (+) | TE4 |
| PE/PPE family | 155 RefSeq CDS with product matching `PE`, `PPE`, `PE-PGRS`, `PE family` or `PPE family` | row 8 |

## Acceptance checks

Verdicts: pass, fail, not testable (or reported where no threshold applies). Shares carry a Wilson
95 % CI; other numbers their GIPPA error estimate or `unestimated`. No threshold is loosened.

| # | Question | Expected (threshold) | Read from |
|---|---|---|---|
| 1 | Mapping (TE1) | ≥ 99.5 % of route-A accepted peptides map exactly in frame | `mapping/`, `mapping_fdr.tsv` |
| 2 | Expression (TE3) | ribosomal proteins (GIPPA's set, Bakta names) ≥ 90 % detected; median abundance in the top quartile | `genes.tsv`, `detection_summary.tsv` |
| 3 | Stoichiometry (TE4) | log2(RplL/RplJ) in [1, 3]; log2(AtpA/AtpD) in [−1, 1]. RplL/RplJ failed in the Nissle pilot and in Phase 1 E1 (L7 N-acetylation not searched): a known risk | `genes.tsv` |
| 4 | Contaminants (TE5) | every accepted trypsin/keratin peptide in the contaminant bucket, 0 placed on the genome | `unmapped.tsv` |
| 5 | Gene-level FDR | decoy-estimated gene FDR of route-A detected genes ≤ 0.01 | `genes.tsv`, `detection_summary.tsv` |
| 6 | Operons | ribosomal-protein and ATP-synthase operon partners concordant: no discordant flag, no within-operon q ≤ 0.05 | `pairs.tsv`, `operon_flags.tsv` |
| 7 | Known biology: ESX-1 | EsxA and EsxB both detected at gene q ≤ 0.01 (route A) | `genes.tsv` |
| 8 | Known biology: PE/PPE | the detected share of PE/PPE CDS is lower than that of all other CDS (two-sided Fisher exact p ≤ 0.05); repetitive, GC-rich PE-PGRS proteins give few unique tryptic peptides | `genes.tsv` (RefSeq product carried by stop) |
| 9 | Genome enrichment | geNomad prophage regions (phiRv1 near 1.78 Mb, phiRv2 near 2.97 Mb, if called): observed < expected and q ≤ 0.05 in every region; not testable if none called. Islands reported | `regions_stats.tsv`, `islands.tsv` |
| 10 | Novel ORFs vs RefSeq | every called route-B ORF compared with RefSeq CDS/pseudogenes (≥ 50 % same-strand overlap); share with Wilson CI (reported); each classed known-RefSeq, known-paper or new | `novel_orfs.tsv` |
| 11 | Class q floor | classes with < ~100 candidates report `novel_orfs.class_q_floor_above_alpha`; alpha not loosened | run warnings |
| 12 | Route B cost | route B keeps ≥ 95 % of route A's detected genes; DB, library, wall times reported | `gippa compare`, provenance |
| 13 | Variants | **not testable** (same-genome search, no de novo arm); `variants.tsv` rows reported | `variants.tsv` |
| 14 | Second round | only if route B calls ORFs: patched proteins detected at gene q ≤ 0.01 | round-2 `genes.tsv` |
| 15 | Circular topology | the run's `replicons` table lists NC_000962.3 as circular (deviation took effect) | run manifest / `replicons` |
