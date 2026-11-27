# GIPPA: genome-integrated analysis of bacterial proteomes with an error estimate for every claim

<!-- Draft manuscript, round 1 (docs agent). Every number is copied from the file named in the
`src` comment after its paragraph or table; paths are relative to the repository root unless they
start with $GIPPA_ROOT. Abbreviations used in the comments:
CAL111 = $GIPPA_ROOT/logs/calibration/full-d2bd17e42dcf-20261125-073329 (calibration_full on
d2bd17e4, Slurm array 2022338, the run readiness 1.1.1 reuses);
CAL110 = $GIPPA_ROOT/logs/calibration/full-4a6fa998300c-20261109-160014 (the 1.1.0 full tier);
RD = validation/release/1.1.1/readiness.json and READINESS.md;
PIL = validation/pilot/metrics/metrics.json; PIL111 = validation/pilot/recheck_1_1_1/metrics_1_1_1.json;
CAMP = validation/campaign/ (not yet on main: Phase 1 in worktree campaign-ecoli, Phase 2 summary,
group A and novel_pc/ in campaign-p2a, group B in campaign-p2b, group C in campaign-p2c);
P2S = CAMP/PHASE2_SUMMARY.md; NPC = CAMP/NOVEL_ORF_PC_THRESHOLD.md and CAMP/novel_pc/.
Records of CAL111 equal those of CAL110 except the ADR-0078 decoy-loss attack streams (4 records,
lower FDP) and the new decoy-pairing attack (3 records); checked record by record, 1,103 shared keys.
Aggregates only: no peptide sequence, no DIA-NN row, no per-gene value except named marker genes. -->

**Authors.** [PLACEHOLDER: Khoa: author list, affiliations, ORCID iDs, corresponding author]

**Venue.** [PLACEHOLDER: Khoa, Q-K11. Draft written for a bioRxiv preprint, then a journal to be chosen]

## Abstract

**Motivation.** Bacterial proteomics is usually analysed against a reference proteome, as lists of
proteins. Questions about genome context (whether prophages, plasmids, operons or biosynthetic
clusters are expressed) need the strain's own genome and an error estimate for every claim, and
existing workflows rarely provide both.

**Results.** GIPPA takes one assembled genome and the output of any search engine, places every
peptide on the genome with a decoy-estimated placement FDR, and builds one gene table keyed by
Bakta locus tags with gene-level FDR, abundance and detectability. Around that table it tests
region enrichment, abundance islands, operon concordance, function, variants and novel ORFs. Each
statistical output carries a q-value, FDR, FWER or confidence interval, or is flagged unestimated
or uncalibrated. In 1,208 passing simulation records of the release calibration, the known
failures are listed, and two methods ship flagged. On an *E. coli* Nissle 1917 pilot, a strain
database detected 3,051 of 4,685 genes (gene q ≤ 0.01) against 2,645 with a K-12 reference. In
eleven public DIA datasets (17 genomes, archaea included), every accepted peptide mapped exactly,
and the six-frame route kept at least 96.8 % of detected genes. No novel ORF was called, because
the small candidate classes could not reach q ≤ 0.01. Waived and deferred release criteria (the
runtime target, one ribosomal stoichiometry check) are reported as such.

**Availability.** MIT licence; Python package `gippa` 1.1.1 on PyPI; source at
https://github.com/khoa-yelo/gippa. Third-party tools run as separate programs under their own
licences.

<!-- src: abstract word count ≤ 250 (checked with wc); 1,208: fig2_calibration.tsv panel b status
pass (CAL111); 3,051 / 2,645 of 4,685: PIL111 › m3.A2/A1.subsets.all; eleven datasets, 17
genomes, 96.8 %: P2S §2.1, §2.4 and CAMP/ecoli/REPORT.md row 1, row 14; 0 novel ORFs: P2S §2.6;
NPC §1; waivers: RD rows R2, R5, R7 (ADR-0080) -->

## 1. Introduction

A bacterial proteome is shaped by the genome that encodes it. Genes in one operon are transcribed
together; prophages and plasmids carry genes that may be silent in a given condition; codon usage,
promoters and terminators affect how much of a protein is made. Proteomics measures the outcome,
but the usual analysis loses the genome: spectra are searched against a reference proteome, often
of a type strain, and the result is a list of protein groups. Genes that the reference lacks
cannot be found, peptides that differ from the reference are lost or misassigned, and the
coordinates needed to ask about genomic context are gone [CITE: proteogenomics of bacteria,
reference-database mismatch].

Proteogenomic workflows restore the genome by searching six-frame translations or strain
databases [CITE: bacterial proteogenomics workflows and six-frame searches]. They raise two
statistical problems. First, a larger search space and a mapping step add sources of error that
the search engine's FDR does not cover: a peptide can be identified correctly and still be placed
at the wrong locus, and a gene can be called detected on peptides that would match a decoy as
well. Second, context questions are multiple-testing problems on dependent data. Neighbouring
genes share transcription units (TUs), abundance is autocorrelated along the replicon, and a test
that treats genes as independent will reject too often [CITE: spatial autocorrelation of gene
expression along bacterial chromosomes].

We built GIPPA (genome integrated prokaryotic proteomics analysis) to answer such questions for a
single sample, on any bacterial or archaeal genome, with an error estimate attached to every claim
or an explicit statement that none exists. Three design choices follow from that aim. GIPPA does
not search spectra: it takes peptides or proteins from any engine and inherits the engine's
identification error, which it reports but does not re-estimate. It estimates the two error
sources it adds itself: the placement FDR of a peptide on the genome and the FDR of calling a gene
detected. And it treats calibration as part of the software: every statistical method has a
simulation suite with pre-set tolerances and attack nulls built to make it fail, and a method that
fails is shipped flagged or with a narrowed claim, not silently.

This paper describes GIPPA 1.1.1 and the evidence behind it: simulation calibration, eight
true-example assertions on real genomes, a pre-registered pilot on *E. coli* Nissle 1917, and a
campaign on public data-independent acquisition (DIA) datasets from eleven organisms. We report
what passed, what failed, what was waived and what is deferred. The Supplementary Methods
(`docs/preprint/methods.md`) and Supplementary Results (`docs/preprint/results.md`) give every
formula, default and record.

<!-- src: docs/preprint/outline.md §2 item 2; methods.md §1 (GIPPA does not search spectra;
inherits identification error); GIPPA_SPEC.md §0.2 R3-R4; CITATION.cff (name) -->

## 2. Design and implementation

### 2.1 Inputs and routes

The required input is a genome FASTA file. A Bakta JSON annotation [Schwengers et al. 2021] is
used as is; a GFF3 or GenBank annotation from another tool is passed to Bakta, which keeps its CDS
and adds the rest. Evidence can be a peptide table, a protein FASTA file with an abundance table,
or both. `gippa convert` writes the peptide table from DIA-NN 2.x reports [Demichev et al. 2020],
Cascadia 0.0.7 de novo output [Sanders et al. 2025] or a generic table; only the peptide column is
required, and q-values, PEPs, decoy flags, intensities and run labels are used when present
(Supplementary Methods §2). Rows above the input q cut-off (0.01) are kept for FDR estimation but
never count as evidence.

For users who search their own spectra, `gippa makedb` writes a strain database from the Bakta
annotation: the CDS and sORF proteins, optionally six-frame open reading frames (ORFs) of at least
20 codons with class-tagged identifiers, optionally proteins patched by a previous run, and a
381-entry contaminant library [Frankenfield et al. 2022]. It adds no decoys; decoys come from the
search engine. We call a search against Bakta proteins route A and a search that adds six-frame
ORFs route B; the campaign (Section 4.4) ran both on the same spectra.

<!-- src: methods.md §2 (Peptides, Converters), §4; gippa makedb --help (cli_help: --six-frame-orfs,
--min-orf-aa default 20, --patch, --no-contaminants); src/gippa/data/db_sources.yaml (contaminants,
381 entries, methods.md §5.5) -->

### 2.2 Annotation layers

Bakta 1.12.2 (database 6.0) supplies genes, CRISPR arrays, resistance genes and the origin of
replication. geNOMAD 1.12.0 [Camargo et al. 2024] calls prophages and plasmids, antiSMASH 8.0.4
[Blin et al. 2025] calls biosynthetic gene clusters (BGCs), and UniOP [Su et al. 2024] gives operon
pair probabilities, with an intergenic-distance model as fallback. Transcription units are maximal
runs of same-strand genes linked with probability ≥ 0.5. Promotech [Chevez-Guardado and
Peña-Castillo 2021], Promoter Calculator [LaFleur et al. 2022], TerminatorNet [Tjaden 2026],
TransTermHP [Kingsford et al. 2007] and OSTIR [Roots et al. 2021] add promoter, terminator and
translation-initiation layers, and codon-usage indices are computed in process. These layers
predict σ70 promoters and intrinsic terminators only, so a missing call is never read as absence;
a missing layer is flagged rather than reported as "no". The operon and terminator cut-offs were
fixed on *E. coli* K-12 RegulonDB 14.5.0 data [Salgado et al. 2024] (Supplementary Methods §3).

<!-- src: methods.md §3 (tool versions, TU rule, cut-offs from RegulonDB 14.5.0); NOTICE.md -->

### 2.3 Mapping peptides to the genome

Every replicon is translated in six frames in the run's genetic code (11, 4 or 25). Peptides of
7–40 residues are matched exactly with an Aho–Corasick automaton [Aho and Corasick 1975], with
I and L folded and circular replicons extended across the origin. Peptides of 10–40 residues
without an exact hit are searched with up to one (10–11 residues) or two (12–40 residues)
differences, where a difference is a substitution, an adjacent swap, a mass-equivalent block
exchange or a ±1-nt frameshift. Proteins longer than 40 residues are aligned with miniprot
[Li 2023] and, failing that, tblastn [Camacho et al. 2009]. Each placement is classed by locus
(in-frame CDS, sORF, N- or C-terminal extension, out-of-frame, antisense, intergenic, pseudogene,
RNA overlap) and by uniqueness. Peptides that place nowhere are matched against the contaminant
library or written to `unmapped.tsv` with a reason (Supplementary Methods §5).

<!-- src: methods.md §2 (genetic codes 11, 4, 25), §5.1-5.5 -->

### 2.4 The gene table

The gene table has one row per annotated CDS, sORF or pseudogene, keyed by its Bakta locus tag. A
peptide counts as evidence for a gene when its input q and placement q are both ≤ 0.01 (or it is
an exact hit from a same-strain search), its locus class assigns it to that gene, and it is unique
to the gene or shared only with identical proteins. The detection tier is chosen once per run from
what the input carries: `decoy_picked` (decoys with protein identifiers that pair to genes),
`mayu` (targets with q or PEP but unusable decoys), `two_peptide` (no q or PEP, as for de novo
input) or `protein_list` (proteins only). The first two estimate a gene-level FDR (Section 3.2); the
last two do not, and say so. Abundance is iBAQ [Schwanhäusser et al. 2011] by default (Top3 as an
alternative) on unique peptides; GIPPA does not normalise across samples. A cross-fitted
detectability model predicts each gene's probability of detection from its observable peptides,
length, hydropathy, pI and transmembrane segments, so that context tests can compare observed with
expected detection (Supplementary Methods §6–8).

<!-- src: methods.md §6 (evidence rule, tiers), §7 (abundance), §8 (detectability) -->

### 2.5 Context, operon, variant, novel-ORF and function modules

Region tests ask whether the genes of a prophage, BGC, CRISPR or resistance region are detected
more or less often than detectability predicts, or are more or less abundant, than genes at
random positions on the same replicon. Autocorrelation statistics (TU-block Moran's I and join
counts along the replicon), an abundance-island scan with family-wise error control, and an
origin-to-terminus gradient describe abundance along each replicon. The operon module tests
discordance of abundance within predicted operons and looks for missed operons among non-operon
pairs. The variant module calls amino-acid differences from tolerant placements and flags
mass-equivalent calls that could be modifications. The novel-ORF module groups peptides that fall
outside annotated genes into candidate ORFs by class (N-terminal extension, intergenic, out of
frame, antisense) and estimates a class-wise FDR from engine decoys of the six-frame database.
Function enrichment runs over-representation analysis (ORA) for detection and fgsea
[Korotkevich et al. 2021] for abundance, over GO, KEGG and COG collections, with TUs as the unit.
A determinants module models detection and abundance from sequence features. `gippa compare`
reports agreement between runs on the same genome (Supplementary Methods §9–17).

<!-- src: methods.md §9-17 -->

### 2.6 Interface, outputs and reproducibility

GIPPA is a Python package (Python ≥ 3.11) with a command-line interface of nine commands (Table 1)
and a bundled Nextflow workflow [Di Tommaso et al. 2017] with profiles for local, Slurm, conda,
Apptainer and Docker execution. Each run writes 22 documented result tables (18 under `tables/`,
4 under `mapping/`), each a TSV file with a JSON schema; `gippa compare` adds three; nine
plot-support tables back the figures (Table 2). A self-contained HTML report holds 15 plots in a
colour-blind-safe palette, each with its data table. The run manifest records input checksums, tool
and database versions, container digests, the resolved configuration and the seed. All randomness
flows from one master seed (20261014 by default) through per-step generators, so the same inputs
and seed give byte-identical tables (TE7, Section 4.2).

<!-- src: gippa --help and gippa <cmd> --help (captured 2026-11-27 at the 1.1.1 tree of this branch);
docs/reference/tables.md (34 documented tables: 18 tables/, 4 mapping/, 3 compare, 9 plots/support);
methods.md §19-20; src/gippa/config/defaults.yaml:2 (seed) -->

**Table 1.** Commands and options of `gippa` 1.1.1 (from `gippa <command> --help`). Every command
also takes `--config FILE`, `--set KEY=VALUE` (any configuration key, repeatable), `--log-level`
and `--version`, except `demo` and `example`.

| Command | Purpose | Options |
|---|---|---|
| `setup` | download databases and models; writes `db_manifest.json` | `--db-dir`; `--bakta light\|full\|skip`; `--genomad/--no-genomad` (non-commercial use); `--antismash/--no-antismash`; `--promotech-models/--no-promotech-models`; `--go/--no-go`; `--kegg/--no-kegg` (academic use); `--cog/--no-cog`; `--contaminants/--no-contaminants`; `--minimal` (no geNOMAD, no antiSMASH); `--check`; `--images [all]`, `--no-images`; `--profile auto\|apptainer\|docker`; `--cache-dir`; `--install-uniop-from-git`, `--uniop`, `--no-uniop`; `--nextflow/--no-nextflow`; `--accept-download-plan` (accepts no licence terms) |
| `preflight` | validate inputs and installed tools | `--genome` (required); `--proteins`; `--peptides`; `--abundance`; `--annotation`; `--outdir`; `--strict` (warnings exit 3); `--profile` (`none`: inputs only) |
| `run` | run the analysis | `--genome` (required); `--outdir` (required); `--peptides`; `--proteins`; `--abundance`; `--annotation`; `--reference-peptides`; `--replicon-depth`; `--db-dir`; `--threads`; `--memory` (default 64.GB); `--modules core,mobile,bgc,function`; `--genome-relation related\|same`; `--abundance-missing median\|exclude\|min\|zero\|mean`; `--abundance-type ibaq\|top3\|maxlfq\|intensity\|unknown`; `--abundance-scale auto\|log2\|linear`; `--promoters targeted\|genome`; `--terminators terminatornet\|transtermhp`; `--add-bactermfinder`; `--publication` (10,000 permutations and bootstrap draws); `--seed`; `--parquet`; `--sensitivity basic\|full\|off`; `--strict`; `--profile` (standard, farmshare, farmshare_local, slurm, apptainer, conda, docker; default auto); `--resume`; `--executor nextflow\|local` |
| `makedb` | build a strain search database (no FDR is built on it) | `--annotation`; `--genome`; `--outdir` (required); `--db-dir`; `--six-frame-orfs`; `--min-orf-aa` (default 20); `--patch`; `--no-contaminants` |
| `convert` | search-engine output to the peptide table | `INPUTS…`; `--from diann\|casanovo\|instanovo\|cascadia\|generic` (required; casanovo and instanovo exit "not implemented"); `--out` (required); `--abundance-out`; `--source database\|denovo`; `--search-db same\|other` |
| `compare` | compare runs on the same genome | `RUNS…`; `--outdir` (required); `--labels`; `--allow-cross-major` |
| `report` | rebuild plots and report from tables | `OUT`; `--no-interactive`; `--plots`; `--strict` (size budget); `--redact` (shareable report) |
| `example fetch` | download, verify and unpack the example bundle | `--outdir`; `--version X.Y` (bundle version); `--from-file` (offline) |
| `demo` | bundled 100-kb simulated demo, no Nextflow, containers or network | `--outdir`; `--full` (through `gippa run` with the real tools); `--open` |

<!-- src: $GIPPA_ROOT/.. scratch capture of `gippa <cmd> --help` for run, setup, preflight, makedb,
convert, compare, report, example fetch, demo (gippa 1.1.1, this branch); methods.md §2
(Converters: Casanovo and InstaNovo exit 5 "not implemented") -->

**Table 2.** Documented output tables (`docs/reference/tables.md`, generated from the schemas).

| Location | Tables |
|---|---|
| `mapping/` (4) | `peptide_hits` (peptide × accepted locus), `protein_hits`, `mapping_fdr` (placement FDR strata, de novo precision bins), `unmapped` |
| `tables/` (18) | `genes` (one row per CDS, sORF or pseudogene), `detection_summary` (tier, FDR curve, O − E, QC), `regions_stats`, `autocorrelation`, `islands`, `operons`, `pairs`, `operon_flags`, `pileup`, `variants`, `novel_orfs`, `category_summary`, `enrichment_detection` (ORA), `enrichment_gsea`, `determinants_features`, `determinants_model`, `determinants_attribution`, `sensitivity` |
| `gippa compare` (3) | `compare_detection`, `compare_abundance` (Spearman ρ with TU-block bootstrap CI, Bland–Altman), `compare_stability` (Jaccard indices) |
| `plots/support/` (9; not result tables) | `replicon_landmarks`, `gc_skew`, `variant_windows`, `gsea_running`, `determinants_shap`, `determinants_pd`, `null_regions`, `null_islands`, `null_determinants` |

<!-- src: docs/reference/tables.md (Path and description lines of all 34 sections) -->

### 2.7 How GIPPA was built

GIPPA was written against a fixed specification, with architecture decision records (ADRs) that
override it once the project owner accepts them; a change to any default, threshold, tolerance or
gate needs such a record. Most of the code, tests and documentation, including this draft, were
written by AI coding agents working in separate roles (specification, tests, implementation,
review, statistics); every commit names its role, and continuous integration checks that no
commit changes both code and its tests and that protected files change only with an accepted ADR.
[PLACEHOLDER: Khoa: AI-assistance statement, wording and scope.]

<!-- src: methods.md §25; AGENTS.md; scripts/ci/check_roles.py; scripts/ci/check_protected.py -->

## 3. Error control

GIPPA reports every p- and q-value, passing or not. Each claim belongs to one family with a stated
error rate (Table 3), and no error rate is controlled across families. Where no valid estimate
exists, the row carries a flag (`ora_uncalibrated`, `null_uncalibrated`, `ci_unestimated`,
`class_error_unestimated` and others) instead of a number that looks like one.

**Table 3.** Error-control families and their status in GIPPA 1.1.1 (condensed from Supplementary
Methods Table M2).

| Claim | Family | Error rate | Level | Calibration status |
|---|---|---|---|---|
| Identification | the input search | search-engine q, inherited | 0.01 | not re-estimated by GIPPA |
| Peptide placement | per decoy type × difference stratum | weighted target–decoy FDR | 0.01 | list-level calibrated; contaminant attack is a known failure at q = 0.01; gate M3 not met |
| Gene detection | eligible genes (+ own units of unpaired decoys) | picked FDR (D + 1)/T, or MAYU | 0.01 | calibrated at 0.01; `two_peptide` and `protein_list` unestimated |
| Novel ORFs | per class and regime | class FDR from engine decoys | 0.01 | within tolerance in the gate; anti-conservative for multi-peptide false ORFs (known failure) |
| Variants | none | placement q of supporting peptides | 0.01 | gate M7 not met (F1 0.874 < 0.90) |
| Region tests | per region type × statistic | BH FDR over circular-shift p | 0.05 | within tolerance |
| Autocorrelation | 20 lags per replicon × statistic | BH FDR | 0.05 | at level with exact TUs only; anti-conservative with predicted TUs (flagged on every run) |
| Abundance islands | all windows, both directions | FWER (single-step maxT) | 0.05 | within a stated dependence regime |
| Operon discordance | both-detected operon pairs | BH FDR | 0.10 | claim narrowed to smooth, monotone abundance noise |
| Missed operons | non-operon pairs | decoy FDR | 0.10 | vacuous with the shipped decoys (no calls) |
| ORA | per collection × direction | BH FDR | 0.05 | shipped flagged `ora_uncalibrated` |
| GSEA | per collection | BH FDR | 0.05 | restricted claim (flags for clustered or confounded sets) |
| Determinants | one p per model | permutation p | – | shipped flagged `null_uncalibrated` |

<!-- src: methods.md Table M2 (§18); results.md §1.2 (status per family; numbers re-checked in
CAL111: identical records); docs/adr/ADR-0034, ADR-0035, ADR-0049, ADR-0055, ADR-0078 -->

### 3.1 Placement FDR

A peptide can be a correct identification and still be placed at the wrong locus, or be a wrong
identification that happens to match the genome. GIPPA estimates the share of accepted placements
that are chance hits of a wrong peptide with genome decoys: every translated frame reversed for
database-search input, and shuffled within stop-free segments for de novo input, because de novo
errors often reverse a stretch of sequence. Peptides are ranked by a chance-hit score, lists are
formed per decoy type and number of differences m ∈ {0, 1, 2}, and within a list
FDR(e) = (W(e) + c)/T(e), where W sums decoy hits weighted by an estimate of the probability that
each is wrong, and c is a +1-type correction [Barber and Candès 2015; He et al. 2015]. With c = 1 a
list's q-values cannot fall below 1/T, so a stratum with fewer than about 100 targets cannot accept
anything at 0.01; since 1.1.1 the run warns when this happens (ADR-0079). Exact hits from a
same-strain search receive q = 0 and are flagged `same_db`; for such input GIPPA adds no
identification error estimate of its own. The genome decoys under-count chance hits of wrong
peptides that resemble real proteins, a limitation recorded in the methods note.

<!-- src: methods.md §5.4; docs/adr/ADR-0010, ADR-0019, ADR-0079; CHANGELOG.md 1.1.1 Features;
docs/methods/placement-fdr.md › Limitations -->

### 3.2 Gene-level FDR

In the `decoy_picked` tier GIPPA adapts the picked protein FDR [Savitski et al. 2015] to genes:
in each gene, the best accepted target peptide competes with the best decoy peptide paired to it,
and FDR(t) = (D(t) + 1)/T(t) over target-won genes T and decoy-won units D. Release 1.1.1 changed
how decoys enter (ADR-0078). In 1.1.0, decoys that could not be paired to one gene were dropped
and each decoy win was rescaled by 1/f, where f is the paired share. That assumes decoys are lost
at random. When protein identifiers are not locus tags and genes are paired only by the votes of
accepted peptides (for example a search against a related strain), a decoy pairs only if its
protein already has an accepted target, so exactly the decoys that would have won their gene are
lost. The campaign found this (Section 4.4): on an NCM3722 search mapped to the MG1655 genome,
1.1.0 reported 2,399 genes at maximum q 0.0051 while the realised FDR, estimated with all decoys,
was about 0.026–0.028. In 1.1.1 every q-passing decoy competes, in its own unit if it maps to no
gene or to several, and the 1/f weight is retired; the same run then reports 2,305 genes at
maximum q 0.0100. The MAYU tier [Reiter et al. 2009] uses a closed-form expectation of false
genes from the peptide-level error, with no +1 term and no finite-sample guarantee.

<!-- src: methods.md §6 (Picked gene FDR, Changes in 1.1.1, MAYU); docs/adr/ADR-0078 (Context:
N2 2,399 genes, max q 0.0051, realised about 0.026-0.028); P2S §4 row #606 (fixed: 2,305 genes at
max q 0.0100); CAMP/ecoli/REPORT_1b.md (N1 baseline, Rerun) -->

### 3.3 Context statistics under dependence

Gene-level tests along a genome are dependent: genes share TUs and abundance is autocorrelated.
GIPPA's region tests therefore use circular shifts. The region's genes are moved together around
the replicon, which keeps its size, its internal structure and the autocorrelation of the
background, and the studentised observed-minus-expected statistic is recomputed at every shift
(Supplementary Methods §9). The p-values of one region type and statistic form a BH family
[Benjamini and Hochberg 1995]. Confidence intervals for region effects use a TU-block bootstrap.
Abundance islands are scanned over windows of TUs with a block bootstrap, and the maximum
statistic over all windows of all replicons gives a family-wise error rate. Its validity depends
on the block length matching the dependence; the claim is restricted to a stated regime
(ADR-0021 Amendment 1) and the run reports the block length and lag-1 Moran's I as diagnostics,
but no rule turns them into an inside-or-outside verdict.

<!-- src: methods.md §9-10; docs/methods/region-tests.md; docs/methods/islands.md; results.md §1.2
(islands domain, no output flags it) -->

### 3.4 Operons, function and determinants

Operon discordance compares the abundances of the two genes of a predicted operon pair with a
per-pair t statistic whose noise model depends on abundance, BH-adjusted at 0.10. Missed operons
are non-operon pairs whose partners behave like operon partners; their FDR is estimated with
convergent and sandwich decoy pairs. ORA uses a TU-level conditional Poisson-binomial null with
detectability as covariate; its p-values are shipped flagged `ora_uncalibrated` because the null
fails under realistic TU errors (ADR-0035). GSEA runs R fgsea 1.36.2 as a separate process over
TU-level ranks (ADR-0038); its claim is restricted to sets that are not genomically clustered and
is conditional on detection (ADR-0025, ADR-0026). The determinants permutation p-value fails the
dependence attacks and is shipped flagged `null_uncalibrated` (ADR-0034). A flag does not remove a
row from its family.

<!-- src: methods.md §11, §14, §15, §18; docs/adr/ADR-0025, ADR-0026, ADR-0034, ADR-0035, ADR-0038 -->

## 4. Validation

### 4.1 Simulation calibration and attack nulls

Every statistical method was calibrated on data from gippasim, a simulator written within the
project that plants genomes, TUs, prophages, abundance with replicon, region, island and AR(1)
effects, detection by a logistic model, and search and de novo errors at known rates. Every claim
in this section is therefore conditional on gippasim's generative model. The criteria are fixed in
a protected file: at nominal levels 0.01, 0.05 and 0.10, the observed rejection rate must not
exceed α + 2√(α(1−α)/n), the mean false discovery proportion (FDP) must not exceed
q + 2·SD(FDP)/√R, the family-wise rate must not exceed α + 2√(α(1−α)/R), and null p-values must
pass a two-sided Kolmogorov–Smirnov test at p > 0.001. A pass rules out excess error larger than
two Monte Carlo standard errors, not excess within that margin; for α = 0.05 the bound is about
1.5α at R = 300 and 1.3α at R = 1,000. A check that exceeds its bound in one block of R replicates
is re-run on a fresh block and fails only if both exceed it (ADR-0015).

<!-- src: tests/calibration/tolerances.yaml (levels, se_multiplier 2.0, ks_min_p 0.001);
methods.md §21 (criteria, two-block confirmation); results.md §1.1 (bounds 1.5α at R 300, 1.3α at
R 1,000) -->

The full tier for release 1.1.1 ran 15 suites on commit `d2bd17e4` (Slurm array 2022338, 25
November 2026), with 100 to 1,000 replicates per suite. All 15 tasks exited 0: of 558 tests, 462
passed, 90 were skipped (shard and merge placeholders) and 6 were expected failures, with none
failed. The 284 attack tests gave 276 passed, 3 skipped and 5 expected failures. Of the harness
records at levels 0.01–0.10 (Fig. 2b), 1,208 passed, 5 were unconfirmed (over the bound in one
block of two), 14 failed in arms that are recorded only or whose rows ship flagged, and 17 failed
as known, strict expected failures. Null p-values were uniform where the claim applies (Fig. 2a),
and observed error stayed at or below the bound in the representative streams of Table 4 and
Fig. 2c. The records equal those of the 1.1.0 full tier except for the ADR-0078 streams: the
decoy-loss attack became more conservative (mean FDP 0.0072 against 0.0098 at q = 0.01 with 20 %
of decoys lost) and the new vote-route attack passed (0.0098 at q = 0.01, bound 0.0103).

<!-- src: RD row R3 (calibration_full: 15 suites, 558 tests, 462/90/6; array 2022338), row R4
(284 attack tests: 276/3/5); local/logs/nightly.tsv:194-195; Fig. 2b counts: fig2_calibration.py
stdout {"pass": 1208, "unconfirmed": 5, "fail_flagged": 14, "fail_xfail": 17} over CAL111
lanes/*/*.jsonl + *.merged.jsonl, last block per (file, stream, criterion, level), criteria level,
fdr, fwer, fdr_counts; decoy loss: CAL111 lanes/gene_fdr-group*/gene_fdr.jsonl
gene_fdr.attack_decoy_loss.L20 q 0.01 observed 0.00717 (CAL110 0.00983); attack_decoy_pairing.reid100
q 0.01 observed 0.00977 bound 0.01033 -->

**Table 4.** Observed error against nominal level, representative streams of the 1.1.1 full tier
(α or q = 0.05 unless stated). Bounds are nominal + 2 Monte Carlo SE (FDR rows: empirical SD of the
FDP). R replicates; m pooled tests per replicate.

| Method | Stream | Criterion | Observed | Bound | R (m) |
|---|---|---|---|---|---|
| placement | read-through, shuffled decoys | FDR | 0.0355 | 0.0510 | 300 |
| gene FDR | picked, π = 0.05, q = 0.01 | FDR | 0.00941 | 0.01031 | 300 |
| gene FDR | MAYU (PEP), π = 0.05, q = 0.01 | FDR | 0.00962 | 0.01022 | 300 |
| region tests | HL shift, chromosome pooled | level | 0.0505 | 0.0569 | 1,000 (4) |
| region tests | prophage HL shift | FDR | 0.063 | 0.0654 | 1,000 |
| autocorrelation | Moran's I lag 1, exact TUs | level | 0.055 | 0.0638 | 1,000 |
| autocorrelation | same, UniOP-like TU errors (recorded only) | level | 0.071 | 0.0638 | 1,000 |
| islands | AR(1) ρ = 0 / 0.5 / 0.9 | FWER | 0.022 / 0.020 / 0.024 | 0.0695 | 500 |
| islands | ρ = 0.9, short lag, 1,800 genes (outside the claim; recorded) | FWER | 0.130 | 0.0695 | 500 |
| operon discordance | gate null, 1,800 genes | level / FDR | 0.0499 / 0.034 | 0.0510 / 0.0615 | 1,000 (200) |
| missed operons | shipped decoys | FDR | no calls (vacuous) | – | 300 |
| GSEA | pooled null, R engine | level | 0.0415 | 0.0525 | 300 (100) |
| novel ORFs | class FDR, database regime | FDR | 0.0453 | 0.0511 | 300 |
| novel ORFs | multi-peptide false ORFs (known failure) | FDR | 0.0788 | 0.0511 | 300 |
| variants | F1 of `high` ALT calls (≥ 0.90; known failure) | F1 | 0.874 | 0.90 | 100 |
| compare | ρ CI non-coverage, TU block | level | 0.048 | 0.0638 | 1,000 |

<!-- src: results.md Table R2 and Table R3b (CAL110 line references), re-read in CAL111 at the
same lines (identical values): lanes/placement-group1/placement.jsonl:2; lanes/gene_fdr-group1/
gene_fdr.jsonl:44,47; lanes/regions-group2/regions.jsonl:21; lanes/regions-group0/regions.jsonl:8;
lanes/autocorrelation-group0/autocorrelation.jsonl:77,44; lanes/islands-group0/islands.jsonl:32,35,38;
lanes/islands-group0/islands_recorded.jsonl:2; lanes/discordance-group1/discordance.jsonl:17,21;
lanes/missed_operons-group0/missed_operons.jsonl:35; lanes/gsea-group0/gsea.jsonl:45;
novel_orfs.merged.jsonl:38, 23; variants.merged.jsonl:3; lanes/compare-group0/compare.jsonl:41 -->

The failures matter as much as the passes, and we state them in the paper rather than only in the
supplement. (i) The contaminant attack `flat_c200` exceeds the placement FDR bound at q = 0.01
(0.0150 against 0.0102). (ii) When a false novel ORF carries several wrong peptides, the class FDR
is anti-conservative (0.0788 at q = 0.05), because an ORF's q is its best peptide's q
(ADR-0050). (iii) Variant calls miss the F1 ≥ 0.90 target (0.874), and the isobaric-swap attacks
give family-wise rates of 0.97 and 1.0; gate M7 was recorded as not met. (iv) TU-block Moran's I
holds level only with exact TUs; with UniOP-like TU errors it rejects at 0.071, and under the
distance-model fallback at 0.354, at α = 0.05; every such row ships flagged `tu_block_tu_error`.
(v) Operon discordance fails under a step-shaped abundance-noise attack; its claim was narrowed to
smooth, monotone noise after that attack failed (ADR-0055 Amendment 1). (vi) The island FWER is
breached at ρ = 0.9 with short block lengths (0.096 and 0.130 at 1,000 and 1,800 genes); the claim
was likewise narrowed after the breach. (vii) The shipped missed-operon decoys made no call in
any of 300 power replicates, so its pass is vacuous. (viii) ORA and the determinants permutation
test pass their gates only because their rows are flagged. All of these were known before
release and are listed in the readiness record as expected failures or flags.

<!-- src: results.md §1.2-1.3 and Table R3b (values identical in CAL111: placement.merged.jsonl:43;
novel_orfs.merged.jsonl:22-24; variants.merged.jsonl:1-3; autocorrelation 44, distance-model arm;
discordance D-1b-step; islands_recorded.jsonl:1-2; missed_operons.jsonl:35); docs/adr/ADR-0021
Amendment 1, ADR-0049, ADR-0050, ADR-0055 Amendment 1, ADR-0034, ADR-0035 -->

### 4.2 True-example assertions on real genomes

Eight true-example assertions (TE1–TE8) check biology that any correct run must reproduce
(Table 5). They are evaluated on the Nissle pilot and on weekly *E. coli* K-12 and Nissle runs.
TE1 passed with 28,129 of 28,129 A2 peptides mapping exactly in frame, but every A2 peptide is a
substring of the searched database, so TE1 checks that database, genome and mapping agree, not
identification error. TE2 passed vacuously: a same-database search cannot yield an amino-acid
variant, and the weekly in-silico proxy also had no variant site. TE4 failed: the ratio of the
ribosomal stalk protein L7/L12 (RplL) to L10 (RplJ), expected near 4:1 (log2 in [1, 3]), was
−0.48 on iBAQ, while AtpA/AtpD passed. The cause was not established; release 1.1.1 waives readiness
row R5 for this pair only (ADR-0080). The campaign later showed the same failure in five of nine lineages where the pair could be
tested (Section 4.4).

<!-- src: results.md Table R4 and §3; PIL › te1, te_real.te2, te_real.te4.pairs; RD row R5
(TE1-TE3, TE5, TE7, TE8 pass; TE4 fail; waived, ADR-0070 via ADR-0080);
local/ci/validation/2026-11-25_te8_2020593.json, 2026-11-25_te2_2020593.json (1.1.1 weekly,
d2bd17e4); P2S §2.5 (TE4: fails in V. natriegens, M. xanthus, C. difficile, M. tuberculosis and
E. coli; passes in C. crescentus, B. theta, Synechococcus, 6 of 7 T. pallidum strains) -->

**Table 5.** True-example assertions (pilot values from the 1.0.0-era freeze; TE6 and TE8 from the
1.1.1 weekly run, commit `d2bd17e4`).

| TE | Assertion | Outcome | Value |
|---|---|---|---|
| TE1 | ≥ 99.5 % of A2 peptides map exactly, in frame | pass | 28,129 / 28,129 (Wilson 95 % lower bound 0.99986) |
| TE2 | no `high` ALT call off start or stop codons | pass, vacuous | 0 ALT sites in A2 and in the weekly proxy |
| TE3 | ribosomal proteins ≥ 90 % detected, median abundance in the top quartile | pass | 53 / 53 detected [0.932, 1.0]; median percentile rank 0.987 |
| TE4 | log2(RplL/RplJ) ∈ [1, 3]; log2(AtpA/AtpD) ∈ [−1, 1] | **fail (waived, R5)** | RplL/RplJ −0.480 (Top3 −0.531); AtpA/AtpD −0.118 |
| TE5 | trypsin and keratin peptides land in the contaminant bucket | pass | A1 59 / 59; A2 61 / 61 |
| TE6 | K-12: Bakta oriC next to mnmG | pass (any-call rule) | 2 oriC calls, 1 at the true origin |
| TE7 | same seed gives byte-identical tables | pass (tables only) | 2,786 gated files, 0 differ; integration test on a simulated genome passes |
| TE8 | UniProt K-12 peptides ≥ 99.5 % exact on K-12; Nissle placement FDR ≤ 1 % | pass | 0.99900 (151,272 / 151,424); Nissle placement FDR 0.00533 |

<!-- src: results.md Table R4; local/ci/validation/2026-11-25_te8_2020593.json › shares.k12_exact
0.99899, counts 151272/151424, shares.nissle_placement_fdr 0.005331; local/logs/gippa-weekly.2020593_0.out:319
("TE6: 2 oriC calls, at the true origin: 1"); RD row R6 -->

### 4.3 The Nissle 1917 pilot

The pilot was specified before it ran. Three wild-type *E. coli* Nissle 1917 DIA runs (PXD078459,
timsTOF HT) were searched with DIA-NN 2.3.1 in three arms: A1 against the UniProt K-12 proteome
(a generic reference; Nissle belongs to phylogroup B2, K-12 to A), A2 against the Nissle strain
database built by `gippa makedb` from the wild-type assembly (GCF_021559835.1), and A3, an
exploratory de novo arm (Cascadia 0.0.7). GIPPA ran end to end on all three arms. The pilot
numbers below were measured at the 1.0.0 report-freeze commit; the values that depend on gene
detection were recomputed with the 1.1.1 detection code on the frozen inputs, and a control with
the unchanged freeze code reproduced the recorded tables byte for byte (Slurm jobs 2002006 and
2001992). Whether the rest of the 1.1.1 code reproduces the remaining pilot values was not
checked.

<!-- src: methods.md §23; results.md §4 (arms, provenance); validation/pilot/recheck_1_1_1/README.md
(control 2001992 byte-identical; fix 2002006); results.md Provenance (124 files changed between
fa49b8f1 and 1.1.0; not re-run) -->

**Detection.** Of 4,685 strain genes, A2 detected 3,051 (0.651 [0.638, 0.665], gene q ≤ 0.01) and
A1 2,645 (0.565 [0.550, 0.579]) (Fig. 3a). The 916 strain-specific genes, which share no tryptic
peptide with K-12, are undetectable in A1 by construction; A2 detected 316 of them (0.345
[0.315, 0.376]). Plasmid genes were detected in A2 at 0.769 [0.497, 0.918] against 0.077 in A1,
and prophage genes at 0.272 [0.204, 0.352] against 0.022. The Wilson intervals treat genes as
independent, whereas genes of one TU co-vary, so they are optimistic; they describe this sample,
not a population. The A2 search's identification error is the search engine's, inherited at
q ≤ 0.01; GIPPA's placement FDR does not apply to A2 because every A2 peptide is an exact hit of the
searched database, so readiness row R2 (mapping FDR ≤ 1 %) was not measured and is waived.

<!-- src: PIL111 › m3.A2/A1.subsets (all, strain_specific, plasmid, prophage); results.md §4.2
and Table R6; RD row R2 (waived, ADR-0071 via ADR-0080) -->

**Pre-registered hypotheses.** H1 stated that A2 detects more strain genes than A1, mostly plasmid,
prophage and BGC genes. The "more" clause holds: 434 genes were detected only in A2 and 28 only in
A1, and a TU-level test gave 290 against 24 units (p = 1.8 × 10⁻⁵⁹). The "mostly" clause fails:
only 106 of the 434 A2-only genes (0.244 [0.206, 0.287]) lie in plasmids, prophages or BGCs, against
the pre-registered requirement of more than half (Fig. 3b). Most A2-only genes (316 of 434) are the
strain-specific genes that A1 could not detect, so the "more" clause was expected by design. H1 is
therefore recorded as not supported. H2 (A1's genome-inconsistent FDP is at least its entrapment
FDP) needs an entrapment benchmark and was deferred. H3 (context conclusions agree between arms)
could not be tested as stated: neither arm called an island or an operon flag, and region tests
were not compared; the recorded "supported" rests on substitute criteria (abundance Spearman
ρ = 0.978 [0.973, 0.983] over 2,613 genes, and a Jaccard index of 0.588 for enriched categories that
include flagged ORA rows), so we do not count H3 as confirmed. H4 (de novo precision agrees with
precision against A2) was consistent with the lower bound in all 10 score bins; the check can
refute but not confirm.

<!-- src: PIL111 › h1 (b 434, c 28, b_tu 290, c_tu 24, p_tu 1.80e-59, mostly 106/434, strain_specific
316/434, verdict not_supported), m4 (rho 0.978173 [0.972762, 0.982741], n_shared 2613), m5
(enriched_categories jaccard 0.587699); results.md §4.4 and Table R7 (H2 deferred; H3 not tested
as stated; H4 consistent with lower bound, 10 bins) -->

**Free entrapment.** K-12 peptides that are absent from the Nissle genome cannot be true if the
genome is complete and the sample is Nissle, so they serve as free entrapment for A1 [Wen et al.
2025]. At precursor level (28,347 targets, 157 entrapment hits) the combined FDP estimate was
0.0236, with a realised lower bound of 0.0055 and a one-sided 95 % upper bound of 0.0265 (Fig. 3c).
The point estimate exceeds the 1 % target, but the lower bound does not, so the assessment is
inconclusive. The upper bound assumes independent, equal-chance false matches; in simulations with
clustered false matches it covered only 0.72–0.75 of the time. Two biases run in opposite
directions and are not quantified: near-homologous entrapment peptides push the estimate up, and
DIA-NN's use of protein information in its q-values pushes it down.

<!-- src: PIL › m8.precursor (N_T 28347, N_E 157, fdp_combined 0.023607, fdp_lower 0.005508,
ci.ucb95 0.026523, assessment inconclusive); results.md §4.5; CAL111 lanes/compare-group0/
compare_attack_free_entrapment_recorded.jsonl:4,7 (D 0.754, E 0.719) -->

**Novel ORFs, variants, de novo.** A2 called no novel ORF and no variant, as a same-database search
allows none. A1 called one novel ORF, which the pilot report attributes to an annotation artefact.
The A3 de novo arm passed its feasibility rule on one run (549 of 863 high-scoring de novo peptides,
0.636 [0.604, 0.668], matched the genome's six-frame translation, against a shuffled-decoy share of
0.006); its variant class error rates are unestimated because DIA data give no spectrum pairing.

<!-- src: results.md §4.6-4.7; validation/pilot/metrics/a3_feasibility.json -->

### 4.4 A multi-organism campaign on public DIA data

To see how GIPPA behaves beyond *E. coli*, we ran it on public DIA datasets with a fixed recipe and
answer keys written before the results. Phase 1 used *E. coli* K-12 MG1655 (PXD058422, four runs,
Orbitrap Eclipse) with route A (Bakta proteins) and route B (plus six-frame ORFs), and a de novo
arm. Phase 1b searched an NCM3722 dataset (PXD060434) against its own genome and against MG1655, to
test strain-variant calls. Phase 2 ran ten datasets from ten lineages (Table 6), from a 1.14-Mb
genome grown in rabbit cells (*Treponema pallidum*, seven strains) to a 9.1-Mb genome
(*Myxococcus xanthus*) and two archaea. Each dataset had a pre-registered answer key, committed
before its report. Datasets were fetched from PRIDE and other ProteomeXchange members under an
accepted ADR (ADR-0077), Thermo raw files were converted with ThermoRawFileParser, and every run
used DIA-NN 2.3.1 with the pilot's settings, the GIPPA 1.1.1 release code (`bfc74ff6`), or, for the
gapped *T. pallidum* assemblies, main `38e1cdba` (1.1.1 plus three assembly-gap fixes released with
it). The campaign reports hold aggregates only.

<!-- src: "GIPPA demo and discovery campaign agent handover plan.md" (goal, recipe, Phase 2 table);
CAMP/datasets.tsv (PXD058422, Orbitrap Eclipse); CAMP/ecoli/REPORT.md (E1/E2/E4; 4 runs);
CAMP/ecoli/REPORT_1b.md (PXD060434 NCM3722); P2S header (Code: bfc74ff6; 38e1cdba for T. pallidum);
docs/adr/ADR-0077-campaign-sources-and-raw-conversion.md; git log of CAMP/*/ANSWER_KEY.md before
REPORT.md (e.g. c16fa3ec 2026-11-19 before 7702782e 2026-11-20; 854c143c before 13f85cc6) -->

**Table 6.** Campaign datasets and headline results (route A unless stated; shares with Wilson
95 % CI, optimistic because genes and peptides co-vary).

| Organism (genome) | Data | Runs, instrument | Genes detected / coding (share) | Route-A gene FDR at q ≤ 0.01 | Route B keeps | Novel-ORF candidates (called) |
|---|---|---|---|---|---|---|
| *E. coli* K-12 MG1655 | PXD058422 | 4, Orbitrap Eclipse | 2,691 / 4,368 (0.616 [0.602, 0.630]) | 0.0093 | 0.984 [0.979, 0.988] | 124 (0) |
| *T. pallidum*, 7 strains | PXD065022 | 3 per strain, timsTOF Pro | 668–880 / 972–976 (0.685–0.903) | 0.0011–0.0060 | 0.973–0.999 | 200 (0) |
| *Caulobacter crescentus* NA1000 | PXD069004 | 3, Exploris 480 | 3,005 / 3,723 (0.807 [0.794, 0.819]) | 0.0093 | 0.971 [0.964, 0.976] | 164 (0) |
| *Bacteroides thetaiotaomicron* | PXD055395 | 16, timsTOF Pro | 2,244 / 4,931 (0.455 [0.441, 0.469]) | 0.0098 | 0.972 [0.964, 0.978] | 99 (0) |
| *Synechococcus* CC9311 | PXD072227 | 5, Q Exactive | 1,897 / 2,859 (0.664 [0.646, 0.681]) | 0.0084 | 0.968 [0.960, 0.975] | 117 (0) |
| *Vibrio natriegens* | PXD049476 | 4, Exploris 480 | 2,515 / 4,551 (0.553 [0.538, 0.567]) | 0.0099 | 0.969 [0.961, 0.975] | 78 (0) |
| *Myxococcus xanthus* DK1622 | PXD060688 | 4, Exploris 480 | 4,864 / 7,318 (0.665 [0.654, 0.675]) | 0.0097 | 0.985 [0.981, 0.988] | 265 (0) |
| *Clostridioides difficile* 630Δerm | PXD065346 | 3, DIA | 1,928 / 3,824 (0.504 [0.488, 0.520]) | 0.0099 | 0.988 [0.982, 0.992] | 53 (0) |
| *Mycobacterium tuberculosis* H37Rv | PXD057335 | 3, Exploris 480 | 2,942 / 4,075 (0.722 [0.708, 0.736]) | 0.0099 | 0.977 [0.971, 0.982] | 247 (0) |
| *Haloferax volcanii* H26 (DS2 genome) | PXD066008 | 5, timsTOF Pro 2 / HT | 2,978 / 4,023 (0.740 [0.726, 0.754]) | 0.0097 | 0.972 [0.966, 0.977] | 142 (0) |
| *Saccharolobus islandicus* REY15A | PXD079213 | 3, Orbitrap Astral | 2,109 / 2,807 (0.751 [0.735, 0.767]) | 0.0085 | 0.983 [0.976, 0.988] | 92 (0) |

<!-- src: P2S §2.2 (gene FDR), §2.3 (runs, instrument, detected/coding, Wilson CIs), §2.4 (route B
kept), §2.6 (candidates, 0 called); CAMP/ecoli/REPORT.md Key numbers (E1 2,691, FDR 0.0093 1.1.1
recheck 0.00929), row 14 (98.4 % [97.9, 98.8]), row 9/15 (124 candidates, 0 called); E. coli
share 2,691/4,368 Wilson computed by fig4_campaign.py (0.616 [0.602, 0.630]); T. pallidum per
strain: PXD065022/metrics.json › T2.per_strain and strains.*.metrics.S9_route_b (fig4_campaign.tsv);
accessions: handover plan Phase 2 table -->

**Mapping and detection.** Across the 16 Phase 2 genomes, all 302,396 accepted route-A peptides
mapped exactly in frame at their own gene, with no miss of any kind; the lowest Wilson lower bound
was 99.86 % (*T. pallidum* Grady, 2,746 peptides). The same held for *E. coli* (31,071 of 31,071).
Every run used the `decoy_picked` tier, and every route-A gene FDR estimate at q ≤ 0.01 was at most
0.0099. Detection depth ranged from 0.455 of coding genes (*B. thetaiotaomicron*, 11.5-min
gradients) to 0.903 (*T. pallidum* Haiti B) (Fig. 4a). Ribosomal proteins (TE3) passed in every
dataset. In the two archaea, Bakta names almost no ribosomal protein, so this check used RefSeq
names carried to Bakta genes by a shared stop codon, as pre-registered. Searching six-frame ORFs
(route B) increased the database 18- to 34-fold in records, but kept at least 96.8 % of route-A
genes in every dataset (Fig. 4b); the extra cost fell on the DIA-NN search (1.5–5.2× longer), not on
GIPPA, whose run times were similar for the two routes except in *M. tuberculosis*.

<!-- src: P2S §2.1 (302,396; 0 misses; 99.860 % Grady 2,746), §2.2 (decoy_picked in 34 runs;
≤ 0.0099), §2.3 (range; archaea TE3 via RefSeq names), §2.4 (records 18.4x-33.5x; kept ≥ 96.8 %;
DIA-NN 1.5x-5.2x; gippa run similar except Mtb); CAMP/ecoli/REPORT.md row 1 (31,071) -->

**Answer keys.** Of 158 pre-registered Phase 2 checks, 91 passed, 10 failed, 34 were not testable
or not run and 23 were reported without a threshold (Fig. 4d). Five of the ten failures were the
TE4 ratio (Fig. 4c). With *E. coli*, RplL/RplJ failed in six of 15 genomes where it could be
evaluated, while AtpA/AtpD passed in every genome where it could be tested. Pass or fail did not
follow peptide counts, and the verdict changed with the abundance measure near the band edges; one
candidate explanation, that the N-terminally acetylated L7 form escapes the search, was not
tested. Of the other five, one was a wrong premise of the answer key (a selenoprotein tail that Bakta
annotates as a separate gene), one a vector caveat stated in the key in advance (pyrE/pyrF
detected in *S. islandicus*), one a test-design limit (a prophage region that spans a whole
plasmid, which a circular-shift null cannot vary, plus small regions), and two were *T. pallidum*
checks: pseudogene distal parts (one half untestable by construction) and a strain abundance
comparison missed on two small differences with error unestimated. Pre-registered checks that passed included the absence of pHV2 proteins (0 of 8 genes) and the
pyrE2 check in *H. volcanii* H26, the ESX-1 and PE/PPE checks in *M. tuberculosis* (PE/PPE genes
detected at 12.7 % against 74.5 % for other genes, Fisher p = 1.1 × 10⁻⁵⁴), under-detection of prophage
regions in *V. natriegens* and *M. xanthus* (all four regions O < E, q ≤ 0.05), and detection
differences between the two *V. natriegens* chromosomes (65.3 % against 38.0 %, difference 27.3
points [24.4, 30.2]). Across the 13 prophage regions tested in four genomes, detection was below
expectation in all 13 and significant in 6.

<!-- src: P2S §1 (verdict table; total 91 / 10 / 34 / 23 of 158; fail list), §2.5 (TE4 table;
AtpA/AtpD 8 of 8 plus E. coli), §2.7 (prophage rule; 13/13 direction, 6/13 significance; V.
natriegens chromosomes), §3 (premise errors and design limits); fig4_campaign.tsv panel c
(15 genomes, 6 fail) -->

**Strain variants.** At default settings GIPPA did not recover substitutions between close strains.
In Phase 1b, 11 of 31 NCM3722–MG1655 substitution sites were covered by peptides, but none was
accepted (recall 0/11 [0, 0.26]); in the *T. pallidum* cross-strain arms, 0 of 159 covered sites
were accepted. The reason is the floor of the placement FDR: tolerant strata held 1 to 51 targets,
so their q could not fall below 1/T > 0.01 even with no decoy hit. Exact placements in the same
arms had pooled placement FDR estimates of 0.35 % (NCM3722) and 0.06–0.10 % (*T. pallidum*).
Release 1.1.1 adds a run warning for this floor (ADR-0079) but does not change acceptance.

<!-- src: CAMP/ecoli/REPORT_1b.md verdicts 1b-1 (0/11 [0.00, 0.26]), 1b-5 (0.35 %), "Why recall is
zero"; P2S §2.9 (0 of 159; strata 1-51 targets; floors 0.030-0.143), §2.2 (T. pallidum placement
FDR 0.059-0.101 %); docs/adr/ADR-0079 -->

**Novel ORFs.** No novel ORF was called in any of the 17 genomes at the default α = 0.01 (1,581
candidates; Fig. 5). The reason is structural: the class FDR is estimated within each class from
engine decoys, and nearly every class held fewer than 100 candidates, so its q could not reach
0.01; every route-B run carried the `class_q_floor_above_alpha` warning. The two classes with more
than 100 candidates, the N-terminal extensions of *M. xanthus* and *M. tuberculosis*, came closest
(minimum class FDR 0.0127 and 0.0118). N-terminal extensions had the lowest class FDR in every
dataset and mostly matched RefSeq starts upstream of Bakta's; antisense candidates had the highest
(0.14–1.0) and rarely matched RefSeq. We did not loosen α.

To see where true signal sits on this q scale, we defined positive controls after the runs but
before looking at other candidates: candidates that known biology confirms, chiefly the N-terminal
part of RF2 (prfB), which is made by a programmed +1 frameshift that Bakta does not model, plus
the infC late start, the PrdB selenocysteine read-through, and RefSeq proteins that Bakta missed
(Table 7). The three primary prfB controls had class FDRs of 0.111 to 0.167, so we used the largest,
τ_prfB = 0.167, as an exploratory reporting tier. It admitted 581 of 1,581 candidates (15
controls). Of the 566 non-control candidates, 423 (74.7 % [71.0, 78.1]) overlapped a RefSeq gene on
the same strand and replicon; GIPPA's own class FDRs imply about 75 false candidates among the 581
(12.9 %), an estimate that assumes the class-wise decoy model holds in these small classes and is
unchecked. This tier is a campaign analysis, not a GIPPA call: the threshold was chosen from the
same data, it is about 17 times α, the controls are biased toward abundant proteins, and the
three prfB controls themselves lack a start codon, so GIPPA's class rule would reject them at any
q. The candidates in this tier are leads for review, not discoveries.

<!-- src: P2S §2.6 (0 called; classes < 100; M. xanthus 0.0127, M. tuberculosis 0.0118; nterm lowest,
antisense 0.14-1.0), §6 (alpha not loosened); NPC §2 (controls, rules defined before passers were
looked at), §3 (τ_prfB 0.1667; 15/23 controls), §4.1 (581 pass; 566 non-control; 423/566 74.7 %
[71.0, 78.1]; 74.8 expected false, 12.9 %), §5 (caveats: same data, 17 × α, controls biased, no
start codon); Khoa 2026-11-24 request quoted in NPC header -->

**Table 7.** Positive controls for the novel-ORF tier (route B; class FDR is GIPPA's `class_fdr`;
none is called at α = 0.01).

| Control | Datasets (candidates) | Class | Class FDR |
|---|---|---|---|
| prfB N-terminal part (programmed +1 frameshift; primary) | *E. coli*, *B. thetaiotaomicron*, *V. natriegens* (3) | N-terminal extension | 0.143, 0.111, 0.167 |
| prfB leader | *E. coli*, *Synechococcus* (2) | intergenic | 0.143, 0.429 |
| infC late start (AUU codon) | *E. coli*, *C. difficile* (2) | N-terminal extension | 0.143 |
| PrdB Sec read-through and N-terminal part | *C. difficile* (2) | N-terminal extension | 0.143 |
| Bakta-missed RefSeq proteins: cdzC/D/I | *C. crescentus* (3) | intergenic | 0.077 |
| Bakta-missed ribosomal proteins L36, L34 | *M. tuberculosis* (2) | intergenic | 0.100 |
| Bakta-missed ribosomal proteins L37ae, L40e | *S. islandicus* (2) | out of frame, intergenic | 0.143, 0.167 |
| Bakta-missed ribosomal protein L28 (rpmB) | *T. pallidum*, 7 strains (7) | intergenic | 0.5–1.0 |

<!-- src: NPC §2 table "Controls found (23)"; CAMP/novel_pc/novel_pc_threshold.json › controls_found -->

**Defects found.** The campaign found four GIPPA defects, all fixed in 1.1.1 or with it: the
anti-conservative picked gene FDR on the vote route (Section 3.2; found in Phase 1b); crashes of
Promotech and Promoter Calculator on genomes with assembly gaps; and TerminatorNet coordinates that
shifted after the first non-ACGT base. The gap fixes change nothing on gap-free genomes. Route A
runs with locus-tag identifiers, such as *E. coli* E1 and E2, kept the same detected genes under
the 1.1.1 FDR rule; only their q-values moved slightly.

<!-- src: P2S §4 (#606/PR #607, #611/#612, #614/#616, #615/#617, ADR-0079); CHANGELOG.md 1.1.1
Fixes; CAMP/ecoli/REPORT.md "Recheck with 1.1.1 code" (same genes; E1 max q 0.00933 → 0.00929) -->

### 4.5 Release readiness

A release script evaluates ten acceptance rows from the specification and six housekeeping rows.
For 1.1.1 (commit `3f8793e0`, 27 November 2026) the overall verdict was pass, with the following
reading (Table 8). Three rows are waived, not passed: R2 because the A2 mapping FDR was not
measured, R5 because TE4 fails (for the RplL/RplJ pair only), and R7 because the measured runtime is
outside the specified range (Section 5). R10 (co-abundance against K-12 operons) is deferred by the
specification to a full validation. R3 and R4 pass without a waiver, but a pass means that no test
failed outside the listed expected failures; it does not mean that every method is calibrated in
every regime (Section 4.1). R9 passes because the pilot report records its five pass criteria as
met, with TE4 reported and explained but its cause not established.

<!-- src: RD (commit 3f8793e0, computed 2026-11-27T01:33:41Z, overall pass; rows R1-R10, H1-H6);
docs/adr/ADR-0080 (waivers R2, R5 TE4 only, R7; R3 not waived); RD rows R3, R4, R9 (reading of R3, R4,
R9) -->

**Table 8.** Release readiness of GIPPA 1.1.1 (acceptance rows).

| Row | Criterion | Status | Value |
|---|---|---|---|
| R1 | ≥ 99.5 % of A2 peptides map exactly | pass | 100.00 % |
| R2 | decoy-estimated mapping FDR ≤ 1 % | **waived** (ADR-0071, ADR-0080) | not measured (same-database search) |
| R3 | every statistical method passes its calibration suite | pass (no waiver; expected failures listed) | full tier 558 tests, 0 failed, 6 xfailed; fast tier 29,049 passed, 17 xfailed |
| R4 | every attack null stays at level | pass (expected failures listed) | 284 attack tests: 276 passed, 3 skipped, 5 xfailed |
| R5 | every true-example assertion the data allow passes | **waived** for TE4 only (ADR-0070, ADR-0080) | TE1–TE3, TE5, TE7, TE8 pass; TE4 fails |
| R6 | same seed gives byte-identical tables | pass | TE7; integration 967 passed; nf-test 144 tests |
| R7 | standard run on 16 threads and 64 GB ≈ 35 min (16–55) | **waived** (ADR-0072, ADR-0080) | 3 h 41 min |
| R8 | core modules ≈ 25 min | reported only | 3 h 01 min |
| R9 | the pilot passes | pass | 5 criteria; TE4 explained, cause not established |
| R10 | co-abundance against K-12 operons | **deferred** to full validation | – |

<!-- src: RD rows R1-R10 (status, value, threshold) -->

## 5. Performance

The standard-run benchmark was measured on the 1.0.0-era code on one FarmShare "wheat" node (Intel
Xeon Gold 6148) with 16 threads and 64 GB, on the K-12 genome with simulated peptide tables. A
standard run with all modules took 3 h 41 min (mean of two repeats; uncertainty unestimated), about
six times the specification's estimate of 35 min; core modules alone took 3 h 01 min. Two
annotation steps, Promoter Calculator and TerminatorNet, used about 92 % of summed task time in the
first repeat. Peak memory reached 62.4 GiB, 99.8 % of the allocation; whether that reflects a need
or a limit was not established. In the campaign, on 16 CPUs in one job, `gippa run` took 35 min
(*T. pallidum*, 1.14 Mb) to 8 h 10 min (*M. xanthus*, route A), and *M. xanthus* peaked at 60.9 GB.
`gippa demo` finishes in about 75 s without containers or databases.

<!-- src: benchmarks/baseline.json › standard_run.wheat.k12_standard, k12_core (13,232.5 s; 10,873.5 s;
2 repeats; peak 62.35 GiB); docs/adr/ADR-0072 (two repeats; not established); results.md §5
(92 % in two steps; demo 75 s); P2S §2.4 (gippa run A/B times; M. xanthus 60.9 GB) -->

## 6. Availability and licences

GIPPA 1.1.1 is MIT-licensed and installs with `pip install gippa==1.1.1`; the source, specification,
methods notes, ADRs, calibration suites, pilot aggregates and readiness records are at
https://github.com/khoa-yelo/gippa. Containers are pinned by SHA-256 digest. `gippa demo` runs a
bundled simulated example without network access, and `gippa example fetch` downloads a 1.09-MB
real example (two 200-kb regions of the Nissle 1917 genome with subsets of the pilot data) for five
input routes. [PLACEHOLDER: Zenodo DOI; documentation URL; Bioconda package.]

GIPPA never imports or vendors a third-party tool; CI fails if the source tree imports a GPL,
AGPL, non-commercial or unlicensed tool (Table 9). Users must respect the following licences:

- **GPL and AGPL tools run as separate programs.** Bakta, Promotech, Promoter Calculator, OSTIR,
  TransTermHP and BacTermFinder (GPL-3.0) and antiSMASH (AGPL-3.0) run as separate processes from
  pinned environments or container images; R (GPL-2) runs fgsea (MIT) as a separate process.
- **geNOMAD is licensed for non-commercial use only.** A commercial user must switch off the
  `mobile` module (or skip the database with `gippa setup --no-genomad`) or obtain a licence.
- **KEGG is for academic use.** Non-academic users need a KEGG licence, or can use the GO and COG
  collections only (`gippa setup --no-kegg`).
- **UniOP has no licence file.** It is fetched from its authors' repository at a pinned commit only
  when the user passes `--install-uniop-from-git` and so accepts its terms; it is never vendored or
  published as an image. Without it, GIPPA uses the distance-model fallback and flags the run.
- **DIA-NN is proprietary** (academic licence). GIPPA does not distribute it; the pilot and campaign
  used DIA-NN 2.3.1 Academia.

<!-- src: NOTICE.md (licence table and separate-process statement); scripts/ci/check_imports.py;
gippa setup --help (--genomad "non-commercial use", --kegg "academic use", --install-uniop-from-git
"UniOP has no licence file ... you accept its terms", --accept-download-plan "does not accept any
licence terms"); methods.md §20 (Tools and licences), §24; docs/preprint/outline.md §5;
validation/example/bundles.tsv (1.09 MB) via methods.md §20 -->

**Table 9.** Third-party tools in GIPPA 1.1.1.

| Tool | Version | Licence | Role | How run |
|---|---|---|---|---|
| Bakta | 1.12.2 (DB 6.0) | GPL-3.0 | genes, CRISPR, AMR, oriC | separate process |
| geNOMAD | 1.12.0 (DB 1.9) | non-commercial | prophages, plasmids | separate process |
| antiSMASH | 8.0.4 | AGPL-3.0 | BGCs | separate process |
| UniOP | git 3ed14dfc | no licence file | operon pairs | separate process; never redistributed |
| Promotech | git 56251ad9 | GPL-3.0 | promoters | separate process |
| Promoter Calculator | 1.2.4 | GPL-3.0 or later | transcription rates | separate process |
| TerminatorNet | git 72180309 | MIT | terminators | separate process |
| TransTermHP | 2.09 | GPL-3.0 or later | terminator cross-check | separate process |
| BacTermFinder (opt-in) | git 02487044 | GPL-3.0 | terminators | separate process |
| OSTIR | 1.1.3 | GPL-3.0 | translation initiation | separate process |
| miniprot | 0.18 | MIT | protein alignment | separate process |
| BLAST+ | 2.16.0 | public domain | protein alignment, contaminants | separate process |
| R / fgsea | 4.5.3 / 1.36.2 | GPL-2 / MIT | GSEA | separate process |
| codon-bias | 0.5.0 | MIT | codon indices | in process |
| DIA-NN | user's own | proprietary (academic) | spectrum search, outside GIPPA | not distributed |

<!-- src: NOTICE.md table; methods.md Table M3 -->

## 7. Limitations

**What the evidence covers.** All calibration rests on gippasim, a simulator built within the
project, so each calibration claim holds only under its generative model; the attack nulls probe
departures from it but do not remove the dependence. The full validation set out in the
specification was not done: no K-12 benchmark dataset, no FDRBench entrapment and no co-abundance
AUROC (readiness R10, deferred), and the sealed EDL933 holdout was not used, so no holdout claim is
made. Real-data evidence is one pilot sample and the campaign, whose answer keys test biology the
authors expected, not error rates. The pilot numbers other than detection were measured on
1.0.0-era code. The campaign evidence is not yet merged into the main branch of the repository.
[PLACEHOLDER: update once the campaign PRs merge.]

**What GIPPA does not estimate.** Identification error is the search engine's; GIPPA inherits it
at q ≤ 0.01. The `two_peptide` and `protein_list` tiers have no gene FDR; in simulation the
two-peptide rule's gene-level FDP was 0.054 at q = 0.01. Variant calls have no family-wise or
false-discovery control, and the variant gate was not met. Abundance has no test or interval. ORA
and the determinants permutation test are shipped flagged as uncalibrated.

**Known weak spots.** The placement FDR cannot accept tolerant placements in small strata, so
cross-strain variants are not recovered at defaults (Section 4.4). Novel-ORF classes in typical
bacterial genomes are too small to reach q ≤ 0.01, so the module will usually call nothing;
multi-peptide false ORFs make the class FDR anti-conservative. TU-block autocorrelation statistics
are anti-conservative with predicted TUs; the island FWER and operon discordance claims hold only
within stated regimes, and no output tells the user whether the data lie inside them. The genome
decoys under-count chance hits of protein-like wrong peptides. The RplL/RplJ stoichiometry check
fails in about half of the genomes tested, for reasons not established. The measured runtime is
several times the specification's target.

**What version 1 does not do.** One sample per run, with no tests across samples; no spectrum
searching; no start-codon correction, signal-peptide prediction or PTM naming; no recoding
catalogue (selenocysteine, programmed frameshifts) beyond reporting where peptides fall; and
deferred features listed in the changelog (assembly checks from frameshift-spanning peptides,
proteome mass fractions, RNA-seq or Ribo-seq inputs, homology evidence for novel ORFs).

<!-- src: methods.md §21 (Scope of the evidence); docs/preprint/outline.md §6; results.md §1.2
(two_peptide FDP 0.054; variants), §4; changelog.d/+known-limitations.misc.md (via outline §6);
P2S §2.6, §2.9, §6 -->

## 8. Discussion

GIPPA makes one claim about itself that we can test: that every number it reports either carries an
error estimate or says it has none. The calibration shows where that holds and where it does not.
Placement, gene-level detection, region tests and GSEA hold their nominal levels within the
tolerance in simulation; autocorrelation, islands and operon discordance hold them only in stated
regimes; ORA and the determinants test do not, and say so on every row. A user can read the flags
and the methods notes and know which column to trust.

The campaign tested something calibration cannot: whether the defaults behave sensibly on real
data from organisms the developers did not tune for. They did in the parts that matter most for
the gene table. Mapping was exact for every accepted peptide in 17 genomes, the six-frame route lost
few genes, and gene-level FDR estimates sat at or below 0.01. The campaign also found the most
consequential defect in the project's history, the decoy-pairing error in the picked gene FDR,
which no simulation suite had covered because the simulator paired decoys through locus tags. The
fix and its calibration attack shipped in 1.1.1. That the defect was found by running real
cross-strain searches, not by reasoning about the estimator, argues for keeping such campaigns part
of release testing.

The campaign's negative results are informative too. Two error-control rules that are correct in
isolation, a +1 correction in the placement FDR and class-wise decoy competition for novel ORFs,
make GIPPA nearly powerless for the two discovery tasks that motivate proteogenomics, strain
variants and new genes, at the sizes typical of bacterial genomes. GIPPA now warns when a stratum
or class cannot reach α, and we report the post hoc positive-control tier only as a guide to where
true signal sits. A principled fix, such as pooling classes or strata with a calibrated mixture,
would change a statistical rule and needs its own calibration and decision record. We think it is
better to ship a method that says "no call possible here" than one that calls at an error rate it
cannot estimate.

[PLACEHOLDER: Khoa: comparison with related tools and positioning; needs citations, e.g.
[CITE: bacterial proteogenomics pipelines], [CITE: DIA-NN-based microbial proteomics workflows],
[CITE: entrapment-based FDR assessment, beyond Wen et al. 2025].]

<!-- src: Sections 3-5 above; docs/adr/ADR-0078 (Context; calibration attack added); P2S §6 (Novel
ORFs: any change needs an ADR and a Khoa decision) -->

## Data and code availability

- **Proteomics data.** Pilot: PXD078459 (iProX). Campaign: PXD058422 (*E. coli* K-12), PXD060434
  (*E. coli* NCM3722), PXD065022, PXD069004, PXD055395, PXD072227, PXD049476, PXD060688, PXD065346,
  PXD057335, PXD066008 and PXD079213 (Table 6). Cited by accession, not redistributed.
- **Genomes.** Nissle 1917 GCF_021559835.1 (pilot); K-12 MG1655 U00096.3; the campaign genomes are
  listed with their NCBI accessions and checksums in `validation/campaign/genomes.tsv`. Proteome:
  UniProt UP000000625 (pilot arm A1).
- **Not redistributed.** RegulonDB (its terms), DIA-NN output rows and peptide sequences (aggregates
  only are committed), DIA-NN itself and the Cascadia checkpoint.
- **Code and evidence.** https://github.com/khoa-yelo/gippa (MIT): source, `GIPPA_SPEC.md`,
  `docs/methods/`, `docs/adr/`, `tests/calibration/`, `validation/pilot/`,
  `validation/release/1.1.1/`, and (once merged) `validation/campaign/`. The figure scripts of this
  paper are in `docs/preprint/figures/`. Calibration harness records are kept on the development
  cluster; [PLACEHOLDER: archive the CAL111 records with the Zenodo release].
- **Publication of the pilot figures.** [PLACEHOLDER: Khoa, Q-K11: which pilot figures may appear.]

<!-- src: methods.md §24; handover plan Phase 2 table (accessions); CAMP/genomes.tsv (worktree
campaign-p2a); docs/preprint/outline.md §5 -->

## Funding, competing interests, acknowledgements

[PLACEHOLDER: Khoa: funding; competing interests; acknowledgements (Stanford FarmShare).]

## Figure legends

**Figure 1. The GIPPA workflow.** Inputs are a genome, an annotation and peptide or protein
evidence from any search engine (`gippa convert` reads DIA-NN, Cascadia and generic tables);
`gippa makedb` optionally writes a strain database, with or without six-frame ORFs, for the user's
own search. Annotation layers are computed by third-party tools run as separate programs. Peptides
are placed on the six-frame genome with exact and tolerant matching; genome decoys give a placement
FDR. The gene table holds one row per Bakta locus tag, with a detection tier, gene-level FDR,
abundance and detectability. Context, operon, variant, novel-ORF, function and determinant modules
build on it. Red boxes mark where GIPPA reports an error estimate; methods without a valid estimate
are flagged. Draft schematic, drawn programmatically.

**Figure 2. Simulation calibration of GIPPA 1.1.1.** Full calibration tier on commit `d2bd17e4`
(Slurm array 2022338). (a) Quantile–quantile plots of null p-values (−log10 scale) for four streams
where the claim applies, and the distance-model operon fallback, a known conservative failure; grey,
pointwise 95 % band for 1,000 uniform p-values. (b) Every harness record at nominal level 0.01,
0.05 or 0.10 (level, FDR, FWER criteria): observed rate against its bound, both divided by the
nominal level; points above the diagonal exceed the bound. Colours: passing records; unconfirmed
(above the bound in one block of two); failing records in arms that are recorded only or whose rows
ship flagged; known failures (strict expected failures, each tied to an ADR or open decision).
Records with zero observed error are drawn at 0.02. (c) Observed rejection rate or mean FDP against
the nominal level for representative streams; dashes mark each stream's bound (nominal + 2 Monte
Carlo SE). Points are offset horizontally for legibility. All results are conditional on the
gippasim simulator.

**Figure 3. The *E. coli* Nissle 1917 pilot.** Three wild-type DIA runs (PXD078459), searched against
the UniProt K-12 proteome (A1, orange) or the Nissle strain database (A2, blue). Values from the
1.1.1 detection code on the frozen pilot inputs. (a) Share of genes detected (gene q ≤ 0.01) by
subset, with the number of genes in each subset; Wilson 95 % intervals treat genes as independent
and are optimistic. A1 cannot detect strain-specific genes by construction. (b) Pre-registered
hypothesis H1: share of discordant genes detected only in A2; share of A2-only genes that are
strain-specific; share that lie in plasmids, prophages or BGCs (vermillion), against the
pre-registered "mostly" requirement (> 0.5, dashed) and the base rate among all genes (grey).
(c) Free entrapment for A1: realised lower bound, combined estimate (dot) and one-sided 95 % upper
bound of the FDP at precursor level (the level DIA-NN's q controls) and peptide level, against the
1 % target. The upper bound assumes independent false matches.

**Figure 4. The campaign on public DIA data.** GIPPA 1.1.1 code on eleven datasets (17 genomes;
*E. coli* K-12 Phase 1 in orange). (a) Route-A genes detected at gene q ≤ 0.01 as a share of coding
genes, with counts. (b) Share of route-A detected genes that route B (database plus six-frame ORFs)
also detects; dashed, the pre-registered prior of 0.95. Intervals are Wilson 95 %, optimistic
because genes co-vary. (c) TE4: log2(RplL/RplJ) on iBAQ (filled, the verdict) and Top3 (open),
against the band [1, 3] (grey); vermillion, fail. The archaea are not testable. (d) Outcomes of the
pre-registered answer-key checks per dataset: pass, fail, not testable or not run, and reported
without a threshold. *E. coli* includes one dropped check counted as not testable.

**Figure 5. Novel ORFs in the campaign.** (a) For each run and candidate class, the smallest class
FDR against the number of candidates in the class (log scales); dashed, α = 0.01; dotted, the post
hoc tier τ_prfB = 0.167; solid, 1/n for reference (the class q cannot fall below 1/T, with T the
class's targets). No class reaches α. (b) Candidates at or below τ_prfB per run: positive controls
(vermillion), non-control candidates that overlap a RefSeq gene on the same strand and replicon
(blue), and those that do not (grey); labels give passers of all candidates. GIPPA called none of
them. The tier was chosen after the runs from the same data and is not an error-controlled call.

## References

- Aho AV, Corasick MJ (1975). Efficient string matching: an aid to bibliographic search. *Commun ACM* 18:333–340.
- Barber RF, Candès EJ (2015). Controlling the false discovery rate via knockoffs. *Ann Stat* 43:2055–2085.
- Benjamini Y, Hochberg Y (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. *J R Stat Soc B* 57:289–300.
- Blin K, et al. (2025). antiSMASH 8.0: extended gene cluster detection capabilities and analyses of chemistry, enzymology, and regulation. *Nucleic Acids Res* 53(W1):W32–W38.
- Camacho C, et al. (2009). BLAST+: architecture and applications. *BMC Bioinformatics* 10:421.
- Camargo AP, et al. (2024). Identification of mobile genetic elements with geNomad. *Nat Biotechnol* 42:1303–1312.
- Chevez-Guardado R, Peña-Castillo L (2021). Promotech: a general tool for bacterial promoter recognition. *Genome Biol* 22:318.
- Demichev V, et al. (2020). DIA-NN: neural networks and interference correction enable deep proteome coverage in high throughput. *Nat Methods* 17:41–44.
- Di Tommaso P, et al. (2017). Nextflow enables reproducible computational workflows. *Nat Biotechnol* 35:316–319.
- Frankenfield AM, et al. (2022). Protein contaminants matter: building universal protein contaminant libraries for DDA and DIA proteomics. *J Proteome Res* 21:2104–2113.
- He K, et al. (2015). A theoretical foundation of the target-decoy search strategy for false discovery rate control in proteomics. arXiv:1501.00537.
- Kingsford CL, Ayanbule K, Salzberg SL (2007). Rapid, accurate, computational discovery of Rho-independent transcription terminators illuminates their relationship to DNA uptake. *Genome Biol* 8:R22.
- Korotkevich G, et al. (2021). Fast gene set enrichment analysis. *bioRxiv* 060012.
- LaFleur TL, Hossain A, Salis HM (2022). Automated model-predictive design of synthetic promoters to control transcriptional profiles in bacteria. *Nat Commun* 13:5159.
- Li H (2023). Protein-to-genome alignment with miniprot. *Bioinformatics* 39:btad014.
- Reiter L, et al. (2009). Protein identification false discovery rates for very large proteomics data sets generated by tandem mass spectrometry. *Mol Cell Proteomics* 8:2405–2417.
- Roots CT, Lukasiewicz AJ, Barrick JE (2021). OSTIR: open source translation initiation rate prediction. *J Open Source Softw* 6:3362.
- Salgado H, et al. (2024). RegulonDB v12.0: a comprehensive resource of transcriptional regulation in *E. coli* K-12. *Nucleic Acids Res* 52:D255–D264.
- Sanders J, et al. (2025). A transformer model for de novo sequencing of data-independent acquisition mass spectrometry data. *Nat Methods* 22:1447–1453.
- Savitski MM, et al. (2015). A scalable approach for protein false discovery rate estimation in large proteomic data sets. *Mol Cell Proteomics* 14:2394–2404.
- Schwanhäusser B, et al. (2011). Global quantification of mammalian gene expression control. *Nature* 473:337–342.
- Schwengers O, et al. (2021). Bakta: rapid and standardized annotation of bacterial genomes via alignment-free sequence identification. *Microb Genom* 7:000685.
- Su H, Zhang R, Söding J (2024). UniOP: a universal operon prediction for high-throughput prokaryotic (meta-)genomic data using intergenic distance. *bioRxiv* doi:10.1101/2024.11.11.623000.
- Tjaden B (2026). TerminatorNet: comprehensive identification of intrinsic transcription terminators in bacteria. *Bioinformatics* 42:btag116.
- Wen B, et al. (2025). Assessment of false discovery rate control in tandem mass spectrometry analysis using entrapment. *Nat Methods* 22:1454–1463.
- Dataset papers of the campaign (PMC identifiers from the campaign plan): PXD065022 PMC12687308; PXD055395 PMC13253324; PXD069004 PMC13046382; PXD049476 PMC11066055; PXD065346 PMC12850539; PXD072227 PMC13753778; PXD057335 PMC12913811; PXD060688 PMC12288472; PXD066008 PMC12949519; PXD079213 PMC13736083. [CITE: full references to be completed from these PMC records; Phase 1 PXD058422 and Phase 1b PXD060434 papers.]

<!-- src: docs/preprint/methods.md › References (methods), verified there; handover plan Phase 2
table (PMC identifiers) -->
