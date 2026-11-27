# GIPPA: genome-integrated analysis of bacterial proteomes with an error estimate for every claim

<!-- Draft manuscript, round 2 (docs agent). Every number is copied from the file named in the
`src` comment after its paragraph or table; paths are relative to the repository root unless they
start with $GIPPA_ROOT. Abbreviations used in the comments:
CAL111 = $GIPPA_ROOT/logs/calibration/full-d2bd17e42dcf-20261125-073329 (calibration_full on
d2bd17e4, Slurm array 2022338, the run readiness 1.1.1 reuses);
CAL110 = $GIPPA_ROOT/logs/calibration/full-4a6fa998300c-20261109-160014 (the 1.1.0 full tier);
RD = validation/release/1.1.1/readiness.json and READINESS.md;
PIL = validation/pilot/metrics/metrics.json; PIL111 = validation/pilot/recheck_1_1_1/metrics_1_1_1.json;
CAMP = validation/campaign/ (not yet on main: Phase 1 in worktree campaign-ecoli, Phase 2 summary,
group A and novel_pc/ in campaign-p2a, group B in campaign-p2b, group C in campaign-p2c);
P2S = CAMP/PHASE2_SUMMARY.md; NPC = CAMP/NOVEL_ORF_PC_THRESHOLD.md and CAMP/novel_pc/;
SUPP = docs/preprint/supplementary_tables.md (Tables S1-S4).
Records of CAL111 equal those of CAL110 except the ADR-0078 decoy-loss attack streams (4 records,
lower FDP) and the new decoy-pairing attack (3 records); checked record by record, 1,103 shared keys.
Aggregates only: no peptide sequence, no DIA-NN row; gene rows only for named marker genes. -->

**Authors.** [PLACEHOLDER: Khoa: author list, affiliations, ORCID iDs, corresponding author]

**Venue.** [PLACEHOLDER: Khoa, Q-K11. Draft written for a bioRxiv preprint, then a journal to be chosen]

## Abstract

**Motivation.** Bacterial proteomes are usually searched against a reference proteome, which
misses strain-specific genes and loses genome context. A strain database fixes the first problem
but adds two errors that search engines do not estimate: placing peptides on the genome and calling
genes detected.

**Results.** GIPPA takes one assembled genome and any search engine's output. It estimates a
placement FDR with genome decoys and a picked gene-level FDR, builds one gene table, and runs
dependence-aware context tests. Each value carries a q-value, FDR, FWER or interval, or a flag
saying none exists. In simulation, 1,208 of 1,244 calibration records passed. Five were
unconfirmed and 31 failed, all known before release or recorded only, including a contaminant
attack on placement, variant calls below target F1 and autocorrelation with predicted operons. In an *E. coli* Nissle 1917 pilot, a strain database detected 3,051 of 4,685
genes at gene q ≤ 0.01 against 2,645 with a K-12 reference, mainly through 316 genes without a
K-12 counterpart. The pre-registered hypothesis that the gain lies mostly in mobile elements was
not supported. In eleven public DIA datasets (17 genomes) the six-frame route kept at least
96.8 % of detected genes, but no strain variant or novel ORF was called at default error levels:
small decoy strata and candidate classes cannot reach q ≤ 0.01. Release 1.1.1 waives three
acceptance criteria (real-data placement FDR, one ribosomal stoichiometry check, runtime) and
defers one.

**Availability.** MIT licence; Python package `gippa` 1.1.1; https://github.com/khoa-yelo/gippa.

<!-- src: abstract 241 words (checked with wc, headings excluded); 1,208 of 1,244, 5, 31 (14 + 17):
fig2_calibration.py stdout over CAL111; contaminant attack, variants F1, autocorrelation: results.md
Table R3b (identical in CAL111); 3,051 / 2,645 of 4,685 and 316 of 434: PIL111 › m3, h1; H1
not_supported: PIL111 › h1.verdict; eleven datasets, 17 genomes, 96.8 %: P2S §2.3-2.4 and
CAMP/ecoli/REPORT.md row 14; 0 variants, 0 novel ORFs: P2S §2.6, §2.9; CAMP/ecoli/REPORT_1b.md;
waivers and deferral: RD rows R2, R5, R7, R10 -->

## 1. Introduction

A bacterial proteome is shaped by the genome that encodes it. Genes in one operon are transcribed
together; prophages and plasmids carry genes that may be silent in a given condition; codon usage,
promoters and terminators affect how much of a protein is made. Proteomics measures the outcome,
but the usual analysis loses the genome. Spectra are searched against a reference proteome, often
of a type strain, and the result is a list of protein groups. Genes that the reference lacks
cannot be found, peptides that differ from the reference are lost or misassigned, and the
coordinates needed to ask about genomic context are gone [CITE: reference-database mismatch in
bacterial proteomics].

Proteogenomic workflows restore the genome by searching six-frame translations or strain
databases [CITE: bacterial proteogenomics workflows and six-frame searches], and several tools map
identified peptides back to genome coordinates [CITE: peptide-to-genome mapping and proteogenomic
annotation tools]. To our knowledge, most stop at the identification step: they report the search
engine's FDR, but not the error of placing a peptide on the genome or of calling a gene detected,
and they test genome-context questions, if at all, as if genes were independent [CITE: survey of
proteogenomics tools and their error control]. Both gaps matter. A peptide can be identified
correctly and still be placed at the wrong locus, and a gene can be called detected on peptides
that match a decoy as well. Neighbouring genes share transcription units (TUs), and abundance is
autocorrelated along the replicon, so a test that treats genes as independent rejects too often
[CITE: spatial autocorrelation of gene expression along bacterial chromosomes].

GIPPA adds the two error estimates that a strain-database search lacks, and says where no
estimate exists. It is built for a single sample on any bacterial or archaeal genome. It does not
search spectra: it takes peptides or proteins from any engine and inherits the engine's
identification error, which it reports but does not re-estimate. Its contributions are:

1. a placement FDR, estimated with decoys built from the genome itself, for every peptide placed
   on the genome;
2. a picked gene-level FDR for strain databases, with decoys that compete gene by gene;
3. context tests that respect dependence along the replicon: circular-shift region tests,
   TU-block autocorrelation, an abundance-island scan with family-wise error control, and operon
   tests;
4. calibration shipped with the software: every statistical method has a simulation suite with
   fixed tolerances and attack nulls built to make it fail, and a method that fails is shipped
   flagged or with a narrowed claim, with the failure disclosed.

What this buys is shown by the pilot (Section 4.3): searching the strain's own database instead
of a K-12 reference detected 406 more genes of *E. coli* Nissle 1917 at a gene-level q of 0.01,
most of them genes K-12 does not have. We describe GIPPA 1.1.1 and the evidence behind it, namely
simulation calibration, eight true-example assertions on real genomes, a pre-registered pilot and a
campaign on public data-independent acquisition (DIA) datasets from eleven organisms. We report
what passed, what failed, what was waived and what was deferred. The Supplementary Methods
(`docs/preprint/methods.md`), Supplementary Results (`docs/preprint/results.md`) and Supplementary
Tables (`docs/preprint/supplementary_tables.md`) give every formula, default and record.

<!-- src: methods.md §1 (GIPPA does not search spectra; inherits identification error); GIPPA_SPEC.md
§0.2 R3-R4; 406 = 3,051 − 2,645 (PIL111 › m3.A2/A1.subsets.all.k); positioning claims are
placeholders until citations are added -->

## 2. Design and implementation

### 2.1 Inputs, routes and arms

The required input is a genome FASTA file. A Bakta JSON annotation [Schwengers et al. 2021] is
used as is; a GFF3 or GenBank annotation from another tool is passed to Bakta, which keeps its CDS
and adds the rest. Evidence can be a peptide table, a protein FASTA file with an abundance table,
or both. `gippa convert` writes the peptide table from DIA-NN 2.x reports [Demichev et al. 2020],
Cascadia 0.0.7 de novo output [Sanders et al. 2025] or a generic table. Only the peptide column is
required; q-values, PEPs, decoy flags, intensities and run labels are used when present. Rows above
the input q cut-off (0.01) are kept for FDR estimation but never count as evidence.

For users who search their own spectra, `gippa makedb` writes a strain database from the Bakta
annotation. It contains the CDS and sORF proteins, optionally six-frame open reading frames (ORFs)
of at least 20 codons with class-tagged identifiers, optionally proteins patched by an earlier
run, and a 381-entry contaminant library [Frankenfield et al. 2022]. It adds no decoys; decoys
come from the search engine. We call a search against the Bakta proteins route A, and a search
that adds the six-frame ORFs route B. Each search configuration of the same spectra is an arm:
the pilot (Section 4.3) has a K-12-reference arm, a strain-database arm and a de novo arm. The
command-line options are listed in Table S1.

<!-- src: methods.md §2 (Peptides, Converters), §4; gippa makedb --help (--six-frame-orfs, --min-orf-aa
default 20, --patch, --no-contaminants); src/gippa/data/db_sources.yaml (contaminants, 381 entries,
methods.md §5.5); SUPP Table S1 -->

### 2.2 Annotation layers

Bakta 1.12.2 (database 6.0) supplies genes, CRISPR arrays, resistance genes and the origin of
replication. geNOMAD 1.12.0 [Camargo et al. 2024] calls prophages and plasmids, antiSMASH 8.0.4
[Blin et al. 2025] calls biosynthetic gene clusters (BGCs), and UniOP [Su et al. 2024] gives operon
pair probabilities, with an intergenic-distance model as fallback. Transcription units are maximal
runs of same-strand genes linked with probability ≥ 0.5. Promotech [Chevez-Guardado and
Peña-Castillo 2021], Promoter Calculator [LaFleur et al. 2022], TerminatorNet [Tjaden 2026],
TransTermHP [Kingsford et al. 2007] and OSTIR [Roots et al. 2021] add promoter, terminator and
translation-initiation layers, and codon-usage indices are computed in process. These layers
predict σ70 promoters and intrinsic terminators only, so a missing call is never read as absence,
and a missing layer is flagged rather than reported as "no". The operon and terminator cut-offs
were fixed on *E. coli* K-12 RegulonDB 14.5.0 data [Salgado et al. 2024].

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
library or written to `unmapped.tsv` with a reason.

<!-- src: methods.md §2 (genetic codes 11, 4, 25), §5.1-5.5 -->

### 2.4 The gene table

The gene table has one row per annotated CDS, sORF or pseudogene, keyed by its Bakta locus tag. A
peptide counts as evidence for a gene when its input q and placement q are both ≤ 0.01 (or it is
an exact hit from a same-strain search), its locus class assigns it to that gene, and it is unique
to the gene or shared only with identical proteins. Fig. 6 shows one locus: the peptides of
*E. coli* prfB, including four that lie upstream of the annotated start, count for the gene. The
detection tier is chosen once per run from what the input carries: `decoy_picked` (decoys with
protein identifiers that pair to genes), `mayu` (targets with q or PEP but unusable decoys),
`two_peptide` (no q or PEP, as for de novo input) or `protein_list` (proteins only). The first two
estimate a gene-level FDR (Section 3.2); the last two do not, and say so. Abundance is iBAQ
[Schwanhäusser et al. 2011] by default (Top3 as an alternative) on unique peptides; GIPPA does not
normalise across samples. A cross-fitted detectability model predicts each gene's probability of
detection from its observable peptides, length, hydropathy, pI and transmembrane segments, so that
context tests can compare observed with expected detection.

<!-- src: methods.md §6 (evidence rule, tiers), §7 (abundance), §8 (detectability);
$GIPPA_ROOT/campaign/ecoli_k12/runs/out_e2_fix/tables/genes.tsv (KOIAIB_03083 prfB:
n_peptides_unique 17, flags nterm_extension_evidence) and mapping/peptide_hits.tsv (13 cds_inframe
and 4 nterm_extension placements on KOIAIB_03083); fig6_prfb_locus.tsv -->

### 2.5 Context, operon, variant, novel-ORF and function modules

Region tests ask whether the genes of a prophage, BGC, CRISPR or resistance region are detected
more or less often than detectability predicts, or are more or less abundant, than genes at
random positions on the same replicon. Autocorrelation statistics (TU-block Moran's I and join
counts along the replicon), an abundance-island scan with family-wise error control, and an
origin-to-terminus gradient describe abundance along each replicon. The operon module tests
discordance of abundance within predicted operons and looks for missed operons among non-operon
pairs. The variant module calls amino-acid differences from tolerant placements and flags
mass-equivalent calls that could be modifications. The novel-ORF module groups peptides outside
annotated genes into candidate ORFs by class (N-terminal extension, intergenic, out of frame,
antisense) and estimates a class-wise FDR (Section 3.4). Function enrichment runs
over-representation analysis (ORA) for detection and fgsea [Korotkevich et al. 2021] for
abundance, over GO, KEGG and COG collections, with TUs as the unit. A determinants module models
detection and abundance from sequence features, and `gippa compare` reports agreement between runs
on the same genome.

<!-- src: methods.md §9-17 -->

### 2.6 Interface, outputs and reproducibility

GIPPA is a Python package (Python ≥ 3.11) with a command-line interface of nine commands (`setup`,
`preflight`, `run`, `makedb`, `convert`, `compare`, `report`, `example fetch` and `demo`; every
option is in Table S1) and a bundled Nextflow workflow [Di Tommaso et al. 2017] with profiles for
local, Slurm, conda, Apptainer and Docker execution. Each run writes 22 documented result tables
(18 under `tables/`, 4 under `mapping/`), each a TSV file with a JSON schema; `gippa compare` adds
three, and nine plot-support tables back the figures (Table 1). A self-contained HTML report holds
15 plots in a colour-blind-safe palette, each with its data table. The run manifest records input
checksums, tool and database versions, container digests, the resolved configuration and the seed.
All randomness flows from one master seed (20261014 by default) through per-step generators, so
the same inputs and seed give byte-identical tables (TE7, Section 4.2).

<!-- src: gippa --help and gippa <cmd> --help (captured 2026-11-27 at the 1.1.1 tree of this branch;
SUPP Table S1); docs/reference/tables.md (34 documented tables: 18 tables/, 4 mapping/, 3 compare,
9 plots/support); methods.md §19-20; src/gippa/config/defaults.yaml:2 (seed) -->

### 2.7 How GIPPA was built

GIPPA was written against a fixed specification, with decision records that override it once the
project owner accepts them; a change to any default, threshold, tolerance or gate needs such a
record. Most of the code, tests, simulator and documentation, including this draft, were written by
AI coding agents working in separate roles (specification, tests, implementation, review,
statistics). Every commit names its role, and continuous integration checks that no commit
changes both code and its tests and that protected files change only with an accepted record. The
validation evidence is therefore not independent of the authors of the code: the simulator, the
tests, the answer keys and the reports came from the same project. [PLACEHOLDER: Khoa:
AI-assistance statement, wording and scope, and a statement on validation independence.]

<!-- src: methods.md §25; AGENTS.md; scripts/ci/check_roles.py; scripts/ci/check_protected.py -->

## 3. Error control

GIPPA reports every p- and q-value, passing or not. Each claim belongs to one family with a stated
error rate (Table 2), and no error rate is controlled across families. Where no valid estimate
exists, the row carries a flag (`ora_uncalibrated`, `null_uncalibrated`, `ci_unestimated`,
`class_error_unestimated` and others) instead of a number that looks like one. The calibration
status in Table 2 is the simulation evidence of Section 4.1, where every failure is listed.

### 3.1 Placement FDR

A peptide can be a correct identification and still be placed at the wrong locus, or be a wrong
identification that happens to match the genome. GIPPA estimates the share of accepted placements
that are chance hits of a wrong peptide with genome decoys: every translated frame reversed for
database-search input, and shuffled within stop-free segments for de novo input, because de novo
errors often reverse a stretch of sequence. Peptides are ranked by a chance-hit score, lists are
formed per decoy type and number of differences m ∈ {0, 1, 2}, and within a list
FDR(e) = (W(e) + c)/T(e). Here W sums decoy hits weighted by an estimate of the probability that
each is wrong, and c is a +1-type correction [Barber and Candès 2015; He et al. 2015]. With c = 1
a list's q-values cannot fall below 1/T, so a stratum with fewer than about 100 targets cannot
accept anything at 0.01; since 1.1.1 the run warns when this happens. Exact hits from a
same-strain search receive q = 0 and are flagged `same_db`; for such input GIPPA adds no error
estimate of its own. The genome decoys under-count chance hits of wrong peptides that resemble
real proteins.

<!-- src: methods.md §5.4; docs/adr/ADR-0010, ADR-0019, ADR-0079; CHANGELOG.md 1.1.1 Features;
docs/methods/placement-fdr.md › Limitations -->

### 3.2 Gene-level FDR

In the `decoy_picked` tier GIPPA adapts the picked protein FDR [Savitski et al. 2015] to genes. In
each gene, the best accepted target peptide competes with the best decoy peptide paired to it, and
FDR(t) = (D(t) + 1)/T(t) over target-won genes T and decoy-won units D. Every decoy that passes the
input q cut-off competes: a decoy that pairs to one gene competes in that gene, and a decoy that
maps to no gene or to several competes in its own unit, which can only add to D. The estimator has
no finite-sample proof; its calibration rests on simulation, including an attack in which genes
are paired only through the votes of accepted peptides (Section 4.1). The MAYU tier [Reiter et al.
2009] uses a closed-form expectation of false genes from the peptide-level error, with no +1 term
and no finite-sample guarantee. Release 1.1.1 introduced the own-unit rule after the campaign
found that the earlier rule could be anti-conservative (Section 4.4).

<!-- src: methods.md §6 (Picked gene FDR, Changes in 1.1.1, MAYU); docs/adr/ADR-0078 (Decision,
Calibration attack); docs/methods/detection-tiers.md "Decoy pairing" -->

### 3.3 Context statistics under dependence

Gene-level tests along a genome are dependent: genes share TUs and abundance is autocorrelated.
GIPPA's region tests therefore use circular shifts. The region's genes are moved together around
the replicon, which keeps the region's size, its internal structure and the autocorrelation of the
background, and the studentised observed-minus-expected statistic is recomputed at every shift. The
p-values of one region type and statistic form a family controlled with the Benjamini–Hochberg (BH)
procedure [Benjamini and Hochberg 1995]. On a small replicon the shift null has few distinct values,
so p-values are coarse, and a region that covers its whole replicon cannot be shifted at all, so its
p-value is 1 by construction. Confidence intervals for region effects use a TU-block bootstrap.
Abundance islands are scanned over windows of TUs with a block bootstrap, and the maximum statistic
over all windows of all replicons gives a family-wise error rate (FWER). Its validity depends on
the block length matching the dependence. The claim is restricted to a stated regime, and the run
reports the block length and lag-1 Moran's I as diagnostics, but no rule turns them into an
inside-or-outside verdict.

<!-- src: methods.md §9-10; docs/methods/region-tests.md; docs/methods/islands.md; P2S §3 (whole-
replicon region: circular shift cannot vary it); results.md §1.2 (islands domain) -->

### 3.4 Novel ORFs, operons, function and determinants

The novel-ORF class FDR is estimated within each class and never pooled. In the database regime,
used whenever the search included `makedb` six-frame entries, each candidate competes with the
engine decoy of its own tagged entry. The class FDR keeps a 1/f rescaling term from an earlier
decision record, which equals 1 in the pipeline because every tagged decoy is paired. In the de novo
regime, peptides of a class are ranked by the chance-hit score as in Section 3.1 and an ORF takes
the minimum q of its peptides. Both regimes have the same 1/T floor, so a class needs about 100
candidates to be callable at 0.01. Operon discordance compares the abundances of the two genes of
a predicted operon pair with a per-pair statistic whose noise model depends on abundance,
BH-adjusted at 0.10. Missed operons are non-operon pairs whose partners behave like operon
partners; their FDR is estimated with convergent and sandwich decoy pairs. ORA uses a TU-level
conditional Poisson-binomial null with detectability as covariate; its p-values ship flagged
`ora_uncalibrated` because the null fails under realistic TU errors. GSEA runs R fgsea 1.36.2 as a
separate process over TU-level ranks. Its claim is restricted to sets that are not clustered on
the genome, and it is conditional on detection. The `detection_confounded` flag is a screen, not a
guarantee: in the simulations behind it (run on a Python port of fgsea), unflagged rows still
rejected at up to about 1.3α, so that check is gated at 1.5α by a statistician decision, a looser
bound than the α + 2 SE of the tolerance file (Section 4.1). The determinants permutation p-value
fails the dependence attacks and ships flagged `null_uncalibrated`. A flag does not remove a row
from its family.

<!-- src: methods.md §13 (Class FDR: picked_decoy and denovo_subgroup regimes; "(D(t)/f + 1)/T(t) ...
ADR-0009's rescaling, which ADR-0078 retired from the gene FDR ... but not from the class FDR"; q
floor), §11, §14 (ORA), §15 (GSEA: 1.3α, 1.15α, 1.07α on the port; statistician decision gates the
unflagged-row check at 1.5α), §18; src/gippa/novel/fdr.py:300,308 and src/gippa/novel/class_call.py:73-90
(makedb_id_map pairs every tagged decoy, f = 1; statistician follow-up, round 2); campaign manifests
paired_share 1.0; docs/adr/ADR-0025, ADR-0026, ADR-0034, ADR-0035, ADR-0038 -->

**Table 2.** Error-control families and their status in GIPPA 1.1.1. Calibration status is from
simulation (Section 4.1); "within tolerance" means the estimate did not exceed nominal + 2 Monte
Carlo SE, not that it equals the nominal level.

| Claim | Family | Error rate | Level | Calibration status |
|---|---|---|---|---|
| Identification | the input search | search-engine q, inherited | 0.01 | not re-estimated; DIA-NN's precursor-level q is used as the peptide cut-off, so peptide-level FDR is not controlled |
| Peptide placement | per decoy type × difference stratum | weighted target–decoy FDR | 0.01 | within tolerance in the gated streams; the contaminant attack fails at q = 0.01; not measured on real data |
| Gene detection | eligible genes plus own units of unpaired decoys | picked FDR (D + 1)/T, or MAYU | 0.01 | within tolerance (no finite-sample proof); `two_peptide` and `protein_list` unestimated |
| Novel ORFs | per class and regime | class FDR from decoys | 0.01 | database regime within tolerance; de novo regime fails for multi-peptide and exact systematic false ORFs |
| Variants | none | placement q of supporting peptides | 0.01 | no family-wise or FDR control; F1 0.874 against the 0.90 target |
| Region tests | per region type × statistic | BH FDR over circular-shift p | 0.05 | within tolerance; prophage FDR 0.063 at q = 0.05 is within tolerance, not at level |
| Autocorrelation | 20 lags per replicon × statistic | BH FDR | 0.05 | at level with exact TUs only; anti-conservative with predicted TUs (flagged on every run) |
| Abundance islands | all windows, both directions | FWER (single-step maxT) | 0.05 | within a stated dependence regime only |
| Operon discordance | both-detected operon pairs | BH FDR | 0.10 | claim narrowed to smooth, monotone abundance noise |
| Missed operons | non-operon pairs | decoy FDR | 0.10 | vacuous with the shipped decoys (no calls) |
| ORA | per collection × direction | BH FDR | 0.05 | ships flagged `ora_uncalibrated` |
| GSEA | per collection | BH FDR | 0.05 | restricted claim; level checks on the R engine; unflagged rows under detection confounding pass only a 1.5α screen (Python port) |
| Determinants | one p per model | permutation p | – | ships flagged `null_uncalibrated` |

<!-- src: methods.md Table M2 (§18); results.md §1.2 and Table R3b (values identical in CAL111);
PIL › m8.peptide.assessment_reason ("q_value is DIA-NN's precursor-level Global.Q.Value"); RD row R2;
CAL111 novel_orfs.merged.jsonl strict-xfail streams; results.md line 199 (prophage 26 % above q) -->

## 4. Validation

### 4.1 Simulation calibration and attack nulls

Every statistical method was calibrated on data from gippasim, a simulator written within the
project. It plants genomes, TUs, prophages, abundance with replicon, region, island and AR(1)
effects, detection by a logistic model, and search and de novo errors at known rates. Every claim in
this section is therefore conditional on gippasim's generative model. The criteria are fixed in a
protected file. At nominal levels 0.01, 0.05 and 0.10, the observed rejection rate must not exceed
α + 2√(α(1−α)/n), the mean false discovery proportion (FDP) must not exceed q + 2·SD(FDP)/√R, the
family-wise rate must not exceed α + 2√(α(1−α)/R), and null p-values must pass a two-sided
Kolmogorov–Smirnov test at p > 0.001. A pass rules out excess error larger than two Monte Carlo
standard errors, not excess within that margin: for α = 0.05 the bound is about 1.5α at R = 300
and 1.3α at R = 1,000. A check that exceeds its bound in one block of R replicates is re-run on a
fresh block and fails only if both blocks exceed it; one excess alone is "unconfirmed".

<!-- src: tests/calibration/tolerances.yaml (levels, se_multiplier 2.0, ks_min_p 0.001);
methods.md §21 (criteria, two-block confirmation, ADR-0015); results.md §1.1 (1.5α at R 300,
1.3α at R 1,000) -->

The release calibration of 1.1.1 ran all 15 suites with 100 to 1,000 replicates each; no test
failed outside the expected failures (Table S4). Fig. 2b shows the 1,244 records with a level, FDR
or FWER criterion at the three nominal levels. Of these, 1,208 passed, 5 were unconfirmed, 14
failed in arms that are recorded only or whose rows ship flagged, and 17 were known failures
declared before the run. (The run holds 1,338 records in all; the rest are uniformity tests,
vacuous checks with too few targets, checks at other nominal levels such as the 1.5α screen, and
other criteria.) Null p-values were uniform where the claim applies (Fig. 2a). Observed error is
compared with its bound for representative streams in Fig. 2c and Table S2, which also lists the
streams that exceed their bounds. Version 1.1.1 changed only the gene-level decoy rule
(Section 3.2): the decoy-loss attack became more conservative (mean FDP 0.0072 against 0.0098 at
q = 0.01 with 20 % of decoys lost), and the new vote-route attack passed (0.0098, bound 0.0103).

<!-- src: RD rows R3, R4 (558 tests, 462/90/6; 284 attack tests 276/3/5; array 2022338); Fig. 2b
counts: fig2_calibration.py stdout {"pass": 1208, "unconfirmed": 5, "fail_flagged": 14, "fail_xfail":
17} over CAL111 lanes/*/*.jsonl + *.merged.jsonl, last block per (file, stream, criterion, level),
criteria level, fdr, fwer, fdr_counts at 0.01/0.05/0.10; 1,338 = all last-block records (uniform 29,
fdr_counts vacuous 26, pooled_p 2, slopes 2, citl 1, pooled_f1 2, the rest at other levels);
results.md §1.1 "about 1,375 harness records" counts the 1.1.0 run differently; decoy loss: CAL111
lanes/gene_fdr-group*/gene_fdr.jsonl attack_decoy_loss.L20 q 0.01 observed 0.00717 (CAL110 0.00983);
attack_decoy_pairing.reid100 q 0.01 observed 0.00977 bound 0.01033 -->

**Known failures and restricted claims.** We list them here once; Table 2 and Section 7 refer
back to this list. All were known before release and ship as expected failures, flags or narrowed
claims.

1. *Placement.* The contaminant attack exceeds the placement FDR bound at q = 0.01 (0.0150 against
   0.0102), and the placement gate of the specification was recorded as not met.
2. *Novel ORFs, de novo regime.* When a false ORF carries several wrong peptides, the class FDR
   is anti-conservative (0.0183, 0.0788 and 0.1139 at q = 0.01, 0.05 and 0.10), because an ORF's q
   is its best peptide's q. Exact systematic de novo errors fail more strongly: a wrong read that
   exactly matches a non-coding translation places at zero cost and is invisible to the decoys. In
   the clustered-error simulation the out-of-frame and antisense classes reached FDPs of 0.0749 and
   0.0719 at q = 0.01, about seven times nominal, and the intergenic class 0.0153; GC-poor
   intergenic reads gave 0.0126 (bound 0.0112). The pilot and campaign calls use the database
   regime, so these failures do not touch them.
3. *Variants.* Calls miss the F1 ≥ 0.90 target (0.874), and the isobaric-swap attacks give
   family-wise rates of 0.97 and 1.0; the variant gate was recorded as not met.
4. *Autocorrelation.* TU-block Moran's I holds level only with exact TUs. With UniOP-like TU errors it
   rejects at 0.071 (bound 0.0638), and under the distance-model fallback at 0.354, at α = 0.05;
   every such row ships flagged `tu_block_tu_error`.
5. *Operon discordance* fails under a step-shaped abundance-noise attack; its claim was narrowed to
   smooth, monotone noise after that attack failed.
6. *Islands.* The FWER is breached at AR(1) ρ = 0.9 with short block lengths (0.096 and 0.130 at
   1,000 and 1,800 genes, bound 0.0695); the claim was narrowed after the breach.
7. *Region tests.* The prophage FDR is 0.063 at q = 0.05 (bound 0.0654): within tolerance, 26 %
   above q, not at level.
8. *GSEA.* Most level claims were checked on the shipped R engine, but the flag checks ran only on
   the Python port. One arm of the detection-confounding attack (`null_b2`) rejected at 0.0833 at
   α = 0.05, above the α + 2 SE bound of 0.0752, and passes only the 1.5α screen (bound 0.1054).
9. *Missed operons.* The shipped decoys made no call in any of 300 power replicates, so the pass is
   vacuous.
10. *ORA and determinants* pass their gates only because their rows are flagged.

<!-- src: results.md §1.2-1.3 and Table R3b (values identical in CAL111: placement.merged.jsonl:43;
variants.merged.jsonl:1-3; autocorrelation 44, distance-model arm; discordance D-1b-step;
islands_recorded.jsonl:1-2; missed_operons.jsonl:35; gsea null_b2); CAL111 novel_orfs.merged.jsonl
strict-xfail records (multi_peptide_orfs.full 0.01831/0.07881/0.11392; g4k3.cds_outframe 0.07488;
g4k3.antisense 0.07188; g4k3.intergenic 0.01529; gc 0.01258 bound 0.01122); docs/methods/novel-orfs.md:265-275
(ADR-0050 mechanism); methods.md §15 (GSEA engines); PROGRESS.md gate records (M3, M7 not met) via
results.md §7 -->

### 4.2 True-example assertions on real genomes

Eight true-example assertions (TE1–TE8) check biology that any correct run must reproduce
(Table 3). They are evaluated on the Nissle pilot and on weekly *E. coli* K-12 and Nissle runs.
Two of them are consistency checks rather than tests of error control. TE1 passed with 28,129 of
28,129 strain-database peptides mapping exactly in frame, but every such peptide is a substring of
the searched database, so TE1 checks that database, genome and mapping agree. TE2 passed
vacuously: a same-database search cannot yield an amino-acid variant, and the weekly in-silico
proxy had no variant site either. TE4 failed: the ratio of the ribosomal stalk protein L7/L12
(RplL) to L10 (RplJ), expected near 4:1 (log2 in [1, 3]), was −0.48 on iBAQ, while AtpA/AtpD
passed. The cause was not established, and release 1.1.1 waives this check for that pair only. The
campaign later showed the same failure in 6 of 15 genomes (Section 4.4).

<!-- src: results.md Table R4 and §3; PIL › te1, te_real.te2, te_real.te4.pairs; RD row R5 (TE4 fail;
waived, ADR-0070 via ADR-0080); local/ci/validation/2026-11-25_te8_2020593.json,
2026-11-25_te2_2020593.json (1.1.1 weekly, d2bd17e4); P2S §2.5 -->

### 4.3 The Nissle 1917 pilot

The pilot was specified before it ran. Three wild-type *E. coli* Nissle 1917 DIA runs (PXD078459,
timsTOF HT) were searched with DIA-NN 2.3.1 in three arms. The K-12-reference arm searched the
UniProt K-12 proteome, a generic reference: Nissle belongs to phylogroup B2 and K-12 to A. The
strain-database arm searched the Nissle database built by `gippa makedb` from the wild-type
assembly (GCF_021559835.1). The de novo arm (Cascadia 0.0.7) was exploratory. GIPPA ran end to end
on all three. The pilot was measured on the 1.0.0-era code; the values that depend on gene
detection were recomputed with the 1.1.1 detection code on the frozen inputs, and a control with the
unchanged code reproduced the recorded tables byte for byte. The other pilot values were not
recomputed (Supplementary Results, Provenance).

<!-- src: methods.md §23; results.md §4 (arms, provenance); validation/pilot/recheck_1_1_1/README.md
(control 2001992 byte-identical; fix 2002006); results.md Provenance -->

**Detection.** Of 4,685 strain genes, the strain-database arm detected 3,051 (0.651 [0.638,
0.665], gene q ≤ 0.01) and the K-12-reference arm 2,645 (0.565 [0.550, 0.579]) (Fig. 3a). The 916
strain-specific genes share no tryptic peptide with K-12, so the reference arm cannot detect them;
the strain-database arm detected 316 of them (0.345 [0.315, 0.376]). Plasmid genes were detected at
0.769 [0.497, 0.918] with the strain database against 0.077 [0.014, 0.333] with the reference, and
prophage genes at 0.272 [0.204, 0.352] against 0.022 [0.008, 0.063]. These Wilson intervals treat
genes as independent, whereas genes of one TU co-vary, so they are optimistic. They describe this
sample, not a population. The strain-database arm's identification error is the search engine's,
inherited at q ≤ 0.01. GIPPA's placement FDR does not apply to it, because every peptide is an
exact hit of the searched database, so the release criterion on the placement FDR of real data was
not measured and is waived.

<!-- src: PIL111 › m3.A2/A1.subsets (all, strain_specific, plasmid 1/13 [0.0137, 0.3331], prophage
3/136 [0.0075, 0.0628]); results.md §4.2 and Table R6; RD row R2 (waived, ADR-0071 via ADR-0080) -->

**Pre-registered hypotheses.** H1 stated that the strain database detects more strain genes than
the reference, mostly plasmid, prophage and BGC genes. The first clause was expected by design:
434 genes were detected only with the strain database and 28 only with the reference, and 316 of
the 434 are strain-specific genes the reference could not detect. At TU level, 290 of 314 discordant
units favoured the strain database (0.924 [0.889, 0.948]; one-sided exact binomial test,
p = 1.8 × 10⁻⁵⁹). The clause that carried information failed: 106 of the 434 genes (0.244
[0.206, 0.287]) lie in plasmids, prophages or BGCs, against the pre-registered requirement of more
than half (Fig. 3b). H1 is therefore recorded as not supported. H2 (the reference arm's
genome-inconsistent FDP is at least its entrapment FDP) needs an entrapment benchmark and was
deferred. H3 (context conclusions agree between arms) could not be tested as stated: neither arm
called an island or an operon flag, and region tests were not compared. Its recorded "supported"
rests on substitute criteria: abundance Spearman ρ = 0.978 [0.973, 0.983] over 2,613 shared genes,
and a Jaccard index of 0.588 (error unestimated) for enriched categories, which include flagged ORA
rows. We do not count H3 as confirmed. H4 (de novo precision agrees with precision against the
strain database) was consistent with the lower bound in all 10 score bins; that check can refute
but not confirm.

<!-- src: PIL111 › h1 (b 434, c 28, b_tu 290, c_tu 24, mixed_tu 3, unit_effect 0.92357 [0.88879,
0.94810], p_tu 1.80e-59, mostly 106/434, strain_specific 316/434, verdict not_supported); validation/
pilot/metrics.py:553 (binomtest(k, n, 0.5, alternative="greater")); PIL111 › m4 (rho 0.978173
[0.972762, 0.982741], n_shared 2613), m5 (enriched_categories jaccard 0.587699, no interval);
results.md §4.4 and Table R7 -->

**Free entrapment.** K-12 peptides that are absent from the Nissle genome cannot be true if the
genome is complete and the sample is Nissle, so they serve as free entrapment for the reference
arm [Wen et al. 2025]. At precursor level, the unit that DIA-NN's q-value controls (28,347 targets,
157 entrapment hits), the combined FDP estimate was 0.0236, with a lower bound of 0.0055 and a
one-sided 95 % upper bound of 0.0265 (Fig. 3c). The point estimate exceeds the 1 % target, but the
lower bound does not, so the assessment is inconclusive. At peptide level, which DIA-NN's
precursor q does not control, the estimate was 0.0273. Of the 152 peptide-level entrapment hits,
136 are near-homologs of Nissle sequences, which push the estimate up; excluding them and short
peptides leaves 7 hits among 17,040 targets and an FDP of 0.0027 (upper bound 0.0049). The upper
bound is an estimated bound in expectation, not a bound on this sample, and it assumes independent,
equal-chance false matches. In simulations with clustered false matches it covered only 0.72–0.75
of the time. DIA-NN's use of protein information in its q-values biases the estimate the other
way, by an unquantified amount.

<!-- src: PIL › m8.precursor (N_T 28347, N_E 157, fdp_combined 0.023607, fdp_lower 0.005508,
ci.ucb95 0.026523, assessment inconclusive); m8.peptide (N_E 152, fdp_combined 0.027344, assessment
not_assessed, reason "peptide-level FDR is not controlled by the search"); m8.peptide.sensitivity.far
(breakdown N_E near_homolog 136, lt10 9, far 7; N_T_10 17040; fdp_combined_far 0.0027438,
ucb95_far 0.0049276); m8.caveats ("an estimated upper bound in expectation ... not a bound on this
sample"); results.md §4.5; CAL111 lanes/compare-group0/compare_attack_free_entrapment_recorded.jsonl:4,7 -->

**Novel ORFs, variants, de novo.** The strain-database arm called no novel ORF and no variant, as a
same-database search allows none. The reference arm called one novel ORF, which the pilot report
attributes to an annotation artefact. The de novo arm passed its feasibility rule on one run: 549
of 863 high-scoring de novo peptides (0.636 [0.604, 0.668]) matched the genome's six-frame
translation. Their shuffled-decoy hit share was 0.006 (5 of 863); all five also match the genome,
so this is a chance-hit rate for short peptides, not a false-match share. Variant class error rates
are unestimated for the de novo arm, because DIA data give no spectrum pairing.

<!-- src: results.md §4.6-4.7; validation/pilot/metrics/a3_feasibility.json (n_top 863, n_match 549,
share 0.63615 [0.60352, 0.66758], decoy_share 0.0057937 = 5/863; contaminant_share is a different
5-peptide set; validation/pilot/a3_feasibility.py:378-379); statistician follow-up, round 2 -->

### 4.4 A multi-organism campaign on public DIA data

**Design.** To see how GIPPA behaves beyond *E. coli*, we ran it on public DIA datasets with a
fixed recipe and answer keys written before each dataset's results. Phase 1 used *E. coli* K-12
MG1655 (PXD058422, four runs, Orbitrap Eclipse) with route A and route B, and a de novo arm. Phase
1b searched an NCM3722 dataset (PXD060434) against its own genome and against MG1655, to test
strain-variant calls. Phase 2 ran ten datasets from ten lineages, from a 1.14-Mb genome grown in
rabbit cells (*Treponema pallidum*, seven strains) to a 9.1-Mb genome (*Myxococcus xanthus*), and
two archaea (Table S3). In this paper "eleven datasets" means Phase 1 and Phase 2; Phase 1b is
counted separately. The datasets were chosen by the campaign plan to span lineages, genome sizes
and instruments from public whole-lysate DIA deposits, with reserves for datasets that failed the
first checks. It is a convenience sample, not a random one. Nine further accessions (reserves and
paired deposits) were inventoried but not run. Runs were selected within each deposit by written
rules (for example wild-type, exponential phase, DIA only). Every run used DIA-NN 2.3.1 with the
pilot's settings and the GIPPA 1.1.1 release code; the gapped *T. pallidum* assemblies used main
branch code with three assembly-gap fixes released with 1.1.1.

<!-- src: "GIPPA demo and discovery campaign agent handover plan.md" (goal, recipe, Phase 2 table,
Reserves); CAMP/datasets.tsv; CAMP/selection_rules.tsv (24 rule rows); CAMP/ directories PXD041262,
PXD065462, PXD065465, PXD066241, PXD069492, PXD073896, PXD074523, PXD075511, PXD080037 (files.tsv
only: inventoried, not run); CAMP/ecoli/REPORT.md; CAMP/ecoli/REPORT_1b.md; P2S header (Code:
bfc74ff6; 38e1cdba for T. pallidum) -->

**What the answer keys can and cannot show.** Each dataset's answer key was committed before its
report (Phase 1 on 19 November 2026, its report on 20 November; Phase 2 keys on 21 November, after
the Phase 1 results, including its TE4 failure, were known). The 158 Phase 2 checks have no
multiplicity adjustment; each is a separate pre-registered expectation. Comparisons of detection
depth across datasets are confounded by run count, gradient length and instrument. Two kinds of
result are consistency checks, not validation, and we report them as such. First, every accepted
route-A peptide in the 16 Phase 2 genomes (302,396 peptides; lowest Wilson lower bound 99.86 %) and
in *E. coli* (31,071) mapped exactly in frame at its own gene. That is expected, because route-A
peptides are substrings of a database built from the same genome. In the cross-strain arms, where
the genome differs from the searched database, 0.10–1.23 % of peptides found no placement. Second,
every route-A gene-level FDR estimate at q ≤ 0.01 was at most 0.0099 (Table S3), which follows from
the cut-off; it is the estimator's own value, not an external measurement. All runs used the
`decoy_picked` tier.

<!-- src: git log of CAMP/*/ANSWER_KEY.md (c16fa3ec 2026-11-19; 18a47b54 2026-11-19; 854c143c,
e3fb1514, 9070f5be 2026-11-21) and reports (7702782e 2026-11-20; 13f85cc6 2026-11-22); P2S §1
(158 checks), §2.1 (302,396; 99.860 % Grady 2,746), §2.2 (decoy_picked in 34 runs; ≤ 0.0099),
§2.9 (T. pallidum cross-strain unplaced 0.10-0.90 %); CAMP/ecoli/REPORT.md row 1 (31,071);
CAMP/ecoli/REPORT_1b.md 1b-4 (NCM3722 on MG1655 unplaced 305/24,842 = 1.23 %) -->

**Detection and the six-frame route.** Detection depth ranged from 0.455 of coding genes
(*Bacteroides thetaiotaomicron*, 11.5-min gradients) to 0.903 (*T. pallidum* Haiti B) (Fig. 4a).
In *T. pallidum*, 65.9–77.3 % of the distinct accepted peptides came from the host or the
contaminant library; we did not analyse how that affects the decoy model. Ribosomal proteins
(TE3) passed in every dataset; in the two archaea, where Bakta names almost no ribosomal protein,
the check used RefSeq names carried to Bakta genes by a shared stop codon, as pre-registered.
Searching six-frame ORFs (route B) increased the database 18- to 34-fold in records, and kept 96.8 %
to 99.9 % of route-A genes (point estimates; the lowest Wilson lower bound was 95.8 %, in
*T. pallidum* Grady) (Fig. 4b). The extra cost fell on the DIA-NN search (1.2–5.2× longer), not on
GIPPA, whose run times were similar for the two routes except in *Mycobacterium tuberculosis*.

<!-- src: P2S §2.3 (range; T. pallidum 65.9-77.3 % CONT_ peptides; archaea TE3 via RefSeq names),
§2.4 (records 18.4x-33.5x; kept 96.8 %-99.9 %; lowest Wilson lower bound T. pallidum Grady 0.9578
from PXD065022/metrics.json S9_route_b k/n via fig4_campaign.tsv; DIA-NN about 1.2x for
T. pallidum to 5.2x; gippa run similar except Mtb) -->

**TE4.** Five of the ten Phase 2 failures were the RplL/RplJ ratio (Fig. 4c). With *E. coli*, it
failed in 6 of 15 genomes, or 5 of 9 lineages counting the seven *T. pallidum* strains as one, while
AtpA/AtpD passed in every genome where it could be tested. Pass or fail did not follow peptide
counts, and near the band edges the verdict changed with the abundance measure (iBAQ or Top3). One
candidate explanation, that the N-terminally acetylated L7 form escapes the search, was not tested.

<!-- src: P2S §1 (fail list), §2.5 (TE4 table: fails T. pallidum Madras, V. natriegens, M. xanthus,
C. difficile, M. tuberculosis, E. coli; AtpA/AtpD 8 of 8 plus E. coli; acetylation not measured);
fig4_campaign.tsv panel c (15 genomes, 6 fail) -->

**Other failures.** Of the other five Phase 2 failures, one was a wrong premise of the answer key (a
selenoprotein tail that Bakta annotates as a separate gene). One was a vector caveat stated in the
key in advance (pyrE/pyrF detected in *Saccharolobus islandicus*). One was a test-design limit: a
*Clostridioides difficile* prophage region that spans a whole plasmid, whose circular-shift p is 1
by construction, plus small regions. Two were *T. pallidum* checks: pseudogene distal parts, one
half untestable by construction, and a strain abundance comparison missed on two small
differences whose error is unestimated.

<!-- src: P2S §1 "What the 10 Phase 2 fails are", §3 -->

**Checks that passed.** Among the pre-registered checks that passed were the absence of pHV2
proteins (0 of 8 genes detected) and the pyrE2 check in *Haloferax volcanii* H26, and the ESX-1
and PE/PPE checks in *M. tuberculosis*. PE/PPE genes were detected at 12.7 % against 74.5 % for
other genes (difference −61.8 points [−66.4, −55.4], Newcombe interval; Fisher p = 1.1 × 10⁻⁵⁴).
Both the interval and the test assume independent genes, but PE/PPE genes cluster on the genome.
Prophage regions in *Vibrio natriegens* and *M. xanthus* were under-detected (all four regions
O < E, q ≤ 0.05), and detection differed between the two *V. natriegens* chromosomes (65.3 %
against 38.0 %, difference 27.3 points [24.4, 30.2]). In the four genomes whose prophage regions
were judged, detection was below expectation in 13 of 13 regions and q ≤ 0.05 within its BH family
in 6. That pooled count has two selection problems: it leaves out *B. thetaiotaomicron*'s 16
regions, which were reported only (0 of 16 at q ≤ 0.05), and it includes the *C. difficile*
whole-plasmid region.

<!-- src: P2S §1 (H. volcanii 8 pHV2 P, 9 pyrE2 P; M. tuberculosis 7 ESX-1 P, 8 PE/PPE P), §2.7
(V. natriegens and M. xanthus prophages; chromosomes; 13/13 and 6/13; B. theta 16 regions, 0/16,
reported); campaign-p2c/validation/campaign/PXD057335/REPORT.md row 8 (PE/PPE 19/150 = 12.7 %
[8.3, 18.9]; other 2,923/3,925 = 74.5 % [73.1, 75.8]; Fisher p 1.1e-54); Newcombe hybrid score
interval computed here from the two Wilson intervals: −0.6180 [−0.6641, −0.5538] -->

**Strain variants.** At default settings GIPPA did not recover substitutions between close strains.
In Phase 1b, 11 of 31 NCM3722–MG1655 substitution sites were covered by peptides, but none was
accepted (recall 0/11 [0, 0.26]). In the *T. pallidum* cross-strain arms, 0 of 159 covered sites
were accepted. The reason is the floor of the placement FDR: the one-difference strata held 7 to 51
targets in the cross-strain arms and 13 in Phase 1b, so their q could not fall below 1/T > 0.01 even
with no decoy hit. Exact placements in the same arms had pooled placement FDR estimates of 0.35 %
(NCM3722; error of the estimate unestimated) and 0.06–0.10 % (*T. pallidum*). Release 1.1.1 adds a
run warning for this floor but does not change acceptance.

<!-- src: CAMP/ecoli/REPORT_1b.md verdicts 1b-1 (0/11 [0.00, 0.26]), 1b-5 ("pass (error unestimated)",
0.35 %), "Why recall is zero" (1-difference stratum 13 targets); P2S §2.9 (0 of 159; m = 1 strata
7-51 targets, m = 2 strata 8-20), §2.2 (T. pallidum placement FDR 0.059-0.101 %); docs/adr/ADR-0079 -->

**Novel ORFs.** No novel ORF was called in any of the 17 genomes at the default α = 0.01 (1,581
candidates; Fig. 5a). The class FDR is estimated within each class, and nearly every class held
fewer than 100 candidates, so its q could not reach 0.01; every route-B run carried the
`class_q_floor_above_alpha` warning. Size is not the whole story: the two classes with more than
100 candidates, the N-terminal extensions of *M. xanthus* (115) and *M. tuberculosis* (109), also
stayed above α (minimum class FDR 0.0127 and 0.0118). N-terminal extensions had the lowest class FDR
in every dataset and mostly matched RefSeq starts upstream of Bakta's. Antisense candidates had the
highest (0.14–1.0) and rarely matched RefSeq. We did not loosen α. Fig. 6 shows the best-known
case, *E. coli* prfB, whose N-terminal part is made by a programmed +1 frameshift that Bakta does
not model: four route-B peptides lie upstream of Bakta's start and one in the 25-codon leader, and
both candidates have class FDR 0.143.

<!-- src: P2S §2.6 (0 called; classes < 100; M. xanthus 0.0127 (115 cand.), M. tuberculosis 0.0118
(109 cand.); nterm lowest, antisense 0.14-1.0), §6 (alpha not loosened); out_e2_fix/tables/novel_orfs.tsv
(NOVEL_U00096.3_3035184_m nterm_extension 4 peptides class_fdr 0.142857 no_start_codon;
NOVEL_U00096.3_3036205_m intergenic 1 peptide 0.142857); CAMP/ecoli/REPORT.md row 8 -->

**A post hoc positive-control tier.** To see where true signal sits on this q scale, we defined
positive controls after the runs but before any passer was looked at. They are candidates that
known biology confirms: chiefly the N-terminal part of prfB, plus the infC late start, the PrdB
selenocysteine read-through and RefSeq proteins that Bakta missed (Table 4). The 23 controls are
leads already listed with their q-values in the campaign reports. The three primary prfB controls
had class FDRs of 0.111 to 0.167, and we used the largest, τ_prfB = 0.167, as an exploratory
reporting tier. Because τ is the largest of their q-values, the primary controls pass by
construction; of all 23 controls, 15 passed (65.2 % [44.9, 81.2]). The tier admitted 581 of 1,581
candidates (Fig. 5b). Of the 566 non-control passers, 423 (74.7 % [71.0, 78.1]) overlapped a RefSeq
gene, meaning at least 50 % of one RefSeq CDS on the same strand and replicon; a candidate without a
start codon counts as no match. That share is not a precision measure. 417 of the 566 are
N-terminal extensions, which overlap their own gene by design, and a RefSeq match is a weak proxy
in both directions. GIPPA's own class FDRs imply about 75 false candidates among the 581 (12.9 %),
but that estimate relies on the class-FDR estimator, which has known failures (Section 4.1), and it
assumes the decoy model holds in these small classes. This tier is a campaign analysis, not a
GIPPA call or a change of α: the threshold was chosen from the same data, it is about 17 times α,
the controls are biased toward abundant proteins, and the three prfB controls lack a start codon,
so GIPPA's class rule would reject them at any q. The candidates in this tier are leads for review,
not discoveries.

<!-- src: NPC §2 ("Positive controls (defined before any passer was looked at)"; 23 controls), §3
(τ_prfB 0.1667; 15/23, 65.2 % [44.9, 81.2]), §4.1 (581 pass; 566 non-control; 423/566 74.7 %
[71.0, 78.1]; class table: nterm_extension 417 non-control passers; 74.8 expected false, 12.9 %),
§5 (caveats: same data, 17 × α, RefSeq match a weak precision proxy, controls biased, no start
codon); P2S §2.6 (RefSeq match rule ≥ 50 % same strand) and NPC §6 (same-replicon rule);
Khoa 2026-11-24 request quoted in NPC header -->

**A defect found by the campaign.** The most consequential finding concerned GIPPA itself. In
1.1.0, decoys that could not be paired to one gene were dropped from the gene-level FDR, and each
decoy win was rescaled by 1/f, where f is the paired share. That assumes decoys are lost at random.
When protein identifiers are not locus tags and genes are paired only by the votes of accepted
peptides, as in a search against a related strain, a decoy pairs only if its protein already has
an accepted target, so exactly the decoys that would have won their gene are lost. On the NCM3722
search mapped to the MG1655 genome, 1.1.0 reported 2,399 genes at maximum q 0.0051, while the FDR
estimated with all decoys was about 0.026–0.028. Under the 1.1.1 rule (Section 3.2) the same run
reports 2,305 genes at maximum q 0.0100, and its detected set equals that of the matched-genome run
(Jaccard 1.000). No simulation suite had covered the case, because the simulator paired decoys
through locus tags; a calibration attack for it shipped with the fix. Runs with locus-tag
identifiers kept the same detected genes under the new rule; only their q-values moved slightly.
The campaign also found three crashes or coordinate shifts on genomes with assembly gaps
(Promotech, Promoter Calculator, TerminatorNet), fixed with 1.1.1 and without effect on gap-free
genomes.

<!-- src: docs/adr/ADR-0078 (Context: N2 2,399 genes, max q 0.0051, FDR estimated with all decoys
about 0.026-0.028; Calibration attack); P2S §4 (#606/#607 fixed: 2,305 genes at max q 0.0100;
#611/#612, #614/#616, #615/#617); CAMP/ecoli/REPORT_1b.md (N2 fixed, Jaccard with N1 1.000);
CAMP/ecoli/REPORT.md "Recheck with 1.1.1 code" (same genes; E1 max q 0.00933 → 0.00929);
CHANGELOG.md 1.1.1 Fixes -->

### 4.5 Release readiness

A release script evaluates ten acceptance criteria from the specification and six housekeeping
checks (Table S4). For 1.1.1 the overall verdict was pass. Three criteria are waived, not passed:
the real-data placement FDR, because the strain-database search makes it unmeasurable; the
true-example criterion, because TE4 fails for the RplL/RplJ pair; and the runtime target, because
the measured runtime is outside the specified range (Section 5). Co-abundance against K-12 operons
is deferred by the specification to a full validation. The two calibration criteria pass without a
waiver, but a pass means that no test failed outside the declared expected failures, not that every
method is calibrated in every regime (Section 4.1). The pilot criterion passes because the pilot's
five pass criteria were met; those criteria do not include the pre-registered hypotheses, and H1
was not supported.

<!-- src: RD (commit 3f8793e0, computed 2026-11-27T01:33:41Z, overall pass; rows R1-R10, H1-H6);
docs/adr/ADR-0080 (waivers R2, R5 TE4 only, R7); validation/pilot/REPORT.md §15.6 table (five
bullets; no hypothesis bullet); SUPP Table S4 -->

## 5. Performance

The standard-run benchmark was measured on the 1.0.0-era code and has not been repeated on 1.1.1.
It ran on one FarmShare "wheat" node (Intel Xeon Gold 6148) with 16 threads and 64 GB, on the K-12
genome with simulated peptide tables. A standard run with all modules took 3 h 41 min (mean of
two repeats; uncertainty unestimated), about six times the specification's estimate of 35 min;
core modules alone took 3 h 01 min. Two annotation steps, Promoter Calculator and TerminatorNet,
used about 92 % of summed task time in the first repeat. Peak memory reached 62.4 GiB, 99.8 % of
the allocation; whether that reflects a need or a limit was not established. In the campaign, on 16
CPUs in one job with 1.1.1 code, `gippa run` took 35 min (*T. pallidum*, 1.14 Mb) to 8 h 10 min
(*M. xanthus*, route A), and *M. xanthus* peaked at 60.9 GB. `gippa demo` finishes in about 75 s
without containers or databases.

<!-- src: benchmarks/baseline.json › standard_run.wheat.k12_standard, k12_core (13,232.5 s; 10,873.5 s;
2 repeats; peak 62.35 GiB); docs/adr/ADR-0072 (two repeats; not established); results.md §5
(92 % in two steps; demo 75 s); P2S §2.4 (gippa run A/B times 35 min to 8 h 10; M. xanthus 60.9 GB) -->

## 6. Availability and licences

GIPPA 1.1.1 is MIT-licensed and installs with `pip install gippa==1.1.1`. The source,
specification, methods notes, decision records, calibration suites, pilot aggregates and readiness
records are at https://github.com/khoa-yelo/gippa. Containers are pinned by SHA-256 digest.
`gippa demo` runs a bundled simulated example without network access, and `gippa example fetch`
downloads a 1.09-MB real example (two 200-kb regions of the Nissle 1917 genome with subsets of the
pilot data) for five input routes. [PLACEHOLDER: Zenodo DOI; documentation URL; Bioconda package.]

GIPPA never imports or vendors a third-party tool; CI fails if the source tree imports a GPL,
AGPL, non-commercial or unlicensed tool (Table 5). Users must respect the following licences:

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
licence terms"); methods.md §20 (Tools and licences), §24; validation/example/bundles.tsv (1.09 MB)
via methods.md §20 -->

## 7. Limitations

**What the evidence covers.** All calibration rests on gippasim, a simulator built within the
project, so each calibration claim holds only under its generative model; the attack nulls probe
departures from it but do not remove the dependence. There is no external real-data ground truth:
every real-data FDR in this paper is GIPPA's or DIA-NN's own estimate. The full validation in the
specification was not done (no K-12 benchmark dataset, no FDRBench entrapment, no co-abundance
AUROC), and the sealed EDL933 holdout was not used, so no holdout claim is made. The real-data
evidence is one pilot sample and the campaign, whose answer keys test biology the authors expected,
not error rates. The code, simulator, tests, answer keys and this paper were written by AI agents
within one project, so the validation is not independent (Section 2.7). Most pilot numbers were
measured on 1.0.0-era code. The campaign evidence is not yet merged into the repository's main
branch. [PLACEHOLDER: update once the campaign PRs merge.]

**What GIPPA does not estimate.** Identification error is the search engine's. With DIA-NN input,
GIPPA uses the precursor-level q as the peptide cut-off, so peptide-level FDR is not controlled;
the pilot's free entrapment estimated it at 0.027 before removing near-homolog hits. The
`two_peptide` and `protein_list` tiers have no gene FDR; in simulation the two-peptide rule's
gene-level FDP was 0.054 at q = 0.01. The picked gene FDR has no finite-sample proof. Variant
calls have no family-wise or false-discovery control. Abundance has no test or interval.

**What this means for users.** Read the flags: ORA, determinants and TU-block autocorrelation with
predicted TUs carry flags because they are not calibrated (Section 4.1, items 4 and 10), and island
and operon-discordance claims hold only in stated regimes that no output checks. Novel-ORF calls
from de novo input can be anti-conservative (item 2). Cross-strain variants will usually not be
accepted at defaults, and novel-ORF classes in typical bacterial genomes are too small to reach
q ≤ 0.01, so both modules will usually call nothing (Section 4.4). Region-test p-values are coarse
on small replicons and meaningless for a region that covers its whole replicon. Detection shares
depend on how many runs are pooled and on gradient length, so they are not comparable across
datasets. In archaea, the ribosomal-protein check depends on RefSeq names. The host-dominated
*T. pallidum* searches were not analysed for effects on the decoy model. The RplL/RplJ stoichiometry
check fails in 6 of 15 genomes for reasons not established. The measured runtime is several times
the specification's target.

**What version 1 does not do.** One sample per run, with no tests across samples; no spectrum
searching; no start-codon correction, signal-peptide prediction or PTM naming; no recoding
catalogue (selenocysteine, programmed frameshifts) beyond reporting where peptides fall; and
deferred features listed in the changelog (assembly checks from frameshift-spanning peptides,
proteome mass fractions, RNA-seq or Ribo-seq inputs, homology evidence for novel ORFs).

<!-- src: methods.md §21 (Scope of the evidence); docs/preprint/outline.md §6; results.md §1.2
(two_peptide FDP 0.054; variants), §4; PIL › m8.peptide; changelog.d/+known-limitations.misc.md
(via outline §6); P2S §2.3, §2.6, §2.9, §3, §6 -->

## 8. Discussion

GIPPA makes one claim about itself that we can test: that every number it reports either carries an
error estimate or says it has none. The calibration shows where that holds. In simulation, the
placement FDR, the gene-level FDR, the region tests and the GSEA level checks stay within
tolerance in their gated streams, with exceptions we list: a contaminant attack on placement, a
prophage FDR within tolerance but not at level, and GSEA's detection-confounding screen, which
passes only at 1.5α on a port of the shipped engine. Autocorrelation, islands and operon
discordance hold only in stated regimes. ORA and the determinants test do not hold, and say so on
every row. A user can read the flags and the methods notes and know which column to trust.

The campaign tested something calibration cannot: how the defaults behave on organisms the
developers did not tune for. Its most useful result was a defect. Real cross-strain searches
exposed the decoy-pairing error in the picked gene FDR, which no simulation suite covered because the
simulator paired decoys through locus tags. The fix and a calibration attack for it shipped in
1.1.1. That the defect was found by running real searches, not by reasoning about the estimator,
argues for keeping such campaigns part of release testing. The campaign's consistency checks
(exact mapping on route A, gene-level FDR estimates at the cut-off) confirm that the pipeline is
wired correctly, not that its error estimates are right.

The negative results are informative too. Two error-control rules that are correct in isolation,
a +1 correction in the placement FDR and class-wise decoy competition for novel ORFs, make GIPPA
nearly powerless at genome sizes typical of bacteria for the two discovery tasks that motivate
proteogenomics, strain variants and new genes. GIPPA now warns when a stratum or class cannot reach
α, and we report the post hoc positive-control tier only as a guide to where true signal sits. A
principled fix, such as pooling classes or strata with a calibrated mixture, would change a
statistical rule and needs its own calibration and decision record. We think it is better to ship a
method that says "no call possible here" than one that calls at an error rate it cannot estimate.

[PLACEHOLDER: Khoa: comparison with related tools and positioning; needs citations, e.g.
[CITE: bacterial proteogenomics pipelines], [CITE: DIA-NN-based microbial proteomics workflows],
[CITE: entrapment-based FDR assessment, beyond Wen et al. 2025].]

<!-- src: Sections 3-5 above; docs/adr/ADR-0078 (Context; calibration attack added); P2S §6 (Novel
ORFs: any change needs an ADR and a Khoa decision) -->

## Data and code availability

- **Proteomics data.** Pilot: PXD078459 (iProX). Campaign: PXD058422 (*E. coli* K-12, Phase 1),
  PXD060434 (*E. coli* NCM3722, Phase 1b), and the ten Phase 2 datasets PXD065022, PXD069004,
  PXD055395, PXD072227, PXD049476, PXD060688, PXD065346, PXD057335, PXD066008 and PXD079213
  (Table S3). Cited by accession, not redistributed.
- **Genomes.** Nissle 1917 GCF_021559835.1 (pilot); K-12 MG1655 U00096.3; the campaign genomes are
  listed with their NCBI accessions and checksums in `validation/campaign/genomes.tsv`. Proteome:
  UniProt UP000000625 (pilot K-12-reference arm).
- **Not redistributed.** RegulonDB (its terms), DIA-NN output rows and peptide sequences (only
  aggregates are committed), DIA-NN itself and the Cascadia checkpoint.
- **Code and evidence.** https://github.com/khoa-yelo/gippa (MIT): source, `GIPPA_SPEC.md`,
  `docs/methods/`, decision records in `docs/adr/`, `tests/calibration/`, `validation/pilot/`,
  `validation/release/1.1.1/`, and (once merged) `validation/campaign/`. The figure scripts of this
  paper are in `docs/preprint/figures/`. Calibration records are kept on the development cluster;
  [PLACEHOLDER: archive the CAL111 records with the Zenodo release].
- **Publication of the pilot figures.** [PLACEHOLDER: Khoa, Q-K11: which pilot figures may appear,
  including Supplementary Figures S1 and S2.]

<!-- src: methods.md §24; handover plan Phase 2 table (accessions); CAMP/genomes.tsv (worktree
campaign-p2a); docs/preprint/outline.md §5 -->

## Funding, competing interests, acknowledgements

[PLACEHOLDER: Khoa: funding; competing interests; acknowledgements (Stanford FarmShare).]

## Figure legends

**Figure 1. The GIPPA workflow.** Inputs are a genome, an annotation and peptide or protein
evidence from any search engine; `gippa convert` reads DIA-NN, Cascadia and generic tables.
`gippa makedb` optionally writes a strain database from the Bakta annotation, with or without
six-frame ORFs, for the user's own search, which runs outside GIPPA. Annotation layers come from
third-party tools run as separate programs. Peptides are placed on the six-frame genome with exact
and tolerant matching, and genome decoys give a placement FDR. The gene table holds one row per
Bakta locus tag with a detection tier, gene-level FDR, abundance and detectability. Context,
operon, variant, novel-ORF, function and determinant modules build on it. Solid vermillion boxes
mark where GIPPA reports an error estimate; dashed grey boxes mark steps with no estimate of their
own (the inherited search q, presence calls of annotation layers), which are flagged. Draft
schematic, drawn programmatically.

**Figure 2. Calibration of GIPPA 1.1.1 in simulation.** Full calibration tier on commit `d2bd17e4`,
the run that readiness 1.1.1 reuses. All results are conditional on the gippasim simulator.
Colours denote methods in panels a and c, and line styles and markers separate streams.
(a) Quantile–quantile plots of null p-values (−log10 scale) for four streams where the claim
applies, and the distance-model operon fallback, a known conservative failure; grey, pointwise
95 % band for 1,000 uniform p-values. (b) Each of the 1,244 records with a level, FDR or FWER
criterion at nominal level 0.01, 0.05 or 0.10: observed rate against its bound, both divided by
the nominal level. Points above the solid curve (observed = bound) exceed the bound; the dotted
line marks observed = nominal. Markers: passing records; unconfirmed (above the bound in one block of
two); failing records in arms that are recorded only or whose rows ship flagged; known failures
declared before the run. Records with zero observed error are drawn at 0.02. (c) Observed rejection
rate or mean FDP against the nominal level for representative streams; short bars mark each stream's
bound (nominal + 2 Monte Carlo SE). Points are offset horizontally for legibility.

**Figure 3. The *E. coli* Nissle 1917 pilot.** Three wild-type DIA runs (PXD078459), searched
against the UniProt K-12 proteome (K-12-reference arm, black squares) or the Nissle strain database
(strain-database arm, blue circles). Detection values come from the 1.1.1 detection code applied to
the frozen pilot inputs. (a) Share of genes detected (gene q ≤ 0.01) by subset, with the number of
genes per subset; Wilson 95 % intervals treat genes as independent and are optimistic. The
reference arm cannot detect strain-specific genes. (b) Pre-registered hypothesis H1: the share of
discordant genes detected only with the strain database; the share of those that are
strain-specific; the share that lie in plasmids, prophages or BGCs (vermillion), against the
pre-registered "mostly" requirement (> 0.5, dashed) and the base rate among all genes (grey).
(c) Free entrapment for the reference arm: lower bound, combined estimate (dot) and one-sided 95 %
upper bound of the FDP at precursor level (the level DIA-NN's q controls), at peptide level (not
controlled by that q), and at peptide level without near-homolog and short entrapment hits; dashed,
the 1 % target. The upper bound holds in expectation and assumes independent false matches.

**Figure 4. The campaign on public DIA data.** GIPPA 1.1.1 on eleven datasets (17 genomes; *E. coli*
K-12 from Phase 1 in black). All panels use the same genome order and names. (a) Route-A genes
detected at gene q ≤ 0.01 as a share of coding genes, with counts. (b) Share of route-A detected
genes that route B (database plus six-frame ORFs) also detects; dashed, the pre-registered prior of
0.95. Intervals are Wilson 95 %, optimistic because genes co-vary. (c) TE4: log2(RplL/RplJ) on iBAQ
(filled, the verdict) and Top3 (open), against the band [1, 3] (grey); vermillion triangles, fail.
The archaea are not testable. (d) Outcomes of the pre-registered answer-key checks per dataset
(*T. pallidum* strains as one dataset): pass, fail, not testable or not run, and reported without a
threshold. One dropped *E. coli* check is not shown. Depth comparisons across datasets are
confounded by run count, gradient length and instrument.

**Figure 5. Novel ORFs in the campaign.** (a) For each run and candidate class, the smallest class
FDR against the number of candidates in the class (log scales). Dashed, α = 0.01; dotted, the post
hoc tier τ_prfB = 0.167; solid, 1/n as a guide only (the class q cannot fall below 1/T, where T
counts the class's targets, which is not the plotted n, so some points lie below the line). No
class reaches α. (b) Candidates at or below τ_prfB per run: positive controls (vermillion),
non-control candidates that overlap a RefSeq gene (≥ 50 % of one same-strand CDS on the same
replicon; blue), and those that do not (grey); labels give passing / all candidates. GIPPA called
none of them. The tier was chosen after the runs from the same data and is not an error-controlled
call.

**Figure 6. Genome context at one locus: *E. coli* K-12 prfB.** Phase 1 route-B run (GIPPA 1.1.1
code), U00096.3 positions 3,034,950–3,036,600. (a) Tracks from top: RefSeq CDS (grey; prfB in two
segments joined by the programmed +1 frameshift at codon 26), Bakta CDS (blue; Bakta annotates prfB
from an internal start), GIPPA's two novel-ORF candidates (outlines; neither called), and peptide
coverage as the number of distinct accepted placements covering each position, by locus class.
Coverage is shown instead of peptides; no sequence is drawn. (b) The `genes.tsv` rows of the three
genes in the window (gene q from the `decoy_picked` tier; the q-value 4.0 × 10⁻⁴ is the tier's
floor in this run).

**Supplementary Figures S1 and S2 (pending Q-K11).** GIPPA report plots V01 (genome overview) and
V07 (observed minus expected detection by replicon, feature class, COG category and region) for the
pilot's strain-database arm. They were re-rendered with the 1.1.1 report code from the 1.1.1
detection rerun of the pilot. GC skew and origin marks are omitted because the rerun outputs lack
their support tables.

## Tables

**Table 1.** Documented output tables (`docs/reference/tables.md`, generated from the schemas).

| Location | Tables |
|---|---|
| `mapping/` (4) | `peptide_hits` (peptide × accepted locus), `protein_hits`, `mapping_fdr` (placement FDR strata, de novo precision bins), `unmapped` |
| `tables/` (18) | `genes` (one row per CDS, sORF or pseudogene), `detection_summary` (tier, FDR curve, O − E, QC), `regions_stats`, `autocorrelation`, `islands`, `operons`, `pairs`, `operon_flags`, `pileup`, `variants`, `novel_orfs`, `category_summary`, `enrichment_detection` (ORA), `enrichment_gsea`, `determinants_features`, `determinants_model`, `determinants_attribution`, `sensitivity` |
| `gippa compare` (3) | `compare_detection`, `compare_abundance` (Spearman ρ with TU-block bootstrap CI, Bland–Altman), `compare_stability` (Jaccard indices) |
| `plots/support/` (9; not result tables) | `replicon_landmarks`, `gc_skew`, `variant_windows`, `gsea_running`, `determinants_shap`, `determinants_pd`, `null_regions`, `null_islands`, `null_determinants` |

<!-- src: docs/reference/tables.md (Path and description lines of all 34 sections) -->

Table 2 (error-control families) is in Section 3.

**Table 3.** True-example assertions (pilot values from the 1.0.0-era freeze; TE6 and TE8 from the
1.1.1 weekly run).

| TE | Assertion | Outcome | Value |
|---|---|---|---|
| TE1 | ≥ 99.5 % of strain-database peptides map exactly, in frame (consistency check) | pass | 28,129 / 28,129 (Wilson 95 % lower bound 0.99986) |
| TE2 | no `high` ALT call off start or stop codons | pass, vacuous | 0 ALT sites in the strain-database arm and in the weekly proxy |
| TE3 | ribosomal proteins ≥ 90 % detected, median abundance in the top quartile | pass | 53 / 53 detected [0.932, 1.0]; median percentile rank 0.987 |
| TE4 | log2(RplL/RplJ) ∈ [1, 3]; log2(AtpA/AtpD) ∈ [−1, 1] | **fail (waived for this pair)** | RplL/RplJ −0.480 (Top3 −0.531); AtpA/AtpD −0.118 |
| TE5 | trypsin and keratin peptides land in the contaminant bucket | pass | K-12-reference arm 59 / 59; strain-database arm 61 / 61 |
| TE6 | K-12: Bakta oriC next to mnmG | pass (any-call rule) | 2 oriC calls, 1 at the true origin |
| TE7 | same seed gives byte-identical tables | pass (tables only) | 2,786 gated files, 0 differ; integration test on a simulated genome passes |
| TE8 | UniProt K-12 peptides ≥ 99.5 % exact on K-12; Nissle placement FDR ≤ 1 % | pass | 0.99900 (151,272 / 151,424); Nissle placement FDR 0.00533 |

<!-- src: results.md Table R4; local/ci/validation/2026-11-25_te8_2020593.json › shares.k12_exact
0.99899, counts 151272/151424, shares.nissle_placement_fdr 0.005331; local/logs/gippa-weekly.2020593_0.out:319
("TE6: 2 oriC calls, at the true origin: 1"); RD rows R5, R6 -->

**Table 4.** Positive controls for the novel-ORF tier (route B; class FDR is GIPPA's `class_fdr`;
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

**Table 5.** Third-party tools in GIPPA 1.1.1.

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

Supplementary Tables S1 (command-line options), S2 (representative calibration streams), S3
(campaign datasets) and S4 (release readiness) are in `docs/preprint/supplementary_tables.md`.

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
