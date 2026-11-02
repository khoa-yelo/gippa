# GIPPA pilot report (arms A1, A2 and A3)

Rendered by `validation/pilot/metrics.py render --final` from `metrics.json` (docs/specs/M10-05-pilot-metrics.md, docs/specs/M10-15-pilot-gate.md).
Every value below is followed by the file and key it comes from; none is typed by hand.
This template holds no numbers. Only aggregates appear here: no sequence, no gene list, no per-gene value.

## 1. Status

| §15.6 bullet | Status | Evidence |
|---|---|---|
| GIPPA runs end to end on A1, A2 and A3 | pass (A3 = Cascadia) (metrics.json › checklist.0) | metrics.json › evidence.exit_codes; reproduce.md › Evidence; reproduce.md › S8 1914948 (COMPLETED, exit 0:0, fa49b8f), S9 1915264 (COMPLETED, exit 0:0, fa49b8f), S11 1916086 (COMPLETED, exit 0:0, fa49b8f) (metrics.json › checklist.0) |
| ≥ 99.5 % of A2 peptides map exactly (TE1) | pass (metrics.json › checklist.1) | metrics.json › te1.share_exact_inframe, te1.passed; reproduce.md › S14 1921926 (COMPLETED, exit 0:0, fa49b8f) (metrics.json › checklist.1) |
| A1 and A2 detections and abundances are compared | pass (metrics.json › checklist.2) | metrics.json › m4, m5; reproduce.md › S12 1916515 (COMPLETED, exit 0:0, fa49b8f) (metrics.json › checklist.2) |
| Measured runtimes replace the estimates | pass (metrics.json › checklist.3) | baseline.json › standard_run.wheat.{k12_standard,k12_core,sim5m_standard,sim5m_core}.median_s; ADR-0072-benchmark-baseline-1-0-0 accepted-by-khoa; docs/user-guide/runtimes.md (metrics.json › checklist.3) |
| TE2–TE5 are reported; failures explained | reported (1 not passed, explained) (metrics.json › checklist.4) | metrics.json › te_real.te2-te5.status; TE4 not passed: report_notes.md › Diagnosis: TE4 (metrics.json › checklist.4) |

A failing bullet is explained in the notes below; no threshold was moved to make a bullet pass (R3).
GIPPA commit: fa49b8f145ea4598a5f896f9d0b1956912fc5df3 (metrics.json › gippa_commit). Spec SHA: 87e381908f8cb9e5020e57cc0dd349641425bcae (metrics.json › spec).
Submitted for Khoa's review (§14.8 item 4).

## 2. Data and arms

PXD078459 wild-type runs (M10-01 manifest), DIA-NN 2.3.1 Academia (ADR-0051). A1 searched UniProt
K-12 (UP000000625) plus contaminants; A2 searched the Nissle strain database (`strain_db.faa`) plus
contaminants. A3 is de novo: Cascadia 0.0.7 on the same runs or, when M10-08's verdict rules it
out, Casanovo 5.2.1 on PXD053727 / UTI89 (fallback, a different strain).
A3 branch of the runbook: cascadia (reproduce.md › evidence.runbook); M10-08 verdict pass (metrics.json › evidence.a3_feasibility) (reason: NA (metrics.json › evidence.a3_feasibility)).
Job IDs are in `reproduce.md` (Evidence table) and in the status table above.

## 3. Metrics

### Metric 1: identifications at 1 % FDR (union over runs)

Only the precursor count is at a controlled 1 % FDR (global precursor q-value). Peptide counts are
not FDR-controlled exactly. Protein groups are groups with ≥ 1 accepted precursor (no protein-level FDR control).
Per-run counts are "accepted at global q ≤ 0.01 among the rows DIA-NN reported for that run", not a per-run FDR.

| | A1 | A2 |
|---|---|---|
| precursors | 28639 (metrics.json › m1.A1.union.precursors) | 33764 (metrics.json › m1.A2.union.precursors) |
| peptides (all) | 23792 (metrics.json › m1.A1.union.peptides) | 28234 (metrics.json › m1.A2.union.peptides) |
| peptides (non-contaminant) | 23704 (metrics.json › m1.A1.union.peptides) | 28129 (metrics.json › m1.A2.union.peptides) |
| peptides (contaminant-only) | 88 (metrics.json › m1.A1.union.peptides) | 105 (metrics.json › m1.A2.union.peptides) |
| protein groups | 2798 (metrics.json › m1.A1.union.protein_groups) | 3228 (metrics.json › m1.A2.union.protein_groups) |
| protein groups with ≥ 2 accepted peptides | 2433 (metrics.json › m1.A1.union.protein_groups_min2) | 2790 (metrics.json › m1.A2.union.protein_groups_min2) |

### Metric 2: mapping shares (non-contaminant query peptides)

| | A1 | A2 |
|---|---|---|
| query peptides (U) | 23792 (metrics.json › m2.A1) | 28234 (metrics.json › m2.A2) |
| contaminant-only (C) | 88 (metrics.json › m2.A1) | 105 (metrics.json › m2.A2) |
| denominator (D) | 23704 (metrics.json › m2.A1) | 28129 (metrics.json › m2.A2) |
| exact | 0.9936 [0.9925, 0.9945] (metrics.json › m2.A1.shares) | 1 [0.9999, 1] (metrics.json › m2.A2.shares) |
| tolerant only | 0.005147 [0.004313, 0.006141] (metrics.json › m2.A1.shares) | 0 [0, 0.0001365] (metrics.json › m2.A2.shares) |
| unmapped | 0.001266 [0.0008867, 0.001806] (metrics.json › m2.A1.shares) | 0 [0, 0.0001365] (metrics.json › m2.A2.shares) |

Genome-inconsistent peptides in A1 (`detection_summary.tsv`): rate 0.058 (metrics.json › m2.A1.genome_inconsistent), count 1450 (metrics.json › m2.A1.genome_inconsistent); equal to tolerant + unmapped: no (metrics.json › m2.A1.genome_inconsistent).

### Metric 3: coverage of the strain CDS

Shares with nominal 95 % Wilson intervals, descriptive only. A1 and A2 are paired on the same genes,
so the arms are not compared by interval overlap; hypothesis 1 is the paired comparison.
"Strain-specific" is defined by peptide identity (no shared in-silico tryptic peptide with the K-12
digest, I/L collapsed), a proxy for "no K-12 counterpart"; it does not mean "undetectable by A1",
which can still detect such genes through tolerant placements.

| subset | A1 | A2 |
|---|---|---|
| all | 0.5705 [0.5563, 0.5847] (metrics.json › m3.A1.subsets.all) | 0.6525 [0.6388, 0.666] (metrics.json › m3.A2.subsets.all) |
| chromosome | 0.5719 [0.5577, 0.586] (metrics.json › m3.A1.subsets.chromosome) | 0.6522 [0.6384, 0.6657] (metrics.json › m3.A2.subsets.chromosome) |
| plasmid | 0.07692 [0.01371, 0.3331] (metrics.json › m3.A1.subsets.plasmid) | 0.7692 [0.4974, 0.9182] (metrics.json › m3.A2.subsets.plasmid) |
| plasmid:NZ_CP058218.1 | 0 [0, 0.4345] (metrics.json › m3.A1.subsets.plasmid:NZ_CP058218.1) | 1 [0.5655, 1] (metrics.json › m3.A2.subsets.plasmid:NZ_CP058218.1) |
| plasmid:NZ_CP058219.1 | 0.125 [0.02242, 0.4709] (metrics.json › m3.A1.subsets.plasmid:NZ_CP058219.1) | 0.625 [0.3057, 0.8632] (metrics.json › m3.A2.subsets.plasmid:NZ_CP058219.1) |
| unknown | NA (empty_denominator) (metrics.json › m3.A1.subsets.unknown) | NA (empty_denominator) (metrics.json › m3.A2.subsets.unknown) |
| prophage | 0.02206 [0.00753, 0.06285] (metrics.json › m3.A1.subsets.prophage) | 0.2721 [0.2043, 0.3523] (metrics.json › m3.A2.subsets.prophage) |
| bgc | 0.3108 [0.2567, 0.3705] (metrics.json › m3.A1.subsets.bgc) | 0.5458 [0.484, 0.6063] (metrics.json › m3.A2.subsets.bgc) |
| strain_specific | 0 [0, 0.004176] (metrics.json › m3.A1.subsets.strain_specific) | 0.3461 [0.316, 0.3775] (metrics.json › m3.A2.subsets.strain_specific) |
| k12_shared | 0.7111 [0.6964, 0.7254] (metrics.json › m3.A1.subsets.k12_shared) | 0.7289 [0.7145, 0.7429] (metrics.json › m3.A2.subsets.k12_shared) |
| no_peptides | 0 [0, 0.2775] (metrics.json › m3.A1.subsets.no_peptides) | 0 [0, 0.2775] (metrics.json › m3.A2.subsets.no_peptides) |
| pseudogene | 0 [0, 0.5615] (metrics.json › m3.A1.subsets.pseudogene) | 0 [0, 0.5615] (metrics.json › m3.A2.subsets.pseudogene) |

Subsets: all of G (cds, sorf), by replicon (`plasmid:<contig>` per plasmid; `unknown`: neither
chromosome nor plasmid), by region (prophage, bgc), by K-12 class (`no_peptides`: no in-silico
tryptic peptide); `pseudogene` is outside G.

### Metric 4: quantitative agreement (shared genes, from `gippa compare`)

Spearman rho 0.978 (metrics.json › m4) [0.9731 (metrics.json › m4), 0.9826 (metrics.json › m4)] on 2635 (metrics.json › m4) shared genes (TU-block bootstrap); Bland–Altman bias -0.195 (metrics.json › m4) log2; flags: NA (metrics.json › m4).

### Metric 5: stability of conclusions (Jaccard, from `gippa compare`)

| family | Jaccard | note |
|---|---|---|
| islands | NA (metrics.json › m5.islands) | none in either arm (metrics.json › m5.islands) |
| operon flags | NA (metrics.json › m5.operon_flags) | none in either arm (metrics.json › m5.operon_flags) |
| enriched categories | 0.563 (metrics.json › m5.enriched_categories) | NA (metrics.json › m5.enriched_categories) |

### Metric 6: novel ORFs and variants by class and confidence, A2/A3 overlap (M10-10)

Descriptive, unestimated: counts and Jaccards carry no error estimate (R4); the two shares carry
nominal 95 % Wilson intervals that describe the sample of calls, not an FDR. Per arm, novel ORFs
(`novel_orfs.tsv`) are counted as called, not called with an unestimated class FDR, and not called.

| class | A1 called | A1 unestimated | A1 not_called | A2 called | A2 unestimated | A2 not_called | A3 called | A3 unestimated | A3 not_called |
|---|---|---|---|---|---|---|---|---|---|
| nterm_extension | 0 (metrics.json › m6.A1.novel_orfs.by_class.nterm_extension) | 15 (metrics.json › m6.A1.novel_orfs.by_class.nterm_extension) | 0 (metrics.json › m6.A1.novel_orfs.by_class.nterm_extension) | 0 (metrics.json › m6.A2.novel_orfs.by_class.nterm_extension) | 0 (metrics.json › m6.A2.novel_orfs.by_class.nterm_extension) | 0 (metrics.json › m6.A2.novel_orfs.by_class.nterm_extension) | 0 (metrics.json › m6.A3.novel_orfs.by_class.nterm_extension) | 0 (metrics.json › m6.A3.novel_orfs.by_class.nterm_extension) | 1 (metrics.json › m6.A3.novel_orfs.by_class.nterm_extension) |
| cterm_extension | 0 (metrics.json › m6.A1.novel_orfs.by_class.cterm_extension) | 0 (metrics.json › m6.A1.novel_orfs.by_class.cterm_extension) | 0 (metrics.json › m6.A1.novel_orfs.by_class.cterm_extension) | 0 (metrics.json › m6.A2.novel_orfs.by_class.cterm_extension) | 0 (metrics.json › m6.A2.novel_orfs.by_class.cterm_extension) | 0 (metrics.json › m6.A2.novel_orfs.by_class.cterm_extension) | 0 (metrics.json › m6.A3.novel_orfs.by_class.cterm_extension) | 0 (metrics.json › m6.A3.novel_orfs.by_class.cterm_extension) | 0 (metrics.json › m6.A3.novel_orfs.by_class.cterm_extension) |
| intergenic | 1 (metrics.json › m6.A1.novel_orfs.by_class.intergenic) | 6 (metrics.json › m6.A1.novel_orfs.by_class.intergenic) | 0 (metrics.json › m6.A1.novel_orfs.by_class.intergenic) | 0 (metrics.json › m6.A2.novel_orfs.by_class.intergenic) | 0 (metrics.json › m6.A2.novel_orfs.by_class.intergenic) | 0 (metrics.json › m6.A2.novel_orfs.by_class.intergenic) | 0 (metrics.json › m6.A3.novel_orfs.by_class.intergenic) | 0 (metrics.json › m6.A3.novel_orfs.by_class.intergenic) | 1 (metrics.json › m6.A3.novel_orfs.by_class.intergenic) |
| cds_outframe | 0 (metrics.json › m6.A1.novel_orfs.by_class.cds_outframe) | 0 (metrics.json › m6.A1.novel_orfs.by_class.cds_outframe) | 0 (metrics.json › m6.A1.novel_orfs.by_class.cds_outframe) | 0 (metrics.json › m6.A2.novel_orfs.by_class.cds_outframe) | 0 (metrics.json › m6.A2.novel_orfs.by_class.cds_outframe) | 0 (metrics.json › m6.A2.novel_orfs.by_class.cds_outframe) | 0 (metrics.json › m6.A3.novel_orfs.by_class.cds_outframe) | 0 (metrics.json › m6.A3.novel_orfs.by_class.cds_outframe) | 1 (metrics.json › m6.A3.novel_orfs.by_class.cds_outframe) |
| antisense | 0 (metrics.json › m6.A1.novel_orfs.by_class.antisense) | 0 (metrics.json › m6.A1.novel_orfs.by_class.antisense) | 0 (metrics.json › m6.A1.novel_orfs.by_class.antisense) | 0 (metrics.json › m6.A2.novel_orfs.by_class.antisense) | 0 (metrics.json › m6.A2.novel_orfs.by_class.antisense) | 0 (metrics.json › m6.A2.novel_orfs.by_class.antisense) | 0 (metrics.json › m6.A3.novel_orfs.by_class.antisense) | 0 (metrics.json › m6.A3.novel_orfs.by_class.antisense) | 2 (metrics.json › m6.A3.novel_orfs.by_class.antisense) |
| pseudogene | 0 (metrics.json › m6.A1.novel_orfs.by_class.pseudogene) | 0 (metrics.json › m6.A1.novel_orfs.by_class.pseudogene) | 0 (metrics.json › m6.A1.novel_orfs.by_class.pseudogene) | 0 (metrics.json › m6.A2.novel_orfs.by_class.pseudogene) | 0 (metrics.json › m6.A2.novel_orfs.by_class.pseudogene) | 0 (metrics.json › m6.A2.novel_orfs.by_class.pseudogene) | 0 (metrics.json › m6.A3.novel_orfs.by_class.pseudogene) | 0 (metrics.json › m6.A3.novel_orfs.by_class.pseudogene) | 0 (metrics.json › m6.A3.novel_orfs.by_class.pseudogene) |
| sorf | 0 (metrics.json › m6.A1.novel_orfs.by_class.sorf) | 0 (metrics.json › m6.A1.novel_orfs.by_class.sorf) | 0 (metrics.json › m6.A1.novel_orfs.by_class.sorf) | 0 (metrics.json › m6.A2.novel_orfs.by_class.sorf) | 0 (metrics.json › m6.A2.novel_orfs.by_class.sorf) | 0 (metrics.json › m6.A2.novel_orfs.by_class.sorf) | 0 (metrics.json › m6.A3.novel_orfs.by_class.sorf) | 0 (metrics.json › m6.A3.novel_orfs.by_class.sorf) | 0 (metrics.json › m6.A3.novel_orfs.by_class.sorf) |
| total | 1 (metrics.json › m6.A1.novel_orfs.totals) | 21 (metrics.json › m6.A1.novel_orfs.totals) | 0 (metrics.json › m6.A1.novel_orfs.totals) | 0 (metrics.json › m6.A2.novel_orfs.totals) | 0 (metrics.json › m6.A2.novel_orfs.totals) | 0 (metrics.json › m6.A2.novel_orfs.totals) | 0 (metrics.json › m6.A3.novel_orfs.totals) | 0 (metrics.json › m6.A3.novel_orfs.totals) | 5 (metrics.json › m6.A3.novel_orfs.totals) |

ALT calls (`variants.tsv`, `call = ALT`) by type and confidence; for A3 every class error rate is
expected unestimated (no spectrum pairing).

| ALT type | A1 high | A1 candidate | A2 high | A2 candidate | A3 high | A3 candidate |
|---|---|---|---|---|---|---|
| substitution | 4 (metrics.json › m6.A1.variants.alt.substitution) | 114 (metrics.json › m6.A1.variants.alt.substitution) | 0 (metrics.json › m6.A2.variants.alt.substitution) | 0 (metrics.json › m6.A2.variants.alt.substitution) | 0 (metrics.json › m6.A3.variants.alt.substitution) | 41 (metrics.json › m6.A3.variants.alt.substitution) |
| swap | 0 (metrics.json › m6.A1.variants.alt.swap) | 0 (metrics.json › m6.A1.variants.alt.swap) | 0 (metrics.json › m6.A2.variants.alt.swap) | 0 (metrics.json › m6.A2.variants.alt.swap) | 1 (metrics.json › m6.A3.variants.alt.swap) | 125 (metrics.json › m6.A3.variants.alt.swap) |
| block | 0 (metrics.json › m6.A1.variants.alt.block) | 0 (metrics.json › m6.A1.variants.alt.block) | 0 (metrics.json › m6.A2.variants.alt.block) | 0 (metrics.json › m6.A2.variants.alt.block) | 0 (metrics.json › m6.A3.variants.alt.block) | 28 (metrics.json › m6.A3.variants.alt.block) |
| frameshift | 0 (metrics.json › m6.A1.variants.alt.frameshift) | 0 (metrics.json › m6.A1.variants.alt.frameshift) | 0 (metrics.json › m6.A2.variants.alt.frameshift) | 0 (metrics.json › m6.A2.variants.alt.frameshift) | 0 (metrics.json › m6.A3.variants.alt.frameshift) | 0 (metrics.json › m6.A3.variants.alt.frameshift) |

A1: mixed calls 1 (metrics.json › m6.A1.variants); class error rate numeric 0 (metrics.json › m6.A1.variants.class_error), unestimated 119 (metrics.json › m6.A1.variants.class_error).
A2: mixed calls 0 (metrics.json › m6.A2.variants); class error rate numeric 0 (metrics.json › m6.A2.variants.class_error), unestimated 0 (metrics.json › m6.A2.variants.class_error).
A3: mixed calls 53 (metrics.json › m6.A3.variants); class error rate numeric 0 (metrics.json › m6.A3.variants.class_error), unestimated 248 (metrics.json › m6.A3.variants.class_error).

A2/A3 overlap. Novel ORFs key on the stop codon (contig, strand, end), which does not depend on the
peptides; ALT sites key on (contig, codon start, strand, ALT residue) as `gippa compare` does. The
`high` set is each arm's own high-confidence ALT rows.

| set | A2 | A3 | shared | union | Jaccard | A3 confirmed by A2 | A2 found by A3 |
|---|---|---|---|---|---|---|---|
| novel ORFs, called | 0 (metrics.json › m6.overlap_a2_a3.novel_orfs.called) | 0 (metrics.json › m6.overlap_a2_a3.novel_orfs.called) | 0 (metrics.json › m6.overlap_a2_a3.novel_orfs.called) | 0 (metrics.json › m6.overlap_a2_a3.novel_orfs.called) | not evaluable (both_empty) | NA (empty_denominator) (metrics.json › m6.overlap_a2_a3.novel_orfs.called) | NA (empty_denominator) (metrics.json › m6.overlap_a2_a3.novel_orfs.called) |
| novel ORFs, any row | 0 (metrics.json › m6.overlap_a2_a3.novel_orfs.any_row) | 5 (metrics.json › m6.overlap_a2_a3.novel_orfs.any_row) | 0 (metrics.json › m6.overlap_a2_a3.novel_orfs.any_row) | 5 (metrics.json › m6.overlap_a2_a3.novel_orfs.any_row) | 0 (metrics.json › m6.overlap_a2_a3.novel_orfs.any_row) | 0 [0, 0.4345] (metrics.json › m6.overlap_a2_a3.novel_orfs.any_row) | NA (empty_denominator) (metrics.json › m6.overlap_a2_a3.novel_orfs.any_row) |
| ALT sites | 0 (metrics.json › m6.overlap_a2_a3.alt_sites.alt) | 195 (metrics.json › m6.overlap_a2_a3.alt_sites.alt) | 0 (metrics.json › m6.overlap_a2_a3.alt_sites.alt) | 195 (metrics.json › m6.overlap_a2_a3.alt_sites.alt) | 0 (metrics.json › m6.overlap_a2_a3.alt_sites.alt) | 0 [0, 0.01932] (metrics.json › m6.overlap_a2_a3.alt_sites.alt) | NA (empty_denominator) (metrics.json › m6.overlap_a2_a3.alt_sites.alt) |
| ALT sites, high confidence | 0 (metrics.json › m6.overlap_a2_a3.alt_sites.high) | 1 (metrics.json › m6.overlap_a2_a3.alt_sites.high) | 0 (metrics.json › m6.overlap_a2_a3.alt_sites.high) | 1 (metrics.json › m6.overlap_a2_a3.alt_sites.high) | 0 (metrics.json › m6.overlap_a2_a3.alt_sites.high) | 0 [0, 0.7935] (metrics.json › m6.overlap_a2_a3.alt_sites.high) | NA (empty_denominator) (metrics.json › m6.overlap_a2_a3.alt_sites.high) |

Novel ORFs called in both arms with the same class: 0 (metrics.json › m6.overlap_a2_a3.novel_orfs.called).
`gippa compare`'s novel-ORF Jaccard (keyed on `orf_id`, which can split one ORF): NA (metrics.json › m6.overlap_a2_a3.novel_orfs).
High-confidence ALT Jaccard: 0 (metrics.json › m6.overlap_a2_a3.alt_sites.high); called novel-ORF Jaccard: not evaluable (both_empty).

### Metric 8: A1's free entrapment (from `entrapment.py`, M10-07)

The entrapment is free: it is the K-12-only part of A1's own database, the K-12 peptides absent from
Nissle's translation after the I/L and mass-equivalent rules. Estimator: Wen et al. (combined
method; reference in the methods note). r = K-12-only / shared in-silico peptides = 0.307 (metrics.json › m8.universe) (35523 (metrics.json › m8.universe) K-12-only, 115731 (metrics.json › m8.universe) shared); at precursor level
r_prec = 0.304 (metrics.json › m8.universe) (from the configured charge and m/z bounds, not the library).
The assessment is made only at precursor level, the unit DIA-NN's global q-value controls; the peptide level is §15.5's number and is not assessed.

| | peptide (§15.5) | precursor (assessed) |
|---|---|---|
| N_T (shared) | 23517 (metrics.json › m8.peptide) | 28347 (metrics.json › m8.precursor) |
| N_E (K-12-only) | 152 (metrics.json › m8.peptide) | 157 (metrics.json › m8.precursor) |
| excluded as contaminant | 120 (metrics.json › m8.peptide) | 133 (metrics.json › m8.precursor) |
| not in the in-silico universe | 0 (metrics.json › m8.peptide) | 0 (metrics.json › m8.precursor) |
| `fdp_combined` (estimated upper bound in expectation) | 0.0273 (metrics.json › m8.peptide) | 0.0236 (metrics.json › m8.precursor) |
| `fdp_lower` (realised lower bound, conditional on the caveats) | 0.00642 (metrics.json › m8.peptide) | 0.00551 (metrics.json › m8.precursor) |
| with 95 % confidence, under equal chance, the FDP is at most | 0.0308 (metrics.json › m8.peptide.ci) | 0.0265 (metrics.json › m8.precursor.ci) |
| assessment | not_assessed (metrics.json › m8.peptide) | inconclusive (metrics.json › m8.precursor) |

Caveats:

- near-homolog entrapment: the estimate errs high (§15.5) (metrics.json › m8.caveats)
- an estimated upper bound in expectation (Wen 2025, Suppl. Note 2.2), not a bound on this sample (metrics.json › m8.caveats)
- assumes the Nissle genome is complete, is the sample's strain, and that no K-12 material is in the sample (metrics.json › m8.caveats)
- assessment is made at precursor level, the unit DIA-NN's global q controls (metrics.json › m8.caveats)

The 95 % bound assumes equal chance and independent false matches; not valid under clustered (protein-level) errors (metrics.json › m8.precursor.ci).
Warning codes of `entrapment.py` (details in `entrapment.json`; a digest mismatch is a lead finding): none (metrics.json › m8).

## 4. TE1 on real data (A2)

Every A2 target identification, true or false, is a substring of an annotated protein, so TE1 checks
that database, genome and mapping agree; it says nothing about identification error.
The expected misses are structural (N-terminal products of non-ATG starts, Sec sites); anything else is a mapping defect.

Share exact in frame at an own gene: 1 (metrics.json › te1) of 28129 (metrics.json › te1) peptides, nominal 95 % Wilson interval [0.9999 (metrics.json › te1), 1 (metrics.json › te1)] (misses cluster by gene, so its coverage is optimistic; the interval is not a test).
Pass (point estimate at or above the threshold): yes (metrics.json › te1); lower end below the threshold: no (metrics.json › te1).

Misses (0 (metrics.json › te1), first matching cause): non-ATG N-terminal 0 (metrics.json › te1.misses), exact at another gene 0 (metrics.json › te1.misses), exact in another class none (metrics.json › te1.misses), tolerant only 0 (metrics.json › te1.misses), unplaced none (metrics.json › te1.misses); masked exact in frame 18 (metrics.json › te1).
Non-ATG N-terminal products in the denominator: 42 (metrics.json › te1.non_atg_nterm). Peptides with no locus tag of `genes.tsv` among their sources: 0 (metrics.json › te1); with only NA q-values: 0 (metrics.json › te1).

Secondary, never the verdict: excluding non-ATG products 1 [0.9999, 1] (metrics.json › te1); without peptides whose accepted rows are all at the cut-off NA (no_crossings) (metrics.json › te1).
Removed not evaluable (no_crossings) peptides (A2 q_crossings 0 (metrics.json › te1.share_without_q_at_cutoff)); a pass that becomes a fail under this superset: no (metrics.json › te1.share_without_q_at_cutoff).
A1 diagnostic (exact in frame at any gene, no pass flag): 0.9921 [0.9909, 0.9932] (metrics.json › te1_a1_diagnostic).

## 5. True examples (TE2–TE5) and determinism (TE7)

From `te_real.json` and `te7.json` (docs/specs/M10-06-te-real-and-te7.md); every definition was fixed
before the real data were read (R3). Shares carry nominal 95 % Wilson intervals, descriptive only.
Every TE that did not pass has a `## Diagnosis: TE<n>` section in the notes below, naming its cause class.

| TE | Status | Headline |
|---|---|---|
| TE2: no unexplained `high` ALT call in A2 | pass (metrics.json › te_real.te2) | violations 0 (metrics.json › te_real.te2.counts) |
| TE3: ribosomal proteins detected (≥ 90 %) and in the top quartile | pass (metrics.json › te_real.te3) | detected 1 [0.9324, 1] (metrics.json › te_real.te3.detection); median percentile rank 0.987 (metrics.json › te_real.te3.abundance) |
| TE4: operon stoichiometry (iBAQ, molar scale) | fail (metrics.json › te_real.te4) | log2 RplL/RplJ -0.48 (metrics.json › te_real.te4.pairs.rplL/rplJ) in [1, 3]; log2 AtpA/AtpD -0.118 (metrics.json › te_real.te4.pairs.atpA/atpD) in [−1, 1] |
| TE5: trypsin and keratin peptides in the contaminant bucket (A1 and A2) | pass (metrics.json › te_real.te5) | in the bucket: A1 1 [0.9389, 1] (metrics.json › te_real.te5.A1), A2 1 [0.9408, 1] (metrics.json › te_real.te5.A2) |
| TE7: a second A2 run gives byte-identical tables | pass (metrics.json › te7.gippa) | gate files compared 2786 (metrics.json › te7.gippa) |

TE7 is decided by GIPPA's own tables only; the DIA-NN layer (a second search of the same spectra) is an
upstream finding, never the verdict: not_run (metrics.json › te7.diann). Raw data to GIPPA tables reproduce byte for byte: no (metrics.json › te7).
A `known_lapack` difference is still a fail (M4-20 decision 9; writer rounding is the decision 10 backlog item).
TE7 is not a §15.6 bullet: a TE7 fail is shown here and in section 9 with its diagnosis and listed as an open item.

## 6. Hypotheses

Test family for this issue: {H1}, size 1, no BH.

**H1** (A2 detects more strain CDS than A1, mostly plasmid, prophage and BGC genes): verdict not_supported (metrics.json › h1); failed clauses: mostly (metrics.json › h1).
`not_supported` means not supported by the pre-registered rule; `not_evaluable` means the rule could not be applied (no discordant TU).
A2-only genes b 417 (metrics.json › h1), A1-only genes c 33 (metrics.json › h1); A2-only share of discordant genes 0.9267 [0.8988, 0.9473] (metrics.json › h1); coverage difference 0.082 (metrics.json › h1).
These are nominal 95 % Wilson intervals; they describe the effects and are not the test (a
two-sided 95 % interval and the one-sided 0.05 test can disagree near the boundary).
TU-level test: 275 (metrics.json › h1) against 31 (metrics.json › h1) units (1 (metrics.json › h1) mixed, left out), p_tu 2.5e-50 (metrics.json › h1), unit effect 0.8987 [0.8598, 0.9277] (metrics.json › h1).
Gene-level p_gene 4.6e-86 (metrics.json › h1) (assumes independent genes; decides nothing).
Share of A2-only genes in plasmid, prophage or BGC 0.2542 [0.2148, 0.2981] (metrics.json › h1) (base rate 0.08538 [0.07771, 0.09372] (metrics.json › h1)); strain-specific share of A2-only genes 0.7602 [0.717, 0.7987] (metrics.json › h1).
The arms are two database searches of the same spectra, so discordance is largely fixed by sequence
divergence rather than sampled, and dependence at the scale of a region remains; p is not a replication probability.

**H3** (context conclusions agree between A1 and A2): verdict supported (metrics.json › h3); rho criterion holds (metrics.json › h3.criteria.rho); islands not_evaluable (metrics.json › h3.criteria.islands); operon flags not_evaluable (metrics.json › h3.criteria.operon_flags); enriched categories holds (metrics.json › h3.criteria.enriched_categories).
The numbers are shown whatever the verdict: rho 0.978 (metrics.json › h3.criteria.rho) [0.9731 (metrics.json › h3.criteria.rho), 0.9826 (metrics.json › h3.criteria.rho)].

| family | n_a | n_b | n_shared | union | Jaccard | status |
|---|---|---|---|---|---|---|
| islands | 0 (metrics.json › h3.criteria.islands) | 0 (metrics.json › h3.criteria.islands) | 0 (metrics.json › h3.criteria.islands) | 0 (metrics.json › h3.criteria.islands) | NA (metrics.json › h3.criteria.islands) | not_evaluable (metrics.json › h3.criteria.islands) |
| operon flags | 0 (metrics.json › h3.criteria.operon_flags) | 0 (metrics.json › h3.criteria.operon_flags) | 0 (metrics.json › h3.criteria.operon_flags) | 0 (metrics.json › h3.criteria.operon_flags) | NA (metrics.json › h3.criteria.operon_flags) | not_evaluable (metrics.json › h3.criteria.operon_flags) |
| enriched categories | 437 (metrics.json › h3.criteria.enriched_categories) | 285 (metrics.json › h3.criteria.enriched_categories) | 260 (metrics.json › h3.criteria.enriched_categories) | 462 (metrics.json › h3.criteria.enriched_categories) | 0.563 (metrics.json › h3.criteria.enriched_categories) | holds (metrics.json › h3.criteria.enriched_categories) |

**H2** (the strain-specific database lowers genome-inconsistent identifications): deferred (full validation, FDRBench); shown side by side only: A1 genome-inconsistent share 0.058 (metrics.json › m2.A1.genome_inconsistent) and the free entrapment estimate (peptide level) 0.0273 (metrics.json › m8.peptide).

**H4** (genome-anchored de novo precision against the precision measured on A2): verdict consistent_with_lower_bound (metrics.json › h4).
DIA gives no spectrum pairing, so the only measure against A2 is sequence-level: the share of a
bin's de novo queries that A2 identified at 1 % FDR (I = L), a lower bound on the bin's precision.
The data can therefore count against H4 but never confirm it. A bin is reportable when GIPPA reports
its precision and it holds at least 20 queries; the verdict uses the Bonferroni-adjusted relation
over the 10 (metrics.json › h4) reportable bins (`relation_adj`), the 95 % relation is descriptive.

| bin | n | A2-confirmed | genome precision | CI low | CI high | c_value | U (rescaled, descriptive) | L (adjusted) | relation 95 % | relation adjusted | nested share | A2 confirmation among searchable peptides (descriptive) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 (metrics.json › h4.bins.0) | 137 (metrics.json › h4.bins.0) | 0.2628 [0.1963, 0.3422] (metrics.json › h4.bins.0) | 0.293 (metrics.json › h4.bins.0.genome) | 0.2194 (metrics.json › h4.bins.0.genome) | 0.3658 (metrics.json › h4.bins.0.genome) | 0.998 (metrics.json › h4.bins.0.genome) | 0.397 (metrics.json › h4.bins.0) | 0.2628 [0.1722, 0.3791] (metrics.json › h4.bins.0) | consistent (metrics.json › h4.bins.0) | consistent (metrics.json › h4.bins.0) | 0.0219 [0.007475, 0.0624] (metrics.json › h4.bins.0) | 1 [0.9036, 1] (metrics.json › h4.bins.0) |
| 1 (metrics.json › h4.bins.1) | 136 (metrics.json › h4.bins.1) | 0.3603 [0.2845, 0.4438] (metrics.json › h4.bins.1) | 0.352 (metrics.json › h4.bins.1.genome) | 0.2768 (metrics.json › h4.bins.1.genome) | 0.4348 (metrics.json › h4.bins.1.genome) | 0.998 (metrics.json › h4.bins.1.genome) | 0.471 (metrics.json › h4.bins.1) | 0.3603 [0.2553, 0.4806] (metrics.json › h4.bins.1) | consistent (metrics.json › h4.bins.1) | consistent (metrics.json › h4.bins.1) | 0 [0, 0.02747] (metrics.json › h4.bins.1) | 1 [0.9273, 1] (metrics.json › h4.bins.1) |
| 2 (metrics.json › h4.bins.2) | 136 (metrics.json › h4.bins.2) | 0.375 [0.2981, 0.4587] (metrics.json › h4.bins.2) | 0.395 (metrics.json › h4.bins.2.genome) | 0.3116 (metrics.json › h4.bins.2.genome) | 0.4741 (metrics.json › h4.bins.2.genome) | 0.998 (metrics.json › h4.bins.2.genome) | 0.508 (metrics.json › h4.bins.2) | 0.375 [0.2683, 0.4953] (metrics.json › h4.bins.2) | consistent (metrics.json › h4.bins.2) | consistent (metrics.json › h4.bins.2) | 0.007353 [0.001299, 0.04047] (metrics.json › h4.bins.2) | 1 [0.93, 1] (metrics.json › h4.bins.2) |
| 3 (metrics.json › h4.bins.3) | 137 (metrics.json › h4.bins.3) | 0.4745 [0.3927, 0.5576] (metrics.json › h4.bins.3) | 0.497 (metrics.json › h4.bins.3.genome) | 0.4097 (metrics.json › h4.bins.3.genome) | 0.578 (metrics.json › h4.bins.3.genome) | 0.998 (metrics.json › h4.bins.3.genome) | 0.613 (metrics.json › h4.bins.3) | 0.4745 [0.3594, 0.5923] (metrics.json › h4.bins.3) | consistent (metrics.json › h4.bins.3) | consistent (metrics.json › h4.bins.3) | 0 [0, 0.02728] (metrics.json › h4.bins.3) | 1 [0.9442, 1] (metrics.json › h4.bins.3) |
| 4 (metrics.json › h4.bins.4) | 136 (metrics.json › h4.bins.4) | 0.5221 [0.4387, 0.6042] (metrics.json › h4.bins.4) | 0.531 (metrics.json › h4.bins.4.genome) | 0.4458 (metrics.json › h4.bins.4.genome) | 0.619 (metrics.json › h4.bins.4.genome) | 0.998 (metrics.json › h4.bins.4.genome) | 0.657 (metrics.json › h4.bins.4) | 0.5221 [0.4039, 0.6378] (metrics.json › h4.bins.4) | consistent (metrics.json › h4.bins.4) | consistent (metrics.json › h4.bins.4) | 0.02206 [0.00753, 0.06285] (metrics.json › h4.bins.4) | 1 [0.9487, 1] (metrics.json › h4.bins.4) |
| 5 (metrics.json › h4.bins.5) | 136 (metrics.json › h4.bins.5) | 0.5074 [0.4243, 0.59] (metrics.json › h4.bins.5) | 0.527 (metrics.json › h4.bins.5.genome) | 0.4349 (metrics.json › h4.bins.5.genome) | 0.6136 (metrics.json › h4.bins.5.genome) | 0.998 (metrics.json › h4.bins.5.genome) | 0.651 (metrics.json › h4.bins.5) | 0.5074 [0.39, 0.6239] (metrics.json › h4.bins.5) | consistent (metrics.json › h4.bins.5) | consistent (metrics.json › h4.bins.5) | 0.01471 [0.004042, 0.05203] (metrics.json › h4.bins.5) | 1 [0.9473, 1] (metrics.json › h4.bins.5) |
| 6 (metrics.json › h4.bins.6) | 137 (metrics.json › h4.bins.6) | 0.635 [0.5518, 0.7109] (metrics.json › h4.bins.6) | 0.651 (metrics.json › h4.bins.6.genome) | 0.5706 (metrics.json › h4.bins.6.genome) | 0.7243 (metrics.json › h4.bins.6.genome) | 0.998 (metrics.json › h4.bins.6.genome) | 0.756 (metrics.json › h4.bins.6) | 0.635 [0.5152, 0.7402] (metrics.json › h4.bins.6) | consistent (metrics.json › h4.bins.6) | consistent (metrics.json › h4.bins.6) | 0 [0, 0.02728] (metrics.json › h4.bins.6) | 1 [0.9577, 1] (metrics.json › h4.bins.6) |
| 7 (metrics.json › h4.bins.7) | 136 (metrics.json › h4.bins.7) | 0.6985 [0.6168, 0.7693] (metrics.json › h4.bins.7) | 0.727 (metrics.json › h4.bins.7.genome) | 0.6485 (metrics.json › h4.bins.7.genome) | 0.8023 (metrics.json › h4.bins.7.genome) | 0.998 (metrics.json › h4.bins.7.genome) | 0.835 (metrics.json › h4.bins.7) | 0.6985 [0.5797, 0.7956] (metrics.json › h4.bins.7) | consistent (metrics.json › h4.bins.7) | consistent (metrics.json › h4.bins.7) | 0.01471 [0.004042, 0.05203] (metrics.json › h4.bins.7) | 1 [0.9611, 1] (metrics.json › h4.bins.7) |
| 8 (metrics.json › h4.bins.8) | 136 (metrics.json › h4.bins.8) | 0.7132 [0.6322, 0.7826] (metrics.json › h4.bins.8) | 0.752 (metrics.json › h4.bins.8.genome) | 0.678 (metrics.json › h4.bins.8.genome) | 0.8181 (metrics.json › h4.bins.8.genome) | 0.998 (metrics.json › h4.bins.8.genome) | 0.847 (metrics.json › h4.bins.8) | 0.7132 [0.5951, 0.808] (metrics.json › h4.bins.8) | consistent (metrics.json › h4.bins.8) | consistent (metrics.json › h4.bins.8) | 0.02206 [0.00753, 0.06285] (metrics.json › h4.bins.8) | 1 [0.9619, 1] (metrics.json › h4.bins.8) |
| 9 (metrics.json › h4.bins.9) | 137 (metrics.json › h4.bins.9) | 0.8832 [0.8187, 0.9268] (metrics.json › h4.bins.9) | 0.893 (metrics.json › h4.bins.9.genome) | 0.834 (metrics.json › h4.bins.9.genome) | 0.9438 (metrics.json › h4.bins.9.genome) | 0.998 (metrics.json › h4.bins.9.genome) | 0.966 (metrics.json › h4.bins.9) | 0.8832 [0.7846, 0.9401] (metrics.json › h4.bins.9) | consistent (metrics.json › h4.bins.9) | consistent (metrics.json › h4.bins.9) | 0.0365 [0.01569, 0.08259] (metrics.json › h4.bins.9) | 1 [0.9692, 1] (metrics.json › h4.bins.9) |

Descriptive, both over the reportable bins: n-weighted mean genome precision 0.5618 (metrics.json › h4); pooled A2-confirmed share 0.5433 [0.5167, 0.5695] (metrics.json › h4).
Searchable peptides: GIPPA's digest (trypsin_p, at most 1 missed cleavage, 7-30 residues, no N-terminal Met excision; not DIA-NN's settings) of --strain-db, I folded to L (metrics.json › h4). A2's contaminant share of its accepted target rows (counted in c_value): 0.002406 [0.002107, 0.002746] (metrics.json › h4).

Caveats:

- The score bins span the de novo scores of (0.8, 1] only (Cascadia's reported range). (metrics.json › h4.caveats)
- c_value comes from A2's own exact matches and the confirmations from A2's identifications, so A2 enters both sides: a shared dependence, not an independent check. (metrics.json › h4.caveats)
- Class error rates of A3's ALT calls are unestimated (no spectrum pairing). (metrics.json › h4.caveats)
- Dependent rows: nested or overlapping de novo peptides of one locus share correctness, so the genome bin CI can be too narrow and a not_supported can come from that undercoverage; nested_share shows how clustered each bin is. (metrics.json › h4.caveats)
- Contaminant rows count in c_value (they never match the genome), which biases the genome precision up and lowers the power to refute; de novo contaminant peptides are not in A2's set and count as unconfirmed. (metrics.json › h4.caveats)
- A2's identifications are at 1 % FDR, so the lower bound can be off by up to that share (not corrected). (metrics.json › h4.caveats)
- A failure confined to a single score bin reads consistent_with_lower_bound: not_supported needs two consecutive reportable bins below the bound (the pre-registered fallback). (metrics.json › h4.caveats)

## 7. Not evaluated

One line per metric or hypothesis whose value is null with a reason: the pilot data do not allow it (§15.7 "reported where the pilot data allow").

- `te1.share_without_q_at_cutoff`: not evaluable (no_crossings)

## 8. Runtimes

### Metric 7: runtime and peak memory (pilot measurements on a shared cluster, not the benchmark)

| | A1 | A2 |
|---|---|---|
| wall time | 2:36:07 (metrics.json › m7.A1) | 2:50:20 (metrics.json › m7.A2) |
| sum of step wall times | 10:05:49 (metrics.json › m7.A1) | 12:05:36 (metrics.json › m7.A2) |
| max step peak RSS (one task, not the arm's total memory) | 7.5 GiB (metrics.json › m7.A1) | 7.5 GiB (metrics.json › m7.A2) |
| complete (no cached step; else totals are lower bounds) | yes (metrics.json › m7.A1) | yes (metrics.json › m7.A2) |

### Benchmark (M10-12, `benchmarks/baseline.json`, wheat nodes)

Medians of the standard run against the estimates of §13.2 (35 min, range 16–55, and 25 min); a median at an estimate is reported, not judged.

| cell | median wall time | commit | CPU |
|---|---|---|---|
| K-12, standard | 3:40:32 (metrics.json › evidence.benchmark.cells.k12_standard) | fa49b8f145ea4598a5f896f9d0b1956912fc5df3 (metrics.json › evidence.benchmark.cells.k12_standard) | Intel Xeon Gold 6148 (CPU_SKU:6148) (metrics.json › evidence.benchmark.cells.k12_standard) |
| K-12, core | 3:01:14 (metrics.json › evidence.benchmark.cells.k12_core) | fa49b8f145ea4598a5f896f9d0b1956912fc5df3 (metrics.json › evidence.benchmark.cells.k12_core) | Intel Xeon Gold 6148 (CPU_SKU:6148) (metrics.json › evidence.benchmark.cells.k12_core) |
| simulated, standard | 3:45:00 (metrics.json › evidence.benchmark.cells.sim5m_standard) | fa49b8f145ea4598a5f896f9d0b1956912fc5df3 (metrics.json › evidence.benchmark.cells.sim5m_standard) | Intel Xeon Gold 6148 (CPU_SKU:6148) (metrics.json › evidence.benchmark.cells.sim5m_standard) |
| simulated, core | 3:29:06 (metrics.json › evidence.benchmark.cells.sim5m_core) | fa49b8f145ea4598a5f896f9d0b1956912fc5df3 (metrics.json › evidence.benchmark.cells.sim5m_core) | Intel Xeon Gold 6148 (CPU_SKU:6148) (metrics.json › evidence.benchmark.cells.sim5m_core) |

Benchmark definition (ADR-M10-D): ADR-0072-benchmark-baseline-1-0-0 (metrics.json › evidence.benchmark), status accepted-by-khoa (metrics.json › evidence.benchmark).

## 9. Reproducibility

Freeze commit of the runbook: fa49b8f145ea4598a5f896f9d0b1956912fc5df3 (reproduce.md › evidence.runbook) (runbook status final (reproduce.md › evidence.runbook)); every arm, the compare and the metrics ran at it (`validation/pilot/reproduce.md`, Evidence table).
The runbook's check writes `verify-<commit>.json` beside the freeze runs; it is part of the review packet.
TE7 at the freeze commit: pass (metrics.json › te7.gippa).
The benchmark was not re-run at the freeze commit (M10-12 L12); its commit is in section 8.
Row-level outputs stay on FarmShare and are removed only by the runbook's guarded clean-up step.

## 10. Deviations from the spec

- DIA-NN 2.3.1 Academia instead of 2.7.0 (ADR-0051).
- Row-level inputs and outputs are not committed (Khoa's rule of 2026-10-23); aggregates and hashes are (2026-10-24).
- Protein groups are groups with at least one accepted precursor, with no protein-level FDR control (metric 1).
- Strain-specific genes are defined by peptide identity, not orthology (metric 3).
- Region tests are not compared, and the Jaccards are not restricted to genes detected in
  both arms (H3; a conservative stand-in, biased against H3).
- H2 is deferred to the full validation; metric 8 stands beside it without a verdict.
- The A3 branch when it is the fallback or absent, every deferral made at the freeze, and any
  accepted ADR that changes a pilot setting are stated in the notes.

## Notes

## Diagnosis: TE4

Cause class: search

TE4 failed on the rplL/rplJ pair and passed on atpA/atpD; the values, the range and the status are in section 5 and are unchanged. The range was not moved (R3).

The cause is not proven. iBAQ and the pre-registered Top3 secondary agree, and the observable-peptide counts do not point to a detection gap. The best-supported explanation is DIA-NN's choice of elution peak for the few peptides that carry the rplL signal; no chromatogram was inspected. The class names that reading; a true ratio below the expected one in this sample is not excluded.

Waiver: te4-rplL-rplJ-waiver-1-0-0 (docs/adr). It covers this pair only, for version 1.0.0; TE4 stays reported as a fail. The evidence is in that record.

## Read these with the numbers

- **Entrapment (metric 8).** The first attempt at the freeze commit failed its guard because of a bug in the tool: its six-frame index had no start cells, so peptides with an initiator M at annotated non-ATG starts were placed only through the Bakta-protein rule, not as exact, and the check against A1's exact route stopped the run. PR #522 fixed the tool. The fix changes which rule places these peptides, not the counts: the pre-fix tool, run without that check, gives the same N_T, N_E and FDP values. Precursor level: combined FDP 2.36 % (an estimated upper bound in expectation, not a bound on this sample); realised lower bound 0.55 % (if the genome is complete, is the sample's strain, and no K-12 material is present); 95 % upper confidence bound 2.65 % (only under equal chance and independent false matches); assessment `inconclusive`. Peptide level: combined FDP 2.73 %, assessment `not_assessed`. In plain terms, at precursor level the 1 % target is neither confirmed nor shown to fail. Known biases run both ways: near-homolog entrapment pushes the estimate high, and DIA-NN's use of protein information in A1's q-values pushes it low; neither is quantified.
- **R2.** A2's mapping FDR was not measured, under r2-mapping-fdr-waiver-1-0-0 (docs/adr). A2's accepted peptides are exempt from the decoy lists, so no value exists. `metrics.json` has no `mapping_fdr` key, and the threshold is unchanged.
- **TE4** stays `fail` (te4-rplL-rplJ-waiver-1-0-0; diagnosis above).
- **TE2** passes with 0 ALT sites in A2, so the pass is vacuous: it checked nothing.
- **TE7** compares the two GIPPA runs and passes on GIPPA's tables. The DIA-NN layer was not run, so reproducibility from raw data to tables is not shown.
- **H1** is not supported: its "mostly" clause fails, with 106 of the 417 A2-only genes in plasmid, prophage or BGC regions.
- **Runtime (ADR-0072).** The benchmark has 2 repeats, not 3, and readiness row R7 is waived for 1.0.0; the range of §13.2 was not met. Each "median" in section 8 is the mean of 2 runs. As ADR-0072 writes them: K-12 standard 3:40:32, K-12 core 3:01:13, simulated standard 3:44:59, simulated core 3:29:05. All four means fall on a half second: ADR-0072 truncates them, and the table in section 8 and `summary.md` round them to the nearest even second, which leaves K-12 standard at 3:40:32 and moves the other three up (3:01:14, 3:45:00, 3:29:06). Simulated standard is marked noisy (spread 16.0 %, above the 15 % threshold). With 2 repeats a spread is a weak estimate.
- **Novel ORFs in A1.** A1 searched a K-12 proteome, so these rows are annotation gaps, not new proteins. Its one called ORF (rank 1, class intergenic in A1's `tables/novel_orfs.tsv`) is the first 25 codons of prfB, K-12 release factor RF2, which has a programmed +1 frameshift at `CTT TGA C`. Its peptides lie at residues 1 and 11, and Bakta's prfB CDS starts downstream of the frameshift. It is not a novel protein. This identification comes from sequence and known biology, with no separate test.
- **Six-frame.** No freeze arm searched a six-frame database. A separate exploratory arm, A4, ran after the freeze and is not part of the gate. It called no novel ORFs, but no call was possible at its class sizes, so this is not a negative result (`docs/user-guide/makedb.md`).
- **Input genome.** The arms use the pass-3 wild-type genome, by Khoa's decision of 25 October.
- **UniProt release.** The release of the A1 database UP000000625 is not in this report's inputs or in A1's `provenance.json`. It is recorded in `$GIPPA_DB_DIR/prefetch_manifest.tsv` against the same sha256 as `metrics.json › inputs.k12_fasta`, and in `$GIPPA_DATA_DIR/uniprot/UP000000625.provenance.json`.

## Deviations at the freeze

- Entrapment: the rerun that counts ran from main (commit db3a95f1), not the pinned freeze commit. `src/` and `workflow/` are identical at the two commits; the tool and its job script differ (#522). The failed first attempt is a superseded row in the runbook's Evidence table.
- Genome and annotation: the arms use the pass-3 wild-type reference and its recorded Bakta files, not the committed chassis genome named in the M10-14 commands (README, "Pass-3 runs"). The recorded GBFF sits next to the JSON, so antiSMASH ran.
- Freeze window: a pinned detached worktree held the commit, not the main checkout (M10-14 question L4 took the other default; passes 2 and 3 used a pinned worktree too).
- A3 feasibility: the input is the pass-3 evaluation on the wild-type genome, not the committed verdict on the chassis genome. Both verdicts are pass.
- The Cascadia table is the pass-2 conversion. De novo search does not depend on the reference genome, and it was not re-run (index decision 4).
- TE7 compares the two GIPPA runs. The DIA-NN replicate layer was not re-run.
- Benchmark: all four cells ran at the freeze commit, so section 9's sentence that it was not re-run there does not apply. Repeat 3 was cancelled on Khoa's decision (ADR-0072). `benchmark_summary.py --rounds 2` is main's script (#523), because the script at the freeze commit has no such option. `metrics.py` reads ADR-0072 and the benchmark's `baseline.json` from the report branch, since neither the record nor its `standard_run` entry is at the freeze commit, so S14's relative paths do not re-run from the pinned worktree as written.
- The clean-up step S15 was not run.
- `metrics.json` has no `mapping_fdr` key: `metrics.py` at this commit does not compute one (see R2 above).
