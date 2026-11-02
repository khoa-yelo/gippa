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
