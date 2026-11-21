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

## Addendum: 1.1.1 detection recheck

This section covers decoy-pairing-independent-of-acceptance (docs/adr), PR #607, released in 1.1.1. In 1.1.0, an unpaired decoy competed in the picked gene FDR only when its own unit had a q-passing target. With the fix, every q-passing decoy competes. Every pilot `decoy_picked` arm was recomputed with the 1.1.0 code and with the fix. Then the gene table, detectability, context, function and compare steps were rerun on the freeze outputs, with only `detection.py` swapped. The control run reproduces the recorded tables byte for byte. The Slurm job IDs, the scripts and the full before and after values are in `validation/pilot/recheck_1_1_1/README.md`. The superseding values are in `validation/pilot/metrics/metrics_1_1_1.json`, under the key paths of `metrics.json`.

The values above are the 1.0.0 record and stay as they are. For a 1.1.1 reading, the following are superseded: metric 3 rows all, chromosome and k12_shared in both arms, A1 bgc and A2 strain_specific; metric 4 rho, its interval, shared genes and bias; the metric 5 enriched-category Jaccard; and in H1 and H3 every count, share, interval and p-value listed below. No threshold changed. No verdict changed: H1 stays not_supported (its "mostly" clause still fails), H3 stays supported, and TE3, TE4 and the plasmid and prophage coverage are unchanged.

| arm | detected genes | maximum gene q | decoys paired | own units | decoy wins |
|---|---|---|---|---|---|
| A1 (pass 3 and freeze) | 28 fewer | rises, stays below 0.01 | 189 | 0 → 33 | 4 → 37 |
| A2 (pass 3, freeze and replicate) | 6 fewer | rises, stays below 0.01 | 323 | 0 → 2 | 49 → 51 |
| A4 (six-frame, exploratory) | unchanged | rises, stays below 0.01 | 168 | 0 → 4 | 15 → 19 |

A1 searched K-12 accessions, so 38 of its decoys have no gene pair, and 33 of them now compete as their own units. Pass 1 changes the same way (A1 28 fewer genes, A2 5 fewer).

- **H1** (1.0.0 → 1.1.1): A2-only genes b 417 → 434, A1-only genes c 33 → 28; TU-level units 275 against 31 (1 mixed) → 290 against 24 (3 mixed); the plasmid, prophage or BGC genes among the A2-only genes stay at 106, now of 434 rather than 417. The "Read these with the numbers" line on H1 reads 106 of the 434 A2-only genes in 1.1.1.
- **H3** enriched categories (n_a / n_b / n_shared / union): 437 / 285 / 260 / 462 → 405 / 292 / 258 / 439; the criterion still holds, and so does rho.
- **Metric 3**: A1 bgc 78 → 76 genes, A2 strain_specific 317 → 316 genes; plasmid, per-plasmid, prophage, A1 strain_specific, A2 bgc, no_peptides and pseudogene are unchanged.
