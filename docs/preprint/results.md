# Results

<!-- Draft results section for the GIPPA 1.1.0 preprint (revision R3, after audits R1–R3). Every
number is copied from the file named in the `src` comment after its paragraph or table. Paths are
relative to the repository root. Abbreviation used in the comments: CAL =
local/logs/calibration/full-4a6fa998300c-20261109-160014 (the release-lineage calibration_full run).
Pilot rows are aggregates only (no peptide sequence, no DIA-NN row, no per-gene value except the
named TE4 genes of the specification). -->

## Provenance of the results

Each subsection names the code state it used.
Calibration on simulated data and the attack nulls ran at commit `4a6fa998`. The release code is
`e43683c7`; the tag v1.1.0 points at `a2a95c76`, which adds only the readiness record. Between
`4a6fa998` and `e43683c7`, only the version string in `src/gippa/__init__.py` changed in `src/`,
`workflow/` and `tests/calibration/`. The fast calibration tier, the integration and nf-test runs,
the K-12 weekly run (TE6, TE8) and the housekeeping rows of the readiness check ran on `e43683c7`
itself.

The Nissle 1917 pilot, the true-example values TE1–TE5 and TE7, and the standard-run benchmark were
measured on the 1.0.0 report-freeze commit `fa49b8f1` (30 October 2026) and were not re-run on 1.1.0
code. Between `fa49b8f1` and `e43683c7`, 124 files in `src/` and `workflow/` changed (+6,484/−1,090
lines), mostly plotting (24 in `viz/`), report (23), table schemas (10) and workflow tests (13).
Statistical code also changed:
- mapping: 3 files, including `loci.py`;
- operons: 5 files, including `discordance.py`, `flags.py` and `missed.py`;
- novel ORFs: 5 files;
- variants: 4 files, including putative alternative starts, which are capped at `candidate`;
- gene table: 2 files;
- determinants: 2 files;
- pipeline steps (`steps/`, 6 files, including `steps/mapping.py` and `steps/detectability.py`),
  input/output (`io/`, 7 files, including `io/run_view.py`, which `gippa compare` reads),
  preflight (3) and `model/genes.py`;
- output columns were renamed for 1.1.

No file in `context/`, `stats/`, `function/` or `compare/` changed, although `compare` reads its
inputs through the changed `io/run_view.py`. Apart from the detection recomputation described
next, the pilot numbers describe 1.0.0-era code; whether 1.1.0 reproduces them (TE1, Tables R5–R7,
operon flags, ALT and novel-ORF calls) was not checked.

GIPPA 1.1.1 corrects the picked gene FDR (ADR-0078, PR #607; Methods Section 6). The pilot values
that depend on gene detection (Table R6, Section 4.3, and the H1 and H3 rows of Table R7) were
therefore recomputed with the 1.1.1 detection code on the frozen pilot inputs. The gene table,
detectability, context-statistics and function steps of the A1 and A2 freeze runs, and
`gippa compare`, were re-run with the freeze code `fa49b8f1`, in which only
`src/gippa/genetable/detection.py` was replaced by its 1.1.1 version. A control with the unchanged
freeze code (Slurm job 2001992) reproduced the recorded tables byte for byte; the reported values
come from the run with the fix (Slurm job 2002006). The fix removed 28 detections from A1
(2,673 to 2,645) and 6 from A2 (3,057 to 3,051) and added none. The outcomes of TE3, TE4, the
plasmid and prophage rows of Table R6, and H1 and H3 did not change. All other pilot values are
from the freeze runs.

<!-- src: local/scratch/a1-recheck/{rerun.sbatch (control, Slurm 2001992), rerun_fix.sbatch (fix, Slurm 2002006), rerun_steps.py, derived.py}; logs/gippa-a1recheck-2001992.out, logs/gippa-a1recheck-fix-2002006.out; diff -rq of the freeze src/gippa against local/scratch/a1-recheck/hybrid/src/gippa (only genetable/detection.py differs; 1.1.1 detection.py from PR #607, merged as 92fbf3af); cmp of runs/ctl/{a1,a2}/out/tables/*.tsv and runs/ctl/cmp/tables/*.tsv against local/runs/pilot/freeze/fa49b8f (0 differ); derived_ctl.json, derived_fix.json › m3.*.all, h1, h3; detected/{v1.1.0,fix7d4c4bfd}/freeze__fa49b8f__out_{a1,a2}.tsv (genes at q ≤ 0.01: 28 and 6 lost, 0 gained; no ribosomal, rpl or atp gene among them) -->

The 1.1.0 readiness check evaluates rows R1, R5 (TE1–TE5), R6 (TE7), R7, R8 and R9 directly on the
1.0.0-era pilot metrics and benchmark file. No decision record covers reusing the passing rows.
ADR-0076 covers only the three waived rows (R2, R5 for TE4 only, and R7), because "Nothing these
rows measure has changed since 1.0.0"; that refers to the data, which were not regenerated, not to
the code.

| Label | Commit(s) | `src/`+`workflow/` difference to `e43683c7` | Results that use it |
|---|---|---|---|
| release code | `e43683c7` | none | fast calibration tier, integration, nf-test, TE6, TE8, readiness H1–H6 |
| release-equivalent | `4a6fa998`, `20283834` | version string only | calibration_full and attacks (`4a6fa998`); weekly Nissle TE2 and K-12 distance-model fallback (`20283834`) |
| 1.1.0 development | `39f9f3a2`; `e6f4aa3c` | 2 files (+7/−3); 20 files (+780/−128) | example routes A–E (`39f9f3a2`); 9 November nightly runs (`e6f4aa3c`) |
| 1.0.0-era | `fa49b8f1` | 124 files (+6,484/−1,090) | Nissle 1917 pilot (A1, A2, A3, compare); TE1–TE5, TE7; standard-run benchmark |
| 1.0.0-era with 1.1.1 detection | `fa49b8f1` with `genetable/detection.py` of `92fbf3af` | as 1.0.0-era, with the 1.1.1 `detection.py` | pilot detection, Table R6, Section 4.3, H1 and H3 (Slurm 2002006; control 2001992) |
| older | named per row | – | UniOP AUROC and cut-offs (ADR-0022, Slurm 1879748); TE2 in-silico proxy (2026-10-23) |

<!-- src: git for-each-ref refs/tags/v1.1.0 (a2a95c76); git diff --stat e43683c7 a2a95c76 (READINESS.md, readiness.json only); validation/release/1.1.0/readiness.json › commit, rows[R1,R5,R6,R7,R8,R9].evidence; docs/adr/ADR-0076-waivers-1-1.md:1-40 (title, Context line 17, Scope);
validation/pilot/metrics/metrics.json › gippa_commit; benchmarks/baseline.json › standard_run.wheat.k12_standard.commit; CAL/task-0.json › commit;
git diff --stat / --shortstat / --name-only fa49b8f1, 4a6fa998, 20283834, 39f9f3a2, e6f4aa3c .. e43683c7 -- src workflow (file counts per module; commit c6b05217 "one name per quantity across tables (1.1 renames)"; 98ff402b putative starts) -->

## 1. Calibration on simulated data

### 1.1 Harness and run

Each statistical method was calibrated on data from gippasim, the project's own simulator, which
the project's implementer role writes and its test-author role tests. Every calibration result in this section is
therefore conditional on gippasim's generative model. The full validation was not done: the K-12
dataset PXD075568, FDRBench entrapment and the K-12 co-abundance AUROC remain open, and the sealed
EDL933 holdout was not used, so no holdout claim is made. The criteria are fixed in
`tests/calibration/tolerances.yaml` and applied by `tests/calibration/_criteria.py`. Three error-rate
criteria are one-sided (the observed rate must not exceed the bound), and one uniformity test is
two-sided:
- **level:** the share of p ≤ α is no more than α + 2·√(α(1−α)/n), where n = R replicates, or R·m
  for rows that pool m p-values per replicate;
- **uniformity:** a two-sided Kolmogorov–Smirnov test against Uniform(0, 1) with p > 0.001; the
  check refuses to run when more than 5 % of the p-values are tied (discrete p);
- **FDR:** mean FDP no more than q + 2·SD(FDP)/√R;
- **FWER:** the share of replicates with at least one false discovery is no more than
  α + 2·√(α(1−α)/R), with α = 0.05 by default.

Fixed minimums or bands also apply:
- F1 ≥ 0.90 for `high` variant calls and novel ORFs;
- planted-peptide mapping recall ≥ 0.999;
- detectability slope within [0.8, 1.2];
- de novo precision |bias| ≤ 0.05;
- determinants ΔR² |bias| ≤ 0.05.

Some coverage floors are coded in the test files, not in `tolerances.yaml`. Examples are the
coverage of the free-entrapment upper bound (≥ 0.95 − k·√(0.95·0.05/R), arms A and B), the de novo
precision interval (0.95 − k·SE) and the class error-rate intervals (0.95 − 3 SE, looser than the k = 2 used elsewhere).

Conservative error is allowed; anti-conservative error beyond the bound fails. The nominal levels
were 0.01, 0.05 and 0.10.

Checks run through the harness's confirmation step use two blocks (ADR-0015). A check that exceeds
its bound on replicates 0 to R−1 is re-run on a fresh block (R to 2R−1), and it fails only if both
blocks exceed the bound. An excess in one block only is logged as `unconfirmed` and does not count as
a failure. The decision record estimates that, for about 46 checks per nightly run, this lowers
false alarms from about 0.75 to about 0.04 per run, and lowers power against a 1.5-fold FDR
inflation at α = 0.05 from 0.92 to 0.85. The full tier judged about 1,375 harness records one by
one, with no adjustment across suites. At the same per-check rate (about 0.023 per block), about 0.7
confirmed false alarms would be expected if every method were exactly at level (our calculation).

**What a pass shows.** A pass means that the observed rate did not exceed the nominal level plus two
Monte Carlo standard errors. It rules out larger excess error, not excess within that margin. For a level check at α = 0.05, the bound is
about 1.9α at R = 100, 1.5α at R = 300 and 1.3α at R = 1,000. At α = 0.01 and R = 300 it is about
2.1α, and for an FWER check at R = 500 it is about 1.4α. Where a gated point estimate exceeds the
nominal level, we report it as within tolerance, not as at level.

The full tier ran 15 suites as 15 Slurm jobs (1966680–1966694, recorded as "array 1966680") on commit `4a6fa998` on 9–10 November 2026.
The islands task was re-run as 1968544 after the first attempt ended with exit 143. Replicates per
suite were:
- 100: detectability, determinants, variants;
- 300: placement, gene_fdr, missed_operons, gsea, novel_orfs, denovo_precision;
- 500: islands;
- 1,000: regions, autocorrelation, discordance, ora, compare.

<!-- src: AGENTS.md (Conventions: gippasim is written by the implementer, tested by the test-author); docs/preprint/outline.md:139-145 (full validation not done; holdout sealed); tests/calibration/tolerances.yaml:1-19; tests/calibration/_criteria.py:219-300 (level_bound, check_level_pooled with R·m, FDR with sd/√R); tests/calibration/_criteria.py:107,243-260 (check_uniform: sps.kstest two-sided, MAX_TIE_SHARE 0.05), :312 (check_fwer uses level_bound(alpha, R)); tests/calibration/attacks/test_attack_free_entrapment.py:256-257,307-313; tests/calibration/test_denovo_precision_calibration.py:109; tests/calibration/attacks/test_attack_class_error_rates.py:187; CAL/task-N.json (commit 4a6fa998, exit 0);
docs/adr/ADR-0015-calibration-two-block-confirmation.md:8-20 (Decision; Evidence: 0.75 → 0.04; power 0.92 → 0.85); bounds computed from level_bound with k = 2;
CAL/failed-1966685/task-5.json; local/logs/nightly.tsv:124 -->

All 15 tasks exited 0. The junit roll-up counted 555 tests: 459 passed, 90 skipped, 6 xfailed,
0 failed, 0 errors and 0 xpassed. Most skips are shard and confirmation placeholders whose verdict is
given by the merge step.

The fast tier (calibration_nightly) passed on the release code `e43683c7` (28,937 passed, 148
skipped, 17 xfailed, 3 xpassed; line coverage 97 %, minimum 90 %). Readiness row R3 records a pass
without a waiver: it reuses the full tier from `4a6fa998` because the intervening diff touches only
non-statistical paths, and it checks only for unexpected failures (Sections 1.3 and 7). For 1.0.0
the row was waived (ADR-0073), because the fast tier had not run on the release commit.

<!-- src: local/logs/nightly.tsv:124, 137 (1972754_2, e43683c7); validation/release/1.1.0/READINESS.md:11 (R3, "reused from 4a6fa99 (reusable: diff ⊆ NON_STATS_PATHS)");
validation/release/1.0.0/READINESS.md (R3 waived, ADR-0073) -->

### 1.2 Results by method family

Table R1 summarises each suite. Table R2 gives representative observed error rates against their
bounds at α = 0.05, with m for pooled rows.

**Placement FDR (mapping).** The placement suite produced 366 gated harness records, of which 365
passed. A further 6 read-through-balance records were vacuous, because no read-through target
occurred. At q = 0.05 the observed mean FDP was 0.000158 in the stressed database arm (R = 300) and
0.0355 in the stressed read-through arm with shuffled decoys.

The one failing record was the contaminant attack `flat_c200` at q = 0.01, with FDP 0.01496 against a
bound of 0.01019. It is a strict xfail whose marker reads "contaminant weight pending (M3-25
decision 10 left it open)". The same attack passed at q = 0.05 and 0.10.

The methods note classes the 100-kb genome-level lists as smoke checks: most replicates contained no
false peptide, so their gates pass trivially. The records do not mark them, so an unknown number of
the 365 passes carry no FDR evidence. The FDR evidence comes from the list-level suites and the stressed stratum-0 arm. Gate M3
was recorded as not met at the 1.0.0 freeze and shipped in 1.0.0 as a known limit; no later gate
record exists.

<!-- src: CAL/lanes/placement-group4/placement.jsonl:20; CAL/lanes/placement-group1/placement.jsonl:2; CAL/placement.merged.jsonl:43-45;
CAL/lanes/placement-group{3,4,5,6,9}/placement_readthrough_balance.jsonl (6 streams, status vacuous, rt_share 0.0); docs/methods/placement-fdr.md:776-787; PROGRESS.md:155-165 -->

**Gene-level FDR (picked, MAYU).** All 132 gene_fdr records passed. At q = 0.05 with π = 0.05, the
observed mean FDP was 0.0486 for the picked FDR (bound 0.0506) and 0.0493 for MAYU with PEP input
(bound 0.0505).

The methods note records two limits:
- the levels 0.05 and 0.10 are weak and vacuous in gippasim, because peptides are capped at
  q ≤ 0.01 and the loosest reachable gene FDR is about 0.053;
- at q = 0.05 the FDP is bimodal, with a mean above q (0.053 ± 0.003; p90 0.135), when 50 % of
  target q-values are missing.

At q = 0.01 the observed FDP was 0.00941 (picked; bound 0.01031) and 0.00962 (MAYU with PEP; bound
0.01022). The simulation places false hits on genes in proportion to their observable peptides,
which is MAYU's own model. The unattacked suite therefore partly assumes what it tests; departures
from that model are covered by the attack files. The `two_peptide` tier carries no FDR estimate; in
the same simulation its gene-level FDP was 0.054 at q = 0.01 (bound 0.0105; recorded, not gated).

<!-- src: CAL/lanes/gene_fdr-group1/gene_fdr.jsonl:44,47 and the q 0.01 rows of the same streams; CAL/lanes/gene_fdr-group1/gene_fdr_two_peptide.jsonl:1 (0.05385, bound 0.01047, ok false); docs/methods/detection-tiers.md:207-213, 296-301 -->

**Detectability and detection limit.** The cross-fitted detectability slope had mean 0.981 over 100
replicates, inside the band [0.8, 1.2]. The detection-limit calibration-in-the-large (CITL) mean was
−0.0078, and single replicates ranged from −0.293 to 0.296. `tolerances.yaml` defines a tolerance
|CITL| ≤ 0.2 (ADR-0030), but no calibration test applies it; the test records the value ungated.

<!-- src: CAL/lanes/detectability-group0/detectability.jsonl:7,12 (citl list: min −0.2933, max 0.2962); tests/calibration/tolerances.yaml:6; tests/calibration/test_detection_limit_citl.py:1-10; docs/adr/ADR-0030*.md:8-22 -->

**Region tests (circular shifts).** All 146 region records passed. The pooled Hodges–Lehmann shift
test on the chromosome at ρ = 0 rejected at 0.0505 at α = 0.05. The bound was 0.0569, from
R = 1,000 replicates × m = 4 pooled tests. The within-type FDR for prophages was 0.063 at
q = 0.05 (bound 0.0654). This point estimate is 26 % above q and within tolerance; it is not at level.

The TU-block bootstrap CIs used for region effects are checked against an expected small-sample
percentile coverage c(n), not against 0.95. For example, in the two-group HL attack the target
non-coverage was 0.0873 and the observed non-coverage 0.081, about 0.92 coverage; with few large TUs
the target was 0.0638 and the observed value 0.051. TU-block percentile intervals with few TUs are
therefore expected to under-cover.

<!-- src: CAL/lanes/regions-group2/regions.jsonl:21 (m 4); CAL/lanes/regions-group0/regions.jsonl:8; CAL/lanes/regions-group*/ (regions.attack_bootstrap.hl_two_group level 0.0873 observed 0.081; few_large_tus level 0.0638 observed 0.051);
tests/calibration/attacks/test_attack_bootstrap.py:110-132 -->

**Autocorrelation.** With exact transcription units (TUs), the TU-block Moran's I at lag 1 rejected
at 0.055 (bound 0.0638, R = 1,000), and the BH family had an any-rejection rate (FWER) of 0.008.

Real runs never have exact TUs, and the statistic was anti-conservative under every TU-error model
tested:
- **UniOP-like TU errors** (recorded, not gated): 0.071 at α = 0.05, above the bound of 0.0638;
- **UniOP's K-12 error rates** (recorded, not gated): 0.432 / 0.659 / 0.758 at α = 0.01 / 0.05 / 0.10;
- **distance-model operon fallback** (ADR-0061, strict xfail): 0.149 / 0.354 / 0.476 at
  α = 0.01 / 0.05 / 0.10, against bounds of 0.0163 / 0.0638 / 0.1190.

Of 174 records, 163 passed, 4 were unconfirmed and 7 failed. The failures are 3 strict-xfail records
for the distance-model arm and 4 recorded-only records (UniOP-like at α = 0.05, and UniOP K-12 rates
at all three levels). The statistic is within tolerance only for exact TUs, where its point
estimates were 0.011 / 0.055 / 0.113 at α = 0.01 / 0.05 / 0.10 (bounds 0.0163 / 0.0638 / 0.1190), all
above nominal, with a pooled excess for large heavy-tailed TUs (below).

Four autocorrelation streams were `unconfirmed` (over the bound in one block of 1,000 replicates,
not both). An exact-TU attack with large, heavy-tailed TUs gave 0.065 and 0.057 at α = 0.05; pooled
(0.061), it exceeds the bound for 2,000 replicates (0.0597). A second is
the UniOP-like arm at α = 0.01 (0.017 and 0.016; pooled 0.0165 against 0.0145 at 2,000 replicates).
The other two pooled within that bound: the UniOP-like join-count lag 20 at α = 0.10 (0.120 and
0.105) and the join-count residual engine at α = 0.01 (0.017 and 0.008).

Because real runs never have exact TUs, every TU-block Moran's I row ships with the flag
`tu_block_tu_error` on every run, and with `tu_source_distance_model` under the fallback. Lag
statistics have no confidence interval (flag `ci_unestimated`), although specification §8.2 asks for
an effect size with a 95 % CI beside every p or q.

<!-- src: CAL/lanes/autocorrelation-group0/autocorrelation.jsonl:76-78, 44 (uniop_like, 0.071, status fail), 47, 49, 51; CAL/lanes/autocorrelation-group2/autocorrelation.jsonl:57, 23, 25, 27;
CAL/lanes/autocorrelation-group{0,2}/autocorrelation.jsonl (status unconfirmed, blocks 1 and 2; pooled bounds computed with level_bound at n = 2,000); docs/methods/autocorrelation.md:141-147, 186;
tests/calibration/attacks/test_attack_operon_distance_fallback.py:54-55 ("recorded, not gated") -->

**Abundance islands.** All 69 gated island records passed. At α = 0.05 the FWER was 0.022, 0.020 and
0.024 at AR(1) ρ = 0, 0.5 and 0.9, and 0.056 for the chromosome-plus-plasmid step null (bound 0.0695,
R = 500).

The FWER claim is limited to a stated dependence regime (ADR-0021 Amendment 1). The full tier at
`4a6fa998` (islands rerun 1968544) also recorded, without gating, the two breaches outside that regime: ρ = 0.9 with short
resolved lag gave 0.096 at 1,000 genes and 0.130 at 1,800 genes, against the bound of 0.0695
(R = 500).

This domain was narrowed after the calibration had exposed failures: Amendment 1 (20 October 2026)
states that the earlier inference "was wrong". No output flags whether a run lies inside the domain.
The report shows the resolved block length ℓ̂, the block length ℓ′ and the lag-1 Moran's I per
replicon as diagnostics, but no rule turns them into an inside/outside verdict. The low
direction of the island test is a structural zero at normal detection rates: a null result in the low
direction is not evidence of no silenced island.

<!-- src: docs/adr/ADR-0021-island-oriter-block-length.md:56-61; docs/methods/islands.md:23, 194, 209-227;
CAL/lanes/islands-group0/islands.jsonl:32,35,38,47; CAL/lanes/islands-group0/islands_recorded.jsonl:1-2 (job 1968544); tests/calibration/test_island_fwer.py:179;
docs/methods/islands.md:166-172, 323-326; docs/adr/ADR-0021-island-oriter-block-length.md (Amendment 1); src/gippa/report/content.py:238-249 (island_diagnostics table) -->

**Operons.** The discordance gate passed at R = 1,000:
- pooled per-pair level 0.0499 at α = 0.05 (bound 0.0510, m = 200);
- BH FDR 0.034 (bound 0.0615).

Its three failing records come from two strict-xfail tests:
- the step-shaped σ(abundance) attack D-1b-step at α = 0.05 and 0.10 (second block 0.0517 and
  0.1031, against bounds of 0.0510 and 0.1013), which lies outside the narrowed claim of ADR-0055
  Amendment 1 (smooth, monotone σ);
- p-value uniformity under the distance-model fallback (KS p 2.5e-06 in both blocks, against the
  minimum of 0.001; ADR-0061). The departure is conservative: the pooled level was 0.0078 at
  α = 0.05.

The claim was narrowed to smooth, monotone σ(abundance) by Amendment 1 (24 October 2026), after the
step attack failed. No output tells a user whether their data have smooth σ.

For missed operons, the convergent-only decoy arm made calls and stayed at level (FDP 0.00065 at
q = 0.05), and the engine's planted arm had FDP 0.0036. The shipped configuration (convergent plus
sandwich decoys) made no call in any of the 300 power-arm replicates (vacuity floor 1/600 calls per
replicate). Of the suite's 48 records, 26 are `vacuous`: 18 from the shipped decoys (6 streams × 3
levels), and 8 from other arms, 7 of them at q = 0.01, including the convergent-only arm. For the
shipped decoys the result is a vacuous pass, not evidence of power or of level; the convergent-only
arm gives level evidence at q ≥ 0.05 only.

<!-- src: CAL/lanes/discordance-group1/discordance.jsonl:17 (m 200), 21; CAL/lanes/discordance-group*/discordance.jsonl (attack_d1b.step.loglinear blocks 1-2; fallback_a.operon_null uniform 2.458e-06 / 2.543e-06);
docs/adr/ADR-0055-operon-discordance-abundance-noise.md:113-118; CAL/lanes/missed_operons-group0/missed_operons.jsonl:35 (all 300 values n_called 0; vacuous floor n_called 0.00167; tests/calibration/_criteria.py:96, 960-962), 38, 44, and status "vacuous" records (26);
CAL/lanes/discordance-group*/discordance.jsonl (fallback_a.operon_null level 0.05 observed 0.00784);
docs/adr/ADR-0055-operon-discordance-abundance-noise.md (Amendment 1); PROGRESS.md:188-197 -->

**Function.** Over-representation analysis (ORA) is shipped uncalibrated: every ORA row carries the
flag `ora_uncalibrated`. In the full tier the fitted ORA null was within tolerance (0.0235 over, 0.0303
under, bound 0.0556 with m = 6). However, 10 of 72 records failed. For example, the large-TU
distance-model arm rejected at 0.0841 against 0.0544 at α = 0.05 (m = 10). Four of the ten failures
are UniOP-like TU-error arms (for example 0.062 against 0.0544 at α = 0.05, m = 10). One UniOP-like
stream was `unconfirmed` (0.0563 and 0.0525; pooled 0.0544 against 0.0531 at 2R). The 10 ORA tests
pass although 10 records fail, because the flag makes those records recorded-only.

The pipeline runs GSEA on R fgsea 1.36.2, but the GSEA function's default engine, used by most
GSEA calibration tests, is the Python port `fgsea_simple`. Streams on the port are evidence about
the method, not about the shipped engine. GSEA was within tolerance under TU exchangeability:
- on R fgsea, the unpooled null 0.0167 / 0.0567 / 0.127 at α = 0.01 / 0.05 / 0.10 (bounds 0.0215 /
  0.0752 / 0.1346), above nominal at every level, and the pooled null over 100 tested sets per replicate from six
  collections, 0.0415 (bound 0.0525, m = 100);
- on the Python port engine, the TU null 0.0333 (bound 0.0752).

Its two strict xfails (Python port engine) are the restricted-claim cases of ADR-0025 and ADR-0026:
detectability-confounded sets and AR(1)-clustered sets.

The suite's "111 pass" includes a detection-confounding attack on rows the detection screen did
*not* flag, gated at 1.5α, not α (statistician decision #174: the flag is a screen). This screen and the
genomic-cluster checks ran on the Python port engine and have not been re-measured on R. One arm
(`null_b2`) rejected at 0.0833 at α = 0.05 and 0.137 at α = 0.10: within the 1.5α bounds (0.105,
0.191) but above the α bounds at R = 300 (0.0752, 0.1346). No α-level claim is made for unflagged
rows under detection confounding. Two further limits are recorded in the methods note:
- BH spillover from flagged rows has not been re-measured on the R engine, so no FWER or FDR claim is
  made for unflagged rows in a collection that contains flagged rows;
- NES has no confidence interval (ADR-0027).

<!-- src: CAL/lanes/ora-group1/ora.jsonl:67,70 (m 6), 42 (m 10), 26-28, 55, 57, 59 (uniop_like; 27-28 unconfirmed; pooled bound level_bound(0.05, 20,000)); CAL/lanes/gsea-group0/gsea.jsonl:45 (gsea.m815.null.r.full, m 100; tests/calibration/test_function_calibration.py:24-25, 52 M_GSEA = 100 over six collections); CAL/lanes/gsea-group1/gsea.jsonl:45-47 (gsea.r.null), 52; src/gippa/function/gsea.py:5-7, 545 (default engine SIMPLE_ENGINE); src/gippa/steps/function.py:17-19 (pipeline uses R fgsea); tests/calibration/test_gsea_r_engine_calibration.py (gsea.r.null on R); tests/calibration/attacks/test_attack_gsea_flags.py, test_attack_gsea_structure.py, test_gsea_tu_calibration.py (no engine argument: port);
CAL/lanes/gsea-group*/ (gsea.attack.flags.detection.null_b2: level 0.075 observed 0.08333 bound 0.1054; level 0.15 observed 0.13667 bound 0.1912); tests/calibration/attacks/test_attack_gsea_flags.py:366-386 (SCREEN = 1.5);
docs/methods/function-gsea.md:54-60;
src/gippa/function/ora.py:52; docs/methods/function.md:316,561,645; CAL/task-9.*.junit.xml -->

**Determinants.** The gated determinants records passed in the full tier, within tolerance: null
rejection rates ranged from 0.00 (LightGBM at α = 0.01) to 0.12 (LightGBM at α = 0.10, bound 0.16),
and at α = 0.05 from 0.02 to 0.07 (bound 0.0936, R = 100). The
permutation p-value nevertheless ships flagged `null_uncalibrated` (ADR-0034), because it fails the
dependence and detection-selection attacks.

In the release full tier, those attack arms rejected at α = 0.05 (R = 100; flagged, so recorded and
not gated) at the following rates:

| Model | AR(1) | AR(1) + positional feature | detection selection |
|---|---|---|---|
| LightGBM | 0.06 | 0.93 | 0.26 |
| elastic net | 0.18 | 0.42 | 0.83 |

A flagged TU attack also recorded univariate-CI non-coverage of 0.07 for AUROC (bound 0.0695) and
0.20 / 0.218 for the gene-level ρ / AUROC intervals (R = 500). Under the distance-model fallback,
the true-TU cross-check recorded non-coverage of 0.074 (ρ) and 0.072 (AUROC) against 0.0695
(flagged). The ΔR² interval is not calibrated: under the null it missed 0 in 0.388 of LightGBM
replicates (elastic net 0.02; bound 0.0695; recorded, not gated). The full-tier pass is therefore
not a level claim for real proteomes.

<!-- src: CAL/determinants.merged.jsonl:1-21 (min 0.0 line 1; max 0.12 line 3; α 0.05 max 0.07 line 17); CAL/lanes/determinants-group0/determinants_m808.jsonl:1 (null_delta_r2 miss_rate 0.388 / 0.02); CAL/lanes/determinants-group1/operon_fallback_arm_b_determinants.jsonl:1 (rho_true_tu_x 0.074, auroc_true_tu_x 0.072, flagged); CAL/lanes/determinants-group1/determinants_attack_tu.jsonl (univariate_ci_level, flagged); CAL/lanes/determinants-group{0,1,2}/determinants_attack_null.jsonl (arm, model, rates["0.05"], flagged true);
docs/adr/ADR-0034-determinants-null-uncalibrated.md (Context, Decision) -->

**Novel ORFs and variants.** Class-level FDR for novel ORFs was within tolerance at q = 0.05 when false ORFs
were generated as in the gate (picked 0.0453, de novo 0.0498, R = 300). In the de novo regime it is
anti-conservative when a false ORF carries several wrong peptides: the multi-peptide attack gave FDP
0.0788 at q = 0.05 (Table R3b). This follows from the shipped rule (ORF q = minimum peptide q).
ADR-0050 kept that rule as a documented limitation for exact systematic multi-peptide errors
(G4-k3). No ADR names the chance-hit multi-peptide attack: ADR-0050 defers the k-chances excess it
measures to the backlog (decision 3), and the attack's marker cites M7-11 decision 3.
The pooled novel-ORF F1 was 0.938 (SE 0.0014), against the minimum of 0.90.

Nine novel-ORF records failed, all strict xfails:
- the multi-peptide false-ORF attack at all three levels (marker "M7-11 decision 3: ORF q = min
  peptide q");
- the GC-poor intergenic attack at q = 0.01;
- five cells with exact systematic de novo errors (ADR-0050), with FDP up to 0.0749 at q = 0.01.

For variants, all three gated results are strict xfails:
- the pooled F1 of `high` ALT calls was 0.874 (SE 0.0018), against 0.90;
- the random isobaric-swap attack gave FWER 0.97 against a bound of 0.0936;
- the systematic isobaric-swap attack gave FWER 1.0 against the same bound.

These records use the scope and criterion from before ADR-0049, which, after G1 failed, redefined
the scope as unflagged `high` calls against identifiable truth. Under that scope:
- truth that is mass-equivalent to the reference leaves the recall denominator (26 % at 100 kb,
  24 % at 1 Mb);
- 37 % of `high` substitution calls are flagged and are not scored.

Under that scope, F1 was 0.932 [0.924, 0.941] at 1 Mb (scope c, code of 23 October), and 0.933
with the decided flags (scope d, no CI at 1 Mb). No full-tier record exists under that scope. Gate
M7 was recorded as not met at the 1.0.0 freeze and shipped in 1.0.0 as a known limit; no later gate
record exists.

<!-- src: CAL/novel_orfs.merged.jsonl:38,41,92,22-24,43,117,120-121,123-124; CAL/variants.merged.jsonl:1-3; PROGRESS.md:199-212 ("pre-ADR scope", "old criterion");
docs/adr/ADR-0049-m7-variant-precision.md:50-64 (Evidence table, scopes c and d; 26 %, 24 %, 37 %); docs/adr/ADR-0050-denovo-orf-q-exact-errors.md:1-10 (Context: G4-k3; Decision option 3); tests/calibration/attacks/test_attack_novel_class_fdr.py:25-34, 288-289 (de novo regime); PROGRESS.md:208; docs/methods/novel-orfs.md:125-140 -->

**Compare, entrapment and pilot-hypothesis attacks.**
- **ρ confidence interval.** The TU-block bootstrap CI for Spearman ρ in `gippa compare` covered at
  0.952 (non-coverage 0.048, bound 0.0638, R = 1,000; mean width 0.201). A gene-level bootstrap,
  recorded for contrast, covered at 0.787.
- **H1 test.** The TU-level H1 test held level under a uniform null (0.0449, bound 0.0544,
  R = 10,000). The gene-level p of H1, recorded only, was anti-conservative: 0.0948 / 0.1737 / 0.227
  at α = 0.01 / 0.05 / 0.10.
- **H4 test.** The H4 false-confirmation FWER was 0.047 (bound 0.0638).
- **Free entrapment.** The free-entrapment ucb95 covered at 0.9736 in arm A (bound 0.9469,
  R = 20,000). The recorded arms D and E, in which false matches are clustered (the regime where
  ucb95 is stated to be invalid), covered at only 0.754 and 0.719. In arm D the combined estimate
  averaged 0.009 against a mean true FDP of 0.016.
- **De novo precision.** This suite is gated by assertion and writes no harness record. The gated
  claim concerns realised precision: the per-bin bias against it was 0.0015 to 0.0094 across
  10 bins (tolerance 0.05). Against the simulator's generative truth table, which is not gated,
  the bias reached −0.066, beyond the 0.05 tolerance of §8.4.

<!-- src: CAL/lanes/compare-group0/compare.jsonl:41,27,38,13; CAL/lanes/compare-group0/compare_coverage.jsonl:1-2; CAL/lanes/compare-group0/compare_attack_h1_recorded.jsonl:1-3;
CAL/lanes/compare-group0/compare_attack_free_entrapment_recorded.jsonl:4 (E 0.71875), 7 (D 0.75415); GIPPA_SPEC.md §8.4 (de novo precision);
CAL/lanes/denovo_precision-group0/placement_m3_21_diagnostics.jsonl:1 › bias_vs_realised, bias_vs_truth_table; tests/calibration/test_denovo_precision_sim.py:113-130 -->

**Table R1.** Full-tier calibration by suite (commit `4a6fa998`, jobs 1966680–1966694). "Records" are harness
records (criterion, level, observed, bound); the status of the last block is counted.

| Suite | R | Tests (pass / skip / xfail) | Harness records | Status for the paper |
|---|---|---|---|---|
| placement | 300 | 87 / 1 / 0 | 366 gated: 365 pass, 1 fail (strict xfail); 6 read-through-balance records vacuous | pass at list level; `flat_c200` at q 0.01 fails (strict xfail); gate M3 not met |
| gene_fdr | 300 | 42 / 0 / 0 | 132 pass | pass; levels 0.05/0.10 weak or vacuous |
| detectability | 100 | 11 / 0 / 0 | 28: 27 pass, 1 recorded | pass (slope); detection-limit CITL recorded |
| regions | 1,000 | 50 / 0 / 0 | 146 pass | pass |
| autocorrelation | 1,000 | 43 / 0 / 1 | 174: 163 pass, 4 unconfirmed (2 exceed the pooled 2R bound: one exact-TU attack; UniOP-like at α 0.01), 7 fail (3 strict xfail, 4 recorded-only) | within tolerance for exact TUs only; anti-conservative with UniOP-like or UniOP K-12 TU errors (recorded) and under the distance-model fallback (strict xfail) |
| islands | 500 | 29 / 0 / 0 | 69 pass; 2 recorded breaches outside the regime | pass within the stated dependence regime |
| discordance | 1,000 | 28 / 0 / 2 | 39: 36 pass, 3 fail (strict xfail) | narrowed claim (smooth, monotone σ) |
| missed_operons | 300 | 27 / 0 / 0 | 48: 22 pass, 26 vacuous | vacuous for the shipped decoys (0 calls); convergent-only arm vacuous at q 0.01 |
| ora | 1,000 | 10 / 0 / 0 | 72: 61 pass, 10 fail (4 UniOP-like arms), 1 unconfirmed (pooled above the 2R bound) | flagged `ora_uncalibrated` |
| gsea | 300 | 55 / 0 / 2 | 111 pass (includes the unflagged-row detection attack gated at 1.5α; one arm 0.0833 at α 0.05, above the α bound 0.0752) | within tolerance under TU exchangeability; restricted claim; no α-level claim for unflagged rows under detection confounding |
| determinants | 100 | 21 / 0 / 0 | 21 pass | flagged `null_uncalibrated` (fails attacks, ADR-0034) |
| novel_orfs | 300 | 21 / 86 / 0 | 125: 116 pass, 9 fail (strict xfail) | class FDR within tolerance in the gate; anti-conservative in the de novo regime when a false ORF carries several wrong peptides; F1 pass; 9 strict xfails |
| denovo_precision | 300 | 2 / 0 / 0 | none (asserted in the test) | pass by assertion for realised precision; bias against generative precision −0.066 (not gated) |
| variants | 100 | 7 / 3 / 1 | 3: 3 fail (strict xfail) | fails F1 ≥ 0.90 (pre-ADR-0049 scope); gate M7 not met |
| compare | 1,000 (ρ CI, H4); 10,000 (H1); 20,000 (free entrapment) | 26 / 0 / 0 | 41: 40 attack records ok, 1 pass | pass; gene-level CI under-covers (recorded) |
| **total** | | **459 / 90 / 6** (555 tests) | | |

<!-- src: CAL/task-0.json … task-14.json; CAL/task-N.*.junit.xml; CAL/<suite>.merged.jsonl and CAL/lanes/<suite>-group*/*.jsonl (status counts; audit R1-A recount); local/logs/nightly.tsv:124 -->

**Table R2.** Observed error against nominal level, representative streams (α or q = 0.05 unless
stated; commit `4a6fa998`). For level and FWER rows the bound uses n = R, or n = R·m when m is shown.
For FDR rows the bound uses the empirical SD of the FDP.

| Suite | Stream | Criterion | Observed | Bound | R (m) |
|---|---|---|---|---|---|
| placement | sim_db_stressed.pooled | FDR | 0.000158 | 0.0500 | 300 |
| placement | readthrough_stressed.shuffle | FDR | 0.0355 | 0.0510 | 300 |
| gene_fdr | picked, π 0.05, q = 0.01 | FDR | 0.00941 | 0.01031 | 300 |
| gene_fdr | MAYU (PEP), π 0.05, q = 0.01 | FDR | 0.00962 | 0.01022 | 300 |
| gene_fdr | picked, π 0.05 (near saturation) | FDR | 0.0486 | 0.0506 | 300 |
| gene_fdr | MAYU (PEP), π 0.05 (near saturation) | FDR | 0.0493 | 0.0505 | 300 |
| regions | HL shift, ρ 0, chromosome pooled | level | 0.0505 | 0.0569 | 1,000 (4) |
| regions | prophage HL shift | FDR | 0.063 | 0.0654 | 1,000 |
| autocorrelation | Moran's I TU-block lag 1 (exact TUs) | level | 0.055 | 0.0638 | 1,000 |
| autocorrelation | same, UniOP-like TU errors (recorded) | level | 0.071 | 0.0638 | 1,000 |
| autocorrelation | same, distance-model TUs (strict xfail) | level | 0.354 | 0.0638 | 1,000 |
| islands | ρ = 0 / 0.5 / 0.9 | FWER | 0.022 / 0.020 / 0.024 | 0.0695 | 500 |
| islands | ρ = 0.9, short lag, 1,000 / 1,800 genes (recorded) | FWER | 0.096 / 0.130 | 0.0695 | 500 |
| discordance | operon null, 1,800 genes | level / FDR | 0.0499 / 0.034 | 0.0510 / 0.0615 | 1,000 (200) |
| missed_operons | convergent-only decoys | FDR (counts) | 0.00065 | 0.0504 | 300 |
| missed_operons | shipped decoys, power arm | FDR (counts) | 0.0 (vacuous) | none | 300 |
| ora | fitted null, over / under | level | 0.0235 / 0.0303 | 0.0556 | 1,000 (6) |
| ora | distance-model large-TU under (fails; flag) | level | 0.0841 | 0.0544 | 1,000 (10) |
| gsea | pooled null, 100 sets per replicate from six collections (R engine) | level | 0.0415 | 0.0525 | 300 (100) |
| gsea | unflagged rows, detection-confounded `null_b2` (gated at 1.5α; Python port engine) | FWER | 0.0833 | 0.1054 (1.5α); α bound 0.0752 | 300 |
| determinants | null, LightGBM / elastic net | level | 0.02 / 0.02 | 0.0936 | 100 |
| novel_orfs | class FDR, picked / de novo | FDR | 0.0453 / 0.0498 | 0.0511 / 0.0507 | 300 |
| novel_orfs | pooled F1 (≥ 0.90) | F1 | 0.938 | 0.90 | 300 |
| variants | pooled F1 of `high` calls (≥ 0.90; strict xfail) | F1 | 0.874 | 0.90 | 100 |
| compare | ρ CI non-coverage (TU block) | level | 0.048 | 0.0638 | 1,000 |
| compare | H4 false confirmation | FWER | 0.047 | 0.0638 | 1,000 |

<!-- src: CAL/lanes/placement-group4/placement.jsonl:20; CAL/lanes/placement-group1/placement.jsonl:2; CAL/lanes/gene_fdr-group1/gene_fdr.jsonl:44,47;
CAL/lanes/regions-group2/regions.jsonl:21; CAL/lanes/regions-group0/regions.jsonl:8; CAL/lanes/autocorrelation-group0/autocorrelation.jsonl:77,44;
CAL/lanes/autocorrelation-group2/autocorrelation.jsonl:25; CAL/lanes/islands-group0/islands.jsonl:32,35,38; CAL/lanes/islands-group0/islands_recorded.jsonl:1-2;
CAL/lanes/discordance-group1/discordance.jsonl:17,21; CAL/lanes/missed_operons-group0/missed_operons.jsonl:38,35; CAL/lanes/ora-group1/ora.jsonl:67,70,42;
CAL/lanes/gsea-group0/gsea.jsonl:45; gsea.attack.flags.detection.null_b2 (CAL/lanes/gsea-group*); CAL/lanes/gene_fdr-group1/gene_fdr.jsonl (q 0.01 rows); CAL/determinants.merged.jsonl:2,5; CAL/novel_orfs.merged.jsonl:38,41,92; CAL/variants.merged.jsonl:3; CAL/lanes/compare-group0/compare.jsonl:41,38 -->

### 1.3 Known failures (expected failures)

A known failure is a test marked xfail: it is expected to fail. At the full tier the marks reported
here are strict, so an unexpected pass turns the run red. Some marks are strict only at the full
tier and non-strict at the fast tier (for example, the distance-model operon fallback and the fast
D-1b-step mode), where R does not resolve the failure; the fast tier recorded 17 xfailed and 3
xpassed tests. The junit roll-up that readiness row R3 reads counted 6 xfailed tests (Table R3a). The merge step of sharded
suites records further strict-xfail verdicts in the merged records, and the roll-up does not count
these (Table R3b).

Three of the merged verdicts appear in the junit files as skipped attack tests, whose verdict is the
merge step:
- `flat_c200`;
- the GC-poor intergenic attack;
- the multi-peptide false-ORF attack.

Readiness R3 therefore reads "6 xfailed", although 13 further merged records (1 placement, 9
novel-ORF, 3 variant) are known failures. Some markers still call ADR-0049 and ADR-0050
"proposed"; both were accepted on 23 October 2026.

**Table R3a.** Strict xfails counted by the junit roll-up (6).

| Suite | Test | Reason recorded in the marker |
|---|---|---|
| autocorrelation | distance-model Moran's I lag 1 (arm B) | ADR-0061: uncalibrated under the distance-model fallback |
| discordance | distance-model p uniformity (arm A) | ADR-0061 |
| discordance | D-1b-step | ADR-0055 Amendment 1: step σ(abundance) outside the smooth-monotone claim |
| gsea | detectability-confounded sets | ADR-0025/ADR-0026: no all-rows level claim |
| gsea | AR(1)-clustered sets (`ar1_run15`) | ADR-0025/ADR-0026 |
| variants | A6a isobaric rows flagged | ADR-0049: isobaric blocks aligned as substitutions |

<!-- src: CAL/task-4.*.junit.xml, task-6.*.junit.xml, task-9.*.junit.xml, task-13.*.junit.xml (xfail entries); local/logs/nightly.tsv:124; tests/calibration/attacks/test_attack_operons_sigma.py:45-52 (fast strict=False, full strict=True); tests/calibration/attacks/test_attack_operon_distance_fallback.py:113-117 (strict=full); tests/calibration/attacks/test_attack_m7_gate.py:295-296 and test_attack_gsea_structure.py:52-53 (strict=True); local/logs/nightly.tsv:137 (fast tier at e43683c7: 17 xfailed, 3 xpassed);
docs/adr/ADR-0049-m7-variant-precision.md:2, docs/adr/ADR-0050-denovo-orf-q-exact-errors.md:2 (Status lines) -->

**Table R3b.** Strict-xfail verdicts in merged records that the roll-up does not count.

| Stream | Criterion, level | Observed | Bound | Marker or ADR |
|---|---|---|---|---|
| placement `flat_c200` contaminant attack | FDR 0.01 | 0.01496 | 0.01019 | contaminant weight pending (M3-25 decision 10) |
| novel ORFs, multi-peptide false ORFs | FDR 0.01 / 0.05 / 0.10 | 0.0183 / 0.0788 / 0.1139 | 0.0109 / 0.0511 / 0.1012 | M7-11 decision 3 (ORF q = min peptide q) |
| novel ORFs, GC-poor intergenic | FDR (counts) 0.01 | 0.01258 | 0.01122 | ADR-0050 |
| novel ORFs, G4-k3 intergenic | FDR (counts) 0.01 | 0.01529 | 0.01136 | ADR-0050 |
| novel ORFs, G4-k3 CDS out-of-frame | FDR (counts) 0.01 / 0.05 | 0.0749 / 0.0749 | 0.0102 / 0.0502 | ADR-0050 |
| novel ORFs, G4-k3 antisense | FDR (counts) 0.01 / 0.05 | 0.0719 / 0.0719 | 0.0102 / 0.0502 | ADR-0050 |
| variants A6b (random isobaric swaps) | FWER 0.05 | 0.97 | 0.0936 | ADR-0049 |
| variants A6c (systematic swaps) | FWER 0.05 | 1.0 | 0.0936 | ADR-0049 |
| variants G1 | pooled F1 ≥ 0.90 | 0.874 (SE 0.0018) | 0.90 | ADR-0049 |

<!-- src: CAL/placement.merged.jsonl:43; CAL/novel_orfs.merged.jsonl:22-24,43,117,120-121,123-124; CAL/variants.merged.jsonl:1-3 (outcome "xfail") -->

## 2. Attack nulls

Several attack nulls fail, and they are known failures. The specification's isobaric-swap attack
is tested in the variants gate, not in the attack files, and it fails: A6b gave FWER 0.97 and A6c
gave 1.0, against 0.0936 (Table R3b). Also failing as strict xfails (Section 1.3) are the
distance-model fallback, the D-1b-step σ attack, two GSEA restricted-claim structures, A6a isobaric
flagging, contaminant weighting (`flat_c200`), multi-peptide or systematic de novo false ORFs, and
GC-poor intergenic exact products (ADR-0050). The flagged methods (ORA,
determinants) pass their attacks only through their flags. The unflagged GSEA rows pass their
detection attack only under a 1.5α screen (Python port engine; Section 1.2).

The attack nulls are simulations built to make a method reject falsely, including:
- autocorrelated abundance without islands;
- detectability confounded with region membership;
- gene sets packed into one operon;
- high-copy plasmids;
- de novo errors that reverse a stretch;
- isobaric swaps everywhere;
- further method-specific attacks, in 55 attack test files.

The high-copy plasmid attack is a deterministic invariance check, not a rejection rate: results at
20 plasmid copies must equal those at one copy.

In the full tier at `4a6fa998`, readiness counts all 14 attack suites as passed: of 281 attack
tests, 273 passed, 3 were skipped and 5 xfailed (0 failed, 0 errors, 0 xpassed).
- The 5 xfails are the five attack tests in Table R3a: two distance-model fallback tests, D-1b-step
  and two GSEA restricted-claim tests.
- The 3 skipped attack tests are those whose verdict is a merged strict xfail (Table R3b).

"14 of 14 suites passed" thus means only that no attack test failed outside the known failures
listed above and in Section 1.3.

<!-- src: local/logs/nightly.tsv:125; validation/release/1.1.0/READINESS.md:12 (R4); GIPPA_SPEC.md:1980-1991; tests/calibration/attacks/test_region_tests_attacks.py:24-28 (A4 invariance); tests/calibration/attacks/ (55 files);
CAL/task-*.junit.xml (attack test outcomes); CAL/variants.merged.jsonl:1-2 (tests/calibration/test_variants_gate.py); docs/adr/ADR-0034-determinants-null-uncalibrated.md (Decision) -->

## 3. True-example assertions and operon calibration on real genomes

The specification defines eight true-example assertions (TE1–TE8) on real genomes. Readiness row R5
lists seven of them. It is recorded as **waived**, because TE4 failed (Table R4). TE6 was measured in
the weekly K-12 run but is not a readiness row.

**Table R4.** True-example assertions.

| TE | Assertion (short) | Outcome | Value | Code state |
|---|---|---|---|---|
| TE1 | ≥ 99.5 % of A2 peptides map exactly, in frame, to annotated CDS | pass | 28,129 of 28,129 (1.0; Wilson 95 % lower bound 0.99986) | 1.0.0-era pilot |
| TE2 | no `high` ALT call off start or stop codon sites | pass, **vacuous by design** | 0 ALT sites, 0 `high`, 0 violations: with a same-database search A2 cannot yield ALT calls, so it checked nothing | 1.0.0-era pilot; release-equivalent weekly in-silico proxy also vacuous (0 ALT sites) |
| TE3 | ribosomal proteins ≥ 90 % detected, median abundance in the top quartile | pass | 53 of 53 detected (1.0 [0.932, 1.0]); median percentile rank 0.987 (iBAQ) | 1.0.0-era pilot |
| TE4 | log2(RplL/RplJ) ∈ [1, 3]; log2(AtpA/AtpD) ∈ [−1, 1] | **fail (waived)** | rplL/rplJ −0.480 (Top3 −0.531); atpA/atpD −0.118 (pass) | 1.0.0-era pilot |
| TE5 | trypsin and keratin peptides land in the contaminant bucket | pass | A1 59 of 59; A2 61 of 61 | 1.0.0-era pilot |
| TE6 | K-12: Bakta oriC next to mnmG (mnmG or mioC nearest, mnmG within 2 kb) | pass by the any-call rule; 1 of 2 oriC calls meets the criterion; not a readiness row | "2 oriC calls, at the true origin: 1" | release code (weekly 1972755_0) |
| TE7 | same seed gives byte-identical tables | pass (GIPPA tables only) | 2,786 gate files compared, 0 differ (13,439 non-gate files excluded; 81 figures identical); DIA-NN layer not run | 1.0.0-era pilot |
| TE8 | UniProt K-12 peptides ≥ 99.5 % exact on K-12; Nissle placement FDR ≤ 1 % | pass | K-12 exact 0.99900 (151,272 of 151,424); Nissle placement FDR 0.00533 | release code (weekly 1972755_0) |

<!-- src: validation/release/1.1.0/READINESS.md:13 (R5); validation/pilot/metrics/metrics.json › te1, te_real.te2, te_real.te3, te_real.te4.pairs, te_real.te5, te7.gippa (n_gate, n_excluded, n_figures, figures_differing);
local/ci/validation/2026-11-10_te2_1969418.json › counts; local/logs/gippa-weekly.1972755_0.out:319-320; scripts/slurm/weekly.sbatch:19-20, 170-192 (TE6 per-call criterion; step passes if any call meets it); GIPPA_SPEC.md:2373-2378; local/ci/validation/2026-11-11_te8_1972755.json › counts, shares, passed -->

**TE1.** Of the 28,129 non-contaminant A2 peptides, 18 exact in-frame peptides were masked, 42
started at non-ATG N-termini, and none missed. The pilot report gives two caveats. Every A2
identification is a substring of an annotated protein, so TE1 checks that database, genome and
mapping agree, and "says nothing about identification error". Misses
cluster by gene, so the Wilson interval's coverage is optimistic, and "the interval is not a test".

<!-- src: validation/pilot/metrics/metrics.json › te1 (n_peptides 28129, exact_inframe_masked 18, non_atg_nterm.peptides 42, n_misses 0), m2.A2.n_denominator 28129, m1.A2.union.peptides.non_contaminant 28129;
validation/pilot/REPORT.md:174-180 -->

**TE2.** TE2 passed in the pilot with zero ALT sites in A2 (`vacuous: true`). Every A2
identification is a substring of an annotated protein of the searched database, so A2 cannot yield
an ALT call or a novel ORF, so TE2 on A2 is uninformative by design. The weekly in-silico proxy at
`20283834` passed again, also vacuously (0 ALT sites). Its note gives two reasons why no
stop-codon ALT site reached high confidence. The selenocysteine peptides of one formate
dehydrogenase subunit were rejected by the stratum-1 placement FDR. The other Sec protein has no
tryptic peptide over its Sec codon. An earlier in-silico proxy
(2026-10-23, before the stop-codon read-through changes) had 233 ALT sites and 80 `high` calls, all
at start codons, with 0 violations. That proxy predates the current code.

<!-- src: validation/pilot/metrics/metrics.json › te_real.te2.vacuous; local/ci/validation/2026-11-10_te2_1969418.json › notes; validation/proxies/2026-10-23_te2.json › counts -->

**TE4.** The waiver review recorded an iBAQ SD of 0.99 among 48 ribosomal proteins expected at 1:1.
It also found that 53 of 2,766 A2 genes have a peptide that swings by more than 3 log2 between runs.
The cause was not proven. Readiness applies the 1.0.0 waiver record ADR-0070. ADR-0076 extends it to
1.1.0 for this pair only.

<!-- src: docs/adr/ADR-0070-te4-rplL-rplJ-waiver-1-0-0.md:39,51; validation/release/1.1.0/readiness.json › rows[R5].waiver; docs/adr/ADR-0076-waivers-1-1.md:28-33 -->

**TE6.** The weekly check evaluated the specification's criterion for each oriC call: mnmG or mioC
is the nearest gene, and mnmG lies within 2 kb. One of the two calls met it, and the step passes
when any call does. The second call does not meet the criterion.

**TE7 and TE8.** TE7 compared two GIPPA runs on the same DIA-NN output at the freeze commit, so
reproducibility from raw spectra to tables was not shown. Readiness row R6 adds the release-code
integration run (967 passed, 84 skipped, 1 xfailed, 1 xpassed) and nf-test (144 tests). The
integration suite contains a two-run byte-identity test on a simulated 100-kb genome
(`test_m10_06_sim_te7.py`). It ran and passed in the release-code run 1973700_0: the run has no
`-m` filter, and none of the 84 skips is that test. The first integration run on the release commit
(1972753_0) had one failure, a known test flake. The rerun cited by R6 (1973700_0) passed. TE8 also
reported that, on the Nissle genome, 0.764 of K-12 peptides mapped exactly, 0.0895 mapped tolerantly
and 0.147 were unmapped.

<!-- src: validation/pilot/REPORT.md (Notes: TE7); tests/integration/test_m10_06_sim_te7.py:1-14; scripts/slurm/nightly.sbatch (task_integration: pytest tests/integration tests/contract -q); local/ci/runs/1973700_0/integration.txt:122 and -ra skip summary (84 skips: test_envs_built 15, test_m7_gate 64, test_m4_gate 2, test_m5_gate 2, test_bakta_regions 1); scripts/slurm/nightly.sbatch:214; local/logs/nightly.tsv:135-136; docs/backlog.md:1474-1475 (flake, test_m11_19 example 21); local/ci/validation/2026-11-11_te8_1972755.json › shares -->

**Operon predictor on K-12.** The UniOP AUROC against 2,739 labelled RegulonDB 14.5.0 pairs (1,730
positive, 1,009 negative) was 0.911, with a directon-bootstrap 95 % percentile CI of [0.898, 0.923].
It applies to K-12 only, and may be partly in-sample, because UniOP's training data were not
verified.

At the shipped pair threshold of 0.5, UniOP had precision 0.896 and recall 0.853 (F1 0.8739).
Precision depends on prevalence; the K-12 labelled π is 0.632. The intergenic-distance model is the
fallback when UniOP is missing:
- its AUROC was 0.907 [0.892, 0.920], and UniOP exceeded it by 0.0041 [0.0009, 0.0077];
- at the same shipped threshold (p ≥ 0.5, d ≤ 102 bp), its precision was 0.828 [0.805, 0.847] and
  its recall 0.916 [0.901, 0.931] (TP 1,584, FP 328), below the 0.90 precision target of
  specification §15.4.

The AUROCs come from Slurm job 1879748 (run 24 October 2026, 01:05 UTC). The distance-model
operating point comes from the ADR-0022 addendum (Slurm 1892896). Both runs predate the 1.0.0
freeze.

<!-- src: docs/adr/ADR-0022-uniop-terminatornet-cutoffs.md:3, 78-104 (AUROC line 84), 120-126, 248-286 (addendum table line 274; below target line 286);
src/gippa/data/calibration.yaml:14-22 (precision_at_pair_threshold 0.8285 [0.8048, 0.8475], recall 0.9156 [0.9008, 0.9305]); src/gippa/config/defaults.yaml:32 (pair_threshold 0.5, fallback distance) -->

The weekly distance-model run (`20283834`) measured a different quantity: TU pairs from a full GIPPA
run with UniOP disabled, against RegulonDB TU pairs, gave precision 0.708 and recall 0.557 (TP 963,
FP 397, FN 767). This pipeline-level agreement is not comparable with ADR-0022's pair-classification
precision. The
matching UniOP comparison in the same weekly run was "not available".

<!-- src: local/logs/gippa-weekly.1969418_2.out:116; local/logs/weekly.tsv:9,12 -->

## 4. The Nissle 1917 pilot (PXD078459)

All pilot results were produced at the 1.0.0 report-freeze commit `fa49b8f1` and were not re-run
on 1.1.0 code (see Provenance). Readiness 1.1.0 reads the same committed `metrics.json`.

The input was the wild-type runs of PXD078459, searched with DIA-NN 2.3.1 Academia, not the 2.7.0
named in the specification (ADR-0051). The three arms were:
- **A1:** the UniProt K-12 proteome UP000000625 plus contaminants;
- **A2:** the Nissle strain database plus contaminants;
- **A3:** Cascadia 0.0.7 de novo sequencing.

Several inputs were not produced at the freeze:
- the freeze runs used the pass-3 wild-type genome and its recorded Bakta files;
- A1's peptide table is the pass-1 conversion, reused because the K-12 search does not depend on the
  Nissle genome;
- the A3 Cascadia table is the pass-2 conversion, which was not re-run;
- the free-entrapment rerun that counts ran from main commit `db3a95f1`. Its `src/` and `workflow/`
  equal those of `fa49b8f1`, but it includes a fix to the entrapment script
  (`validation/pilot/entrapment.py`: initiator M at annotated non-ATG starts).

<!-- src: validation/pilot/REPORT.md:1-30, 306-354; validation/pilot/metrics/metrics.json › gippa_commit, inputs.peptides_a1; validation/pilot/README.md:603 (A1 pass-1 table); validation/pilot/REPORT.md:346 (db3a95f1); git diff --stat fa49b8f1 db3a95f1 -- validation/pilot/entrapment.py (+50/−15);
docs/adr/ADR-0076-waivers-1-1.md (Decision, Scope); git diff --stat fa49b8f1 e43683c7 -- src workflow -->

### 4.1 End-to-end runs, identifications and mapping

GIPPA completed on all three arms at the freeze commit (A2 job 1914948, A1 1915264, A3 1916086, all
exit 0). Table R5 gives identification and mapping aggregates. Precursors were filtered at DIA-NN's
global q ≤ 0.01. That 1 % is the search engine's own estimate, which GIPPA inherits and does not
verify; on A1 the free-entrapment point estimate is 0.024 (Section 4.5). Peptide and protein-group
counts are not FDR-controlled.

**Table R5.** Identifications and mapping, A1 and A2 (union over runs).

| Quantity | A1 (K-12 reference) | A2 (strain database) |
|---|---|---|
| precursors (global q ≤ 0.01) | 28,639 | 33,764 |
| peptides: all / non-contaminant / contaminant-only | 23,792 / 23,704 / 88 | 28,234 / 28,129 / 105 |
| protein groups / with ≥ 2 peptides | 2,798 / 2,433 | 3,228 / 2,790 |
| exact mapping share [Wilson 95 %] | 0.9936 [0.9925, 0.9945] | 1.0 [0.99986, 1.0] |
| tolerant-only share | 0.00515 | 0 |
| unmapped share | 0.00127 | 0 |
| genome-inconsistent peptides (from `detection_summary.tsv`) | 1,450 of 25,002 (0.058) | not applicable |

The genome-inconsistent count uses a different denominator from the share rows. It is not the sum of
the tolerant-only and unmapped peptides.

<!-- src: validation/pilot/metrics/metrics.json › m1.A1/A2.union; m2.A1/A2.shares; m2.A1.genome_inconsistent (count 1450, n 25002, value 0.0580, consistent false);
validation/pilot/REPORT.md:11 (job IDs), 58 -->

The decoy-estimated mapping FDR of A2 (readiness R2) was not measured: all 28,129 A2 peptides come
from the searched database (`same_db`) and are exempt from both decoy lists, so the row is waived. A
context-only value on A1 was 7.27e-05 over 23,675 accepted peptides. The waiver records that this
value cannot fail by construction and does not validate A2. Identification error in A2 is inherited
from DIA-NN (input q ≤ 0.01) and is not estimated by GIPPA.

<!-- src: docs/adr/ADR-0071-r2-mapping-fdr-waiver-1-0-0.md:15-17, 30-41; validation/release/1.1.0/READINESS.md:10 (R2) -->

### 4.2 Detection and detection tiers

The gene set has 4,685 strain genes (CDS and sORF). A1 detected 2,645 of them
(0.5646 [0.5503, 0.5787]) and A2 detected 3,051 (0.6512 [0.6375, 0.6647]). Of the 916
strain-specific genes (no tryptic peptide shared with K-12), A1 cannot detect any by construction,
because it searched only the K-12 proteome; it detected 0. A2 detected 316 (0.3450 [0.3149,
0.3764]). Table R6 gives shares by region. The freeze run manifests record the
detection tier `decoy_picked` for A1 and A2, and `two_peptide` for A3, because de novo rows carry no
decoys; that tier has no FDR estimate (Section 1.2: simulated gene FDP 0.054 at q = 0.01).

**Table R6.** Genes detected by subset (share [Wilson 95 %]). The Wilson intervals treat genes as
independent. Genes of one TU, and peptides of one protein, co-vary, so the intervals (here and in
Table R5) are optimistic. They describe this single sample (three wild-type runs of one condition,
union over runs), not a population quantity.

| Subset (genes) | A1 | A2 |
|---|---|---|
| all (4,685) | 0.5646 [0.5503, 0.5787] | 0.6512 [0.6375, 0.6647] |
| strain-specific (916) | 0 (structural; not an estimate) | 0.3450 [0.3149, 0.3764] |
| plasmid (13) | 0.0769 [0.0137, 0.3331] | 0.7692 [0.4974, 0.9182] |
| prophage (136) | 0.0221 [0.0075, 0.0628] | 0.2721 [0.2043, 0.3523] |
| BGC (251) | 0.3028 [0.2493, 0.3623] | 0.5458 [0.4840, 0.6063] |
| K-12-shared (3,759) | 0.7036 [0.6888, 0.7180] | 0.7276 [0.7131, 0.7416] |

<!-- src: local/scratch/a1-recheck/derived_fix.json › m3.A1/A2 (1.1.1 detection recompute, Slurm 2002006; the 1.0.0-era values in validation/pilot/metrics/metrics.json › m3 are reproduced byte-identically by the control, Slurm 2001992, derived_ctl.json); local/scratch/a1-recheck/out_fix.jsonl (freeze/fa49b8f/out_a1, out_a2: genes_q01); validation/pilot/REPORT.md:81-84 (G = cds, sorf);
local/runs/pilot/freeze/fa49b8f/out_{a1,a2,a3}/run_manifest.json › detection_tier, git_commit -->

### 4.3 Comparison of A1 and A2

On the 2,613 genes with measured abundance in both arms (1,674 TUs), the Spearman correlation of
abundance was 0.978, with a TU-block bootstrap 95 % CI of [0.973, 0.983] (1,000 resamples). The
Bland–Altman bias was −0.194 log2, with limits of agreement [−1.362, 0.975]. The Jaccard index of
enriched functional categories was 0.588 (405 in A1, 292 in A2, 258 shared).

Both arms quantify the same spectra. The correlation therefore measures how far the choice of
search database perturbs abundance on genes quantified in both arms; it does not measure
reproducibility or biological agreement. The enriched-category sets combine ORA rows, which ship
flagged `ora_uncalibrated`, with GSEA rows at q ≤ α, so this Jaccard is descriptive.

No islands and no operon flags were called in either arm, so their Jaccard indices are undefined. The
Jaccard indices for ALT sites (118 in A1, 0 in A2) and novel ORFs (1 in A1, 0 in A2) were 0. Jaccard
indices carry no uncertainty estimate.

<!-- src: local/scratch/a1-recheck/derived_fix.json › m4, m5.enriched_categories (1.1.1 detection recompute, Slurm 2002006); validation/pilot/metrics/metrics.json › m5.islands, m5.operon_flags, m5.other_families (unchanged in local/scratch/a1-recheck/derived_fix.json); m5.islands; m5.operon_flags; m5.other_families; src/gippa/compare/stability.py:40-52 (enriched = ORA ∪ GSEA at q ≤ alpha); validation/pilot/REPORT.md:174-176 (same spectra, same database) -->

### 4.4 Pre-registered hypotheses

Specification §15.7 fixed four hypotheses before the pilot, to be "reported where the pilot data
allow". Table R7 quotes them and gives their recorded outcomes.

**Table R7.** Hypotheses H1–H4 (wording of specification §15.7).

| Hypothesis | Recorded outcome | What was evaluated |
|---|---|---|
| H1: "A2 detects more strain CDS than A1, mostly plasmid, prophage and BGC genes." | **not supported** (failed clause: "mostly") | 434 genes detected only in A2, 28 only in A1; A2-only share 0.9394 [0.9138, 0.9577]; TU-level test 290 vs 24 units, p_tu 1.80e-59, unit effect 0.9236 [0.8888, 0.9481]; "mostly" share 106 of 434 = 0.2442 [0.2062, 0.2868] against a base rate of 0.0854; 316 of the 434 A2-only genes (0.728 [0.684, 0.768]) are strain-specific and undetectable in A1 by construction, so the "more" clause is expected by design |
| H2: "A1's genome-inconsistent FDP is at least its entrapment FDP (full validation, FDRBench), since it errs high." | **deferred** to the full validation | not tested; shown side by side only: A1 genome-inconsistent share 0.058 and A1 free-entrapment estimate 0.0273 (peptide level) |
| H3: "Context conclusions (islands, region tests, operon concordance) agree between A1 and A2 on shared genes." | recorded verdict "supported" under the implemented rule (ρ and enriched-category Jaccard); **none of the three named conclusions was evaluated** | islands and operon flags not evaluable (none in either arm); region tests not compared; Jaccards not restricted to shared genes; evaluated instead: ρ 0.978 [0.973, 0.983] and enriched-category Jaccard 0.588 |
| H4: "The genome-anchored de novo precision agrees with precision measured against A2 identifications." | **consistent with lower bound** (the data can refute but never confirm) | 10 reportable bins, 1,364 queries; all 10 bins consistent; n-weighted mean genome precision 0.562; pooled A2-confirmed share 0.543 [0.517, 0.570] |

<!-- src: GIPPA_SPEC.md:2696-2703 (§15.7 wording); local/scratch/a1-recheck/derived_fix.json › h1 (b 434, c 28, b_tu 290, c_tu 24, p_tu, unit_effect, mostly, strain_specific_share 316/434; 1.1.1 detection recompute, Slurm 2002006), h3 (enriched_categories jaccard 0.587699, verdict supported); validation/pilot/REPORT.md:233 (H2 shown side by side), 312-313;
metrics.json › h3.criteria (rho holds, enriched_categories holds, islands/operon_flags not_evaluable), h3.region_tests ("not compared in pass 1"), h3.verdict; metrics.json › h4 -->

H3 as pre-registered was therefore not tested. The recorded "supported" rests only on criteria
added by the implementation. One of them, the enriched-category Jaccard, depends on uncalibrated ORA
q-values; the other, ρ, compares two searches of the same spectra (Section 4.3).

Specification §15.7 fixed the wording of H1, not its test. The test unit (the TU), the TU-level test
and the rule that the "mostly" share must exceed 0.5 were specified in `docs/specs/M10-05`. That
file was committed at 2026-10-24 04:22 PDT. At that time the DIA-NN reports had been converted
(2026-10-23, 23:02), but no GIPPA pilot run had yet produced gene-level detections.

For H1, the gene-level p (5.20e-95) decides nothing: the compare attack
in Section 1.2 recorded this p as anti-conservative (0.1737 at α = 0.05), and the verdict rests on
the TU-level test and the "mostly" clause. The table uses the specification's wording of H2,
not the pilot report's.

For H4, A2 enters both sides (its matches set the c-value; its identifications set the
confirmations), so the check is not independent. The score bins span de novo scores in (0.8, 1]
only, and DIA data give no spectrum pairing.

<!-- src: local/scratch/a1-recheck/derived_fix.json › h1.p_gene (5.202e-95; 1.1.0 value 4.65e-86 in derived_ctl.json); CAL/lanes/compare-group0/compare_attack_h1_recorded.jsonl:2; docs/specs/M10-05-pilot-metrics.md:266-298 (commit e81c1559, 2026-10-24T04:22-07:00); sacct 1880208 start 2026-10-23T23:02 (DIA-NN convert, validation/pilot/README.md:217-242); first GIPPA pilot runs 1892236 ff. (README.md:326-329); validation/pilot/REPORT.md:233 (H2 wording), H4 caveats -->

### 4.5 Free entrapment (A1)

K-12 peptides absent from the Nissle genome served as free entrapment for A1.
- **Precursor level** (N_T 28,347, N_E 157): the combined FDP estimate was 0.0236, the realised lower
  bound `fdp_lower` 0.00551, and the one-sided 95 % upper bound `ucb95` 0.0265. The assessment was
  **inconclusive**: the point estimate exceeds the 1 % target, and the interval does not exclude
  1 %, so the target was neither confirmed nor shown to fail. Under clustered false matches
  (simulation arms D and E, Section 1.2), ucb95 covered only 0.719–0.754.
- **Peptide level** (N_T 23,517, N_E 152): 0.0273, 0.00642 and 0.0308; not assessed, because the
  search does not control peptide-level FDR.
- The ratio r of K-12-only to shared in-silico entries was 0.307 for peptides and 0.304 for
  precursors.

The estimates rest on stated assumptions:
- `ucb95` assumes equal-chance, independent false matches and is not valid under clustered errors;
- the combined estimate is an estimated upper bound in expectation, not a bound on this sample;
  under clustered false matches (simulation arm D) it averaged 0.009 against a true FDP of 0.016, so
  it is not an upper bound in that regime;
- the lower bound assumes that the Nissle genome is complete, that it is the sample's strain, and that
  no K-12 material is present.

Known biases run both ways, and neither is quantified: near-homolog entrapment pushes the estimate
high, and DIA-NN's use of protein information in A1's q-values pushes it low. Restricted to
entrapment peptides far from any Nissle homolog (N_E 7), the peptide-level combined estimate was
0.00274 (ucb95 0.00493).

<!-- src: validation/pilot/metrics/metrics.json › m8.precursor, m8.peptide, m8.universe.r, m8.universe.r_prec, m8.peptide.sensitivity.far, m8.caveats; validation/pilot/REPORT.md:332; CAL/lanes/compare-group0/compare_attack_free_entrapment_recorded.jsonl:4,7 -->

### 4.6 Novel ORFs and variants

A2 called no novel ORF and no ALT variant, as its same-database search allows none (Section 3). The
Jaccard indices of 0 against A1 for these families (Section 4.3) therefore reflect the design. A1 called one novel ORF (21 more were unestimated), which the pilot
report attributes to an annotation artefact, not a novel protein. A1 had 4 `high` and 114 candidate substitution ALT calls.

For A3:
- 0 novel ORFs were called (5 not called);
- there was 1 `high` and 125 candidate swap calls, 41 candidate substitutions and 28 candidate blocks;
- the class error rates of A3's ALT calls are unestimated.

A2 and A3 shared no ALT site (0 in A2, 195 in A3).

<!-- src: validation/pilot/metrics/metrics.json › m6.A1/A2/A3.novel_orfs.totals, m6.*.variants.alt, m6.overlap_a2_a3.alt_sites; validation/pilot/REPORT.md:99-143 and Notes (Novel ORFs in A1) -->

### 4.7 A3 feasibility

Before arm A3 was run, a feasibility rule required at least 500 distinct 7–30-mers with a de novo
score above 0.8. It also required that at least 0.20 of the top 1,000 match the six-frame translation
of the genome.

The committed verdict file covers one run, WT_1: 863 peptides passed, 549 matched, a share of 0.636
(Wilson 95 % CI 0.604–0.668), verdict pass. The shuffled-decoy share was 0.0058. This CI is
optimistic, because peptides of one protein co-vary.

The committed file was evaluated on the engineered chassis genome; the freeze used an uncommitted
re-evaluation of WT_1 on the wild-type genome (also a pass, equal aggregates). The README lists WT_2
(0.644, 0.611–0.675) and WT_3 (0.623, 0.589–0.654) as passing, without committed verdict files. The check is a precondition for A3, not a GIPPA A3 result.

<!-- src: validation/pilot/metrics/a3_feasibility.json › run, n_pass, n_match, share, share_ci95, decoy_share, rule, matches_by_replicon;
validation/pilot/README.md:443-452, 604-606; validation/pilot/REPORT.md:349 -->

### 4.8 Pilot runtimes

On the shared cluster with 4 allocated CPUs, the pipeline wall time recorded in the run manifest was
2:36:07 for A1 (sum of step wall times 10:05:49) and 2:50:20 for A2 (12:05:36). The corresponding
sacct elapsed times were 2:36:36 and 2:50:55. The largest single-step peak RSS was 7.5 GiB in both
arms. These are pilot measurements on a shared cluster, not the benchmark.

<!-- src: validation/pilot/metrics/metrics.json › m7.A1/A2.wall_time_s, sum_task_wall_s, sacct.Elapsed, sacct.AllocCPUS, max_step_peak_rss_bytes; validation/pilot/REPORT.md §8 -->

## 5. Performance

**Standard-run benchmark.** The benchmark ran at the 1.0.0-era commit `fa49b8f1` on the K-12 MG1655
genome and on a simulated 5-Mb genome (`sim5m`). Both inputs used simulated peptide tables from
gippasim, so "K-12" here means a real genome with simulated proteomics. The run settings were:
- one FarmShare "wheat" node (wheat-01, Intel Xeon Gold 6148);
- 16 threads and 64 GB, with 56 GB per task;
- profile `farmshare_local,apptainer`;
- array 1916711, on 31 October 2026.

Each cell ran twice; repeat 3 was cancelled (ADR-0072). The "median" is the mean of 2 runs, with
unestimated uncertainty. The node was busy, with a median allocation of 60–80 of
80 CPUs. Table R8 gives the results.

**Table R8.** Standard-run benchmark (1.0.0-era, 2 repeats per cell; simulated peptide tables).

| Cell | Modules; Bakta | Repeats (s) | Mean (s) | Mean (h:mm:ss) | Spread | Peak RSS (sacct max, GiB) | Estimate (s) | Ratio to estimate |
|---|---|---|---|---|---|---|---|---|
| K-12 standard | core, mobile, bgc, function; live, full DB | 12,502; 13,963 | 13,232.5 | 3:40:32 | 0.110 | 62.35 | 2,100 | 6.30 |
| K-12 core | core; recorded JSON | 10,896; 10,851 | 10,873.5 | 3:01:13 | 0.004 | 37.87 | 1,500 | 7.25 |
| sim5m standard | core, mobile, bgc, function; live, full DB | 14,578; 12,421 | 13,499.5 | 3:44:59 | 0.160 (noisy) | 62.38 | 2,100 | 6.43 |
| sim5m core | core; recorded JSON | 12,146; 12,945 | 12,545.5 | 3:29:05 | 0.064 | 41.20 | 1,500 | 8.36 |

The peak RSS of 62.35–62.38 GiB in the standard cells is 99.8 % of the 62.5 GiB (64,000 MB)
allocation. Peaks
differed widely between repeats (K-12 standard 62.4 and 49.2 GiB; K-12 core 26.4 and 37.9 GiB).
ADR-0072 states that it is not established why, or whether the 62 GiB peaks reflect a limit, and that
64 GB is "the tested allocation, not … a measured requirement".

<!-- src: benchmarks/baseline.json › standard_run.wheat.* (median_s, runs_s, spread, noisy, peak_rss_sacct_max_bytes, estimate_s, ratio_to_estimate, nodes, cpu, mem_mb 64000, task_memory "56.GB", rounds, slurm_jobs, date, commit);
docs/adr/ADR-0072-benchmark-baseline-1-0-0.md:15-22 (Decision 1: "Both use simulated peptide tables (gippasim)"), 25-43 (Evidence, "Not established") -->

The specification estimated a standard run at about 35 min (range 16–55 min). The measured K-12
standard run took 3:41 (13,232 s). Readiness row R7 therefore failed its range and is **waived**;
the runtime claim is not met. Row R8 (core modules, "about 25 min", reported only) reads 3:01
(10,874 s). Two steps, `promoter_calculator` and `terminatornet_chunk`, took about 92 % of summed
task time in repeat 1 of the K-12 standard cell. The h:mm:ss means are as ADR-0072 writes them; the
pilot report rounds three of them up by one second.

<!-- src: validation/release/1.1.0/READINESS.md:15-16; docs/adr/ADR-0072-benchmark-baseline-1-0-0.md:39-40; validation/pilot/REPORT.md (Notes: Runtime) -->

**Other timings.**
- **Demo.** `gippa demo` from the 1.1.0 PyPI wheel finished in 75 s on a wheat node with Python
  3.12.3. The smoke script measured 77 s, including start-up.
- **Example fetch.** `gippa example fetch` from the v1.1.0 GitHub release exited 0.
- **Micro-benchmarks.** The nightly suite passed 59 tests at `4a6fa998` (580.09 s, wheat-05), after a
  first run at the same commit failed two timing budgets. It passed 59 at `a7562487`. It failed three
  timing budgets on the tagged commit `a2a95c76` (code-identical to `e43683c7`), then passed 59 on
  each night from 13 to 16 November at `37a01337`.
- **End-to-end on a simulated genome.** The nightly run on the 100-kb simulated genome exited 0, with
  both the light Bakta database and recorded full-DB JSON, on every night from 10 to 12 November.
  Those commits differ from the release at most in the version string. The run on 9 November, at
  development commit `e6f4aa3c`, also exited 0.

<!-- src: local/logs/gippa-pypi-smoke.1977437.out:3,9,14; local/logs/gippa-pypi-smoke.1977441.out:30; local/logs/nightly.tsv:120,126,131,139,142,147,151,155 (benchmarks);
local/logs/nightly.tsv:117,121,133,138 (sim100k); git diff --shortstat <commit> e43683c7 -- src workflow -->

## 6. Example bundle routes A–E

The example bundle `nissle_r200k` (bundle version 1.1) holds the pilot's data restricted to one
200-kb region (routes A–D). Route E instead uses a separate 200-kb region of the post-freeze exploratory
six-frame arm A4 (DIA-NN job 1918156), which was not part of the pilot gate.

Each route ran once (Table R9) at the development commit `39f9f3a2` on 9 November 2026, with modules core and
bgc, the light Bakta database, TransTermHP terminators, Promoter Calculator off and 4 threads. For
every route, all measured outputs equalled the bundle's expected values on every key except wall
time. GIPPA itself produced those expected values, so the match is a regression and reproducibility
check of the routes, not validation against external truth. The bundle was built at commit
`52ddd95f`, which is not an ancestor of the release commit.

**Table R9.** Example routes (one run each, commit `39f9f3a2`).

| Route | Input | Detection tier | Genes | Detected | Detected (≥ 2 peptides) | With abundance | Peptides placed | Novel ORFs (candidates / called) | Wall time (s) |
|---|---|---|---|---|---|---|---|---|---|
| A | DIA-NN report, `--search-db same` (2,334 rows) | mayu | 158 | 99 | 91 | 102 | 774 | – | 339.2 |
| B | protein FASTA and abundance table (106 records) | protein_list | 158 | 106 | 0 | 106 | 0 | – | 500.6 |
| C | Cascadia `.ssl` × 3 (76 rows; format demonstration) | two_peptide | 158 | 4 | 4 | – | 21 | – | 270.0 |
| D | generic peptide table (2,088 rows) | mayu | 158 | 99 | 91 | 102 | 774 | – | 273.8 |
| E | DIA-NN report against a six-frame database (arm A4; 3,650 rows) | decoy_picked | 178 | 114 | 100 | 114 | 1,041 | 8 / 0 | 641.8 |

<!-- TODO ADR-0078 route E: was run at `39f9f3a2` (1.1.0 development code) in the `decoy_picked` tier, so its Detected, Detected (≥ 2 peptides) and With abundance values use the 1.1.0 picked gene FDR. On the full arm A4 the fix leaves 2,965 detected genes unchanged (maximum detected q 0.0055 to 0.0067; local/scratch/a1-recheck/out_v110.jsonl, out_fix.jsonl › pilot/sixframe-fa49b8f/out_a4), but the 200-kb route E subset was not recomputed. Re-run route E with 1.1.1 (and update the bundle's expected/E.json if it changes) before quoting it as 1.1.1 output. -->

<!-- src: local/example/routes/1966530/gippa-example/runs/measured/{A..E}.json; runs/wall_times.tsv; expected/{A..E}.json; bundle.json › bundle_version, build_commit, counts;
runs/A/run_manifest.json › git_commit; gippa-example/PROVENANCE.md:36-45 (route E from arm A4, Slurm 1918156); docs/adr/ADR-0074-example-dataset.md (decision 3; Amendment 2026-11-06) -->

Wall times include conversion and preflight and exceeded the bundle's informational expected
wall times (247.3, 409.0, 210.0, 188.9 and 439.5 s for A–E) by 22–46 %. Route E's zero called novel
ORFs is not a negative result: in all four decoy-tagged classes the smallest reachable q was at least
0.2, above α, so no call was possible at these class sizes.

<!-- src: local/example/routes/1966530/gippa-example/expected/{A..E}.json › wall_time_s; runs/measured/E.json › novel_classes (q_floor, floor_above_alpha) -->

## 7. Release readiness of GIPPA 1.1.0

The readiness check was computed on the release code `e43683c7` at 2026-11-12T05:50 UTC and reported
overall **pass** (Table R10). Of the ten acceptance rows:
- five passed under readiness's rules (R1, R3, R4, R6, R9);
- one (R8) is reported only and measured about 7× its estimate;
- three were waived;
- one was deferred.

The three waivers cover:
- a mapping FDR that was not measured (R2);
- the failing TE4 pair (R5);
- a runtime outside the specified range (R7).

A readiness pass means only that the row's check succeeded; for R3 and R4, "no *unexpected*
failure". Specification §13.2 reads "Every statistical method passes its calibration
suite … and every attack null stays at level", and §8.5 says that a failing attack blocks the merge.
That wording is not met. Behind the two passes lie:
- known failures: 6 xfailed tests and 13 merged strict-xfail records (Tables R3a and R3b), which the
  calibration merge counts as non-failures (ADR-0057);
- two methods that pass only through their flags (ORA, determinants);
- GSEA unflagged rows that pass only under a 1.5α screen (Python port engine);
- recorded-only or vacuous results that would fail an α-level reading (Section 1.2): 4 autocorrelation
  records under UniOP-like and UniOP K-12 TU errors, 2 island breaches outside the regime, the
  gene-level H1 p (0.1737 at α = 0.05), the de novo truth-table bias (−0.066), the vacuous
  missed-operon pass, 3 unconfirmed streams whose pooled value exceeds the bound, the determinants
  ΔR² interval (0.388 miss rate) and distance-fallback CIs (0.074 / 0.072), and the `two_peptide`
  gene FDP (0.054 at q = 0.01).

<!-- src: validation/release/1.1.0/READINESS.md:9-18; scripts/release/readiness.py:373 (R4), R3 threshold "exit 0 and 0 failed for R or an ancestor C with reusable(C, R)"; GIPPA_SPEC.md §13.2, §8.5:1991;
docs/adr/ADR-0057-calibration-merge-honours-strict-xfail.md:1-8; benchmarks/baseline.json › standard_run.wheat.k12_core.ratio_to_estimate 7.249 -->

**Table R10.** Readiness 1.1.0.

| Row | Criterion | Status | Value (abridged) | Where the evidence was produced |
|---|---|---|---|---|
| R1 | ≥ 99.5 % of A2 peptides map exactly | pass | 100.00 % | 1.0.0-era pilot |
| R2 | decoy-estimated mapping FDR ≤ 1 % | **waived**ᵃ | "no mapping_fdr in validation/pilot/metrics/metrics.json and no --freeze-run" | 1.0.0-era pilot |
| R3 | every method passes its calibration suite | readiness pass (no waiver)ᵇ | full tier: 15 suites, 0 failed; 555 tests: 459 passed, 90 skipped, 6 xfailed, 0 failed; fast tier: 28,937 passed, 148 skipped, 17 xfailed, 3 xpassed; coverage 97 % | full tier at `4a6fa998` (reused); fast tier at `e43683c7` |
| R4 | every attack null stays at level | readiness passᵇ | 14 suites, 0 failed; 281 attack tests: 273 passed, 3 skipped, 5 xfailed | `4a6fa998` (reused) |
| R5 | every true-example assertion the data allow passes | **waived**ᵃ (TE4 only) | TE1, TE2, TE3, TE5, TE7, TE8 pass; TE4 fail | pilot 1.0.0-era; TE8 release code |
| R6 | same seed gives byte-identical tables | pass | TE7 pass; integration 967 passed; nf-test 144 tests | TE7 1.0.0-era; integration release code |
| R7 | standard run ≈ 35 min (16–55) | **waived**ᵃ | 3:41 > 0:55 (13,232 s); wheat | 1.0.0-era benchmark |
| R8 | core modules ≈ 25 min | reported only (readiness marks it pass) | 3:01 (10,874 s), about 7× the estimate; wheat | 1.0.0-era benchmark |
| R9 | the pilot passes (§15.6) | pass | 5 bullets; not passed: 0 | 1.0.0-era pilot |
| R10 | co-abundance (moves to the full validation) | **deferred** | by specification | none |
| H1 | milestone gates M1–M10 recorded in PROGRESS.md | pass | 10 gate blocks | 1.0.0-era gate records |
| H2 | version agreement; no release tag yet | pass | 1.1.0 everywhere; no tag v1.1.0 | release code |
| H3 | changelog built | pass | built | release code |
| H4 | containers pinned by digest | pass | 17 digests; 3 conda-only or in-process | release code |
| H5 | licences listed; GIPPA is MIT | pass | all tools listed; MIT | release code |
| H6 | CI check runs on the commit | pass | 11 check runs succeeded; summary, deploy, guards skipped | release code |

ᵃ The readiness file applies the 1.0.0 waiver records ADR-0071 (R2), ADR-0070 (R5) and ADR-0072 (R7).
The waiver mechanism does not check a record's version. ADR-0076 is the decision that extends these
records to 1.1.0, for R5 only for TE4's rplL/rplJ pair.
ᵇ Readiness checks exit 0 and "0 failed" for the release commit or a reusable ancestor; strict
xfails do not count as failures. The §13.2 wording is not met: see Tables R3a and R3b (6 xfailed tests,
13 merged records), the flagged methods (ORA, determinants) and the GSEA 1.5α screen (Python port engine).

<!-- src: validation/release/1.1.0/READINESS.md:1-24; validation/release/1.1.0/readiness.json › rows[*].status, value, evidence, waiver (ADR-0071, ADR-0070, ADR-0072);
docs/adr/ADR-0076-waivers-1-1.md:10-14 (apply_waivers does not read the version), 28-33 -->

Three rows need a precise reading.
- **R5** reads waived as a whole row, because the waiver mechanism works on whole rows. TE4 stays
  reported as a fail, and TE2's pass within R5 is vacuous (Section 3).
- **R9** passed because PROGRESS.md has a Khoa review line under Gate M10 and no REPORT.md status
  starts with "fail" or "pending". The specification's fifth bullet reads "TE2–TE5 are reported. A
  failure is investigated and explained in the report before Khoa's review." Its recorded status is
  "reported (1 not passed, explained)", for TE4. Here "explained" means investigated and reported:
  the cause of the TE4 failure was not established (Section 3).
- **H1** checks only that ten gate headings (`### Gate M1` … `M10`) exist in PROGRESS.md. It does not
  check whether each gate was met.

The gate records were compiled from the repository and logs on 1 November 2026 (1.0.0 freeze) and
are not reviewed verdicts. Only M9 and M10 carry a recorded review acceptance. The
records give:
- **met:** M1, M4, M5, M9 and M10;
- **met with exceptions:** M2, M6 and M8;
- **not met:** M3 and M7, which were shipped in 1.0.0 as disclosed known limits by a decision of
  2 November 2026.

H1 therefore reads pass although gates M3 and M7 were not met.

<!-- src: docs/adr/ADR-0076-waivers-1-1.md:28-33; scripts/release/readiness.py:470, 478, 494 (R9: Khoa review line under Gate M10; fail/pending statuses); GIPPA_SPEC.md:2694; validation/pilot/metrics/metrics.json › checklist[4];
validation/release/1.1.0/READINESS.md:19 (H1 threshold "### Gate M1 … M10"); PROGRESS.md:127-131, 133-288 (gate blocks; M9 line 226, M10 line 237) -->

Readiness 1.0.0 (commit `1aae3b67`, 3 November 2026) differs from 1.1.0 as follows:
- **R3** changed from waived to pass, because the fast-tier calibration ran on the 1.1.0 release
  commit, whereas for 1.0.0 it had not (ADR-0073).
- **R4** has the same counts (281 / 273 / 3 / 5), re-measured at `4a6fa998` ("array 1966680", jobs 1966680–1966694).
  Readiness 1.0.0 used array 1924911, reused from `1099769`.
- **R2, R5 and R7** remained waived (pilot and benchmark not re-run); R5 and R7 values are unchanged. R5's TE8 evidence moved from the weekly run 1914510_0 at `fa49b8f1` to
  1972755_0 at `e43683c7`. R2 is still unmeasured, but its value text changed: 1.0.0 read "missing:
  <freeze-run>/mapping/mapping_fdr.tsv", and 1.1.0 reads it from `metrics.json`.
- **R6:** the integration count rose from 964 to 967 passed tests, and nf-test from 130 to 144 tests.

<!-- src: validation/release/1.0.0/READINESS.md:9-24 (R2 line 10); validation/release/1.1.0/READINESS.md:10,13; docs/adr/ADR-0073-r3-nightly-calibration-waiver-1-0-0.md -->

The nightly runs around the release also produced failures and non-passing exits:
- **Time limits.** Of the seven fast-tier calibration rows from 9 to 12 November, four ended at the
  shard time limit (exit 143). The backlog attributes this to slow nodes with 3 shards and a 6-h
  limit; 4 shards with a 10-h limit completed.
- **Test failures.** Two rows at development commit `e6f4aa3c` each recorded one failing test (a
  report-render property test in 1956100_2, a plot-style test in 1965524_2); the records do not show
  whether they were investigated. Full passing runs followed on
  `20283834` and `e43683c7`.
- **Integration flake.** The release-commit integration run failed once on a known test flake (the
  fake container log records of example 21 interleave under concurrent tasks) and passed on rerun.
- **Islands task.** The first islands full-tier task ended with exit 143 and passed on rerun
  (1968544).

<!-- src: local/logs/nightly.tsv:115,119,123,130,134,137,141 (calibration_nightly rows), 135-136 (integration);
local/ci/runs/1956100_2-s2/calibration_nightly.txt (FAILED tests/unit/test_m11_03_render.py::test_property_overall); local/ci/runs/1965524_2-s3/calibration_nightly.txt (FAILED tests/unit/viz/test_m12_04_pr5_locus_style.py::test_any_status_mix_gives_the_right_fill_per_gene);
docs/backlog.md:1470-1475; CAL/failed-1966685/task-5.json -->

## Figures

The outline's figures are listed with the status of their source. All pilot figures exist only as
1.0.0-era renderings from the freeze runs (`$GIPPA_ROOT/runs/pilot/freeze/fa49b8f/`). Since then 24
files in `src/gippa/viz/` changed, including V13, which now draws every candidate, so the figures
need re-rendering or a label. Which pilot figures may appear is still an open question (Q-K11).

- **Figure 1** (new workflow schematic): not available; it is still to be drawn (outline: pending,
  M11-08).
- **Figure 2** (V01 `circos`, from the A2 freeze run `out_a2/plots/V01_circos`): where the A2
  proteome lies along the Nissle 1917 genome; descriptive.
- **Figure 3** (new panels, from the CAL harness records): a QQ plot of null p-values and observed FDP
  against nominal q. This figure **needs new plotting code**; no code for these panels exists.
- **Figure 4** (V07 `obs_vs_exp` and V08 `cog_abundance`, from the A1 and A2 freeze runs
  `out_a1/plots` and `out_a2/plots`): detection against detectability, by region and functional
  category.
- **Figure 5** (V11 `operon_pairs`, from the A2 freeze run `out_a2/plots/V11_operon_pairs`): operon
  partners on the pilot. This figure is descriptive and is not evidence for the simulation
  calibration of operon tests.
- **Figure 6** (V13 `novel_orfs`, from the outline's A2 and A3 freeze runs `out_a2/plots` and
  `out_a3/plots`): novel ORFs by class and call status. A2 has no novel-ORF rows, so its panel would
  be empty. A3 has 5 rows, none of them called. A1 has the only called novel ORF and its own V13
  rendering (`out_a1/plots/V13_novel_orfs`), so it is the candidate source; the choice needs review.
- **Supplementary Figure S1** (V16 `compare`, from `cmp_a1_a2/` and `cmp_a1_a2_a3/` of the freeze
  runs): detection and abundance agreement between arms.

<!-- src: docs/preprint/outline.md §3; directory listing of local/runs/pilot/freeze/fa49b8f/{out_a1,out_a2,out_a3}/plots and cmp_a1_a2*/; out_a3/run_manifest.json › git_commit;
validation/pilot/metrics/metrics.json › m6.A2.novel_orfs.n_rows (0), m6.A3.novel_orfs.totals; git diff --stat fa49b8f1 e43683c7 -- src/gippa/viz (24 files); PROGRESS.md:38 (M12-05a V13) -->
