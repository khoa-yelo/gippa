# Results

<!-- Draft results section for the GIPPA 1.1.0 preprint. Every number is copied from the file named in
the `src` comment after its paragraph or table. Paths are relative to the repository root.
Abbreviation used in the comments: CAL = local/logs/calibration/full-4a6fa998300c-20261109-160014
(the release-lineage calibration_full run). Pilot rows are aggregates only (no peptide sequence, no
DIA-NN row, no per-gene value except the named TE4 genes). -->

## Provenance of the results

The results below come from four code states, and each subsection says which one it used.
Calibration on simulated data, the attack nulls, the integration and nf-test runs, the K-12 weekly
run (TE6, TE8) and the housekeeping rows of the readiness check were run on the 1.1.0 release
code. The calibration evidence was produced at commit `4a6fa998`; between that commit and the
release commit `e43683c7` only the version string in `src/gippa/__init__.py` changed in `src/`,
`workflow/` and `tests/calibration/`. The Nissle 1917 pilot, the true-example values TE1–TE5 and
TE7, and the standard-run benchmark were measured on the 1.0.0 report-freeze commit `fa49b8f1`
(30 October 2026). Between `fa49b8f1` and `e43683c7`, 124 files in `src/` and `workflow/` changed
(+6,484/−1,090 lines). The 1.1.0 readiness check reuses these 1.0.0-era results under ADR-0076,
which records that "Nothing these rows measure has changed since 1.0.0" for the three waived rows
(R2, R5 for TE4 only, R7). The pilot and the benchmark were not re-run on 1.1.0 code. The example
routes and the weekly Nissle and distance-model runs used 1.1.0 release-candidate commits.

| Label | Commit(s) | Results that use it |
|---|---|---|
| release code (1.1.0) | `e43683c7`; `4a6fa998` (differs only in the version string) | calibration_full, attacks, calibration_nightly, integration, nf-test, TE6, TE8, readiness H1–H6, `gippa demo` from the PyPI wheel |
| 1.1.0rc1-era | `39f9f3a2`, `20283834`, `e6f4aa3c` | example routes A–E; weekly Nissle TE2; K-12 distance-model fallback |
| 1.0.0-era | `fa49b8f1` | Nissle 1917 pilot (A1, A2, A3, compare); TE1–TE5, TE7; standard-run benchmark |
| older | named per row | UniOP AUROC and cut-offs (ADR-0022, HEAD `3db69d9b`); TE2 in-silico proxy (2026-10-23) |

<!-- src: validation/release/1.1.0/readiness.json › commit, computed; docs/adr/ADR-0076-waivers-1-1.md:17-24;
validation/pilot/metrics/metrics.json › gippa_commit; benchmarks/baseline.json › standard_run.wheat.k12_standard.commit;
CAL/task-0.json › commit; digest_E_results.md provenance table (git diff --stat counts) -->

## 1. Calibration on simulated data

### 1.1 Harness and run

Each statistical method was calibrated on data from the gippasim simulator against four one-sided
criteria fixed in `tests/calibration/tolerances.yaml`: the level of a p-value (share of p ≤ α no
more than α + 2 SE), uniformity (Kolmogorov–Smirnov p > 0.001), the false discovery rate (mean FDP
no more than q + 2 SE over replicates) and the family-wise error rate (share of replicates with at
least one false discovery no more than 0.05 + 2 SE). Conservative error is allowed; anti-conservative
error beyond the bound fails. The nominal levels were 0.01, 0.05 and 0.10. The full tier ran 15
suites as one Slurm array (1966680; the islands task was re-run as 1968544 after the first attempt
ended with exit 143) on commit `4a6fa998` on 9–10 November 2026. Replicates per suite were 100
(detectability, determinants, variants), 300 (placement, gene_fdr, missed_operons, gsea, novel_orfs,
denovo_precision), 500 (islands) and 1,000 (regions, autocorrelation, discordance, ora, compare).

<!-- src: tests/calibration/tolerances.yaml:1-19 (levels, se_multiplier 2.0, ks_min_p 0.001, replicates fast/full);
tests/calibration/_criteria.py:219-273; CAL/task-N.json (commit 4a6fa998, exit 0); local/logs/nightly.tsv:124 -->

All 15 tasks exited 0. The junit roll-up counted 555 tests: 459 passed, 90 skipped, 6 xfailed,
0 failed, 0 errors, 0 xpassed. Most skips are shard and confirmation placeholders whose verdict is
given by the merge step. The fast tier (calibration_nightly) on the release commit itself passed
with 28,937 passed, 148 skipped, 17 xfailed and 3 xpassed tests, with line coverage 97 % (minimum
90 %). Readiness row R3 records a pass without a waiver; for 1.0.0 the same row had been waived
(ADR-0073) because the fast tier had not run on the release commit.

<!-- src: local/logs/nightly.tsv:124 (calibration_full roll-up); local/logs/nightly.tsv:137 (1972754_2, e43683c7, calibration_nightly);
validation/release/1.1.0/READINESS.md:11 (R3); validation/release/1.0.0/READINESS.md (R3 waived, ADR-0073) -->

### 1.2 Results by method family

Table 1 summarises each suite; Table 2 gives representative observed error rates against their
bounds at α = 0.05. Every rate in Table 2 is a rejection or FDP rate over the stated R; its bound is
the nominal level plus 2 SE at that R.

**Placement FDR (mapping).** The placement suite produced 366 harness records, of which 365 passed.
At q = 0.05 the observed mean FDP was 0.000158 in the stressed database arm (R = 300) and 0.0355 in
the stressed read-through arm with shuffled decoys. The one failing record was the contaminant
attack `flat_c200` at q = 0.01 (FDP 0.01496 against a bound of 0.01019), a strict xfail whose marker
reads "contaminant weight pending (M3-25 decision 10 left it open)"; the same attack passed at
q = 0.05 and 0.10. The methods note classes the 100-kb genome-level lists as smoke checks, because
they contained no false peptides in most replicates; list-level suites and the stressed stratum-0
arm carry the FDR evidence. Gate M3 is recorded as not met (shipped as a known limit).

<!-- src: CAL/lanes/placement-group4/placement.jsonl:20; CAL/lanes/placement-group1/placement.jsonl:2;
CAL/placement.merged.jsonl:43-45; docs/methods/placement-fdr.md:776-787; PROGRESS.md:155-165 (Gate M3) -->

**Gene-level FDR (picked, MAYU).** All 132 gene_fdr records passed. At q = 0.05 with π = 0.05 the
observed mean FDP was 0.0486 for the picked FDR (bound 0.0506) and 0.0493 for MAYU with PEP input
(bound 0.0505). The methods note records that the levels 0.05 and 0.10 are weak and vacuous in
gippasim, because peptides are capped at q ≤ 0.01 and the loosest reachable gene FDR is about
0.053; it also records a bimodal FDP (mean 0.053 ± 0.003, p90 0.135) when 50 % of target q-values
are missing.

<!-- src: CAL/lanes/gene_fdr-group1/gene_fdr.jsonl:44,47; docs/methods/detection-tiers.md:207-209, 296-301 -->

**Detectability and detection limit.** The cross-fitted detectability slope had mean 0.981 over 100
replicates, inside the band [0.8, 1.2]. The detection-limit calibration-in-the-large mean was
−0.0078; this record is "recorded" and not gated.

<!-- src: CAL/lanes/detectability-group0/detectability.jsonl:7,12; tests/calibration/tolerances.yaml (detectability_slope) -->

**Region tests (circular shifts).** All 146 region records passed. The pooled Hodges–Lehmann shift
test on the chromosome at ρ = 0 rejected at 0.0505 at α = 0.05 (bound 0.0569, R = 1,000); the
within-type FDR for prophages was 0.063 (bound 0.0654).

<!-- src: CAL/lanes/regions-group2/regions.jsonl:21; CAL/lanes/regions-group0/regions.jsonl:8 -->

**Autocorrelation.** With exact transcription units the TU-block Moran's I at lag 1 rejected at
0.055 (bound 0.0638) and the BH family at 0.008. Under the distance-model operon fallback (ADR-0061),
the same statistic rejected at 0.149 / 0.354 / 0.476 at α = 0.01 / 0.05 / 0.10 (bounds 0.0163 /
0.0638 / 0.1190): a strict xfail. With UniOP's K-12 error rates applied to the TUs (recorded, not
gated) it rejected at 0.432 / 0.659 / 0.758. Of 174 records, 163 passed, 4 were unconfirmed and 7
failed (3 strict xfail, 4 recorded-only arms).

<!-- src: CAL/lanes/autocorrelation-group0/autocorrelation.jsonl:77,47,49,51; CAL/lanes/autocorrelation-group2/autocorrelation.jsonl:57,23,25,27;
digest_E_results.md §2.1 (record status counts) -->

**Abundance islands.** All 69 island records passed. The FWER at α = 0.05 was 0.022, 0.020 and
0.024 at AR(1) ρ = 0, 0.5 and 0.9 (bound 0.0695, R = 500), and 0.056 for the chromosome-plus-plasmid
step null. The FWER claim is limited to a stated dependence regime (ADR-0021 Amendment 1); outside
it, at ρ = 0.9 with short resolved lag, the methods note records breaches (for example 0.096 at
1,000 genes and 0.130 at 1,800 genes against the bound 0.0695). These breaches were measured before
the release run and are not part of the 1.1.0 full tier.

<!-- src: CAL/lanes/islands-group0/islands.jsonl:32,35,38,47; docs/methods/islands.md:166-172, 323-326; docs/adr/ADR-0021-island-oriter-block-length.md (Amendment 1) -->

**Operons.** The discordance gate held the pooled per-pair level (0.0499 at α = 0.05, bound 0.0510)
and the BH FDR (0.034, bound 0.0615) at R = 1,000. Its three failing records are strict xfails: the
step-shaped σ(abundance) attack D-1b-step, which lies outside the narrowed claim of ADR-0055
Amendment 1 (smooth, monotone σ), and the distance-model fallback (ADR-0061). For missed operons,
the convergent-only decoy arm made calls and stayed at level (FDP 0.00065 at q = 0.05), and the
engine's planted arm had FDP 0.0036. The shipped configuration (convergent plus sandwich decoys)
made no call: 26 of the suite's 48 records are `vacuous`. The missed-operon result is therefore
reported as a vacuous pass for the shipped decoys, not as evidence of power or of level.

<!-- src: CAL/lanes/discordance-group1/discordance.jsonl:17,21; CAL/lanes/missed_operons-group0/missed_operons.jsonl:35,38,44;
docs/adr/ADR-0055-operon-discordance-abundance-noise.md (Amendment 1); PROGRESS.md:188-197 (Gate M6) -->

**Function.** Over-representation analysis (ORA) is shipped uncalibrated: every ORA row carries the
flag `ora_uncalibrated`. In the full tier the fitted ORA null was at level (0.0235 over, 0.0303
under, bound 0.0556), but 10 of 72 records failed, for example the large-TU distance-model arm at
0.0841 against 0.0544 at α = 0.05; these tests pass only because of the flag. GSEA, run through R
fgsea 1.36.2, held level under TU exchangeability (single-set null 0.0415, bound 0.0525; TU null
0.0333, bound 0.0752); its two strict xfails are the restricted-claim cases of ADR-0025 and ADR-0026
(detectability-confounded sets and AR(1)-clustered sets).

<!-- src: CAL/lanes/ora-group1/ora.jsonl:67,70,42; CAL/lanes/gsea-group0/gsea.jsonl:45; CAL/lanes/gsea-group1/gsea.jsonl:52;
docs/methods/function.md:316,561,645; CAL junit (task-9 xfails) -->

**Determinants.** The determinants records in the full tier passed (for example LightGBM and elastic
net null level 0.02 at α = 0.05, bound 0.0936, R = 100). The permutation p-value nevertheless ships
flagged `null_uncalibrated` (ADR-0034), because it fails the dependence and detection-selection
attacks (for LightGBM, 0.12, 0.95 and 0.53 at α = 0.05 in three attack arms, measured at R = 100
before the release); the attack test records these rates and does not gate on flagged results. The
full-tier pass is therefore not a level claim for real proteomes.

<!-- src: CAL/determinants.merged.jsonl:2,5,17; docs/adr/ADR-0034-determinants-null-uncalibrated.md (Context, Decision); docs/methods/determinants.md:100-112 -->

**Novel ORFs and variants.** Class-level FDR for novel ORFs was at level (picked 0.0453, de novo
0.0498 at q = 0.05, R = 300), and the pooled novel-ORF F1 was 0.938 (SE 0.0014) against the minimum
0.90. Nine novel-ORF records failed, all strict xfails: the multi-peptide false-ORF attack at all
three levels (marker "M7-11 decision 3: ORF q = min peptide q"), the GC-poor intergenic attack at
q = 0.01, and five cells with exact systematic de novo errors (ADR-0050; FDP up to 0.0749 at
q = 0.01). For variants, the pooled F1 of `high` ALT calls was 0.874 (SE 0.0018) against 0.90, and
the random and systematic isobaric-swap attacks gave FWER 0.97 and 1.0 against a bound of 0.094; all
three are strict xfails. Gate M7 is recorded as not met and was shipped as a known limit. A higher F1,
0.933 [0.924, 0.941], was measured under ADR-0049's narrowed scope on code of 23 October and has no
full-tier record.

<!-- src: CAL/novel_orfs.merged.jsonl:38,41,92,22-24,43,117,120-121,123-124; CAL/variants.merged.jsonl:1-3;
PROGRESS.md:199-212 (Gate M7); docs/adr/ADR-0049-m7-variant-precision.md; docs/adr/ADR-0050-denovo-orf-q-exact-errors.md -->

**Compare, entrapment and pilot-hypothesis attacks.** The TU-block bootstrap CI for Spearman ρ in
`gippa compare` covered at 0.952 (non-coverage 0.048, bound 0.0638, R = 1,000; mean width 0.201).
A gene-level bootstrap, recorded for contrast, covered at 0.787. The TU-level H1 test held level
under a uniform null (0.0449, bound 0.0544, R = 10,000), while the gene-level p of H1, recorded only,
was anti-conservative (0.0948 / 0.1737 / 0.227 at α = 0.01 / 0.05 / 0.10). The H4 false-confirmation
FWER was 0.047 (bound 0.0638), and the free-entrapment ucb95 covered at 0.9736 (bound 0.9469,
R = 20,000) in arm A; a recorded arm E covered at 0.71875. The de novo precision suite produced
diagnostics only: the bias against realised precision ranged from 0.0015 to 0.0094 across 10 bins
(tolerance 0.05).

<!-- src: CAL/lanes/compare-group0/compare.jsonl:41,27,38,13; CAL/lanes/compare-group0/compare_coverage.jsonl:1-2;
CAL/lanes/compare-group0/compare_attack_h1_recorded.jsonl:1-3; CAL/lanes/compare-group0/compare_attack_free_entrapment_recorded.jsonl:4;
CAL/lanes/denovo_precision-group0/placement_m3_21_diagnostics.jsonl:1 › bias_vs_realised -->

**Table 1.** Full-tier calibration by suite (release lineage, commit `4a6fa998`, array 1966680).
"Records" are harness records (criterion, level, observed, bound); the status of the last block is
counted.

| Suite | R | Tests (pass / skip / xfail) | Harness records | Status for the paper |
|---|---|---|---|---|
| placement | 300 | 87 / 1 / 0 | 366: 365 pass, 1 fail (strict xfail) | pass at list level; `flat_c200` at q 0.01 flagged; gate M3 not met |
| gene_fdr | 300 | 42 / 0 / 0 | 132 pass | pass; levels 0.05/0.10 weak or vacuous |
| detectability | 100 | 11 / 0 / 0 | 28: 27 pass, 1 recorded | pass (slope); detection-limit CITL recorded |
| regions | 1,000 | 50 / 0 / 0 | 146 pass | pass |
| autocorrelation | 1,000 | 43 / 0 / 1 | 174: 163 pass, 4 unconfirmed, 7 fail (3 strict xfail, 4 recorded-only) | pass for exact TUs; uncalibrated under the distance-model fallback |
| islands | 500 | 29 / 0 / 0 | 69 pass | pass within the stated dependence regime |
| discordance | 1,000 | 28 / 0 / 2 | 39: 36 pass, 3 fail (strict xfail) | narrowed claim (smooth, monotone σ) |
| missed_operons | 300 | 27 / 0 / 0 | 48: 22 pass, 26 vacuous | vacuous for the shipped decoys |
| ora | 1,000 | 10 / 0 / 0 | 72: 61 pass, 10 fail, 1 unconfirmed | flagged `ora_uncalibrated` |
| gsea | 300 | 55 / 0 / 2 | 111 pass | pass under TU exchangeability; restricted claim |
| determinants | 100 | 21 / 0 / 0 | 21 pass | flagged `null_uncalibrated` (fails attacks, ADR-0034) |
| novel_orfs | 300 | 21 / 86 / 0 | 125: 116 pass, 9 fail (strict xfail) | class FDR and F1 pass; 9 strict xfails |
| denovo_precision | 300 | 2 / 0 / 0 | diagnostics only | diagnostic |
| variants | 100 | 7 / 3 / 1 | 3: 3 fail (strict xfail) | fails F1 ≥ 0.90; gate M7 not met |
| compare | 1,000 (ρ CI); 20,000 and 10,000 (attacks) | 26 / 0 / 0 | 41: 40 attack records ok, 1 pass | pass; gene-level CI under-covers (recorded) |
| **total** | | **459 / 90 / 6** (555 tests) | | |

<!-- src: CAL/task-0.json … task-14.json; CAL/task-N.*.junit.xml; CAL/<suite>.merged.jsonl and CAL/lanes/<suite>-group*/*.jsonl
(status counts as tabulated in digest_E_results.md §2.1; junit totals re-checked: 459/90/6 over 555 unique tests); local/logs/nightly.tsv:124 -->

**Table 2.** Observed error against nominal level, representative gated streams (α or q = 0.05
unless stated; release lineage).

| Suite | Stream | Criterion | Observed | Bound | R |
|---|---|---|---|---|---|
| placement | sim_db_stressed.pooled | FDR | 0.000158 | 0.0500 | 300 |
| placement | readthrough_stressed.shuffle | FDR | 0.0355 | 0.0510 | 300 |
| gene_fdr | picked, π 0.05 | FDR | 0.0486 | 0.0506 | 300 |
| gene_fdr | MAYU (PEP), π 0.05 | FDR | 0.0493 | 0.0505 | 300 |
| regions | HL shift, ρ 0, chromosome pooled | level | 0.0505 | 0.0569 | 1,000 |
| regions | prophage HL shift | FDR | 0.063 | 0.0654 | 1,000 |
| autocorrelation | Moran's I TU-block lag 1 (exact TUs) | level | 0.055 | 0.0638 | 1,000 |
| autocorrelation | same, distance-model TUs (strict xfail) | level | 0.354 | 0.0638 | 1,000 |
| islands | ρ = 0 / 0.5 / 0.9 | FWER | 0.022 / 0.020 / 0.024 | 0.0695 | 500 |
| discordance | operon null, 1,800 genes | level / FDR | 0.0499 / 0.034 | 0.0510 / 0.0615 | 1,000 |
| missed_operons | convergent-only decoys | FDR (counts) | 0.00065 | 0.0504 | 300 |
| missed_operons | shipped decoys, power arm | FDR (counts) | 0.0 (vacuous) | none | 300 |
| ora | fitted null, over / under | level | 0.0235 / 0.0303 | 0.0556 | 1,000 |
| ora | distance-model large-TU under (fails; flag) | level | 0.0841 | 0.0544 | 1,000 |
| gsea | single-set null (R engine) | level | 0.0415 | 0.0525 | 300 |
| determinants | null, LightGBM / elastic net | level | 0.02 / 0.02 | 0.0936 | 100 |
| novel_orfs | class FDR, picked / de novo | FDR | 0.0453 / 0.0498 | 0.0511 / 0.0507 | 300 |
| novel_orfs | pooled F1 (≥ 0.90) | F1 | 0.938 | 0.90 | 300 |
| variants | pooled F1 of `high` calls (≥ 0.90; strict xfail) | F1 | 0.874 | 0.90 | 100 |
| compare | ρ CI non-coverage (TU block) | level | 0.048 | 0.0638 | 1,000 |
| compare | H4 false confirmation | FWER | 0.047 | 0.0638 | 1,000 |

<!-- src: CAL/lanes/placement-group4/placement.jsonl:20; CAL/lanes/placement-group1/placement.jsonl:2; CAL/lanes/gene_fdr-group1/gene_fdr.jsonl:44,47;
CAL/lanes/regions-group2/regions.jsonl:21; CAL/lanes/regions-group0/regions.jsonl:8; CAL/lanes/autocorrelation-group0/autocorrelation.jsonl:77;
CAL/lanes/autocorrelation-group2/autocorrelation.jsonl:25; CAL/lanes/islands-group0/islands.jsonl:32,35,38; CAL/lanes/discordance-group1/discordance.jsonl:17,21;
CAL/lanes/missed_operons-group0/missed_operons.jsonl:38,35; CAL/lanes/ora-group1/ora.jsonl:67,70,42; CAL/lanes/gsea-group0/gsea.jsonl:45;
CAL/determinants.merged.jsonl:2,5; CAL/novel_orfs.merged.jsonl:38,41,92; CAL/variants.merged.jsonl:3; CAL/lanes/compare-group0/compare.jsonl:41,38 -->

### 1.3 Known failures (strict xfails)

A strict xfail is a test that is expected to fail and turns the run red if it passes. The junit
roll-up that readiness row R3 reads counted 6 xfailed tests (Table 3a). The merge step of sharded
suites records further strict-xfail verdicts in the merged records; these are not in the roll-up
count (Table 3b). Three of them (`flat_c200`, the GC-poor intergenic attack and the multi-peptide
false-ORF attack) appear in the junit files as skipped attack tests, whose verdict is the merge
step. Readiness R3 therefore reads "6 xfailed" although 13 further merged records (1 placement,
9 novel-ORF, 3 variant) are known failures.

**Table 3a.** Strict xfails counted by the junit roll-up (6).

| Suite | Test | Reason recorded in the marker |
|---|---|---|
| autocorrelation | distance-model Moran's I lag 1 (arm B) | ADR-0061: uncalibrated under the distance-model fallback |
| discordance | distance-model p uniformity (arm A) | ADR-0061 |
| discordance | D-1b-step | ADR-0055 Amendment 1: step σ(abundance) outside the smooth-monotone claim |
| gsea | detectability-confounded sets | ADR-0025/ADR-0026: no all-rows level claim |
| gsea | AR(1)-clustered sets (`ar1_run15`) | ADR-0025/ADR-0026 |
| variants | A6a isobaric rows flagged | ADR-0049: isobaric blocks aligned as substitutions |

<!-- src: CAL/task-4.*.junit.xml, task-6.*.junit.xml, task-9.*.junit.xml, task-13.*.junit.xml (xfail entries); local/logs/nightly.tsv:124 -->

**Table 3b.** Strict-xfail verdicts in merged records that the roll-up does not count.

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

The attack nulls are simulations built to make a method reject falsely: autocorrelated abundance
without islands, detectability confounded with region membership, gene sets packed into one operon,
high-copy plasmids, de novo errors that reverse a stretch, isobaric swaps everywhere, and further
method-specific attacks (55 attack test files). In the release-lineage full tier, attack tests ran in
14 suites. All 14 suites passed: 281 attack tests, of which 273 passed, 3 were skipped and 5 were
xfailed, with 0 failed, 0 errors and 0 xpassed. The 5 xfails are the five attack tests in Table 3a
(two distance-model fallback tests, D-1b-step, two GSEA restricted-claim tests). The 3 skipped
attack tests are the ones whose verdict is a merged strict xfail (`flat_c200`, the GC-poor
intergenic attack and the multi-peptide false-ORF attack; Table 3b). "14 of 14 suites passed"
thus means that no attack failed outside the known failures listed in Section 1.3, and that
flagged methods (ORA, determinants) pass their attacks only through their flags.

<!-- src: local/logs/nightly.tsv:125; validation/release/1.1.0/READINESS.md:12 (R4); GIPPA_SPEC.md:1980-1991 (attack list);
tests/calibration/attacks/ (55 files); CAL/task-*.junit.xml (attack test outcomes); docs/adr/ADR-0034-determinants-null-uncalibrated.md (Decision) -->

## 3. True-example assertions and operon calibration on real genomes

The specification defines eight true-example assertions (TE1–TE8) on real genomes. Readiness row R5
lists seven of them and records it as **waived**, because TE4 failed (Table 4). TE6 was measured in
the weekly K-12 run but is not a readiness row.

**Table 4.** True-example assertions.

| TE | Assertion (short) | Outcome | Value | Code state |
|---|---|---|---|---|
| TE1 | ≥ 99.5 % of A2 peptides map exactly in frame | pass | 28,129 of 28,129 (1.0; Wilson 95 % lower bound 0.99986) | 1.0.0-era pilot |
| TE2 | no `high` ALT call off start or stop codon sites | pass, **vacuous** | 0 ALT sites, 0 `high`, 0 violations: it checked nothing | 1.0.0-era pilot; 1.1.0rc1 weekly Nissle identical (0 sites) |
| TE3 | ribosomal proteins ≥ 90 % detected, median abundance in the top quartile | pass | 53 of 53 detected (1.0 [0.932, 1.0]); median percentile rank 0.987 (iBAQ) | 1.0.0-era pilot |
| TE4 | log2(RplL/RplJ) ∈ [1, 3]; log2(AtpA/AtpD) ∈ [−1, 1] | **fail (waived, ADR-0070 via ADR-0076)** | rplL/rplJ −0.480 (Top3 −0.531); atpA/atpD −0.118 (pass) | 1.0.0-era pilot |
| TE5 | trypsin and keratin peptides land in the contaminant bucket | pass | A1 59 of 59; A2 61 of 61 | 1.0.0-era pilot |
| TE6 | K-12: Bakta oriC next to mnmG | measured, not a readiness row | "2 oriC calls, at the true origin: 1" | release code (weekly 1972755_0) |
| TE7 | same seed gives byte-identical tables | pass (GIPPA tables only) | 2,786 gate files compared, 0 differ; DIA-NN layer not run | 1.0.0-era pilot |
| TE8 | UniProt K-12 peptides ≥ 99.5 % exact on K-12; Nissle placement FDR ≤ 1 % | pass | K-12 exact 0.99900 (151,272 of 151,424); Nissle placement FDR 0.00533 | release code (weekly 1972755_0) |

<!-- src: validation/release/1.1.0/READINESS.md:13 (R5); validation/pilot/metrics/metrics.json › te1, te_real.te2, te_real.te3, te_real.te4.pairs, te_real.te5, te7.gippa;
local/ci/validation/2026-11-10_te2_1969418.json › counts; local/logs/gippa-weekly.1972755_0.out:319-320;
local/ci/validation/2026-11-11_te8_1972755.json › counts, shares, passed -->

TE1 was evaluated on the 28,129 non-contaminant A2 peptides; 18 exact in-frame peptides were masked,
42 started at non-ATG N-termini and no peptide missed. TE2 passed in the pilot with zero ALT sites
in A2 (`vacuous: true`) and again with zero sites in the weekly Nissle run at `20283834`, whose note
records that the fdnG selenocysteine peptides were rejected by the stratum-1 placement FDR (4
targets, no decoys, q = 0.25) and that fdhF has no tryptic peptide over its Sec codon. An earlier
in-silico proxy (2026-10-23, before the stop-codon read-through changes) had 233 ALT sites and 80
`high` calls, all at start codons, with 0 violations; it predates the current code. For TE4, the
waiver review recorded an iBAQ SD of 0.99 among 48 ribosomal proteins expected at 1:1 and 53 of
2,766 A2 genes with a peptide that swings by more than 3 log2 between runs; the cause was not proven.
TE7 compared two GIPPA runs on the same DIA-NN output at the freeze commit, so reproducibility from
raw spectra to tables was not shown. Readiness row R6 adds the release-code integration run (967
passed, 84 skipped, 1 xfailed, 1 xpassed) and nf-test (144 tests). TE8 also reported that on the
Nissle genome 0.764 of K-12 peptides mapped exactly, 0.0895 mapped tolerantly and 0.147 were unmapped.

<!-- src: validation/pilot/metrics/metrics.json › te1 (exact_inframe_masked 18, non_atg_nterm.peptides 42, n_misses 0), te_real.te2.vacuous;
local/ci/validation/2026-11-10_te2_1969418.json › notes; validation/proxies/2026-10-23_te2.json › counts;
docs/adr/ADR-0070-te4-rplL-rplJ-waiver-1-0-0.md:39,51; validation/pilot/REPORT.md (Notes: TE7, Diagnosis: TE4);
local/logs/nightly.tsv:136; local/ci/validation/2026-11-11_te8_1972755.json › shares -->

**Operon predictor on K-12.** The UniOP AUROC against 2,739 labelled RegulonDB 14.5.0 pairs (1,730
positive, 1,009 negative) was 0.911, with a directon-bootstrap 95 % percentile CI of [0.898, 0.923].
This is a number for K-12 only, not a transfer claim, and it may be partly in-sample, because
UniOP's training data were not verified. At the shipped pair threshold of 0.5, precision was 0.896
and recall 0.853 (F1 0.8739); precision depends on prevalence (K-12 labelled π = 0.632). The
intergenic-distance model, the fallback when UniOP is missing, had AUROC 0.907 [0.892, 0.920];
UniOP exceeded it by 0.0041 [0.0009, 0.0077]. On the labelled pairs, the distance-model pair rule
at p ≥ 0.529415 (d ≤ 97.1 bp) called 1,878 pairs with precision 0.839 and recall 0.910. These
numbers come from Slurm job 1879748 on 24 October 2026, before the 1.0.0 freeze.

<!-- src: docs/adr/ADR-0022-uniop-terminatornet-cutoffs.md:78-136 (AUROC table line 84, distance-model table lines 120-136, decision on 0.5);
digest_E_results.md §4 (job 1879748, HEAD 3db69d9b) -->

A different quantity was measured end to end in the weekly distance-model run (`20283834`): the TU
pairs produced by a full GIPPA run with UniOP disabled, compared with the RegulonDB TU pairs, gave
precision 0.708 and recall 0.557 (TP 963, FP 397, FN 767). This is pipeline-level TU-pair agreement,
not the pair-classification precision of ADR-0022, and the two are not comparable. The matching
UniOP comparison in the same weekly run was "not available".

<!-- src: local/logs/gippa-weekly.1969418_2.out:116; local/logs/weekly.tsv:9,12 -->

## 4. The Nissle 1917 pilot (PXD078459)

All pilot results in this section were produced at the 1.0.0 report-freeze commit `fa49b8f1` and are
reused for 1.1.0 under ADR-0076; the pilot was not re-run on 1.1.0 code. The input was the
wild-type runs of PXD078459, searched with DIA-NN 2.3.1 Academia (not the 2.7.0 named in the
specification; ADR-0051). Arm A1 searched the UniProt K-12 proteome UP000000625 plus contaminants,
arm A2 the Nissle strain database plus contaminants, and arm A3 used Cascadia 0.0.7 de novo
sequencing. The freeze runs used the pass-3 wild-type genome and its recorded Bakta files; the A3
Cascadia table was the pass-2 conversion, which was not re-run. The free-entrapment rerun that
counts ran from main commit `db3a95f1`, whose `src/` and `workflow/` equal those of `fa49b8f1`.

<!-- src: validation/pilot/REPORT.md:1-30 (status, data and arms), 306-354 (deviations); validation/pilot/metrics/metrics.json › gippa_commit;
docs/adr/ADR-0076-waivers-1-1.md:17-24 -->

### 4.1 End-to-end runs, identifications and mapping

GIPPA completed on all three arms at the freeze commit (A2 job 1914948, A1 1915264, A3 1916086, all
exit 0). Table 5 gives identification and mapping aggregates. Only the precursor counts are at a
controlled 1 % FDR (global precursor q); peptide and protein-group counts are not FDR-controlled.

**Table 5.** Identifications and mapping, A1 and A2 (union over runs).

| Quantity | A1 (K-12 reference) | A2 (strain database) |
|---|---|---|
| precursors (global q ≤ 0.01) | 28,639 | 33,764 |
| peptides: all / non-contaminant / contaminant-only | 23,792 / 23,704 / 88 | 28,234 / 28,129 / 105 |
| protein groups / with ≥ 2 peptides | 2,798 / 2,433 | 3,228 / 2,790 |
| exact mapping share [Wilson 95 %] | 0.9936 [0.9925, 0.9945] | 1.0 [0.99986, 1.0] |
| tolerant-only share | 0.00515 | 0 |
| unmapped share | 0.00127 | 0 |
| genome-inconsistent peptides (rate, count) | 0.058, 1,450 | not applicable |

<!-- src: validation/pilot/metrics/metrics.json › m1.A1/A2.union (precursors.all, peptides, protein_groups.all, protein_groups_min2.all);
m2.A1/A2.shares; m2.A1.genome_inconsistent; validation/pilot/REPORT.md:11 (job IDs) -->

The decoy-estimated mapping FDR of A2 (readiness R2) was not measured. All 28,129 A2 peptides with
an exact target hit come from the same database as the search (`same_db`) and are exempt from both
decoy lists, so no accepted A2 peptide was decoy-estimated; the row is waived (ADR-0071, carried by
ADR-0076). A context-only value on A1 was 7.27e-05 over 23,675 accepted peptides; the waiver records
that this value cannot fail by construction and does not validate A2. Identification error in A2 is
inherited from DIA-NN (input q ≤ 0.01) and is not estimated by GIPPA.

<!-- src: docs/adr/ADR-0071-r2-mapping-fdr-waiver-1-0-0.md:15-17, 30-41; validation/release/1.1.0/READINESS.md:10 (R2) -->

### 4.2 Detection and detection tiers

Of the 4,685 strain CDS, A1 detected 2,673 (0.5705 [0.5563, 0.5847]) and A2 detected 3,057
(0.6525 [0.6388, 0.6660]). Of the 916 strain-specific CDS (no tryptic peptide shared with K-12), A1
detected none and A2 detected 317 (0.3461 [0.3160, 0.3775]). Shares by region are in Table 6. The
detection tier of the freeze arms is not recorded in `metrics.json` or the report; the pilot README
records for the pass-2 runs that A1 and A2 used `decoy_picked` and A3 used `two_peptide`, because
de novo rows carry no decoys.

**Table 6.** Genes detected by subset (share [Wilson 95 %]).

| Subset (genes) | A1 | A2 |
|---|---|---|
| all (4,685) | 0.5705 [0.5563, 0.5847] | 0.6525 [0.6388, 0.6660] |
| strain-specific (916) | 0 [0, 0.0042] | 0.3461 [0.3160, 0.3775] |
| plasmid (13) | 0.0769 [0.0137, 0.3331] | 0.7692 [0.4974, 0.9182] |
| prophage (136) | 0.0221 [0.0075, 0.0628] | 0.2721 [0.2043, 0.3523] |
| BGC (251) | 0.3108 [0.2567, 0.3705] | 0.5458 [0.4840, 0.6063] |
| K-12-shared (3,759) | 0.7111 [0.6964, 0.7254] | 0.7289 [0.7145, 0.7429] |

<!-- src: validation/pilot/metrics/metrics.json › m3.A1/A2.subsets.{all,strain_specific,plasmid,prophage,bgc,k12_shared}.share;
validation/pilot/README.md:434-435 (pass-2 detection tiers) -->

### 4.3 Comparison of A1 and A2

On the 2,635 genes with measured abundance in both arms (1,683 TUs), the Spearman correlation of
abundance was 0.978 with a TU-block bootstrap 95 % CI of [0.973, 0.983] (1,000 resamples). The
Bland–Altman bias was −0.195 log2 with limits of agreement [−1.376, 0.987]. The Jaccard index of
enriched functional categories was 0.563 (437 in A1, 285 in A2, 260 shared). No islands and no
operon flags were called in either arm, so their Jaccard indices are undefined. Jaccard indices for
ALT sites (118 in A1, 0 in A2) and novel ORFs (1 in A1, 0 in A2) were 0. Jaccard indices carry no
uncertainty estimate.

<!-- src: validation/pilot/metrics/metrics.json › m4 (spearman_rho, ci_low, ci_high, n_shared, n_units, n_boot, ba_*); m5.enriched_categories; m5.islands; m5.operon_flags; m5.other_families -->

### 4.4 Pre-registered hypotheses

The four hypotheses of specification §15.7 were fixed before the pilot. Table 7 gives their outcomes.

**Table 7.** Hypotheses H1–H4.

| Hypothesis | Verdict | Key numbers |
|---|---|---|
| H1: A2 detects more strain CDS than A1, mostly plasmid, prophage and BGC genes | **not supported** (failed clause: "mostly") | 417 genes detected only in A2, 33 only in A1; A2-only share 0.9267 [0.8988, 0.9473]; TU-level test 275 vs 31 units, p_tu 2.49e-50, unit effect 0.8987 [0.8598, 0.9277]; "mostly" share 106 of 417 = 0.2542 [0.2148, 0.2981] against a base rate of 0.0854 |
| H2: the strain database lowers genome-inconsistent identifications | **deferred** (full validation, FDRBench) | shown side by side only: A1 genome-inconsistent share 0.058; A1 free-entrapment estimate 0.0273 (peptide level) |
| H3: context conclusions agree between A1 and A2 | **supported** on ρ and enriched categories only; islands and operon flags not evaluable | ρ 0.978 [0.973, 0.983]; Jaccard 0.563; region tests not compared |
| H4: genome-anchored de novo precision agrees with precision against A2 | **consistent with lower bound** (the data can refute but never confirm) | 10 reportable bins, 1,364 queries; all 10 bins consistent; n-weighted mean genome precision 0.562; pooled A2-confirmed share 0.543 [0.517, 0.570] |

<!-- src: validation/pilot/metrics/metrics.json › h1 (b, c, prop_a2_only, b_tu, c_tu, p_tu, unit_effect, mostly, mostly_base_rate, failed_clauses, verdict);
validation/pilot/REPORT.md:233 (H2); metrics.json › h3.criteria, h3.region_tests, h3.verdict; metrics.json › h4 (B, n_queries, mean_genome_precision, pooled_L, verdict) -->

For H1, the gene-level p (4.65e-86) decides nothing: the compare attack in Section 1.2 recorded this
p as anti-conservative (0.1737 at α = 0.05), and the verdict rests on the TU-level test and the
"mostly" clause. H3 compared Jaccard indices without restricting them to genes detected in both
arms, which the report records as a conservative stand-in biased against H3. For H4, A2 enters both
sides of the comparison (A2's own matches set the c-value and A2's identifications set the
confirmations), so the check is not independent; the score bins span de novo scores in (0.8, 1]
only, and with DIA data there is no spectrum pairing.

<!-- src: validation/pilot/metrics/metrics.json › h1.p_gene; CAL/lanes/compare-group0/compare_attack_h1_recorded.jsonl:2; validation/pilot/REPORT.md (H1 notes; §10 deviations; H4 caveats) -->

### 4.5 Free entrapment (A1)

K-12 peptides absent from the Nissle genome served as free entrapment for A1. At precursor level
(N_T 28,347, N_E 157), the combined FDP estimate was 0.0236, the realised lower bound `fdp_lower`
0.00551, and the one-sided 95 % upper bound `ucb95` 0.0265; the assessment was **inconclusive**. At
peptide level (N_T 23,517, N_E 152) the values were 0.0273, 0.00642 and 0.0308, not assessed because
peptide-level FDR is not controlled by the search. The ratio r of K-12-only to shared in-silico
peptides was 0.307. The `ucb95` assumes equal-chance, independent false matches and is not valid
under clustered errors; the combined estimate is an estimated upper bound in expectation, not a
bound on this sample. Restricted to entrapment peptides far from any Nissle homolog (N_E 7), the
peptide-level combined estimate was 0.00274 (ucb95 0.00493).

<!-- src: validation/pilot/metrics/metrics.json › m8.precursor, m8.peptide (assessment, assessment_reason, ci.ucb95, ci.assumes), m8.universe.r, m8.peptide.sensitivity.far, m8.caveats -->

### 4.6 Novel ORFs and variants

A2 called no novel ORF and no ALT variant. A1 called one novel ORF (21 more unestimated); the report
identifies it as the first 25 codons of prfB, upstream of RF2's programmed +1 frameshift, which is
not a novel protein. A1 had 4 `high` and 114 candidate substitution ALT calls. A3 had 0 called novel
ORFs (5 not called), 1 `high` and 125 candidate swap calls, 41 candidate substitutions and 28
candidate blocks; the class error rates of A3's ALT calls are unestimated. A2 and A3 shared no ALT
site (0 in A2, 195 in A3).

<!-- src: validation/pilot/metrics/metrics.json › m6.A1/A2/A3.novel_orfs.totals, m6.*.variants.alt, m6.overlap_a2_a3.alt_sites; validation/pilot/REPORT.md:99-143 and Notes (Novel ORFs in A1) -->

### 4.7 A3 feasibility

Before arm A3 was run, a feasibility rule required at least 500 distinct 7–30-mers with de novo
score above 0.8 and a share of at least 0.20 of the top 1,000 matching the six-frame translation of
the genome. The committed verdict file covers one run, WT_1: 863 peptides passed, 549 matched, a
share of 0.636 (Wilson 95 % CI 0.604–0.668), verdict pass; the shuffled-decoy share was 0.0058. This
CI is optimistic because peptides of one protein co-vary. The committed file was evaluated on the
engineered chassis genome; the freeze used a re-evaluation of WT_1 on the wild-type genome, also a
pass with equal aggregates (not committed). The pilot README also lists WT_2 (0.644, 0.611–0.675) and
WT_3 (0.623, 0.589–0.654) as passing; only WT_1 has a committed verdict file. The check is a
precondition for A3, not a GIPPA A3 result.

<!-- src: validation/pilot/metrics/a3_feasibility.json › run, n_pass, n_match, share, share_ci95, decoy_share, rule, matches_by_replicon;
validation/pilot/README.md:443-452 (per-run L6 table), 604-606 (pass-3 re-evaluation); validation/pilot/REPORT.md:349; docs/preprint/outline.md C21 -->

### 4.8 Pilot runtimes

On the shared cluster with 4 allocated CPUs, the A1 run took 2:36:07 (sum of step wall times
10:05:49) and the A2 run 2:50:20 (12:05:36); the largest single-step peak RSS was 7.5 GiB in both
arms. These are pilot measurements on a shared cluster, not the benchmark.

<!-- src: validation/pilot/metrics/metrics.json › m7.A1/A2.wall_time_s, max_step_peak_rss_bytes, sacct.AllocCPUS; validation/pilot/REPORT.md §8 -->

## 5. Performance

**Standard-run benchmark.** The benchmark ran at the 1.0.0-era commit `fa49b8f1` on one FarmShare
"wheat" node (wheat-01, Intel Xeon Gold 6148), with 16 threads, 64 GB, profile
`farmshare_local,apptainer`, array 1916711, on 31 October 2026. Each cell ran twice, not three times
(repeat 3 was cancelled, ADR-0072); the "median" is the mean of 2 runs, and its uncertainty is
unestimated. The node was busy (median allocation 60–80 of 80 CPUs). Table 8 gives the results.

**Table 8.** Standard-run benchmark (1.0.0-era, 2 repeats per cell).

| Cell | Modules; Bakta | Repeats (s) | Mean (s) | Mean (h:mm:ss) | Spread | Peak RSS (sacct max, GiB) | Estimate (s) | Ratio to estimate |
|---|---|---|---|---|---|---|---|---|
| K-12 standard | core, mobile, bgc, function; live, full DB | 12,502; 13,963 | 13,232.5 | 3:40:32 | 0.110 | 62.35 | 2,100 | 6.30 |
| K-12 core | core; recorded JSON | 10,896; 10,851 | 10,873.5 | 3:01:13 | 0.004 | 37.87 | 1,500 | 7.25 |
| sim5m standard | core, mobile, bgc, function; live, full DB | 14,578; 12,421 | 13,499.5 | 3:44:59 | 0.160 (noisy) | 62.38 | 2,100 | 6.43 |
| sim5m core | core; recorded JSON | 12,146; 12,945 | 12,545.5 | 3:29:05 | 0.064 | 41.20 | 1,500 | 8.36 |

<!-- src: benchmarks/baseline.json › standard_run.wheat.{k12_standard,k12_core,sim5m_standard,sim5m_core} (median_s, runs_s, spread, noisy, peak_rss_sacct_max_bytes, estimate_s, ratio_to_estimate, nodes, cpu, mem_mb, rounds, slurm_jobs, date, commit);
docs/adr/ADR-0072-benchmark-baseline-1-0-0.md:25-40 (h:mm:ss means, busy node) -->

The specification estimated a standard run at about 35 min (range 16–55 min). The measured K-12
standard run took 3:41 (13,232 s), so readiness row R7 failed its range and is **waived** (ADR-0072,
carried by ADR-0076); the runtime claim is not met. Row R8 (core modules, "about 25 min", reported
only) reads 3:01 (10,874 s). Two steps, `promoter_calculator` and `terminatornet_chunk`, took about
92 % of summed task time in repeat 1 of the K-12 standard cell. Means in the h:mm:ss column are as
ADR-0072 writes them; the pilot report rounds three of them up by one second.

<!-- src: validation/release/1.1.0/READINESS.md:15-16 (R7, R8); docs/adr/ADR-0072-benchmark-baseline-1-0-0.md:39-40; validation/pilot/REPORT.md (Notes: Runtime) -->

**Other timings (release code).** `gippa demo` from the 1.1.0 PyPI wheel finished in 75 s (the
smoke script measured 77 s including start-up) on a wheat node with Python 3.12.3, and
`gippa example fetch` from the v1.1.0 GitHub release exited 0. The nightly micro-benchmark suite
passed 59 tests at `4a6fa998` (580.09 s, wheat-05) after a first run at the same commit failed two
timing budgets; a run on a post-release documentation commit failed three. The nightly end-to-end
run on the simulated 100-kb genome exited 0 with both the light Bakta database and recorded full-DB
JSON on every night from 9 to 12 November.

<!-- src: local/logs/gippa-pypi-smoke.1977437.out:3,9,14; local/logs/gippa-pypi-smoke.1977441.out:30; local/logs/nightly.tsv:120,126,139; local/logs/nightly.tsv:117,121,138 -->

## 6. Example bundle routes A–E

The example bundle `nissle_r200k` (bundle version 1.1) holds the pilot's data restricted to one
200-kb region (a separate 200-kb region for route E). Each route ran once at the 1.1.0rc1-era commit
`39f9f3a2` on 9 November 2026 with modules core and bgc, the light Bakta database, TransTermHP
terminators, Promoter Calculator off and 4 threads. For every route, all measured outputs equalled
the bundle's expected values on every key except wall time. The bundle was built at commit
`52ddd95f`, which is not an ancestor of the release commit.

**Table 9.** Example routes (one run each, 1.1.0rc1-era).

| Route | Input | Detection tier | Genes | Detected | Detected (≥ 2 peptides) | With abundance | Peptides placed | Novel ORFs (candidates / called) | Wall time (s) |
|---|---|---|---|---|---|---|---|---|---|
| A | DIA-NN report, `--search-db same` (2,334 rows) | mayu | 158 | 99 | 91 | 102 | 774 | – | 339.2 |
| B | protein FASTA and abundance table (106 records) | protein_list | 158 | 106 | 0 | 106 | 0 | – | 500.6 |
| C | Cascadia `.ssl` × 3 (76 rows; format demonstration) | two_peptide | 158 | 4 | 4 | – | 21 | – | 270.0 |
| D | generic peptide table (2,088 rows) | mayu | 158 | 99 | 91 | 102 | 774 | – | 273.8 |
| E | DIA-NN report against a six-frame database (3,650 rows) | decoy_picked | 178 | 114 | 100 | 114 | 1,041 | 8 / 0 | 641.8 |

<!-- src: local/example/routes/1966530/gippa-example/runs/measured/{A,B,C,D,E}.json (convert_rows, detection_tier, n_genes, n_detected, n_detected_2pep, n_abundance,
n_peptides_placed, n_novel_candidates, n_novel_called, wall_time_s); runs/wall_times.tsv; expected/{A..E}.json; bundle.json › bundle_version, build_commit, counts;
runs/A/run_manifest.json › git_commit; docs/adr/ADR-0074-example-dataset.md (decision 3, settings) -->

Wall times include conversion and preflight. The bundle's expected wall times (247.3, 409.0,
210.0, 188.9 and 439.5 s for A–E) are informational and were not met in this run. Route E's zero
called novel ORFs is not a negative result: in all four decoy-tagged classes the smallest reachable
q was at least 0.2, above α, so no call was possible at these class sizes.

<!-- src: local/example/routes/1966530/gippa-example/expected/{A..E}.json › wall_time_s; runs/measured/E.json › novel_classes (q_floor, floor_above_alpha) -->

## 7. Release readiness of GIPPA 1.1.0

The readiness check computed on the release commit `e43683c7` (12 November 2026) reported overall
**pass** (Table 10). Of the ten acceptance rows, six passed, three were waived and one was deferred.
The waivers cover a mapping FDR that was not measured (R2), the failing TE4 pair (R5) and a runtime
outside the specified range (R7). R8 is a reported-only row.

**Table 10.** Readiness 1.1.0.

| Row | Criterion | Status | Value (abridged) | Code state of the evidence |
|---|---|---|---|---|
| R1 | ≥ 99.5 % of A2 peptides map exactly | pass | 100.00 % | 1.0.0-era pilot |
| R2 | decoy-estimated mapping FDR ≤ 1 % | **waived** (ADR-0071 via ADR-0076) | not measured: no `mapping_fdr` | 1.0.0-era pilot |
| R3 | every method passes its calibration suite | pass (no waiver) | full tier 15 suites, 0 failed; 555 tests, 6 xfailed; fast tier 28,937 passed; coverage 97 % | release code |
| R4 | every attack null stays at level | pass | 14 suites, 0 failed; 281 attack tests, 5 xfailed | release code |
| R5 | every true-example assertion the data allow passes | **waived** (ADR-0070 via ADR-0076; TE4 only) | TE1, TE2, TE3, TE5, TE7, TE8 pass; TE4 fail | pilot 1.0.0-era; TE8 release code |
| R6 | same seed gives byte-identical tables | pass | TE7 pass; integration 967 passed; nf-test 144 tests | TE7 1.0.0-era; integration release code |
| R7 | standard run ≈ 35 min (16–55) | **waived** (ADR-0072 via ADR-0076) | 3:41 > 0:55 (13,232 s); wheat | 1.0.0-era benchmark |
| R8 | core modules ≈ 25 min | pass (reported only) | 3:01 (10,874 s); wheat | 1.0.0-era benchmark |
| R9 | the pilot passes (§15.6) | pass | 5 bullets; not passed: 0 | 1.0.0-era pilot |
| R10 | co-abundance (moves to the full validation) | **deferred** | by specification | none |
| H1 | milestone gates M1–M10 recorded in PROGRESS.md | pass | 10 gate blocks | 1.0.0-era gate records |
| H2 | version agreement; no release tag yet | pass | 1.1.0 everywhere; no tag v1.1.0 | release code |
| H3 | changelog built | pass | built | release code |
| H4 | containers pinned by digest | pass | 17 digests; 3 conda-only or in-process | release code |
| H5 | licences listed; GIPPA is MIT | pass | all tools listed; MIT | release code |
| H6 | CI check runs on the commit | pass | 11 check runs succeeded; summary, deploy, guards skipped | release code |

<!-- src: validation/release/1.1.0/READINESS.md:1-24; validation/release/1.1.0/readiness.json › rows[*].status, value, evidence, waiver -->

Three rows need a precise reading. **R5** reads waived as a whole row because the waiver mechanism
works on whole rows; ADR-0076 limits the waiver to TE4's rplL/rplJ pair, and TE4 stays reported as a
fail. TE2's pass within R5 is vacuous (Section 3). **R9** passed because all five §15.6 pilot bullets
were recorded as passed or reported; the fifth bullet ("TE2–TE5 are reported; failures explained")
reads "reported (1 not passed, explained)", referring to TE4. **H1** checks only that ten gate
headings (`### Gate M1` … `M10`) exist in PROGRESS.md. It does not check whether each gate was met.
The gate records themselves, compiled on 1 November 2026 at the 1.0.0 freeze and stated there not to
be Khoa's verdicts, record M1, M4, M5, M9 and M10 as met, M2, M6 and M8 as met with exceptions, and
M3 and M7 as **not met**; M3 and M7 were shipped as disclosed known limits by a decision of
2 November 2026. H1 therefore reads pass although gates M3 and M7 were not met.

<!-- src: docs/adr/ADR-0076-waivers-1-1.md:28-33 (scope); validation/pilot/metrics/metrics.json › checklist; validation/release/1.1.0/READINESS.md:19 (H1 threshold "### Gate M1 … M10");
PROGRESS.md:127-131 (H1 only counts headings; M3 and M7 not met), 133-288 (gate blocks) -->

Compared with readiness 1.0.0 (commit `1aae3b67`, 3 November 2026), R3 changed from waived to pass:
the fast-tier calibration ran on the 1.1.0 release commit, whereas for 1.0.0 it had not (ADR-0073).
R2, R5 and R7 remained waived with identical values, because the pilot and the benchmark were not
re-run. The integration count in R6 rose from 964 to 967 passed tests and nf-test from 130 to 144
tests. Around the release, infrastructure faults occurred and were followed by passing reruns:
calibration_nightly shards hit their time limit (exit 143) on 4 of 7 nightly calibration rows from
9 to 12 November, one integration run on the release commit failed a container-pull test before the
passing rerun cited by R6, and the first islands full-tier task ended with exit 143.

<!-- src: validation/release/1.0.0/READINESS.md:9-24; docs/adr/ADR-0073-r3-nightly-calibration-waiver-1-0-0.md; local/logs/nightly.tsv:115,119,123,130,134,137,141 (calibration_nightly rows), 135-136 (integration);
CAL/failed-1966685/task-5.json -->

## Figures

Figure placeholders are listed only where the outline's figure plan names a source that exists. All
pilot figures come from the 1.0.0-era freeze runs (`fa49b8f1`) under
`$GIPPA_ROOT/runs/pilot/freeze/fa49b8f/`; which pilot figures may appear is still an open question
(Q-K11).

- Figure 1 (new workflow schematic): not available; to be drawn (outline: pending, M11-08).
- Figure 2 (V01 `circos`, A2 freeze run `out_a2/plots/V01_circos`): where the A2 proteome lies along
  the Nissle 1917 genome; descriptive.
- Figure 3 (new panels: QQ plot of null p-values, observed FDP against nominal q; source CAL
  harness records): **needs new plotting code**; no code for these panels exists.
- Figure 4 (V07 `obs_vs_exp` and V08 `cog_abundance`, A1 and A2 freeze runs `out_a1/plots`,
  `out_a2/plots`): detection against detectability, by region and functional category.
- Figure 5 (V11 `operon_pairs`, A2 freeze run `out_a2/plots/V11_operon_pairs`): operon partners on
  the pilot; descriptive, and not evidence for the simulation calibration of operon tests.
- Figure 6 (V13 `novel_orfs`, A2 and A3 freeze runs `out_a2/plots`, `out_a3/plots`): novel ORFs by
  class and call status.
- Supplementary Figure S1 (V16 `compare`, `cmp_a1_a2/` and `cmp_a1_a2_a3/` of the freeze runs):
  detection and abundance agreement between arms.

<!-- src: docs/preprint/outline.md §3 (figure and table plan); directory listing of local/runs/pilot/freeze/fa49b8f/{out_a1,out_a2,out_a3}/plots and cmp_a1_a2*/ (V01, V07, V08, V11, V13, V16 files present);
out_a3/run_manifest.json › git_commit fa49b8f1 -->
