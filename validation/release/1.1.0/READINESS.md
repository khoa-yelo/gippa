# GIPPA 1.1.0rc1 release check

- Commit: `202838346def772d64822d5998eaa88e49619707`
- Computed: 2026-11-11T07:45:21+00:00
- Overall: **pass**

| id | criterion | status | value | threshold | evidence | source |
|---|---|---|---|---|---|---|
| R1 | ≥ 99.5 % of A2 peptides map exactly | pass | 100.00 % | ≥ 99.5 % and te1.passed | validation/pilot/metrics/metrics.json | repo |
| R2 | the decoy-estimated mapping FDR is ≤ 1 % at defaults | waived | no mapping_fdr in validation/pilot/metrics/metrics.json and no --freeze-run | ≤ 1 % | docs/adr/ADR-0071-r2-mapping-fdr-waiver-1-0-0.md (accepted-by-khoa; waived from pending)<br>waiver ADR-0071 | farmshare |
| R3 | Every statistical method passes its calibration suite (§8.4) | pass | calibration_full: 15 suites, 0 failed suites; 555 tests: 0 failed, 0 errors, 459 passed, 90 skipped, 6 xfailed, 0 xpassed; array 1966680; calibration_nightly: 28937 passed, 148 skipped, 17 xfailed, 3 xpassed, 19 warnings in 29214.17s (8:06:54, longest of 4 shards); coverage 97% (min 90%); reused from 4a6fa99 (reusable: diff ⊆ NON_STATS_PATHS) | exit 0 and 0 failed for R or an ancestor C with reusable(C, R) | nightly.tsv 2026-11-10 1966680_calibration 4a6fa998300c0fa0ea6fb01bc9b46f181e24ae35 exit 0<br>nightly.tsv 2026-11-10 1969445_2 20283834 exit 0 | farmshare |
| R4 | every attack null (§8.5) stays at level | pass | calibration_attacks: 14 suites, 0 failed suites; 281 attack tests: 0 failed, 0 errors, 273 passed, 3 skipped, 5 xfailed, 0 xpassed; array 1966680; reused from 4a6fa99 (reusable: diff ⊆ NON_STATS_PATHS) | exit 0 and 0 failed for R or an ancestor C with reusable(C, R) | nightly.tsv 2026-11-10 1966680_calibration 4a6fa998300c0fa0ea6fb01bc9b46f181e24ae35 exit 0 | farmshare |
| R5 | Every true-example assertion that the available data allow (§11.6) passes | waived | TE1 pass; TE2 pass; TE3 pass; TE4 fail (fail); TE5 pass; TE7 pass; TE8 pass | each TE pass; no data → deferred; reported → fail (Q-K5) | validation/pilot/metrics/metrics.json<br>weekly.tsv 2026-11-10 1969418_0 20283834 exit 0<br>docs/adr/ADR-0070-te4-rplL-rplJ-waiver-1-0-0.md (accepted-by-khoa; waived from fail)<br>waiver ADR-0070 | both |
| R6 | The same seed gives byte-identical tables | pass | TE7 pass; integration: integration: 967 passed, 84 skipped, 1 xfailed, 1 xpassed in 2099.81s (0:34:59); nf-test test,stub: SUCCESS: Executed 144 tests in 1228.371s | te7.gippa.status pass and the integration run | validation/pilot/metrics/metrics.json<br>nightly.tsv 2026-11-10 1967956_0 4a6fa998 exit 0 | both |
| R7 | A standard run on 16 threads and 64 GB takes about 35 min (range 16-55 …) | waived | 3:41 > 0:55 (13232 s); node class wheat (16 threads) | 0:16-0:55 (inclusive) | benchmarks/baseline.json > standard_run.wheat.k12_standard<br>docs/adr/ADR-0072-benchmark-baseline-1-0-0.md (accepted-by-khoa; waived from fail)<br>waiver ADR-0072 | repo |
| R8 | Core modules alone take about 25 min | pass | 3:01 (10874 s); node class wheat | reported only | benchmarks/baseline.json > standard_run.wheat.k12_core | repo |
| R9 | The pilot passes (§15.6) | pass | 5 bullets; not passed: 0 | Gate M10 with a Khoa review line; no fail/pending in REPORT.md | PROGRESS.md<br>validation/pilot/REPORT.md | repo |
| R10 | co-abundance … moves to the full validation | deferred | spec: co-abundance moves to the full validation | none | none | repo |
| H1 | milestone gates M1-M10 recorded in PROGRESS.md | pass | 10 gate blocks | ### Gate M1 … M10 | PROGRESS.md | repo |
| H2 | version agreement; release tag not yet created | pass | 1.1.0rc1 everywhere; no tag v1.1.0rc1 | = 1.1.0rc1 | src/gippa/__init__.py<br>CITATION.cff<br>git tag -l | repo |
| H3 | changelog built for the version | pass | built | section present, no fragments | CHANGELOG.md<br>changelog.d/ | repo |
| H4 | containers pinned by digest | pass | 17 digests; 3 conda-only/in-process | <image>@sha256:<64 hex>, no docker:// prefix | workflow/conf/containers.config<br>src/gippa/annotate/tools.py<br>workflow/conf/images.json | repo |
| H5 | licences listed; GIPPA is MIT | pass | all tools listed; MIT | registry ⊆ NOTICE.md table; MIT | NOTICE.md<br>LICENSE | repo |
| H6 | CI check runs on the commit | pass | 11 check runs success; skipped: summary, deploy, guards | every check run success or skipped | gh api …/commits/202838346def772d64822d5998eaa88e49619707/check-runs | farmshare |

## Notes

- Peptide guard (M10-05): runs of 7+ amino-acid letters are shown as [redacted], except names on the fixed allow-list.
