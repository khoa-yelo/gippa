# ADR-0080: 1.1.1 carries forward the R2, R5 (TE4 only) and R7 waivers
- Status: accepted-by-khoa (2026-11-26, in session: "New ADR-0080 (Recommended)")
- Date: 2026-11-26
- Spec sections: §13.2, §11.6, §7.2.5; `scripts/release/readiness.py` rows R2, R5, R7; ADR-0070,
  ADR-0071, ADR-0072, ADR-0073, ADR-0076
- Needs Khoa: yes (release criteria that fail or cannot be evaluated)

Waives readiness R2
Waives readiness R5
Waives readiness R7

## Context
ADR-0070 (R5, TE4's rplL/rplJ pair), ADR-0071 (R2, A2's mapping FDR) and ADR-0072 (R7, the
§13.2 runtime range) each waive their row for 1.0.0 only, and ADR-0076 carried them forward for
1.1.0rc1 and 1.1.0 only. `apply_waivers` in `readiness.py` still does not read a waiver's version,
so readiness on the 1.1.1 release commit R = `d2bd17e4` (#609) says "overall: pass" only because of
these older records.

Nothing these rows measure has changed since `v1.1.0`: `git diff v1.1.0 d2bd17e4` is empty for
`validation/pilot/metrics/` and `benchmarks/baseline.json`.

- **R2:** A2's mapping FDR is still not measured (ADR-0071; no `--freeze-run`).
- **R5, TE4:** the rplL/rplJ ratio comes from the same pilot metrics, so it still fails
  (ADR-0070). Every other TE R5 reads passes: TE1, TE2, TE3, TE5, TE7, and TE8 on the weekly
  run on R (`2020593_0`, exit 0), as in 1.1.0.
- **R7:** a K-12 standard run on wheat takes 3:41 (13,232 s) against the §13.2 range of
  0:16-0:55 (ADR-0072).

R3 is not waived. On R the full tier passes (`2022338_calibration`: 558 tests, 0 failed) and so
does the fast tier: `2023930_2` (4 shards, 10 h limit) gave 29,049 passed, 0 failed, coverage 97 %.
The scheduled nightly's fast tier (`2021519_2`) had hit the 6 h limit at 99 % of shard 1 with
0 failures. ADR-0073, the 1.0.0 R3 waiver, which `apply_waivers` would also have applied, is not
relied on.

## Decision
Release 1.1.1 with R2, R5 (TE4 only) and R7 waived, disclosed in the release notes as in 1.1.0.

**Scope.** This ADR covers 1.1.1 only. For R5 it covers only TE4's rplL/rplJ pair.

## Consequences
- The 1.1.1 release notes say that the 1.1.0 limitations for these rows still apply: A2's mapping
  FDR is not measured, TE4's rplL/rplJ pair fails, and the runtime is the measured one.
- The backlog items stay open: `apply_waivers` should honour a waiver's version scope, and the
  nightly fast tier needs a time limit or shard count that fits it. Neither is changed before the
  release: a change to `readiness.py` or `scripts/slurm/` is not paperwork.
