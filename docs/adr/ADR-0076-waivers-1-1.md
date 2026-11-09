# ADR-0076: 1.1 carries forward the 1.0.0 waivers of R2, R5 (TE4 only) and R7
- Status: accepted-by-khoa (2026-11-09, in session: "New 1.1 waiver ADRs (Recommended)"; "accept ADR-0076 and merge it")
- Date: 2026-11-09
- Spec sections: §13.2, §11.6, §7.2.5; `scripts/release/readiness.py` rows R2, R5, R7; ADR-0070, ADR-0071, ADR-0072
- Needs Khoa: yes (release criteria that fail or cannot be evaluated)

Waives readiness R2
Waives readiness R5
Waives readiness R7

## Context
ADR-0070 (R5, TE4's rplL/rplJ pair), ADR-0071 (R2, A2's mapping FDR) and ADR-0072 (R7, the
§13.2 runtime range) each waive their row for 1.0.0 only. `apply_waivers` in `readiness.py`
does not read the version, so it still applies them to 1.1.0rc1. Readiness on `e6f4aa3c`
(2026-11-09) says "overall: pass" only because of these 1.0.0 records.

Nothing these rows measure has changed since 1.0.0:

- **R2:** A2's mapping FDR is still not measured. The pilot data and the `same_db` exemption
  (§7.2.5 step 5) are unchanged, so 0 accepted A2 peptides are decoy-estimated (ADR-0071).
- **R5, TE4:** the rplL/rplJ ratio comes from the same pilot metrics
  (`validation/pilot/metrics/metrics.json`), so it still fails (ADR-0070).
- **R7:** `benchmarks/baseline.json` is unchanged: a K-12 standard run on wheat takes 3:41
  (13,232 s), against the §13.2 range of 0:16-0:55 (ADR-0072).

## Decision
Release 1.1.0rc1 and 1.1.0 with R2, R5 (TE4 only) and R7 waived, disclosed in the release notes
as in 1.0.0.

**Scope.** This ADR covers 1.1.0rc1 and 1.1.0 only. For R5 it covers only TE4's rplL/rplJ pair.
`apply_waivers` works on whole rows, so before step 8's PR the lead reads R5's value: every other
TE must be `pass` or `deferred`, and TE8 needs a weekly run on the readiness commit. R3 is not
waived. Its full tier passes on `e6f4aa3c` (array 1948919); its fast-tier nightly must pass
without a waiver.

## Consequences
- The release notes of 1.1.0rc1 and 1.1.0 say that A2's mapping FDR is not measured, that TE4's
  rplL/rplJ pair fails, and that the runtime is the measured one.
- After 1.1.0, `apply_waivers` should honour a waiver's version scope, so that no record covers a
  release it does not name (backlog). It is not changed before the release: a change to
  `readiness.py` makes no evidence reusable.
