# ADR-0073: 1.0.0 ships without a fast-tier calibration nightly on the release commit
- Status: accepted-by-khoa (2026-11-03, in session: "Option 2")
- Date: 2026-11-03
- Spec sections: §8.4, §12, §13.2; `docs/specs/M11-03-release-readiness.md` row R3;
  `scripts/release/readiness.py` (`r3`, `reusable`)
- Needs Khoa: yes (a release criterion is not evaluated on the release commit)

Waives readiness R3

## Context
R3 reads "every statistical method passes its calibration suite (§8.4)". `readiness.py` asks for the
newest `calibration_full` row and the newest `calibration_nightly` row (the fast tier) for the release
commit R, or for an ancestor C whose diff to R is reusable. On `870042eb` (main after #526):

- `calibration_full` passes, reused from `10997697` (#524), whose diff to R is in reusable paths only:
  15 suites, 0 failed suites, 552 tests, 0 failed, 0 errors; attacks 14 suites, 281 attack tests,
  0 failed (`logs/nightly.tsv` row `1924911_calibration`, 2026-11-02).
- `calibration_nightly` has no row for R or a reusable ancestor, so R3 reads `pending`.

The fast tier did not run on R for reasons outside the code:
- the 01:30 nightly of 2026-11-03 (array 1927650) failed in its first seconds on barley-04 (exit
  0:53), a node class that does not mount the home directory;
- a rerun of task 2 (1928572) failed the same way on barley-03;
- the rerun restricted to wheat nodes (1928574) is queued; Slurm estimates its start at
  2026-11-04 10:18.

## Evidence
- The fast tier passed at `fa49b8f1`, the report-freeze commit: nightly `1914426_2`, exit 0,
  25,988 passed, 0 failed, coverage 97 % (2026-10-31).
- Between `fa49b8f1` and R, nothing changed under `src/`, `workflow/`, `tests/calibration/` or `ci/`.
  The changes are documentation, the changelog, `PROGRESS.md`, `benchmarks/baseline.json`,
  `validation/pilot/` (the report and the entrapment tool), `scripts/dev/`, `scripts/slurm/entrapment.sbatch`,
  `tests/helpers/entrapment_m10_07.py` and nine `tests/unit/` files. `readiness.py` does not reuse a
  row across `tests/unit/` or `validation/` changes, which is why the `fa49b8f1` row does not count.
- The one fast-tier run in between, at `a629525a` (`1924210_2`), exited 1 with 2 failed tests. Both
  were `tests/unit/test_m6_11_regulondb_fetch.py::test_main_is_idempotent_and_appends_the_second_date`,
  whose fixed "later" date became the current date on 2026-11-02 (UTC). #524 fixed the test; no
  calibration test failed in that run (26,085 passed).
- CI on R's ancestors `bbef9718` and on #526 passed every check, macOS included.

## Decision
Release 1.0.0 with R3 waived for the fast tier on the release commit. The full tier stays the gating
evidence, and no calibration tolerance or test changes. The queued run 1928574 stays queued; its
result is recorded in `PROGRESS.md` when it ends, and a failure there is handled as a 1.0.1 issue
stated in the release notes.

## Consequences
- The fast tier was last observed passing at `fa49b8f1`, not at the release commit. The difference is
  test and tooling code outside the statistics, as listed above.
- The nightly launcher and the calibration launcher need a node exclusion for the barley and rye
  nodes that cannot run these jobs (`docs/backlog.md`).
- The waiver applies to 1.0.0 only.
