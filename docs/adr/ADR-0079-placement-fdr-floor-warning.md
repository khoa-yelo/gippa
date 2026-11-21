# ADR-0079: Warn when a placement stratum's q floor is above α
- Status: proposed
- Date: 2026-11-21
- Spec sections: §7.2.5 (placement FDR, G3), §8.1 rows "Peptide placed at a locus (exact)" and
  "Peptide placed with ≥ 1 difference", R3; mirrors M7-11b rule 8
  (`novel_orfs.class_q_floor_above_alpha`) and the gene-level `q_floor_above_alpha`
- Needs Khoa: yes (adds a run warning to the manifest `warnings` list and to the report; Khoa has
  not chosen between a warning, documentation only, or a statistics change)

## Context
Placement FDR (§7.2.5, docs/methods/placement-fdr.md) builds one list per decoy type × difference
stratum m ∈ {0, 1, 2}. Each list is one §8.1 family. At every target E = e:

  FDR_m(e) = (W_m(e) + c) / T_m(e),  q_i = min over e ≥ E_i of FDR_m(e).

For m ≥ 1, π̂ = 1 and c = 1, so this is (D + 1)/T. W ≥ 0 and T_m(e) ≤ T_m, the stratum's total
targets, so **no q in the stratum can fall below c/T_m** (1/T_m for m ≥ 1). When c/T_m > α
(`stats.alpha.placement`, 0.01), the stratum cannot accept anything, even with zero decoys. A
tolerant stratum with fewer than about 1/α = 100 targets accepts nothing.

The run says nothing about this. The analogous floor in the novel-ORF class FDR warns
(`novel_orfs.class_q_floor_above_alpha`, M7-11b rule 8, Khoa 2026-10-22: "q and α are unchanged").
The gene-level tier warns too (`q_floor_above_alpha`). A user sees zero tolerant placements and
cannot tell "no decoy-supported substitutions" apart from "this stratum could never accept one".

**Campaign NCM3722 N2** (`validation/campaign/ecoli/REPORT_1b.md`, "Why recall is zero"; MG1655
genome, NCM3722 spectra, `--search-db other`; `$GIPPA_ROOT/campaign/ncm3722/out_n2/mapping/mapping_fdr.tsv`):
- reverse m = 1: T = 13 (11 of them truth substitution peptides), D = 7, FDR 0.62, 0 accepted.
  The floor is 1/13 = 0.077 > 0.01. With D = 0 the list would still accept nothing.
- reverse m = 2: T = 1, D = 4, floor 1.0, 0 accepted.
- Truth recall of substitutions: 0/11, all `rejected_placement_fdr`. The only run warning was
  `novel_orfs.decoys_not_makedb`.

## Options

### Option 1: run warning, acceptance unchanged (recommended)
In `placement`, for each list (decoy type × m) with T_m > 0 and q_floor = c/T_m **strictly** > α,
log one WARNING `placement.stratum_q_floor_above_alpha` with data
`{decoy_type, n_diff, n_targets, n_decoys, q_floor, alpha}`. Message: "placement: stratum
<decoy_type> m=<m> has <T> targets; its q floor <q_floor> > alpha <alpha>, so none can be
accepted (a stratum with fewer than about c/alpha = <c/α> targets cannot accept anything at
alpha)".
- For m ≥ 1, c = 1 and q_floor = 1/T_m. For stratum 0 the same rule uses that list's c (de novo:
  max π̂; database: the largest ADR-0019 class π̂). It rarely fires for database stratum 0
  (c ≈ 0.016–0.019 in N2 needs T ≤ 1). It does fire for small de novo inputs (B.6: T < 90 at
  c = 0.9).
- The warning reaches the manifest `warnings` list as `placement:placement.stratum_q_floor_above_alpha`
  through the existing events route. The report renders it wherever it renders the class floor.
  It is not a run-flag token (as M7-11b Q4).
- q, `accepted`, `mapping_fdr.tsv` and α do not change. R3 forbids loosening a threshold to make
  the stratum accept, so this option does not.
- The floor is a **sufficient** condition for zero acceptances. It is not necessary: a stratum at
  or below the floor can still accept nothing because of its decoys. The message says "cannot",
  not "did not".

### Option 2: documentation only
Document the floor in the methods note and the user docs, and add no warning. The cost is that
every run hits the same trap silently. N2 needed a manual read of `mapping_fdr.tsv` to explain a
recall of zero. The class FDR and the gene tier already set the precedent of warning for this.

### Option 3: pool the strata (discussed, not proposed)
Pool m = 1 and m = 2 (or all strata) into one list ranked by E, so that small tolerant strata
borrow targets. **This changes the statistics, and this ADR does not propose it.**
- E already carries the difference count through the Poisson-binomial, so a pooled list is
  well-defined. But π̂ differs (stratum 0 uses π̂ < 1, m ≥ 1 uses π̂ = 1). Exact targets would
  dominate the list, so the first tolerant target would get a q computed mostly from exact
  placements. That moves the error control from "FDR ≤ α among placements with ≥ 1 difference"
  (§8.1 row 2) to a pooled FDR in which the tolerant share can have a much higher FDP. This is the
  same small-family dilution M3-10 rejected for per-class lists in the other direction.
- Pooling only m = 1 with m = 2 keeps π̂ = 1 but helps little in N2 (T = 14 pooled).
- Any such change alters §8.1's families and needs its own ADR, calibration under §8.5 and a
  Khoa decision under R3.

## Recommendation
Option 1, together with the Option 2 documentation (methods paragraph added to
docs/methods/placement-fdr.md "Strata and ranking"). It is reversible: it adds one log event and
changes no number. It matches the two existing floor warnings.

**Where Option 1 would have fired on the campaign runs** (`placement` rows of `mapping_fdr.tsv`,
length class `all`):

| run | list | T | D | FDR | accepted | q floor | warns |
|---|---|---|---|---|---|---|---|
| N2 | reverse m = 1 | 13 | 7 | 0.62 | 0 | 0.077 | yes |
| N2 | reverse m = 2 | 1 | 4 | 5 | 0 | 1.0 | yes |
| N2 | reverse m = 0 | 24,539 | 86 | 5.6e-5 | 24,539 | 7.6e-7 | no |
| N1 | reverse m = 1 | 1 | 7 | 8 | 0 | 1.0 | yes |
| E1, E2, E1-TE7 | reverse m = 1 | 1 | 6 | 7 | 0 | 1.0 | yes |
| E1, E2 | reverse m = 0 | 0 (all `same_db`) | 0 | NA | 0 | none | no |
| E4 (de novo) | shuffle m = 2 | **68** | **0** | 0.0147 | **0** | 0.0147 | yes |
| E4 | shuffle m = 1 | 314 | 1 | 0.0064 | 314 | 0.0032 | no |
| E4 | shuffle m = 0 | 700 | 5 | 0.0063 | 700 | ≤ 0.0014 | no |

E2 does have a tolerant stratum, but it is one target against six decoys. The floor and the decoys
both reject it, so the warning there is information only. E4's m = 2 list is the clearest case:
68 targets, **no** decoy, and still nothing accepted. That is the floor alone (1/68 = 0.0147).

## Spec for tests (only if Option 1 is accepted)
Binding semantics for the test-author. Unit tests on `placement_tables` (or a pure helper the
implementer exposes) and one step test on the event and manifest token. α = 0.01 unless stated.

| # | Input (one list) | Expected |
|---|---|---|
| 1 | reverse, m = 1, T = 13, D = 7 (N2) | one WARNING, data `{decoy_type: reverse, n_diff: 1, n_targets: 13, n_decoys: 7, q_floor: 1/13, alpha: 0.01}`; 0 accepted |
| 2 | reverse, m = 1, T = 100, D = 0 | q_floor = 0.01 = α: **no** warning (strict >); all 100 accepted (B.6 boundary) |
| 3 | reverse, m = 1, T = 99, D = 0 | q_floor 0.010101 > α: warning; 0 accepted |
| 4 | reverse, m = 2, T = 0, D = 5 | no warning (no targets) |
| 5 | shuffle, m = 2, T = 68, D = 0 (E4) | warning, q_floor 1/68; 0 accepted |
| 6 | de novo stratum 0, c = 0.9, T = 89 / 90 | warning (0.01011) / no warning (0.01) |
| 7 | database stratum 0, c = 0.0186, T = 1 / 2 | warning (0.0186) / no warning (0.0093) |
| 8 | α = 0.05 (config), reverse m = 1, T = 13 / T = 25 | warning, `alpha: 0.05` in data (1/13 = 0.077 > 0.05) / no warning (0.04); the rule reads `stats.alpha.placement` |
| 9 | both decoy types and three strata in one run, two lists below the floor | exactly two WARNINGs, one per list, ordered by decoy type then m |

Invariants (property tests, Hypothesis):
- **Acceptance unchanged.** `per_peptide` and `mapping_fdr` are identical to the pre-ADR output for
  the same costs table (golden comparison on the existing fixtures).
- **Soundness of the floor.** For random lists: every target q ≥ c/T_m; and when the warning
  fires, the list's `n_accepted` = 0.
- **Manifest.** A run whose `placement` step logs the event lists
  `placement:placement.stratum_q_floor_above_alpha` in `warnings`. The MVP-04 gate-2 exact-token
  match for `q_floor_above_alpha` must **not** match it (the `_` before `q_floor` is a word
  character, as for the class floor; extend `tests/unit/test_mvp04_gate2_token_match.py`).

Calibration: none needed. No q changes, so the existing §8.5 placement attacks and the M3-21
genome-level calibration stay valid as they are.

## Consequences
- Most database runs with a handful of tolerant hits will carry this warning (E1, E2, N1, N2 all
  do). This is intended: none of those runs could report a tolerant placement at α = 0.01, and the
  user should know that before reading zero substitutions as biology.
- Tests and fixtures that pin the exact manifest `warnings` set of a run with small tolerant
  strata will change (a new entry only). Baselines change only in the warnings list.
- The event name follows the brief (`placement.` prefix). The `placement` step's existing event is
  `mapping.input_prefiltered`; the implementer keeps the name given here unless Khoa prefers
  `mapping.stratum_q_floor_above_alpha` at acceptance.
- If Khoa chooses Option 2, the methods paragraph stays and the warning is not implemented. Option 3
  needs a separate ADR with calibration evidence.
- GIPPA at defaults still cannot report substitutions between closely related strains (N2 recall
  0/11). This ADR makes that visible; it does not fix it.
