"""Calibration of the M12-04 PR 10 ``putative_start`` cap on the gippasim scenario
``alt_start_initiation`` (statistician sign-off ``local/handoff/pr10_statistician.md`` C6,
"Simulation scenario"; docs/specs/M12-04-outputs-revision.md "PR 10" Statistical evidence).

Fast tier: ``replicates_fast.variants`` (R = 20) at 100 kb, nightly. Full tier:
``replicates_full.variants`` (R = 100) at 1 Mb (``calibration_full``). Three arms: genome
relation ``related`` and ``same`` (the sign-off's two runs), and ``table4`` (``related``,
translation table 4 with ``TTA`` in ``annotation.start_codons``: the arm that kills attack
mutant 3). No tolerance, bound or R is new (R3): the bounds are the sign-off's.

Checks (``helpers.putative_start_records.record``), in every replicate unless noted:

- A sites (A1-A3) that form ALT rows: ``putative_start`` share = 1.0 and 0 ``high`` rows;
- B sites (B1, B2, B3, B_origin): 0 ``putative_start`` rows;
- controls: no ``putative_start`` on the out-of-gene K→M rows, and the annotated non-ATG starts
  read by N-terminal M peptides give no ALT row (M11-20 credit);
- ``id_mixups = (1, 1)``: only C1's planted same-db site;
- the independent oracle flags exactly the chain's rows;
- pooled over the replicates: the rule without its K/R guard flags > 0 B2 rows (the guard has
  teeth), and (full tier) the Wilson 95 % upper bound of the B false-cap share is <= 0.02.

Vacuity (no A, B or B2 row in a replicate) raises ``RuntimeError``, never an assertion.

Scenario interface (gippasim tooling; the implementer writes it): ``gippasim.alt_start``
``alt_start_initiation(seed, *, size="100k" | "1M", genome_relation="related" | "same",
translation_table=11, start_codons=("ATG", "GTG", "TTG")) -> AltStartScenario`` with ``genome``
(``{contig: sequence}``), ``gm`` (``GeneModel``), ``cfg`` (``GippaConfig`` with that relation,
table and start codons, ``sample_name`` set), ``peptide_hits`` (C11 rows, 0-based half-open,
unwrapped), ``peptides`` (``peptide``, ``search_db``) and ``truth`` (one row per planted site;
columns ``putative_start_oracle.TRUTH_COLUMNS``); see the smoke test for the per-replicate
counts.
"""

from __future__ import annotations

import json
import os
from typing import Any

import _criteria as C
import pytest
from helpers import putative_start_oracle as O
from helpers import putative_start_records as P

pytestmark = [pytest.mark.calibration, pytest.mark.calsuite("variants")]
NIGHTLY = pytest.mark.calibration_nightly
SUITE = "variants"
ARMS = tuple(P.ARMS)
FALSE_CAP_UPPER = 0.02  # sign-off: pooled B false-cap share, Wilson 95 % upper bound (full tier)


def diagnostics(record: dict[str, Any]) -> None:
    d = C.log_dir()
    d.mkdir(parents=True, exist_ok=True)
    with (d / "variants_putative_start.jsonl").open("a") as fh:
        fh.write(json.dumps(record, default=str) + "\n")


# ---------------------------------------------------------------- per-PR smoke (one replicate)


def test_smoke_one_replicate_scenario_and_record() -> None:
    """The scenario's planted truth (counts per replicate) and one complete, passing record."""
    scn = P.scenario(0, "related", False)
    t = scn.truth
    assert set(O.TRUTH_COLUMNS) <= set(t.columns), t.columns
    counts = t["kind"].value_counts().to_dict()
    assert {k: counts.get(k, 0) for k in ("A1", "A2", "A3", "B1", "B2")} == {
        "A1": 20, "A2": 20, "A3": 10, "B1": 30, "B2": 30}, counts  # fmt: skip
    assert counts.get("B3", 0) >= 5 and counts.get("B_origin", 0) >= 1, counts
    assert counts.get("control_start", 0) >= 5 and counts.get("control_kM", 0) >= 5, counts
    assert counts.get("same_db", 0) == 1, counts
    a = t[t["kind"].isin(O.A_KINDS)]
    assert a["all_first"].all() and not a["prev_kr"][a["kind"] != "A2"].any()
    assert int(a["prev_kr"][a["kind"] == "A2"].sum()) == 5  # 25 % of A2
    assert abs(int((a["strand"] == "-").sum()) - len(a) // 2) <= 1  # half on the minus strand
    assert int((a["crosses_origin"] & (a["kind"] == "A1")).sum()) == 1
    assert not a[a["kind"] == "A2"]["in_gene"].any() and a[a["kind"] == "A3"]["in_gene"].all()
    b = t[t["kind"].isin(("B1", "B2"))]
    assert b["in_gene"].all() and set(b["ref_codon"]) <= {"GTG", "TTG", "CTG"}
    assert not t[t["kind"] == "B1"]["prev_kr"].any() and t[t["kind"] == "B2"]["prev_kr"].all()
    assert int(t[t["kind"] == "B2"]["all_first"].sum()) >= 3  # the guard must have something
    assert not t[t["kind"].isin(("B3", "B_origin"))]["all_first"].any()
    assert (t[t["kind"] == "control_kM"]["ref_codon"] == "AAG").all()
    rec = P.record(0, "related", False)
    diagnostics({"check": "smoke", **rec})
    assert rec["oracle_agrees"], rec
    assert P.failed(rec) == [], rec
    assert rec["n_a_rows"] == rec["n_a_planted"], rec  # every planted A site forms a row
    assert rec["b2_unguarded"] is not None and rec["b2_unguarded"] >= 1, rec


def test_smoke_table4_arm_plants_tta_sites() -> None:
    scn = P.scenario(0, "table4", False)
    a = scn.truth[scn.truth["kind"].isin(O.A_KINDS)]
    assert int((a["ref_codon"] == "TTA").sum()) >= 5, a["ref_codon"].value_counts()
    assert scn.cfg.annotation.translation_table == 4
    assert "TTA" in scn.cfg.annotation.start_codons


# ---------------------------------------------------------------- fast and full tiers


def _run(arm: str, full: bool) -> None:
    R = C.replicates(SUITE, full)
    got = [P.record(i, arm, full) for i in C.replicate_range(R)]
    n_b = sum(r["n_b_rows"] for r in got)
    k_b = sum(r["b_flagged"] for r in got)
    b2_ung = sum(int(r["b2_unguarded"] or 0) for r in got)
    upper = O.wilson_upper(k_b, n_b)
    bad = {r["replicate"]: P.failed(r) for r in got if P.failed(r)}
    disagree = [r["replicate"] for r in got if not r["oracle_agrees"]]
    diagnostics({"check": "putative_start", "arm": arm, "full": full, "R_run": len(got),
                 "n_b_rows": n_b, "b_flagged": k_b, "b_false_cap_wilson_upper": upper,
                 "b2_unguarded": b2_ung, "failed": bad, "oracle_disagrees": disagree})  # fmt: skip
    for r in got:
        if not (r["n_a_rows"] and r["n_b_rows"] and r["n_b2_rows"]):
            raise RuntimeError(f"alt_start_initiation vacuous in replicate {r['replicate']}: {r}")
    assert not bad, f"putative_start thresholds failed ({arm}): {bad}"
    assert not disagree, f"oracle and chain disagree ({arm}) in replicates {disagree}"
    if got and not os.environ.get(C.SHARD_ENV):
        assert b2_ung > 0, f"{arm}: the unguarded rule flags no B2 row (the guard has no teeth)"
    if full:
        assert upper <= FALSE_CAP_UPPER, (arm, k_b, n_b, upper)


@NIGHTLY
@pytest.mark.parametrize("arm", ARMS)
def test_putative_start_thresholds(arm: str) -> None:
    _run(arm, False)


@pytest.mark.calibration_full
@pytest.mark.parametrize("arm", ARMS)
def test_putative_start_thresholds_full(arm: str) -> None:
    _run(arm, True)
