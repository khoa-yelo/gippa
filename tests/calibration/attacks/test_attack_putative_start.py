"""Teeth of the M12-04 PR 10 ``putative_start`` checks (statistician sign-off
``local/handoff/pr10_statistician.md`` C6 "Attack"; GIPPA_SPEC.md §8.5 "A failing attack blocks the
merge"). Each mutant must fail at least one threshold of
tests/calibration/test_putative_start_calibration.py on the same scenario
(``alt_start_initiation``), and the expected one:

1. ``any_first``: "any peptide first" in place of "every" -> B sites flagged (B3, B_origin);
2. ``genomic_order``: positions from the lowest coordinate on the minus strand -> A share < 1;
3. ``hardcoded_starts``: ``{GTG, TTG}`` -> A share < 1 in the ``table4`` (TTA) arm;
4. ``no_modulo``: codons matched without ``mod L`` -> the origin A1 site or B_origin fails;
5. ``no_kr_guard``: condition 5 dropped -> B2 sites flagged;
6. ``no_cap``: the flag without the cap -> A rows reach ``high``;
7. ``exempt_in_id_mixups``: ``putative_start`` in ``EXEMPT_FLAGS`` -> ``id_mixups != (1, 1)``.

Mutants 1-5 are applied to the oracle (``helpers.putative_start_oracle``) on the chain's own
alignments; the calibration test pins the oracle to the chain (``oracle_agrees``), so a mutant of
the oracle stands for the same mutant of the implementation. Mutants 6-7 act on the chain's output
(``helpers.putative_start_records.outcome``). The mutants reuse the calibration replicates
(R = ``replicates_fast.variants``, 100 kb), evaluated over the replicates (pooled: a mutant is
caught when its check fails in at least one replicate). Nightly; records go to
``<log_dir>/variants_putative_start_attack.jsonl``.
"""

from __future__ import annotations

import json
from typing import Any

import _criteria as C
import pytest
from helpers import putative_start_records as P

pytestmark = [pytest.mark.calibration, pytest.mark.calsuite("variants")]
NIGHTLY = pytest.mark.calibration_nightly
SUITE = "variants"

#: mutant -> (arms it runs on, the check it must fail)
EXPECTED: dict[str, tuple[tuple[str, ...], frozenset[str]]] = {
    "any_first": (("related", "same"), frozenset({"b_none"})),
    "genomic_order": (("related", "same"), frozenset({"a_share_1"})),
    "hardcoded_starts": (("table4",), frozenset({"a_share_1"})),
    "no_modulo": (("related", "same"), frozenset({"a_share_1", "b_none"})),
    "no_kr_guard": (("related", "same"), frozenset({"b_none"})),
    "no_cap": (("related", "same"), frozenset({"a_no_high"})),
    "exempt_in_id_mixups": (("related", "same"), frozenset({"id_mixups"})),
}


def _log(rec: dict[str, Any]) -> None:
    d = C.log_dir()
    d.mkdir(parents=True, exist_ok=True)
    with (d / "variants_putative_start_attack.jsonl").open("a") as fh:
        fh.write(json.dumps(rec, default=str) + "\n")


def test_smoke_every_mutant_fails_on_one_replicate() -> None:
    """Per-PR: one ``related`` replicate (and one ``table4`` for mutant 3) already catches each
    mutant; the correct chain passes the same replicate."""
    assert P.failed(P.record(0, "related", False)) == []
    for mutant, (arms, want) in EXPECTED.items():
        arm = arms[0]
        got = set(P.failed(P.record(0, arm, False, mutant)))
        assert got & want, (mutant, arm, got)


@NIGHTLY
@pytest.mark.parametrize("mutant", list(EXPECTED))
def test_mutant_fails_its_check(mutant: str) -> None:
    arms, want = EXPECTED[mutant]
    R = C.replicates(SUITE, False)
    caught: dict[str, int] = {}
    for arm in arms:
        for i in C.replicate_range(R):
            ok = P.record(i, arm, False)
            assert P.failed(ok) == [], (arm, i, ok)  # the correct chain passes the replicate
            rec = P.record(i, arm, False, mutant)
            for check in set(P.failed(rec)) & want:
                caught[check] = caught.get(check, 0) + 1
    _log({"mutant": mutant, "arms": list(arms), "R": R, "caught": caught})
    assert caught, f"mutant {mutant} passed every {sorted(want)} check over R = {R} replicates"
