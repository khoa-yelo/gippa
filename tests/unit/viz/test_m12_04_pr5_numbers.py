# ruff: noqa: RUF001, RUF002, RUF003  (figure text is Unicode: −, ≤)
"""M12-04 PR 5, P-18: numbers in figure text go through ``fmt_q`` / ``fmt_num``
(docs/specs/M12-04-outputs-revision.md § 1 rule 1, § PR 5, AC 5.2).

AC 5.2: no figure text of any plot on fixtures P and K contains a q-value written with more than
2 significant digits, the ASCII ``-`` before a digit, or ``e-0``. Every visible text artist of
every drawn instance is walked after a draw (titles, suptitles, axis and tick labels, legends,
annotations, colour-bar labels), plus each instance's finding title.

Fixtures:

- ``sim_plots``: the committed rich gippasim snapshot run (peptides, iBAQ, every table); P-like.
- ``K``: ``helpers.m12_04_pr5.k_like``, sim_plots carrying the numbers K printed badly: q
  ``4.9e-07`` / ``1.234e-05``, q ``0.432``, a negative Moran's I ``−0.147806``. There is no
  large K-12-like fixture; this one exercises the same writers.
- ``P``, ``R``, ``D``: the PR 1 fixtures (``tests/unit/report/test_m12_04_pr1_fixtures.py``).

Spec examples pinned: V13 "best q 4.9e-07" -> "best q < 0.001"; V03 "q = 0.432" -> "q = 0.43".

A q is a number written right after ``q``, ``q =`` or ``q:``; its significant digits are counted
on the mantissa (``"0.050"`` -> 2). An ASCII minus is ``-`` before a digit that does not join two
word characters (``chr1_1-2`` style identifiers are not signs).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest
from helpers import m12_04_pr5 as h

from gippa.viz import style


@pytest.fixture(autouse=True)
def _style() -> None:
    style.apply()


def _violations(results: dict[str, Any]) -> list[tuple[str, object, str, list[str]]]:
    found = []
    for pid, inst in h.instances(results):
        for text in h.figure_texts(inst):
            bad = h.style_violations(text)
            if bad:
                found.append((pid, inst.suffix, text, bad))
    return found


# --- the checker itself (pure; guards the regexes) --------------------------------------------


@pytest.mark.parametrize(
    ("text", "ok"),
    [
        ("q = 0.43", True),
        ("(q = 0.050, TU-block null)", True),
        ("best q < 0.001", True),
        ("q ≤ 0.05", True),
        ("NES 1.2, q 0.63", True),
        ("Lag-1 Moran's I −0.15", True),
        ("TU-block null", True),
        ("NZ_CP058217_1_1965000-2165000", True),
        ("q = 0.432", False),
        ("q 0.0197", False),
        ("best q 4.9e-07", False),
        ("Moran's I -0.15", False),
        ("(-2.5)", False),
    ],
)
def test_the_checker(text: str, ok: bool) -> None:
    assert (h.style_violations(text) == []) is ok, h.style_violations(text)


# --- AC 5.2 ---------------------------------------------------------------------------------------


@pytest.fixture(scope="module")
def k_results(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Any]:
    return h.render(h.k_like(tmp_path_factory.mktemp("k")))


def test_no_style_violation_on_sim_plots() -> None:
    found = _violations(h.render(h.SIM_OUT))
    assert found == [], "\n".join(map(repr, found))


def test_no_style_violation_on_fixture_k(k_results: dict[str, Any]) -> None:
    found = _violations(k_results)
    assert found == [], "\n".join(map(repr, found))


@pytest.mark.parametrize("name", ["P", "R", "D"])
def test_no_style_violation_on_the_pr1_fixtures(tmp_path: Path, name: str) -> None:
    from unit.report.test_m12_04_pr1_fixtures import build

    run = build(name, tmp_path)
    found = _violations(h.render(run.out))
    assert found == [], "\n".join(map(repr, found))


def test_v13_tiny_best_q_reads_below_0_001(k_results: dict[str, Any]) -> None:
    titles = [i.title for p, i in h.instances(k_results)
              if p == "V13" and str(i.suffix).startswith("locus_")]  # fmt: skip
    assert titles
    for t in titles:
        assert "best q < 0.001" in t, t


def test_v03_region_q_has_two_significant_digits(k_results: dict[str, Any]) -> None:
    """The BGC o_minus_e q of fixture K is 0.432: V03 writes "q = 0.43" (as V07 does)."""
    bgc = [i for p, i in h.instances(k_results) if p == "V03" and str(i.suffix).startswith("bgc")]
    v07 = h.only(k_results, "V07")
    assert any("q = 0.43" in t and "0.432" not in t for t in h.figure_texts(v07))
    for inst in bgc:
        assert "q = 0.43" in inst.title and "0.432" not in inst.title, inst.title
    prophage = [i for p, i in h.instances(k_results)
                if p == "V03" and str(i.suffix).startswith("prophage")]  # fmt: skip
    assert prophage
    for inst in prophage:
        assert "q < 0.001" in inst.title, inst.title
