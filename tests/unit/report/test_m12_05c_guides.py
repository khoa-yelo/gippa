"""M12-05c AC 1, 3, 8 and 11 on ``gippa.report.content.GUIDES`` itself
(docs/specs/M12-05c-plot-guides.md § 2 decision 7, § 4, § 7).

The guide text is static: one ``Guide(question, panels, read, caveats=None)`` per registered plot.
Only the V02 and V07 guides are asserted word for word (AC 11); the others by AC 1-3 and the
§ 4 panel counts. Interface: tests/helpers/m12_05c.py.
"""

from __future__ import annotations

import dataclasses

import pytest
from helpers.m12_04_pr2 import QUESTIONS, flag_tokens, internal_refs
from helpers.m12_05c import LETTERS, NO_CAVEATS, PANEL_COUNTS, V02_GUIDE, V07_GUIDE, guides
from hypothesis import given
from hypothesis import strategies as st


def _registered() -> set[str]:
    from gippa.steps.plots import load_plot_modules
    from gippa.viz.base import registered

    load_plot_modules()
    return set(registered())


def test_guides_keys_are_the_registered_plot_ids() -> None:
    """AC 1 (and § 6: tests keep the keys equal to the registry)."""
    assert set(guides()) == _registered()
    assert len(guides()) == 16


def test_spec_tables_cover_the_registry() -> None:
    """Guards the pinned literals of this test module."""
    assert set(PANEL_COUNTS) == _registered() == set(QUESTIONS)


def test_guide_is_a_frozen_dataclass_with_the_decision_7_fields() -> None:
    from gippa.report.content import Guide

    assert dataclasses.is_dataclass(Guide)
    assert [f.name for f in dataclasses.fields(Guide)] == ["question", "panels", "read", "caveats"]
    g = Guide("Q?", ("one",), "Read.")
    assert g.caveats is None
    with pytest.raises(dataclasses.FrozenInstanceError):
        g.read = "x"  # type: ignore[misc]


@pytest.mark.parametrize("pid", sorted(QUESTIONS))
def test_every_guide_has_its_required_fields(pid: str) -> None:
    """AC 1: non-empty question and read, at least one panel, caveats None or non-empty."""
    g = guides()[pid]
    assert g.question == QUESTIONS[pid]
    assert g.read.strip(), pid
    assert isinstance(g.panels, tuple) and g.panels, pid
    assert all(isinstance(p, str) and p.strip() for p in g.panels), g.panels
    assert g.caveats is None or g.caveats.strip(), pid
    assert (g.caveats is None) == (pid in NO_CAVEATS), (pid, g.caveats)


@pytest.mark.parametrize("pid", sorted(QUESTIONS))
def test_question_helper_returns_the_guide_question(pid: str) -> None:
    from gippa.report.content import question

    assert question(pid) == guides()[pid].question == QUESTIONS[pid]


@pytest.mark.parametrize("pid", sorted(QUESTIONS))
def test_panels_are_lettered_as_the_plot_prints_them(pid: str) -> None:
    """AC 3 and § 4: a one-panel plot has one unlettered item; a lettered plot's k-th item starts
    "(a) ", "(b) ", ... in order."""
    panels = guides()[pid].panels
    assert len(panels) == PANEL_COUNTS[pid], (pid, panels)
    if len(panels) == 1:
        assert not panels[0].startswith("("), panels[0]
        return
    for k, item in enumerate(panels):
        assert item.startswith(f"({LETTERS[k]}) "), (pid, k, item)


@pytest.mark.parametrize("pid", sorted(QUESTIONS))
def test_every_field_passes_the_m12_04_wording_rules(pid: str) -> None:
    """AC 8 / decision 7: no ADR-, §, GIPPA_SPEC, docs/methods and no flag token in any field."""
    g = guides()[pid]
    texts = [g.question, *g.panels, g.read, *([g.caveats] if g.caveats else [])]
    for text in texts:
        assert not internal_refs(text), (pid, internal_refs(text), text)
        assert not flag_tokens(text), (pid, flag_tokens(text), text)


@pytest.mark.parametrize(("pid", "want"), [("V02", V02_GUIDE), ("V07", V07_GUIDE)])
def test_worked_examples_verbatim(pid: str, want: tuple[object, ...]) -> None:
    """AC 11: the § 4 V02 and V07 guides, word for word."""
    g = guides()[pid]
    assert (g.question, g.panels, g.read, g.caveats) == want


def test_v17_has_five_lettered_panels_and_no_caveats() -> None:
    """Example 2."""
    g = guides()["V17"]
    assert g.caveats is None
    assert [p[:4] for p in g.panels] == ["(a) ", "(b) ", "(c) ", "(d) ", "(e) "]


@given(st.sampled_from(sorted(QUESTIONS)))
def test_guides_are_static(pid: str) -> None:
    """§ 6: the guide text is fixed per plot: its entries are immutable."""
    from gippa.report.content import GUIDES

    g = GUIDES[pid]
    with pytest.raises(dataclasses.FrozenInstanceError):
        g.question = "changed?"  # type: ignore[misc]
