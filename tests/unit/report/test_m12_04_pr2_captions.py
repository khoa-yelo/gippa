"""M12-04 PR 2 behaviour 1, AC 2.2, as restated by M12-05c AC 8: every figure caption starts
with "Question " and the plot's §10.2 Question, then the rest of its structured guide
(docs/specs/M12-04-outputs-revision.md § PR 2; docs/specs/M12-05c-plot-guides.md).

Interface assumed (PR 2 test author; updated by the M12-05c PR 1 test author):
- ``gippa.report.content.GUIDES: Mapping[str, Guide]`` replaces ``CAPTIONS`` (M12-05c decision
  7); its keys are the registered plot IDs (lead decision 2);
- the ``<figcaption>`` holds the guide (Question, Panels, How to read, [Caveats,] Data); the file
  links sit in its Data field (M12-05c § 3). The M12-04 one-sentence captions (V02, V17) are
  superseded by the M12-05c § 4 guides, asserted in tests/unit/report/test_m12_05c_*.py.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest
from helpers.m12_04_pr2 import (
    QUESTIONS,
    dress,
    flag_tokens,
    internal_refs,
    sim_report,
    words,
)
from helpers.m12_05c import field, guide_of
from helpers.report_fixtures import Node, Page
from unit.report.test_m12_04_pr1_fixtures import build

#: the index instance of a multi-instance plot (review r1 L2): V13's ``locus`` instance (its
#: suffix-``None`` instance is panel (a)); V03's index is an image-free table since M12-04 PR 4
INDEX_STEMS = {"V13": "V13_novel_orfs_locus"}


def _guides() -> dict[str, Any]:
    from gippa.report.content import GUIDES

    return dict(GUIDES)


def test_guides_keys_are_the_registered_plot_ids() -> None:
    from gippa.steps.plots import load_plot_modules
    from gippa.viz.base import registered

    load_plot_modules()
    assert set(_guides()) == set(registered())
    assert len(_guides()) == 16


def test_question_literal_matches_the_registry() -> None:
    """The pinned §10.2 Questions cover exactly the registered plots (guards the literal)."""
    from gippa.steps.plots import load_plot_modules
    from gippa.viz.base import registered

    load_plot_modules()
    assert set(QUESTIONS) == set(registered())


@pytest.mark.parametrize("pid", sorted(QUESTIONS))
def test_guide_question_then_how_to_read(pid: str) -> None:
    g = _guides()[pid]
    assert g.question == QUESTIONS[pid], g.question
    assert g.read.strip() and g.read.strip().endswith("."), f"{pid}: no how-to-read sentence"
    for text in (g.question, g.read, *g.panels, *([g.caveats] if g.caveats else [])):
        assert not internal_refs(text) and not flag_tokens(text), text


def _plot_of(node: Node) -> str | None:
    cur: Node | None = node
    while cur is not None:
        if "data-plot" in cur.attrs:
            return cur.attrs["data-plot"]
        cur = cur.parent
    return None


def _check_figures(page: Page) -> None:
    guides = _guides()
    figures = page.root.by_tag("figure")
    assert figures, "no <figure> in the report"
    for fig in figures:
        pid = _plot_of(fig)
        assert pid is not None, "a figure outside a plot slot"
        found = fig.by_tag("figcaption")
        assert len(found) == 1, f"{pid}: {len(found)} figcaptions"
        text = words(found[0].text())
        assert text.startswith(f"Question {QUESTIONS[pid]} "), f"{pid}: caption {text!r}"
        links = [a.attrs.get("href", "") for a in found[0].by_tag("a")]
        assert any(h.startswith("../plots/") for h in links), f"{pid}: links lost: {links}"
        read = words(field(guide_of(fig), "How to read").text())
        if found[0].attrs.get("data-stem") == INDEX_STEMS.get(pid):
            # review r1 L2: the index instance has its own how-to-read
            assert read != words(guides[pid].read), f"{pid}: index caption {text!r}"
            continue
        # the links follow the question
        assert words(found[0].by_tag("a")[-1].text()) in text[len(QUESTIONS[pid]) :]
        if guide_of(fig).attrs.get("id") == f"guide-{pid}":
            assert read == words(guides[pid].read), f"{pid}: caption {text!r}"


@pytest.mark.parametrize("name", ["P", "R", "D"])
def test_every_figure_caption_on_the_unit_fixtures(name: str, tmp_path: Path) -> None:
    run = dress(build(name, tmp_path), tmp_path)
    _check_figures(run.page())


def test_every_figure_caption_on_sim_plots() -> None:
    page, _ = sim_report("sim_plots")
    _check_figures(page)


def test_v17_caption_in_the_page(tmp_path: Path) -> None:
    """The whole V17 caption: the guide fields in order, the file links last (M12-05c § 3)."""
    run = dress(build("P", tmp_path), tmp_path)
    cap = run.page().slot("V17").by_tag("figcaption")
    assert cap, "no V17 figcaption"
    g = _guides()["V17"]
    labels = " ".join(words(a.text()) for a in cap[0].by_tag("a"))
    assert labels.startswith("plots/data/V17_qc.tsv SVG"), labels
    want = " ".join(["Question", g.question, "Panels", *g.panels, "How to read", g.read, "Data",
                     "plots/data/V17_qc.tsv; figure files: SVG"])  # fmt: skip
    assert words(cap[0].text()) == words(want), words(cap[0].text())
