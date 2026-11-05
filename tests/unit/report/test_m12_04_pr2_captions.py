"""M12-04 PR 2 behaviour 1, AC 2.2: every figure caption starts with the plot's §10.2 Question,
then one "how to read" sentence from ``gippa.report.content.CAPTIONS``, then the file links
(docs/specs/M12-04-outputs-revision.md § PR 2).

Interface assumed (PR 2 test author):
- ``gippa.report.content.CAPTIONS: Mapping[str, str]``, plot ID (``V01`` ... ``V18``) -> the
  whole caption sentence pair; its keys are the registered plot IDs (lead decision 2);
- the ``<figcaption>`` text is the caption followed by the link labels (``SVG``, ``PNG``, ...);
  the instance stem no longer opens it (spec example "V17_qc SVG PNG PDF data TSV" -> "Are the
  inputs sane? ... SVG PNG PDF data TSV").
"""

from __future__ import annotations

from pathlib import Path

import pytest
from helpers.m12_04_pr2 import (
    QUESTIONS,
    dress,
    flag_tokens,
    internal_refs,
    sim_report,
    words,
)
from helpers.report_fixtures import Node, Page
from unit.report.test_m12_04_pr1_fixtures import build

#: the spec's two worked examples (behaviour 1 and the examples table), verbatim
V02_CAPTION = (
    "Is abundance clustered along the chromosome? Each point is a detected gene's log2 abundance "
    "at its position; the line is a LOESS trend; shaded bands are abundance islands."
)
V17_CAPTION = (
    "Are the inputs sane? Panels show the preflight summary, contig lengths, peptide lengths, "
    "abundances and missing values per feature."
)


#: the index instance of a multi-instance plot (review r1 L2): V03's suffix-``None`` instance,
#: V13's ``locus`` instance (its suffix-``None`` instance is panel (a))
INDEX_STEMS = {"V03": "V03_locus", "V13": "V13_novel_orfs_locus"}


def _captions() -> dict[str, str]:
    from gippa.report.content import CAPTIONS

    return dict(CAPTIONS)


def test_captions_keys_are_the_registered_plot_ids() -> None:
    from gippa.steps.plots import load_plot_modules
    from gippa.viz.base import registered

    load_plot_modules()
    assert set(_captions()) == set(registered())
    assert len(_captions()) == 16


def test_question_literal_matches_the_registry() -> None:
    """The pinned §10.2 Questions cover exactly the registered plots (guards the literal)."""
    from gippa.steps.plots import load_plot_modules
    from gippa.viz.base import registered

    load_plot_modules()
    assert set(QUESTIONS) == set(registered())


@pytest.mark.parametrize("pid", sorted(QUESTIONS))
def test_caption_starts_with_question_then_one_reading_sentence(pid: str) -> None:
    caption = _captions()[pid]
    assert caption.startswith(QUESTIONS[pid] + " "), caption
    rest = caption[len(QUESTIONS[pid]) :].strip()
    assert rest and rest.endswith("."), f"{pid}: no how-to-read sentence after the question"
    assert not internal_refs(caption) and not flag_tokens(caption), caption


def test_caption_examples_from_the_spec() -> None:
    caps = _captions()
    assert caps["V02"] == V02_CAPTION
    assert caps["V17"] == V17_CAPTION


def _plot_of(node: Node) -> str | None:
    cur: Node | None = node
    while cur is not None:
        if "data-plot" in cur.attrs:
            return cur.attrs["data-plot"]
        cur = cur.parent
    return None


def _check_figures(page: Page) -> None:
    caps = _captions()
    figures = page.root.by_tag("figure")
    assert figures, "no <figure> in the report"
    for fig in figures:
        pid = _plot_of(fig)
        assert pid is not None, "a figure outside a plot slot"
        found = fig.by_tag("figcaption")
        assert len(found) == 1, f"{pid}: {len(found)} figcaptions"
        text = words(found[0].text())
        assert text.startswith(QUESTIONS[pid]), f"{pid}: caption {text!r}"
        links = [a.attrs.get("href", "") for a in found[0].by_tag("a")]
        assert any(h.startswith("../plots/") for h in links), f"{pid}: links lost: {links}"
        if found[0].attrs.get("data-stem") == INDEX_STEMS.get(pid):
            # review r1 L2: the index instance has its own caption
            assert not text.startswith(words(caps[pid])), f"{pid}: index caption {text!r}"
            continue
        assert text.startswith(words(caps[pid])), f"{pid}: caption {text!r}"
        # the links follow the caption
        assert words(found[0].by_tag("a")[0].text()) in text[len(words(caps[pid])) :]


@pytest.mark.parametrize("name", ["P", "R", "D"])
def test_every_figure_caption_on_the_unit_fixtures(name: str, tmp_path: Path) -> None:
    run = dress(build(name, tmp_path), tmp_path)
    _check_figures(run.page())


def test_every_figure_caption_on_sim_plots() -> None:
    page, _ = sim_report("sim_plots")
    _check_figures(page)


def test_v17_caption_in_the_page(tmp_path: Path) -> None:
    run = dress(build("P", tmp_path), tmp_path)
    cap = run.page().slot("V17").by_tag("figcaption")
    assert cap, "no V17 figcaption"
    labels = " ".join(words(a.text()) for a in cap[0].by_tag("a"))
    assert labels.startswith("SVG"), labels
    assert words(cap[0].text()) == f"{V17_CAPTION} {labels}", words(cap[0].text())
