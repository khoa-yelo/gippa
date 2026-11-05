"""M12-04 PR 8 (F-01): V08 and V14 when no function term was assigned
(docs/specs/M12-04-outputs-revision.md § PR 8 behaviour, AC 8.1, AC 8.4; PR 4 leaves V14 at
coverage 0 to this PR).

Interface pinned here (PR 8 test author):
- V08 (``cog_abundance``) returns ``Skip("no_terms:function")`` when the run's last
  ``function.collections`` event has at least one collection and every collection has coverage
  0 or 0 terms; the reason is decided inside the plot (it is not a declared skip of
  ``base.skip_reason``), so ``render_plots`` logs it as ``plots.skipped`` and writes no file;
- ``gippa.viz.words.SKIP_TEXT`` has the key ``no_terms:function`` (or the template
  ``no_terms:``) and ``plain("no_terms:function")`` is "no function terms in this run";
- V14 is one ``Instance`` with ``empty=True`` whose title follows the P-14 pattern ("Nothing to
  show: ...") and names the collections like the report notice ("COG, GO or KEGG term", or
  "COG or GO term" without a ``kegg`` key); its data TSV is header-only with ``COLUMNS``;
- with a ``function.collections`` event that has terms, no event at all, or ``function`` off,
  V08 and V14 behave as before.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from helpers.m12_04_pr8 import GO, NO_TERMS, NO_TERMS_TEXT, function_run, partial, render, zero
from unit.report.test_m12_04_pr1_fixtures import MODULES_BC, fixture_p

from gippa.viz import base

NOTHING = "Nothing to show:"


def test_skip_text_has_the_no_terms_entry() -> None:
    from gippa.viz.words import SKIP_TEXT, plain

    assert NO_TERMS in SKIP_TEXT or "no_terms:" in SKIP_TEXT, sorted(SKIP_TEXT)
    assert plain(NO_TERMS) == NO_TERMS_TEXT


@pytest.mark.parametrize("light", [True, False])
def test_v08_skips_when_every_coverage_is_zero(tmp_path: Path, light: bool) -> None:
    run = function_run(tmp_path, collections=zero(), light=light)
    result = render(run, ["V08"])["V08"]
    assert isinstance(result, base.Skip), f"V08 drawn: {[i.title for i in result]}"
    assert result.reason == NO_TERMS


def test_v08_skips_without_kegg_when_cog_and_go_are_zero(tmp_path: Path) -> None:
    run = function_run(tmp_path, collections=zero("cog", *GO))
    result = render(run, ["V08"])["V08"]
    assert isinstance(result, base.Skip) and result.reason == NO_TERMS, result


def test_v08_is_not_a_declared_skip(tmp_path: Path) -> None:
    """The no-terms reason comes from the plot; the declared-skip check is unchanged."""
    run = function_run(tmp_path, collections=zero())
    assert base.skip_reason("V08", run.view()) is None


def test_render_plots_logs_the_skip_and_writes_no_v08_file(tmp_path: Path) -> None:
    from gippa.steps.plots import load_plot_modules, render_plots

    load_plot_modules()
    run = function_run(tmp_path, collections=zero())
    plots = tmp_path / "plots"
    events = render_plots(run.view(), run.cfg(), plots, ["V08"])
    assert [(e.event, e.data) for e in events] == [
        ("plots.skipped", {"plot": "V08", "reason": NO_TERMS})
    ]
    assert not list(plots.rglob("V08_*")), sorted(p.name for p in plots.rglob("*"))


@pytest.mark.parametrize(
    ("colls", "names", "absent"),
    [
        (zero(), "COG, GO or KEGG term", None),
        (zero("cog", *GO), "COG or GO term", "KEGG"),
    ],
    ids=["with_kegg", "no_kegg"],
)
def test_v14_is_an_empty_state_naming_the_collections(
    tmp_path: Path, colls: dict[str, tuple[float, int]], names: str, absent: str | None
) -> None:
    from gippa.viz.enrichment import COLUMNS, EMPTY

    run = function_run(tmp_path, collections=colls)
    result = render(run, ["V14"])["V14"]
    assert not isinstance(result, base.Skip), result
    assert len(result) == 1
    inst = result[0]
    assert inst.empty is True
    assert inst.title.startswith(NOTHING), inst.title
    assert names in inst.title, inst.title
    assert EMPTY not in inst.title, inst.title
    if absent is not None:
        assert absent not in inst.title, inst.title
    assert list(inst.data.columns) == list(COLUMNS) and len(inst.data) == 0


def test_v14_example_wording_with_kegg(tmp_path: Path) -> None:
    """The examples table, run A: "Nothing to show: no gene has a COG, GO or KEGG term."."""
    from helpers.m12_04_pr8 import V14_EMPTY

    run = function_run(tmp_path, collections=zero())
    inst = render(run, ["V14"])["V14"][0]
    assert inst.title.startswith(V14_EMPTY), inst.title


# ---------------------------------------------------------------- unchanged


def test_v08_drawn_at_partial_coverage_without_kegg(tmp_path: Path) -> None:
    """AC 8.4: no KEGG, COG coverage > 0: V08 is drawn."""
    run = function_run(tmp_path, collections=partial(kegg=None), n_cog=60)
    result = render(run, ["V08"])["V08"]
    assert not isinstance(result, base.Skip), result
    assert len(result) == 1 and not result[0].empty


def test_v08_and_v14_unchanged_with_terms(tmp_path: Path) -> None:
    run = function_run(tmp_path, collections=partial(), n_cog=60)
    out = render(run)
    assert not isinstance(out["V08"], base.Skip) and not out["V08"][0].empty
    assert not out["V14"][0].title.startswith(NOTHING + " no gene has"), out["V14"][0].title


def test_v08_unchanged_without_the_event(tmp_path: Path) -> None:
    """An old run (no ``function.collections`` event): coverage unknown, V08 drawn."""
    run = function_run(tmp_path, collections=None)
    result = render(run, ["V08"])["V08"]
    assert not isinstance(result, base.Skip), result


def test_v08_unchanged_with_function_off(tmp_path: Path) -> None:
    run = fixture_p(tmp_path, "B", modules=MODULES_BC)
    result = render(run, ["V08"])["V08"]
    assert not isinstance(result, base.Skip), result
