"""M12-04 PR 4, V03: no text-only index image; the index is a report table (P-12; AC 4.6;
docs/specs/M12-04-outputs-revision.md § PR 4).

``registered()["V03"]`` returns no ``Instance`` with ``suffix is None`` (no ``V03_locus.{svg,png,
pdf}``, on a run with loci or without), ``plots/data/V03_locus.tsv`` is still written by the
plots step, and the report shows it as a table "Locus views, in display order" with the columns
rank, locus, family, estimate, q.

Interface assumed (PR 4 test author): the plots step (``gippa.steps.plots.render_plots``) writes
``data/V03_locus.tsv`` however the plot hands it over (V03 may return a non-``Instance`` item or
the step may ask the plot module for it); the TSV keeps a ``rank`` column; the report table's
header cells read exactly ``rank, locus, family, estimate, q`` (case-insensitive) and the
caption "Locus views, in display order" is the table's ``<caption>`` or text of an enclosing
element; one body row per TSV row, ranks in TSV order.
"""

from __future__ import annotations

import csv
from pathlib import Path

import pytest
from helpers import snapshot_fixture as sf
from helpers.m12_04_pr5 import sig_digits
from helpers.report_fixtures import Node, parse, table_rows

from gippa.viz import base

CAPTION = "Locus views, in display order"
HEADER = ["rank", "locus", "family", "estimate", "q"]
IMAGES = ("V03_locus.svg", "V03_locus.png", "V03_locus.pdf")


def _need(fixture: str) -> None:
    why = sf.missing(fixture)
    if why is not None:
        pytest.fail(why)


def _result(fixture: str) -> object:
    from gippa.io.run_view import RunView
    from gippa.steps.plots import load_plot_modules

    _need(fixture)
    load_plot_modules()
    view = RunView.from_outdir(sf.fixture_out(fixture))
    cfg = sf.fixture_cfg(fixture)
    return sf.styled(lambda: base.registered()["V03"].fn(view, cfg))


def _render(fixture: str, plots: Path) -> None:
    from gippa.io.run_view import RunView
    from gippa.steps.plots import load_plot_modules, render_plots

    _need(fixture)
    load_plot_modules()
    view = RunView.from_outdir(sf.fixture_out(fixture))
    render_plots(view, sf.fixture_cfg(fixture), plots, ["V03"])


def _tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8") as fh:
        body = [ln for ln in fh if not ln.startswith("#")]
    return list(csv.DictReader(body, delimiter="\t"))


@pytest.mark.parametrize("fixture", sf.FIXTURES)
def test_v03_returns_no_index_instance(fixture: str) -> None:
    result = _result(fixture)
    if isinstance(result, base.Skip):
        return
    index = [i for i in result if isinstance(i, base.Instance) and i.suffix is None]  # type: ignore[attr-defined]
    assert index == [], [i.title for i in index]


@pytest.mark.parametrize("fixture", sf.FIXTURES)
def test_no_v03_locus_image_is_written(fixture: str, tmp_path: Path) -> None:
    plots = tmp_path / "plots"
    _render(fixture, plots)
    assert [n for n in IMAGES if (plots / n).exists()] == []


def test_v03_locus_tsv_is_still_written(tmp_path: Path) -> None:
    plots = tmp_path / "plots"
    _render(sf.RICH, plots)
    assert [n for n in IMAGES if (plots / n).exists()] == []
    tsv = plots / "data" / "V03_locus.tsv"
    assert tsv.is_file()
    rows = _tsv(tsv)
    assert rows and "rank" in rows[0]
    assert [int(r["rank"]) for r in rows] == list(range(1, len(rows) + 1))
    # the locus views themselves are still drawn
    assert sorted(p.name for p in plots.glob("V03_locus_*.png"))


def _ancestors(node: Node) -> list[Node]:
    out, cur = [], node.parent
    while cur is not None and len(out) < 4:
        out.append(cur)
        cur = cur.parent
    return out


def test_report_shows_the_locus_index_as_a_table(tmp_path: Path) -> None:
    from gippa.io.plot_files import read_plot_files
    from gippa.io.run_view import RunView
    from gippa.report.build import build_report
    from gippa.steps.plots import skipped_plots

    plots = tmp_path / "plots"
    _render(sf.RICH, plots)
    view = RunView.from_outdir(sf.fixture_out(sf.RICH))
    page = parse(
        build_report(view, read_plot_files(plots), skipped_plots(view), sf.fixture_cfg(sf.RICH))
    )
    tables = [
        t for t in page.root.by_tag("table") if [c.lower() for c in table_rows(t)[0]] == HEADER
    ]
    assert len(tables) == 1, [table_rows(t)[0] for t in page.root.by_tag("table")]
    table = tables[0]
    near = [" ".join(n.text().split()) for n in [table, *_ancestors(table)]]
    assert any(CAPTION in t for t in near), near[:2]
    _, body = table_rows(table)
    rows = _tsv(plots / "data" / "V03_locus.tsv")
    assert [r[0] for r in body] == [str(int(r["rank"])) for r in rows]
    for r in body:
        q = r[4]
        assert q in ("NA", "< 0.001") or sig_digits(q) <= 2, q
