"""M12-04 PR 4 (#565) review round 1: twins that match the static plots, the surviving mutants,
the V06 assumed text with NA values and the V10 panel note's extent
(docs/specs/M12-04-outputs-revision.md § PR 4).

- The Plotly twins in report.html read the same data TSV as the static plot and must show what
  it shows:
  - V06 on D: no (a)/(b) subplots, and the "placed by exact match" line instead;
  - V06: "(CI)" on the (c) y axis only when a CI is drawn;
  - V06: a 1.0 TSV without ``drawn`` still builds the twin (mutant M9);
  - V11: the (c) title names only the drawn classes, legend names are ``class_label``
    ("OP (same operon)"), OP is drawn over NO with opacity 0.6;
  - V10: the text of a partial panel goes through ``plain()``.
- Mutants:
  - M4: a ``flags`` column mixing NA and text renders, the NA cells count as no flag;
  - M5: one placement row with ``targets = 1`` is not exact-only;
  - M2b: a degenerate CI (1 to 1) gets no "(CI)" label;
  - M6b/M6c: the V03 report table matches its TSV in every column;
  - M13: V15 with no features table says "TUs not recorded".
- The 0-row V03 index is a report note (``data-index-none``).
- Finding 4: measured de novo bins whose values are all NA are not "assumed c = 0.99".
- The V10 "Not drawn: …" note stays inside its axes.

Interface assumed (PR 4 review test author):
- a twin is checked on the static instance's data read back as its TSV
  (``helpers.twins.instance_frame``) with the static finding title;
- a twin subplot title is a ``layout.annotations`` text that starts with "(a)", "(b)" or "(c)";
- a twin trace's ``name`` is its legend entry; the later trace of two is drawn on top.
"""

from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Any

import pytest
from helpers import m12_04_pr4 as h
from helpers import snapshot_fixture as sf
from helpers import twins as tw
from helpers.m12_04_pr5 import render
from helpers.report_fixtures import parse, table_rows
from helpers.viz_mapping import EXAMPLE4_TOTALS, denovo_row, full_placement, mapping_out
from helpers.viz_out_tree import make_out
from unit.report.test_m12_04_pr1_fixtures import write_not_fitted

from gippa.viz import base, style
from gippa.viz.words import plain

pytestmark = pytest.mark.filterwarnings("ignore::DeprecationWarning")

PANEL = re.compile(r"^\(([abc])\)")


def _one(out: Path, pid: str) -> base.Instance:
    result = render(out, [pid])[pid]
    assert not isinstance(result, base.Skip), result
    assert len(result) == 1, [i.title for i in result]
    return result[0]


def _twin(pid: str, inst: base.Instance, drop: tuple[str, ...] = ()) -> dict[str, Any]:
    data = inst.data.drop(columns=list(drop))
    return tw.spec(pid, tw.instance_frame(data, inst.coords), inst.title)


def _subplot_titles(s: dict[str, Any]) -> list[str]:
    return [a for a in tw.annotations(s) if PANEL.match(a)]


def _ylabels(inst: base.Instance) -> list[str]:
    return [ax.get_ylabel() for ax in inst.figure.axes if ax.get_ylabel()]


# ---------------------------------------------------------------- V06 twin


def test_v06_twin_on_d_has_no_a_b_subplots_and_says_exact_match(tmp_path: Path) -> None:
    inst = _one(h.d_full(tmp_path).out, "V06")
    s = _twin("V06", inst)
    panels = {PANEL.match(t).group(1) for t in _subplot_titles(s)}  # type: ignore[union-attr]
    assert not panels & {"a", "b"}, _subplot_titles(s)
    assert tw.traces(s, "bar") == []
    assert any("placed by exact match" in a for a in tw.annotations(s)), tw.annotations(s)


def _no_ci_bins(tmp: Path) -> base.Instance:
    rows = [denovo_row(i, 0.9, None, None) for i in (1, 2, 3)]
    return _one(mapping_out(tmp, fdr=[*full_placement(EXAMPLE4_TOTALS), *rows]), "V06")


def test_v06_twin_ci_label_only_with_a_drawn_ci(tmp_path: Path) -> None:
    inst = _no_ci_bins(tmp_path)
    assert not any("(CI)" in t for t in _ylabels(inst))  # the static plot: no CI, no "(CI)"
    s = _twin("V06", inst)
    assert not any("(CI)" in t for t in tw.axis_titles(s)), tw.axis_titles(s)


def test_v06_twin_ci_label_kept_with_a_ci(tmp_path: Path) -> None:
    rows = [denovo_row(i, 0.9, 0.85, 0.95) for i in (1, 2, 3)]
    inst = _one(mapping_out(tmp_path, fdr=[*full_placement(EXAMPLE4_TOTALS), *rows]), "V06")
    assert any("(CI)" in t for t in _ylabels(inst))
    assert any("(CI)" in t for t in tw.axis_titles(_twin("V06", inst)))


def test_v06_twin_builds_from_a_tsv_without_drawn(tmp_path: Path) -> None:
    """Mutant M9: a 1.0 data TSV has no ``drawn`` column; the twin is still built."""
    rows = [denovo_row(i, 0.9, 0.85, 0.95) for i in (1, 2, 3)]
    inst = _one(mapping_out(tmp_path, fdr=[*full_placement(EXAMPLE4_TOTALS), *rows]), "V06")
    s = _twin("V06", inst, drop=("drawn",))
    assert tw.traces(s, "bar")
    assert {PANEL.match(t).group(1) for t in _subplot_titles(s)} == {"a", "b", "c"}  # type: ignore[union-attr]


# ---------------------------------------------------------------- V06 static (M5, M2b, finding 4)


def test_one_tolerant_target_is_not_exact_only(tmp_path: Path) -> None:
    """Mutant M5: ``targets = 1`` in one placement row: the full plot, not the empty state."""
    inst = _one(mapping_out(tmp_path, fdr=full_placement({0: (1, 0, 0.0, 1)})), "V06")
    assert inst.empty is False
    assert "placed by exact match" not in inst.title
    on = sorted(ax.get_label() for ax in inst.figure.axes if ax.axison)
    assert {"a", "b"} <= set(on), on


def test_degenerate_ci_gets_no_ci_label(tmp_path: Path) -> None:
    """Mutant M2b: a CI of zero width (1 to 1) is no drawn CI."""
    rows = [denovo_row(i, 1.0, 1.0, 1.0) for i in (1, 2, 3)]
    inst = _one(mapping_out(tmp_path, fdr=[*full_placement(EXAMPLE4_TOTALS), *rows]), "V06")
    assert not any("(CI)" in t for t in _ylabels(inst)), _ylabels(inst)


def test_measured_bins_with_na_values_are_not_assumed(tmp_path: Path) -> None:
    """Finding 4: bins not ``c_assumed`` whose values are all NA: no "assumed c" text."""
    rows = [denovo_row(i, None, None, None) for i in (1, 2, 3)]
    inst = _one(mapping_out(tmp_path, fdr=[*full_placement(EXAMPLE4_TOTALS), *rows]), "V06")
    texts = h.joined(inst)
    assert not any("assumed c" in t for t in texts), texts


# ---------------------------------------------------------------- V11 twin


def test_v11_twin_c_title_names_only_the_drawn_classes(tmp_path: Path) -> None:
    inst = _one(h.p_pairs(tmp_path).out, "V11")
    s = _twin("V11", inst)
    (c,) = [t for t in _subplot_titles(s) if t.startswith("(c)")]
    assert re.search(r"\bOP\b", c) and not re.search(r"\bNO\b", c), c


def test_v11_twin_legend_uses_class_label(tmp_path: Path) -> None:
    from gippa.viz.operon_pairs import class_label

    inst = _one(h.p_pairs(tmp_path, measured_no=True).out, "V11")
    s = _twin("V11", inst)
    names = [str(t.get("name", "")) for t in tw.traces(s)]
    for c in ("OP", "NO"):
        assert any(n.startswith(class_label(c)) for n in names), (c, names)
    assert not any(re.match(r"^(OP|NO)\b(?! \()", n) for n in names), names


def test_v11_twin_draws_op_over_no_at_opacity_0_6(tmp_path: Path) -> None:
    inst = _one(h.p_pairs(tmp_path, measured_no=True).out, "V11")
    s = _twin("V11", inst)
    # the (c) class traces, named by their code (legend wording is the test above)
    scatter = [t for t in tw.traces(s, "scatter", "scattergl")
               if re.match(r"^(OP|NO)\b", str(t.get("name", "")))]  # fmt: skip
    order = [str(t["name"]).split(" ")[0] for t in scatter]
    assert "OP" in order and "NO" in order, order
    assert order.index("OP") > order.index("NO"), order
    op = [t for t in scatter if str(t["name"]).startswith("OP")]
    assert all(t.get("opacity") == pytest.approx(0.6) for t in op), op


# ---------------------------------------------------------------- V10 twin and static


def _bb_na_rows(contigs: tuple[str, ...] = (h.CONTIG,)) -> list[dict[str, Any]]:
    """Moran's I estimable; every BB row NA and flagged (panel (b) is not drawn)."""
    rows = []
    for contig in contigs:
        for lag in range(1, 6):
            for null in ("gene_label", "tu_block"):
                rows.append(h.ac_row("morans_i", null, lag, 0.05 * lag, contig=contig))
                rows.append(h.ac_row("join_count_bb", null, lag, None,
                                     "few_values;not_estimable", contig=contig))  # fmt: skip
    return rows


def test_v10_twin_partial_panel_text_is_plain(tmp_path: Path) -> None:
    run = h.d_full(tmp_path)
    h.write_ac(run.out, _bb_na_rows())
    inst = _one(run.out, "V10")
    s = _twin("V10", inst)
    notes = tw.annotations(s)
    assert h.raw_tokens(notes) == [], notes
    assert any(plain("not_estimable") in n for n in notes), notes


def test_v10_flags_mixing_na_and_text(tmp_path: Path) -> None:
    """Mutant M4: NA flag cells beside text ones render; an NA cell is no flag."""
    rows = []
    for lag in range(1, 6):
        for null in ("gene_label", "tu_block"):
            rows.append(h.ac_row("morans_i", null, lag, 0.05 * lag, None))
            rows.append(h.ac_row("join_count_bb", null, lag, 0.2, "ci_unestimated"))
    run = h.d_full(tmp_path)
    h.write_ac(run.out, rows)
    inst = _one(run.out, "V10")
    assert inst.empty is False
    a = next(ax for ax in inst.figure.axes if ax.get_label() == f"a:{h.CONTIG}")
    xs = sorted({float(x) for line in a.get_lines() if line.get_marker() not in ("", "None")
                 for x in line.get_xdata()})  # fmt: skip
    assert xs == [1.0, 2.0, 3.0, 4.0, 5.0]


def test_v10_na_flag_cells_are_no_flag_in_the_not_drawn_note(tmp_path: Path) -> None:
    """Mutant M4: an NA ``flags`` cell is no flag, not a token ("nan") in the panel note."""
    rows = []
    for lag in range(1, 6):
        for null in ("gene_label", "tu_block"):
            rows.append(h.ac_row("morans_i", null, lag, 0.05 * lag, "ci_unestimated"))
            rows.append(h.ac_row("join_count_bb", null, lag, None, None))
    run = h.d_full(tmp_path)
    h.write_ac(run.out, rows)
    inst = _one(run.out, "V10")
    b = next(ax for ax in inst.figure.axes if ax.get_label() == f"b:{h.CONTIG}")
    notes = [t.get_text() for t in b.texts if t.get_text().startswith("Not drawn")]
    assert notes, [t.get_text() for t in b.texts]
    assert not any(re.search(r"\bnan\b|\bNA\b|None", n) for n in notes), notes


@pytest.mark.parametrize("contigs", [(h.CONTIG,), (h.CONTIG, "chr2", "chr3")])
def test_v10_not_drawn_note_stays_inside_its_axes(tmp_path: Path, contigs: tuple[str, ...]) -> None:
    """The note "Not drawn: too few values to scan; statistic cannot be estimated" is not
    clipped at the edges of its axes, also in the narrow columns of three replicons."""
    run = h.d_full(tmp_path)
    h.write_ac(run.out, _bb_na_rows(contigs))
    inst = _one(run.out, "V10")
    from matplotlib.backends.backend_agg import FigureCanvasAgg

    fig = inst.figure
    canvas = FigureCanvasAgg(fig)
    canvas.draw()
    renderer = canvas.get_renderer()
    for contig in contigs:
        b = next(ax for ax in fig.axes if ax.get_label() == f"b:{contig}")
        notes = [t for t in b.texts if t.get_text().startswith("Not drawn")]
        assert notes, [t.get_text() for t in b.texts]
        box = b.get_window_extent(renderer)
        for t in notes:
            ext = t.get_window_extent(renderer)
            assert ext.x0 >= box.x0 - 1 and ext.x1 <= box.x1 + 1, (contig, t.get_text(), ext, box)


# ---------------------------------------------------------------- V03 report table (M6b, M6c)


def _v03_page(fixture: str, tmp: Path) -> tuple[Any, list[dict[str, str]]]:
    from gippa.io.plot_files import read_plot_files
    from gippa.io.run_view import RunView
    from gippa.report.build import build_report
    from gippa.steps.plots import load_plot_modules, render_plots, skipped_plots

    why = sf.missing(fixture)
    if why is not None:
        pytest.fail(why)
    load_plot_modules()
    plots = tmp / "plots"
    view = RunView.from_outdir(sf.fixture_out(fixture))
    cfg = sf.fixture_cfg(fixture)
    sf.styled(lambda: render_plots(view, cfg, plots, ["V03"]))
    page = parse(build_report(view, read_plot_files(plots), skipped_plots(view), cfg))
    with (plots / "data" / "V03_locus.tsv").open(encoding="utf-8") as fh:
        rows = list(csv.DictReader([ln for ln in fh if not ln.startswith("#")], delimiter="\t"))
    return page, rows


ESTIMATE = {"island": "FWER p", "missed_operon": "decoy-FDR q"}


def test_v03_report_table_round_trips_every_column(tmp_path: Path) -> None:
    page, rows = _v03_page(sf.RICH, tmp_path)
    tables = [t for t in page.root.by_tag("table") if t.attrs.get("data-index") == "V03"]
    assert len(tables) == 1
    header, body = table_rows(tables[0])
    assert [c.lower() for c in header] == ["rank", "locus", "family", "estimate", "q"]
    assert len(body) == len(rows) > 0
    for cells, r in zip(body, rows, strict=True):
        rank, locus, family, estimate, q = cells
        assert rank == str(int(r["rank"]))
        assert locus == r["locus_id"]
        assert "_" not in family and family == r["family"].replace("_", " "), (family, r)
        assert estimate == ESTIMATE.get(r["family"], "BH q"), (estimate, r)
        assert q == style.fmt_q(None if r["error"] in ("", "NA") else float(r["error"])), (q, r)


FAMILIES = ("prophage", "island", "bgc", "operon_discordance", "missed_operon")


def test_v03_report_rows_name_each_family_in_words() -> None:
    """Mutants M6b/M6c on families the rich fixture lacks (underscores, decoy-FDR q)."""
    from gippa.report.content import locus_index
    from gippa.viz.locus import INDEX_COLUMNS

    header = list(INDEX_COLUMNS)
    rows = []
    for i, fam in enumerate(FAMILIES, start=1):
        est = ESTIMATE.get(fam, "BH q")
        cells = {"rank": str(i), "locus_id": f"L{i}", "kind": fam, "family": fam,
                 "contig": "chr1", "start": "1", "end": "10",
                 "estimate": est, "error": "0.0123", "rendered": "true"}  # fmt: skip
        rows.append([cells[c] for c in header])
    out = locus_index(header, rows)
    assert [r[2] for r in out] == [f.replace("_", " ") for f in FAMILIES]
    assert [r[3] for r in out] == [ESTIMATE.get(f, "BH q") for f in FAMILIES]
    assert [r[4] for r in out] == [style.fmt_q(0.0123)] * len(FAMILIES)
    assert [r[:2] for r in out] == [[str(i), f"L{i}"] for i in range(1, len(FAMILIES) + 1)]


def test_v03_zero_row_index_is_a_report_note(tmp_path: Path) -> None:
    page, rows = _v03_page(sf.EMPTY, tmp_path)
    assert rows == []
    notes = page.root.by_attr("data-index-none", "V03")
    assert len(notes) == 1
    assert "Locus views, in display order" in notes[0].text()
    assert [t for t in page.root.by_tag("table") if t.attrs.get("data-index") == "V03"] == []


# ---------------------------------------------------------------- V15 (M13)


def test_v15_without_features_says_tus_not_recorded(tmp_path: Path) -> None:
    out = make_out(tmp_path)
    write_not_fitted(out, 99)
    result = render(out, ["V15"])["V15"]
    assert not isinstance(result, base.Skip), result
    (inst,) = result
    assert inst.empty is True
    assert inst.title.startswith("Nothing to show: ")
    assert "TUs not recorded" in inst.title, inst.title
    assert "NA TUs" not in inst.title
