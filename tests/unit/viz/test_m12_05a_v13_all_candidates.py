"""M12-05a: V13 draws every novel ORF candidate, lowest q first
(docs/specs/M12-05a-novel-orf-plot.md AC 1-9 and 13, examples E1-E10, and the lead decisions
of 2026-11-09).

Lead decision Q1 ("class q, then peptide q"): the order is ``class_fdr`` ascending, ties by
``best_peptide_q`` ascending, then called first, ``rank``, ``rank_key``. The order oracle and the
fixtures live in ``tests/helpers/m12_05a.py``. Lead decision (from the M12-05c review): no raw
class identifier (``nterm_extension``, ``cterm_extension``, ``cds_outframe``, ``sorf``) appears in
the text of any V13 SVG.

Interface pinned here beyond the spec text:

- the ORF arrow is the one ``Polygon`` on the ``"genes"`` axes whose edge is ``#CC79A7``; its
  dash pattern is read from ``Patch._unscaled_dash_pattern`` (``(0, None)`` solid, ``(0, (3,
  2))`` for ``(0, (3, 2))``), since a legend proxy reports ``get_linestyle() == "solid"`` whatever
  it copied;
- the index image has one ``Text`` per shown row starting with ``"<order>. <orf_id>"``; the
  other three columns are ``Text`` objects on the same line (or the same ``Text``).
"""

from __future__ import annotations

import importlib
import math
import tempfile
import xml.etree.ElementTree as ET
from collections.abc import Iterator, Sequence
from pathlib import Path
from types import ModuleType
from typing import Any

import matplotlib.colors as mcolors
import pandas as pd
import pytest
from helpers import m12_05a as h
from helpers import snapshot_fixture as sf
from helpers import viz_locus as vl
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st
from matplotlib.axes import Axes
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from matplotlib.legend import Legend
from matplotlib.patches import Patch, Polygon
from matplotlib.text import Text

import gippa.viz.locus as locus_mod
from gippa.config import GippaConfig, resolve_config
from gippa.io.plot_files import write_instance
from gippa.io.run_view import RunView
from gippa.viz import base, style
from gippa.viz.base import Instance, Skip

LOCUS_COLUMNS = ["track", "id", "contig", "start", "end", "strand", "value", "category", "fill"]
COUNT_COLUMNS = ["class", "status", "n", "class_fdr_max_called"]


@pytest.fixture(autouse=True)
def _style() -> None:
    style.apply()


@pytest.fixture
def cfg() -> GippaConfig:
    return resolve_config()


@pytest.fixture
def nov() -> ModuleType:
    return importlib.import_module("gippa.viz.novel_orfs")


# --- helpers ------------------------------------------------------------------------------------


def _run_view(view: RunView, cfg: GippaConfig) -> list[Instance]:
    importlib.import_module("gippa.viz.novel_orfs")
    result = base.registered()["V13"].fn(view, cfg)
    assert not isinstance(result, Skip), result
    return list(result)


def _run(
    tmp: Path, rows: Sequence[dict[str, object]], cfg: GippaConfig, **kw: Any
) -> list[Instance]:
    return _run_view(RunView.from_outdir(vl.build_out(tmp, novel_orfs=rows, **kw)), cfg)


def _main(insts: Sequence[Instance]) -> Instance:
    (m,) = [i for i in insts if i.suffix is None]
    return m


def _views(insts: Sequence[Instance]) -> list[Instance]:
    return [i for i in insts if i.suffix is not None and i.suffix.startswith("locus_")]


def _view_ids(insts: Sequence[Instance]) -> list[str]:
    return [str(v.suffix).removeprefix("locus_") for v in _views(insts)]


def _index(insts: Sequence[Instance]) -> Instance:
    found = [i for i in insts if i.suffix == "locus"]
    assert len(found) == 1, [i.suffix for i in insts]
    return found[0]


def _axes(fig: Figure, label: str) -> Axes:
    found = [ax for ax in fig.axes if ax.get_label() == label]
    assert len(found) == 1, [ax.get_label() for ax in fig.axes]
    return found[0]


def _agg(fig: Figure) -> Any:
    renderer = FigureCanvasAgg(fig).get_renderer()
    fig.draw(renderer)
    return renderer


def _hex(c: Any) -> str | None:
    try:
        return mcolors.to_hex(c).upper()
    except (ValueError, TypeError):
        return None


def _collapse(s: str) -> str:
    return " ".join(s.split())


def _texts(fig: Figure) -> list[Text]:
    return [t for t in fig.findobj(Text) if t.get_visible() and t.get_text().strip()]


def _arrow(fig: Figure) -> Polygon:
    ax = _axes(fig, "genes")
    found = [
        p for p in ax.patches if isinstance(p, Polygon) and _hex(p.get_edgecolor()) == h.ORF_EDGE
    ]
    assert len(found) == 1, len(found)
    return found[0]


def _dashes(p: Patch) -> tuple[float, tuple[float, ...] | None]:
    off, seq = p._unscaled_dash_pattern  # type: ignore[attr-defined]
    return float(off), None if seq is None else tuple(float(x) for x in seq)


def _legends(fig: Figure) -> list[Legend]:
    out = list(fig.legends)
    out += [lg for ax in fig.axes if (lg := ax.get_legend()) is not None]
    return out


class _Stub:
    """A cheap ``render_locus``: records the loci, draws nothing (order tests only)."""

    def __init__(self) -> None:
        self.loci: list[locus_mod.Locus] = []

    def __call__(self, lo: locus_mod.Locus, *args: Any, **kw: Any) -> Instance:
        self.loci.append(lo)
        data = pd.DataFrame({c: pd.Series([], dtype=object) for c in LOCUS_COLUMNS})
        return Instance(lo.locus_id, style.new_figure("single", 1.0), data, lo.label, ("start",))


def _patch_stub(mp: pytest.MonkeyPatch, nov: ModuleType) -> _Stub:
    s = _Stub()
    mp.setattr(locus_mod, "render_locus", s)
    if hasattr(nov, "render_locus"):
        mp.setattr(nov, "render_locus", s)
    return s


@pytest.fixture
def stub(nov: ModuleType, monkeypatch: pytest.MonkeyPatch) -> Iterator[_Stub]:
    yield _patch_stub(monkeypatch, nov)


def _index_tsv(ix: Instance, tmp: Path) -> list[list[str]]:
    paths = write_instance(ix, "V13", "novel_orfs", tmp, ["svg"], 100)
    (tsv,) = [p for p in paths if p.suffix == ".tsv"]
    assert tsv.name == "V13_novel_orfs_locus.tsv"
    return [ln.split("\t") for ln in tsv.read_text(encoding="utf-8").splitlines()]


# --- AC 1: every candidate ------------------------------------------------------------------------


def test_ac1_e1_route_e_gives_panel_a_8_views_and_the_index(
    tmp_path: Path, cfg: GippaConfig, stub: _Stub
) -> None:
    insts = _run(tmp_path, h.route_e_rows(), cfg)
    assert len(insts) == 10
    assert insts[0].suffix is None
    assert len(_views(insts)) == 8
    _index(insts)


def test_ac1_e8_25_uncalled_rows_give_22_instances(
    tmp_path: Path, cfg: GippaConfig, stub: _Stub
) -> None:
    insts = _run(tmp_path, h.e8_rows(), cfg)
    assert len(insts) == 22
    assert len(_views(insts)) == 20
    _index(insts)


def test_ac1_e10_header_only_is_the_unchanged_empty_state(tmp_path: Path, cfg: GippaConfig) -> None:
    insts = _run(tmp_path / "run", [], cfg)
    assert len(insts) == 1
    (main,) = insts
    assert main.empty and main.suffix is None
    assert main.title.startswith("Nothing to show: no novel ORF candidates in this run. ")
    paths = write_instance(main, "V13", "novel_orfs", tmp_path / "plots", ["svg"], 100)
    assert not [p for p in paths if p.name.startswith("V13_novel_orfs_locus")]


@pytest.mark.parametrize("called", [None, False])
def test_ac1_called_na_or_false_still_gets_a_view(
    tmp_path: Path, cfg: GippaConfig, called: bool | None, stub: _Stub
) -> None:
    rows = [h.orf("NOVEL_chr1_5001_p", called=called, rank=None, class_fdr="0.3")]
    insts = _run(tmp_path, rows, cfg)
    assert _view_ids(insts) == ["NOVEL_chr1_5001_p"]
    assert len(insts) == 3


# --- AC 2: order ----------------------------------------------------------------------------------


def test_ac2_e1_order(tmp_path: Path, cfg: GippaConfig, stub: _Stub) -> None:
    assert _view_ids(_run(tmp_path, h.route_e_rows(), cfg)) == h.E1_ORDER


@pytest.mark.parametrize("example", ["e4", "e5", "e6", "e7"])
@pytest.mark.parametrize("reverse", [False, True])
def test_ac2_tie_examples(
    tmp_path: Path, cfg: GippaConfig, stub: _Stub, example: str, reverse: bool
) -> None:
    rows, want = getattr(h, f"{example}_rows")()
    rows = rows[::-1] if reverse else rows
    assert _view_ids(_run(tmp_path, rows, cfg)) == want


def _rich() -> RunView:
    why = sf.missing(sf.RICH)
    if why:
        pytest.fail(why)
    return RunView.from_outdir(sf.fixture_out(sf.RICH))


#: spec E2: 13 pseudogenes (q 0.00332862) in rank order, then 7 antisense (q 0.00383386)
E2_FIRST_20 = [
    "NOVEL_chr1_45143_m", "NOVEL_chr1_44225_m", "NOVEL_chr1_269313_p", "NOVEL_chr1_249144_m",
    "NOVEL_chr1_106664_m", "NOVEL_chr1_6789_m", "NOVEL_chr1_268776_p", "NOVEL_chr1_7293_m",
    "NOVEL_chr1_249558_m", "NOVEL_chr1_106406_m", "NOVEL_chr1_92882_p", "NOVEL_chr1_93092_p",
    "NOVEL_chr1_269094_p", "NOVEL_chr1_15774_p", "NOVEL_chr1_176135_p", "NOVEL_chr1_208659_p",
    "NOVEL_chr1_276721_m", "NOVEL_chr1_78659_m", "NOVEL_chr1_213777_m", "NOVEL_chr1_110417_m",
]  # fmt: skip


def test_ac2_e2_sim_plots_first_20_views(cfg: GippaConfig, stub: _Stub) -> None:
    view = _rich()
    insts = _run_view(view, cfg)
    ids = _view_ids(insts)
    assert ids == E2_FIRST_20
    table = view.table("novel_orfs")
    assert table is not None
    ranks = {str(r.orf_id): r.rank for r in table.itertuples(index=False)}
    assert [int(ranks[i]) for i in ids[:13]] == [1, 2, 3, 4, 7, 9, 13, 17, 18, 20, 23, 24, 32]
    assert [int(ranks[i]) for i in ids[13:]] == [5, 6, 8, 10, 11, 12, 14]


def test_ac2_e3_rank_32_is_drawn_and_ranks_15_16_19_are_not(
    cfg: GippaConfig, stub: _Stub, tmp_path: Path
) -> None:
    insts = _run_view(_rich(), cfg)
    ids = _view_ids(insts)
    assert ids.index("NOVEL_chr1_269094_p") == 12  # view 13
    lines = _index_tsv(_index(insts), tmp_path)
    col = {c: i for i, c in enumerate(lines[0])}
    rendered = {r[col["orf_id"]]: r[col["rendered"]] for r in lines[1:]}
    for oid in ("NOVEL_chr1_185480_p", "NOVEL_chr1_247263_p", "NOVEL_chr1_130084_m"):
        assert rendered[oid] == "false", oid


def _row_strategy() -> st.SearchStrategy[dict[str, object]]:
    return st.fixed_dictionaries(
        {
            "class_fdr": st.sampled_from(["0.01", "0.02", "0.4", None, "unestimated"]),
            "best_q": st.sampled_from([0.001, 0.002, None]),
            "called": st.sampled_from([True, False, None]),
            "rank": st.sampled_from([1, 2, 3, None]),
            "n_peptides": st.sampled_from([1, 2]),
        }
    )


@settings(
    max_examples=12, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture]
)
@given(spec=st.lists(_row_strategy(), min_size=1, max_size=30), seed=st.randoms())
def test_ac2_property_order_is_the_oracle_and_ignores_row_order(
    cfg: GippaConfig, nov: ModuleType, spec: list[dict[str, Any]], seed: Any
) -> None:
    rows = []
    for i, s in enumerate(spec):
        start = 1_001 + 3_000 * i
        rows.append(h.orf(f"NOVEL_chr1_{start}_p", start_longest=start, end=start + 299,
                          class_fdr=s["class_fdr"], best_q=s["best_q"], called=s["called"],
                          rank=s["rank"], n_peptides=s["n_peptides"]))  # fmt: skip
    shuffled = list(rows)
    seed.shuffle(shuffled)
    want = h.expected_order(rows)
    got = []
    with tempfile.TemporaryDirectory() as tmp, pytest.MonkeyPatch.context() as mp:
        _patch_stub(mp, nov)
        for k, table in enumerate((rows, shuffled)):
            insts = _run(Path(tmp) / str(k), table, cfg)
            ids = _view_ids(insts)
            lines = _index_tsv(_index(insts), Path(tmp) / f"plots{k}")
            col = {c: i for i, c in enumerate(lines[0])}
            got.append((ids, [r[col["orf_id"]] for r in lines[1:]]))
    assert got[0] == got[1]
    ids, listed = got[0]
    assert listed == want
    assert ids == want[: h.TOP_N]
    by_id = {str(r["orf_id"]): r for r in rows}
    keys = [h.order_key(by_id[i])[:2] for i in ids]
    assert keys == sorted(keys)  # non-decreasing in (class q, peptide q), NA as +inf


# --- AC 3: style ----------------------------------------------------------------------------------


def _orf_entries(fig: Figure) -> list[tuple[Legend, int, str]]:
    return [
        (lg, i, _collapse(t.get_text()))
        for lg in _legends(fig)
        for i, t in enumerate(lg.get_texts())
        if "novel ORF" in t.get_text()
    ]


@pytest.mark.parametrize("called", [True, False])
def test_ac3_arrow_style_follows_the_call(tmp_path: Path, cfg: GippaConfig, called: bool) -> None:
    rows = [h.orf("NOVEL_chr1_5001_p", called=called, rank=1 if called else None)]
    (view,) = _views(_run(tmp_path, rows, cfg))
    _agg(view.figure)
    arrow = _arrow(view.figure)
    face = mcolors.to_rgba(arrow.get_facecolor())
    if called:
        assert face[3] == pytest.approx(0.35)
        assert _hex(face[:3]) == h.ORF_EDGE
        assert _dashes(arrow)[1] is None
    else:
        assert face[3] == 0
        assert _dashes(arrow) == (0.0, (3.0, 2.0))
    assert arrow.get_linewidth() == pytest.approx(1.0)


@pytest.mark.parametrize("called", [True, False])
def test_ac3_one_key_entry_matching_the_arrow(
    tmp_path: Path, cfg: GippaConfig, called: bool
) -> None:
    rows = [h.orf("NOVEL_chr1_5001_p", called=called, rank=1 if called else None)]
    (view,) = _views(_run(tmp_path, rows, cfg))
    _agg(view.figure)
    entries = _orf_entries(view.figure)
    assert len(entries) == 1, [e[2] for e in entries]
    lg, i, text = entries[0]
    assert text == (h.ORF_KEY_CALLED if called else h.ORF_KEY_UNCALLED)
    handle = lg.legend_handles[i]
    assert isinstance(handle, Patch)
    arrow = _arrow(view.figure)
    assert mcolors.to_rgba(handle.get_facecolor()) == pytest.approx(
        mcolors.to_rgba(arrow.get_facecolor())
    )
    assert _hex(handle.get_edgecolor()) == h.ORF_EDGE
    assert _dashes(handle) == _dashes(arrow)


def test_ac3_key_texts_are_module_constants() -> None:
    assert getattr(locus_mod, "ORF_KEY_CALLED", None) == h.ORF_KEY_CALLED
    assert getattr(locus_mod, "ORF_KEY_UNCALLED", None) == h.ORF_KEY_UNCALLED


def test_ac3_e1_every_view_is_dashed_and_unfilled(tmp_path: Path, cfg: GippaConfig) -> None:
    views = _views(_run(tmp_path, h.route_e_rows(), cfg))
    assert len(views) == 8
    for v in views:
        arrow = _arrow(v.figure)
        assert mcolors.to_rgba(arrow.get_facecolor())[3] == 0, v.suffix
        assert _dashes(arrow) == (0.0, (3.0, 2.0)), v.suffix


def test_ac3_e7_na_q_called_is_filled_and_unestimated_uncalled_is_dashed(
    tmp_path: Path, cfg: GippaConfig
) -> None:
    rows, _ = h.e7_rows()
    views = {str(v.suffix).removeprefix("locus_"): v for v in _views(_run(tmp_path, rows, cfg))}
    assert set(views) == {h.E7_E, h.E7_F, h.E7_G}
    assert "class FDR unestimated" in _collapse(views[h.E7_E].title)
    assert mcolors.to_rgba(_arrow(views[h.E7_E].figure).get_facecolor())[3] == pytest.approx(0.35)
    assert _dashes(_arrow(views[h.E7_G].figure)) == (0.0, (3.0, 2.0))


# --- AC 4: title ----------------------------------------------------------------------------------


def test_ac4_e9_title_says_not_called(tmp_path: Path, cfg: GippaConfig) -> None:
    (view,) = _views(_run(tmp_path, [h.e9_row()], cfg, replicons=h.C10))
    want = (f"{h.E9_ID} (antisense, not called): 1 peptide, best q {style.fmt_q(0.002)}, "
            f"class FDR {style.fmt_num(0.4)}")  # fmt: skip
    assert _collapse(view.title) == want
    assert want in _collapse(" ".join(t.get_text() for t in _texts(view.figure)))


def test_ac4_called_title_says_called(tmp_path: Path, cfg: GippaConfig) -> None:
    (view,) = _views(_run(tmp_path, [h.e9_row(called=True)], cfg, replicons=h.C10))
    assert _collapse(view.title).startswith(f"{h.E9_ID} (antisense, called): 1 peptide, ")


# --- AC 5: the index image and its cap note -------------------------------------------------------


def test_ac5_e2_index_title(cfg: GippaConfig, stub: _Stub) -> None:
    ix = _index(_run_view(_rich(), cfg))
    want = f"{h.INDEX_TITLE} (top 20 of 57; all in tables/novel_orfs.tsv)"
    assert _collapse(ix.title) == want
    assert any(want in _collapse(t.get_text()) for t in _texts(ix.figure))


def test_ac5_e1_index_title(tmp_path: Path, cfg: GippaConfig, stub: _Stub) -> None:
    ix = _index(_run(tmp_path, h.route_e_rows(), cfg))
    assert _collapse(ix.title) == f"{h.INDEX_TITLE} (8 candidates)"


def test_ac5_e8_index_title(tmp_path: Path, cfg: GippaConfig, stub: _Stub) -> None:
    ix = _index(_run(tmp_path, h.e8_rows(), cfg))
    assert _collapse(ix.title) == f"{h.INDEX_TITLE} (top 20 of 25; all in tables/novel_orfs.tsv)"


def test_ac5_e9_index_title_is_singular(tmp_path: Path, cfg: GippaConfig) -> None:
    insts = _run(tmp_path, [h.e9_row()], cfg, replicons=h.C10)
    assert len(insts) == 3
    assert _collapse(_index(insts).title) == f"{h.INDEX_TITLE} (1 candidate)"


def _line_texts(fig: Figure) -> dict[int, str]:
    """order number -> every text on that index line, joined."""
    r = _agg(fig)
    texts = _texts(fig)
    out: dict[int, str] = {}
    for t in texts:
        head = t.get_text().strip().split(".", 1)[0]
        if not head.isdigit() or not t.get_text().strip().startswith(f"{head}. "):
            continue
        ext = t.get_window_extent(r)
        y = (ext.y0 + ext.y1) / 2
        same = [u for u in texts
                if (e := u.get_window_extent(r)).y0 - 0.5 <= y <= e.y1 + 0.5]  # fmt: skip
        out[int(head)] = " | ".join(_collapse(u.get_text()) for u in same)
    return out


def test_ac5_e1_index_lines(tmp_path: Path, cfg: GippaConfig, stub: _Stub) -> None:
    rows = {str(r["orf_id"]): r for r in h.route_e_rows()}
    ix = _index(_run(tmp_path, list(rows.values()), cfg))
    assert ix.figure.get_figwidth() == pytest.approx(style.DOUBLE_IN)
    lines = _line_texts(ix.figure)
    assert sorted(lines) == list(range(1, 9))
    for k, oid in enumerate(h.E1_ORDER, start=1):
        line = lines[k]
        assert f"{k}. {oid}" in line, line
        assert "not called" in line, line
        assert f"class FDR {style.fmt_num(float(str(rows[oid]['class_fdr'])))}" in line, line
        cls = str(rows[oid]["class"])
        if cls in ("antisense", "intergenic"):
            assert cls in line, line
        assert not h.raw_class_tokens(line), line


def test_ac5_index_lines_are_capped_at_20(tmp_path: Path, cfg: GippaConfig, stub: _Stub) -> None:
    ix = _index(_run(tmp_path, h.e8_rows(), cfg))
    assert sorted(_line_texts(ix.figure)) == list(range(1, 21))


def test_ac5_index_line_of_a_called_unestimated_row(
    tmp_path: Path, cfg: GippaConfig, stub: _Stub
) -> None:
    rows, _ = h.e7_rows()
    lines = _line_texts(_index(_run(tmp_path, rows, cfg)).figure)
    assert f"2. {h.E7_E}" in lines[2]
    assert "class FDR unestimated" in lines[2]
    assert "called" in lines[2] and "not called" not in lines[2], lines[2]
    assert "not called" in lines[1] and "not called" in lines[3]


# --- AC 6: the index TSV --------------------------------------------------------------------------


def test_ac6_e1_index_tsv(tmp_path: Path, cfg: GippaConfig, stub: _Stub) -> None:
    lines = _index_tsv(_index(_run(tmp_path / "run", h.route_e_rows(), cfg)), tmp_path / "plots")
    assert lines[0] == h.INDEX_COLUMNS
    body = [dict(zip(lines[0], r, strict=True)) for r in lines[1:]]
    assert [r["order"] for r in body] == [str(i) for i in range(1, 9)]
    assert [r["orf_id"] for r in body] == h.E1_ORDER
    assert {r["called"] for r in body} == {"false"}
    assert {r["rendered"] for r in body} == {"true"}
    assert {r["rank"] for r in body} == {"NA"}
    assert body[0]["class"] == "antisense" and body[5]["class"] == "nterm_extension"
    assert body[0]["class_fdr"] == "0.4" and body[5]["class_fdr"] == "1"
    assert body[0]["best_peptide_q"] == "0.000874751"
    assert body[0]["n_peptides"] == "1"


def test_ac6_e2_index_tsv_lists_all_57(cfg: GippaConfig, stub: _Stub, tmp_path: Path) -> None:
    lines = _index_tsv(_index(_run_view(_rich(), cfg)), tmp_path)
    assert lines[0] == h.INDEX_COLUMNS
    body = [dict(zip(lines[0], r, strict=True)) for r in lines[1:]]
    assert len(body) == 57
    assert [r["order"] for r in body] == [str(i) for i in range(1, 58)]
    assert [r["rendered"] for r in body] == ["true"] * 20 + ["false"] * 37
    assert [r["orf_id"] for r in body[:20]] == E2_FIRST_20
    r39 = body[38]
    assert (r39["orf_id"], r39["class"], r39["called"]) == ("NOVEL_chr1_48669_p", "antisense",
                                                            "false")  # fmt: skip
    assert float(r39["class_fdr"]) == pytest.approx(0.00383386)
    assert all(r["called"] == "true" for r in body[:38])
    assert all(r["called"] == "false" for r in body[38:])
    assert body[56]["orf_id"] == "NOVEL_chr1_74837_p"


def test_ac6_e7_class_fdr_is_na_for_no_q(tmp_path: Path, cfg: GippaConfig, stub: _Stub) -> None:
    rows, want = h.e7_rows()
    lines = _index_tsv(_index(_run(tmp_path / "run", rows, cfg)), tmp_path / "plots")
    body = [dict(zip(lines[0], r, strict=True)) for r in lines[1:]]
    assert [r["orf_id"] for r in body] == want
    fdr = {r["orf_id"]: r["class_fdr"] for r in body}
    assert fdr[h.E7_E] == "NA" and fdr[h.E7_G] == "NA"
    assert fdr[h.E7_F] == "0.9"
    called = {r["orf_id"]: r["called"] for r in body}
    assert called == {h.E7_E: "true", h.E7_F: "false", h.E7_G: "false"}
    rank = {r["orf_id"]: r["rank"] for r in body}
    assert rank == {h.E7_E: "1", h.E7_F: "NA", h.E7_G: "NA"}


def test_ac6_e8_index_tsv_has_25_rows_5_not_rendered(
    tmp_path: Path, cfg: GippaConfig, stub: _Stub
) -> None:
    lines = _index_tsv(_index(_run(tmp_path / "run", h.e8_rows(), cfg)), tmp_path / "plots")
    body = [dict(zip(lines[0], r, strict=True)) for r in lines[1:]]
    assert len(body) == 25
    assert [r["rendered"] for r in body] == ["true"] * 20 + ["false"] * 5
    assert [r["orf_id"] for r in body] == [f"NOVEL_chr1_{1_001 + 3_000 * i}_p" for i in range(25)]


def test_ac6_called_na_is_false_in_the_index(tmp_path: Path, cfg: GippaConfig, stub: _Stub) -> None:
    rows = [h.orf("NOVEL_chr1_5001_p", called=None, rank=None)]
    lines = _index_tsv(_index(_run(tmp_path / "run", rows, cfg)), tmp_path / "plots")
    assert dict(zip(lines[0], lines[1], strict=True))["called"] == "false"


# --- AC 7: the locus TSV's category cell ----------------------------------------------------------


def _locus_tsv(tmp: Path, cfg: GippaConfig, row: dict[str, object], **kw: Any) -> list[list[str]]:
    (view,) = _views(_run(tmp / "run", [row], cfg, **kw))
    paths = write_instance(view, "V13", "novel_orfs", tmp / "plots", ["svg"], 100)
    (tsv,) = [p for p in paths if p.suffix == ".tsv"]
    return [ln.split("\t") for ln in tsv.read_text(encoding="utf-8").splitlines()]


@pytest.mark.parametrize("example", ["e9", "example3"])
def test_ac7_only_the_novel_orf_category_differs(
    tmp_path: Path, cfg: GippaConfig, example: str
) -> None:
    def row(called: bool) -> tuple[dict[str, object], dict[str, Any]]:
        if example == "e9":
            return h.e9_row(called=called), {"replicons": h.C10}
        r = h.orf("NOVEL_chr1_1201_m", cls="antisense", strand="-", start_longest=1_500,
                  end=1_201, called=called, rank=1 if called else None)  # fmt: skip
        hits = [vl.hit("PEPTWDEK", "chr1", 1_301, 1_324, "-")]
        genes = [vl.gene("g1", "chr1", 1_001, 1_700, strand="+", log2=20.0)]
        return r, {"hits": hits, "genes": genes}

    r_yes, kw_yes = row(True)
    r_no, kw_no = row(False)
    yes = _locus_tsv(tmp_path / "yes", cfg, r_yes, **kw_yes)
    no = _locus_tsv(tmp_path / "no", cfg, r_no, **kw_no)
    assert yes[0] == LOCUS_COLUMNS == no[0]
    cat = LOCUS_COLUMNS.index("category")
    nv_yes = [r for r in yes[1:] if r[0] == "novel_orf"]
    nv_no = [r for r in no[1:] if r[0] == "novel_orf"]
    assert [r[cat] for r in nv_yes] == ["called"]
    assert [r[cat] for r in nv_no] == ["not_called"]
    assert [r for r in yes[1:] if r[0] != "novel_orf"] == [r for r in no[1:] if r[0] != "novel_orf"]
    strip = [r[:cat] + r[cat + 1 :] for r in nv_yes]
    assert strip == [r[:cat] + r[cat + 1 :] for r in nv_no]


# --- AC 8: geometry does not depend on the call ---------------------------------------------------


@pytest.mark.parametrize("example", ["e9", "example3"])
def test_ac8_geometry_equal_called_and_not(tmp_path: Path, cfg: GippaConfig, example: str) -> None:
    def view(called: bool) -> Instance:
        if example == "e9":
            rows, kw = [h.e9_row(called=called)], {"replicons": h.C10}
        else:
            rows = [h.orf("NOVEL_chr1_1201_m", cls="antisense", strand="-", start_longest=1_500,
                          end=1_201, start_best_ostir=1_401, called=called,
                          rank=1 if called else None)]  # fmt: skip
            kw = {"genes": [vl.gene("g1", "chr1", 1_001, 1_700, strand="+", log2=20.0)]}
        (v,) = _views(_run(tmp_path / str(called), rows, cfg, **kw))
        _agg(v.figure)
        return v

    a, b = view(True), view(False)
    assert _arrow(a.figure).get_xy().tolist() == _arrow(b.figure).get_xy().tolist()
    ga, gb = _axes(a.figure, "genes"), _axes(b.figure, "genes")
    assert ga.get_xlim() == gb.get_xlim()
    assert [t.get_text() for t in ga.get_xticklabels()] == [
        t.get_text() for t in gb.get_xticklabels()
    ]
    assert _dashes(_arrow(a.figure)) != _dashes(_arrow(b.figure))  # only the style differs


def test_ac8_e9_window_is_unwrapped_across_the_origin(tmp_path: Path, cfg: GippaConfig) -> None:
    (v,) = _views(_run(tmp_path, [h.e9_row()], cfg, replicons=h.C10))
    nv = v.data[v.data["track"] == "novel_orf"]
    assert [(int(r.start), int(r.end), r.strand) for r in nv.itertuples()] == [(9_899, 10_150, "-")]


# --- AC 9: panel (a) values of e6f4aa3c; the e1 title drops "0 candidates unestimated" (Khoa) ---

_NAN = math.nan
PANEL_A = {
    "e1": ("No novel ORF called",
           [("nterm_extension", "not_called", 1, _NAN), ("intergenic", "not_called", 1, _NAN),
            ("cds_outframe", "not_called", 1, _NAN), ("antisense", "not_called", 5, _NAN)]),
    "e7": ("1 novel ORF called (class FDR unestimated), 1 unestimated",
           [("intergenic", "called", 1, _NAN), ("intergenic", "not_called", 1, _NAN),
            ("intergenic", "unestimated", 1, _NAN)]),
    "e2": (f"38 novel ORFs called (class FDR ≤ {style.fmt_num(0.00383386)}), 0 unestimated",
           [("nterm_extension", "not_called", 2, _NAN), ("intergenic", "not_called", 5, _NAN),
            ("cds_outframe", "not_called", 7, _NAN), ("antisense", "called", 25, 0.00383386),
            ("antisense", "not_called", 2, _NAN), ("pseudogene", "called", 13, 0.00332862),
            ("sorf", "not_called", 3, _NAN)]),
}  # fmt: skip


@pytest.mark.parametrize("example", ["e1", "e2", "e7"])
def test_ac9_panel_a_data_and_title_unchanged(
    tmp_path: Path, cfg: GippaConfig, stub: _Stub, example: str
) -> None:
    if example == "e2":
        insts = _run_view(_rich(), cfg)
    else:
        rows = h.route_e_rows() if example == "e1" else h.e7_rows()[0]
        insts = _run(tmp_path, rows, cfg)
    main = _main(insts)
    title, want = PANEL_A[example]
    assert _collapse(main.title) == title
    assert list(main.data.columns) == COUNT_COLUMNS
    got = [(str(r["class"]), str(r["status"]), int(r["n"]), float(r["class_fdr_max_called"]))
           for _, r in main.data.iterrows() if r["n"]]  # fmt: skip
    assert [g[:3] for g in got] == [w[:3] for w in want]
    for g, w in zip(got, want, strict=True):
        assert (math.isnan(g[3]) and math.isnan(w[3])) or g[3] == pytest.approx(w[3])


# --- lead decision 2026-11-09: no raw class identifier in V13 SVG text ----------------------------


def _svg_text(path: Path) -> str:
    return " ".join(" ".join(ET.parse(path).getroot().itertext()).split())


def test_no_raw_class_identifier_in_any_v13_svg(tmp_path: Path, cfg: GippaConfig) -> None:
    classes = ("nterm_extension", "cterm_extension", "intergenic", "cds_outframe", "antisense",
               "pseudogene", "sorf")  # fmt: skip
    rows = []
    for i, cls in enumerate(classes):
        for called in (True, False):
            s = 1_001 + 6_000 * i + (3_000 if called else 0)
            rows.append(h.orf(f"NOVEL_chr1_{s}_p", cls=cls, start_longest=s, end=s + 299,
                              called=called, rank=i + 1 if called else None))  # fmt: skip
    insts = _run(tmp_path / "run", rows, cfg)
    assert len(_views(insts)) == len(rows)
    found = {}
    for inst in insts:
        for p in write_instance(inst, "V13", "novel_orfs", tmp_path / "plots", ["svg"], 100):
            if p.suffix == ".svg":
                tokens = h.raw_class_tokens(_svg_text(p))
                if tokens:
                    found[p.name] = sorted(set(tokens))
    assert found == {}


def test_panel_a_class_ticks_are_plain_words_in_enum_order(
    tmp_path: Path, cfg: GippaConfig
) -> None:
    rows = [
        h.orf("NOVEL_chr1_1001_p", cls="sorf", start_longest=1001, end=1100),
        h.orf("NOVEL_chr1_2001_p", cls="nterm_extension", start_longest=2001, end=2300),
        h.orf("NOVEL_chr1_4001_p", cls="antisense", start_longest=4001, end=4300),
    ]
    fig = _main(_run(tmp_path, rows, cfg)).figure
    r = _agg(fig)
    ax = _axes(fig, "counts")
    labels = [t.get_text() for t in sorted(ax.get_yticklabels(),
                                           key=lambda t: -t.get_window_extent(r).y0)]  # fmt: skip
    labels = [s for s in labels if s.strip()]
    assert len(labels) == 7, labels
    assert not [s for s in labels if h.raw_class_tokens(s) or "_" in s], labels
    low = [s.lower() for s in labels]
    assert "n-terminal" in low[0] and "c-terminal" in low[1], labels
    assert "intergenic" in low[2] and "frame" in low[3], labels
    assert "antisense" in low[4] and "pseudogene" in low[5], labels


# --- AC 13: determinism and purity ----------------------------------------------------------------


def test_ac13_e1_files_are_byte_identical_across_renders(tmp_path: Path, cfg: GippaConfig) -> None:
    out = vl.build_out(tmp_path / "run", novel_orfs=h.route_e_rows())
    written = []
    for name in ("x", "y"):
        files: dict[Path, bytes] = {}
        for inst in _run_view(RunView.from_outdir(out), cfg):
            for p in write_instance(
                inst, "V13", "novel_orfs", tmp_path / name, ["svg", "png"], 100
            ):
                files[p.relative_to(tmp_path / name)] = p.read_bytes()
        written.append(files)
    assert Path("V13_novel_orfs_locus.svg") in written[0]
    assert Path("data/V13_novel_orfs_locus.tsv") in written[0]
    assert len([p for p in written[0] if p.name.startswith("V13_novel_orfs_locus_")
                and p.suffix == ".svg"]) == 8  # fmt: skip
    assert written[0] == written[1]


def test_ac13_no_pyplot_figure_and_every_text_at_least_6_pt(
    tmp_path: Path, cfg: GippaConfig
) -> None:
    import matplotlib.pyplot as plt

    before = plt.get_fignums()
    insts = _run(tmp_path, h.route_e_rows(), cfg)
    assert plt.get_fignums() == before
    ix = _index(insts)
    for inst in insts:
        assert isinstance(inst.figure, Figure)
        _agg(inst.figure)
        sizes = [t.get_fontsize() for t in _texts(inst.figure)]
        assert sizes, inst.suffix
        assert min(sizes) >= 6, inst.suffix
    assert len(_texts(ix.figure)) >= 8
