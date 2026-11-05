# ruff: noqa: RUF001, RUF002, RUF003  (figure text is Unicode: Δ, ², ±, −, ≤)
"""M12-04 PR 7b: V15 gating and labels (P-07; docs/specs/M12-04-outputs-revision.md § PR 7b,
ADR-0075 item 4, Q6).

Gating rule (m = elastic-net ``cv_r2_mean``, s = ``cv_r2_sd``): the model has **no predictive
power** when m is NA, m ≤ 0, or m − s ≤ 0 (s NA: only m ≤ 0 counts). Then only panels (a) and (d)
are drawn, with the banner ``BANNER``, and the title reads "Elastic net CV R² <m> ± <s>: no
predictive power". Labels: panel titles wrap at the panel width, PD labels use
``words.FEATURE_LABEL``, the missing-layer PD note is plain ("CAI: not computed (codon usage
unavailable)") and sits inside the panel, footers cite no ADR, the (d) legend is outside the
axes, and no text leaves the figure.

V15 presents the M8 statistics and computes none, so there is no calibration test here (§8.4
covers statistical claims; the m ± s display rule is a display rule, Q6). The gating rule gets
Hypothesis properties instead: one on the pure helper (boundary at 0, an interval touching 0,
NA), one rendered on a coarse grid (panels, data TSV, banner and title follow the helper).

Interface assumed beyond the spec text (PR 7b test author):

- ``gippa.viz.determinants.no_predictive_power(model: pd.DataFrame) -> bool``: pure, takes the
  ``determinants_model`` table as ``RunView.table`` returns it and applies the rule to the
  ``elastic_net`` / ``all`` rows ``cv_r2_mean`` and ``cv_r2_sd`` (a missing row counts as NA).
  When ``cv_r2_mean``'s ``ci_low``/``ci_high`` are both finite (never so in the runs we have), an
  interval that includes 0 (``ci_low ≤ 0``) also counts as no predictive power (ADR-0075 item 4,
  "its CI includes 0").
- Axes labels stay as M9-20 set them (``"a"``, ``"b"`` / ``"b_<feature>"``, ``"c"``, ``"d"``);
  with no predictive power there is no ``b``/``b_*``/``c`` axes at all (not an empty one).
- The banner is a ``Text`` artist anywhere in the figure; whitespace (line wrapping) is ignored
  when matching texts.
- A run without abundance (univariate-only ``determinants_model``, no attribution rows, no
  support files) returns an empty state or a ``Skip``; either way no "not available (no
  determinants_pd support file)" text (the empty-state wording itself is PR 4's).
"""

from __future__ import annotations

import math
import re
import tempfile
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from helpers import m12_04_pr5 as p5
from helpers import viz_determinants as vd
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st
from matplotlib.axes import Axes
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from matplotlib.legend import Legend
from matplotlib.text import Text
from matplotlib.transforms import Bbox

from gippa.config import GippaConfig, resolve_config
from gippa.io.plot_files import write_instance
from gippa.io.run_view import RunView
from gippa.viz import base, style, words
from gippa.viz.base import Instance, Skip

BANNER = "The model does not predict abundance in this run; feature effects are not shown."
NO_POWER = "no predictive power"
B_MEAN, B_SD = -0.0174037, 0.0656274  # run B (AC 7b.1)
CAI_NOTE = "CAI: not computed (codon usage unavailable)"
RIBO_LABEL = "ribosomal reference (not in R²)"
FORBIDDEN = re.compile(r"ADR-\d|layer_missing|promotech_")  # AC 7b.3
PD_FEATURES = ("cai", "tai", "tir_ostir", "promoter_score")
_n = style.fmt_num


@pytest.fixture(autouse=True)
def _style() -> None:
    style.apply()


@pytest.fixture
def cfg() -> GippaConfig:
    return resolve_config()


# ---------------------------------------------------------------- helpers


def _v15() -> base.PlotSpec:
    from gippa.steps.plots import load_plot_modules

    load_plot_modules()
    return base.registered()["V15"]


def _result(out: Path, cfg: GippaConfig) -> list[Instance] | Skip:
    return _v15().fn(RunView.from_outdir(out), cfg)


def _render(out: Path, cfg: GippaConfig) -> Instance:
    result = _result(out, cfg)
    assert not isinstance(result, Skip), result
    assert len(result) == 1
    return result[0]


def _with_cv(
    tmp: Path, cfg: GippaConfig, mean: float | None, sd: float | None, **kw: Any
) -> Instance:
    """The main V15 fixture with the elastic net's CV R² mean and SD replaced (None: NA)."""
    return _render(vd.build_out(tmp, model=_model(mean, sd), **kw), cfg)


def _model(mean: float | None, sd: float | None, lo: float | None = None,
           hi: float | None = None) -> list[dict[str, object]]:  # fmt: skip
    rows = vd.model_table()
    for r in rows:
        if r["model"] == "elastic_net" and r["metric"] == "cv_r2_mean":
            r["value"], r["ci_low"], r["ci_high"] = mean, lo, hi
            if lo is not None and hi is not None:
                r["flags"] = None
        elif r["model"] == "elastic_net" and r["metric"] == "cv_r2_sd":
            r["value"] = sd
    return rows


def _flat(text: str) -> str:
    return " ".join(text.split())


def _renderer(fig: Figure) -> Any:
    """An Agg renderer after a full draw (the plot figures carry a base canvas)."""
    canvas = FigureCanvasAgg(fig)
    canvas.draw()
    return canvas.get_renderer()


def _texts(fig: Figure) -> list[str]:
    _renderer(fig)
    return [_flat(t.get_text()) for t in fig.findobj(Text) if t.get_visible() and t.get_text()]


def _labels(inst: Instance) -> list[str]:
    return [ax.get_label() for ax in inst.figure.axes]


def _axes(inst: Instance, label: str) -> Axes:
    found = [ax for ax in inst.figure.axes if ax.get_label() == label]
    assert len(found) == 1, (label, _labels(inst))
    return found[0]


def _hidden_panels(inst: Instance) -> list[str]:
    return [lab for lab in _labels(inst) if lab in ("b", "c") or lab.startswith("b_")]


def _tsv_panels(inst: Instance, tmp: Path) -> set[str]:
    """The ``panel`` values of the written data TSV."""
    files = write_instance(inst, "V15", "determinants", tmp, ["svg"], 50)
    tsv = next(p for p in files if p.suffix == ".tsv")
    back = pd.read_csv(tsv, sep="\t", na_values=["NA"], keep_default_na=False, dtype=str)
    return set(back["panel"])


def _assert_gated(inst: Instance, tmp: Path) -> None:
    """AC 7b.1: exactly panels (a) and (d), no (b)/(c) data rows, the banner, the title."""
    assert _hidden_panels(inst) == [], _labels(inst)
    _axes(inst, "a")
    _axes(inst, "d")
    assert set(inst.data["panel"].astype(str)) == {"a", "d"}
    assert _tsv_panels(inst, tmp) == {"a", "d"}
    assert BANNER in _texts(inst.figure), _texts(inst.figure)
    assert inst.title.startswith("Elastic net CV R² "), inst.title
    assert inst.title.endswith(f": {NO_POWER}"), inst.title


def _assert_full(inst: Instance, tmp: Path) -> None:
    labels = _labels(inst)
    assert {"a", "c", "d"} <= set(labels), labels
    assert {f"b_{f}" for f in PD_FEATURES} <= set(labels), labels
    assert set(inst.data["panel"].astype(str)) == {"a", "b", "c", "d"}
    assert _tsv_panels(inst, tmp) == {"a", "b", "c", "d"}
    assert BANNER not in _texts(inst.figure)
    assert NO_POWER not in inst.title, inst.title


def _model_frame(mean: float | None, sd: float | None, lo: float | None = None,
                 hi: float | None = None) -> pd.DataFrame:  # fmt: skip
    """A ``determinants_model`` frame as ``RunView.table`` returns it (NA as NaN)."""
    frame = pd.DataFrame(_model(mean, sd, lo, hi))
    for c in ("value", "ci_low", "ci_high"):
        frame[c] = pd.to_numeric(frame[c], errors="coerce").astype(float)
    return frame


def _expected(mean: float | None, sd: float | None) -> bool:
    """The spec's rule, written out (§ PR 7b, Q6)."""
    m = math.nan if mean is None else mean
    s = math.nan if sd is None else sd
    if math.isnan(m) or m <= 0:
        return True
    return (not math.isnan(s)) and m - s <= 0


# ---------------------------------------------------------------- gating: the pure helper


def test_helper_examples_from_the_spec() -> None:
    from gippa.viz.determinants import no_predictive_power as npp

    assert npp(_model_frame(B_MEAN, B_SD)) is True  # B
    assert npp(_model_frame(0.05, 0.06)) is True  # m − s < 0 (AC 7b.2)
    assert npp(_model_frame(0.133, 0.035)) is False  # K
    assert npp(_model_frame(0.0359, 0.0043)) is False  # the snapshot fixture value
    assert npp(_model_frame(0.21, 0.03)) is False  # the M9-20 main fixture


def test_helper_boundaries_at_zero() -> None:
    from gippa.viz.determinants import no_predictive_power as npp

    assert npp(_model_frame(0.0, 0.0)) is True  # m = 0
    assert npp(_model_frame(0.0, None)) is True  # m = 0, s NA
    assert npp(_model_frame(0.25, 0.25)) is True  # m − s = 0: the interval touches 0
    assert npp(_model_frame(0.25, 0.24)) is False
    assert npp(_model_frame(None, 0.03)) is True  # m NA
    assert npp(_model_frame(None, None)) is True
    assert npp(_model_frame(0.05, None)) is False  # edge case: s NA with m > 0 is drawn
    assert npp(_model_frame(-0.05, None)) is True


def test_helper_missing_rows_count_as_na() -> None:
    from gippa.viz.determinants import no_predictive_power as npp

    frame = _model_frame(0.1, 0.02)
    no_sd = frame[~((frame["model"] == "elastic_net") & (frame["metric"] == "cv_r2_sd"))]
    assert npp(no_sd.reset_index(drop=True)) is False  # s missing, m > 0
    no_mean = frame[~((frame["model"] == "elastic_net") & (frame["metric"] == "cv_r2_mean"))]
    assert npp(no_mean.reset_index(drop=True)) is True  # m missing


def test_helper_reads_the_elastic_net_not_lightgbm() -> None:
    """ADR-0033: the elastic net is primary; a good LightGBM does not lift the gate."""
    from gippa.viz.determinants import no_predictive_power as npp

    rows = vd.model_table(en_cv=(B_MEAN, B_SD, 0.4), lgb_cv=(0.5, 0.01, 0.01))
    frame = pd.DataFrame(rows)
    frame["value"] = pd.to_numeric(frame["value"], errors="coerce").astype(float)
    assert npp(frame) is True
    rows = vd.model_table(en_cv=(0.2, 0.01, 0.01), lgb_cv=(-0.3, 0.1, 0.9))
    frame = pd.DataFrame(rows)
    frame["value"] = pd.to_numeric(frame["value"], errors="coerce").astype(float)
    assert npp(frame) is False


def test_helper_an_estimated_ci_that_includes_zero_hides() -> None:
    """ADR-0075 item 4 ("or its CI includes 0"): used when the interval is estimated."""
    from gippa.viz.determinants import no_predictive_power as npp

    assert npp(_model_frame(0.1, 0.02, -0.01, 0.2)) is True
    assert npp(_model_frame(0.1, 0.02, 0.0, 0.2)) is True  # touches 0
    assert npp(_model_frame(0.1, 0.02, 0.05, 0.15)) is False


_MEAN = st.one_of(
    st.none(),
    st.just(0.0),
    st.just(-0.0),
    st.floats(-1.0, 1.0, allow_nan=False),
    st.floats(-1e-9, 1e-9, allow_nan=False),
)
_SD = st.one_of(st.none(), st.just(0.0), st.floats(0.0, 1.0, allow_nan=False))


@settings(max_examples=300, deadline=None)
@given(_MEAN, _SD)
def test_property_helper_matches_the_rule(mean: float | None, sd: float | None) -> None:
    """Any (m, s), with the boundary at 0 and NA drawn often; the CV interval unestimated."""
    from gippa.viz.determinants import no_predictive_power as npp

    assert npp(_model_frame(mean, sd)) is _expected(mean, sd), (mean, sd)


@settings(max_examples=200, deadline=None)
@given(st.floats(1e-6, 1.0), st.floats(0.0, 1.0))
def test_property_interval_touching_zero_hides(mean: float, frac: float) -> None:
    """s = m (the interval m ± s touches 0) hides; any s < m with m > 0 shows."""
    from gippa.viz.determinants import no_predictive_power as npp

    assert npp(_model_frame(mean, mean)) is True
    sd = mean * frac
    if sd < mean:
        assert npp(_model_frame(mean, sd)) is False, (mean, sd)


# ---------------------------------------------------------------- gating: the figure


def test_ac_7b_1_run_b_draws_only_a_and_d(tmp_path: Path, cfg: GippaConfig) -> None:
    inst = _with_cv(tmp_path, cfg, B_MEAN, B_SD)
    _assert_gated(inst, tmp_path / "plots")
    assert inst.title == f"Elastic net CV R² {_n(B_MEAN)} ± {_n(B_SD)}: {NO_POWER}", inst.title
    assert inst.title == "Elastic net CV R² −0.017 ± 0.066: no predictive power"


def test_ac_7b_1_no_shap_or_pd_artists(tmp_path: Path, cfg: GippaConfig) -> None:
    """The support files exist, but nothing from them is drawn: no SHAP colour bar either."""
    inst = _with_cv(tmp_path, cfg, B_MEAN, B_SD)
    text = "\n".join(_texts(inst.figure))
    for gone in ("SHAP", "partial dependence", "feature value rank", "not available"):
        assert gone not in text, (gone, text)


def test_ac_7b_2_interval_below_zero_is_gated(tmp_path: Path, cfg: GippaConfig) -> None:
    _assert_gated(_with_cv(tmp_path, cfg, 0.05, 0.06), tmp_path / "plots")


def test_ac_7b_2_snapshot_value_draws_all_four(tmp_path: Path, cfg: GippaConfig) -> None:
    _assert_full(_with_cv(tmp_path, cfg, 0.0359, 0.0043), tmp_path / "plots")


def test_example_k_draws_all_four(tmp_path: Path, cfg: GippaConfig) -> None:
    _assert_full(_with_cv(tmp_path, cfg, 0.133, 0.035), tmp_path / "plots")


def test_snapshot_fixture_draws_all_four(tmp_path: Path) -> None:
    """The committed ``sim_plots`` run (elastic net 0.132854 ± 0.0351372): gating unchanged."""
    inst = p5.only(p5.render(p5.SIM_OUT, ["V15"]), "V15")
    _assert_full(inst, tmp_path / "plots")


def test_edge_sd_na_with_positive_mean_is_drawn(tmp_path: Path, cfg: GippaConfig) -> None:
    _assert_full(_with_cv(tmp_path, cfg, 0.05, None), tmp_path / "plots")


def test_edge_mean_na_is_gated(tmp_path: Path, cfg: GippaConfig) -> None:
    _assert_gated(_with_cv(tmp_path, cfg, None, 0.03), tmp_path / "plots")


def test_edge_interval_touching_zero_is_gated(tmp_path: Path, cfg: GippaConfig) -> None:
    _assert_gated(_with_cv(tmp_path, cfg, 0.25, 0.25), tmp_path / "plots")


def test_edge_mean_exactly_zero_is_gated(tmp_path: Path, cfg: GippaConfig) -> None:
    _assert_gated(_with_cv(tmp_path, cfg, 0.0, 0.0), tmp_path / "plots")


#: (m, s) in hundredths: integers, so m − s is exact at the boundary after the %.6g round trip
_GRID = st.tuples(st.integers(-10, 10), st.one_of(st.none(), st.integers(0, 10)))


@settings(
    max_examples=12,
    deadline=None,
    suppress_health_check=[HealthCheck.function_scoped_fixture, HealthCheck.too_slow],
)
@given(_GRID)
def test_property_rendered_panels_follow_the_rule(ms: tuple[int, int | None]) -> None:
    m100, s100 = ms
    mean, sd = m100 / 100, (None if s100 is None else s100 / 100)
    cfg = resolve_config()
    with tempfile.TemporaryDirectory() as tmp:
        inst = _render(vd.build_out(Path(tmp), model=_model(mean, sd)), cfg)
        try:
            if _expected(mean, sd):
                _assert_gated(inst, Path(tmp) / "plots")
            else:
                _assert_full(inst, Path(tmp) / "plots")
        finally:
            inst.figure.clear()


def test_edge_no_abundance_has_no_not_available_panels(tmp_path: Path, cfg: GippaConfig) -> None:
    """Run C (no abundance): univariate-only rows, no predictions, no support files."""
    out = vd.build_out(tmp_path, model=vd.univariate_rows(), attribution=[], support=False)
    result = _result(out, cfg)
    if isinstance(result, Skip):
        return
    assert len(result) == 1
    inst = result[0]
    assert inst.empty, _texts(inst.figure)
    text = "\n".join(_texts(inst.figure))
    assert "not available" not in text, text
    assert "support file" not in text, text


# ---------------------------------------------------------------- labels


def _pd_xlabel(inst: Instance, feature: str) -> str:
    return _flat(_axes(inst, f"b_{feature}").get_xlabel())


def test_pd_labels_use_feature_label(tmp_path: Path, cfg: GippaConfig) -> None:
    inst = _render(vd.build_out(tmp_path), cfg)
    for f in PD_FEATURES:
        lab = _pd_xlabel(inst, f)
        assert words.FEATURE_LABEL[f] in lab, (f, lab)
        assert not re.search(rf"\b{re.escape(f)}\b", lab), (f, lab)  # no column name
        assert "_" not in lab, (f, lab)


def _layer_missing(tmp: Path, cfg: GippaConfig) -> Instance:
    pdf = vd.pd_frame(("tir_ostir", "promoter_score"))
    out = vd.build_out(tmp, pd_=pdf, flags=["layer_missing:codon"])
    return _render(out, cfg)


def test_cai_not_computed_note_is_plain(tmp_path: Path, cfg: GippaConfig) -> None:
    inst = _layer_missing(tmp_path, cfg)
    ax = _axes(inst, "b_cai")
    inside = [_flat(t.get_text()) for t in ax.texts if t.get_visible()]
    assert CAI_NOTE in inside, inside
    tai = " ".join(_flat(t.get_text()) for t in _axes(inst, "b_tai").texts)
    assert "not computed" in tai, tai
    every = "\n".join(_texts(inst.figure))
    assert "layer_missing" not in every, every
    assert "not available (layer_missing" not in every, every


def test_cai_note_is_inside_the_panel_not_on_a_tick(tmp_path: Path, cfg: GippaConfig) -> None:
    inst = _layer_missing(tmp_path, cfg)
    ax = _axes(inst, "b_cai")
    renderer = _renderer(inst.figure)
    note = [t for t in ax.texts if CAI_NOTE in _flat(t.get_text())]
    assert len(note) == 1, [t.get_text() for t in ax.texts]
    box, bb = ax.get_window_extent(renderer), note[0].get_window_extent(renderer)
    assert box.x0 - 1 <= bb.x0 and bb.x1 <= box.x1 + 1, (box, bb)
    assert box.y0 - 1 <= bb.y0 and bb.y1 <= box.y1 + 1, (box, bb)
    ticks = [_flat(t.get_text()) for t in ax.get_xticklabels() + ax.get_yticklabels()]
    assert not any("not computed" in t or "CAI" in t for t in ticks), ticks


def _all_fixtures(tmp: Path, cfg: GippaConfig) -> list[Instance]:
    return [
        _render(vd.build_out(tmp / "main"), cfg),
        _with_cv(tmp / "b", cfg, B_MEAN, B_SD),
        _layer_missing(tmp / "codon", cfg),
        _render(
            vd.build_out(
                tmp / "ostir",
                pd_=vd.pd_frame(("cai", "tai", "promoter_score")),
                flags=["layer_missing:ostir"],
            ),
            cfg,
        ),
        _render(vd.build_out(tmp / "nosup", support=False), cfg),
    ]


def test_ac_7b_3_no_internal_references(tmp_path: Path, cfg: GippaConfig) -> None:
    for inst in _all_fixtures(tmp_path, cfg):
        for t in [inst.title, *_texts(inst.figure)]:
            assert not FORBIDDEN.search(t), t


def test_ac_7b_3_no_internal_references_in_the_twin(tmp_path: Path, cfg: GippaConfig) -> None:
    import json

    from gippa.viz.interactive import figure_json

    inst = _render(vd.build_out(tmp_path), cfg)
    text = figure_json("V15", inst.data, inst.title, None)
    assert text is not None
    fig = json.loads(text)
    shown = [str(a.get("text", "")) for a in fig["layout"].get("annotations", [])]
    shown += [str(t.get("name", "")) for t in fig["data"]]
    for s in shown:
        assert not FORBIDDEN.search(s), s


def test_footers_are_plain_caveats(tmp_path: Path, cfg: GippaConfig) -> None:
    """The ΔR² footer stays a caveat (correlated blocks; ΔR² = the loss on removing a block)."""
    text = "\n".join(_texts(_render(vd.build_out(tmp_path), cfg).figure))
    assert "ADR" not in text, text
    assert re.search(r"block is removed|removing a block", text), text
    assert "correlated" in text, text


def test_d_legend_is_outside_the_axes(tmp_path: Path, cfg: GippaConfig) -> None:
    """Example K: all four panels; the (d) legend names the ribosomal reference, outside."""
    inst = _with_cv(tmp_path, cfg, 0.133, 0.035)
    fig = inst.figure
    renderer = _renderer(fig)
    legends = [lg for lg in fig.findobj(Legend)
               if RIBO_LABEL in [t.get_text() for t in lg.get_texts()]]  # fmt: skip
    assert len(legends) == 1, [[t.get_text() for t in lg.get_texts()] for lg in
                               fig.findobj(Legend)]  # fmt: skip
    box = _axes(inst, "d").get_window_extent(renderer)
    leg = legends[0].get_window_extent(renderer)
    inter = Bbox.intersection(box, leg)
    area = 0.0 if inter is None else max(inter.width, 0.0) * max(inter.height, 0.0)
    assert area < 1.0, (box, leg)


def _title_texts(ax: Axes) -> list[Text]:
    return [t for t in (ax.title, ax._left_title, ax._right_title) if t.get_text()]


def _union(boxes: list[Bbox]) -> Bbox:
    return Bbox.union(boxes)


def test_panel_titles_wrap_at_the_panel_width(tmp_path: Path, cfg: GippaConfig) -> None:
    for inst in (_render(vd.build_out(tmp_path / "m"), cfg),
                 _with_cv(tmp_path / "b", cfg, B_MEAN, B_SD)):  # fmt: skip
        fig = inst.figure
        renderer = _renderer(fig)
        b_axes = [ax for ax in fig.axes if ax.get_label().startswith("b")]
        b_box = _union([ax.get_window_extent(renderer) for ax in b_axes]) if b_axes else None
        for ax in fig.axes:
            lab = ax.get_label()
            if lab not in ("a", "c", "d") and not lab.startswith("b"):
                continue
            width = (b_box if lab.startswith("b") and b_box is not None
                     else ax.get_window_extent(renderer)).width  # fmt: skip
            for t in _title_texts(ax):
                w = t.get_window_extent(renderer).width
                assert w <= width + 2.0, (lab, t.get_text(), w, width)


def _outside(fig: Figure) -> list[tuple[str, Any]]:
    """Visible texts drawn outside the figure. Tick labels count only for the ticks drawn
    (inside the view interval); matplotlib keeps stale tick objects beyond it."""
    renderer = _renderer(fig)
    fb = fig.bbox
    ticks, drawn = set(), set()
    for ax in fig.axes:
        for axis in (ax.xaxis, ax.yaxis):
            for tick in [*axis.majorTicks, *axis.minorTicks]:
                ticks |= {id(tick.label1), id(tick.label2)}
            for tick in axis._update_ticks():
                drawn |= {id(tick.label1), id(tick.label2)}
    bad = []
    for t in fig.findobj(Text):
        if not (t.get_visible() and t.get_text()):
            continue
        if id(t) in ticks and id(t) not in drawn:
            continue
        bb = t.get_window_extent(renderer)
        if not np.all(np.isfinite([bb.x0, bb.x1, bb.y0, bb.y1])):
            continue
        if bb.x0 < fb.x0 - 1 or bb.x1 > fb.x1 + 1 or bb.y0 < fb.y0 - 1 or bb.y1 > fb.y1 + 1:
            bad.append((t.get_text(), bb.bounds))
    return bad


def test_ac_7b_4_no_text_leaves_the_figure(tmp_path: Path, cfg: GippaConfig) -> None:
    for inst in _all_fixtures(tmp_path, cfg):
        assert _outside(inst.figure) == [], inst.title


def test_ac_7b_4_snapshot_fixture_text_inside(tmp_path: Path) -> None:
    inst = p5.only(p5.render(p5.SIM_OUT, ["V15"]), "V15")
    assert p5.sf.styled(lambda: _outside(inst.figure)) == []
