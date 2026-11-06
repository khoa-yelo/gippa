# ruff: noqa: RUF002, RUF003  (figure text is Unicode: Δ, ², ±, −, ≤)
"""M12-04 PR 7b review round 1 (#562): V15's twin follows the gate, run C is an empty state,
the (d) legend clears its x label and the footers, and the gate agrees with the report card.

Items (lead, review round 1):

1. The Plotly twin of a gated V15 (run B) has only the (a) and (d) panels, no "(b)"/"(c)"
   subplot title or annotation, and the banner as an annotation; an ungated run keeps all four.
2. Run C (univariate rows only, no abundance, run flag ``layer_missing:abundance``) is an empty
   state, not a skip: an ``Instance`` with ``empty=True`` whose title starts "Nothing to show: ",
   says the run has no abundance, ends with ".", holds no raw token and never says "no
   determinants model". (The PR 7b test ``test_edge_no_abundance_has_no_not_available_panels`` is
   tightened the same way.)
3. On the K and B fixtures, the (d) legend's window extent is at least 4 px (at 100 dpi) clear of
   the (d) x label and overlaps no footer text (the figure-level texts other than the banner).
4. For a model whose ``cv_r2_mean`` has an estimated interval with ``ci_low ≤ 0`` while m − s > 0,
   V15 says "no predictive power" and the report's determinants card treats the model as not
   predictive. The card and the plot use the one helper,
   ``gippa.viz.determinants.no_predictive_power``; the test pins agreement (the card shows
   ``gippa.report.build.NO_POWER`` exactly when the helper says so), not the card's wording.
"""

from __future__ import annotations

import json
import re
import tempfile
from pathlib import Path
from typing import Any

import pytest
from helpers import viz_determinants as vd
from helpers.report_content import det_rows, write_det_model
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from matplotlib.legend import Legend
from matplotlib.text import Text
from matplotlib.transforms import Bbox
from unit.report.test_m12_04_pr1_fixtures import Run, fixture_d, fixture_r, text

from gippa.config import GippaConfig, resolve_config
from gippa.io.plot_files import write_instance
from gippa.io.run_view import RunView
from gippa.viz import base, style
from gippa.viz.base import Instance, Skip

BANNER = "The model does not predict abundance in this run; feature effects are not shown."
NO_POWER = "no predictive power"
B_MEAN, B_SD = -0.0174037, 0.0656274
K_MEAN, K_SD = 0.133, 0.035
RIBO_LABEL = "ribosomal reference (not in R²)"
RAW_TOKEN = re.compile(r"\w+_\w+|\w+:\w+|ADR-\d|§\s?\d")
GAP_PX = 4.0


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


def _with_cv(tmp: Path, cfg: GippaConfig, mean: float, sd: float, **kw: Any) -> Instance:
    return _render(vd.build_out(tmp, model=_model(mean, sd, **kw)), cfg)


def _flat(s: str) -> str:
    return " ".join(s.split())


def _renderer(fig: Figure, dpi: float | None = None) -> Any:
    if dpi is not None:
        fig.set_dpi(dpi)
    canvas = FigureCanvasAgg(fig)
    canvas.draw()
    return canvas.get_renderer()


def _gap(a: Bbox, b: Bbox) -> float:
    """The clearance between two boxes in pixels (0 when they touch or overlap)."""
    dx = max(a.x0, b.x0) - min(a.x1, b.x1)
    dy = max(a.y0, b.y0) - min(a.y1, b.y1)
    return max(dx, dy, 0.0)


def _overlap(a: Bbox, b: Bbox) -> float:
    w = min(a.x1, b.x1) - max(a.x0, b.x0)
    h = min(a.y1, b.y1) - max(a.y0, b.y0)
    return max(w, 0.0) * max(h, 0.0)


# ---------------------------------------------------------------- item 1: the twin


def _twin(inst: Instance, tmp: Path) -> dict[str, Any]:
    """The twin built as the plots step builds it: from the written data TSV."""
    import pandas as pd

    from gippa.steps.plots import TWIN_ID_DTYPES
    from gippa.viz.interactive import figure_json

    write_instance(inst, "V15", "determinants", tmp, ["svg"], 50)
    data = pd.read_csv(tmp / "data" / "V15_determinants.tsv", sep="\t", na_values=["NA"],
                       keep_default_na=False, dtype=TWIN_ID_DTYPES)  # fmt: skip
    blob = figure_json("V15", data, inst.title, None)
    assert blob is not None
    fig: dict[str, Any] = json.loads(blob)
    return fig


def _annotations(fig: dict[str, Any]) -> list[str]:
    return [_flat(re.sub(r"<br\s*/?>", " ", str(a.get("text", ""))))
            for a in fig["layout"].get("annotations", [])]  # fmt: skip


def _xaxes(fig: dict[str, Any]) -> list[str]:
    return sorted(k for k in fig["layout"] if re.fullmatch(r"xaxis\d*", k))


def test_item1_gated_twin_has_only_a_and_d(tmp_path: Path, cfg: GippaConfig) -> None:
    inst = _with_cv(tmp_path / "b", cfg, B_MEAN, B_SD)
    fig = _twin(inst, tmp_path / "plots")
    notes = _annotations(fig)
    assert not [n for n in notes if n.startswith(("(b)", "(c)"))], notes
    assert not [n for n in notes if "(b)" in n or "(c)" in n], notes
    assert any(n.startswith("(a)") for n in notes), notes
    assert any(n.startswith("(d)") for n in notes), notes
    assert len(_xaxes(fig)) == 2, _xaxes(fig)


def test_item1_gated_twin_has_the_banner(tmp_path: Path, cfg: GippaConfig) -> None:
    inst = _with_cv(tmp_path / "b", cfg, B_MEAN, B_SD)
    assert BANNER in _annotations(_twin(inst, tmp_path / "plots"))


def test_item1_ungated_twin_keeps_four_panels(tmp_path: Path, cfg: GippaConfig) -> None:
    inst = _with_cv(tmp_path / "k", cfg, K_MEAN, K_SD)
    fig = _twin(inst, tmp_path / "plots")
    notes = _annotations(fig)
    for head in ("(a)", "(b)", "(c)", "(d)"):
        assert any(n.startswith(head) for n in notes), (head, notes)
    assert BANNER not in notes, notes
    assert len(_xaxes(fig)) > 2, _xaxes(fig)


# ---------------------------------------------------------------- item 2: run C


def _assert_c_empty_state(result: list[Instance] | Skip) -> None:
    assert not isinstance(result, Skip), result
    assert len(result) == 1
    inst = result[0]
    assert inst.empty, inst.title
    t = inst.title
    assert t.startswith("Nothing to show: "), t
    assert t.endswith("."), t
    assert "abundance" in t.lower(), t
    assert re.search(r"\bno\b.*\babundance\b", t.lower()), t
    assert "no determinants model" not in t.lower(), t
    assert not RAW_TOKEN.search(t), t
    shown = [_flat(x.get_text()) for x in inst.figure.findobj(Text) if x.get_text()]
    assert not any("not available" in s or "support file" in s for s in shown), shown


def test_item2_run_c_vd_fixture_is_an_empty_state(tmp_path: Path, cfg: GippaConfig) -> None:
    """Univariate rows only, no attribution rows, no support files."""
    out = vd.build_out(tmp_path, model=vd.univariate_rows(), attribution=[], support=False)
    _assert_c_empty_state(_result(out, cfg))


def test_item2_run_c_with_the_abundance_run_flag_is_an_empty_state(
    tmp_path: Path, cfg: GippaConfig
) -> None:
    """The PR 1 fixture D (run C): univariate AUROC rows and ``layer_missing:abundance``."""
    run = fixture_d(tmp_path)
    _assert_c_empty_state(_result(run.out, cfg))


# ---------------------------------------------------------------- item 3: the (d) legend


def _d_legend(fig: Figure) -> Legend:
    found = [lg for lg in fig.findobj(Legend)
             if RIBO_LABEL in [t.get_text() for t in lg.get_texts()]]  # fmt: skip
    assert len(found) == 1, [[t.get_text() for t in lg.get_texts()] for lg in fig.findobj(Legend)]
    return found[0]


def _footers(fig: Figure) -> list[Text]:
    return [t for t in fig.texts if t.get_visible() and t.get_text()
            and _flat(t.get_text()) != BANNER]  # fmt: skip


@pytest.mark.parametrize(("name", "mean", "sd"), [("K", K_MEAN, K_SD), ("B", B_MEAN, B_SD)])
def test_item3_d_legend_clears_the_x_label(
    tmp_path: Path, cfg: GippaConfig, name: str, mean: float, sd: float
) -> None:
    inst = _with_cv(tmp_path, cfg, mean, sd)
    fig = inst.figure
    renderer = _renderer(fig, dpi=100)
    (ax,) = [a for a in fig.axes if a.get_label() == "d"]
    leg = _d_legend(fig).get_window_extent(renderer)
    xlab = ax.xaxis.label.get_window_extent(renderer)
    assert ax.get_xlabel(), name
    assert _overlap(leg, xlab) == 0.0, (name, leg, xlab)
    assert _gap(leg, xlab) >= GAP_PX, (name, _gap(leg, xlab), leg, xlab)


@pytest.mark.parametrize(("name", "mean", "sd"), [("K", K_MEAN, K_SD), ("B", B_MEAN, B_SD)])
def test_item3_d_legend_overlaps_no_footer(
    tmp_path: Path, cfg: GippaConfig, name: str, mean: float, sd: float
) -> None:
    inst = _with_cv(tmp_path, cfg, mean, sd)
    fig = inst.figure
    renderer = _renderer(fig, dpi=100)
    leg = _d_legend(fig).get_window_extent(renderer)
    feet = _footers(fig)
    assert feet, name
    for t in feet:
        bb = t.get_window_extent(renderer)
        assert _overlap(leg, bb) == 0.0, (name, t.get_text()[:60], leg, bb)


# ---------------------------------------------------------------- item 4: plot and card agree


def _card(run: Run) -> str:
    cards = {k: text(v) for k, v in run.page().cards().items()}
    assert "model_r2" in cards, sorted(cards)
    return cards["model_r2"]


def _card_rows(mean: float, sd: float, lo: float | None, hi: float | None) -> list[dict]:
    rows = det_rows("elastic_net", mean, sd, 0.3, "null_uncalibrated")
    for r in rows:
        if r["metric"] == "cv_r2_mean":
            r["ci_low"], r["ci_high"] = lo, hi
    return rows


def _helper(out: Path) -> bool:
    from gippa.viz.determinants import no_predictive_power

    model = RunView.from_outdir(out).table("determinants_model")
    assert model is not None
    return no_predictive_power(model)


def test_item4_estimated_ci_including_zero_plot_says_no_power(
    tmp_path: Path, cfg: GippaConfig
) -> None:
    inst = _with_cv(tmp_path, cfg, 0.1, 0.02, lo=-0.01, hi=0.2)
    assert NO_POWER in inst.title, inst.title


def test_item4_estimated_ci_including_zero_card_says_no_power(tmp_path: Path) -> None:
    from gippa.report.build import NO_POWER as CARD_NO_POWER

    run = fixture_r(tmp_path)
    write_det_model(run.out, _card_rows(0.1, 0.02, -0.01, 0.2))
    assert _helper(run.out) is True
    card = _card(run)
    assert CARD_NO_POWER in card, card


@settings(
    max_examples=15,
    deadline=None,
    suppress_health_check=[HealthCheck.function_scoped_fixture, HealthCheck.too_slow],
)
@given(
    st.integers(-20, 30),
    st.integers(0, 20),
    st.one_of(st.none(), st.tuples(st.integers(-20, 20), st.integers(0, 20))),
)
def test_property_item4_card_follows_the_helper(
    m100: int, s100: int, ci: tuple[int, int] | None
) -> None:
    """In hundredths (exact after the %.6g round trip); the CI is (lo, lo + width) or NA."""
    from gippa.report.build import NO_POWER as CARD_NO_POWER

    mean, sd = m100 / 100, s100 / 100
    lo, hi = (None, None) if ci is None else (ci[0] / 100, (ci[0] + ci[1]) / 100)
    with tempfile.TemporaryDirectory() as tmp:
        run = fixture_r(Path(tmp))
        write_det_model(run.out, _card_rows(mean, sd, lo, hi))
        want = _helper(run.out)
        card = _card(run)
    assert (CARD_NO_POWER in card) is want, (mean, sd, lo, hi, card)
