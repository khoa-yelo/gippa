# ruff: noqa: RUF001, RUF002, RUF003  (figure text is Unicode: −, ≤)
"""M12-04 PR 6, V07 layout (P-08; docs/specs/M12-04-outputs-revision.md § PR 6, AC 6.1, 6.2,
6.4).

Interface assumed where the spec leaves it open (test author; reported to the lead):

- The figure is always double width (``style.DOUBLE_IN``). The forest is the one axes whose x
  label is ``XLABEL_PP`` (the spec's "full label").
- The data TSV keeps the M9 columns and adds ``o_minus_e_pp``, ``ci_low_pp`` and ``ci_high_pp``
  = 100 · value / ``n``. Points and CI ends are drawn at the ``_pp`` values.
- Group headers are bold text artists: "Replicon", "Feature class", "COG category" for the basic
  blocks; the region block's header is "Region" (or "Region type").
- A COG row label starts with "<letter> " and a word of the letter's NCBI description.
- With regions, the no-hit title reads "No group departs from detectability at q ≤ <alpha> (…
  tested)" and names the region count with a correct plural.
- Above the height cap, the forest draws the rows with the largest |pp| and one text reads
  "<k> more rows in the data TSV", k = data rows − drawn rows; the data TSV keeps every row.
"""

from __future__ import annotations

import math
import re
import tempfile
from pathlib import Path

import pytest
from helpers import m12_04_pr5 as h5
from helpers import m12_04_pr6 as h
from helpers.viz_artists import markers
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st
from matplotlib.axes import Axes

from gippa.config import GippaConfig, resolve_config
from gippa.io.cog_categories import cog_descriptions
from gippa.io.run_view import RunView
from gippa.viz import base, style
from gippa.viz.base import Instance, Skip

XLABEL_PP = "observed − expected detected (percentage points)"
PP_COLUMNS = ("o_minus_e_pp", "ci_low_pp", "ci_high_pp")
HEADERS = {"replicon": {"Replicon"}, "feature_class": {"Feature class"},
           "cog_category": {"COG category"}, "region": {"Region", "Region type"}}  # fmt: skip
MORE_ROWS = re.compile(r"^(\d[\d,]*) more rows in the data TSV$")


@pytest.fixture(autouse=True)
def _style() -> None:
    style.apply()


@pytest.fixture
def cfg() -> GippaConfig:
    return resolve_config()


def _render(out: Path, cfg: GippaConfig) -> Instance:
    import gippa.viz.detection  # noqa: F401  (registers V07)

    result = base.registered()["V07"].fn(RunView.from_outdir(out), cfg)
    assert not isinstance(result, Skip), result
    assert len(result) == 1
    return result[0]


def _forest(inst: Instance) -> Axes:
    found = [ax for ax in inst.figure.axes if ax.get_xlabel() == XLABEL_PP]
    assert len(found) == 1, [ax.get_xlabel() for ax in inst.figure.axes]
    return found[0]


def _g6(x: float) -> float:
    return float(f"{x:.6g}")


# --- AC 6.1: width and the pp columns -----------------------------------------------------------


@pytest.mark.parametrize("regions", ["few", "none"])
def test_figure_is_double_width(tmp_path: Path, cfg: GippaConfig, regions: str) -> None:
    """A small forest (6 rows) is double width too (example C: overlapping ticks)."""
    out = h.v07_out(tmp_path, regions=h.A_REGIONS if regions == "few" else None)
    fig = _render(out, cfg).figure
    assert fig.get_figwidth() == pytest.approx(style.DOUBLE_IN)


def test_data_has_the_pp_columns_and_their_values(tmp_path: Path, cfg: GippaConfig) -> None:
    data = _render(h.v07_out(tmp_path), cfg).data
    assert set(PP_COLUMNS) <= set(data.columns), list(data.columns)
    for r in data.to_dict("records"):
        n = h.num(r["n"])
        for raw, pp in zip(("o_minus_e", "ci_low", "ci_high"), PP_COLUMNS, strict=True):
            if h.isna(r[raw]):
                assert h.isna(r[pp]), r
            else:
                assert h.num(r[pp]) == pytest.approx(_g6(100 * h.num(r[raw]) / n), rel=1e-5), r


def test_example_replicon_row_in_pp(tmp_path: Path, cfg: GippaConfig) -> None:
    """chr1: O − E 2.1 of 40 genes = 5.25 pp, CI −1.0..5.2 -> −2.5..13 pp."""
    data = _render(h.v07_out(tmp_path), cfg).data
    row = data[data["group"].astype(str) == "chr1"].iloc[0]
    assert h.num(row["o_minus_e_pp"]) == pytest.approx(5.25)
    assert h.num(row["ci_low_pp"]) == pytest.approx(-2.5)
    assert h.num(row["ci_high_pp"]) == pytest.approx(13.0)


def test_x_axis_is_percentage_points(tmp_path: Path, cfg: GippaConfig) -> None:
    inst = _render(h.v07_out(tmp_path), cfg)
    ax = _forest(inst)
    xs = sorted(round(m.x, 6) for m in markers(ax))
    want = sorted(round(h.num(v), 6) for v in inst.data["o_minus_e_pp"] if not h.isna(v))
    assert xs == want


def test_ci_ends_are_drawn_in_pp(tmp_path: Path, cfg: GippaConfig) -> None:
    inst = _render(h.v07_out(tmp_path), cfg)
    ax = _forest(inst)
    xs: set[float] = set()
    for line in ax.get_lines():
        xs.update(round(float(x), 6) for x in line.get_xdata())
    for coll in ax.collections:
        for seg in getattr(coll, "get_segments", list)():
            xs.update(round(float(x), 6) for x, _ in seg)
    for lo, hi in zip(inst.data["ci_low_pp"], inst.data["ci_high_pp"], strict=True):
        if not h.isna(lo):
            assert round(h.num(lo), 6) in xs and round(h.num(hi), 6) in xs, (lo, hi)


_ROWS = st.lists(
    st.tuples(st.integers(1, 5000), st.floats(-0.9, 0.9), st.floats(0.0, 0.5)),
    min_size=1, max_size=8,
)  # fmt: skip


@settings(max_examples=20, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(rows=_ROWS)
def test_property_pp_is_100_ome_over_n(rows: list[tuple[int, float, float]]) -> None:
    """AC 6.1 for any n (1 to K-like 5,000) and any O − E inside ±n."""
    style.apply()
    regions = []
    for i, (n, frac, half) in enumerate(rows):
        ome = round(frac * n, 4)
        regions.append(h.region(f"prophage_{i}", "prophage", n_eligible=n, observed=1.0,
                                expected=1.0 - ome, effect=ome,
                                ci=(ome - half * n, ome + half * n), q=0.5))  # fmt: skip
    with tempfile.TemporaryDirectory() as tmp:
        data = _render(h.v07_out(Path(tmp), regions=regions), resolve_config()).data
    for r in data.to_dict("records"):
        n = h.num(r["n"])
        for raw, pp in zip(("o_minus_e", "ci_low", "ci_high"), PP_COLUMNS, strict=True):
            assert h.num(r[pp]) == pytest.approx(_g6(100 * h.num(r[raw]) / n), rel=1e-5)


# --- AC 6.2: no n = 0 rows; group headers ---------------------------------------------------------


def test_rows_with_n_zero_are_dropped(tmp_path: Path, cfg: GippaConfig) -> None:
    """K: "crispr_2 (n = 0)" was drawn."""
    inst = _render(h.v07_out(tmp_path), cfg)
    assert not any(h.num(n) == 0 for n in inst.data["n"]), inst.data
    assert "crispr_2" not in set(inst.data["group"].astype(str))
    assert not any("crispr_2" in t for t in h.texts(inst.figure))


def test_sim_plots_has_no_n_zero_row(tmp_path: Path) -> None:
    """The sim_plots run has ``crispr_1`` with ``n_eligible = 0``."""
    inst = h5.only(h5.render(h5.sim_copy(tmp_path), ["V07"]), "V07")
    assert not any(h.num(n) == 0 for n in inst.data["n"])


def _headers(inst: Instance) -> dict[str, list[bool]]:
    found: dict[str, list[bool]] = {}
    for t in h.drawn_texts(inst.figure):
        found.setdefault(t.get_text(), []).append(h.is_bold(t))
    return found


def test_every_group_present_has_a_bold_header(tmp_path: Path, cfg: GippaConfig) -> None:
    inst = _render(h.v07_out(tmp_path), cfg)
    found = _headers(inst)
    for block in sorted(set(inst.data["block"].astype(str))):
        names = HEADERS[block] & set(found)
        assert names, f"no header for {block}: {sorted(found)}"
        assert any(any(found[n]) for n in names), f"header of {block} is not bold"


def test_no_cog_group_header_without_function(tmp_path: Path, cfg: GippaConfig) -> None:
    """Edge case: no function (B, C) -> no COG group header."""
    out = h.v07_out(tmp_path, flags=["layer_missing:cog"])
    inst = _render(out, cfg)
    found = _headers(inst)
    assert "COG category" not in found
    assert {"Replicon", "Feature class"} <= set(found)


def test_cog_rows_show_letter_and_short_name(tmp_path: Path, cfg: GippaConfig) -> None:
    fig = _render(h.v07_out(tmp_path), cfg).figure
    words = {w.lower() for w in re.findall(r"[A-Za-z]{4,}", cog_descriptions()["J"])}
    labels = [t for t in h.texts(fig) if t.startswith("J ")]
    assert labels, h.texts(fig)
    assert any(set(re.findall(r"[a-z]{4,}", t.lower())) & words for t in labels), labels


# --- title ---------------------------------------------------------------------------------------


def test_title_counts_the_tested_items_by_kind(tmp_path: Path, cfg: GippaConfig) -> None:
    """Example A: "No group departs from detectability at q ≤ 0.05 (2 regions, 3 classes
    tested)"; the counts come from the data."""
    title = _render(h.v07_out(tmp_path), cfg).title
    assert title.startswith("No group departs from detectability at q ≤ 0.05 ("), title
    assert title.endswith(" tested)"), title
    assert re.search(r"\b2 regions\b", title), title


def test_title_singular_region(tmp_path: Path, cfg: GippaConfig) -> None:
    title = _render(h.v07_out(tmp_path, regions=h.A_REGIONS[:1]), cfg).title
    assert re.search(r"\(1 region\b", title) or re.search(r"\b1 region[,) ]", title), title
    assert "1 regions" not in title


# --- many rows: height cap and the "more rows" note ----------------------------------------------


def _drawn_rows(inst: Instance) -> int:
    return len(markers(_forest(inst)))


def test_height_grows_with_rows_up_to_a_cap(tmp_path: Path, cfg: GippaConfig) -> None:
    heights = []
    for n in (2, 30, 400, 800):
        out = h.v07_out(tmp_path / str(n), summary=[], regions=h.many_regions(n))
        heights.append(_render(out, cfg).figure.get_figheight())
    assert heights[0] < heights[1], heights
    assert heights[2] == pytest.approx(heights[3]), heights


def test_rows_beyond_the_cap_are_counted_in_a_note(tmp_path: Path, cfg: GippaConfig) -> None:
    out = h.v07_out(tmp_path, summary=[], regions=h.many_regions(400))
    inst = _render(out, cfg)
    assert len(inst.data) == 400  # the data TSV keeps every row
    drawn = _drawn_rows(inst)
    assert 0 < drawn < 400
    notes = [m for t in h.texts(inst.figure) if (m := MORE_ROWS.match(t))]
    assert len(notes) == 1, [t for t in h.texts(inst.figure) if "more" in t]
    assert int(notes[0].group(1).replace(",", "")) == 400 - drawn


def test_rows_beyond_the_cap_are_the_smallest_pp(tmp_path: Path, cfg: GippaConfig) -> None:
    out = h.v07_out(tmp_path, summary=[], regions=h.many_regions(400))
    inst = _render(out, cfg)
    shown = sorted(abs(m.x) for m in markers(_forest(inst)))
    every = sorted((abs(h.num(v)) for v in inst.data["o_minus_e_pp"]), reverse=True)
    assert shown == pytest.approx(sorted(every[: len(shown)]))


# --- AC 6.4: no text leaves the figure ------------------------------------------------------------


def test_texts_stay_in_the_figure_on_the_a_like_run(tmp_path: Path, cfg: GippaConfig) -> None:
    """Example A: the x label was clipped "observed − expected detected ge…"."""
    inst = _render(h.v07_out(tmp_path), cfg)
    assert h.out_of_bounds(inst.figure) == []
    assert XLABEL_PP in h.texts(inst.figure)


def test_texts_stay_in_the_figure_with_many_rows(tmp_path: Path, cfg: GippaConfig) -> None:
    out = h.v07_out(tmp_path, summary=h.A_BASIC, regions=h.many_regions(400))
    assert h.out_of_bounds(_render(out, cfg).figure) == []


def test_texts_stay_in_the_figure_on_sim_plots(tmp_path: Path) -> None:
    inst = h5.only(h5.render(h5.sim_copy(tmp_path), ["V07"]), "V07")
    assert h.out_of_bounds(inst.figure) == []
    assert XLABEL_PP in h.texts(inst.figure)


def test_k_like_wide_ci_does_not_squash_the_categories(tmp_path: Path, cfg: GippaConfig) -> None:
    """K: a replicon row n = 4,367 with CI ±70 genes is ±1.6 pp: every row on one pp scale, so
    the x range is set by the percentages, not by gene counts."""
    summary = [
        *h.basic("replicon", "chr1", n_genes=4367, n_detected=2792, expected=2780.0,
                 ci=(-58.0, 82.0)),
        *h.basic("cog_category", "J", n_genes=12, n_detected=3, expected=6.4, ci=(-5.0, -1.5)),
    ]  # fmt: skip
    inst = _render(h.v07_out(tmp_path, summary=summary, regions=None), cfg)
    lo, hi = _forest(inst).get_xlim()
    assert hi - lo < 100, (lo, hi)  # J spans −42..−12.5 pp; chr1 ±2 pp
    assert lo <= -41.6 and math.isfinite(hi)
