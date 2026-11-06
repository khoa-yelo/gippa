# ruff: noqa: RUF001, RUF002, RUF003  (figure text is Unicode: −, …)
"""M12-04 PR 6, V18 small-n (P-10; docs/specs/M12-04-outputs-revision.md § PR 6, AC 6.3, 6.4).

Interface assumed where the spec leaves it open (test author; reported to the lead):

- Panel (b) keeps the M9 axes labels ``"b:<table> / <family>"``; a family with fewer than 10
  p-values gets no such axes. One drawn text starts "Not shown (fewer than 10 tests): " and lists
  each such family as "<human name> (<n>)"; the regions ``o_minus_e`` family alone reads exactly
  "Not shown (fewer than 10 tests): region O − E (2)".
- Panel titles use human names: the regions ``o_minus_e`` histogram title contains "Region
  tests (O − E)", the autocorrelation ``join_count_bb / tu_block`` one "Join count, TU-block
  null"; no panel-(b) title holds an underscore (a raw table or key token).
- When every family has fewer than 10 p-values and nothing else is drawn (constant
  detectability, no null draws: run C), V18 is an empty-state ``Instance`` (``empty=True``)
  whose message is ``ALL_SMALL``.
- With no calibration curve, the reliability panel is one text line and the figure is shorter
  than the same run with a fitted curve.
- Panel letters read "(a)"; no text uses the M9 "a  " style. Contig names in panel titles are
  shortened with one middle ellipsis to at most 24 characters.
"""

from __future__ import annotations

import re
import tempfile
from pathlib import Path

import pytest
from helpers import m12_04_pr5 as h5
from helpers import m12_04_pr6 as h
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st
from matplotlib.axes import Axes

from gippa.config import GippaConfig, resolve_config
from gippa.io.run_view import RunView
from gippa.viz import base, style
from gippa.viz.base import Instance, Skip

NOTE_HEAD = "Not shown (fewer than 10 tests): "
REGION_NOTE = "Not shown (fewer than 10 tests): region O − E (2)"
ALL_SMALL = "Nothing to show: no test family has 10 or more tests in this run."
REGION_TITLE = "Region tests (O − E)"
JOIN_TITLE = "Join count, TU-block null"
MIN_TESTS = 10


@pytest.fixture(autouse=True)
def _style() -> None:
    style.apply()


@pytest.fixture
def cfg() -> GippaConfig:
    return resolve_config()


def _result(out: Path, cfg: GippaConfig) -> list[Instance]:
    import gippa.viz.diagnostics  # noqa: F401  (registers V18)

    result = base.registered()["V18"].fn(RunView.from_outdir(out), cfg)
    assert not isinstance(result, Skip), result
    assert len(result) == 1
    return result


def _render(out: Path, cfg: GippaConfig) -> Instance:
    return _result(out, cfg)[0]


def _hist(inst: Instance) -> dict[str, Axes]:
    return {ax.get_label()[2:]: ax for ax in inst.figure.axes if ax.get_label().startswith("b:")}


def _title(ax: Axes) -> str:
    return " ".join(t for t in (ax.get_title("left"), ax.get_title(), ax.get_title("right")) if t)


# --- AC 6.3: small families are listed, not drawn ------------------------------------------------


def test_two_region_tests_are_listed_not_drawn(tmp_path: Path, cfg: GippaConfig) -> None:
    """Example A: "regions_stats hl_shift, n = 2" with ticks 0.25, 0.75 -> the note."""
    inst = _render(h.v18_out(tmp_path, constant=False, n_regions=2, n_join_tu=12), cfg)
    assert "regions_stats / o_minus_e" not in _hist(inst)
    assert "autocorrelation / join_count_bb / tu_block" in _hist(inst)
    assert REGION_NOTE in h.texts(inst.figure), [t for t in h.texts(inst.figure) if "Not" in t]


def test_family_of_exactly_ten_is_drawn(tmp_path: Path, cfg: GippaConfig) -> None:
    inst = _render(h.v18_out(tmp_path, constant=False, n_regions=10, n_join_tu=9), cfg)
    assert "regions_stats / o_minus_e" in _hist(inst)
    assert "autocorrelation / join_count_bb / tu_block" not in _hist(inst)
    notes = [t for t in h.texts(inst.figure) if t.startswith(NOTE_HEAD)]
    assert len(notes) == 1 and notes[0].endswith("(9)"), notes


_SIZES = st.tuples(*(st.integers(0, 25) for _ in range(5)))


@settings(max_examples=15, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(sizes=_SIZES)
def test_property_drawn_iff_ten_or_more(sizes: tuple[int, int, int, int, int]) -> None:
    """Five families of 0..25 p-values: a histogram exactly for n ≥ 10, every other non-empty
    family in the one note with its count."""
    style.apply()
    keys = ("regions_stats / o_minus_e", "regions_stats / hl_shift",
            "autocorrelation / join_count_bb / tu_block",
            "autocorrelation / morans_i / gene_label", "pairs / both_detected_t")  # fmt: skip
    n_regions, n_hl, n_join, n_moran, n_pairs = sizes
    with tempfile.TemporaryDirectory() as tmp:
        out = h.v18_out(Path(tmp), constant=False, n_regions=n_regions, n_hl=n_hl,
                        n_join_tu=n_join, n_moran_gl=n_moran, n_pairs=n_pairs)  # fmt: skip
        inst = _render(out, resolve_config())
    drawn = set(_hist(inst))
    assert drawn == {k for k, n in zip(keys, sizes, strict=True) if n >= MIN_TESTS}
    small = sorted(n for n in sizes if 0 < n < MIN_TESTS)
    notes = [t for t in h.texts(inst.figure) if t.startswith(NOTE_HEAD)]
    if small:
        assert len(notes) == 1, notes
        assert sorted(int(c) for c in re.findall(r"\((\d+)\)", notes[0])) == small, notes
    else:
        assert notes == []


def test_count_ticks_are_integers(tmp_path: Path, cfg: GippaConfig) -> None:
    """10 p-values over 20 bins: counts 0..2, so default ticks would be 0.25, 0.5, …"""
    inst = _render(h.v18_out(tmp_path, constant=False, n_regions=10, n_join_tu=11), cfg)
    h.draw(inst.figure)
    assert _hist(inst)
    for key, ax in _hist(inst).items():
        lo, hi = ax.get_ylim()
        ticks = [float(y) for y in ax.get_yticks() if lo - 1e-9 <= y <= hi + 1e-9]
        assert ticks and all(y == int(y) for y in ticks), (key, ticks)


def test_every_family_small_is_an_empty_state(tmp_path: Path, cfg: GippaConfig) -> None:
    """Edge case C: constant detectability, two families of 2 and 4 tests, no null draws."""
    inst = _render(h.v18_out(tmp_path, constant=True, n_regions=2, n_join_tu=4), cfg)
    assert inst.empty, inst.title
    assert ALL_SMALL in [inst.title, *h.texts(inst.figure)]


# --- human names, letters, contigs ---------------------------------------------------------------


def test_panel_titles_use_human_names(tmp_path: Path, cfg: GippaConfig) -> None:
    out = h.v18_out(tmp_path, constant=False, n_regions=12, n_join_tu=12, n_moran_gl=12,
                    n_pairs=12)  # fmt: skip
    hist = _hist(_render(out, cfg))
    assert REGION_TITLE in _title(hist["regions_stats / o_minus_e"])
    assert JOIN_TITLE in _title(hist["autocorrelation / join_count_bb / tu_block"])
    for key, ax in hist.items():
        assert "_" not in _title(ax), (key, _title(ax))


def test_sim_plots_panel_b_titles_have_no_raw_tokens(tmp_path: Path) -> None:
    inst = h5.only(h5.render(h5.sim_copy(tmp_path), ["V18"]), "V18")
    for key, ax in _hist(inst).items():
        assert "_" not in _title(ax), (key, _title(ax))


def test_panel_letters_use_parentheses(tmp_path: Path, cfg: GippaConfig) -> None:
    inst = _render(h.v18_out(tmp_path, constant=False, n_regions=12), cfg)
    found = h.texts(inst.figure)
    assert any(t.startswith("(a) ") for t in found), found
    assert not any(re.match(r"^[a-c] {2}", t) for t in found), found


def test_no_calibration_curve_collapses_the_reliability_panel(
    tmp_path: Path, cfg: GippaConfig
) -> None:
    """Example C: a large empty top half for "detectability constant: no calibration curve"."""
    fit = _render(h.v18_out(tmp_path / "fit", constant=False, n_regions=12), cfg)
    const = _render(h.v18_out(tmp_path / "const", constant=True, n_regions=12), cfg)
    assert const.figure.get_figheight() < fit.figure.get_figheight()
    lines = [t for t in h.texts(const.figure) if "calibration curve" in t]
    assert len(lines) == 1 and "\n" not in lines[0], lines


def test_long_contig_is_shortened_in_titles(tmp_path: Path, cfg: GippaConfig) -> None:
    """Example A: "nissle1917_NZ_CP058217_1_1965000_2165…" clipped -> a middle ellipsis."""
    out = h.v18_out(tmp_path, constant=False, n_regions=12)
    h.add_island_nulls(out)
    inst = _render(out, cfg)
    titles = [_title(ax) for ax in inst.figure.axes if ax.get_label().startswith("c:islands")]
    per_contig = [t for t in titles if any(h.ELLIPSIS in w for w in t.split())]
    assert per_contig, titles
    for t in titles:
        assert h.LONG_CONTIG not in t, t
    for t in per_contig:
        tokens = [w.strip("()/,") for w in t.split() if h.ELLIPSIS in w]
        assert all(h.middle_shortened(w, h.LONG_CONTIG) for w in tokens), tokens


# --- AC 6.4: no text leaves the figure ------------------------------------------------------------


def test_texts_stay_in_the_figure_with_a_long_contig(tmp_path: Path, cfg: GippaConfig) -> None:
    out = h.v18_out(tmp_path, constant=False, n_regions=2, n_join_tu=12)
    h.add_island_nulls(out)
    assert h.out_of_bounds(_render(out, cfg).figure) == []


def test_texts_stay_in_the_figure_on_sim_plots(tmp_path: Path) -> None:
    inst = h5.only(h5.render(h5.sim_copy(tmp_path), ["V18"]), "V18")
    assert h.out_of_bounds(inst.figure) == []
