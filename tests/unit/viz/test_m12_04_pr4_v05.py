"""M12-04 PR 4, V05 ``mapping_sankey`` with one route and one outcome (P-15; AC 4.7;
docs/specs/M12-04-outputs-revision.md § PR 4).

When exactly one route and one outcome carry flow (A: 774 peptides, all exact and placed; B: 106
proteins, all miniprot and placed), V05 is one horizontal bar at ⅓ of the default height, titled
"All 774 peptides placed by exact match" / "All 106 proteins placed (miniprot)". With peptides
and proteins (two routes) the full alluvial is drawn.

The default V05 height is the height of the full alluvial (on the rich snapshot fixture, several
routes and outcomes). Counts are written with ``style.fmt_count``.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from helpers import m12_04_pr4 as h
from helpers import snapshot_fixture as sf
from helpers.m12_04_pr5 import render
from helpers.viz_mapping import mapping_out, peptide_hit
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from gippa.viz import base, style


def _v05(out: Path) -> base.Instance:
    result = render(out, ["V05"])["V05"]
    assert not isinstance(result, base.Skip), result
    assert len(result) == 1, [i.title for i in result]
    return result[0]


@pytest.fixture(scope="module")
def default_height() -> float:
    why = sf.missing(sf.RICH)
    if why is not None:
        pytest.fail(why)
    return float(_v05(sf.fixture_out(sf.RICH)).figure.get_figheight())


def test_a_like_one_bar_titled_exact_match(tmp_path: Path, default_height: float) -> None:
    inst = _v05(h.p_exact(tmp_path).out)
    assert inst.title == "All 774 peptides placed by exact match"
    assert inst.figure.get_figheight() <= default_height / 3 + 1e-9


def test_c_like_de_novo_peptides_one_bar(tmp_path: Path, default_height: float) -> None:
    inst = _v05(h.d_full(tmp_path).out)
    assert inst.title == "All 21 peptides placed by exact match"
    assert inst.figure.get_figheight() <= default_height / 3 + 1e-9


def test_b_like_one_bar_titled_miniprot(tmp_path: Path, default_height: float) -> None:
    inst = _v05(h.r_placed(tmp_path).out)
    assert inst.title == "All 106 proteins placed (miniprot)"
    assert inst.figure.get_figheight() <= default_height / 3 + 1e-9


def test_peptides_and_proteins_draw_the_full_plot(tmp_path: Path, default_height: float) -> None:
    """Edge case: two input classes, so two routes: the full alluvial, at the default height."""
    inst = _v05(h.both_placed(tmp_path).out)
    assert not inst.title.startswith("All ")
    assert inst.figure.get_figheight() == pytest.approx(default_height)


def test_one_route_two_outcomes_is_the_full_plot(tmp_path: Path, default_height: float) -> None:
    """Exact matches plus unexplained unmapped peptides: two outcomes, so no single bar."""
    from helpers.viz_mapping import unmapped

    hits = [peptide_hit(f"PEPT{i:04d}K", start=101 + 60 * i) for i in range(30)]
    gone = [unmapped(f"LOST{i:04d}K", "no_hit") for i in range(5)]
    inst = _v05(mapping_out(tmp_path, hits=hits, unmapped_rows=gone))
    assert not inst.title.startswith("All ")
    assert inst.figure.get_figheight() == pytest.approx(default_height)


@settings(
    max_examples=4,
    deadline=None,
    suppress_health_check=[HealthCheck.too_slow, HealthCheck.function_scoped_fixture],
)
@given(n=st.integers(min_value=2, max_value=12_000))
def test_property_single_route_title_and_height(
    tmp_path_factory: pytest.TempPathFactory, default_height: float, n: int
) -> None:
    hits = [peptide_hit(f"PEPT{i:05d}K", start=101 + 60 * i) for i in range(n)]
    inst = _v05(mapping_out(tmp_path_factory.mktemp("v05"), hits=hits, unmapped_rows=[]))
    assert inst.title == f"All {style.fmt_count(n)} peptides placed by exact match"
    assert inst.figure.get_figheight() <= default_height / 3 + 1e-9
