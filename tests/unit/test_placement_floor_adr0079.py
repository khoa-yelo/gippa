"""ADR-0079 Option 1: the placement stratum q-floor warning (docs/adr/ADR-0079, "Spec for tests").

For each placement list (decoy type x difference stratum m) with T_m > 0 targets, the q floor is
c/T_m (c = 1 for m >= 1; de novo stratum 0: max pi_denovo; database stratum 0: the largest ADR-0019
class pi-hat). No q of the list can fall below it. When the floor is **strictly** above
``stats.alpha.placement`` the list cannot accept anything and ``placement`` warns
``placement.stratum_q_floor_above_alpha`` with data
``{decoy_type, n_diff, n_targets, n_decoys, q_floor, alpha}``. q, ``accepted``, ``mapping_fdr``
and alpha are unchanged.

Seam assumed by the test author (the implementer may choose names inside it freely otherwise):

* ``placement_tables(costs, cfg)`` returns ``PlacementTables`` with a new field
  ``floor_breaches``: a sequence, one item per breaching list, ordered by decoy type then m
  (the order of the ``mapping_fdr`` rows). An item is a mapping or a dataclass whose fields are
  exactly the event data keys ``decoy_type, n_diff, n_targets, n_decoys, q_floor, alpha``.
* a pure rule helper ``stratum_floor_breach(*, decoy_type, n_diff, n_targets, n_decoys, c,
  alpha)`` in ``gippa.mapping.placement_fdr`` returning such an item, or None when T = 0 or
  c/T <= alpha. Row 7 (database stratum 0, c = 0.0186, T = 1 / 2) needs it: a database
  stratum-0 list with T <= 2 always has c close to 1 (every exact target is in the pi-hat
  population, so U95(k, n) >= (n - 1)/n), so that row cannot be built through costs.

The step test (the logged WARNING) is in tests/unit/test_placement_floor_step_adr0079.py.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from helpers.placement_floor import (
    DATA_KEYS,
    Ids,
    as_dict,
    cfg,
    costs_frame,
    pep_of,
    row9_costs,
    stratum,
)
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from gippa.mapping.placement_fdr import placement_tables

GOLDEN = Path(__file__).resolve().parents[1] / "data" / "placement_floor_adr0079"


# ------------------------------------------------------------------------------------ helpers


def breaches(costs: pd.DataFrame, alpha: float = 0.01) -> tuple[Any, list[dict[str, Any]]]:
    out = placement_tables(costs, cfg(alpha))
    assert hasattr(out, "floor_breaches"), "PlacementTables lacks the ADR-0079 floor_breaches"
    items = [as_dict(b) for b in out.floor_breaches]
    for b in items:
        assert set(b) == DATA_KEYS, b
    return out, items


def all_row(fdr: pd.DataFrame, dt: str, m: int) -> dict[str, Any]:
    sel = fdr[(fdr["decoy_type"] == dt) & (fdr["n_diff"] == m) & (fdr["length_class"] == "all")]
    assert len(sel) == 1, (dt, m)
    return dict(sel.iloc[0])


def check_breach(b: Mapping[str, Any], dt: str, m: int, t: int, d: int, floor: float,
                 alpha: float = 0.01) -> None:  # fmt: skip
    assert (b["decoy_type"], int(b["n_diff"]), int(b["n_targets"]), int(b["n_decoys"])) == (
        dt, m, t, d,
    ), b  # fmt: skip
    assert isinstance(b["decoy_type"], str)
    assert float(b["q_floor"]) == pytest.approx(floor, rel=1e-12), b
    assert float(b["alpha"]) == alpha, b


def rule(**kw: Any) -> dict[str, Any] | None:
    from gippa.mapping.placement_fdr import stratum_floor_breach

    res = stratum_floor_breach(**kw)
    return None if res is None else as_dict(res)


# --------------------------------------------------------------------- spec table rows 1 - 9


def test_row1_n2_reverse_m1_t13_d7_warns() -> None:
    out, bs = breaches(costs_frame(stratum(Ids(), m=1, t=13, d=7)))
    assert len(bs) == 1, bs
    check_breach(bs[0], "reverse", 1, 13, 7, 1 / 13)
    assert all_row(out.mapping_fdr, "reverse", 1)["n_accepted"] == 0
    assert not out.per_peptide["accepted"].any()


def test_row2_t100_d0_floor_equals_alpha_no_warning_all_accepted() -> None:
    out, bs = breaches(costs_frame(stratum(Ids(), m=1, t=100, d=0)))
    assert bs == []
    assert all_row(out.mapping_fdr, "reverse", 1)["n_accepted"] == 100
    assert int(out.per_peptide["accepted"].sum()) == 100


def test_row3_t99_d0_warns_none_accepted() -> None:
    out, bs = breaches(costs_frame(stratum(Ids(), m=1, t=99, d=0)))
    assert len(bs) == 1, bs
    check_breach(bs[0], "reverse", 1, 99, 0, 1 / 99)
    assert float(bs[0]["q_floor"]) == pytest.approx(0.010101, abs=1e-6)
    assert all_row(out.mapping_fdr, "reverse", 1)["n_accepted"] == 0


def test_row4_no_targets_no_warning() -> None:
    out, bs = breaches(costs_frame(stratum(Ids(), m=2, t=0, d=5)))
    assert bs == []
    assert all_row(out.mapping_fdr, "reverse", 2)["targets"] == 0


def test_row4_rule_no_targets_is_none() -> None:
    assert rule(decoy_type="reverse", n_diff=2, n_targets=0, n_decoys=5, c=1.0, alpha=0.01) is None


def test_row5_e4_shuffle_m2_t68_d0_warns() -> None:
    rows = stratum(Ids(), m=2, t=68, d=0, source="denovo", pi=0.5)
    out, bs = breaches(costs_frame(rows))
    assert len(bs) == 1, bs
    check_breach(bs[0], "shuffle", 2, 68, 0, 1 / 68)  # m >= 1: c = 1 whatever pi_denovo is
    assert all_row(out.mapping_fdr, "shuffle", 2)["n_accepted"] == 0


@pytest.mark.parametrize(("t", "warns"), [(89, True), (90, False)])
def test_row6_denovo_stratum0_c09(t: int, warns: bool) -> None:
    out, bs = breaches(costs_frame(stratum(Ids(), m=0, t=t, d=0, source="denovo", pi=0.9)))
    acc = all_row(out.mapping_fdr, "shuffle", 0)["n_accepted"]
    if warns:
        assert len(bs) == 1, bs
        check_breach(bs[0], "shuffle", 0, 89, 0, 0.9 / 89)
        assert float(bs[0]["q_floor"]) == pytest.approx(0.01011, abs=1e-5)
        assert acc == 0
    else:
        assert bs == []
        assert acc == 90


@pytest.mark.parametrize(("t", "c", "warns"), [(89, 0.9, True), (90, 0.9, False)])
def test_row6_rule(t: int, c: float, warns: bool) -> None:
    b = rule(decoy_type="shuffle", n_diff=0, n_targets=t, n_decoys=0, c=c, alpha=0.01)
    if warns:
        assert b is not None
        check_breach(b, "shuffle", 0, t, 0, c / t)
    else:
        assert b is None


def test_row7_database_stratum0_c00186_t1_warns() -> None:
    b = rule(decoy_type="reverse", n_diff=0, n_targets=1, n_decoys=0, c=0.0186, alpha=0.01)
    assert b is not None
    check_breach(b, "reverse", 0, 1, 0, 0.0186)


def test_row7_database_stratum0_c00186_t2_no_warning() -> None:
    assert rule(decoy_type="reverse", n_diff=0, n_targets=2, n_decoys=3, c=0.0186,
                alpha=0.01) is None  # fmt: skip  # 0.0093 <= 0.01


def test_row7_database_stratum0_through_placement_tables_uses_the_class_c() -> None:
    """A tiny database stratum-0 list: c is the largest class pi-hat (here U95(0, 1) = 0.95 for
    one exact target alone in its class), q floor = c/T = 0.95 > alpha."""
    out, bs = breaches(costs_frame(stratum(Ids(), m=0, t=1, d=0)))
    c = float(out.stratum0_classes.loc[out.stratum0_classes["n"] > 0, "pi_wrong"].max())
    assert c == pytest.approx(0.95)
    assert len(bs) == 1, bs
    check_breach(bs[0], "reverse", 0, 1, 0, c)


@pytest.mark.parametrize(("t", "warns"), [(13, True), (25, False)])
def test_row8_alpha_from_config(t: int, warns: bool) -> None:
    out, bs = breaches(costs_frame(stratum(Ids(), m=1, t=t, d=0)), alpha=0.05)
    if warns:
        assert len(bs) == 1, bs
        check_breach(bs[0], "reverse", 1, 13, 0, 1 / 13, alpha=0.05)
    else:
        assert bs == []  # 1/25 = 0.04 <= 0.05
        assert all_row(out.mapping_fdr, "reverse", 1)["n_accepted"] == 25


def test_row8_the_same_list_at_alpha_001_warns() -> None:
    _, bs = breaches(costs_frame(stratum(Ids(), m=1, t=25, d=0)), alpha=0.01)
    assert len(bs) == 1
    check_breach(bs[0], "reverse", 1, 25, 0, 1 / 25)


def test_row9_two_breaches_ordered_by_decoy_type_then_m() -> None:
    out, bs = breaches(row9_costs())
    assert len(bs) == 2, bs
    check_breach(bs[0], "reverse", 1, 13, 7, 1 / 13)
    check_breach(bs[1], "shuffle", 2, 68, 0, 1 / 68)
    for b in bs:
        assert all_row(out.mapping_fdr, b["decoy_type"], int(b["n_diff"]))["n_accepted"] == 0
    assert all_row(out.mapping_fdr, "reverse", 0)["n_accepted"] > 0  # the large lists accept


def test_breach_counts_match_the_mapping_fdr_all_row() -> None:
    out, bs = breaches(row9_costs())
    for b in bs:
        row = all_row(out.mapping_fdr, b["decoy_type"], int(b["n_diff"]))
        assert (int(b["n_targets"]), int(b["n_decoys"])) == (row["targets"], row["decoys"])


def test_empty_costs_no_breach() -> None:
    _, bs = breaches(costs_frame([]))
    assert bs == []


# ------------------------------------------------------------ invariant: acceptance unchanged


def golden_cases() -> dict[str, tuple[pd.DataFrame, float]]:
    """The costs tables of the golden comparison; goldens written by the pre-ADR code (main
    38e1cdba) to tests/data/placement_floor_adr0079/<name>.{per,fdr}.parquet."""
    ids = Ids()
    return {
        "row1": (costs_frame(stratum(Ids(), m=1, t=13, d=7)), 0.01),
        "row2": (costs_frame(stratum(Ids(), m=1, t=100, d=0)), 0.01),
        "row6": (costs_frame(stratum(Ids(), m=0, t=89, d=0, source="denovo", pi=0.9)), 0.01),
        "row8": (costs_frame(stratum(Ids(), m=1, t=13, d=0)), 0.05),
        "row9": (row9_costs(), 0.01),
        "mixed": (
            costs_frame(
                stratum(ids, m=0, t=60, d=3)
                + stratum(ids, m=1, t=5, d=1)
                + stratum(ids, m=0, t=40, d=2, source="denovo", pi=0.3)
                + stratum(ids, m=2, t=120, d=0, source="denovo", pi=0.3)
            ),
            0.01,
        ),
    }


@pytest.mark.parametrize("name", sorted(golden_cases()))
def test_per_peptide_and_mapping_fdr_equal_the_pre_adr_golden(name: str) -> None:
    costs, alpha = golden_cases()[name]
    out = placement_tables(costs, cfg(alpha))
    per = pd.read_parquet(GOLDEN / f"{name}.per.parquet")
    fdr = pd.read_parquet(GOLDEN / f"{name}.fdr.parquet")
    pd.testing.assert_frame_equal(out.per_peptide.reset_index(drop=True), per, rtol=1e-12)
    pd.testing.assert_frame_equal(out.mapping_fdr.reset_index(drop=True), fdr, rtol=1e-12)
    assert hasattr(out, "floor_breaches")


# ---------------------------------------------------------- invariant: soundness of the floor

LISTS: tuple[tuple[str, int], ...] = (
    ("reverse", 1), ("reverse", 2), ("shuffle", 0), ("shuffle", 1), ("shuffle", 2),
)  # fmt: skip
SOURCE = {"reverse": "database", "shuffle": "denovo"}


@st.composite
def random_costs(draw: st.DrawFn) -> tuple[pd.DataFrame, float, float]:
    """Random lists for every decoy type x m except database stratum 0 (whose c is the ADR-0019
    class estimate, pinned by row 7); de novo rows carry pi_denovo in [0.05, 1]."""
    alpha = draw(st.sampled_from([0.01, 0.02, 0.05, 0.1]))
    ids = Ids()
    rows: list[dict[str, Any]] = []
    pis: list[float] = []
    for dt, m in LISTS:
        t = draw(st.integers(0, 130))
        d = draw(st.integers(0, 8))
        src = SOURCE[dt]
        for i in ids.take(t):
            pi = draw(st.floats(0.05, 1.0)) if src == "denovo" else None
            pis += [] if pi is None else [pi]
            rows.append({"peptide": pep_of(i), "m_target": m, "source": src, "pi_denovo": pi,
                         "e_target": draw(st.floats(1e-6, 10.0))})  # fmt: skip
        for i in ids.take(d):
            pi = draw(st.floats(0.05, 1.0)) if src == "denovo" else None
            pis += [] if pi is None else [pi]
            rows.append({"peptide": pep_of(i), "m_decoy": m, "source": src, "pi_denovo": pi,
                         "e_decoy": draw(st.floats(1e-6, 10.0))})  # fmt: skip
    c_dn = max(pis) if pis else 1.0
    return costs_frame(rows), alpha, c_dn


def c_of(dt: str, m: int, c_dn: float) -> float:
    return c_dn if (dt, m) == ("shuffle", 0) else 1.0


@settings(max_examples=40, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(random_costs())
def test_floor_is_sound_and_the_breaches_are_exactly_the_lists_above_alpha(
    case: tuple[pd.DataFrame, float, float],
) -> None:
    costs, alpha, c_dn = case
    out, bs = breaches(costs, alpha)
    per, fdr = out.per_peptide, out.mapping_fdr
    want = []
    for dt, m in LISTS:
        row = all_row(fdr, dt, m)
        t = int(row["targets"])
        floor = c_of(dt, m, c_dn) / t if t else math.nan
        q = per.loc[(per["decoy_type"] == dt) & (per["m_target"] == m), "placement_q"]
        q = q.to_numpy(dtype=np.float64)
        assert len(q) == t
        if t:
            assert np.all(q >= min(1.0, floor) * (1 - 1e-12)), (dt, m, floor, q.min())
        if t and floor > alpha:
            want.append((dt, m, t, int(row["decoys"]), floor))
            assert row["n_accepted"] == 0, (dt, m)
    assert [(b["decoy_type"], int(b["n_diff"]), int(b["n_targets"]), int(b["n_decoys"]))
            for b in bs] == [w[:4] for w in want]  # fmt: skip
    for b, w in zip(bs, want, strict=True):
        assert float(b["q_floor"]) == pytest.approx(w[4], rel=1e-12)
        assert float(b["alpha"]) == alpha


@settings(max_examples=200, deadline=None)
@given(
    t=st.integers(0, 10_000),
    d=st.integers(0, 50),
    m=st.integers(0, 2),
    c=st.floats(1e-4, 1.0),
    alpha=st.floats(1e-4, 0.5),
)
def test_rule_fires_iff_targets_and_c_over_t_strictly_above_alpha(
    t: int, d: int, m: int, c: float, alpha: float
) -> None:
    b = rule(decoy_type="reverse", n_diff=m, n_targets=t, n_decoys=d, c=c, alpha=alpha)
    if t > 0 and c / t > alpha:
        assert b is not None
        check_breach(b, "reverse", m, t, d, c / t, alpha=alpha)
    else:
        assert b is None
