"""M12-04 PR 4, V06 ``mapping_fdr`` degenerate states (P-02, P-03; AC 4.1, 4.2;
docs/specs/M12-04-outputs-revision.md § PR 4).

- P-02: when every ``mapping_fdr`` placement row has ``targets = 0`` (no tolerant targets) and
  there are no de novo rows, V06 is one empty state: "All 774 peptides placed by exact match; the
  tolerant-route FDR is not needed." With de novo rows (C) panels (a) and (b) collapse to one text
  line and (c) is drawn. When (b) does draw, its NA strata are listed in one stacked note.
- P-03: when every ``denovo_precision`` row is ``c_assumed``, panel (c) draws no points and holds
  "De novo precision not measured (no reference peptides); assumed c = 0.99"; an axis label says
  "(CI)" only when a CI is drawn.

Interface assumed (PR 4 test author):
- the V06 data TSV gains a column ``drawn``; panel (c) rows that are not drawn have ``drawn``
  false (``False``/``"false"``) and keep ``value``, ``ci_low``, ``ci_high``;
- V06 reads ``peptide_hits`` (optional input) for the peptide count of the empty state, which is
  written with ``style.fmt_count``;
- the panel axes keep their labels ``a``, ``b``, ``c``; a collapsed (a)/(b) leaves no axes
  labelled ``a`` or ``b`` with its axis on.
"""

from __future__ import annotations

import math
import re
from pathlib import Path

import pandas as pd
import pytest
from helpers import m12_04_pr4 as h
from helpers.m12_04_pr5 import render
from helpers.viz_mapping import (
    EXAMPLE4_TOTALS,
    denovo_row,
    full_placement,
    mapping_out,
    peptide_hit,
)
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st
from matplotlib.axes import Axes
from matplotlib.collections import LineCollection, PathCollection

from gippa.viz import base, style

ASSUMED_TEXT = "De novo precision not measured (no reference peptides); assumed c = 0.99"


def _v06(out: Path) -> list[base.Instance]:
    result = render(out, ["V06"])["V06"]
    assert not isinstance(result, base.Skip), result
    return list(result)


def _one(out: Path) -> base.Instance:
    found = _v06(out)
    assert len(found) == 1, [i.title for i in found]
    return found[0]


def _ax(inst: base.Instance, label: str) -> Axes | None:
    found = [ax for ax in inst.figure.axes if ax.get_label() == label]
    return found[0] if found else None


def _markers(ax: Axes) -> int:
    n = 0
    for line in ax.get_lines():
        if line.get_marker() not in (None, "", "None", " ", "none"):
            n += len(line.get_xdata())
    for coll in ax.collections:
        if isinstance(coll, PathCollection):
            n += len(coll.get_offsets())
    return n


def _ci_drawn(ax: Axes) -> bool:
    """An error bar: a vertical segment of non-zero length (LineCollection or 2-point line)."""
    for coll in ax.collections:
        if isinstance(coll, LineCollection):
            for seg in coll.get_segments():
                ys = [float(p[1]) for p in seg]
                xs = {float(p[0]) for p in seg}
                if len(ys) == 2 and len(xs) == 1 and ys[0] != ys[1]:
                    return True
    return False


def _labels(inst: base.Instance) -> list[str]:
    return [t for ax in inst.figure.axes for t in (ax.get_xlabel(), ax.get_ylabel()) if t]


def _false(v: object) -> bool:
    return v is False or str(v).strip().lower() in ("false", "0")


def _exact_hits(n: int) -> list[dict]:
    return [peptide_hit(f"PEPT{i:05d}K", start=101 + 60 * i) for i in range(n)]


# ---------------------------------------------------------------- AC 4.1


def test_exact_only_without_denovo_is_one_empty_state(tmp_path: Path) -> None:
    """A: 774 peptides placed by exact match, every placement row ``targets = 0``."""
    inst = _one(h.p_exact(tmp_path).out)
    assert inst.empty is True
    assert inst.suffix is None
    assert "placed by exact match" in inst.title
    assert re.search(r"\b774\b", inst.title), inst.title
    assert h.raw_tokens([inst.title]) == []


@settings(
    max_examples=8,
    deadline=None,
    suppress_health_check=[HealthCheck.too_slow, HealthCheck.function_scoped_fixture],
)
@given(n=st.integers(min_value=2, max_value=12_000))
def test_property_exact_only_title_holds_the_peptide_count(
    tmp_path_factory: pytest.TempPathFactory, n: int
) -> None:
    tmp = tmp_path_factory.mktemp("exact")
    out = mapping_out(tmp, hits=_exact_hits(n), unmapped_rows=[], fdr=h.exact_only_rows())
    inst = _one(out)
    assert inst.empty is True
    assert "placed by exact match" in inst.title
    assert re.search(rf"(?<![\d,]){re.escape(style.fmt_count(n))}(?![\d,])", inst.title), (
        n,
        inst.title,
    )


# ---------------------------------------------------------------- AC 4.2 (fixture D)


def test_d_panel_c_draws_no_markers_and_keeps_the_values(tmp_path: Path) -> None:
    inst = _one(h.d_full(tmp_path).out)
    assert inst.empty is False  # de novo rows: (c) is drawn
    ax = _ax(inst, "c")
    assert ax is not None, [a.get_label() for a in inst.figure.axes]
    assert _markers(ax) == 0
    assert "drawn" in inst.data.columns, list(inst.data.columns)
    c = inst.data[inst.data["panel"] == "c"]
    assert len(c) == h.DENOVO_BINS
    assert [float(v) for v in c["value"]] == pytest.approx([1.0] * h.DENOVO_BINS)
    assert [float(v) for v in c["ci_low"]] == pytest.approx([1.0] * h.DENOVO_BINS)
    assert [float(v) for v in c["ci_high"]] == pytest.approx([1.0] * h.DENOVO_BINS)
    assert all(_false(v) for v in c["drawn"]), list(c["drawn"])


def test_d_panel_c_says_precision_was_assumed(tmp_path: Path) -> None:
    inst = _one(h.d_full(tmp_path).out)
    assert ASSUMED_TEXT in h.joined(inst)


def test_d_no_axis_label_says_ci_when_no_ci_is_drawn(tmp_path: Path) -> None:
    inst = _one(h.d_full(tmp_path).out)
    assert not any(_ci_drawn(ax) for ax in inst.figure.axes)
    assert [t for t in _labels(inst) if "(CI)" in t] == []


def test_d_panels_a_and_b_collapse_to_one_text_line(tmp_path: Path) -> None:
    """C: no tolerant targets but de novo rows: (a) and (b) are one text line, not two empty
    panels with log axes."""
    inst = _one(h.d_full(tmp_path).out)
    drawn = [
        ax.get_label() for ax in inst.figure.axes if ax.get_label() in ("a", "b") and ax.axison
    ]
    assert drawn == []
    texts = h.joined(inst)
    assert "no peptide placements" not in texts  # the 1.0 per-panel text, once per panel
    assert "NA (no targets)" not in texts


_PREC = st.floats(min_value=0.5, max_value=0.99, allow_nan=False).map(lambda x: round(x, 3))


@settings(
    max_examples=10,
    deadline=None,
    suppress_health_check=[HealthCheck.too_slow, HealthCheck.function_scoped_fixture],
)
@given(rows=st.lists(st.tuples(_PREC, st.booleans(), st.booleans()), min_size=1, max_size=6))
def test_property_ci_label_only_with_a_drawn_ci(
    tmp_path_factory: pytest.TempPathFactory, rows: list[tuple[float, bool, bool]]
) -> None:
    """Any mix of de novo bins (with or without a CI, assumed or measured): "(CI)" is in an axis
    label only if an error bar is drawn."""
    denovo = []
    for i, (p, has_ci, assumed) in enumerate(rows):
        lo, hi = (round(p - 0.05, 3), round(min(p + 0.05, 1.0), 3)) if has_ci else (None, None)
        denovo.append(denovo_row(i + 1, p, lo, hi, "c_assumed" if assumed else None))
    out = mapping_out(
        tmp_path_factory.mktemp("ci"), fdr=[*full_placement(EXAMPLE4_TOTALS), *denovo]
    )
    inst = _one(out)
    ci = any(_ci_drawn(ax) for ax in inst.figure.axes)
    if not ci:
        assert [t for t in _labels(inst) if "(CI)" in t] == [], rows


# ---------------------------------------------------------------- the stacked NA note


def _na_texts(ax: Axes) -> list[str]:
    return [
        " ".join(t.get_text().split())
        for t in ax.texts
        if "no targets" in t.get_text() or "NA" in t.get_text().split()
    ]


def test_k_like_na_strata_are_one_stacked_note(tmp_path: Path) -> None:
    """K-like: tolerant targets in stratum 0, none in strata 1 and 2: (a) and (b) are drawn as in
    1.0, and the two NA strata share one note instead of two overlapping labels."""
    totals = {0: (900, 3, 0.0033, 880), 1: (0, 0, None, 0), 2: (0, 0, None, 0)}
    rows = [*full_placement(totals), *full_placement({}, decoy_type="shuffle")]
    inst = _one(mapping_out(tmp_path, fdr=rows))
    assert inst.empty is False
    a, b = _ax(inst, "a"), _ax(inst, "b")
    assert a is not None and a.axison and b is not None and b.axison
    assert a.get_yscale() == "log"
    assert inst.title.startswith("Placement FDR ≤ 0.01 in strata 0")
    title = " ".join(inst.title.split())
    # every NA text artist (one per stratum in 1.0), the suptitle left out
    notes = [t for ax in inst.figure.axes for t in _na_texts(ax)]
    notes += [
        " ".join(t.get_text().split())
        for t in inst.figure.texts
        if "no targets" in t.get_text() and " ".join(t.get_text().split()) != title
    ]
    assert len(notes) == 1, notes
    assert re.search(r"\b1\b", notes[0]) and re.search(r"\b2\b", notes[0]), notes[0]


def test_k_like_with_denovo_keeps_panels_a_b_and_c(tmp_path: Path) -> None:
    """Edge case: tolerant targets present: V06 is unchanged apart from the NA note."""
    rows = [*full_placement(EXAMPLE4_TOTALS), *[denovo_row(i, 0.9, 0.85, 0.95) for i in (1, 2, 3)]]
    inst = _one(mapping_out(tmp_path, fdr=rows))
    assert inst.empty is False
    assert [ax.get_label() for ax in inst.figure.axes if ax.get_label() in "abc"] == ["a", "b", "c"]
    assert _markers(_ax(inst, "c")) == 3  # type: ignore[arg-type]
    assert "drawn" in inst.data.columns
    c = inst.data[inst.data["panel"] == "c"]
    assert not any(_false(v) for v in c["drawn"]), list(c["drawn"])
    na = inst.data[(inst.data["panel"] == "b") & (inst.data["metric"] == "fdr_stratum")]
    assert (
        sum(
            1
            for v in na["value"]
            if v is None or (isinstance(v, float) and math.isnan(v)) or v is pd.NA
        )
        == 1
    )
