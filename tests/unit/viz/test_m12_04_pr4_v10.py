"""M12-04 PR 4, V10 ``correlogram`` when Moran's I cannot be estimated (P-06; AC 4.4;
docs/specs/M12-04-outputs-revision.md § PR 4).

When every ``morans_i`` row of ``autocorrelation.tsv`` is NA or flagged ``not_estimable`` /
``few_values``, the whole plot is one empty state: "Too few detected genes (4 of 158) to estimate
autocorrelation." (the BB panel is not drawn from 4 genes either). A partly estimable correlogram
keeps both panels, and its NA lags and NA panels are described through ``gippa.viz.words.plain``,
never as raw flag tokens.

Interface assumed (PR 4 test author): V10 reads ``genes_tsv`` (optional) for the detected count;
the count is written ``<detected> of <genes>``.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest
from helpers import m12_04_pr4 as h
from helpers.m12_04_pr5 import render
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from gippa.viz import base
from gippa.viz.words import plain


def _v10(out: Path) -> list[base.Instance]:
    result = render(out, ["V10"])["V10"]
    assert not isinstance(result, base.Skip), result
    return list(result)


def test_d_is_one_empty_instance_with_the_detected_count(tmp_path: Path) -> None:
    found = _v10(h.d_full(tmp_path).out)
    assert len(found) == 1
    inst = found[0]
    assert inst.empty is True
    assert re.search(rf"\b{h.D_DETECTED} of {h.N_GENES}\b", inst.title), inst.title
    assert "autocorrelation" in inst.title.lower()


def test_d_has_no_raw_flag_token(tmp_path: Path) -> None:
    found = _v10(h.d_full(tmp_path).out)
    texts = [t for inst in found for t in h.joined(inst)]
    assert not any("tu_block_tu_error" in t for t in texts)
    assert h.raw_tokens(texts) == [], texts
    assert [t for t in texts if h.SNAKE.search(t)] == []
    assert len(found) == 1 and found[0].empty is True


@settings(
    max_examples=8,
    deadline=None,
    suppress_health_check=[HealthCheck.too_slow, HealthCheck.function_scoped_fixture],
)
@given(
    flags=st.sampled_from(
        ["not_estimable", "few_values", "few_values;not_estimable", h.C_MORAN_FLAGS]
    ),
    na_value=st.booleans(),
)
def test_property_no_estimable_moran_row_is_one_empty_state(
    tmp_path_factory: pytest.TempPathFactory, flags: str, na_value: bool
) -> None:
    """Every Moran's I row NA, or flagged whatever its value: one empty state, BB rows aside."""
    rows = []
    for lag in range(1, 6):
        for null in ("gene_label", "tu_block"):
            rows.append(h.ac_row("morans_i", null, lag, None if na_value else 0.3, flags))
            rows.append(h.ac_row("join_count_bb", null, lag, 2.0))
    run = h.d_full(tmp_path_factory.mktemp("ac"))
    h.write_ac(run.out, rows)
    found = _v10(run.out)
    assert len(found) == 1 and found[0].empty is True, [i.title for i in found]
    assert f"{h.D_DETECTED} of {h.N_GENES}" in found[0].title


def test_partly_estimable_keeps_both_panels_and_notes_na_lags_in_plain_words(
    tmp_path: Path,
) -> None:
    run = h.d_full(tmp_path)
    h.write_ac(run.out, h.partial_autocorrelation(na_lags=(4, 5)))
    found = _v10(run.out)
    assert len(found) == 1
    inst = found[0]
    assert inst.empty is False
    on = sorted(
        ax.get_label()
        for ax in inst.figure.axes
        if ax.axison and ax.get_label() in (f"a:{h.CONTIG}", f"b:{h.CONTIG}")
    )
    assert on == [f"a:{h.CONTIG}", f"b:{h.CONTIG}"]
    text = " ".join(h.joined(inst))
    assert plain("few_values") in text or plain("not_estimable") in text, text
    assert h.raw_tokens(h.joined(inst)) == []


def test_partial_panel_flags_go_through_plain(tmp_path: Path) -> None:
    """Moran's I estimable, every BB row NA and flagged: panel (b) names the flags in words."""
    rows = []
    for lag in range(1, 6):
        for null in ("gene_label", "tu_block"):
            rows.append(h.ac_row("morans_i", null, lag, 0.05 * lag))
            rows.append(
                h.ac_row(
                    "join_count_bb",
                    null,
                    lag,
                    None,
                    "ci_unestimated;not_estimable;tu_block_tu_error",
                )
            )
    run = h.d_full(tmp_path)
    h.write_ac(run.out, rows)
    (inst,) = _v10(run.out)
    texts = h.joined(inst)
    assert h.raw_tokens(texts) == [], texts
    assert any(plain("not_estimable") in t for t in texts), texts
