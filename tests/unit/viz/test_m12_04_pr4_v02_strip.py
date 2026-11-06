"""M12-04 PR 4, V02 without abundance: the detection strip (P-05; AC 4.3;
docs/specs/M12-04-outputs-revision.md § PR 4).

When ``abundance_type = none`` or no gene has a finite ``abundance_log2``, V02 draws one row of
ticks per replicon coloured by ``detection_status`` with no y axis, titled "No abundance in this
run: detection along the replicon (4 of 158 detected)". Edge cases: a linear contig runs
1…length with no wrap; an origin-crossing gene's tick sits at its start (the M9 convention);
minus-strand genes share the one row.

Interface assumed (PR 4 test author): the strip rows keep the 1.0 V02 data columns (``kind``
``rug``, ``position_bp`` 0-based in the ``Instance``, ``category`` = ``detection_status``) with
``value`` NA; one panel axes per replicon, labelled with the contig; the ticks are drawn as
vertical segments (``vlines``: a ``LineCollection``).
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import matplotlib.colors as mcolors
import pandas as pd
import pytest
from helpers import m12_04_pr4 as h
from helpers.m12_04_pr5 import render
from helpers.viz_out_tree import gene_row, make_out, write_schema_table
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st
from matplotlib.axes import Axes
from matplotlib.collections import LineCollection

from gippa.viz import base, style

STRIP = re.compile(
    r"No abundance in this run: detection along the replicon \((\d+) of (\d+) detected\)"
)
LINEAR, CIRC = "plasmid1", "chr1"
L_LIN, L_CIRC = 10_000, 20_000


def _v02(out: Path) -> base.Instance:
    result = render(out, ["V02"])["V02"]
    assert not isinstance(result, base.Skip), result
    assert len(result) == 1
    return result[0]


def _segments(ax: Axes) -> list[tuple[float, float, float, str]]:
    """``(x, y0, y1, colour)`` of every vertical tick segment of ``ax``."""
    out = []
    for coll in ax.collections:
        if not isinstance(coll, LineCollection):
            continue
        colours = coll.get_colors()
        for i, seg in enumerate(coll.get_segments()):
            xs = {round(float(p[0]), 9) for p in seg}
            ys = sorted(float(p[1]) for p in seg)
            if len(xs) == 1 and len(ys) == 2:
                colour = mcolors.to_hex(colours[i % len(colours)]).upper()
                out.append((xs.pop(), round(ys[0], 6), round(ys[1], 6), colour))
    return out


def _no_y_axis(ax: Axes) -> bool:
    ticks = [t.get_text() for t in ax.get_yticklabels() if t.get_visible() and t.get_text()]
    return ax.get_ylabel() == "" and (not ax.yaxis.get_visible() or ticks == [])


# ---------------------------------------------------------------- fixture D


def test_d_title_is_the_detection_strip(tmp_path: Path) -> None:
    inst = _v02(h.d_full(tmp_path).out)
    m = STRIP.fullmatch(inst.title)
    assert m is not None, inst.title
    assert (int(m.group(1)), int(m.group(2))) == (h.D_DETECTED, h.N_GENES)


def test_d_data_has_no_abundance_values(tmp_path: Path) -> None:
    inst = _v02(h.d_full(tmp_path).out)
    assert "abundance_log2" not in inst.data.columns or inst.data["abundance_log2"].isna().all()
    assert pd.to_numeric(inst.data["value"]).isna().all()
    assert set(inst.data["kind"]) <= {"rug", "island", "oric", "ter"}
    assert STRIP.fullmatch(inst.title), inst.title


def test_d_strip_has_no_y_axis_and_status_colours(tmp_path: Path) -> None:
    inst = _v02(h.d_full(tmp_path).out)
    assert STRIP.fullmatch(inst.title), inst.title
    panels = [ax for ax in inst.figure.axes if ax.get_label() == h.CONTIG]
    assert len(panels) == 1
    ax = panels[0]
    assert _no_y_axis(ax), (ax.get_ylabel(), [t.get_text() for t in ax.get_yticklabels()])
    segs = _segments(ax)
    assert len(segs) == h.N_GENES
    colours = {c for *_, c in segs}
    assert colours == {style.DETECTION["detected"].upper(), style.DETECTION["not_detected"].upper()}


# ---------------------------------------------------------------- edge cases


def _strip_tree(tmp: Path, genes: list[dict], name: str = "OUT") -> Path:
    out = make_out(tmp, name)
    pre = out / "preflight" / "preflight.json"
    pre.parent.mkdir(parents=True, exist_ok=True)
    pre.write_text(
        json.dumps(
            {
                "schema_version": "1",
                "gippa_version": "0.0.0",
                "inputs": [],
                "replicons": [
                    {"id": CIRC, "length": L_CIRC, "topology": "circular"},
                    {"id": LINEAR, "length": L_LIN, "topology": "linear"},
                ],
            }
        ),
        encoding="utf-8",
    )
    write_schema_table(out / "tables" / "genes.tsv", "genes", genes)
    return out


def _gene(
    tag: str, contig: str, start: int, end: int, strand: str = "+", status: str = "not_detected"
) -> dict:
    return gene_row(
        locus_tag=tag,
        contig=contig,
        start=start,
        end=end,
        strand=strand,
        detection_status=status,
        detected=status == "detected",
        abundance_log2=None,
        abundance_type="none",
        abundance_scale=None,
        crosses_origin=end > L_CIRC and contig == CIRC,
    )


def _edge_genes() -> list[dict]:
    return [
        _gene("c1", CIRC, 1001, 1600, "+", "detected"),
        _gene("c2", CIRC, 5001, 5600, "-", "detected"),
        _gene("c3", CIRC, 9001, 9600, "-", "single_peptide_unestimated"),
        _gene("cori", CIRC, L_CIRC - 99, L_CIRC + 200, "+", "not_detected"),  # crosses the origin
        _gene("l1", LINEAR, 1, 300, "+", "detected"),
        _gene("l2", LINEAR, 4001, 4600, "-", "not_detected"),
        _gene("lend", LINEAR, L_LIN - 299, L_LIN, "-", "not_detected"),
    ]


def _row(inst: base.Instance, tag: str) -> pd.Series:
    hit = inst.data[(inst.data["locus_tag"] == tag) & (inst.data["kind"] == "rug")]
    assert len(hit) == 1, (tag, inst.data)
    return hit.iloc[0]


def test_origin_crossing_gene_tick_is_at_its_start(tmp_path: Path) -> None:
    inst = _v02(_strip_tree(tmp_path, _edge_genes()))
    assert STRIP.fullmatch(inst.title), inst.title
    assert int(_row(inst, "cori")["position_bp"]) == L_CIRC - 100  # 0-based start


def test_linear_contig_runs_from_1_to_length_without_wrap(tmp_path: Path) -> None:
    inst = _v02(_strip_tree(tmp_path, _edge_genes()))
    assert STRIP.fullmatch(inst.title), inst.title
    lin = inst.data[(inst.data["contig"] == LINEAR) & (inst.data["kind"] == "rug")]
    pos = [int(v) for v in lin["position_bp"]]
    assert pos and all(0 <= p < L_LIN for p in pos), pos
    ax = next(a for a in inst.figure.axes if a.get_label() == LINEAR)
    lo, hi = ax.get_xlim()
    assert lo == pytest.approx(0.0, abs=1e-9)
    xs = [x for x, *_ in _segments(ax)]
    assert xs and all(lo <= x <= hi for x in xs)


def test_minus_strand_genes_share_the_one_row(tmp_path: Path) -> None:
    inst = _v02(_strip_tree(tmp_path, _edge_genes()))
    assert STRIP.fullmatch(inst.title), inst.title
    for contig in (CIRC, LINEAR):
        ax = next(a for a in inst.figure.axes if a.get_label() == contig)
        spans = {(y0, y1) for _, y0, y1, _ in _segments(ax)}
        assert len(spans) == 1, (contig, spans)
        assert _no_y_axis(ax)


def test_edge_tree_title_counts_and_status_colours(tmp_path: Path) -> None:
    inst = _v02(_strip_tree(tmp_path, _edge_genes()))
    m = STRIP.fullmatch(inst.title)
    assert m is not None, inst.title
    assert (int(m.group(1)), int(m.group(2))) == (3, 7)
    ax = next(a for a in inst.figure.axes if a.get_label() == CIRC)
    colours = {c for *_, c in _segments(ax)}
    want = {
        style.STATUS_COLOUR[s].upper()
        for s in ("detected", "single_peptide_unestimated", "not_detected")
    }
    assert colours == want


@settings(
    max_examples=10,
    deadline=None,
    suppress_health_check=[HealthCheck.too_slow, HealthCheck.function_scoped_fixture],
)
@given(
    statuses=st.lists(
        st.sampled_from(
            ["detected", "not_detected", "partial_excluded", "single_peptide_unestimated"]
        ),
        min_size=1,
        max_size=15,
    )
)
def test_property_strip_title_counts_the_detected_genes(
    tmp_path_factory: pytest.TempPathFactory, statuses: list[str]
) -> None:
    genes = [
        _gene(f"g{i}", CIRC, 101 + 1000 * i, 400 + 1000 * i, "+-"[i % 2], s)
        for i, s in enumerate(statuses)
    ]
    inst = _v02(_strip_tree(tmp_path_factory.mktemp("strip"), genes))
    m = STRIP.fullmatch(inst.title)
    assert m is not None, inst.title
    assert (int(m.group(1)), int(m.group(2))) == (statuses.count("detected"), len(statuses))
