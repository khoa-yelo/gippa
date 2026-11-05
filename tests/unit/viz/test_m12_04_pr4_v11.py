"""M12-04 PR 4, V11 ``operon_pairs`` empty and partial states (P-11; AC 4.5;
docs/specs/M12-04-outputs-revision.md § PR 4).

- Every panel would be empty (C: no abundance, so no measured pair, no lfdr, no partner
  scatter): one empty state, "No measured operon pairs: partner comparisons need abundance."
- The title is built from what is drawn: no "too few measured pairs for densities" while (a)
  draws densities; (c)'s title names only the classes drawn; the legend entry for OP is
  "OP (same operon)"; NO is drawn in ``#009E73`` (not the detection "unestimated" orange), and OP
  on top with alpha 0.6.

Interface assumed (PR 4 test author): legends may sit on the figure or on an axes; a legend
label names a pair class by starting with its code (``OP …``, ``NO …``); panel axes keep the
labels ``a``, ``b``, ``c``; the (c) scatter of a class is one ``PathCollection``.
"""

from __future__ import annotations

import re
from pathlib import Path

from helpers import m12_04_pr4 as h
from helpers.m12_04_pr5 import hexc, legend_entries, render
from matplotlib.axes import Axes
from matplotlib.collections import PathCollection

from gippa.viz import base, style

EMPTY = "No measured operon pairs: partner comparisons need abundance."
CLASS = re.compile(r"^(OP|NO|CV|SW|DV)\b")
NO_COLOUR = "#009E73"


def _v11(out: Path) -> list[base.Instance]:
    result = render(out, ["V11"])["V11"]
    assert not isinstance(result, base.Skip), result
    return list(result)


def _ax(inst: base.Instance, label: str) -> Axes:
    found = [ax for ax in inst.figure.axes if ax.get_label() == label]
    assert len(found) == 1, [ax.get_label() for ax in inst.figure.axes]
    return found[0]


def _legend_classes(inst: base.Instance) -> list[str]:
    inst.figure.canvas.draw()
    return [lab for lab, _ in legend_entries(inst.figure)]


def _scatter(ax: Axes, colour: str) -> list[PathCollection]:
    return [
        c
        for c in ax.collections
        if isinstance(c, PathCollection)
        and len(c.get_offsets())
        and hexc(c.get_facecolors()[0]) == colour.upper()
    ]


def _drawn_classes(inst: base.Instance) -> set[str]:
    """Classes with a density in (a) or a point in (c) (data TSV)."""
    d = inst.data
    a = d[(d["panel"] == "a") & d["bin"].notna()]
    c = d[d["panel"] == "c"]
    return {str(x) for x in a["class"]} | {str(x) for x in c["class"]}


# ---------------------------------------------------------------- fixture D


def test_d_is_one_empty_instance(tmp_path: Path) -> None:
    found = _v11(h.d_full(tmp_path).out)
    assert len(found) == 1
    assert found[0].empty is True
    assert found[0].title == EMPTY


# ---------------------------------------------------------------- fixture P (A-like)


def test_p_title_is_built_from_what_is_drawn(tmp_path: Path) -> None:
    (inst,) = _v11(h.p_pairs(tmp_path).out)
    assert _drawn_classes(inst) == {"OP"}
    assert len(_ax(inst, "a").patches) > 0  # (a) draws the OP density
    assert "too few measured pairs" not in inst.title
    assert re.search(r"\bOP\b", inst.title), inst.title
    assert not re.search(r"\bNO\b", inst.title), inst.title


def test_p_panel_c_title_names_only_the_drawn_classes(tmp_path: Path) -> None:
    (inst,) = _v11(h.p_pairs(tmp_path).out)
    title = _ax(inst, "c").get_title(loc="left")
    assert re.search(r"\bOP\b", title), title
    assert not re.search(r"\bNO\b", title), title


def test_p_legend_labels_are_exactly_the_classes_drawn(tmp_path: Path) -> None:
    (inst,) = _v11(h.p_pairs(tmp_path).out)
    labels = _legend_classes(inst)
    assert [lab for lab in labels if lab.startswith("(c)")] == [], labels
    codes = {m.group(1) for lab in labels if (m := CLASS.match(lab))}
    assert codes == _drawn_classes(inst) == {"OP"}, labels
    assert any(lab.startswith("OP (same operon)") for lab in labels), labels


def test_no_colour_is_not_a_detection_colour(tmp_path: Path) -> None:
    (inst,) = _v11(h.p_pairs(tmp_path, measured_no=True).out)
    assert _drawn_classes(inst) == {"OP", "NO"}
    detection = {c.upper() for c in style.DETECTION.values()}
    c = _ax(inst, "c")
    assert _scatter(c, NO_COLOUR), [
        hexc(x.get_facecolors()[0]) for x in c.collections if isinstance(x, PathCollection)
    ]
    assert NO_COLOUR not in detection
    for coll in c.collections:
        if isinstance(coll, PathCollection) and len(coll.get_offsets()):
            assert hexc(coll.get_facecolors()[0]) != style.DETECTION["unestimated"].upper()
    no_lines = [p for p in _ax(inst, "a").patches if (p.get_label() or "").startswith("NO")]
    assert no_lines and all(hexc(p.get_edgecolor()) not in detection for p in no_lines)


def test_op_is_drawn_on_top_with_alpha_0_6(tmp_path: Path) -> None:
    (inst,) = _v11(h.p_pairs(tmp_path, measured_no=True).out)
    c = _ax(inst, "c")
    op = _scatter(c, style.OKABE_ITO["blue"])
    no = _scatter(c, NO_COLOUR)
    assert op and no
    assert min(x.get_zorder() for x in op) > max(x.get_zorder() for x in no)
    assert all(x.get_alpha() == 0.6 for x in op)


def test_measured_legend_has_op_and_no_labels(tmp_path: Path) -> None:
    (inst,) = _v11(h.p_pairs(tmp_path, measured_no=True).out)
    labels = _legend_classes(inst)
    codes = {m.group(1) for lab in labels if (m := CLASS.match(lab))}
    assert codes == {"OP", "NO"}, labels
    assert any(lab.startswith("OP (same operon)") for lab in labels), labels
    title = _ax(inst, "c").get_title(loc="left")
    assert re.search(r"\bOP\b", title) and re.search(r"\bNO\b", title), title
