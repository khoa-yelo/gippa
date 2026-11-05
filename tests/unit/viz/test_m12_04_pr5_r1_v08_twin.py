"""M12-04 PR 5, review round 1 (L2): the V08 interactive twin matches the static V08.

- the twin carries the "detected / total" note next to its per-column counts;
- the "none" category (genes without a COG letter) is hatched in every trace type the twin draws:
  bars (nothing measured), a violin (more than 50 measured values) and a strip (50 or fewer).

Hover text with exact TSV values (``%.6g``) is data and is not checked here.

Interface assumed (test-author, round 1):

- the note is a layout annotation whose text contains exactly ``"detected / total"``;
- "hatched" means a Plotly pattern with a non-empty ``shape``: ``marker.pattern.shape`` (bars) or
  ``fillpattern.shape`` (a filled scatter, e.g. a hatched backdrop or outline). A violin and a
  marker strip have no pattern of their own, so the twin adds a patterned trace at the "none"
  column: every x of that trace lies within 0.5 of the column's x (its index in the category
  order). No trace of another column is patterned.
"""

from __future__ import annotations

from typing import Any

import pytest
from helpers import twins as tw

from gippa.viz import detection as det

NOTE = det.COUNT_NOTE


def _rows(n_none: int, measured: bool) -> list[tuple[object, ...]]:
    """J: 5 measured (a strip); none: ``n_none`` measured plus one undetected."""
    cats = [("J", "information", [(f"j{i}", 8.0 + 0.5 * i) for i in range(5)]),
            ("none", "none", [(f"n{i:03d}", 1.0 + 0.01 * i) for i in range(n_none)])]  # fmt: skip
    rows: list[tuple[object, ...]] = []
    for letter, sup, genes in cats:
        everything = [*genes, ("u1", None)] if letter == "none" else genes
        vals = sorted(v for _, v in genes)
        med = vals[len(vals) // 2] if measured else None
        rows.append((letter, sup, None, None, None, med, len(everything), len(genes), 0))
        for tag, v in everything:
            rows.append((letter, sup, tag, v if measured else None, v is not None, None, None,
                         None, None))  # fmt: skip
    return rows


def _spec(n_none: int, measured: bool) -> dict[str, Any]:
    return tw.spec("V08", tw.frame(det.V08_COLUMNS, _rows(n_none, measured)))


def _xs(t: dict[str, Any]) -> list[float]:
    return [float(x) for x in tw.values(t.get("x")) if x is not None]


def _shape(t: dict[str, Any]) -> str:
    pattern = (t.get("marker") or {}).get("pattern") or {}
    fill = t.get("fillpattern") or {}
    shape = pattern.get("shape") or fill.get("shape") or ""
    return "".join(shape) if isinstance(shape, list) else str(shape)


def _column(s: dict[str, Any], letter: str) -> float:
    ticks = s["layout"]["xaxis"]
    return float(tw.values(ticks["tickvals"])[tw.values(ticks["ticktext"]).index(letter)])


MODES = {"bars": (3, False), "violin": (60, True), "strip": (3, True)}


@pytest.mark.parametrize("mode", list(MODES))
def test_twin_has_the_detected_total_note(mode: str) -> None:
    s = _spec(*MODES[mode])
    texts = [str(a.get("text", "")) for a in s["layout"].get("annotations", [])]
    assert any(NOTE in t for t in texts), texts


@pytest.mark.parametrize("mode", list(MODES))
def test_twin_hatches_the_none_category_in_every_trace_type(mode: str) -> None:
    s = _spec(*MODES[mode])
    x_none, x_j = _column(s, "none"), _column(s, "J")
    traces = tw.traces(s)
    if mode == "violin":
        assert any(tw.ttype(t) == "violin" for t in traces), [tw.ttype(t) for t in traces]

    def at(t: dict[str, Any], x0: float) -> bool:
        xs = _xs(t)
        return bool(xs) and all(abs(x - x0) <= 0.5 for x in xs)

    hatched = [t for t in traces if _shape(t)]
    assert any(at(t, x_none) for t in hatched), [(tw.ttype(t), t.get("name")) for t in traces]
    assert not any(at(t, x_j) for t in hatched), [(tw.ttype(t), t.get("name")) for t in hatched]
