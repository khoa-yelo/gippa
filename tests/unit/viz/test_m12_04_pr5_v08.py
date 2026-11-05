"""M12-04 PR 5, P-19: V08 ``cog_abundance`` palette and counts (docs/specs/
M12-04-outputs-revision.md § PR 5, AC 5.5).

- the "none" category (genes without a COG letter) is ``#E0E0E0`` with a hatch, so it no longer
  shares ``#999999`` with "poorly characterised" (S, R), which keeps the COG grey;
- the rotated or flat counts over the columns get an axis note "detected / total";
- the title's n equals the count printed over its column.

Interface assumed (test-author, PR 5):

- The "none" category is drawn by patches (bars in the detected-fraction mode, the legend swatch
  in both modes) whose facecolor is ``#E0E0E0`` and whose hatch is non-empty. Exactly one legend
  entry has that style; the legend entry labelled ``"poorly"`` (or containing "poorly") stays
  ``#999999`` and unhatched.
- The counts over the columns stay one text per column, ``"<a>/<b>"`` (spaces allowed around the
  slash), centred on the column's x; the note is a figure or axes text containing exactly
  ``"detected / total"``.
- The title of the "highest median" form ends ``"(<m> log2, <n> genes)"``; ``<n>`` equals ``<a>``
  of that category's column label (the implementer may change either the title or the label).
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import pytest
from helpers import m12_04_pr5 as h
from helpers.viz_out_tree import gene_row, make_out, write_schema_table
from matplotlib.patches import Rectangle

from gippa.config import GippaConfig, resolve_config
from gippa.io.run_view import RunView
from gippa.viz import base, style
from gippa.viz.base import Instance, Skip

NONE_FILL = "#E0E0E0"
COG_GREY = "#999999"
COUNT = re.compile(r"^\s*(\d[\d,]*)\s*/\s*(\d[\d,]*)\s*$")


@pytest.fixture(autouse=True)
def _style() -> None:
    style.apply()


@pytest.fixture
def cfg() -> GippaConfig:
    return resolve_config()


def g(tag: str, cog: str | None, log2: float | None, *, detected: bool = True,
      atype: str = "ibaq", scale: str | None = "molar") -> dict[str, object]:  # fmt: skip
    return gene_row(locus_tag=tag, feature_type="cds", cog_category=cog, abundance_log2=log2,
                    abundance_imputed=False, abundance_type=atype, abundance_scale=scale,
                    detected=detected,
                    detection_status="detected" if detected else "not_detected")  # fmt: skip


def _render(tmp: Path, genes: list[dict[str, object]], cfg: GippaConfig) -> Instance:
    import gippa.viz.detection  # noqa: F401  (registers V08)

    out = make_out(tmp)
    write_schema_table(out / "tables" / "genes.tsv", "genes", genes)
    result = base.registered()["V08"].fn(RunView.from_outdir(out), cfg)
    assert not isinstance(result, Skip), result
    (inst,) = result
    return inst


#: J: 5 measured + 2 detected without a value + 2 undetected (title n 5 vs label "7/9" on main);
#: S: 3 measured; none: 3 measured
MEASURED = (
    [g(f"j{i}", "J", 20.0 + i) for i in range(5)]
    + [g("jd1", "J", None), g("jd2", "J", None)]
    + [g("ju1", "J", None, detected=False), g("ju2", "J", None, detected=False)]
    + [g(f"s{i}", "S", 18.0 + i) for i in range(3)]
    + [g(f"n{i}", None, 17.0 + i) for i in range(3)]
)
#: no abundance: detected fraction bars
FRACTION = (
    [g(f"j{i}", "J", None, atype="none", scale=None, detected=i < 2) for i in range(4)]
    + [g(f"s{i}", "S", None, atype="none", scale=None, detected=i < 1) for i in range(3)]
    + [g(f"n{i}", None, None, atype="none", scale=None, detected=i < 2) for i in range(3)]
)


def _plot_axes(inst: Instance) -> Any:
    (ax,) = [a for a in inst.figure.axes if a.get_xlabel() == "COG category"]
    return ax


def _column_counts(inst: Instance) -> dict[str, tuple[int, int]]:
    """``tick letter -> (a, b)`` of the count text over each column."""
    ax = _plot_axes(inst)
    inst.figure.canvas.draw()
    ticks = {round(float(t.get_position()[0])): t.get_text() for t in ax.get_xticklabels()}
    out: dict[str, tuple[int, int]] = {}
    for t in ax.texts:
        m = COUNT.match(t.get_text())
        if m:
            x = getattr(t, "xy", t.get_position())[0]
            out[ticks[round(float(x))]] = (int(m.group(1).replace(",", "")),
                                           int(m.group(2).replace(",", "")))  # fmt: skip
    return out


def test_title_n_equals_the_count_over_its_column(tmp_path: Path, cfg: GippaConfig) -> None:
    """AC 5.5."""
    inst = _render(tmp_path, MEASURED, cfg)
    m = re.search(r"\(\S+ log2, (\d[\d,]*) genes?\)", inst.title)
    assert m, inst.title
    assert inst.title.startswith("J ")
    counts = _column_counts(inst)
    assert int(m.group(1).replace(",", "")) == counts["J"][0], (inst.title, counts)


def test_counts_get_the_detected_total_note(tmp_path: Path, cfg: GippaConfig) -> None:
    inst = _render(tmp_path, MEASURED, cfg)
    texts = h.figure_texts(inst)
    assert any("detected / total" in t for t in texts), texts


def _legend(inst: Instance) -> list[tuple[str, tuple[str, str, str | None]]]:
    return [(label, h.patch_style(hd)) for label, hd in h.legend_entries(inst.figure)]


@pytest.mark.parametrize("genes", [MEASURED, FRACTION], ids=["abundance", "fraction"])
def test_none_is_light_grey_and_hatched_poorly_keeps_the_cog_grey(
    tmp_path: Path, cfg: GippaConfig, genes: list[dict[str, object]]
) -> None:
    inst = _render(tmp_path, genes, cfg)
    entries = _legend(inst)
    none = [e for e in entries if e[1][0] == NONE_FILL]
    assert len(none) == 1, entries
    assert none[0][1][2], f"the 'none' legend swatch has no hatch: {none}"
    poorly = [e for e in entries if "poorly" in e[0]]
    assert len(poorly) == 1, entries
    assert poorly[0][1][0] == COG_GREY and not poorly[0][1][2]
    assert style.cog_colour("S") == COG_GREY


def test_none_bar_is_light_grey_and_hatched_in_the_fraction_mode(
    tmp_path: Path, cfg: GippaConfig
) -> None:
    inst = _render(tmp_path, FRACTION, cfg)
    ax = _plot_axes(inst)
    ticks = [t.get_text() for t in ax.get_xticklabels()]
    bars = {round(p.get_x() + p.get_width() / 2): p for p in ax.patches
            if isinstance(p, Rectangle) and p.get_height() > 0}  # fmt: skip
    none_bar = bars[ticks.index("none")]
    s_bar = bars[ticks.index("S")]
    assert h.hexc(none_bar.get_facecolor()) == NONE_FILL
    assert none_bar.get_hatch()
    assert h.hexc(s_bar.get_facecolor()) == COG_GREY
    assert not s_bar.get_hatch()
