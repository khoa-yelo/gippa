# ruff: noqa: RUF001, RUF002, RUF003  (figure text is Unicode: −)
"""M12-04 PR 5, review round 1 (L1): numbers in figure text go through the ``style`` formatters.

- V05: the title's totals go through ``fmt_count`` (45231 -> "45,231").
- V15: ``determinants._f2`` writes zero as ``fmt_num`` does, never "0.00" or "−0.00".
- V02: the legend's LOESS span goes through ``fmt_num`` or ``fmt_alpha``, not ``:g``.

Interface assumed (test-author, round 1): ``mapping_sankey.title(counts)`` keeps its signature;
``determinants._f2(v)`` stays the V15 two-decimal writer; the V02 legend text stays
"trend (LOESS, span <frac>)".
"""

from __future__ import annotations

from pathlib import Path

import pytest
from helpers import m12_04_pr5 as h
from helpers import viz_landscape as vls

from gippa.config import resolve_config
from gippa.io.run_view import RunView
from gippa.viz import style
from gippa.viz.style import fmt_alpha, fmt_count, fmt_num


@pytest.fixture(autouse=True)
def _style() -> None:
    style.apply()


def test_v05_title_counts_are_comma_grouped() -> None:
    from gippa.viz import mapping_sankey as ms

    counts = dict.fromkeys(ms.GRID, 0)
    counts["peptide", "exact", "placed"] = 40_000
    unmapped = [k for k in ms.GRID if k[0] == "peptide" and k[1] == "unmapped"]
    counts[unmapped[0]] = 5_231
    counts["protein", "miniprot", "placed"] = 12_345
    text = ms.title(counts)
    assert f"of {fmt_count(45_231)} input peptides" in text, text
    assert "of 45,231 input peptides" in text, text
    assert f"of {fmt_count(12_345)} proteins" in text, text


@pytest.mark.parametrize("v", [0.0, -0.0])
def test_v15_zero_is_written_as_fmt_num_writes_it(v: float) -> None:
    from gippa.viz.determinants import _f2

    got = _f2(v)
    assert got in (fmt_num(0), fmt_num(0.0)), got
    assert "−" not in got and "-" not in got, got


def test_v02_legend_span_uses_the_style_formatter(tmp_path: Path) -> None:
    frac = 0.1234567  # ":g" writes 0.123457; fmt_num 0.12; fmt_alpha 0.1234567
    cfg = resolve_config(overrides=[f"context.islands.detrend_span_frac={frac}"])
    genes = [vls.gene(f"g{i}", 2_000 * i + 1_000, float(i % 5)) for i in range(40)]
    out = vls.landscape_out(tmp_path, genes=genes, islands=[])
    (inst,) = vls.plot_spec("V02").fn(RunView.from_outdir(out), cfg)
    texts = h.figure_texts(inst)
    allowed = {f"trend (LOESS, span {fmt_num(frac)})", f"trend (LOESS, span {fmt_alpha(frac)})"}
    spans = [t for t in texts if t.startswith("trend (LOESS")]
    assert spans and set(spans) <= allowed, spans
