# ruff: noqa: RUF003  (figure text is Unicode: −)
"""M12-04 PR 5, review round 1 (M1): V02 and V10 agree on the lag-1 Moran's I title.

Lead decision (binding): both plots use the TU-block lag-1 ``morans_i`` row only when its
``value``, ``q_value`` and envelope (``null_lo``, ``null_hi``) are all finite; otherwise both fall
back to the gene-label row and say "gene-label null". Both pick the same contig.

Interface assumed (test-author, round 1): the null is tried first, then the contigs in FASTA
order (the docstrings' "first contig with a finite TU-block row"): a TU-block row that is usable
on chr2 wins over the gene-label row on chr1. The gene-label fallback row obeys the same
finiteness rule. Titles keep the PR 5 form "Moran's I <fmt_num(value)>" and
"(<q_equals(q)>, <null> null)".
"""

from __future__ import annotations

import re
import tempfile
from collections.abc import Mapping
from pathlib import Path

import pytest
from helpers import viz_landscape as vls
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from gippa.config import resolve_config
from gippa.io.run_view import RunView
from gippa.viz import style
from gippa.viz.style import fmt_num, q_equals

REPLICONS = (("chr1", 100_000, "circular", "chromosome"), ("chr2", 60_000, "circular", "plasmid"))
VALUE = re.compile(r"Moran's I ([^\s:]+)")
CLAUSE = re.compile(r"\((q [^,]+), (TU-block|gene-label) null\)")

TU, GL = ("morans_i", "tu_block", 1), ("morans_i", "gene_label", 1)
GL_OK = (-0.12, -0.002, -0.052, 0.048, 0.13)
TU_CHR2_OK = (0.31, -0.004, -0.094, 0.086, 0.04)


@pytest.fixture(autouse=True)
def _style() -> None:
    style.apply()


Over = Mapping[tuple[str, str, int], vls.Cell]


def _titles(tmp: Path, chr1: Over, chr2: Over | None = None) -> dict[str, str]:
    genes = [vls.gene(f"a{i}", 2_000 + 3_000 * i, 18.0 + (i % 7) * 0.5) for i in range(30)]
    rows = vls.correlogram_rows("chr1", overrides=chr1)
    replicons = REPLICONS[:1]
    if chr2 is not None:
        genes += [vls.gene(f"b{i}", 2_000 + 3_000 * i, 18.0 + (i % 5) * 0.5, contig="chr2")
                  for i in range(18)]  # fmt: skip
        rows += vls.correlogram_rows("chr2", overrides=chr2, oriter=False)
        replicons = REPLICONS
    out = vls.landscape_out(tmp, replicons=replicons, genes=genes, autocorr=rows)
    view, cfg = RunView.from_outdir(out), resolve_config()
    titles = {}
    for pid in ("V02", "V10"):
        (inst,) = vls.plot_spec(pid).fn(view, cfg)
        titles[pid] = inst.title
        inst.figure.clear()
    return titles


def _parse(title: str) -> tuple[str, str, str]:
    v, c = VALUE.search(title), CLAUSE.search(title)
    assert v is not None and c is not None, title
    return v.group(1), c.group(1), c.group(2)


def _want(cell: vls.Cell, null: str) -> tuple[str, str, str]:
    return fmt_num(cell[0]), q_equals(cell[4]), null


CASES = {
    # TU-block q NA on the only contig -> gene-label row (V10 printed "q = NA, TU-block null")
    "tu_q_na": ({TU: (-0.15, -0.004, -0.094, 0.086, None), GL: GL_OK}, None, (GL_OK, "gene-label")),
    # TU-block envelope NA, q finite -> gene-label row (V02 printed the TU-block q)
    "tu_envelope_na": ({TU: (-0.15, -0.004, None, None, 0.92), GL: GL_OK}, None,
                       (GL_OK, "gene-label")),
    "tu_value_na": ({TU: (None, -0.004, -0.094, 0.086, 0.92), GL: GL_OK}, None,
                    (GL_OK, "gene-label")),
    # chr1's TU-block row unusable, chr2's fine -> chr2's TU-block row in both
    "tu_na_chr1_fine_chr2": ({TU: (-0.15, -0.004, -0.094, 0.086, None), GL: GL_OK},
                             {TU: TU_CHR2_OK}, (TU_CHR2_OK, "TU-block")),
    "tu_envelope_na_chr1_fine_chr2": ({TU: (-0.15, -0.004, None, None, 0.92), GL: GL_OK},
                                      {TU: TU_CHR2_OK}, (TU_CHR2_OK, "TU-block")),
}  # fmt: skip


@pytest.mark.parametrize("case", list(CASES))
def test_v02_and_v10_title_the_same_usable_row(tmp_path: Path, case: str) -> None:
    chr1, chr2, (cell, null) = CASES[case]
    titles = _titles(tmp_path, chr1, chr2)
    got = {pid: _parse(t) for pid, t in titles.items()}
    assert got["V02"] == got["V10"], titles
    assert got["V02"] == _want(cell, null), titles


@settings(max_examples=15, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(
    q=st.one_of(st.none(), st.floats(1e-6, 1.0)),
    lo=st.one_of(st.none(), st.just(-0.094)),
    hi=st.one_of(st.none(), st.just(0.086)),
    value=st.one_of(st.none(), st.floats(-0.6, 0.6).map(lambda v: round(v, 4))),
    two=st.booleans(),
)
def test_any_na_pattern_v02_and_v10_agree(
    q: float | None, lo: float | None, hi: float | None, value: float | None, two: bool
) -> None:
    tu = (value, -0.004, lo, hi, q)
    chr1 = {TU: tu, GL: GL_OK}
    chr2 = {TU: TU_CHR2_OK} if two else None
    with tempfile.TemporaryDirectory() as tmp:
        titles = _titles(Path(tmp), chr1, chr2)
    got = {pid: _parse(t) for pid, t in titles.items()}
    assert got["V02"] == got["V10"], titles
    usable = None not in (value, lo, hi, q)
    if usable:
        want = _want(tu, "TU-block")
    elif two:
        want = _want(TU_CHR2_OK, "TU-block")
    else:
        want = _want(GL_OK, "gene-label")
    assert got["V02"] == want, titles
