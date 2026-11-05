# ruff: noqa: RUF002, RUF003  (figure text is Unicode: −, ≤)
"""M12-04 PR 5, P-04: the V02 and V10 titles report the lag-1 Moran's I with the TU-block null
and name it (docs/specs/M12-04-outputs-revision.md § PR 5, AC 5.3; ADR-0075 item 4).

Spec examples (run A): V02 "lag-1 Moran's I −0.15 (q = 0.92, TU-block null)"; V10 "Lag-1
Moran's I −0.15: within the TU-block envelope (q = 0.92, TU-block null)". Edge case: an old run
without TU-block rows falls back to the gene-label row and names "gene-label null". The
gene-label rows stay in ``autocorrelation.tsv`` and in V10's envelopes (data rows).

Interface assumed (test-author, PR 5): the clause is ``"(q = <fmt_q(q)>, TU-block null)"`` (or
``"gene-label null"`` in the fallback), with the Moran's I written ``fmt_num(value)``; when
``fmt_q`` gives ``"< 0.001"`` the clause reads ``"(q < 0.001, TU-block null)"``. The q is the
``q_value`` of the ``morans_i`` lag-1 row of that null for the first replicon in FASTA order.
"""

from __future__ import annotations

import tempfile
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import pytest
from helpers import m12_04_pr5 as h
from helpers import viz_landscape as vls
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from gippa.config import resolve_config
from gippa.io.run_view import RunView
from gippa.viz import style
from gippa.viz.style import fmt_num, fmt_q


@pytest.fixture(autouse=True)
def _style() -> None:
    style.apply()


def clause(q: float, null: str) -> str:
    text = fmt_q(q)
    return f"(q {text}, {null})" if text.startswith("<") else f"(q = {text}, {null})"


def _lag1(rows: list[dict[str, str]], null: str) -> tuple[float, float]:
    (r,) = [r for r in rows if r["statistic"] == "morans_i" and r["lag"] == "1"
            and r["null_model"] == null and r["contig"] == "chr1"]  # fmt: skip
    return float(r["value"]), float(r["q_value"])


# --- sim_plots (the P-like rich fixture) ---------------------------------------------------


@pytest.fixture(scope="module")
def sim() -> dict[str, Any]:
    return h.render(h.SIM_OUT, ["V02", "V10"])


@pytest.mark.parametrize("pid", ["V02", "V10"])
def test_titles_use_the_tu_block_null_and_name_it(sim: dict[str, Any], pid: str) -> None:
    """AC 5.3 on sim_plots: TU-block lag-1 q 0.97028 -> "0.97" (gene-label: 0.13986)."""
    value, q = _lag1(h.read_table(h.SIM_OUT, "tables/autocorrelation.tsv"), "tu_block")
    inst = h.only(sim, pid)
    assert clause(q, "TU-block null") in inst.title, inst.title
    assert f"Moran's I {fmt_num(value)}" in inst.title, inst.title
    assert "gene-label" not in inst.title
    assert " ".join(inst.title.split()) in [" ".join(t.split()) for t in h.figure_texts(inst)]


def test_v10_keeps_the_gene_label_envelope(sim: dict[str, Any]) -> None:
    inst = h.only(sim, "V10")
    assert "TU-block null" in inst.title, inst.title
    assert {"gene_label", "tu_block"} <= set(inst.data["null_model"])


# --- the fallback (old runs without TU-block rows) ------------------------------------------------


@pytest.fixture(scope="module")
def no_tu_block(tmp_path_factory: pytest.TempPathFactory) -> tuple[Path, dict[str, Any]]:
    out = h.sim_copy(tmp_path_factory.mktemp("no_tu"))
    h.rewrite_table(out, "tables/autocorrelation.tsv",
                    lambda r: None if r["null_model"] == "tu_block" else r)  # fmt: skip
    return out, h.render(out, ["V02", "V10"])


@pytest.mark.parametrize("pid", ["V02", "V10"])
def test_without_tu_block_rows_the_title_names_the_gene_label_null(
    no_tu_block: tuple[Path, dict[str, Any]], pid: str
) -> None:
    out, results = no_tu_block
    value, q = _lag1(h.read_table(out, "tables/autocorrelation.tsv"), "gene_label")
    inst = h.only(results, pid)
    assert clause(q, "gene-label null") in inst.title, inst.title
    assert f"Moran's I {fmt_num(value)}" in inst.title, inst.title
    assert "TU-block null" not in inst.title


# --- property: the title's q is the TU-block q, whatever the gene-label q -------------------------


def _out(tmp: Path, overrides: Mapping[tuple[str, str, int], vls.Cell]) -> Path:
    genes = [vls.gene(f"g{i}", 2_000 + 3_000 * i, 18.0 + (i % 7) * 0.5) for i in range(30)]
    for gr in genes:
        gr["abundance_scale"] = "molar"
    rows = vls.correlogram_rows(overrides=overrides)
    return vls.landscape_out(tmp, genes=genes, autocorr=rows)


@settings(max_examples=8, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(
    value=st.floats(-0.6, 0.6).map(lambda v: round(v, 6)),
    q_tu=st.one_of(st.floats(1e-6, 1.0), st.sampled_from([0.0004, 0.05, 0.92, 1.0])),
    q_gl=st.floats(1e-6, 1.0),
)
def test_any_q_pair_titles_the_tu_block_q(value: float, q_tu: float, q_gl: float) -> None:
    cfg = resolve_config()
    over = {("morans_i", "tu_block", 1): (value, -0.004, -0.094, 0.086, q_tu),
            ("morans_i", "gene_label", 1): (value, -0.002, -0.052, 0.048, q_gl)}  # fmt: skip
    with tempfile.TemporaryDirectory() as tmp:
        view = RunView.from_outdir(_out(Path(tmp), over))
        for pid in ("V02", "V10"):
            (inst,) = vls.plot_spec(pid).fn(view, cfg)
            assert clause(q_tu, "TU-block null") in inst.title, (pid, inst.title)
            assert f"Moran's I {fmt_num(value)}" in inst.title, (pid, inst.title)
            inst.figure.clear()
