# ruff: noqa: RUF003  (figure text is Unicode: −, ≤, ≥, Δ)
"""Shared fixtures for the M12-04 PR 4 tests: degenerate and empty plot states (P-02, P-03,
P-05, P-06, P-11, P-12, P-15 and the plot half of P-14; docs/specs/M12-04-outputs-revision.md
§ PR 4).

The PR 1 fixtures P, R and D (``tests/unit/report/test_m12_04_pr1_fixtures.py``) carry no
``autocorrelation.tsv`` or ``pairs.tsv``, no de novo precision rows and no tolerant-route
targets, and R's protein hits are not assigned to a gene (so V05 counts none as placed). The
builders here add what the evidence runs A, B and C had, so the spec's "Before" column is what
the 1.0 plots draw on them:

- ``p_exact(tmp)``: P with every ``mapping_fdr`` placement row at ``targets = 0`` (both decoy
  types, every length class and stratum 0-2) and no de novo rows: A (774 peptides, exact only);
- ``d_full(tmp)``: D with the same exact-only placement rows plus 10 ``denovo_precision`` rows
  at precision 1, CI 1-1, flag ``c_assumed`` (C's panel c); a C-like ``autocorrelation.tsv``
  (every ``morans_i`` row NA and flagged ``ci_unestimated;few_values;not_estimable;
  tu_block_tu_error``, the BB join-count rows finite with spikes at lags 10 and 14); and a
  ``pairs.tsv`` of 20 OP and 10 NO pairs with no abundance (no ``delta_log2``, no ``lfdr``);
- ``p_pairs(tmp, measured_no)``: P with 20 measured OP pairs and 10 NO pairs, measured only
  when ``measured_no`` (A: NO unmeasured, so (a) and (c) draw OP only);
- ``r_placed(tmp)``: R with its 106 proteins assigned to genes (B: "106 proteins placed");
- ``both_placed(tmp)``: P plus 12 proteins assigned to genes (two input classes, two routes).

Text helpers: ``joined`` (every figure text, whitespace collapsed), ``raw_tokens`` (flag and
skip tokens of ``gippa.viz.words`` found in a text, as whole words) and ``NOTHING`` (the P-14
empty-state prefix).
"""

from __future__ import annotations

import re
from collections.abc import Iterable, Sequence
from pathlib import Path
from typing import Any

from unit.report.test_m12_04_pr1_fixtures import (
    CONTIG,
    Run,
    fixture_d,
    fixture_p,
    fixture_r,
    write_proteins,
)
from unit.report.test_m12_04_pr1_fixtures import (
    N_GENES as N_GENES,  # re-exported: the PR 4 tests count against it
)

from helpers.m12_04_pr2 import flag_tokens, internal_refs
from helpers.m12_04_pr5 import figure_texts
from helpers.viz_mapping import denovo_row, full_placement
from helpers.viz_out_tree import write_schema_table

#: the P-14 empty-state prefix (AC 4.8)
NOTHING = "Nothing to show:"
#: D: genes detected of N_GENES (fixture D detects the first 4)
D_DETECTED = 4
#: C's de novo precision bins (10 bins at precision 1, CI 1-1, c assumed 0.99)
DENOVO_BINS = 10
C_VALUE = 0.99
#: the C flags of every ``morans_i`` row (the run-C V10 panel (a) text)
C_MORAN_FLAGS = "ci_unestimated;few_values;not_estimable;tu_block_tu_error"
N_OP, N_NO = 20, 10


# ---------------------------------------------------------------- mapping_fdr


def exact_only_rows() -> list[dict[str, Any]]:
    """``placement`` rows with ``targets = 0`` everywhere, both decoy types (A and C)."""
    return [*full_placement({}), *full_placement({}, decoy_type="shuffle")]


def c_assumed_rows(n: int = DENOVO_BINS) -> list[dict[str, Any]]:
    return [denovo_row(i + 1, 1.0, 1.0, 1.0, "c_assumed") for i in range(n)]


def write_fdr(out: Path, rows: Sequence[dict[str, Any]]) -> None:
    write_schema_table(out / "mapping" / "mapping_fdr.tsv", "mapping_fdr", rows)


# ---------------------------------------------------------------- autocorrelation


def ac_row(
    statistic: str,
    null: str,
    lag: int,
    value: float | None,
    flags: str | None = "ci_unestimated",
    contig: str = CONTIG,
    q: float | None = 0.5,
) -> dict[str, Any]:
    finite = value is not None
    return {
        "contig": contig,
        "statistic": statistic,
        "lag": lag,
        "null_model": null,
        "value": value,
        "null_mean": 0.0 if finite else None,
        "null_lo": -0.2 if finite else None,
        "null_hi": 0.2 if finite else None,
        "stat": statistic,
        "effect": value,
        "ci_low": None,
        "ci_high": None,
        "p_value": q if finite else None,
        "q_value": q if finite else None,
        "method": "permutation",
        "n_perm": 199,
        "flags": flags,
    }


def c_autocorrelation(lags: Iterable[int] = range(1, 21)) -> list[dict[str, Any]]:
    """C: every Moran's I row NA and flagged; the BB rows finite, spikes at lags 10 and 14."""
    rows = []
    for lag in lags:
        for null in ("gene_label", "tu_block"):
            rows.append(ac_row("morans_i", null, lag, None, C_MORAN_FLAGS))
            spike = 3.0 if lag in (10, 14) else 0.1
            rows.append(ac_row("join_count_bb", null, lag, spike, "ci_unestimated"))
    return rows


def partial_autocorrelation(
    na_lags: Sequence[int] = (4, 5), n_lags: int = 5, na_flags: str = "few_values;not_estimable"
) -> list[dict[str, Any]]:
    """Moran's I estimable at some lags, NA (flagged) at ``na_lags``; BB finite throughout.

    Every row carries a flag: a ``flags`` column mixing NA and text is read back as ``float``
    NaN cells, which 1.0 ``correlogram.drawn`` cannot split (a separate bug, see the PR 4
    test-author report)."""
    rows = []
    for lag in range(1, n_lags + 1):
        for null in ("gene_label", "tu_block"):
            if lag in na_lags:
                rows.append(ac_row("morans_i", null, lag, None, na_flags))
            else:
                rows.append(ac_row("morans_i", null, lag, 0.05 * lag))
            rows.append(ac_row("join_count_bb", null, lag, 0.2))
    return rows


def write_ac(out: Path, rows: Sequence[dict[str, Any]]) -> None:
    write_schema_table(out / "tables" / "autocorrelation.tsv", "autocorrelation", rows)


# ---------------------------------------------------------------- pairs


def pair_row(i: int, cls: str, a: float | None, b: float | None) -> dict[str, Any]:
    measured = a is not None and b is not None
    return {
        "pair_id": f"P{i:04d}",
        "pair_class": cls,
        "gene_a": f"g{2 * i:05d}",
        "gene_b": f"g{2 * i + 1:05d}",
        "strand": "+",
        "intergenic_bp": 40,
        "detected_a": a is not None,
        "detected_b": b is not None,
        "abundance_log2_a": a,
        "abundance_log2_b": b,
        "delta_log2": (b - a) if measured else None,  # type: ignore[operator]
        "lfdr": 0.5 if measured else None,
        "flags": None,
    }


def write_pairs(out: Path, *, measured_op: bool, measured_no: bool) -> None:
    rows = []
    for i in range(N_OP):
        a = 20.0 + 0.1 * i if measured_op else None
        rows.append(pair_row(i, "OP", a, None if a is None else a + 0.05 * (i % 7)))
    for j in range(N_NO):
        a = 18.0 + 0.2 * j if measured_no else None
        rows.append(pair_row(N_OP + j, "NO", a, None if a is None else a + 0.4 * (j % 5)))
    write_schema_table(out / "tables" / "pairs.tsv", "pairs", rows)


# ---------------------------------------------------------------- fixtures


def p_exact(tmp: Path) -> Run:
    run = fixture_p(tmp, "P_exact")
    write_fdr(run.out, exact_only_rows())
    return run


def d_full(tmp: Path) -> Run:
    run = fixture_d(tmp, "D_full")
    write_fdr(run.out, [*exact_only_rows(), *c_assumed_rows()])
    write_ac(run.out, c_autocorrelation())
    write_pairs(run.out, measured_op=False, measured_no=False)
    return run


def p_pairs(tmp: Path, *, measured_no: bool = False) -> Run:
    run = fixture_p(tmp, "P_pairs_no" if measured_no else "P_pairs")
    write_pairs(run.out, measured_op=True, measured_no=measured_no)
    return run


def write_assigned_proteins(out: Path, n: int, route: str = "miniprot") -> None:
    """``n`` proteins assigned to genes (``gene`` set, as M3 writes placed proteins)."""
    rows = [
        {
            "protein_id": f"JNOIEE_{i:05d}",
            "protein_len_aa": 100,
            "route": route,
            "contig": CONTIG,
            "start": 101 + 1000 * i,
            "end": 400 + 1000 * i,
            "strand": "+",
            "rank": 1,
            "identity": 1.0,
            "query_coverage": 1.0,
            "class": "identical",
            "gene": f"g{i:05d}",
        }
        for i in range(n)
    ]
    write_schema_table(out / "mapping" / "protein_hits.tsv", "protein_hits", rows)


def r_placed(tmp: Path) -> Run:
    run = fixture_r(tmp, "R_placed")
    write_assigned_proteins(run.out, 106)
    return run


def both_placed(tmp: Path) -> Run:
    run = fixture_p(tmp, "both_placed", roles=("genome", "peptides", "proteins"))
    write_proteins(run.out, 0)
    write_assigned_proteins(run.out, 12)
    return run


# ---------------------------------------------------------------- text


def joined(inst: Any) -> list[str]:
    """``figure_texts`` with whitespace (line breaks too) collapsed to single spaces."""
    return [" ".join(t.split()) for t in figure_texts(inst)]


def raw_tokens(texts: Iterable[str]) -> list[str]:
    """Flag/skip tokens (``gippa.viz.words`` keys with ``_`` or ``:``) and internal references."""
    found: list[str] = []
    for t in texts:
        found += flag_tokens(t) + internal_refs(t)
    return sorted(set(found))


#: any snake_case word (a raw token the dictionary may not know yet)
SNAKE = re.compile(r"\b[a-z0-9]+_[a-z0-9_]+\b")
