"""M12-05b review gaps (docs/specs/M12-05b-novel-orf-genes.md rules 3, 4 and 8): three cases the
first tests left to mutants.

- Rule 4 on a LINEAR replicon: a gene that abuts the span (no base strictly between) is the
  flank at distance 0; only genes wholly on the other side of the span are skipped.
- Rule 3: the span is ``start_longest`` through the stop codon, not the best start, so a gene in
  the longest-start extension overlaps and a flank distance is measured from the longest start.
- Rule 8: an older tree read through ``RunView`` gets the nine columns as NA with the schema's
  types (``integer`` -> nullable ``Int64``; ``string`` -> object), the same dtypes a current tree
  reads with, never float NaN.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from helpers.m12_05b_context import (
    DISTANCES,
    NEW_COLUMNS,
    RC1_COLUMNS,
    cell,
    frame,
    model,
    oracle,
    orf_row,
)
from helpers.report_fixtures import run_out
from helpers.viz_out_tree import write_schema_table

from gippa.model.genes import GeneModel

L1 = 10_000


def _gc(table: pd.DataFrame, gm: GeneModel) -> pd.DataFrame:
    from gippa.novel.gene_context import gene_context

    return gene_context(table, gm)


def _row_cells(out: pd.DataFrame, i: int) -> dict[str, object]:
    return {c: cell(out[c].iloc[i]) for c in NEW_COLUMNS}


# ---------------------------------------------------------------- gap 1: distance 0, linear


@pytest.mark.parametrize(
    ("strand", "lo", "hi", "expect"),
    [
        # [400, 499): T_0001 [100, 400) ends where the span starts; T_0002 [500, 1100) is 1 away
        ("+", 400, 499, {"upstream_gene": "T_0001", "upstream_distance_nt": 0,
                         "downstream_gene": "T_0002", "downstream_distance_nt": 1}),
        # [401, 500): the span ends where T_0002 starts
        ("+", 401, 500, {"upstream_gene": "T_0001", "upstream_distance_nt": 1,
                         "downstream_gene": "T_0002", "downstream_distance_nt": 0}),
        # on - the sides swap: upstream is the higher-coordinate side
        ("-", 400, 499, {"upstream_gene": "T_0002", "upstream_distance_nt": 1,
                         "downstream_gene": "T_0001", "downstream_distance_nt": 0}),
        ("-", 401, 500, {"upstream_gene": "T_0002", "upstream_distance_nt": 0,
                         "downstream_gene": "T_0001", "downstream_distance_nt": 1}),
    ],
)  # fmt: skip
def test_abutting_gene_on_a_linear_replicon_is_the_flank_at_distance_0(
    strand: str, lo: int, hi: int, expect: dict[str, object]
) -> None:
    gm = model(c1="linear")
    row = orf_row("c1", strand, lo, hi, "intergenic", None, L1)
    got = _row_cells(_gc(frame([row]), gm), 0)
    for c, v in expect.items():
        assert got[c] == v, (c, got[c], v)
    assert got["other_overlapping_genes"] is None  # abutting is not overlapping (rule 3)
    assert got == oracle(row, gm)


# ---------------------------------------------------------------- gap 2: longest span counts


def _best_differs(strand: str, lo: int, hi: int, best_lo: int, best_hi: int) -> dict[str, Any]:
    """A row with span ``[lo, hi)`` from ``start_longest`` and the shorter best-start span
    ``[best_lo, best_hi)`` (same stop codon), both in the ``orf_table`` representation."""
    row = orf_row("c1", strand, lo, hi, "intergenic", None, L1)
    n_best = best_hi - best_lo
    assert n_best % 3 == 0 and (best_hi == hi if strand == "+" else best_lo == lo)
    row["start_best_ostir"] = best_lo if strand == "+" else best_hi - 1
    row["length_aa_best"] = n_best // 3 - 1
    assert row["start_best_ostir"] != row["start_longest"]
    assert row["length_aa_best"] < row["length_aa_longest"]
    return row


CASES = [
    # + : the longest span [390, 498) reaches into T_0001 [100, 400); the best one [402, 498)
    # does not. So T_0001 overlaps and the upstream flank wraps to T_0004 [9500, 9900) (490).
    (_best_differs("+", 390, 498, 402, 498),
     {"other_overlapping_genes": "T_0001", "upstream_gene": "T_0004",
      "upstream_distance_nt": 490, "downstream_gene": "T_0002", "downstream_distance_nt": 2}),
    # + : longest [3000, 3300), best [3150, 3300): upstream T_0003 ends at 2500, 500 not 650
    (_best_differs("+", 3000, 3300, 3150, 3300),
     {"other_overlapping_genes": None, "upstream_gene": "T_0003",
      "upstream_distance_nt": 500, "downstream_gene": "T_0004", "downstream_distance_nt": 6200}),
    # - : longest [6000, 6300), best [6000, 6150): upstream (higher side) T_0004 at 9500,
    # 3200 not 3350
    (_best_differs("-", 6000, 6300, 6000, 6150),
     {"other_overlapping_genes": None, "upstream_gene": "T_0004",
      "upstream_distance_nt": 3200, "downstream_gene": "T_0003", "downstream_distance_nt": 3500}),
    # - : longest [1100, 1304) reaches into T_0003 [1300, 2500); best [1100, 1298) does not
    (_best_differs("-", 1100, 1304, 1100, 1298),
     {"other_overlapping_genes": "T_0003", "upstream_gene": "T_0004",
      "upstream_distance_nt": 8196, "downstream_gene": "T_0002", "downstream_distance_nt": 0}),
]  # fmt: skip


@pytest.mark.parametrize(("row", "expect"), CASES)
def test_overlaps_and_flanks_follow_the_longest_start(
    row: dict[str, Any], expect: dict[str, object]
) -> None:
    gm = model()
    got = _row_cells(_gc(frame([row]), gm), 0)
    for c, v in expect.items():
        assert got[c] == v, (c, got[c], v)
    assert got == oracle(row, gm)  # the oracle reads only start_longest / length_aa_longest


def test_longest_start_rows_together_keep_their_index() -> None:
    gm = model()
    rows = [r for r, _ in CASES]
    out = _gc(frame(rows, index=[7, 3, 11, 5]), gm)
    assert list(out.index) == [7, 3, 11, 5]
    for i, (_row, expect) in enumerate(CASES):
        got = _row_cells(out, i)
        assert {c: got[c] for c in expect} == expect, i


# ---------------------------------------------------------------- gap 3: old-tree dtypes


ROW = {"orf_id": "NOVEL_chr1_101_p", "class": "intergenic", "contig": "chr1", "strand": "+",
       "frame": 3, "start_longest": 101, "end": 400, "length_aa_longest": 99,
       "length_aa_best": 99, "n_peptides": 2, "class_fdr": "0.004", "called": True,
       "overlapping_gene": "g3"}  # fmt: skip
VALUES = {"overlapping_gene_label": "lacZ", "overlap_relation": "antisense",
          "other_overlapping_genes": "g2", "upstream_gene": "g4", "upstream_gene_label": "g4",
          "upstream_distance_nt": 7, "downstream_gene": "g1", "downstream_gene_label": "dnaA",
          "downstream_distance_nt": 0}  # fmt: skip


def _read(out: Path) -> pd.DataFrame:
    from gippa.io.run_view import RunView

    frame_: pd.DataFrame = RunView.from_outdir(out).table("novel_orfs")
    return frame_


def _tree(tmp: Path, rows: list[dict[str, Any]], drop: tuple[str, ...]) -> Path:
    out = run_out(tmp)
    write_schema_table(out / "tables" / "novel_orfs.tsv", "novel_orfs", rows, drop=drop)
    return out


@pytest.mark.parametrize("n_rows", [0, 1, 2])
def test_old_tree_fill_has_the_schema_dtypes(tmp_path: Path, n_rows: int) -> None:
    """A 1.1.0rc1 header (the nine columns absent): distances read as nullable ``Int64`` and the
    string columns as object, all NA; never float NaN (``float64``)."""
    old = _read(_tree(tmp_path / "old", [ROW] * n_rows, NEW_COLUMNS))
    new = _read(_tree(tmp_path / "new", [ROW | VALUES] * n_rows, ()))
    assert [c for c in RC1_COLUMNS if c not in old.columns] == []
    for c in NEW_COLUMNS:
        want = "Int64" if c in DISTANCES else object
        assert old[c].dtype == want, (c, old[c].dtype)
        assert old[c].dtype == new[c].dtype, (c, old[c].dtype, new[c].dtype)
        assert old[c].dtype != "float64", c
        assert old[c].isna().all(), c
        assert all(v is pd.NA or v is None for v in old[c]), (c, list(old[c]))
    if n_rows:
        assert int(new["upstream_distance_nt"].iloc[0]) == 7
        assert int(new["downstream_distance_nt"].iloc[0]) == 0
