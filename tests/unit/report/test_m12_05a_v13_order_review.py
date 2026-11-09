"""M12-05a AC 10, review round 1 test gap (docs/specs/M12-05a-novel-orf-plot.md § 6.3):
``gippa.report.content.display_order`` with the index rows not in ``order`` order and the
stale stems given in reverse, so both the sort of the index rows by ``order`` and the
alphabetical sort of the stale stems are exercised (the earlier tests passed both pre-sorted).
"""

from __future__ import annotations

import random

import pytest
from helpers import m12_05a as h

MAIN = "V13_novel_orfs"
INDEX = "V13_novel_orfs_locus"
PREFIX = "V13_novel_orfs_locus_"
STALE = [
    "V13_novel_orfs_locus_NOVEL_gone_1_p",
    "V13_novel_orfs_locus_AAA_old",
    "V13_novel_orfs_extra",
    "V13_novel_orfs_locus_zzz_old",
]


def _rows(ids: list[str]) -> list[dict[str, str | None]]:
    out: list[dict[str, str | None]] = []
    for k, oid in enumerate(ids, start=1):
        out.append({"order": str(k), "orf_id": oid, "class": "antisense", "called": "false",
                    "class_fdr": "0.4", "best_peptide_q": "0.001", "n_peptides": "1",
                    "rank": None, "rendered": "true"})  # fmt: skip
    return out


def _shuffled(rows: list[dict[str, str | None]], seed: int) -> list[dict[str, str | None]]:
    out = list(rows)
    random.Random(seed).shuffle(out)
    return out


@pytest.mark.parametrize("seed", [0, 1, 2, 3])
def test_display_order_sorts_shuffled_index_rows_and_stale_stems(seed: int) -> None:
    from gippa.report.content import display_order

    ids = list(h.E1_ORDER)
    rows = _shuffled(_rows(ids), seed)
    assert [r["orf_id"] for r in rows] != ids  # the TSV order is not the display order
    stems = sorted([MAIN, INDEX, *(PREFIX + i for i in ids), *STALE], reverse=True)
    want = [MAIN, INDEX, *(PREFIX + i for i in ids), *sorted(STALE)]
    assert display_order(stems, rows) == want


def test_display_order_reversed_index_rows() -> None:
    from gippa.report.content import display_order

    ids = list(h.E1_ORDER)
    rows = _rows(ids)[::-1]
    stems = [*reversed(STALE), *(PREFIX + i for i in ids), INDEX, MAIN]
    want = [MAIN, INDEX, *(PREFIX + i for i in ids), *sorted(STALE)]
    assert display_order(stems, rows) == want


def test_display_order_numeric_not_text_order() -> None:
    """``order`` 10 sorts after 9 (a text sort would put "10" before "2")."""
    from gippa.report.content import display_order

    ids = [f"NOVEL_chr1_{1_001 + 3_000 * i}_p" for i in range(12)]
    rows = _shuffled(_rows(ids), 7)
    stems = sorted([MAIN, INDEX, *(PREFIX + i for i in ids)], reverse=True)
    assert display_order(stems, rows) == [MAIN, INDEX, *(PREFIX + i for i in ids)]
