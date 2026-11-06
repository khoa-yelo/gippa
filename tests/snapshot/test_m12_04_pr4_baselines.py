"""M12-04 PR 4 baselines (AC 4.6; docs/specs/M12-04-outputs-revision.md § 2 and § PR 4).

PR 4 removes the text-only V03 index image, so its two baselines go: no snapshot case of either
fixture is ``V03_locus`` (the plot no longer returns the index instance), and the files
``sim_plots-V03_locus.png`` and ``sim_empty-V03_locus.png`` are deleted from
``tests/snapshot/baseline``. The locus views themselves (``V03_locus_<id>``) keep their cases.

The deletion is the lead's step after Khoa approves the PR's baselines (§ 1 rule 4), so
``test_v03_index_baselines_are_deleted`` fails until then.
"""

from __future__ import annotations

import pytest
from helpers import snapshot_fixture as sf

REMOVED = ("sim_plots-V03_locus.png", "sim_empty-V03_locus.png")


@pytest.mark.parametrize("fixture", sf.FIXTURES)
def test_v03_index_is_no_snapshot_case(fixture: str) -> None:
    why = sf.missing(fixture)
    if why is not None:
        pytest.fail(why)
    stems = [stem for _, stem in sf._rendered_pid(fixture, "V03")]
    assert "V03_locus" not in stems, stems
    if fixture == sf.RICH:  # the locus views stay
        assert len([s for s in stems if s.startswith("V03_locus_")]) == sf.LOCUS_KEEP, stems


def test_v03_index_baselines_are_deleted() -> None:
    assert [n for n in REMOVED if (sf.BASELINE_DIR / n).exists()] == []
