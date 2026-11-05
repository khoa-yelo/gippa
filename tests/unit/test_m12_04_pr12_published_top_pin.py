"""M12-04 PR 12 review round 1, L1: the published top-level names of ``OUT/`` are listed twice,
``publishedTop()`` in workflow/modules/local/gippa_step/main.nf (Nextflow's ``publishAs``) and
``gippa.runner.local.PUBLISHED_TOP`` (the local executor, the README); this pins them together
(pattern: tests/unit/test_local_publish_names.py)."""

from __future__ import annotations

import re
from pathlib import Path

from gippa.runner.local import PUBLISHED_TOP

ROOT = Path(__file__).resolve().parents[2]
MAIN_NF = ROOT / "workflow" / "modules" / "local" / "gippa_step" / "main.nf"


def _published_top_nf() -> list[str]:
    m = re.search(r"def\s+publishedTop\s*\(\s*\)\s*\{\s*\[(.*?)\]\s*as\s+Set\s*\}",
                  MAIN_NF.read_text(encoding="utf-8"), re.DOTALL)  # fmt: skip
    assert m, f"publishedTop() not found in {MAIN_NF}"
    return re.findall(r"'([^']+)'", m.group(1))


def test_published_top_equals_the_nextflow_list() -> None:
    names = _published_top_nf()
    assert len(names) == len(set(names)), f"duplicate in publishedTop(): {names}"
    assert set(names) == PUBLISHED_TOP
