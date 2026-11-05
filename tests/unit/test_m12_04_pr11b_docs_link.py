"""M12-04 PR 11b: ``docs/reference/tables.md`` links to the 1.0 -> 1.1 migration section.

docs/specs/M12-04-outputs-revision.md § 5: the migration note goes into
``docs/user-guide/outputs.md`` > "Changes in 1.1" "with a link from ``tables.md``".
``tables.md`` is generated (``scripts/docs/gen_reference.py`` ``render_tables``; ``--check``
in CI), so the link must come from the generator: the exact line is pinned in
``test_m11_07_gen_reference.CHANGES_LINK`` / ``TABLES_GOLDEN``.
"""

from __future__ import annotations

import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
TABLES_MD = REPO / "docs" / "reference" / "tables.md"
OUTPUTS_MD = REPO / "docs" / "user-guide" / "outputs.md"
TARGET = "../user-guide/outputs.md#changes-in-11"


def _slug(heading: str) -> str:
    """Python-Markdown's default toc slug (mkdocs ``toc``): drop non-word characters, lower
    case, spaces to hyphens."""
    text = re.sub(r"[^\w\s-]", "", heading).strip().lower()
    return re.sub(r"[-\s]+", "-", text)


def test_the_anchor_is_the_changes_in_1_1_heading() -> None:
    headings = re.findall(r"^#{1,6} (.+?)\s*$", OUTPUTS_MD.read_text(encoding="utf-8"), re.M)
    assert "Changes in 1.1" in headings
    assert _slug("Changes in 1.1") == TARGET.split("#", 1)[1]


def test_tables_md_links_to_the_changes_in_1_1_section() -> None:
    text = TABLES_MD.read_text(encoding="utf-8")
    links = re.findall(r"\[([^\]]+)\]\(([^)\s]+)\)", text)
    assert any(url == TARGET for _, url in links), (
        f"docs/reference/tables.md has no link to {TARGET}"
    )
    # before the first table section, so a reader sees it first
    assert text.index(TARGET) < text.index("\n## ")
