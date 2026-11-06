"""The 1.1 breaking-change notes (M12-04 PR 11b/11c), wherever they are now.

Before release prep they are ``changelog.d/*.removal.md`` fragments; the release-prep towncrier
build (docs/release.md step 6) consumes them into the ``### Removals and breaking changes``
part of ``## GIPPA 1.1.0`` in CHANGELOG.md.
"""

from __future__ import annotations

import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
SECTION = re.compile(r"(?ms)^## GIPPA 1\.1\.0 \(.*?(?=^## GIPPA |\Z)")
REMOVALS = re.compile(r"(?ms)^### Removals and breaking changes\n(.*?)(?=^### |\Z)")


def removal_texts() -> list[str]:
    """One text per fragment; after the build, the built removals section as one text."""
    frags = sorted((REPO / "changelog.d").glob("*.removal.md"))
    if frags:
        return [p.read_text(encoding="utf-8") for p in frags]
    section = SECTION.search((REPO / "CHANGELOG.md").read_text(encoding="utf-8"))
    removals = REMOVALS.search(section.group(0)) if section else None
    return [removals.group(1)] if removals else []
