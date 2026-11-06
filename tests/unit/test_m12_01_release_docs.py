"""M12-01 AC 19: ``docs/release.md`` has the build, attach and verify steps of the example bundle
(spec docs/specs/M12-01-example-dataset.md, section 5 "Release procedure"). Written in PR 8a; the
prose lands in PR 9 (``docs(release): build and attach the example bundle``), not PR 8b.

Step 13b is Khoa's sub-step of step 13: it is looked for between the ``### 13.`` and ``### 14.``
headings, so ``test_m11_05_runbook.py``'s steps 1-17 stay as they are.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
RUNBOOK = REPO / "docs" / "release.md"
ASSET = "gippa-example-<X.Y>.tar.gz"


@pytest.fixture(scope="module")
def text() -> str:
    return RUNBOOK.read_text(encoding="utf-8")


def _between(text: str, step: int) -> str:
    """Text from ``### <step>.`` to ``### <step + 1>.`` (sub-steps such as 13b included)."""
    start = re.search(rf"^### {step}\. ", text, re.M)
    end = re.search(rf"^### {step + 1}\. ", text, re.M)
    assert start and end, step
    return text[start.end() : end.start()]


def test_step_6_builds_measures_and_pins_the_bundle(text: str) -> None:
    body = _between(text, 6)
    for needle in ("scripts/slurm/build_example.sbatch", "scripts/slurm/example_routes.sbatch",
                   "expected/", "db_sources.yaml", "validation/example/bundles.tsv",
                   "asset_digest"):  # fmt: skip
        assert needle in body, needle
    assert re.search(r"[Pp]atch", body) and "tolerance" in body  # patch: re-run, no new bundle


def test_step_13b_khoa_attaches_and_verifies_the_asset(text: str) -> None:
    body = _between(text, 13)
    assert re.search(r"13b", body), "no step 13b"
    sub = body[body.index("13b") :]
    assert "**Khoa**" in sub or "(Khoa)" in sub
    assert f"gh release upload v<X.Y>.0 $GIPPA_ROOT/example/{ASSET}" in sub
    assert "--repo khoa-yelo/gippa" in sub
    assert "gh release view v<X.Y>.0 --json assets" in sub and "digest" in sub
    assert "Zenodo" in sub  # the Zenodo record holds no bundle file


def test_step_14_smoke_fetches_the_example(text: str) -> None:
    body = _between(text, 14)
    assert 'gippa example fetch --outdir "$smoke/example"' in body
    assert "needs_net" in body


def test_step_15_records_the_asset_digest(text: str) -> None:
    body = _between(text, 15)
    assert "asset_digest" in body and "bundles.tsv" in body
