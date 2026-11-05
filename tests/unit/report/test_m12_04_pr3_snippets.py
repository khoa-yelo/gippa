"""M12-04 PR 3, behaviour 1 and AC 3.1: keyed ``<!-- report:begin [key] -->`` blocks in
``docs/methods`` and the snippet generator (docs/specs/M12-04-outputs-revision.md § PR 3).

Interface pinned here (PR 3 test author, beyond the spec text):
- a block opens with ``<!-- report:begin -->`` or ``<!-- report:begin <key> -->`` (key
  ``[A-Za-z0-9_-]+``, one space each side) and closes with ``<!-- report:end -->``; a page may
  hold several keyed blocks;
- ``scripts/dev/sync_report_snippets.py`` writes ``<topic>.html`` for the unkeyed block and
  ``<topic>.<key>.html`` for each keyed one; ``--check`` names ``<topic>.<key>`` when that
  snippet differs, is missing or is an orphan; ``main(argv) -> int`` stays importable;
- every block of the PR 3 pages has ≤ 120 words (whitespace-separated tokens of the Markdown
  between the markers); the long ``islands`` block predates PR 3 and stays as it is;
- ``detection-tiers.md`` holds one keyed block per tier of ``gippa.preflight.evidence.Tier``;
  ``denovo-precision`` and ``abundance`` are unkeyed (AC 3.2 names them so);
- snippet text is user-facing prose (rule 2, 3): no ``ADR-NNNN``, ``§``, ``GIPPA_SPEC``,
  ``docs/methods``, ``/home/``, ``farmshare``, ``local/`` or ``$GIPPA_ROOT``.
"""

from __future__ import annotations

import importlib.util
import re
import shutil
import subprocess
import sys
import tempfile
import typing
from pathlib import Path
from types import ModuleType

import pytest
from helpers.m12_04_pr2 import INTERNAL_REFS
from hypothesis import given, settings
from hypothesis import strategies as st

ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "scripts" / "dev" / "sync_report_snippets.py"
DOCS = ROOT / "docs" / "methods"
SNIPPETS = ROOT / "src" / "gippa" / "report" / "snippets"
END = "<!-- report:end -->"
BLOCK = re.compile(r"<!-- report:begin(?: ([A-Za-z0-9_-]+))? -->(.*?)<!-- report:end -->", re.S)

#: the pages PR 3 marks (spec § PR 3 behaviour 1)
PR3_PAGES = ("detection-tiers", "placement-fdr", "denovo-precision", "abundance", "operons",
             "variants", "novel-orfs", "function", "determinants")  # fmt: skip
MAX_WORDS = 120
SITE_LEAKS = (re.compile("farmshare", re.I), re.compile("/home/"), re.compile(r"(?<![\w.])local/"),
              re.compile(r"\$GIPPA_ROOT"))  # fmt: skip


def tiers() -> tuple[str, ...]:
    from gippa.preflight.evidence import Tier

    return typing.get_args(Tier)


def blocks(page: str) -> list[tuple[str | None, str]]:
    text = (DOCS / f"{page}.md").read_text(encoding="utf-8")
    return [(m.group(1), m.group(2)) for m in BLOCK.finditer(text)]


def _run(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run([sys.executable, str(SCRIPT), *args], capture_output=True, text=True,
                          cwd=ROOT, check=False)  # fmt: skip


def _module() -> ModuleType:
    spec = importlib.util.spec_from_file_location("sync_report_snippets_pr3", SCRIPT)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ---------------------------------------------------------------- the docs (AC 3.1)


@pytest.mark.parametrize("page", PR3_PAGES)
def test_pr3_page_has_a_report_block(page: str) -> None:
    assert blocks(page), f"docs/methods/{page}.md has no <!-- report:begin [key] --> block"


@pytest.mark.parametrize("page", PR3_PAGES)
def test_every_pr3_block_has_at_most_120_words(page: str) -> None:
    for key, body in blocks(page):
        n = len(body.split())
        assert 0 < n <= MAX_WORDS, f"{page}.{key or ''}: {n} words"


def test_islands_block_stays() -> None:
    assert [k for k, _ in blocks("islands")] == [None]
    assert (SNIPPETS / "islands.html").is_file()


def test_detection_tiers_has_one_keyed_block_per_tier() -> None:
    keys = [k for k, _ in blocks("detection-tiers")]
    assert None not in keys, "detection-tiers.md: blocks there are keyed by tier"
    assert sorted(keys) == sorted(tiers()), keys
    assert "decoy_picked" in keys


def test_tier_vocabulary_is_the_four_tiers() -> None:
    """Guard: if the vocabulary grows, detection-tiers.md needs a block for the new tier."""
    assert set(tiers()) == {"decoy_picked", "mayu", "two_peptide", "protein_list"}


@pytest.mark.parametrize("page", ["denovo-precision", "abundance"])
def test_ac32_pages_are_unkeyed(page: str) -> None:
    assert [k for k, _ in blocks(page)] == [None]


def test_snippet_files_for_every_tier_and_page() -> None:
    have = {p.stem for p in SNIPPETS.glob("*.html")}
    for tier in tiers():
        assert f"detection-tiers.{tier}" in have, tier
    for page in PR3_PAGES[1:]:
        assert page in have or any(h.startswith(f"{page}.") for h in have), page
    assert {"denovo-precision", "abundance", "islands"} <= have


def test_repo_snippets_in_sync_with_keyed_blocks() -> None:
    res = _run("--check")
    assert res.returncode == 0, res.stdout + res.stderr
    assert (SNIPPETS / "detection-tiers.mayu.html").is_file()


@pytest.mark.parametrize("path", sorted(SNIPPETS.glob("*.html")), ids=lambda p: p.name)
def test_snippets_hold_no_internal_refs_or_site_paths(path: Path) -> None:
    text = path.read_text(encoding="utf-8")
    for pat in (*INTERNAL_REFS, *SITE_LEAKS):
        assert not pat.search(text), (path.name, pat.pattern)


# ---------------------------------------------------------------- the generator


def _toy(docs: Path) -> None:
    docs.mkdir(parents=True)
    (docs / "toy.md").write_text(
        "# Toy\n\nOutside.\n\n"
        "<!-- report:begin alpha -->\nAlpha text.\n<!-- report:end -->\n\nBetween.\n\n"
        "<!-- report:begin two_peptide -->\nTwo peptide text,\nsecond line.\n\n- item\n"
        f"{END}\n",
        encoding="utf-8",
    )
    (docs / "plain.md").write_text(f"# Plain\n\n<!-- report:begin -->\nPlain text.\n{END}\n",
                                   encoding="utf-8")  # fmt: skip


def test_keyed_blocks_write_topic_key_html(tmp_path: Path) -> None:
    docs, out = tmp_path / "docs", tmp_path / "out"
    _toy(docs)
    res = _run("--docs", str(docs), "--out", str(out))
    assert res.returncode == 0, res.stdout + res.stderr
    assert sorted(p.name for p in out.iterdir()) == [
        "plain.html", "toy.alpha.html", "toy.two_peptide.html"]  # fmt: skip
    assert (out / "toy.alpha.html").read_text("utf-8") == "<p>Alpha text.</p>\n"
    assert (out / "toy.two_peptide.html").read_text("utf-8") == (
        "<p>Two peptide text, second line.</p>\n<ul><li>item</li></ul>\n")  # fmt: skip
    assert (out / "plain.html").read_text("utf-8") == "<p>Plain text.</p>\n"
    assert _run("--check", "--docs", str(docs), "--out", str(out)).returncode == 0


def test_check_names_an_edited_keyed_snippet(tmp_path: Path) -> None:
    docs, out = tmp_path / "docs", tmp_path / "out"
    _toy(docs)
    assert _run("--docs", str(docs), "--out", str(out)).returncode == 0
    (out / "toy.alpha.html").write_text("<p>edited</p>\n", encoding="utf-8")
    res = _run("--check", "--docs", str(docs), "--out", str(out))
    assert res.returncode == 1
    assert "toy.alpha" in res.stdout + res.stderr
    assert "toy.two_peptide" not in res.stdout + res.stderr


def test_check_flags_missing_and_orphan_keyed_snippets(tmp_path: Path) -> None:
    docs, out = tmp_path / "docs", tmp_path / "out"
    _toy(docs)
    assert _run("--docs", str(docs), "--out", str(out)).returncode == 0
    (out / "toy.two_peptide.html").unlink()
    res = _run("--check", "--docs", str(docs), "--out", str(out))
    assert res.returncode == 1 and "toy.two_peptide" in res.stdout
    assert _run("--docs", str(docs), "--out", str(out)).returncode == 0
    (out / "toy.gone.html").write_text("<p>x</p>\n", encoding="utf-8")
    res = _run("--check", "--docs", str(docs), "--out", str(out))
    assert res.returncode == 1 and "toy.gone" in res.stdout
    assert _run("--docs", str(docs), "--out", str(out)).returncode == 0
    assert not (out / "toy.gone.html").exists()


def test_check_on_a_copy_of_the_repo_docs_catches_a_tier_edit(tmp_path: Path) -> None:
    docs, out = tmp_path / "docs", tmp_path / "out"
    shutil.copytree(DOCS, docs)
    shutil.copytree(SNIPPETS, out)
    md = docs / "detection-tiers.md"
    text = md.read_text(encoding="utf-8")
    i = text.index("<!-- report:begin mayu -->")
    j = text.index(END, i)
    md.write_text(text[:j] + "Extra line.\n" + text[j:], encoding="utf-8")
    res = _run("--check", "--docs", str(docs), "--out", str(out))
    assert res.returncode == 1 and "detection-tiers.mayu" in res.stdout


KEYS = st.lists(st.from_regex(r"[a-z][a-z0-9_]{0,10}", fullmatch=True), min_size=1, max_size=5,
                unique=True)  # fmt: skip


@settings(max_examples=25, deadline=None)
@given(keys=KEYS, unkeyed=st.booleans())
def test_property_one_file_per_block_and_check_round_trips(keys: list[str], unkeyed: bool) -> None:
    """Every block gives exactly one file, named by its key; a fresh write passes ``--check``."""
    mod = _module()
    with tempfile.TemporaryDirectory() as tmp:
        docs, out = Path(tmp) / "docs", Path(tmp) / "out"
        docs.mkdir()
        parts = ["# T", ""]
        if unkeyed:
            parts += ["<!-- report:begin -->", "Plain.", END, ""]
        parts += [f"<!-- report:begin {k} -->\nText {k}.\n{END}\n" for k in keys]
        (docs / "topic.md").write_text("\n".join(parts), encoding="utf-8")
        assert mod.main(["--docs", str(docs), "--out", str(out)]) == 0
        want = {f"topic.{k}.html" for k in keys} | ({"topic.html"} if unkeyed else set())
        assert {p.name for p in out.iterdir()} == want
        for k in keys:
            assert (out / f"topic.{k}.html").read_text("utf-8") == f"<p>Text {k}.</p>\n"
        assert mod.main(["--check", "--docs", str(docs), "--out", str(out)]) == 0
