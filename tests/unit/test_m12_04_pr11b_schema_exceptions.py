"""M12-04 PR 11b: the ADR-0075 schema-diff exception file (AC 11b.2, 11b.3).

docs/specs/M12-04-outputs-revision.md > PR 11b "Schema-diff exception":

- ``scripts/release/schema_exceptions.json`` maps ``{"v1.0.0": {"adr": "ADR-0075", "renamed":
  {"<table>": {"<old>": "<new>", ...}}, "retyped": {...}, "reordered": ["genes"]}}``;
- ``schema_diff.check`` accepts a removed column without a major bump only if the column is in
  ``renamed`` for the tag compared against, its new name exists, and ``docs/adr/ADR-0075-*.md``
  exists;
- ``retyped`` and ``reordered`` are documentation: the check reads names only.

Test-author readings (beyond the spec text):

- ``check(root)`` reads the exception file and the ADR from ``root``'s working tree (the file is
  written after the tag it excuses, so it cannot come from the tag). A missing exception file
  means no exceptions. The ADR is found as ``docs/adr/<adr>-*.md``, ``<adr>`` taken from the
  entry.
- "its new name exists" means in the same table of the working tree.
- The exception excuses the removal only: the new name is still an added column and still needs
  the ``output_schema_version`` minor bump.
- The spec names the new test file ``test_m12_04_schema_exceptions.py``; the lead's brief fixes
  ``test_m12_04_pr11b_*.py``, used here.
"""

from __future__ import annotations

import importlib.util
import json
import os
import subprocess
import sys
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest
from helpers.m12_04_renames import RENAMES, TRIPLES

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "release" / "schema_diff.py"
EXCEPTIONS = REPO / "scripts" / "release" / "schema_exceptions.json"
EXC_REL = "scripts/release/schema_exceptions.json"
TABLES_REL = "src/gippa/schemas/tables"
COPIED = (
    "src/gippa/__init__.py",
    "src/gippa/schemas/run_manifest.schema.json",
    "src/gippa/schemas/datapackage.json",
)

_GIT_ENV = {
    "GIT_AUTHOR_NAME": "Fixture",
    "GIT_AUTHOR_EMAIL": "fixture@example.invalid",
    "GIT_COMMITTER_NAME": "Fixture",
    "GIT_COMMITTER_EMAIL": "fixture@example.invalid",
    "GIT_CONFIG_GLOBAL": os.devnull,
    "GIT_CONFIG_SYSTEM": os.devnull,
}

PAIRS_1_0 = ["pair_id", "class", "gene_a", "promoter_between", "flags"]
PAIRS_1_1 = ["pair_id", "pair_class", "gene_a", "promoter_score_between", "flags"]
GENES = ["locus_tag", "contig", "start", "end"]
ENTRY = {"pairs": {"class": "pair_class", "promoter_between": "promoter_score_between"}}


def _load() -> ModuleType:
    spec = importlib.util.spec_from_file_location("gippa_release_schema_diff_m1204", SCRIPT)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def git(repo: Path, *args: str) -> str:
    out = subprocess.run(
        ["git", "-C", str(repo), *args],
        env={**os.environ, **_GIT_ENV},
        check=True,
        capture_output=True,
        text=True,
    )
    return out.stdout.strip()


def _write(root: Path, rel: str, text: str) -> None:
    path = root / rel
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def write_tree(root: Path, *, version: str, output: str, tables: dict[str, list[str]]) -> None:
    _write(root, "src/gippa/__init__.py", f'"""GIPPA."""\n\n__version__ = "{version}"\n')
    tdir = root / TABLES_REL
    if tdir.is_dir():
        for old in tdir.glob("*.schema.json"):
            old.unlink()
    for stem, cols in tables.items():
        doc = {"name": stem, "fields": [{"name": c, "type": "string"} for c in cols]}
        _write(root, f"{TABLES_REL}/{stem}.schema.json", json.dumps(doc, indent=2))
    props = {"schema_version": {"type": "string"}}
    _write(root, "src/gippa/schemas/run_manifest.schema.json", json.dumps({"properties": props}))
    package = {"name": "gippa-run", "output_schema_version": output, "resources": []}
    _write(root, "src/gippa/schemas/datapackage.json", json.dumps(package))


def write_exceptions(root: Path, doc: dict[str, Any] | None, *, adr_file: bool = True) -> None:
    if doc is not None:
        _write(root, EXC_REL, json.dumps(doc, indent=2))
    if adr_file:
        _write(root, "docs/adr/ADR-0075-outputs-revision-1-1.md", "# ADR-0075\n")


def commit(root: Path, msg: str, tag: str | None = None) -> None:
    git(root, "add", "-A")
    git(root, "commit", "-q", "--allow-empty", "-m", msg)
    if tag:
        git(root, "tag", tag)


def released(tmp: Path) -> Path:
    root = tmp / "repo"
    root.mkdir()
    git(root, "init", "-q", "-b", "main")
    write_tree(root, version="1.0.0", output="1.2", tables={"genes": GENES, "pairs": PAIRS_1_0})
    commit(root, "release 1.0.0", tag="v1.0.0")
    return root


def head(
    root: Path,
    doc: dict[str, Any] | None,
    *,
    adr_file: bool = True,
    pairs: list[str] | None = None,
    genes: list[str] | None = None,
    output: str = "1.3",
) -> list[str]:
    tables = {
        "genes": GENES if genes is None else genes,
        "pairs": PAIRS_1_1 if pairs is None else pairs,
    }
    write_tree(root, version="1.1.0", output=output, tables=tables)
    write_exceptions(root, doc, adr_file=adr_file)
    commit(root, "1.1 renames")
    return [str(p) for p in _load().check(root)]


def entry(renamed: dict[str, dict[str, str]] | None = None, **kw: Any) -> dict[str, Any]:
    return {"adr": "ADR-0075", "renamed": ENTRY if renamed is None else renamed,
            "retyped": {}, "reordered": ["genes"], **kw}  # fmt: skip


def _about(problems: list[str], table: str, col: str) -> bool:
    return any(table in p and col in p for p in problems)


# ---------------------------------------------------------------- AC 11b.3: fixture repos


def test_listed_renames_pass_without_a_major_bump(tmp_path: Path) -> None:
    root = released(tmp_path)
    assert head(root, {"v1.0.0": entry()}) == []


def test_no_exception_file_still_fails(tmp_path: Path) -> None:
    root = released(tmp_path)
    out = head(root, None)
    assert _about(out, "pairs", "class") and _about(out, "pairs", "promoter_between"), out


def test_an_exception_naming_a_missing_adr_file_fails(tmp_path: Path) -> None:
    root = released(tmp_path)
    out = head(root, {"v1.0.0": entry()}, adr_file=False)
    assert _about(out, "pairs", "class"), out


def test_an_exception_naming_another_adr_id_needs_that_adr_file(tmp_path: Path) -> None:
    """ADR-0075 exists, the entry names ADR-0099, which does not."""
    root = released(tmp_path)
    out = head(root, {"v1.0.0": entry(adr="ADR-0099")})
    assert _about(out, "pairs", "class"), out


def test_an_exception_listed_for_another_tag_does_not_apply(tmp_path: Path) -> None:
    root = released(tmp_path)
    out = head(root, {"v0.9.0": entry()})
    assert _about(out, "pairs", "class") and _about(out, "pairs", "promoter_between"), out


def test_an_unlisted_removed_column_still_fails(tmp_path: Path) -> None:
    """The listed renames pass; the dropped genes.end, not listed, fails alone."""
    root = released(tmp_path)
    out = head(root, {"v1.0.0": entry()}, genes=GENES[:-1])
    assert _about(out, "genes", "end"), out
    assert not any("class" in p or "promoter_between" in p for p in out), out


def test_an_unlisted_rename_in_a_listed_table_still_fails(tmp_path: Path) -> None:
    root = released(tmp_path)
    out = head(root, {"v1.0.0": entry(renamed={"pairs": {"class": "pair_class"}})})
    assert _about(out, "pairs", "promoter_between"), out
    assert not any("column class " in p for p in out), out


def test_a_rename_whose_new_name_is_missing_fails(tmp_path: Path) -> None:
    """``promoter_between`` is listed, but HEAD has no ``promoter_score_between``."""
    root = released(tmp_path)
    out = head(root, {"v1.0.0": entry()}, pairs=["pair_id", "pair_class", "gene_a", "flags"])
    assert _about(out, "pairs", "promoter_between"), out


def test_a_rename_listed_under_another_table_does_not_apply(tmp_path: Path) -> None:
    root = released(tmp_path)
    out = head(root, {"v1.0.0": entry(renamed={"genes": dict(ENTRY["pairs"])})})
    assert _about(out, "pairs", "class"), out


def test_retyped_and_reordered_do_not_excuse_a_removal(tmp_path: Path) -> None:
    root = released(tmp_path)
    doc = {"v1.0.0": entry(renamed={}, retyped={"pairs": {"class": "number"}},
                           reordered=["pairs"])}  # fmt: skip
    out = head(root, doc)
    assert _about(out, "pairs", "class"), out


def test_the_new_name_still_needs_the_schema_minor_bump(tmp_path: Path) -> None:
    root = released(tmp_path)
    out = head(root, {"v1.0.0": entry()}, output="1.2")
    assert _about(out, "pairs", "pair_class") and "added" in " ".join(out), out


# ---------------------------------------------------------------- the file in this repository


def _doc() -> dict[str, Any]:
    assert EXCEPTIONS.is_file(), f"{EXCEPTIONS} missing"
    doc = json.loads(EXCEPTIONS.read_text(encoding="utf-8"))
    assert isinstance(doc, dict)
    return doc


def test_exception_file_lists_exactly_the_renames() -> None:
    doc = _doc()
    assert set(doc) == {"v1.0.0"}, sorted(doc)
    e = doc["v1.0.0"]
    assert e["adr"] == "ADR-0075"
    assert e["renamed"] == {t: dict(m) for t, m in RENAMES.items()}
    assert isinstance(e["retyped"], dict)
    assert e["reordered"] == ["genes"]
    assert list((REPO / "docs" / "adr").glob("ADR-0075-*.md")), "ADR-0075 file missing"


def _tag_text(rel: str) -> str | None:
    res = subprocess.run(["git", "-C", str(REPO), "show", f"v1.0.0:{rel}"],
                         capture_output=True, text=True, check=False)  # fmt: skip
    return res.stdout if res.returncode == 0 else None


@pytest.fixture
def real_release(tmp_path: Path) -> Path:
    """A fixture repo whose ``v1.0.0`` holds this repository's v1.0.0 schema files and whose
    HEAD holds the working tree's schema files, the real exception file and ADR-0075."""
    if _tag_text("src/gippa/__init__.py") is None:
        pytest.skip("no v1.0.0 tag in this clone")
    root = tmp_path / "real"
    root.mkdir()
    git(root, "init", "-q", "-b", "main")
    listing = subprocess.run(["git", "-C", str(REPO), "ls-tree", "--name-only",
                              f"v1.0.0:{TABLES_REL}"], capture_output=True, text=True,
                             check=True).stdout.split()  # fmt: skip
    for rel in (*COPIED, *(f"{TABLES_REL}/{n}" for n in listing if n.endswith(".json"))):
        text = _tag_text(rel)
        assert text is not None, rel
        _write(root, rel, text)
    commit(root, "release 1.0.0", tag="v1.0.0")
    for old in (root / TABLES_REL).glob("*.schema.json"):
        old.unlink()
    for path in (REPO / TABLES_REL).glob("*.schema.json"):
        _write(root, f"{TABLES_REL}/{path.name}", path.read_text(encoding="utf-8"))
    for rel in COPIED:
        _write(root, rel, (REPO / rel).read_text(encoding="utf-8"))
    for adr in (REPO / "docs" / "adr").glob("ADR-0075-*.md"):
        _write(root, f"docs/adr/{adr.name}", adr.read_text(encoding="utf-8"))
    return root


def _run(root: Path, doc: dict[str, Any]) -> list[str]:
    _write(root, EXC_REL, json.dumps(doc, indent=2))
    commit(root, "head")
    return [str(p) for p in _load().check(root)]


def test_the_real_exception_file_makes_the_real_schemas_pass(real_release: Path) -> None:
    assert _run(real_release, _doc()) == []


@pytest.mark.parametrize(("table", "old", "new"), TRIPLES, ids=[f"{t}.{o}" for t, o, _ in TRIPLES])
def test_the_real_file_minus_one_rename_fails_for_that_column(
    real_release: Path, table: str, old: str, new: str
) -> None:
    doc = _doc()
    renamed = {t: dict(m) for t, m in doc["v1.0.0"]["renamed"].items()}
    del renamed[table][old]
    doc["v1.0.0"]["renamed"] = renamed
    out = _run(real_release, doc)
    assert any(f"{table}: column {old} removed" in p for p in out), out
    assert len(out) == 1, out


def test_this_repository_has_the_renames_and_passes_the_release_gate() -> None:
    """AC 11b.2: the M11-14h repository check passes through the exception file, and is not
    passing merely because nothing was renamed."""
    for table, renames in RENAMES.items():
        doc = json.loads((REPO / TABLES_REL / f"{table}.schema.json").read_text("utf-8"))
        names = {f["name"] for f in doc["fields"]}
        assert set(renames.values()) <= names, table
    mod = _load()
    try:
        problems = mod.check(REPO)
    except mod.NoPreviousRelease as exc:
        pytest.skip(str(exc))
    assert problems == []
