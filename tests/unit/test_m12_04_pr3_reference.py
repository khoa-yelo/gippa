"""M12-04 PR 3, behaviour 2 (D-01) and AC 3.3, 3.4: every table in the reference page has a
path, descriptions lose their "(GIPPA_SPEC.md §…)" parentheticals, and the user guide no longer
says "GIPPA 1.0 builds" (docs/specs/M12-04-outputs-revision.md § PR 3).

Interface pinned here (PR 3 test author, beyond the spec text):
- ``render_tables(schemas_dir)`` keeps its signature. The Path line of a table is, in order:
  the datapackage resource ``path``; for ``compare_*`` tables
  ``<compare outdir>/tables/<name>.tsv``; else ``gippa.io.run_view.PATHS[<name>]``, with the
  ``genes`` schema read as the ``genes_tsv`` key (``tables/genes.tsv``);
- a parenthetical that starts with ``GIPPA_SPEC.md`` is dropped from table and field
  descriptions together with the space before it; other parentheticals stay; the schema files
  are not edited.
"""

from __future__ import annotations

import importlib.util
import json
import re
import sys
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest

REPO = Path(__file__).resolve().parents[1].parent
SCRIPT = REPO / "scripts" / "docs" / "gen_reference.py"
SCHEMAS = REPO / "src" / "gippa" / "schemas"
TABLES_MD = REPO / "docs" / "reference" / "tables.md"
USER_GUIDE = REPO / "docs" / "user-guide"
COMPARE_DIR = "<compare outdir>/tables"


@pytest.fixture(scope="module")
def gr() -> ModuleType:
    spec = importlib.util.spec_from_file_location("gen_reference_pr3", SCRIPT)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    sys.modules["gen_reference_pr3"] = mod
    spec.loader.exec_module(mod)
    return mod


def _sections(md: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for chunk in re.split(r"^## ", md, flags=re.M)[1:]:
        head, _, body = chunk.partition("\n")
        out[head.strip().strip("`")] = body
    return out


def _path_line(body: str) -> str:
    lines = [ln for ln in body.splitlines() if ln.startswith("- Path:")]
    assert len(lines) == 1, lines
    return lines[0]


def _expected_path(name: str) -> str:
    from gippa.io.run_view import PATHS

    package = json.loads((SCHEMAS / "datapackage.json").read_text(encoding="utf-8"))
    for res in package["resources"]:
        if res["name"] == name:
            return str(res["path"])
    if name.startswith("compare_"):
        return f"{COMPARE_DIR}/{name}.tsv"
    return PATHS["genes_tsv" if name == "genes" else name]


def _table_names() -> list[str]:
    return sorted(json.loads(p.read_text(encoding="utf-8"))["name"]
                  for p in (SCHEMAS / "tables").glob("*.schema.json"))  # fmt: skip


# ---------------------------------------------------------------- AC 3.3 on the committed page


def test_ac33_tables_md_has_no_missing_path() -> None:
    text = TABLES_MD.read_text(encoding="utf-8")
    assert "Path: —" not in text, [ln for ln in text.splitlines() if "Path: —" in ln][:5]


def test_ac33_tables_md_has_no_spec_reference() -> None:
    text = TABLES_MD.read_text(encoding="utf-8")
    assert "GIPPA_SPEC.md" not in text


def test_example_genes_path() -> None:
    body = _sections(TABLES_MD.read_text(encoding="utf-8"))["genes"]
    assert _path_line(body) == "- Path: `tables/genes.tsv`"


@pytest.mark.parametrize("name", _table_names())
def test_every_table_has_its_path(name: str) -> None:
    body = _sections(TABLES_MD.read_text(encoding="utf-8"))[name]
    assert _path_line(body) == f"- Path: `{_expected_path(name)}`"


@pytest.mark.parametrize("name", ["compare_abundance", "compare_detection", "compare_stability"])
def test_compare_tables_live_under_the_compare_outdir(name: str) -> None:
    body = _sections(TABLES_MD.read_text(encoding="utf-8"))[name]
    assert _path_line(body) == f"- Path: `<compare outdir>/tables/{name}.tsv`"


def test_generated_tables_page_matches_the_committed_one(gr: ModuleType) -> None:
    assert gr.render_tables(SCHEMAS) == TABLES_MD.read_text(encoding="utf-8")


# ---------------------------------------------------------------- the generator on toy schemas


def _write(path: Path, obj: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def _toy(root: Path) -> Path:
    _write(root / "datapackage.json", {"name": "toy", "resources": [
        {"name": "islands", "path": "custom/islands.tsv", "schema": "tables/islands.schema.json"},
    ]})  # fmt: skip
    _write(
        root / "tables" / "islands.schema.json",
        {
            "name": "islands",
            "description": "One row per island.",
            "fields": [{"name": "island_id", "type": "string", "description": "Island ID"}],
        },
    )
    _write(root / "tables" / "genes.schema.json", {
        "name": "genes",
        "description": "One row per gene (GIPPA_SPEC.md §9.3; M6-07, M6-08).",
        "fields": [
            {"name": "locus_tag", "type": "string",
             "description": "Locus tag (GIPPA_SPEC.md §9.3) of the gene."},
            {"name": "mw_da", "type": "number", "description": "Molecular weight (Da)."},
        ]})  # fmt: skip
    _write(root / "tables" / "compare_toy.schema.json", {
        "name": "compare_toy", "description": "Toy comparison (GIPPA_SPEC.md §7.10).",
        "fields": [{"name": "x", "type": "string", "description": "X"}]})  # fmt: skip
    _write(root / "tables" / "variants.schema.json", {
        "name": "variants", "description": "Variants.",
        "fields": [{"name": "x", "type": "string", "description": "X"}]})  # fmt: skip
    return root


def test_toy_datapackage_path_wins(gr: ModuleType, tmp_path: Path) -> None:
    secs = _sections(gr.render_tables(_toy(tmp_path)))
    assert _path_line(secs["islands"]) == "- Path: `custom/islands.tsv`"


def test_toy_run_view_paths_fill_unlisted_tables(gr: ModuleType, tmp_path: Path) -> None:
    secs = _sections(gr.render_tables(_toy(tmp_path)))
    assert _path_line(secs["genes"]) == "- Path: `tables/genes.tsv`"
    assert _path_line(secs["variants"]) == "- Path: `tables/variants.tsv`"


def test_toy_compare_table_path(gr: ModuleType, tmp_path: Path) -> None:
    secs = _sections(gr.render_tables(_toy(tmp_path)))
    assert _path_line(secs["compare_toy"]) == "- Path: `<compare outdir>/tables/compare_toy.tsv`"


def test_toy_spec_parentheticals_removed(gr: ModuleType, tmp_path: Path) -> None:
    text = gr.render_tables(_toy(tmp_path))
    assert "GIPPA_SPEC" not in text
    secs = _sections(text)
    assert "One row per gene." in secs["genes"]
    assert "Toy comparison." in secs["compare_toy"]
    assert "| Locus tag of the gene. |" in secs["genes"]
    assert "Molecular weight (Da)." in secs["genes"], "other parentheticals stay"


def test_schema_files_untouched_by_rendering(gr: ModuleType, tmp_path: Path) -> None:
    root = _toy(tmp_path)
    before = (root / "tables" / "genes.schema.json").read_bytes()
    gr.render_tables(root)
    assert (root / "tables" / "genes.schema.json").read_bytes() == before


# ---------------------------------------------------------------- AC 3.4


def test_ac34_no_user_guide_page_says_gippa_1_0_builds() -> None:
    hits = [p.relative_to(REPO).as_posix() for p in sorted(USER_GUIDE.rglob("*.md"))
            if re.search(r"GIPPA 1\.0 builds", p.read_text(encoding="utf-8"))]  # fmt: skip
    assert hits == []


def test_plots_page_says_gippa_builds_16_plots() -> None:
    text = (USER_GUIDE / "plots-and-report.md").read_text(encoding="utf-8")
    assert "GIPPA builds 16 plots" in text
