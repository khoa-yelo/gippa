"""M12-04 PR 11a, T-06 (first half): the five ``genes.tsv`` columns whose value is the same on every
row by construction are marked run-level (docs/specs/M12-04-outputs-revision.md § PR 11a).

AC 11a.3: ``detection_tier``, ``abundance_type``, ``abundance_scale``, ``imputation_policy`` and
``annotation_source`` carry ``"gippa:run_level": true`` and a description starting
"Run-level:"; ``docs/reference/tables.md`` renders "run-level" for them (and only for them).

Interface assumed (test-author reading): the generator puts the lower-case word ``run-level`` in
the field's row of the table (any cell); a field without the property never shows it. The
description prefix is capitalised ("Run-level:"), so the lower-case check sees the property, not
the description.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
SCHEMAS = REPO / "src" / "gippa" / "schemas"
GENES_SCHEMA = SCHEMAS / "tables" / "genes.schema.json"
TABLES_MD = REPO / "docs" / "reference" / "tables.md"
RUN_LEVEL = ("detection_tier", "abundance_type", "abundance_scale", "imputation_policy",
             "annotation_source")  # fmt: skip
PROP = "gippa:run_level"


def _fields() -> dict[str, dict]:
    data = json.loads(GENES_SCHEMA.read_text(encoding="utf-8"))
    return {f["name"]: f for f in data["fields"]}


def _gen_reference():  # type: ignore[no-untyped-def]
    sys.path.insert(0, str(REPO / "scripts" / "docs"))
    import gen_reference

    return gen_reference


def _genes_rows(text: str) -> dict[str, str]:
    m = re.search(r"^## `genes`\n(.*?)(?=^## |\Z)", text, re.S | re.M)
    assert m, "no genes section"
    return {r.group(1): r.group(0) for r in re.finditer(r"^\| `([a-z0-9_]+)` \|.*$", m.group(1),
                                                         re.M)}  # fmt: skip


@pytest.mark.parametrize("name", RUN_LEVEL)
def test_run_level_columns_carry_the_property(name: str) -> None:
    assert _fields()[name].get(PROP) is True, _fields()[name]


@pytest.mark.parametrize("name", RUN_LEVEL)
def test_run_level_descriptions_start_with_run_level(name: str) -> None:
    assert str(_fields()[name].get("description", "")).startswith("Run-level:")


def test_no_other_genes_column_is_run_level() -> None:
    marked = {n for n, f in _fields().items() if f.get(PROP)}
    assert marked == set(RUN_LEVEL)


def test_property_is_a_boolean_true_when_present() -> None:
    for name, f in _fields().items():
        if PROP in f:
            assert f[PROP] is True, name


def test_schema_with_the_property_is_still_a_frictionless_schema() -> None:
    frictionless = pytest.importorskip("frictionless")
    data = json.loads(GENES_SCHEMA.read_text(encoding="utf-8"))
    descriptor = {k: v for k, v in data.items() if k != "coordinates"}
    frictionless.Schema.from_descriptor(descriptor)


@pytest.mark.parametrize("name", RUN_LEVEL)
def test_tables_md_renders_run_level(name: str) -> None:
    row = _genes_rows(TABLES_MD.read_text(encoding="utf-8"))[name]
    assert "run-level" in row, row


def test_tables_md_does_not_render_run_level_elsewhere_in_genes() -> None:
    rows = _genes_rows(TABLES_MD.read_text(encoding="utf-8"))
    assert {n for n, r in rows.items() if "run-level" in r} == set(RUN_LEVEL)


def test_generator_renders_the_property_on_a_synthetic_schema(tmp_path: Path) -> None:
    """``render_tables`` reads ``gippa:run_level`` itself (not the description text)."""
    gen = _gen_reference()
    d = tmp_path / "schemas"
    (d / "tables").mkdir(parents=True)
    (d / "datapackage.json").write_text(json.dumps({"resources": []}), encoding="utf-8")
    schema = {
        "name": "toy",
        "fields": [
            {"name": "a", "type": "string", "description": "same on every row", PROP: True},
            {"name": "b", "type": "string", "description": "varies"},
        ],
    }
    (d / "tables" / "toy.schema.json").write_text(json.dumps(schema), encoding="utf-8")
    text = gen.render_tables(d)
    rows = {m.group(1): m.group(0) for m in re.finditer(r"^\| `([a-z])` \|.*$", text, re.M)}
    assert "run-level" in rows["a"], rows["a"]
    assert "run-level" not in rows["b"], rows["b"]


def test_generated_tables_md_is_current() -> None:
    gen = _gen_reference()
    assert gen.render_tables(SCHEMAS) == TABLES_MD.read_text(encoding="utf-8")
