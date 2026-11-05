"""M12-04 PR 11b: reading 1.0 tables (Q3, AC 11b.5) and the user migration note (§ 5.1).

docs/specs/M12-04-outputs-revision.md > PR 11b edge cases and § 5. Q3 (default, binding): the
RunView renames 1.0 columns on read with one migration map, so ``gippa report`` and ``gippa
compare`` (both read through ``RunView.from_outdir``) work on 1.0 output trees. This is not an
output alias: no writer emits an old name.

Interface pinned here (test-author reading, beyond the spec text):

- ``gippa.io.migration`` (the allow-listed module of AC 11b.4) holds ``COLUMN_RENAMES:
  Mapping[str, Mapping[str, str]]``, table (schema stem) -> {1.0 name: 1.1 name}, equal to
  § 5.1, and ``migrate_columns(table: str, columns: Sequence[str]) -> list[str]``: each 1.0 name
  of ``table`` replaced by its 1.1 name, every other column unchanged and in place; a header
  holding both a 1.0 name and its 1.1 name is ambiguous and raises ``ValueError``.
- ``RunView`` applies it to every schema-typed TSV before the schema check, so a 1.0 file
  returns a frame with the 1.1 names and the same values.
- The changelog note is a ``changelog.d/*.removal.md`` fragment (spec § 1 rule 5) whose
  migration table has one line per rename naming the table, the 1.0 name and the 1.1 name in
  backticks (old before new). The hand-written docs section is in ``docs/reference/tables.md``
  or ``docs/user-guide/outputs.md`` (spec § 5: the implementer picks).
"""

from __future__ import annotations

import json
import re
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import pytest
from helpers.m12_04_renames import PATHS, RENAMES, TRIPLES, to_1_0
from hypothesis import given
from hypothesis import strategies as st

REPO = Path(__file__).resolve().parents[2]
TABLES = REPO / "src" / "gippa" / "schemas" / "tables"


def _schema(table: str) -> list[dict[str, Any]]:
    doc = json.loads((TABLES / f"{table}.schema.json").read_text(encoding="utf-8"))
    return list(doc["fields"])


def _names_1_1(table: str) -> list[str]:
    """The 1.1 header (the working tree's schema with the 1.0 names mapped forward)."""
    fwd = RENAMES[table]
    return [fwd.get(f["name"], f["name"]) for f in _schema(table)]


def _names_1_0(table: str) -> list[str]:
    return to_1_0(table, [f["name"] for f in _schema(table)])


# ---------------------------------------------------------------- the map


def test_migration_map_equals_the_spec_rename_table() -> None:
    from gippa.io.migration import COLUMN_RENAMES

    assert {t: dict(m) for t, m in COLUMN_RENAMES.items()} == {
        t: dict(m) for t, m in RENAMES.items()
    }


@pytest.mark.parametrize("table", sorted(RENAMES))
def test_rename_map_is_a_bijection_without_collisions(table: str) -> None:
    """Injective; no 1.1 name is a 1.0 name of the table; no 1.1 name collides with a column
    that keeps its name; so the renamed header has as many distinct names as the old one."""
    renames = RENAMES[table]
    assert len(set(renames.values())) == len(renames)
    assert not set(renames.values()) & set(renames)
    kept = set(_names_1_0(table)) - set(renames)
    assert not set(renames.values()) & kept, (table, set(renames.values()) & kept)
    assert set(renames) <= set(_names_1_0(table))


_EXTRA = st.text(alphabet="abcdefghij_", min_size=1, max_size=8).map(lambda s: f"x_{s}")


@st.composite
def _headers(draw: st.DrawFn) -> tuple[str, list[str]]:
    table = draw(st.sampled_from(sorted(RENAMES)))
    cols = draw(st.permutations(_names_1_0(table)))
    keep = draw(st.lists(st.booleans(), min_size=len(cols), max_size=len(cols)))
    extra = draw(st.lists(_EXTRA, max_size=3, unique=True))
    header = [c for c, k in zip(cols, keep, strict=True) if k] + extra
    return table, draw(st.permutations(header))


@given(_headers())
def test_migrate_columns_renames_in_place_and_is_idempotent(case: tuple[str, list[str]]) -> None:
    from gippa.io.migration import migrate_columns

    table, header = case
    out = migrate_columns(table, header)
    assert out == [RENAMES[table].get(c, c) for c in header]
    assert len(set(out)) == len(out) == len(header)
    assert migrate_columns(table, out) == out
    assert to_1_0(table, out) == header


@given(st.lists(_EXTRA, max_size=6, unique=True), st.sampled_from(["genes", "variants", "x"]))
def test_migrate_columns_leaves_other_tables_alone(header: list[str], table: str) -> None:
    from gippa.io.migration import migrate_columns

    cols = [*header, "class", "flag", "mw"]
    assert migrate_columns(table, cols) == cols


@pytest.mark.parametrize(("table", "old", "new"), TRIPLES, ids=[f"{t}.{o}" for t, o, _ in TRIPLES])
def test_a_header_with_both_names_is_refused(table: str, old: str, new: str) -> None:
    from gippa.io.migration import migrate_columns

    with pytest.raises(ValueError, match=re.escape(old)):
        migrate_columns(table, ["pair_id", old, new])


# ---------------------------------------------------------------- RunView on 1.0 trees

_VALUE = {"integer": "1", "number": "0.5", "boolean": "true"}


def _row(table: str, overrides: Mapping[str, str] | None = None) -> dict[str, str]:
    """One valid cell per 1.1 column, by schema type; ``overrides`` by 1.1 name."""
    fwd = RENAMES[table]
    row = {fwd.get(f["name"], f["name"]): _VALUE.get(str(f.get("type")), "x")
           for f in _schema(table)}  # fmt: skip
    row.update(overrides or {})
    return row


def _write_table(out: Path, table: str, rows: list[dict[str, str]], version: str) -> None:
    header = _names_1_1(table)
    written = header if version == "1.1" else to_1_0(table, header)
    path = out / PATHS[table][0]
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["\t".join(written)] + ["\t".join(r[c] for c in header) for r in rows]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _tree(root: Path, version: str, rows: Mapping[str, list[dict[str, str]]]) -> Path:
    out = root / f"OUT_{version}"
    (out / "tables").mkdir(parents=True)
    for table, rs in rows.items():
        _write_table(out, table, rs, version)
    return out


def _distinct_row(table: str) -> dict[str, str]:
    """A row whose renamed cells hold values that identify them after the read."""
    marks = {"pair_class": "OP", "flags": "internal_promoter_supported"}
    over: dict[str, str] = {}
    for f in _schema(table):
        new = RENAMES[table].get(f["name"], f["name"])
        if new in RENAMES[table].values():
            kind = str(f.get("type"))
            over[new] = marks.get(new, {"integer": "17", "number": "0.25"}.get(kind, "x"))
    return _row(table, over)


@pytest.mark.parametrize("version", ["1.0", "1.1"])
@pytest.mark.parametrize("table", sorted(RENAMES))
def test_runview_reads_both_layouts_with_the_1_1_names(
    tmp_path: Path, table: str, version: str
) -> None:
    from gippa.io.run_view import RunView

    row = _distinct_row(table)
    out = _tree(tmp_path, version, {table: [row]})
    frame = RunView.from_outdir(out).table(PATHS[table][1])
    assert frame is not None
    assert not set(RENAMES[table]) & set(frame.columns), list(frame.columns)
    for new in RENAMES[table].values():
        assert new in frame.columns, (new, list(frame.columns))
        got = frame[new].iloc[0]
        want = row[new]
        if isinstance(got, float):
            assert got == float(want), (new, got)
        else:
            assert str(got) == want, (new, got)


def _cfg() -> Any:
    from gippa.config import GippaConfig

    return GippaConfig.model_validate({"seed": 1})


def _compare_rows() -> dict[str, list[dict[str, str]]]:
    flags = [
        _row("operon_flags", {"pair_id": f"PAIR_{i}", "gene_a": f"G{i}", "gene_b": f"G{i + 1}",
                              "pair_class": "OP", "flags": "internal_promoter_supported",
                              "q_value": q, "lfdr": "0.1"})
        for i, q in ((1, "0.001"), (2, "0.002"), (3, "0.9"))
    ]  # fmt: skip
    orfs = [
        _row("novel_orfs", {"orf_id": f"NOVEL_{i}", "called": c, "tir_ostir": "120.5"})
        for i, c in ((1, "true"), (2, "false"), (3, "true"))
    ]
    return {"operon_flags": flags, "novel_orfs": orfs}


@pytest.mark.parametrize("family", ["operon_flags", "novel_orfs"])
def test_compare_sees_the_same_sets_in_a_1_0_and_a_1_1_run(tmp_path: Path, family: str) -> None:
    """``gippa compare`` of a 1.0 run against a 1.1 run (edge case of PR 11b, Q3)."""
    from gippa.compare.stability import family_set
    from gippa.io.run_view import RunView

    rows = _compare_rows()
    old = RunView.from_outdir(_tree(tmp_path, "1.0", rows))
    new = RunView.from_outdir(_tree(tmp_path, "1.1", rows))
    a, b = family_set(old, family, _cfg()), family_set(new, family, _cfg())
    assert a is not None and b is not None
    assert a == b and len(a) == 2, (a, b)


# ---------------------------------------------------------------- the migration note

_FORBIDDEN = re.compile(r"farmshare|/home/|local/|\$GIPPA_ROOT", re.IGNORECASE)


def _pair_line(text: str, table: str, old: str, new: str) -> bool:
    for line in text.splitlines():
        i, j = line.find(f"`{old}`"), line.find(f"`{new}`")
        if table in line and 0 <= i < j:
            return True
    return False


def _removal_text() -> str:
    frags = sorted((REPO / "changelog.d").glob("*.removal.md"))
    return "\n".join(p.read_text(encoding="utf-8") for p in frags)


@pytest.mark.parametrize(("table", "old", "new"), TRIPLES, ids=[f"{t}.{o}" for t, o, _ in TRIPLES])
def test_changelog_removal_fragment_lists_the_rename(table: str, old: str, new: str) -> None:
    text = _removal_text()
    assert text, "no changelog.d/*.removal.md fragment"
    assert _pair_line(text, table, old, new), f"{table}: `{old}` -> `{new}` not in the note"


def test_changelog_note_says_no_aliases_and_names_the_pairs_join() -> None:
    """§ 5 wording: old names are not kept as aliases; § PR 11b: a pairs/operon_flags join on
    pair_id needs suffixes now that both have pair_class and flags."""
    text = _removal_text()
    assert re.search(r"not kept as aliases|no aliases", text, re.IGNORECASE), text
    assert "suffix" in text.lower(), text
    assert not _FORBIDDEN.search(text), _FORBIDDEN.search(text)


@pytest.mark.parametrize(("table", "old", "new"), TRIPLES, ids=[f"{t}.{o}" for t, o, _ in TRIPLES])
def test_docs_migration_section_lists_the_rename(table: str, old: str, new: str) -> None:
    docs = [REPO / "docs" / "reference" / "tables.md", REPO / "docs" / "user-guide" / "outputs.md"]
    text = "\n".join(p.read_text(encoding="utf-8") for p in docs if p.is_file())
    assert _pair_line(text, table, old, new), f"{table}: `{old}` -> `{new}` not in the docs"
