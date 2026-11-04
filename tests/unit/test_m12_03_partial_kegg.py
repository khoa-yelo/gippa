"""M12-03 lead decisions 3-5 (implementer tests): partial KEGG, ``--strict`` and provenance.

- **Partial KEGG (decision 3).** KEGG counts as installed only when both KEGG items are ``ok``:
  with one of them, preflight warns ``environment.kegg``; ``functionResources()`` in
  workflow/main.nf then stages neither (workflow/tests/main_function_partial_kegg.nf.test).
- **Level (decision 4).** ``environment.kegg`` is a warning: exit 0, and 3 under ``--strict``.
- **Provenance (decision 5).** Without KEGG the provenance line reads "KEGG not installed", not
  "not recorded in this run".
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from gippa.exitcodes import ExitCode
from gippa.preflight.environment import GIB, EnvProbe, environment_checks, requirements
from gippa.preflight.report import PreflightReport
from gippa.preflight.run import exit_code
from gippa.setupdb.installed import FUNCTION_RESOURCES
from gippa.setupdb.items import SOURCES

REPO = Path(__file__).resolve().parents[2]
SETUP_MANIFEST = REPO / "tests" / "data" / "db_manifest" / "db_manifest.json"
RULE = "environment.kegg"
KEGG_RESOURCES = ("kegg_ko_pathway", "kegg_pathway_names")


def _db(root: Path, drop: tuple[str, ...]) -> Path:
    """A DB directory: the committed setup manifest without ``drop``, every listed function
    resource present on disk."""
    doc = json.loads(SETUP_MANIFEST.read_text(encoding="utf-8"))
    doc["items"] = [i for i in doc["items"] if i["name"] not in set(drop)]
    root.mkdir(parents=True)
    (root / "db_manifest.json").write_text(json.dumps(doc), encoding="utf-8")
    for name in FUNCTION_RESOURCES:
        path = root / SOURCES[name].path
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("", encoding="utf-8")
    return root


@pytest.mark.parametrize("item", KEGG_RESOURCES)
def test_partial_kegg_is_the_kegg_warning(tmp_path: Path, item: str) -> None:
    from gippa.preflight.environment import read_manifest

    found = read_manifest(str(_db(tmp_path / "db", (item,))))
    assert found is not None and found.get("kegg") is None, found
    reqs = requirements(["core", "function"], "targeted")
    probe = EnvProbe(found=found, free_disk_bytes=500 * GIB, total_ram_bytes=64 * GIB)
    checks = {c.rule_id: c for c in environment_checks(reqs, probe, 0, 0)}
    assert checks["environment.tools"].level == "ok", checks["environment.tools"]
    assert checks[RULE].level == "warn", checks[RULE]


def test_kegg_requirement_is_optional_go_and_cog_are_not() -> None:
    by_name = {r.name: r for r in requirements(["core", "function"], "targeted")}
    assert by_name["kegg"].optional
    assert not by_name["go"].optional and not by_name["cog"].optional
    assert not any(r.optional for r in requirements(["core", "mobile", "bgc"], "genome"))


@pytest.mark.parametrize(("strict", "code"), [(False, ExitCode.SUCCESS),
                                              (True, ExitCode.WARNINGS_STRICT)])  # fmt: skip
def test_kegg_warning_exits_3_only_under_strict(strict: bool, code: ExitCode) -> None:
    reqs = requirements(["core", "function"], "targeted")
    found: dict[str, str | None] = {r.name: "1.0" for r in reqs if r.name != "kegg"}
    found["bakta_db"] = "6.0"
    probe = EnvProbe(found=found, free_disk_bytes=500 * GIB, total_ram_bytes=64 * GIB)
    checks = environment_checks(reqs, probe, 0, 0)
    assert [c.rule_id for c in checks if c.level == "warn"] == [RULE]
    report = PreflightReport("0", [], checks, {}, "", [], [])
    assert exit_code(report, strict=strict) == code


# ---------------------------------------------------------------- provenance


def _provenance_lines(tmp_path: Path, flags: list[str], *, kegg_db: bool) -> dict[str, str]:
    from helpers.report_content import manifest, write_flags, write_gsea
    from helpers.report_fixtures import run_out

    from gippa.io.run_view import RunView
    from gippa.report.content import provenance

    out = run_out(tmp_path)
    write_gsea(out, [None])
    if flags:
        write_flags(out, "function", flags)
    doc = manifest()
    if not kegg_db:
        doc["databases"] = [d for d in doc["databases"] if not str(d["name"]).startswith("kegg")]
    prov = provenance(RunView.from_outdir(out), doc, None)
    return {str(k): str(v) for k, v in prov["lines"]}


def test_provenance_says_kegg_not_installed(tmp_path: Path) -> None:
    lines = _provenance_lines(tmp_path, ["layer_missing:kegg"], kegg_db=False)
    assert lines["KEGG retrieval date"] == "KEGG not installed", lines


def test_provenance_without_the_flag_still_says_not_recorded(tmp_path: Path) -> None:
    lines = _provenance_lines(tmp_path, [], kegg_db=False)
    assert lines["KEGG retrieval date"] == "not recorded in this run", lines


def test_provenance_with_kegg_keeps_the_date(tmp_path: Path) -> None:
    from helpers.report_content import KEGG_DATE

    lines = _provenance_lines(tmp_path, [], kegg_db=True)
    assert lines["KEGG retrieval date"] == KEGG_DATE, lines


# ---------------------------------------------------------------- #541 review M1, M2


@pytest.mark.parametrize("item", KEGG_RESOURCES)
def test_lone_kegg_item_is_not_staged(tmp_path: Path, item: str) -> None:
    """M2: a lone KEGG item is never staged by the local executor's resource lookup."""
    from gippa.setupdb.installed import function_resources

    other = next(k for k in KEGG_RESOURCES if k != item)
    got = function_resources(_db(tmp_path / "db", (other,)))
    assert not set(KEGG_RESOURCES) & set(got), got
    assert {"go_basic", "goslim_generic", "cog2020_fun"} <= set(got), got


def test_both_kegg_items_are_staged(tmp_path: Path) -> None:
    from gippa.setupdb.installed import function_resources

    assert set(function_resources(_db(tmp_path / "db", ()))) == set(FUNCTION_RESOURCES)


def test_provenance_partial_kegg_says_not_installed(tmp_path: Path) -> None:
    """M1: with the flag, a lone KEGG item's retrieval date is not shown."""
    lines = _provenance_lines(tmp_path, ["layer_missing:kegg"], kegg_db=True)
    assert lines["KEGG retrieval date"] == "KEGG not installed", lines


@pytest.mark.parametrize("items", [(), ("kegg_ko_pathway",), ("kegg_pathway_names",),
                                   ("kegg_ko_pathway", "kegg_pathway_names")])  # fmt: skip
def test_kegg_licence_note_only_with_both_items(items: tuple[str, ...]) -> None:
    """M1: the KEGG academic-use note is recorded only when both KEGG items are installed."""
    from types import SimpleNamespace
    from typing import Any

    from gippa.manifest import KEGG_NOTE, licence_notes

    db = {"items": [{"name": n, "status": "ok"} for n in ("go_basic", "cog2020_fun", *items)]}
    src: Any = SimpleNamespace(cfg=SimpleNamespace(modules=["core", "function"]),
                               tool_reports=[], db_manifest=db, operon_meta=None)  # fmt: skip
    notes = licence_notes(src, {"function"})
    assert (KEGG_NOTE in notes) == (len(items) == 2), notes
