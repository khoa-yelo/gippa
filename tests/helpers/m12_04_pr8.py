"""Fixtures of the M12-04 PR 8 tests: one notice when no function term was assigned (F-01;
docs/specs/M12-04-outputs-revision.md § PR 8).

Every fixture is the PR 1 fixture P (``tests/unit/report/test_m12_04_pr1_fixtures.py``: an
A-like peptides + iBAQ run with the modules ``core, bgc, function``) plus what the ``function``
step leaves behind, written by hand (R11: nothing row-level from the example bundle):

- ``logs/steps/function.events.jsonl``: ``function.run_flags`` and, unless ``collections`` is
  ``None`` (an old run), ``function.collections`` as ``gippa.steps.function`` logs it:
  collection -> ``{"n_terms", "coverage", "n_unknown_ids"}``;
- ``tables/category_summary.tsv``: header-only (``category=HEADER``, as run A), one row per
  collection with its coverage (a mapping), or absent (``None``);
- ``tables/enrichment_gsea.tsv`` and ``enrichment_detection.tsv``: header-only (no term, no set);
- the run flag ``bakta_db_light`` (``logs/steps/annotate.events.jsonl``) when ``light``;
- ``layer_missing:kegg`` among the function run flags when the event has no ``kegg`` key
  (#543: ``gippa setup --no-kegg``), unless ``kegg_flag`` says otherwise;
- ``genes.tsv`` with COG letters (``cog_category``) on the first ``n_cog`` genes and none on
  the others (``n_cog = 0``: no gene has a letter, as in run A).

``report(run, tmp)`` renders V08 and V14 for real (``render_plots``) and builds the page the way
the ``plots`` step does (``write_report``: its plot files, ``skipped_plots(view)``).
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from unit.report.test_m12_04_pr1_fixtures import Run, fixture_p

from helpers.m12_04_pr2 import words
from helpers.m12_04_pr5 import rewrite_table
from helpers.report_fixtures import Node, Page, parse
from helpers.viz_out_tree import write_events, write_schema_table

GO = ("go_bp", "go_mf", "go_cc", "go_slim")
ALL = ("cog", "kegg", *GO)
HEADER = "header"

#: § PR 8 behaviour: the notice, its first sentence and its two second sentences
FIRST = "No COG, GO or KEGG terms were assigned to any gene."
FIRST_NO_KEGG = "No COG or GO terms were assigned to any gene."
LIGHT = (
    "With the Bakta light database this is expected: rerun with the full database to get "
    "function results."
)
CHECK = "Check the annotation source."
NOTICE_LIGHT = f"{FIRST} {LIGHT}"
NOTICE_CHECK = f"{FIRST} {CHECK}"
#: the plain text of the V08 skip (``SKIP_TEXT``)
NO_TERMS = "no_terms:function"
NO_TERMS_TEXT = "no function terms in this run"
#: the V14 empty state of the examples table
V14_EMPTY = "Nothing to show: no gene has a COG, GO or KEGG term."
#: the function-section KEGG licence note (``gippa.report.content.KEGG_TERMS``)
ACADEMIC = "academic use only"
RETRIEVAL = "KEGG retrieval date"
NOT_RECORDED = "not recorded"


def zero(*colls: str) -> dict[str, tuple[float, int]]:
    """Every collection in ``colls`` (default: all six) at coverage 0 with 0 terms (run A)."""
    return {c: (0.0, 0) for c in (colls or ALL)}


def partial(cog: float = 0.72, go: float = 0.41, kegg: float | None = 0.66
            ) -> dict[str, tuple[float, int]]:  # fmt: skip
    """COG, every GO collection and (unless ``None``, not installed) KEGG with terms."""
    out = {"cog": (cog, 9)} | {g: (go, 16) for g in GO}
    if kegg is not None:
        out["kegg"] = (kegg, 5)
    return out


def collections_event(colls: Mapping[str, tuple[float | None, int]]) -> dict[str, Any]:
    return {c: {"coverage": cov, "n_terms": n, "n_unknown_ids": 0}
            for c, (cov, n) in colls.items()}  # fmt: skip


def write_function(
    out: Path,
    *,
    collections: Mapping[str, tuple[float | None, int]] | None,
    category: Mapping[str, float] | str | None = HEADER,
    light: bool = True,
    kegg_flag: bool | None = None,
) -> None:
    """The ``function`` step's events and tables (see the module docstring)."""
    no_kegg = collections is not None and "kegg" not in collections
    flags = ["layer_missing:kegg"] if (no_kegg if kegg_flag is None else kegg_flag) else []
    records: list[dict[str, Any]] = [
        {"event": "function.run_flags", "level": "INFO", "step": "function",
         "data": {"flags": flags}},
    ]  # fmt: skip
    if collections is not None:
        records.append({"event": "function.collections", "level": "INFO", "step": "function",
                        "data": collections_event(collections)})  # fmt: skip
    write_events(out / "logs" / "steps" / "function.events.jsonl", records)
    if light:
        write_events(out / "logs" / "steps" / "annotate.events.jsonl",
                     [{"event": "annotate.run_flags", "level": "INFO", "step": "annotate",
                       "data": {"flags": ["bakta_db_light"]}}])  # fmt: skip
    tables = out / "tables"
    if category == HEADER:
        write_schema_table(tables / "category_summary.tsv", "category_summary", [])
    elif isinstance(category, Mapping):
        rows = [{"collection": c, "category": f"{c}_1", "name": f"{c} one", "n_genes": 10,
                 "n_detected": 5, "expected": 5.0, "o_minus_e": 0.0, "coverage": v}
                for c, v in category.items()]  # fmt: skip
        write_schema_table(tables / "category_summary.tsv", "category_summary", rows)
    write_schema_table(tables / "enrichment_gsea.tsv", "enrichment_gsea", [])
    write_schema_table(tables / "enrichment_detection.tsv", "enrichment_detection", [])


def set_cog(out: Path, n_cog: int) -> None:
    """COG letters on the first ``n_cog`` genes (``J``, ``C``, ``E`` in turn), none after."""
    letters = "JCE"
    count = {"i": 0}

    def fn(r: dict[str, str]) -> dict[str, str]:
        i = count["i"]
        count["i"] += 1
        r["cog_category"] = letters[i % 3] if i < n_cog else "NA"
        return r

    rewrite_table(out, "tables/genes.tsv", fn)


def function_run(
    tmp: Path,
    name: str = "F",
    *,
    collections: Mapping[str, tuple[float | None, int]] | None,
    category: Mapping[str, float] | str | None = HEADER,
    light: bool = True,
    kegg_flag: bool | None = None,
    n_cog: int = 0,
) -> Run:
    """Fixture P with the function outputs of ``write_function``."""
    run = fixture_p(tmp, name)
    write_function(run.out, collections=collections, category=category, light=light,
                   kegg_flag=kegg_flag)  # fmt: skip
    set_cog(run.out, n_cog)
    return run


# ---------------------------------------------------------------- rendering


def render(run: Run, pids: Sequence[str] = ("V08", "V14")) -> dict[str, Any]:
    """``pid -> list[Instance] | Skip`` of ``pids`` on ``run`` (the snapshot style)."""
    from gippa.steps.plots import load_plot_modules
    from gippa.viz import base
    from helpers.snapshot_fixture import styled

    def go() -> dict[str, Any]:
        load_plot_modules()
        view, cfg = run.view(), run.cfg()
        return {pid: base.registered()[pid].fn(view, cfg) for pid in pids}

    result: dict[str, Any] = styled(go)
    return result


def report(run: Run, tmp: Path) -> tuple[Page, Path, list[Any]]:
    """The page the ``plots`` step writes for ``run`` (V08 and V14 rendered for real), its
    ``plots/`` directory and the plot events."""
    from gippa.steps.plots import load_plot_modules, render_plots, write_report

    load_plot_modules()
    view, cfg = run.view(), run.cfg()
    plots = tmp / f"plots-{run.name}"
    events = render_plots(view, cfg, plots, ["V08", "V14"])
    page = write_report(view, plots, tmp / f"report-{run.name}", cfg)
    return parse(page.read_text(encoding="utf-8")), plots, events


def function_section(page: Page) -> Node:
    return page.section("function")


def notices(node: Node) -> list[str]:
    return [words(n.text()) for n in node.by_class("notice")]


def section_text(page: Page) -> str:
    return words(function_section(page).text())
