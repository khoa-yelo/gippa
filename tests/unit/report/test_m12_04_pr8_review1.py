"""M12-04 PR 8 review round 1 (#561): which ``function.collections`` event counts, malformed
events, tiny coverages and the collections named in the notice
(docs/specs/M12-04-outputs-revision.md § PR 8).

Rules pinned here (review r1, beyond the PR 8 tests):
1. the run's coverage is the **last** ``function.collections`` event with at least one
   well-formed collection (a mapping holding ``coverage`` or ``n_terms``);
2. an event with no well-formed collection is ignored; when no event qualifies, coverage is
   unknown: no notice, V08 drawn, the coverage note says "not recorded" (header-only
   ``category_summary``), and nothing raises;
3. a coverage strictly between 0 and 1 never reads "0 %" or "100 %" in the coverage note;
4. the notice names the collections whose keys the event holds: COG (``cog``), GO (any
   ``go_*``), KEGG (``kegg``), joined "A or B" / "A, B or C".
"""

from __future__ import annotations

import re
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import pytest
from helpers.m12_04_pr2 import words
from helpers.m12_04_pr8 import (
    FIRST,
    GO,
    LIGHT,
    NOT_RECORDED,
    collections_event,
    function_run,
    function_section,
    notices,
    partial,
    render,
    report,
    zero,
)
from helpers.report_fixtures import Page
from helpers.viz_out_tree import write_events
from unit.report.test_m12_04_pr1_fixtures import Run

from gippa.viz import base

FIRST_COG_KEGG = "No COG or KEGG terms were assigned to any gene."
COVERAGE_PARA = "coverage per collection"
#: the F-01 notice (not V14's own empty-state notice)
NOTICE_START = re.compile(r"No (COG|GO|KEGG)[A-Za-z, ]* terms were assigned to any gene\.")


def _events(run: Run, datas: Sequence[Mapping[str, Any]]) -> None:
    """Rewrite the function events: the run flags, then one ``function.collections`` event per
    item of ``datas`` (written verbatim), in order."""
    records: list[dict[str, Any]] = [
        {"event": "function.run_flags", "level": "INFO", "step": "function",
         "data": {"flags": []}},
    ]  # fmt: skip
    records += [{"event": "function.collections", "level": "INFO", "step": "function",
                 "data": dict(d)} for d in datas]  # fmt: skip
    write_events(run.out / "logs" / "steps" / "function.events.jsonl", records)


def _first_sentences(page: Page) -> list[str]:
    return [n for n in notices(function_section(page)) if NOTICE_START.match(n)]


def _coverage_para(page: Page) -> str:
    found = [words(p.text()) for p in function_section(page).by_tag("p")
             if COVERAGE_PARA in words(p.text())]  # fmt: skip
    assert len(found) == 1, [words(p.text()) for p in function_section(page).by_tag("p")]
    return found[0]


def _run(tmp: Path, datas: Sequence[Mapping[str, Any]]) -> Run:
    run = function_run(tmp, collections=zero(), n_cog=60)  # kegg present: no layer_missing:kegg
    _events(run, datas)
    return run


# ---------------------------------------------------------------- 1. last event wins


def test_partial_then_zero_shows_the_notice(tmp_path: Path) -> None:
    run = _run(tmp_path, [collections_event(partial()), collections_event(zero())])
    page, _, _ = report(run, tmp_path)
    found = _first_sentences(page)
    assert len(found) == 1 and found[0].startswith(f"{FIRST} {LIGHT}"), found


def test_zero_then_partial_keeps_the_section(tmp_path: Path) -> None:
    run = _run(tmp_path, [collections_event(zero()), collections_event(partial())])
    page, _, _ = report(run, tmp_path)
    assert not _first_sentences(page), notices(function_section(page))
    note = _coverage_para(page)
    for want in ("COG 72 %", "GO 41 %", "KEGG 66 %"):
        assert want in note, (want, note)


# ---------------------------------------------------------------- 2. malformed events


def test_malformed_last_event_does_not_hide_a_valid_zero_event(tmp_path: Path) -> None:
    run = _run(tmp_path, [collections_event(zero()), {"cog": "x"}])
    page, _, _ = report(run, tmp_path)
    found = _first_sentences(page)
    assert len(found) == 1 and found[0].startswith(FIRST), notices(function_section(page))
    v08 = render(run, ["V08"])["V08"]
    assert isinstance(v08, base.Skip) and v08.reason == "no_terms:function", v08


@pytest.mark.parametrize("data", [{"cog": {}}, {"cog": {"coverage": "0"}}],
                         ids=["empty_mapping", "string_coverage"])  # fmt: skip
def test_only_malformed_events_mean_coverage_unknown(tmp_path: Path, data: dict[str, Any]) -> None:
    run = _run(tmp_path, [data])
    page, _, _ = report(run, tmp_path)
    assert not _first_sentences(page), notices(function_section(page))
    note = _coverage_para(page)
    assert re.search(r"[Cc]overage[^.]{0,60}" + NOT_RECORDED, note), note
    assert "NA" not in note.split(":", 1)[-1], note
    v08 = render(run, ["V08"])["V08"]
    assert not isinstance(v08, base.Skip), v08


# ---------------------------------------------------------------- 3. tiny coverage


@pytest.mark.parametrize("cov", [1e-5, 0.99999])
def test_partial_coverage_never_reads_0_or_100_percent(tmp_path: Path, cov: float) -> None:
    colls = {c: (cov, 3) for c in ("cog", "kegg", *GO)}
    run = _run(tmp_path, [collections_event(colls)])
    page, _, _ = report(run, tmp_path)
    note = _coverage_para(page)
    assert not re.search(r"(?<![\d.])0 %", note), note
    assert "100 %" not in note, note
    for name in ("COG", "GO", "KEGG"):
        assert re.search(rf"\b{name} \d+ %", note), (name, note)


# ---------------------------------------------------------------- 4. names from the keys


def test_notice_names_only_cog_and_kegg_without_go_keys(tmp_path: Path) -> None:
    run = _run(tmp_path, [collections_event(zero("cog", "kegg"))])
    page, _, _ = report(run, tmp_path)
    found = _first_sentences(page)
    assert len(found) == 1 and found[0].startswith(FIRST_COG_KEGG), found


def test_notice_names_all_three_when_all_keys_are_present(tmp_path: Path) -> None:
    run = _run(tmp_path, [collections_event(zero())])
    page, _, _ = report(run, tmp_path)
    found = _first_sentences(page)
    assert len(found) == 1 and found[0].startswith(FIRST), found
