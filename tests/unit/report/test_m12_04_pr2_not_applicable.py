"""M12-04 PR 2 behaviour 4 (last bullet), AC 2.7: inapplicable sections move to one
``<section id="sec-not-applicable">`` at the end (docs/specs/M12-04-outputs-revision.md § PR 2).

Interface assumed (PR 2 test author):
- a section is inapplicable when its module is off or an input it needs is absent: here
  Variants without ``variants.tsv``, Novel ORFs without ``novel_orfs.tsv``, Operons without
  ``pairs.tsv`` and ``operon_flags.tsv``, Function with the ``function`` module off and no
  drawn instance in any of its plot slots (review r1 M1);
- ``section#sec-not-applicable`` is the last ``<section>``, headed "Not applicable in this run"
  (the spec's example wording); each entry is an ``<li>`` "<Section heading>: <plain reason>";
- an inapplicable section is moved, not repeated: no ``<section id="<its id>">`` of its own;
- the property runs with ``mobile`` and ``bgc`` on, so module-level entries without a section
  (Mobile, BGC content; edge cases 1 and 3) cannot appear in it.
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import pytest
from helpers.m12_04_pr2 import (
    dress,
    flag_tokens,
    internal_refs,
    not_applicable,
    not_applicable_entries,
    section,
    words,
)
from helpers.report_content import (
    write_novel,
    write_operon_flags,
    write_pairs,
    write_sensitivity,
    write_variants,
)
from helpers.report_fixtures import Page
from hypothesis import HealthCheck, example, given, settings
from hypothesis import strategies as st
from unit.report.test_m12_04_pr1_fixtures import Run, build, fixture_p

HEADING = "Not applicable in this run"
#: optional section -> (section id, heading)
OPTIONAL = {"variants": ("variants", "Variants"), "novel": ("novel-orfs", "Novel ORFs"),
            "operons": ("operons", "Operons"), "function": ("function", "Function")}  # fmt: skip


def _check_entries(page: Page) -> list[tuple[str, str]]:
    entries = not_applicable_entries(page)
    for name, reason in entries:
        assert name, f"entry without a section name: {reason!r}"
        assert reason.strip(), f"{name}: empty reason"
        assert not flag_tokens(reason), f"{name}: flag token in {reason!r}"
        assert not internal_refs(reason), f"{name}: internal reference in {reason!r}"
    return entries


def _check_section_shape(page: Page) -> None:
    sec = not_applicable(page)
    assert sec is not None
    assert page.sections()[-1] is sec, "#sec-not-applicable is not the last section"
    heads = [words(h.text()) for h in sec.by_tag("h2")]
    assert heads == [HEADING], heads


# ---------------------------------------------------------------- examples


def test_proteins_only_run_lists_variants(tmp_path: Path) -> None:
    """R (B-like): no variant table, so "Variants: no variant table in this run"."""
    page = dress(build("R", tmp_path), tmp_path).page()
    _check_section_shape(page)
    entries = dict(_check_entries(page))
    assert "Variants" in entries, entries
    assert entries["Variants"].startswith("no variant table in this run"), entries
    assert section(page, "variants") is None, "Variants is listed and still rendered"
    assert not page.root.by_attr("data-plot", "V12"), "V12 still has an inline slot"


def test_core_only_run_lists_bgc_and_keeps_function(tmp_path: Path) -> None:
    """``--modules core`` (edge case 3): the BGC content is not applicable; Function stays,
    because its V08 (COG from the annotation) is drawn (review r1 M1)."""
    page = dress(build("P_core", tmp_path), tmp_path).page()
    entries = _check_entries(page)
    names = [n for n, _ in entries]
    assert "Function" not in names, names
    assert section(page, "function") is not None
    assert any("BGC" in n or "iosynthetic" in n for n, _ in entries), entries


def test_core_only_run_without_a_drawn_function_plot_lists_function(tmp_path: Path) -> None:
    """``--modules core`` with V08 an empty state: nothing of Function is drawn, so it is
    listed (review r1 M1)."""
    run = dress(build("P_core", tmp_path), tmp_path, empty={"V08": "No genes"})
    page = run.page()
    names = [n for n, _ in _check_entries(page)]
    assert "Function" in names, names
    assert section(page, "function") is None


@pytest.mark.parametrize("name", ["P", "R", "D", "P_core"])
def test_entries_are_plain(name: str, tmp_path: Path) -> None:
    page = dress(build(name, tmp_path), tmp_path).page()
    assert _check_entries(page), f"{name}: no not-applicable entry"


# ---------------------------------------------------------------- property


def _run(tmp: Path, *, variants: bool, novel: bool, operons: bool, function: bool,
         v08_drawn: bool = True) -> Run:  # fmt: skip
    modules = ("core", "mobile", "bgc", "function") if function else ("core", "mobile", "bgc")
    run = fixture_p(tmp, "NA", modules=modules)
    out = run.out
    (out / "tables" / "variants.tsv").unlink()
    if variants:
        write_variants(out, [("v1", "high", "0.01")])
    if novel:
        write_novel(out, [("N1", True, "0.01")])
    if operons:
        write_pairs(out)
        write_operon_flags(out, [0.01])
    write_sensitivity(out, [("islands", "block_x2", 2, 3, 2, 0.6667)])
    return dress(run, tmp, None if v08_drawn else {"V08": "No genes"})


_FLAGS = st.fixed_dictionaries({k: st.booleans() for k in OPTIONAL})


@settings(max_examples=16, deadline=None,
          suppress_health_check=[HealthCheck.too_slow])  # fmt: skip
@given(_FLAGS, st.booleans())
@example(dict.fromkeys(OPTIONAL, True), True)
@example(dict.fromkeys(OPTIONAL, False), True)
@example(dict.fromkeys(OPTIONAL, False), False)
def test_not_applicable_iff_some_section_is_inapplicable(
    present: dict[str, bool], v08_drawn: bool
) -> None:
    """Review r1 M1: Function is inapplicable only with its module off and no drawn plot (V08
    drawn from the annotation's COG categories keeps it)."""
    with tempfile.TemporaryDirectory() as tmp:
        page = _run(Path(tmp), **present, v08_drawn=v08_drawn).page()
    applicable = dict(present, function=present["function"] or v08_drawn)
    want = {OPTIONAL[k][1] for k, on in applicable.items() if not on}
    entries = _check_entries(page)
    assert {n for n, _ in entries} == want, (present, entries)
    assert len(entries) == len(want), entries
    if want:
        _check_section_shape(page)
    else:
        assert not_applicable(page) is None
        assert not page.root.by_attr("id", "sec-not-applicable")
    for key, (sid, _) in OPTIONAL.items():
        rendered = section(page, sid) is not None
        on = applicable[key]
        assert rendered == on, f"{sid}: rendered={rendered}, applicable={on}"
