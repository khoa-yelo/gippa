"""M12-04 PR 2 behaviour 4, AC 2.4-2.6: conditional sections
(docs/specs/M12-04-outputs-revision.md § PR 2).

Interface assumed (PR 2 test author):
- the mapping section keeps ``id="mapping"`` and its ``<h2>`` reads "Mapping";
- the Determinants caveats are ``.caveat`` elements, printed only when ``determinants_model``
  has a finite elastic-net (or any model's) ``cv_r2_mean``;
- a section table preview with 0 body rows is not rendered; an element (``p``, ``li`` or
  ``h3``) of its section instead holds the word "none" and an ``<a href="../<PATHS[key]>">``
  (spec example "Islands: none (full table)"); a 0-row ``sensitivity.tsv`` may show the
  sensitivity sentence instead;
- AC 2.6 covers the previews (``table[data-table]``), not the Provenance tables;
- with ``mobile`` off the Genome context section says "prophage and plasmid regions not
  searched"; with ``layer_missing:kegg`` the Function section stays (COG and GO only, #543).
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest
from helpers.m12_04_pr2 import dress, not_applicable_entries, section, sim_report, words
from helpers.report_content import (
    write_islands,
    write_novel,
    write_operon_flags,
    write_pairs,
    write_sensitivity,
)
from helpers.report_fixtures import Node, Page, table_rows
from helpers.snapshot_fixture import fixture_out
from unit.report.test_m12_04_pr1_fixtures import build

SENSITIVITY_NONE = "No flagged results, so there is nothing to compare across sensitivity reruns."


def _section(page: Page, sid: str) -> Node:
    found = section(page, sid)
    assert found is not None, f"no section#{sid}"
    return found


# ---------------------------------------------------------------- AC 2.4 heading


@pytest.mark.parametrize("name", ["R", "P"])
def test_mapping_heading(name: str, tmp_path: Path) -> None:
    page = build(name, tmp_path).page()
    h2s = [words(h.text()) for h in page.root.by_tag("h2")]
    assert "Peptide mapping" not in h2s, h2s
    heads = [words(h.text()) for h in _section(page, "mapping").by_tag("h2")]
    assert heads == ["Mapping"], heads


# ---------------------------------------------------------------- AC 2.5 caveats


@pytest.mark.parametrize("name", ["P", "D"])
def test_no_determinants_caveat_without_a_fitted_model(name: str, tmp_path: Path) -> None:
    """P: both whole models ``too_few_rows`` (cv_r2_mean NA); D: univariate rows only."""
    sec = _section(build(name, tmp_path).page(), "determinants")
    caveats = [words(c.text()) for c in sec.by_class("caveat")]
    assert not caveats, caveats


def test_determinants_caveats_print_with_a_fitted_model(tmp_path: Path) -> None:
    """The control: R has a fitted elastic-net model."""
    sec = _section(build("R", tmp_path).page(), "determinants")
    assert sec.by_class("caveat")


# ---------------------------------------------------------------- sensitivity sentence


def test_sensitivity_sentence_when_nothing_was_flagged(tmp_path: Path) -> None:
    run = build("P", tmp_path)
    write_sensitivity(run.out, [("islands", "block_x2", 0, 0, 0, None),
                                ("operons", "alpha_0.1", 0, 0, 0, None)])  # fmt: skip
    page = run.page()
    sec = _section(page, "sensitivity")
    assert SENSITIVITY_NONE in words(sec.text())
    assert page.table("sensitivity") is None


def test_sensitivity_tables_when_something_was_flagged(tmp_path: Path) -> None:
    run = build("P", tmp_path)
    write_sensitivity(run.out, [("islands", "block_x2", 2, 3, 2, 0.6667),
                                ("operons", "alpha_0.1", 0, 0, 0, None)])  # fmt: skip
    page = run.page()
    assert page.table("sensitivity") is not None
    assert SENSITIVITY_NONE not in words(page.html)


# ---------------------------------------------------------------- AC 2.6 empty previews


def _empty_table_keys(page: Page) -> list[str]:
    return [t.attrs["data-table"] for t in page.root.by_tag("table")
            if "data-table" in t.attrs and not table_rows(t)[1]]  # fmt: skip


def _none_with_link(sec: Node, href: str, *, alt: str | None = None) -> bool:
    for n in sec.iter():
        if n.tag not in ("p", "li", "h3", "div", "span"):
            continue
        if href not in [a.attrs.get("href") for a in n.by_tag("a")]:
            continue
        text = words(n.text())
        if re.search(r"\bnone\b", text) or (alt is not None and alt in text):
            return True
    return False


def test_empty_islands_preview_is_none_with_link(tmp_path: Path) -> None:
    run = build("P", tmp_path)
    write_islands(run.out, [])
    page = run.page()
    assert page.table("islands") is None, "a 0-row islands table is rendered"
    assert _none_with_link(_section(page, "genome-context"), "../tables/islands.tsv")


def test_no_zero_row_preview_unit(tmp_path: Path) -> None:
    """P plus header-only islands, pairs, operon flags and novel ORFs."""
    run = build("P", tmp_path)
    write_islands(run.out, [])
    write_pairs(run.out)
    write_operon_flags(run.out, [])
    write_novel(run.out, [])
    page = dress(run, tmp_path).page()
    assert not _empty_table_keys(page), _empty_table_keys(page)


def test_no_zero_row_preview_sim_empty() -> None:
    """sim_empty: every section table is header-only; each shows "none" and a link to a file
    that exists."""
    from gippa.io.run_view import PATHS, RunView
    from gippa.report.build import SECTIONS

    page, _ = sim_report("sim_empty")
    assert not _empty_table_keys(page), _empty_table_keys(page)
    out = fixture_out("sim_empty")
    view = RunView.from_outdir(out)
    checked = 0
    for s in SECTIONS:
        sec = section(page, s.id)
        for key in s.tables:
            frame = view.table(key) if view.has(key) else None
            if frame is None or len(frame) or sec is None:
                continue
            href = f"../{PATHS[key]}"
            alt = SENSITIVITY_NONE if key == "sensitivity" else None
            assert _none_with_link(sec, href, alt=alt), f"{key}: no 'none' with {href}"
            assert (out / "report" / href).resolve().is_file(), href
            checked += 1
    assert checked >= 5, checked


# ---------------------------------------------------------------- edge cases


def test_mobile_off_is_a_genome_context_note(tmp_path: Path) -> None:
    """P runs core, bgc, function: mobile is off (runs A, B, C)."""
    page = dress(build("P", tmp_path), tmp_path).page()
    sec = _section(page, "genome-context")
    assert "prophage and plasmid regions not searched" in words(sec.text())


def test_no_kegg_keeps_the_function_section(tmp_path: Path) -> None:
    """P_nokegg has no novel ORF table, so a not-applicable list exists; Function is not in
    it."""
    page = dress(build("P_nokegg", tmp_path), tmp_path).page()
    assert section(page, "function") is not None
    names = [n for n, _ in not_applicable_entries(page)]
    assert names, "no not-applicable list although Novel ORFs has no input"
    assert "Function" not in [n for n, _ in not_applicable_entries(page)]
