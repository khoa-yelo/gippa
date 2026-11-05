# ruff: noqa: RUF001, RUF002, RUF003  (the report's wording is Unicode: −, π₀, ℓ̂, Î₁)
"""M12-04 PR 2 behaviour 2 and 6, AC 2.1 and 2.9: plain wording in the report
(docs/specs/M12-04-outputs-revision.md § PR 2, § 1 rule 2).

Interface assumed (PR 2 test author):
- AC 2.1 runs on every PR 1 fixture (dressed with the real skip reasons and drawn plots,
  ``helpers.m12_04_pr2.dress``) and on the ``sim_plots`` and ``sim_empty`` snapshot reports; the
  flag check covers vocabulary tokens with an underscore or a colon only (lead decision 1);
- an ADR or § reference becomes a caveat plus an ``<a href>`` to
  ``https://khoa-yelo.github.io/gippa/methods/<topic>/`` whose ``<topic>`` is a page under
  ``docs/methods`` (lead decision 3); ``gippa.report.checks.external_loads`` does not count it
  (a link is not a load; Q12);
- the missed-operon decoy table (``data-table="missed_operons_bins"``) heads its counts
  "non-operon pairs", "convergent pairs", "sandwich pairs";
- the island diagnostics table (``data-table="island_diagnostics"``) heads the Moran column
  "lag-1 Moran's I of the detrended TU scores (used to choose the block length)" and prints it
  through ``fmt_num`` (−0.399 -> "−0.40", the spec example);
- the Methods section holds ``<dl id="glossary">`` with one ``<dt>`` per term (text exactly the
  spec's term) followed by a non-empty ``<dd>``.
"""

from __future__ import annotations

import unicodedata
from pathlib import Path

import pytest
from helpers.m12_04_pr2 import (
    METHODS_URL,
    dress,
    methods_links,
    methods_topics,
    prose_violations,
    section,
    sim_report,
    words,
)
from helpers.report_content import (
    write_diagnostics,
    write_islands,
    write_missed_event,
    write_pairs,
)
from helpers.report_fixtures import Node, Page, table_rows
from unit.report.test_m12_04_pr1_fixtures import FIXTURES, build

SIM = ("sim_plots", "sim_empty")
GLOSSARY = ("FDR", "FWER", "q-value", "local FDR (lfdr)", "O − E", "π₀", "block length (ℓ̂)",
            "BB join count")  # fmt: skip
MORAN_LABEL = "lag-1 Moran's I of the detrended TU scores (used to choose the block length)"


def _unit_page(name: str, tmp: Path) -> Page:
    return dress(build(name, tmp), tmp).page()


# ---------------------------------------------------------------- AC 2.1


@pytest.mark.parametrize("name", FIXTURES)
def test_no_internal_reference_or_flag_token_in_prose_unit(name: str, tmp_path: Path) -> None:
    bad = prose_violations(_unit_page(name, tmp_path))
    assert not bad, "\n".join(bad)


@pytest.mark.parametrize("fixture", SIM)
def test_no_internal_reference_or_flag_token_in_prose_sim(fixture: str) -> None:
    page, _ = sim_report(fixture)
    bad = prose_violations(page)
    assert not bad, "\n".join(bad)


def test_rich_caveats_have_no_internal_reference(tmp_path: Path) -> None:
    """P plus islands, diagnostics, pairs and the missed-operon event: the genome-context and
    operon caveats are plain too."""
    run = build("P", tmp_path)
    write_islands(run.out, [0.001, 0.2])
    write_diagnostics(run.out, [_diag(-0.399)])
    write_pairs(run.out)
    write_missed_event(run.out)
    bad = prose_violations(dress(run, tmp_path).page())
    assert not bad, "\n".join(bad)


@pytest.mark.parametrize("fixture", SIM)
def test_methods_sentence_naming_docs_methods_is_gone(fixture: str) -> None:
    page, _ = sim_report(fixture)
    assert "documentation (docs/methods)" not in page.html
    assert "docs/methods" not in words(_section(page, "methods").text())


# ---------------------------------------------------------------- methods links


def _section(page: Page, sid: str) -> Node:
    found = section(page, sid)
    assert found is not None, f"no section#{sid}"
    return found


def _check_link_shape(page: Page) -> list[str]:
    links = methods_links(page.root)
    topics = methods_topics()
    for href in links:
        m = METHODS_URL.match(href)
        assert m is not None, f"methods link of the wrong shape: {href!r}"
        assert m.group(1) in topics, f"{href!r}: no docs/methods/{m.group(1)}.md"
    return links


def test_determinants_caveats_link_the_methods_page(tmp_path: Path) -> None:
    """R has a fitted model, so its determinants caveats print; the ADR-0032/33/34 references
    become a link to the docs methods page."""
    from gippa.report.checks import external_loads

    page = _unit_page("R", tmp_path)
    sec = _section(page, "determinants")
    assert sec.by_class("caveat"), "R: no determinants caveat"
    assert methods_links(sec), "R: no methods link in the determinants section"
    assert _check_link_shape(page)
    assert external_loads(page.html) == [], "a methods link counted as an external load"


def test_methods_links_have_the_docs_site_shape() -> None:
    page, _ = sim_report("sim_plots")
    assert _check_link_shape(page), "sim_plots: no methods link at all"


def test_function_caveats_link_the_methods_page() -> None:
    """sim_plots has GSEA rows (ADR-0026, ADR-0027) and the coverage note (§7.7.1)."""
    page, _ = sim_report("sim_plots")
    assert methods_links(_section(page, "function")), "no methods link in Function"


# ---------------------------------------------------------------- decoy headers, diagnostics


def test_missed_operon_table_headers_are_plain(tmp_path: Path) -> None:
    run = build("P", tmp_path)
    write_pairs(run.out)
    write_missed_event(run.out)
    page = run.page()
    table = page.table("missed_operons_bins")
    assert table is not None, "no missed-operon decoy table"
    header, body = table_rows(table)
    assert body
    for plain in ("non-operon pairs", "convergent pairs", "sandwich pairs"):
        assert plain in header, header
    for raw in ("N_NO", "N_cv", "N_sw"):
        assert raw not in header, header


def _diag(moran: float) -> dict[str, object]:
    return {"contig": "chr1", "scanned": True, "block_rule": "auto", "l_hat": 4,
            "block_len": 9, "morans_i_lag1": moran}  # fmt: skip


def test_island_diagnostics_label_and_number(tmp_path: Path) -> None:
    run = build("P", tmp_path)
    write_islands(run.out, [0.001])
    write_diagnostics(run.out, [_diag(-0.399)])
    table = run.page().table("island_diagnostics")
    assert table is not None, "no island diagnostics table"
    header, body = table_rows(table)
    assert MORAN_LABEL in header, header
    assert not [h for h in header if "Î₁" in h], header
    col = header.index(MORAN_LABEL)
    assert body[0][col] == "−0.40", body[0]


# ---------------------------------------------------------------- AC 2.9 glossary


def _norm(text: str) -> str:
    return unicodedata.normalize("NFC", words(text))


def _glossary(page: Page) -> dict[str, str]:
    found = page.root.by_attr("id", "glossary")
    assert len(found) == 1, f"#glossary: {len(found)} found"
    dl = found[0]
    assert dl.tag == "dl", dl.tag
    assert dl in list(_section(page, "methods").iter()), "#glossary is not in Methods"
    out: dict[str, str] = {}
    term: str | None = None
    for child in dl.children:
        if isinstance(child, str):
            continue
        if child.tag == "dt":
            term = _norm(child.text())
        elif child.tag == "dd" and term is not None:
            out[term] = words(child.text())
    return out


@pytest.mark.parametrize("name", ["P", "D"])
def test_glossary_defines_the_eight_terms(name: str, tmp_path: Path) -> None:
    gloss = _glossary(build(name, tmp_path).page())
    for term in GLOSSARY:
        key = _norm(term)
        assert key in gloss, f"glossary lacks {term!r}; has {sorted(gloss)}"
        assert gloss[key], f"{term!r} has an empty definition"


def test_glossary_on_sim_plots() -> None:
    page, _ = sim_report("sim_plots")
    gloss = _glossary(page)
    assert {_norm(t) for t in GLOSSARY} <= set(gloss), sorted(gloss)
