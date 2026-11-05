"""M12-04 PR 8 (F-01): one notice when no function term was assigned
(docs/specs/M12-04-outputs-revision.md § PR 8, AC 8.1-8.4).

The fixtures are fixture P (A-like, ``core, bgc, function``) with the ``function`` step's
events and tables written by hand (``helpers.m12_04_pr8``); V08 and V14 are rendered for real
and the page is built as the ``plots`` step builds it.

Interface pinned here (PR 8 test author, beyond the spec text):
- coverage per collection is read from the last ``function.collections`` event, then from
  ``category_summary.coverage``, then "not recorded" (behaviour 1). A collection counts as
  "no term" when its ``coverage`` is 0 **or** its ``n_terms`` is 0 (a ``null`` coverage with
  0 terms counts);
- "every collection at coverage 0" needs a ``function.collections`` event with at least one
  collection; with no event (an old run) nothing changes;
- the notice is the section's only ``p.notice`` and starts with the spec's two sentences (file
  links after it are allowed). Its first sentence names COG, GO and KEGG, or only COG and GO
  when the event has no ``kegg`` key; it holds no raw flag or skip token;
- the coverage note names collections in upper case with a whole percent (``fmt_pct``): "COG
  72 %", "GO 41 %" (every GO collection is at 41 % in these fixtures, so how GO's four
  collections are combined is not pinned), "KEGG 66 %", or "KEGG not installed" when the event
  has no ``kegg`` key;
- the KEGG licence note of the section ("academic use only", "KEGG retrieval date") prints only
  when the event has ``kegg`` with ``n_terms > 0``; with no event it prints as before (unless
  ``layer_missing:kegg``), so the M9-16 and M12-03 report tests keep passing.
"""

from __future__ import annotations

import functools
import json
import re
import tempfile
from pathlib import Path

from helpers.m12_04_pr2 import flag_tokens, internal_refs, prose_nodes, sim_report, words
from helpers.m12_04_pr8 import (
    ACADEMIC,
    CHECK,
    FIRST,
    FIRST_NO_KEGG,
    GO,
    LIGHT,
    NO_TERMS,
    NOT_RECORDED,
    NOTICE_CHECK,
    NOTICE_LIGHT,
    RETRIEVAL,
    collections_event,
    function_run,
    function_section,
    notices,
    partial,
    report,
    section_text,
    write_function,
    zero,
)
from helpers.report_fixtures import Page, parse
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st
from unit.report.test_m12_04_pr1_fixtures import MODULES_BC, Run, fixture_p


def _f01(page: Page) -> list[str]:
    """The F-01 notices of the Function section (by their first sentence)."""
    return [n for n in notices(function_section(page))
            if n.startswith(FIRST) or n.startswith(FIRST_NO_KEGG)]  # fmt: skip


def _coverage_unknown(text: str) -> bool:
    """Whether the coverage note says "not recorded" (the KEGG retrieval date's own "not
    recorded in this run" is a different note and does not count)."""
    return re.search(r"[Cc]overage[^.]{0,60}" + NOT_RECORDED, text) is not None


# ---------------------------------------------------------------- AC 8.1


def test_ac_8_1_every_coverage_zero_is_one_notice_and_no_figure(tmp_path: Path) -> None:
    """Run A: every collection at coverage 0 with 0 terms, Bakta light."""
    run = function_run(tmp_path, collections=zero())
    page, plots, events = report(run, tmp_path)
    sec = function_section(page)
    found = notices(sec)
    assert len(found) == 1, f"the Function section holds {len(found)} notices: {found}"
    assert found[0].startswith(NOTICE_LIGHT), found[0]
    assert not sec.by_tag("figure"), "a figure is still drawn in the Function section"
    assert not sec.by_tag("svg") and not sec.by_tag("img"), "an image in the Function section"
    assert not sec.by_class("twin"), "an interactive twin in the Function section"
    assert not list(plots.glob("V08_cog_abundance.*")), sorted(p.name for p in plots.iterdir())
    skipped = [e.data for e in events if e.event == "plots.skipped"]
    assert {"plot": "V08", "reason": NO_TERMS} in skipped, events


def test_ac_8_1_the_section_is_replaced_not_annotated(tmp_path: Path) -> None:
    """Behaviour: the section is replaced by the notice, so no caveat, coverage note or KEGG
    licence note is printed beside it (examples rows 2, 3)."""
    run = function_run(tmp_path, collections=zero())
    page, _, _ = report(run, tmp_path)
    text = section_text(page)
    for gone in (ACADEMIC, RETRIEVAL, "coverage per collection", NOT_RECORDED,
                 "GSEA results are conditional", "No GSEA or detection ORA rows"):  # fmt: skip
        assert gone not in text, (gone, text)
    assert not function_section(page).by_class("caveat"), text


def test_ac_8_1_notice_prose_has_no_internal_token(tmp_path: Path) -> None:
    run = function_run(tmp_path, collections=zero())
    page, _, _ = report(run, tmp_path)
    for node in prose_nodes(function_section(page)):
        text = words(node.text())
        assert not flag_tokens(text), (flag_tokens(text), text)
        assert not internal_refs(text), (internal_refs(text), text)
        assert "bakta_db_light" not in text and "no_terms" not in text, text


def test_second_sentence_without_bakta_light_says_check_the_source(tmp_path: Path) -> None:
    run = function_run(tmp_path, collections=zero(), light=False)
    page, _, _ = report(run, tmp_path)
    found = _f01(page)
    assert len(found) == 1 and found[0].startswith(NOTICE_CHECK), notices(function_section(page))
    assert LIGHT not in section_text(page)


def test_light_sentence_only_with_the_bakta_light_flag(tmp_path: Path) -> None:
    run = function_run(tmp_path, collections=zero(), light=True)
    page, _, _ = report(run, tmp_path)
    assert CHECK not in section_text(page)


def test_no_kegg_installed_notice_names_only_cog_and_go(tmp_path: Path) -> None:
    """Edge case #543: no ``kegg`` key; COG and GO at coverage 0: the notice applies and its
    first sentence names only COG and GO."""
    run = function_run(tmp_path, collections=zero("cog", *GO))
    page, plots, _ = report(run, tmp_path)
    found = notices(function_section(page))
    assert len(found) == 1, found
    assert found[0].startswith(f"{FIRST_NO_KEGG} {LIGHT}"), found[0]
    assert not list(plots.glob("V08_cog_abundance.*"))


def test_zero_terms_with_null_coverage_counts_as_no_term(tmp_path: Path) -> None:
    """``coverage 0 (or n_terms = 0)``: a null coverage (no eligible gene) with 0 terms."""
    colls: dict[str, tuple[float | None, int]] = {c: (None, 0) for c in ("cog", "kegg", *GO)}
    run = function_run(tmp_path, collections=colls)
    page, _, _ = report(run, tmp_path)
    assert len(_f01(page)) == 1, notices(function_section(page))


# ---------------------------------------------------------------- AC 8.2


def test_ac_8_2_coverage_from_the_event_when_category_summary_is_header_only(
    tmp_path: Path,
) -> None:
    run = function_run(tmp_path, collections=partial(), n_cog=60)
    page, _, _ = report(run, tmp_path)
    text = section_text(page)
    for want in ("COG 72 %", "GO 41 %", "KEGG 66 %"):
        assert want in text, (want, text)
    assert not _coverage_unknown(text), text
    assert not _f01(page), "partial coverage must keep the section"
    assert page.slot("V08").by_tag("figure"), "V08 is not drawn at partial coverage"


def test_ac_8_2_event_comes_before_category_summary(tmp_path: Path) -> None:
    run = function_run(tmp_path, collections=partial(), category={"cog": 0.3, "kegg": 0.2},
                       n_cog=60)  # fmt: skip
    text = section_text(report(run, tmp_path)[0])
    assert "COG 72 %" in text, text
    for stale in ("0.3", "30 %", "0.2", "20 %"):
        assert stale not in text, (stale, text)


def test_ac_8_2_category_summary_when_the_event_is_missing(tmp_path: Path) -> None:
    """Old run, no event: the category summary's coverage is shown, never "not recorded"."""
    run = function_run(tmp_path, collections=None, category={"cog": 0.3}, n_cog=60)
    page, _, _ = report(run, tmp_path)
    text = section_text(page)
    assert not _coverage_unknown(text), text
    assert "0.3" in text or "30 %" in text, text
    assert not _f01(page)


def test_ac_8_2_not_recorded_only_when_both_sources_are_absent(tmp_path: Path) -> None:
    """Edge case: no event and a header-only category summary: coverage is unknown, so the
    section says "not recorded" and shows no notice."""
    run = function_run(tmp_path, collections=None)
    page, _, _ = report(run, tmp_path)
    text = section_text(page)
    assert _coverage_unknown(text), text
    assert not _f01(page), notices(function_section(page))
    assert page.slot("V08").by_tag("figure"), "an old run's V08 must still be drawn"


def test_ac_8_2_sim_empty_reads_coverage_from_its_event() -> None:
    """The ``sim_empty`` snapshot run: header-only ``category_summary.tsv`` and a
    ``function.collections`` event with real coverages (cog 0.988, go_slim 0.658, kegg 0.658)."""
    page, _ = sim_report("sim_empty")
    text = section_text(page)
    assert not _coverage_unknown(text), text
    assert "COG 99 %" in text and "KEGG 66 %" in text, text


# ---------------------------------------------------------------- AC 8.3


def test_ac_8_3_kegg_note_absent_when_kegg_has_no_terms(tmp_path: Path) -> None:
    colls = partial() | {"kegg": (0.0, 0)}
    run = function_run(tmp_path, collections=colls, n_cog=60)
    text = section_text(report(run, tmp_path)[0])
    assert ACADEMIC not in text and RETRIEVAL not in text, text


def test_ac_8_3_kegg_note_absent_when_kegg_is_not_in_the_event(tmp_path: Path) -> None:
    """KEGG absent from the event, even without the ``layer_missing:kegg`` flag."""
    run = function_run(tmp_path, collections=partial(kegg=None), kegg_flag=False, n_cog=60)
    text = section_text(report(run, tmp_path)[0])
    assert ACADEMIC not in text and RETRIEVAL not in text, text
    assert "KEGG not installed" in text, text


def test_ac_8_3_kegg_note_present_when_kegg_has_terms(tmp_path: Path) -> None:
    run = function_run(tmp_path, collections=partial(), n_cog=60)
    text = section_text(report(run, tmp_path)[0])
    assert ACADEMIC in text and "KEGG retrieval date: 2026-10-15" in text, text


# ---------------------------------------------------------------- AC 8.4


def test_ac_8_4_no_kegg_with_cog_terms_draws_v08_and_says_kegg_not_installed(
    tmp_path: Path,
) -> None:
    run = function_run(tmp_path, collections=partial(kegg=None), n_cog=60)
    page, plots, _ = report(run, tmp_path)
    text = section_text(page)
    assert "KEGG not installed" in text, text
    assert re.search(r"COG 72 %.{0,40}KEGG not installed|KEGG not installed.{0,40}COG 72 %",
                     text), text  # fmt: skip
    assert page.slot("V08").by_tag("figure"), "V08 must be drawn"
    assert list(plots.glob("V08_cog_abundance.svg"))
    assert not _f01(page)


# ---------------------------------------------------------------- unchanged cases


def test_function_off_is_left_to_not_applicable(tmp_path: Path) -> None:
    """Runs B, C: ``function`` off; this PR does nothing (V08 is a core plot and stays)."""
    run = fixture_p(tmp_path, "B", modules=MODULES_BC)
    page, plots, events = report(run, tmp_path)
    assert not any(e.data.get("reason") == NO_TERMS for e in events), events
    assert list(plots.glob("V08_cog_abundance.svg"))
    assert NO_TERMS not in page.html and FIRST not in page.html


# ---------------------------------------------------------------- property


def _write_events(out: Path, colls: dict[str, tuple[float, int]], light: bool) -> None:
    write_function(out, collections=colls, category="header", light=light)
    if not light:
        (out / "logs" / "steps" / "annotate.events.jsonl").unlink(missing_ok=True)


@functools.cache
def _base_run() -> Run:
    """One fixture P for every Hypothesis example (only the events are rewritten)."""
    tmp = Path(tempfile.mkdtemp(prefix="gippa-m1204-pr8-"))
    return function_run(tmp, collections=zero(), n_cog=60)


_COV = st.one_of(st.just((0.0, 0)),
                 st.tuples(st.floats(0.01, 1.0), st.integers(1, 500)))  # fmt: skip


@settings(max_examples=30, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(cog=_COV, go=st.lists(_COV, min_size=4, max_size=4), kegg=st.none() | _COV,
       light=st.booleans())  # fmt: skip
def test_property_notice_iff_no_collection_has_a_term(
    cog: tuple[float, int], go: list[tuple[float, int]], kegg: tuple[float, int] | None,
    light: bool,
) -> None:  # fmt: skip
    """For any ``function.collections`` event: exactly one F-01 notice iff every collection has
    coverage 0; its first sentence names KEGG iff the event has ``kegg``; its second sentence
    follows ``bakta_db_light``; the KEGG licence note prints iff ``kegg`` has terms (and the
    section is kept); "KEGG not installed" iff the event has no ``kegg``."""
    from gippa.report.build import build_report

    colls = {"cog": cog} | dict(zip(GO, go, strict=True))
    if kegg is not None:
        colls["kegg"] = kegg
    run = _base_run()
    _write_events(run.out, colls, light)
    assert json.loads((run.out / "logs/steps/function.events.jsonl").read_text().splitlines()[1])[
        "data"] == collections_event(colls)  # fmt: skip
    page = parse(build_report(run.view(), {}, {}, run.cfg()))
    none = all(n == 0 for _, n in colls.values())
    found = _f01(page)
    text = section_text(page)
    if none:
        assert len(found) == 1, (colls, notices(function_section(page)))
        first = FIRST if kegg is not None else FIRST_NO_KEGG
        assert found[0].startswith(f"{first} {LIGHT if light else CHECK}"), found[0]
        assert ACADEMIC not in text, text
    else:
        assert not found, (colls, found)
        kegg_terms = kegg is not None and kegg[1] > 0
        assert (ACADEMIC in text) == kegg_terms, (colls, text)
        assert ("KEGG not installed" in text) == (kegg is None), (colls, text)
