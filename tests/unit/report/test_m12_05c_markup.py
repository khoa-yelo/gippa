"""M12-05c AC 2-10 and examples 1-10: the guide markup in ``report.html`` and ``compare.html``
(docs/specs/M12-05c-plot-guides.md § 2-7; lead decisions 2026-11-09: the five labels, always
shown, no extra heading, no ``<details>``).

Fixtures: the M12-04 PR 1 unit fixtures dressed with one drawn instance per plot
(``helpers.m12_04_pr2.dress``), a 3-locus V03 slot, a V13 slot with its index, the M12-04 PR 13
linked-view ``OUT/`` and the M9-25 ``sim_plots`` / ``sim_empty`` snapshot fixtures (report and
compare). Interface: tests/helpers/m12_05c.py.
"""

from __future__ import annotations

import dataclasses
import functools
import re
import tempfile
from pathlib import Path
from typing import Any

import pytest
from helpers.m12_04_pr2 import QUESTIONS, dress, flag_tokens, internal_refs, sim_report, words
from helpers.m12_05c import (
    INDEX_STEMS,
    LETTERS,
    SHORT_LABELS,
    data_links,
    data_text,
    field,
    fields,
    figcaption_of,
    full_labels,
    guide_dls,
    guide_html,
    guide_of,
    guides,
    labels,
    links,
    panel_items,
    plot_of,
    slot_figures,
)
from helpers.report_fixtures import NAMES, Page, figure, frame, parse
from unit.report.test_m12_04_pr1_fixtures import build

V03_SUFFIXES = ["bgc_1", "bgc_2", "island_1"]


# ---------------------------------------------------------------- builders


def _write(pid: str, suffixes: list[str | None], tmp: Path, formats: list[str] | None = None,
           empty: frozenset[str | None] = frozenset()) -> Any:  # fmt: skip
    """``PlotFiles`` of ``pid`` with one drawn (or empty-state) instance per suffix."""
    from gippa.io.plot_files import read_plot_files, write_instance
    from gippa.viz import base

    for suffix in suffixes:
        if suffix in empty:
            inst = dataclasses.replace(
                base.empty_state(pid, NAMES[pid], "Nothing here", ["x"], "single"), suffix=suffix
            )
        else:
            inst = base.Instance(suffix=suffix, figure=figure(), data=frame(),
                                 title=f"{pid} {suffix}: a drawn finding")  # fmt: skip
        write_instance(inst, pid, NAMES[pid], tmp, formats or ["svg"], 50)
    return read_plot_files(tmp)[pid]


def _dressed(name: str, tmp: Path, **kw: Any) -> Any:
    return dress(build(name, tmp), tmp, **kw)


def _v13_run(tmp: Path) -> Any:
    from helpers.report_content import write_novel

    run = build("P", tmp)
    write_novel(run.out, [("N1", True, "0.01"), ("N2", False, "0.4")])
    dress(run, tmp)
    run.plots["V13"] = _write("V13", [None, "locus", "locus_N1", "locus_N2"], tmp / "v13")
    return run


@functools.cache
def _sim_compare(fixture: str) -> Page:
    from helpers.snapshot_fixture import build_compare, missing

    why = missing(fixture)
    if why is not None:
        raise RuntimeError(why)
    tmp = Path(tempfile.mkdtemp(prefix=f"gippa-m1205c-cmp-{fixture}-"))
    return parse(build_compare(fixture, tmp).read_text(encoding="utf-8"))


def _sim_plots(fixture: str) -> tuple[Page, Path, dict[str, Any]]:
    from gippa.io.plot_files import read_plot_files

    page, tmp = sim_report(fixture)
    return page, tmp, dict(read_plot_files(tmp / "plots"))


# ---------------------------------------------------------------- the page checker (AC 2-5, 8)


def check_page(page: Page, plots: dict[str, Any]) -> None:
    """Every figure of the report carries one guide: full on the first figure of its plot, short
    on the later ones; fields, Panels, ids, links and Data as in AC 2-5 and AC 8."""
    figures = [f for f in page.root.by_tag("figure") if plot_of(f) is not None]
    assert figures, "no <figure> in a plot slot"
    seen: set[str] = set()
    for fig in figures:
        pid = plot_of(fig)
        assert pid is not None
        g = guides()[pid]
        cap = figcaption_of(fig)
        stem = cap.attrs.get("data-stem", "")
        dl = guide_of(fig)
        assert dl.attrs.get("data-plot") == pid, (stem, dl.attrs)
        got = labels(dl)
        for lab, dd in fields(dl):  # AC 2: every dd has text; AC 8: wording
            text = words(dd.text())
            assert text, f"{stem}: empty dd for {lab!r}"
            for t in (lab, text):
                assert not internal_refs(t) and not flag_tokens(t), (stem, t)
        q = words(field(dl, "Question").text())
        assert q == g.question == QUESTIONS[pid], (stem, q)
        assert words(cap.text()).startswith(f"Question {QUESTIONS[pid]} "), words(cap.text())
        if pid not in seen:  # the first drawn figure: the full guide (AC 2, 4)
            seen.add(pid)
            assert got == full_labels(g), (stem, got)
            assert dl.attrs.get("id") == f"guide-{pid}", (stem, dl.attrs)
            assert not dl.has_class("guide-short"), stem
            items = panel_items(dl)  # AC 3
            assert items == [words(p) for p in g.panels], (stem, items)
            if len(items) > 1:
                assert [i[:4] for i in items] == [f"({c}) " for c in LETTERS[: len(items)]]
            assert words(field(dl, "How to read").text()) == words(g.read), stem
            if g.caveats is not None:
                assert words(field(dl, "Caveats").text()) == words(g.caveats), stem
        else:  # a later figure: the short guide, linking the full one (AC 2, 4)
            assert got == SHORT_LABELS, (stem, got)
            assert "id" not in dl.attrs, (stem, dl.attrs)
            assert dl.has_class("guide-short"), (stem, dl.attrs)
            read = field(dl, "How to read")
            if stem != INDEX_STEMS.get(pid):
                assert words(read.text()).startswith("As for the first"), words(read.text())
                assert "see the guide above" in words(read.text()), words(read.text())
                assert ("#guide-" + pid) in [h for _, h in links(read)], links(read)
        inst = next(i for i in plots[pid].instances if i.stem == stem)  # AC 5
        data = field(dl, "Data")
        assert words(data.text()) == data_text(stem, inst), (stem, words(data.text()))
        assert links(data) == data_links(inst), (stem, links(data))
    for pid in seen:  # AC 4: one full guide per plot per page
        assert len(page.root.by_attr("id", f"guide-{pid}")) == 1, pid


@pytest.mark.parametrize("name", ["P", "R", "D"])
def test_every_figure_has_its_guide_on_the_unit_fixtures(name: str, tmp_path: Path) -> None:
    run = _dressed(name, tmp_path)
    check_page(run.page(), run.plots)


def test_every_figure_has_its_guide_on_sim_plots() -> None:
    page, _, plots = _sim_plots("sim_plots")
    check_page(page, plots)


def test_guide_links_resolve_on_sim_plots() -> None:
    """AC 5: every guide link names a file next to the report, or an id on the page."""
    page, tmp, _ = _sim_plots("sim_plots")
    ids = set(page.ids())
    dls = [dl for cap in page.root.by_tag("figcaption") for dl in guide_dls(cap)]
    assert dls, "no guide on the page"
    for dl in dls:
        for _, href in links(dl):
            if href.startswith("#"):
                assert href[1:] in ids, href
            else:
                assert (tmp / "report" / href).resolve().is_file(), href


# ---------------------------------------------------------------- examples


def test_example_1_v02_full_guide(tmp_path: Path) -> None:
    run = _dressed("P", tmp_path)
    run.plots["V02"] = _write("V02", [None], tmp_path / "v02", ["svg", "png", "pdf"])
    (fig,) = slot_figures(run.page(), "V02")
    dl = guide_of(fig)
    assert dl.attrs.get("id") == "guide-V02" and dl.attrs.get("data-plot") == "V02"
    assert labels(dl) == ("Question", "Panels", "How to read", "Caveats", "Data")
    assert words(field(dl, "Question").text()) == "Is abundance clustered along the chromosome?"
    assert len(panel_items(dl)) == 1
    found = links(field(dl, "Data"))
    assert len(found) == 4, found
    assert found[0] == ("plots/data/V02_abundance_landscape.tsv",
                        "../plots/data/V02_abundance_landscape.tsv")  # fmt: skip
    assert sorted(lab for lab, _ in found[1:]) == ["PDF", "PNG", "SVG"], found
    assert figcaption_of(fig).attrs.get("data-stem") == "V02_abundance_landscape"


def test_example_2_v17_no_caveats_five_lettered_panels(tmp_path: Path) -> None:
    (fig,) = slot_figures(_dressed("P", tmp_path).page(), "V17")
    dl = guide_of(fig)
    assert labels(dl) == ("Question", "Panels", "How to read", "Data")
    assert [i[:4] for i in panel_items(dl)] == ["(a) ", "(b) ", "(c) ", "(d) ", "(e) "]


def test_example_3_v03_three_locus_views(tmp_path: Path) -> None:
    run = _dressed("P", tmp_path)
    run.plots["V03"] = _write("V03", list(V03_SUFFIXES), tmp_path / "v03")
    page = run.page()
    figs = slot_figures(page, "V03")
    assert [figcaption_of(f).attrs.get("data-stem") for f in figs] == [
        f"V03_locus_{s}" for s in V03_SUFFIXES]  # fmt: skip
    first = guide_of(figs[0])
    assert first.attrs.get("id") == "guide-V03"
    assert labels(first) == full_labels(guides()["V03"])
    for fig in figs[1:]:
        dl = guide_of(fig)
        stem = figcaption_of(fig).attrs["data-stem"]
        assert dl.has_class("guide") and dl.has_class("guide-short")
        assert labels(dl) == SHORT_LABELS and "id" not in dl.attrs
        assert ("#guide-V03") in [h for _, h in links(dl)]
        assert (f"plots/data/{stem}.tsv", f"../plots/data/{stem}.tsv") in links(dl)
    assert len(page.root.by_attr("id", "guide-V03")) == 1
    check_page(page, run.plots)


def test_first_drawn_instance_carries_the_full_guide(tmp_path: Path) -> None:
    """§ 6: the first V03 instance is an empty state (a notice), so the second carries the full
    guide."""
    run = _dressed("P", tmp_path)
    run.plots["V03"] = _write("V03", list(V03_SUFFIXES), tmp_path / "v03",
                              empty=frozenset({"bgc_1"}))  # fmt: skip
    page = run.page()
    figs = slot_figures(page, "V03")
    stems = [figcaption_of(f).attrs.get("data-stem") for f in figs]
    assert stems == ["V03_locus_bgc_2", "V03_locus_island_1"], stems
    assert guide_of(figs[0]).attrs.get("id") == "guide-V03"
    assert guide_of(figs[1]).has_class("guide-short")
    assert not page.slot("V03").by_class("notice")[0].by_tag("dl")


def test_example_4_v13_index_and_views(tmp_path: Path) -> None:
    run = _v13_run(tmp_path)
    page = run.page()
    figs = {figcaption_of(f).attrs.get("data-stem"): f for f in slot_figures(page, "V13")}
    assert set(figs) == {"V13_novel_orfs", "V13_novel_orfs_locus",
                         "V13_novel_orfs_locus_N1", "V13_novel_orfs_locus_N2"}  # fmt: skip
    main = guide_of(figs["V13_novel_orfs"])
    assert main.attrs.get("id") == "guide-V13"
    assert labels(main) == full_labels(guides()["V13"])
    index = guide_of(figs["V13_novel_orfs_locus"])
    assert labels(index) == SHORT_LABELS
    read = words(field(index, "How to read").text())
    assert re.search(r"\bindex\b", read, re.I), read
    assert not read.startswith("As for the first"), read
    for stem in ("V13_novel_orfs_locus_N1", "V13_novel_orfs_locus_N2"):
        dl = guide_of(figs[stem])
        assert labels(dl) == SHORT_LABELS and dl.has_class("guide-short")
        assert "#guide-V13" in [h for _, h in links(dl)]
    check_page(page, run.plots)


@pytest.mark.parametrize("name", ["P", "D"])
def test_example_5_redact_keeps_the_guides(name: str, tmp_path: Path) -> None:
    """AC 6: ``redact=True`` gives byte-identical ``dl.guide`` markup."""
    from gippa.report.build import build_report

    run = _v13_run(tmp_path) if name == "P" else _dressed(name, tmp_path)
    plain = guide_html(run.html())
    red = guide_html(build_report(run.view(), run.plots, run.skipped, run.cfg(), redact=True))
    assert plain, "no dl.guide in the report"
    assert red == plain


def test_example_6_linked_views_keep_their_guides(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """AC 6 and example 6: in ``linked`` mode each linked V13 locus view keeps its figcaption
    and guide; the guide markup equals the inlined report's."""
    from unit.report.test_m12_04_pr13_report_size import LOCUS, _html, _out

    from gippa.report import build

    out = _out(tmp_path, 200_000)
    monkeypatch.setattr(build, "REPORT_LINK_THRESHOLD", 10**12)
    inlined = _html(out)
    monkeypatch.setattr(build, "REPORT_LINK_THRESHOLD", 300_000)
    linked = _html(out)
    assert linked != inlined, "the fixture did not switch to linked mode"
    assert guide_html(inlined), "no dl.guide in the report"
    assert guide_html(linked) == guide_html(inlined)
    page = parse(linked)
    figs = {figcaption_of(f).attrs.get("data-stem"): f for f in slot_figures(page, "V13")}
    for stem in LOCUS:
        assert figs[stem].by_tag("img"), f"{stem} is not linked"
        dl = guide_of(figs[stem])
        assert labels(dl) == SHORT_LABELS, (stem, labels(dl))
    assert labels(guide_of(figs["V13_novel_orfs"]))[0] == "Question"


def test_example_7_compare_v16_full_guide() -> None:
    """AC 7."""
    page = _sim_compare("sim_plots")
    caps = [c for c in page.root.by_tag("figcaption") if guide_dls(c)]
    assert len(caps) == 1, f"{len(caps)} figcaptions with a guide"
    cap = caps[0]
    assert cap.attrs.get("data-stem") == "V16_compare", cap.attrs
    (dl,) = guide_dls(cap)
    g = guides()["V16"]
    assert dl.attrs.get("id") == "guide-V16" and dl.attrs.get("data-plot") == "V16"
    assert labels(dl) == full_labels(g)
    items = panel_items(dl)
    assert len(items) == 2 and items[0].startswith("(a) ") and items[1].startswith("(b) ")
    assert words(field(dl, "Question").text()) == QUESTIONS["V16"]
    assert words(cap.text()).startswith(f"Question {QUESTIONS['V16']} Panels (a) ")
    assert ("plots/data/V16_compare.tsv", "../plots/data/V16_compare.tsv") in links(dl)


def test_compare_without_v16_has_no_guide() -> None:
    """§ 6 and decision 6: ``sim_empty``'s V16 is an empty state (``data-empty`` SVG), so its
    compare report has no guide; the empty-state figure itself is unchanged."""
    page = _sim_compare("sim_empty")
    caps = page.root.by_tag("figcaption")
    assert [words(c.text()).split()[0] for c in caps] == ["V16_compare"], [
        words(c.text()) for c in caps]  # fmt: skip
    assert not [d for d in page.root.by_tag("dl") if d.has_class("guide")]


def test_example_8_notice_only_plot_has_no_guide(tmp_path: Path) -> None:
    run = _dressed("P", tmp_path, empty={"V12": "No variant sites in this run"})
    page = run.page()
    slot = page.slot("V12")
    assert slot.by_class("notice"), "the V12 notice is gone"
    assert not slot.by_tag("figure")
    assert not [d for d in page.root.by_tag("dl") if d.attrs.get("data-plot") == "V12"]
    check_page(page, run.plots)


def test_example_9_instance_without_data_tsv(tmp_path: Path) -> None:
    from gippa.io.plot_files import InstanceFiles, PlotFiles

    run = _dressed("P", tmp_path)
    stem = "V02_abundance_landscape"
    svg = (
        '<svg xmlns="http://www.w3.org/2000/svg" width="10" height="10">'
        f'<title>{stem}: a finding</title><rect id="r" width="1" height="1"/></svg>'
    )
    images = (f"{stem}.svg", f"{stem}.png", f"{stem}.pdf")
    run.plots["V02"] = PlotFiles("V02", (InstanceFiles(stem, svg, images, None),))
    (fig,) = slot_figures(run.page(), "V02")
    data = field(guide_of(fig), "Data")
    assert words(data.text()) == "no data table for this figure; figure files: SVG PNG PDF"
    assert links(data) == [(p.rsplit(".", 1)[1].upper(), f"../plots/{p}") for p in images]


def test_example_10_v07_figcaption_text(tmp_path: Path) -> None:
    (fig,) = slot_figures(_dressed("P", tmp_path).page(), "V07")
    text = words(figcaption_of(fig).text())
    assert text.startswith(
        "Question Which regions and categories are detected less than detectability predicts? "
        "Panels One row per group"
    ), text


def test_v17_figcaption_reads_top_to_bottom(tmp_path: Path) -> None:
    """§ 3: without CSS the caption reads "Question / text / Panels / list / ..." in order."""
    run = _dressed("P", tmp_path)
    (fig,) = slot_figures(run.page(), "V17")
    g = guides()["V17"]
    inst = run.plots["V17"].instances[0]
    want = " ".join(["Question", g.question, "Panels", *g.panels, "How to read", g.read,
                     "Data", data_text(inst.stem, inst)])  # fmt: skip
    assert words(figcaption_of(fig).text()) == words(want)


# ---------------------------------------------------------------- edge cases


def test_v03_slot_with_only_the_index_table_has_no_guide(tmp_path: Path) -> None:
    from gippa.io.plot_files import InstanceFiles, PlotFiles

    run = _dressed("P", tmp_path)
    header = ["rank", "locus_id", "family", "estimate", "error"]
    rows: list[list[str | None]] = [["1", "BGC_1", "bgc", "q", "0.01"]]
    inst = InstanceFiles("V03_locus", None, (), "data/V03_locus.tsv", table=(header, rows))
    run.plots["V03"] = PlotFiles("V03", (inst,))
    slot = run.page().slot("V03")
    assert slot.by_tag("table"), "the V03 index table is gone"
    assert not slot.by_tag("figure") and not slot.by_tag("dl")


def test_unknown_plot_id_raises_key_error(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """§ 6: a drawn plot without a ``GUIDES`` entry makes ``build_report`` raise ``KeyError``
    naming the ID."""
    from gippa.report import content

    run = _dressed("P", tmp_path)
    monkeypatch.setattr(content, "GUIDES", {k: v for k, v in guides().items() if k != "V02"})
    with pytest.raises(KeyError, match="V02"):
        run.html()


def test_no_extra_heading_and_no_details(tmp_path: Path) -> None:
    """Lead decisions Q1, Q2: no "Instruction" heading, the guide never inside ``<details>``."""
    page = _v13_run(tmp_path).page()
    dls = [d for d in page.root.by_tag("dl") if d.has_class("guide")]
    assert dls
    for dl in dls:
        cur = dl.parent
        while cur is not None:
            assert cur.tag != "details", "a guide inside <details>"
            cur = cur.parent
    assert "Instruction" not in words(page.root.text())


def test_rendering_leaves_the_plot_files_alone() -> None:
    """AC 9: the guide is in the page only, never in a plot file."""
    page, tmp, _ = _sim_plots("sim_plots")
    assert guide_html(page.html), "no dl.guide on the page"
    for path in (tmp / "plots").rglob("*"):
        if path.is_file() and path.suffix in (".svg", ".tsv"):
            assert "How to read" not in path.read_text(encoding="utf-8"), path


# ---------------------------------------------------------------- AC 10 (after M12-05a)


def _m12_05a_landed() -> bool:
    """M12-05a adds ``content.display_order`` (its § 6.3)."""
    from gippa.report import content

    return hasattr(content, "display_order")


@pytest.mark.xfail(not _m12_05a_landed(), strict=True,
                   reason="M12-05a (V13 lists every candidate by q) has not merged")  # fmt: skip
def test_ac_10_v13_guide_follows_m12_05a(tmp_path: Path) -> None:
    g = guides()["V13"]
    assert "every candidate" in g.panels[1], g.panels[1]
    assert "smallest q first" in g.panels[1], g.panels[1]
    page = _v13_run(tmp_path).page()
    figs = {figcaption_of(f).attrs.get("data-stem"): f for f in slot_figures(page, "V13")}
    index = words(field(guide_of(figs["V13_novel_orfs_locus"]), "How to read").text())
    assert "every candidate" in index, index
