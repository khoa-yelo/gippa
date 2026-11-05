# ruff: noqa: RUF003  (the report's wording is Unicode: ≤, ≥, –, ·)
"""M12-04 PR 1, behaviours 1, 2 and 4: the "This run" box, the version line and the table of
contents (AC 1.6, 1.7; docs/specs/M12-04-outputs-revision.md § PR 1, R-01, R-03).

Interface assumed beyond the spec text (PR 1 test author):
- ``#this-run`` is the first element child of ``<body>`` after the ``<h1>``; its whole text
  (whitespace-normalised) holds every line, so the tests do not pin the markup of a line;
- the tier sentence is ``TIER_TEXT[tier]`` with its ``{...}`` placeholders (if any) filled in;
  every literal fragment between placeholders is in the box, and so is the run's alpha
  (``stats.alpha.identification``) for the tiers that have one (mayu, decoy_picked);
- "#sec-<key>": ``<key>`` is the ``SECTIONS`` id (``summary``, ``detection``, …). The target
  ``sec-<key>`` may be the ``<section>`` itself or an element inside it; either way it exists
  once and lies in that section;
- the TOC count is the row count of the section's main table: ``genes.tsv`` for Gene
  detection, ``variants.tsv`` for Variants, ``novel_orfs.tsv`` for Novel ORFs; it is written
  "(<n>)" in the link text.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import pytest
from helpers.report_fixtures import Node, Page, local_refs
from helpers.viz_out_tree import write_schema_table
from unit.report.test_m12_04_pr1_fixtures import (
    FIXTURES,
    GENOME_STEM,
    REPLICON_BP,
    TITLE_V01,
    TITLE_V02,
    TITLE_V07,
    TITLES,
    Run,
    build,
    fake_plots,
    fixture_d,
    fixture_p,
    fixture_r,
    text,
    write_preflight,
    write_run_manifest,
)

import gippa
from gippa.viz import style

INTERNAL = re.compile(r"ADR-\d{4}|§\s?\d|GIPPA_SPEC|docs/methods")
PLACEHOLDER = re.compile(r"\{[^}]*\}")


def _tier_text(tier: str) -> str:
    from gippa.viz.words import TIER_TEXT

    return TIER_TEXT[tier]


def _box(page: Page) -> Node:
    found = page.root.by_attr("id", "this-run")
    assert len(found) == 1, f"#this-run: {len(found)} found"
    return found[0]


def _box_text(run: Run) -> str:
    return text(_box(run.page()))


def _elements(node: Node) -> list[Node]:
    return [c for c in node.children if isinstance(c, Node)]


def _order(page: Page) -> dict[int, int]:
    return {id(n): i for i, n in enumerate(page.root.iter())}


# ---------------------------------------------------------------- AC 1.6: the box exists, first


@pytest.mark.parametrize("name", FIXTURES)
def test_ac_1_6_this_run_box_on_every_fixture(name: str, tmp_path: Path) -> None:
    _box(build(name, tmp_path).page())


def test_this_run_box_on_a_bare_run(tmp_path: Path) -> None:
    """No manifest, no genes: the box still exists (with "not recorded" lines)."""
    from helpers.report_fixtures import cfg, parse, run_out

    from gippa.io.run_view import RunView
    from gippa.report.build import build_report

    out = run_out(tmp_path)
    _box(parse(build_report(RunView.from_outdir(out), {}, {}, cfg())))


def test_this_run_box_is_the_first_element_after_h1(tmp_path: Path) -> None:
    page = fixture_p(tmp_path).page()
    body = page.root.by_tag("body")[0]
    kids = _elements(body)
    h1 = next(i for i, n in enumerate(kids) if n.tag == "h1")
    assert kids[h1 + 1].attrs.get("id") == "this-run", [(n.tag, n.attrs) for n in kids[:4]]


# ---------------------------------------------------------------- the lines of the box


def test_box_example_a_genome_replicons_genes(tmp_path: Path) -> None:
    box = _box_text(fixture_p(tmp_path))
    assert GENOME_STEM in box, box
    assert not re.search(r"nissle_r200k\.(fna|fa|fasta|gz)", box), box
    assert re.search(r"\b1 replicons?\b", box), box
    assert f"{style.fmt_num(REPLICON_BP)} bp" in box, box
    assert re.search(r"\b158 genes\b", box), box


@pytest.mark.parametrize("ext", ["fna", "fa", "fasta", "fna.gz", "fasta.gz", "fa.gz"])
def test_box_genome_stem_drops_every_suffix(ext: str, tmp_path: Path) -> None:
    import json

    run = fixture_p(tmp_path)
    path = run.out / "run_manifest.json"
    doc = json.loads(path.read_text(encoding="utf-8"))
    for i in doc["inputs"]:
        if i["role"] == "genome":
            i["file"] = i["path"] = f"mygenome.{ext}"
    path.write_text(json.dumps(doc, sort_keys=True), encoding="utf-8")
    box = _box_text(run)
    assert "mygenome" in box and f"mygenome.{ext}" not in box and "mygenome." not in box, box


def test_box_replicon_count_and_total_length(tmp_path: Path) -> None:
    run = fixture_p(tmp_path)
    write_preflight(run.out, [], n_replicons=3)
    box = _box_text(run)
    assert re.search(r"\b3 replicons\b", box), box
    assert f"{style.fmt_num(REPLICON_BP // 3 * 3)} bp" in box, box


@pytest.mark.parametrize(
    ("name", "kind"),
    [("P", "peptides (774)"), ("R", "proteins (106)"), ("D", "de novo peptides (21)")],
)
def test_box_input_kind_and_count(name: str, kind: str, tmp_path: Path) -> None:
    box = _box_text(build(name, tmp_path))
    assert kind in box, box


def test_box_peptide_run_does_not_say_proteins(tmp_path: Path) -> None:
    box = _box_text(fixture_p(tmp_path))
    assert not re.search(r"\bproteins \(\d", box), box
    assert "de novo" not in box, box


def test_box_both_inputs_listed(tmp_path: Path) -> None:
    box = _box_text(build("both", tmp_path))
    assert "peptides (774)" in box and re.search(r"\bproteins \(\d+\)", box), box


def test_box_manifest_without_inputs_says_not_recorded(tmp_path: Path) -> None:
    run = fixture_p(tmp_path)
    write_run_manifest(run.out, roles=None, tier="mayu", modules=run.modules)
    write_preflight(run.out, [])
    box = _box_text(run)
    assert "input: not recorded" in box, box


@pytest.mark.parametrize(
    ("name", "abundance"),
    [("P", ("ibaq", "molar")), ("R", ("maxlfq", "relative")), ("D", None)],
)
def test_box_abundance_label(name: str, abundance: tuple[str, str] | None, tmp_path: Path) -> None:
    box = _box_text(build(name, tmp_path))
    if abundance is None:
        assert re.search(r"abundance:?\s+none\b", box, re.I), box
        assert "log2" not in box, box
    else:
        label = style.abundance_label(*abundance)
        assert label and label in box, (label, box)


@pytest.mark.parametrize(
    ("name", "on", "off"),
    [
        ("P", ("core", "bgc", "function"), ("mobile",)),
        ("R", ("core", "bgc"), ("mobile", "function")),
        ("P_nokegg", ("core", "bgc", "function"), ("mobile",)),
        ("P_core", ("core",), ("bgc", "function", "mobile")),
    ],
)
def test_box_lists_the_modules_run(
    name: str, on: tuple[str, ...], off: tuple[str, ...], tmp_path: Path
) -> None:
    box = _box_text(build(name, tmp_path))
    assert re.search(r"modules?:?\s+" + r",\s*".join(on) + r"(?![\w,])", box), box
    for m in off:
        assert not re.search(rf"\b{m}\b", box), (m, box)


@pytest.mark.parametrize("name", FIXTURES)
def test_ac_1_6_box_holds_the_tier_sentence(name: str, tmp_path: Path) -> None:
    run = build(name, tmp_path)
    box = _box_text(run)
    sentence = " ".join(_tier_text(run.tier).split())
    for fragment in PLACEHOLDER.split(sentence):
        assert fragment.strip() in box, (run.tier, fragment, box)
    if run.tier in ("mayu", "decoy_picked"):
        assert "≤ 0.01" in box, box


def test_box_tier_sentence_follows_alpha(tmp_path: Path) -> None:
    run = fixture_p(tmp_path)
    run.overrides = ("stats.alpha.identification=0.05",)
    box = _box_text(run)
    assert "≤ 0.05" in box, box


@pytest.mark.parametrize("name", FIXTURES)
def test_box_holds_no_internal_reference_or_path(name: str, tmp_path: Path) -> None:
    box = _box_text(build(name, tmp_path))
    assert not INTERNAL.search(box), box
    assert "/home/" not in box and "local/" not in box and "farmshare" not in box.lower(), box


# ---------------------------------------------------------------- key findings


def test_findings_are_the_v01_v07_v02_titles_in_order(tmp_path: Path) -> None:
    run = fixture_p(tmp_path)
    run.plots = fake_plots(["V01", "V07", "V02"])
    box = _box_text(run)
    where = [box.find(t) for t in (TITLE_V01, TITLE_V07, TITLE_V02)]
    assert all(i >= 0 for i in where), box
    assert where == sorted(where), where


def test_at_most_three_findings(tmp_path: Path) -> None:
    run = fixture_p(tmp_path)
    run.plots = fake_plots(["V01", "V07", "V02"])
    box = _box(run.page())
    items = [li for li in box.by_tag("li") if any(text(li) == t for t in TITLES.values())]
    assert len(items) <= 3


def test_skipped_v01_gives_v07_and_v02_only(tmp_path: Path) -> None:
    run = fixture_p(tmp_path)
    run.plots = fake_plots(["V07", "V02"])
    run.skipped = {"V01": "input_absent:tables/islands.tsv"}
    box = _box_text(run)
    assert TITLE_V07 in box and TITLE_V02 in box, box
    assert TITLE_V01 not in box, box


def test_no_findings_list_when_no_title_is_available(tmp_path: Path) -> None:
    run = fixture_p(tmp_path)
    run.skipped = {"V01": "module_off:mobile", "V07": "input_absent:x", "V02": "input_absent:y"}
    box = _box(run.page())
    for lst in box.by_tag("ul") + box.by_tag("ol"):
        assert lst.by_tag("li"), "an empty findings list is left in the box"
    assert "finding" not in text(box).lower(), text(box)


# ---------------------------------------------------------------- behaviour 2: the version line


def test_version_line_same_version(tmp_path: Path) -> None:
    box = _box_text(fixture_p(tmp_path))
    assert f"Report built by GIPPA {gippa.__version__}" in box, box
    assert "from a run made with" not in box, box


def test_ac_1_6_version_line_shows_both_versions_when_they_differ(tmp_path: Path) -> None:
    """Spec example K: "Report built by GIPPA 1.1.0 from a run made with GIPPA 0.0.1"."""
    run = fixture_r(tmp_path)
    write_run_manifest(run.out, roles=("genome", "proteins"), tier="protein_list",
                       modules=run.modules, version="0.0.1")  # fmt: skip
    box = _box_text(run)
    assert f"Report built by GIPPA {gippa.__version__} from a run made with GIPPA 0.0.1" in box


# ---------------------------------------------------------------- AC 1.7: the TOC


def _toc(page: Page) -> Node:
    found = [n for n in page.root.by_tag("nav") if n.attrs.get("id") == "toc"]
    assert len(found) == 1, f"nav#toc: {len(found)} found"
    return found[0]


def _key(section: Node) -> str:
    sid = section.attrs["id"]
    return sid[len("sec-") :] if sid.startswith("sec-") else sid


def _toc_links(page: Page) -> list[Node]:
    return [a for a in _toc(page).by_tag("a")]


@pytest.mark.parametrize("name", ["P", "R", "D", "P_core"])
def test_ac_1_7_one_toc_link_per_section_in_order(name: str, tmp_path: Path) -> None:
    page = build(name, tmp_path).page()
    hrefs = [a.attrs.get("href") for a in _toc_links(page)]
    want = [f"#sec-{_key(s)}" for s in page.sections()]
    assert hrefs == want, (hrefs, want)


def test_ac_1_7_toc_keys_are_the_section_ids(tmp_path: Path) -> None:
    from gippa.report.build import SECTIONS

    page = fixture_p(tmp_path).page()
    assert [_key(s) for s in page.sections()] == [s.id for s in SECTIONS]


@pytest.mark.parametrize("name", ["P", "R", "D"])
def test_ac_1_7_every_toc_target_exists_once_in_its_section(name: str, tmp_path: Path) -> None:
    page = build(name, tmp_path).page()
    for section in page.sections():
        target = f"sec-{_key(section)}"
        found = page.root.by_attr("id", target)
        assert len(found) == 1, (target, len(found))
        node: Any = found[0]
        while node is not None and node is not section:
            node = node.parent
        assert node is section, f"#{target} is not in section {section.attrs['id']}"


@pytest.mark.parametrize("name", FIXTURES)
def test_ac_1_7_link_checker_every_local_ref_resolves(name: str, tmp_path: Path) -> None:
    page = build(name, tmp_path).page()
    ids = page.ids()
    assert len(ids) == len(set(ids)), "duplicate ids"
    missing = sorted(set(local_refs(page)) - set(ids))
    assert not missing, missing


def test_toc_comes_after_the_box_and_before_the_first_section(tmp_path: Path) -> None:
    page = fixture_p(tmp_path).page()
    order = _order(page)
    box, toc, first = _box(page), _toc(page), page.sections()[0]
    assert order[id(box)] < order[id(toc)] < order[id(first)]
    assert toc.parent is not None and toc.parent.tag == "body"


def test_toc_counts_the_main_table_rows(tmp_path: Path) -> None:
    """Spec: "Gene detection (158)", "Variants (0)"."""
    page = fixture_p(tmp_path).page()
    by = {a.attrs.get("href"): text(a) for a in _toc_links(page)}
    assert re.search(r"\(158\)$", by["#sec-detection"]), by
    assert re.search(r"\(0\)$", by["#sec-variants"]), by


def test_toc_counts_novel_orfs(tmp_path: Path) -> None:
    run = fixture_d(tmp_path)
    rows = [{"orf_id": f"N{i}", "class": "intergenic", "contig": "chr1", "strand": "+",
             "end": 500 + i, "n_peptides": 1, "class_fdr": "0.01", "called": True}
            for i in range(3)]  # fmt: skip
    write_schema_table(run.out / "tables" / "novel_orfs.tsv", "novel_orfs", rows)
    by = {a.attrs.get("href"): text(a) for a in _toc_links(run.page())}
    assert re.search(r"\(3\)$", by["#sec-novel-orfs"]), by


def test_toc_absent_table_has_no_count(tmp_path: Path) -> None:
    """R has no variants.tsv: "Variants" without "(0)"."""
    by = {a.attrs.get("href"): text(a) for a in _toc_links(fixture_r(tmp_path).page())}
    assert "Variants" in by["#sec-variants"], by
    assert not re.search(r"\(\d+\)", by["#sec-variants"]), by
    for key in ("summary", "methods", "provenance"):
        assert not re.search(r"\(\d+\)", by[f"#sec-{key}"]), by


def test_toc_link_text_is_the_section_heading(tmp_path: Path) -> None:
    page = fixture_p(tmp_path).page()
    links = _toc_links(page)
    for a, section in zip(links, page.sections(), strict=True):
        heading = text(section.by_tag("h2")[0])
        assert text(a).startswith(heading), (text(a), heading)


def test_toc_is_sticky_by_css_only(tmp_path: Path) -> None:
    page = fixture_p(tmp_path).page()
    toc = _toc(page)
    assert not toc.by_tag("script"), "the TOC needs no script"
    assert not any(k.startswith("on") for n in toc.iter() for k in n.attrs), "no inline handlers"
    css = " ".join(n.text() for n in page.root.by_tag("style"))
    rule = re.search(r"(nav#toc|#toc)[^{]*\{[^}]*\}", css)
    assert rule and re.search(r"position:\s*sticky", rule.group(0)), "#toc is not sticky in CSS"
