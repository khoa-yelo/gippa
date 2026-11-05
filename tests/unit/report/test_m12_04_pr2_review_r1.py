# ruff: noqa: RUF001, RUF002, RUF003  (the report's wording is Unicode: ℓ̂, ℓ′, α)
"""M12-04 PR 2 review round 1: the lead decisions M1-M3 and L1-L6 (PR #555).

Interface pinned here (PR 2 review r1 test author, beyond the spec text):
- M1: a section with a module (Function) is listed in ``#sec-not-applicable`` only when no plot
  slot of the section holds a drawn (non-empty) instance. With the ``function`` module off and a
  drawn V08 (COG categories come from the annotation), ``section#function`` renders, V08 keeps
  its figure and the V14 slot carries an inline note "... the function module was not run";
- M2: an empty-state instance is ``<p class="notice">`` holding its title followed by its links
  (``SVG``, ``PNG``, ``PDF``, ``data TSV``, whichever exist, each ``../plots/<file>``); still no
  ``<svg>`` or ``<img>``. Every data TSV under ``plots/data/`` is linked from ``report.html``;
- M3: the glossary term is ``block length (ℓ′)`` (U+2113 U+2032); its ``<dd>`` names the
  estimated correlation length ℓ̂ and ``min(20, 2ℓ̂+1, number of TUs)`` (spaces free);
- L1: each autocorrelation caveat paragraph says the rows it applies to: "(these rows carry the
  flag '<plain(flag)>' in the table)" with the text of ``gippa.viz.words.plain`` (straight or
  curly single quotes), and no raw flag token;
- L2: the index instance of a multi-instance plot (V03: suffix ``None``, stem ``V03_locus``;
  V13: suffix ``locus``, stem ``V13_novel_orfs_locus``, its suffix-``None`` instance is panel
  (a)) has its own caption: it starts with the plot's Question, says "index" and differs from
  the caption of the drawn instances;
- L3: the skip reason ``not_registered`` reads "this plot is not available in this GIPPA
  version" (case-insensitive) in the slot note;
- L4: the sensitivity sentence appears only when every row has non-NA ``n_primary`` and
  ``n_variant``, both 0; NA counts keep the tables; a table without those columns never crashes
  the page;
- L5: the proteins-only Novel ORFs reason says "novel ORF calling needs peptides";
- L6: the exact listed sets per fixture, derived from the rules (Operons: no ``pairs.tsv`` and
  no ``operon_flags.tsv``; Variants: no ``variants.tsv``; Novel ORFs: no ``novel_orfs.tsv``;
  Function: M1; the BGC entry when ``bgc`` is off; Mobile is never listed).
"""

from __future__ import annotations

import re
import tempfile
import unicodedata
from pathlib import Path
from typing import Any

import pytest
from helpers.m12_04_pr2 import (
    QUESTIONS,
    dress,
    flag_tokens,
    not_applicable_entries,
    section,
    sim_report,
    words,
)
from helpers.report_content import (
    manifest,
    page,
    write_manifest,
    write_pairs,
    write_sensitivity,
)
from helpers.report_fixtures import NAMES, Node, Page, run_out
from helpers.viz_out_tree import write_schema_table
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st
from unit.report.test_m12_04_pr1_fixtures import FIXTURES, TITLE_V02, Run, build, fixture_p

SENSITIVITY_NONE = "No flagged results, so there is nothing to compare across sensitivity reruns."
MODULES_NO_FUNCTION = ("core", "mobile", "bgc")
BGC = "BGC"


def _listed(page_: Page) -> set[str]:
    """The listed names, the BGC content entry as ``BGC``."""
    out = set()
    for name, _ in not_applicable_entries(page_):
        out.add(BGC if ("BGC" in name or "iosynthetic" in name) else name)
    return out


def _outer_svgs(node: Node) -> list[Node]:
    out = []
    for n in node.iter():
        if n.tag != "svg":
            continue
        cur, inner = n.parent, False
        while cur is not None and cur is not node:
            inner = inner or cur.tag == "svg"
            cur = cur.parent
        if not inner:
            out.append(n)
    return out


def _render_real(run: Run, pid: str, plots_dir: Path) -> Any:
    """``pid`` rendered for real on ``run`` (SVG only), read back as ``PlotFiles``."""
    from gippa.io.plot_files import read_plot_files, write_instance
    from gippa.steps.plots import load_plot_modules
    from gippa.viz.base import registered

    load_plot_modules()
    spec = registered()[pid]
    found = spec.fn(run.view(), run.cfg())
    assert isinstance(found, list), f"{pid} skipped: {found}"
    for inst in found:
        write_instance(inst, pid, spec.name, plots_dir, ["svg"], 50)
    return read_plot_files(plots_dir)[pid]


# ---------------------------------------------------------------- M1


def test_function_off_with_cog_categories_keeps_the_function_section(tmp_path: Path) -> None:
    """Function off, ``genes.tsv`` with real COG categories: V08 is drawn (it is a core plot),
    so Function renders with V08 and a V14 note, and is not listed."""
    run = dress(fixture_p(tmp_path, "nofn", modules=MODULES_NO_FUNCTION), tmp_path)
    genes = run.view().table("genes_tsv")
    assert genes is not None and set(genes["cog_category"].dropna()) - {"none"}
    run.plots["V08"] = _render_real(run, "V08", tmp_path / "real-v08")
    assert any(not i.empty for i in run.plots["V08"].instances)
    assert run.skipped.get("V14") == "module_off:function", run.skipped
    page_ = run.page()
    sec = section(page_, "function")
    assert sec is not None, f"Function listed: {not_applicable_entries(page_)}"
    assert "Function" not in {n for n, _ in not_applicable_entries(page_)}
    v08 = page_.slot("V08")
    assert v08.by_tag("figure") and _outer_svgs(v08), "V08 not drawn"
    v14 = words(page_.slot("V14").text())
    assert "the function module was not run" in v14, v14
    assert page_.slot("V08") in sec.iter() and page_.slot("V14") in sec.iter()


@pytest.mark.parametrize("v08", ["empty", "absent"])
def test_function_off_without_a_drawn_plot_is_listed(v08: str, tmp_path: Path) -> None:
    """The converse: function off and V08 an empty state (or no V08 file at all): nothing of
    the section is drawn, so Function is listed and not rendered."""
    run = fixture_p(tmp_path, "nofn", modules=MODULES_NO_FUNCTION)
    dress(run, tmp_path, empty={"V08": "No non-partial genes in this run"})
    if v08 == "absent":
        run.plots.pop("V08")
    page_ = run.page()
    entries = dict(not_applicable_entries(page_))
    assert "Function" in entries, entries
    assert "the function module was not run" in entries["Function"], entries
    assert section(page_, "function") is None


def test_function_on_with_v08_empty_is_rendered(tmp_path: Path) -> None:
    """With the module on the section is applicable whatever its plots hold."""
    run = dress(build("P", tmp_path), tmp_path, empty={"V08": "No non-partial genes"})
    page_ = run.page()
    assert section(page_, "function") is not None
    assert "Function" not in _listed(page_)


# ---------------------------------------------------------------- M2


def _notice_links(notice: Node) -> list[tuple[str, str]]:
    return [(words(a.text()), a.attrs.get("href", "")) for a in notice.by_tag("a")]


def test_empty_state_notice_keeps_its_file_links(tmp_path: Path) -> None:
    """V02 empty (SVG, PNG, data TSV written): the notice is the title, then the links."""
    from gippa.io.plot_files import read_plot_files, write_instance
    from gippa.viz import base

    run = dress(build("P", tmp_path), tmp_path)
    plots_dir = tmp_path / "v02"
    inst = base.empty_state("V02", NAMES["V02"], TITLE_V02, ["x"], "single")
    write_instance(inst, "V02", NAMES["V02"], plots_dir, ["svg", "png"], 50)
    run.plots["V02"] = read_plot_files(plots_dir)["V02"]
    files = run.plots["V02"].instances[0]
    assert files.empty and files.data is not None
    want = [("SVG", f"../plots/{files.stem}.svg"), ("PNG", f"../plots/{files.stem}.png"),
            ("data TSV", f"../plots/{files.data}")]  # fmt: skip
    slot = run.page().slot("V02")
    notices = [n for n in slot.by_class("notice") if n.tag == "p"]
    assert len(notices) == 1, [words(n.text()) for n in notices]
    text = words(notices[0].text())
    assert text.startswith(TITLE_V02), text
    links = _notice_links(notices[0])
    assert sorted(links) == sorted(want), links
    assert [lab for lab, _ in links][-1] == "data TSV", links
    assert text.endswith(" ".join(lab for lab, _ in links)), text
    assert not _outer_svgs(slot) and not slot.by_tag("img"), "an image for an empty state"


@pytest.mark.parametrize("fixture", ["sim_plots", "sim_empty"])
def test_every_data_tsv_is_linked_from_the_report(fixture: str) -> None:
    """Every ``plots/data/*.tsv`` (empty states included) is linked from ``report.html``."""
    page_, tmp = sim_report(fixture)
    links = {a.attrs.get("href", "") for a in page_.root.by_tag("a")}
    tsvs = sorted((tmp / "plots" / "data").glob("*.tsv"))
    assert tsvs, f"{fixture}: no data TSV rendered"
    missing = [p.name for p in tsvs if f"../plots/data/{p.name}" not in links]
    assert not missing, missing


# ---------------------------------------------------------------- M3


def _nfc(text: str) -> str:
    return unicodedata.normalize("NFC", words(text))


def _glossary(page_: Page) -> dict[str, str]:
    found = page_.root.by_attr("id", "glossary")
    assert len(found) == 1, f"{len(found)} #glossary"
    out: dict[str, str] = {}
    term: str | None = None
    for child in found[0].children:
        if isinstance(child, str):
            continue
        if child.tag == "dt":
            term = _nfc(child.text())
        elif child.tag == "dd" and term is not None:
            out[term] = _nfc(child.text())
    return out


def test_glossary_block_length_is_l_prime(tmp_path: Path) -> None:
    gloss = _glossary(build("P", tmp_path).page())
    term = _nfc("block length (ℓ′)")
    assert term in gloss, sorted(gloss)
    assert _nfc("block length (ℓ̂)") not in gloss, "the old term (ℓ̂) is still listed"
    dd = gloss[term]
    l_hat = _nfc("ℓ̂")
    assert l_hat in dd, dd
    assert "estimated correlation length" in dd, dd
    assert _nfc(f"min(20,2{l_hat}+1,numberofTUs)") in dd.replace(" ", ""), dd


# ---------------------------------------------------------------- L1

_TU_ERROR = "tu_block_tu_error"
_DISTANCE = "tu_source_distance_model"


def _autocorr_report(tmp: Path, source: str) -> Page:
    out = run_out(tmp)
    ac_flags = _TU_ERROR if source != "distance_model" else f"{_TU_ERROR};{_DISTANCE}"
    ac = [{"contig": "chr1", "statistic": "morans_i", "lag": lag, "null_model": "tu_block",
           "value": 0.1, "stat": "morans_i", "p_value": 0.01, "q_value": 0.02,
           "method": "tu_block_permutation", "n_perm": 999, "flags": ac_flags}
          for lag in (1, 2)]  # fmt: skip
    write_schema_table(out / "tables" / "autocorrelation.tsv", "autocorrelation", ac)
    write_pairs(out)
    write_manifest(out, manifest(operon_source=source, operon_source_reason="UniOP not installed"))
    return page(out)


def _rows_clause(flag: str) -> re.Pattern[str]:
    from gippa.viz.words import plain

    text = re.escape(plain(flag))
    return re.compile(rf"\(these rows carry the flag ['‘]{text}['’] in the table\)")


def _caveat(page_: Page, start: str) -> str:
    sec = page_.section("genome-context")
    found = [words(p.text()) for p in sec.by_tag("p") if start in words(p.text())]
    assert len(found) == 1, (start, found)
    return found[0]


@pytest.mark.parametrize("source", ["distance_model", "uniop"])
def test_autocorrelation_caveats_name_their_rows(source: str, tmp_path: Path) -> None:
    page_ = _autocorr_report(tmp_path, source)
    tu_error = _caveat(page_, "assumes exact TU boundaries")
    assert _rows_clause(_TU_ERROR).search(tu_error), tu_error
    checked = [tu_error]
    if source == "distance_model":
        distance = _caveat(page_, "when TUs come from the distance model")
        assert _rows_clause(_DISTANCE).search(distance), distance
        checked.append(distance)
    for text in checked:
        assert not flag_tokens(text), flag_tokens(text)
        assert _TU_ERROR not in text and _DISTANCE not in text, text


# ---------------------------------------------------------------- L2

#: plot -> (index stem, a drawn-instance stem prefix)
INDEX_STEMS = {"V03": "V03_locus", "V13": "V13_novel_orfs_locus"}


def _captions_by_stem(page_: Page, pid: str) -> dict[str, str]:
    out = {}
    for cap in page_.slot(pid).by_tag("figcaption"):
        labels = " ".join(words(a.text()) for a in cap.by_tag("a"))
        text = words(cap.text())
        out[cap.attrs.get("data-stem", "")] = text[: len(text) - len(labels)].strip()
    return out


def _plot_files(pid: str, suffixes: list[str | None], tmp: Path) -> Any:
    from helpers.report_fixtures import figure, frame

    from gippa.io.plot_files import read_plot_files, write_instance
    from gippa.viz import base

    for suffix in suffixes:
        inst = base.Instance(suffix=suffix, figure=figure(), data=frame(),
                             title=f"{pid} {suffix}")  # fmt: skip
        write_instance(inst, pid, NAMES[pid], tmp, ["svg"], 50)
    return read_plot_files(tmp)[pid]


@pytest.mark.parametrize(("pid", "suffixes"), [("V03", ["ISL0001", "PRO0002", None]),
                                               ("V13", [None, "locus_N1", "locus"])])  # fmt: skip
def test_index_instance_has_its_own_caption(pid: str, suffixes: list[str | None],
                                            tmp_path: Path) -> None:  # fmt: skip
    from gippa.report.content import CAPTIONS

    run = build("P", tmp_path)
    if pid == "V13":
        from helpers.report_content import write_novel

        write_novel(run.out, [("N1", True, "0.01")])
    dress(run, tmp_path)
    run.plots[pid] = _plot_files(pid, suffixes, tmp_path / f"multi-{pid}")
    caps = _captions_by_stem(run.page(), pid)
    index = INDEX_STEMS[pid]
    assert index in caps, sorted(caps)
    others = {s: c for s, c in caps.items() if s != index}
    assert others, caps
    for stem, cap in others.items():
        assert cap == words(CAPTIONS[pid]), (stem, cap)
    idx = caps[index]
    assert idx.startswith(QUESTIONS[pid]), idx
    assert re.search(r"\bindex\b", idx, re.I), idx
    assert idx != words(CAPTIONS[pid]), "the index instance reuses the drawn caption"
    assert not flag_tokens(idx), idx


def test_sim_plots_v03_index_caption() -> None:
    page_, _ = sim_report("sim_plots")
    caps = _captions_by_stem(page_, "V03")
    assert INDEX_STEMS["V03"] in caps, sorted(caps)
    idx = caps.pop(INDEX_STEMS["V03"])
    assert caps, "sim_plots has no drawn V03 locus"
    assert idx.startswith(QUESTIONS["V03"]) and re.search(r"\bindex\b", idx, re.I), idx
    assert all(c != idx for c in caps.values()), caps


# ---------------------------------------------------------------- L3


def test_not_registered_reads_as_plain_words(tmp_path: Path) -> None:
    run = dress(build("P", tmp_path), tmp_path)
    run.plots.pop("V15", None)
    run.skipped["V15"] = "not_registered"
    page_ = run.page()
    note = words(page_.slot("V15").text())
    assert "this plot is not available in this gippa version" in note.lower(), note
    assert "not_registered" not in words(page_.root.text()), "raw token in the page"


def test_not_registered_is_in_the_ac21_vocabulary() -> None:
    from helpers.m12_04_pr2 import flag_vocabulary

    assert "not_registered" in flag_vocabulary()
    assert flag_tokens("V15: not shown: not_registered") == ["not_registered"]


# ---------------------------------------------------------------- L4


def _write_sensitivity_raw(out: Path, header: list[str], rows: list[list[str]]) -> None:
    path = out / "tables" / "sensitivity.tsv"
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["\t".join(header)] + ["\t".join(r) for r in rows]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _sens_page(tmp: Path, rows: list[tuple[str, str, Any, Any, Any, Any]]) -> Page:
    run = build("P", tmp)
    write_sensitivity(run.out, rows)  # type: ignore[arg-type]
    return run.page()


def test_sensitivity_na_counts_keep_the_tables(tmp_path: Path) -> None:
    page_ = _sens_page(tmp_path, [("islands", "block_x2", None, None, None, None),
                                  ("operons", "pi0_x", 0, 0, 0, None)])  # fmt: skip
    assert SENSITIVITY_NONE not in words(page_.html)
    assert page_.table("sensitivity") is not None


def test_sensitivity_without_count_columns_does_not_crash(tmp_path: Path) -> None:
    run = build("P", tmp_path)
    _write_sensitivity_raw(run.out, ["family", "variant", "jaccard"],
                           [["islands", "block_x2", "NA"]])  # fmt: skip
    page_ = run.page()
    assert page_.section("sensitivity") is not None
    assert SENSITIVITY_NONE not in words(page_.html)


_COUNT = st.one_of(st.none(), st.integers(0, 3))


@settings(max_examples=25, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(st.lists(st.tuples(st.sampled_from(["islands", "operons", "enrichment"]), _COUNT, _COUNT),
                min_size=1, max_size=4))  # fmt: skip
def test_sensitivity_sentence_iff_every_count_is_zero(
    rows: list[tuple[str, int | None, int | None]],
) -> None:
    body = []
    for k, (fam, a, b) in enumerate(rows):
        shared = None if a is None or b is None else min(a, b)
        jac = None if not a or not b or shared is None else shared / (a + b - shared)
        body.append((fam, f"v{k}", a, b, shared, jac))
    with tempfile.TemporaryDirectory() as tmp:
        page_ = _sens_page(Path(tmp), body)
    want = all(a == 0 and b == 0 for _, a, b in rows)
    assert (SENSITIVITY_NONE in words(page_.html)) == want, rows
    assert (page_.table("sensitivity") is None) == want, rows


# ---------------------------------------------------------------- L5


def test_proteins_only_novel_orfs_reason_needs_peptides(tmp_path: Path) -> None:
    entries = dict(not_applicable_entries(dress(build("R", tmp_path), tmp_path).page()))
    assert "Novel ORFs" in entries, entries
    assert "novel ORF calling needs peptides" in entries["Novel ORFs"], entries


def test_peptide_run_novel_orfs_reason_does_not_say_needs_peptides(tmp_path: Path) -> None:
    entries = dict(not_applicable_entries(dress(build("P", tmp_path), tmp_path).page()))
    assert "Novel ORFs" in entries, entries
    assert "needs peptides" not in entries["Novel ORFs"], entries


# ---------------------------------------------------------------- L6

#: fixture -> the sections listed, from the rules (module docstring). No PR 1 fixture writes
#: pairs, operon flags or novel ORFs; R and D have no variants table; P_core has bgc off; every
#: fixture's V08 is drawn (``dress``: genes carry COG ``KL``), so Function is never listed.
LISTED: dict[str, set[str]] = {
    "P": {"Operons", "Novel ORFs"},
    "R": {"Operons", "Variants", "Novel ORFs"},
    "D": {"Operons", "Variants", "Novel ORFs"},
    "P_nokegg": {"Operons", "Novel ORFs"},
    "P_core": {"Operons", "Novel ORFs", BGC},
    "decoy": {"Operons", "Novel ORFs"},
    "both": {"Operons", "Novel ORFs"},
}


def test_listed_table_covers_every_fixture() -> None:
    assert set(LISTED) == set(FIXTURES)


@pytest.mark.parametrize("name", FIXTURES)
def test_exact_listed_set_per_fixture(name: str, tmp_path: Path) -> None:
    from gippa.report.build import SECTIONS

    page_ = dress(build(name, tmp_path), tmp_path).page()
    assert _listed(page_) == LISTED[name], (name, not_applicable_entries(page_))
    heads = {s.heading: s.id for s in SECTIONS}
    rendered = {s.attrs.get("id") for s in page_.sections()}
    missing = {s.id for s in SECTIONS} - rendered
    assert missing == {heads[n] for n in LISTED[name] if n in heads}, (name, missing)


def test_exact_listed_set_run_out(tmp_path: Path) -> None:
    """A bare ``run_out`` (default modules, no tables, no plots)."""
    assert _listed(page(run_out(tmp_path))) == {"Operons", "Variants", "Novel ORFs"}


def test_sim_reports_list_exactly_the_missing_sections() -> None:
    from gippa.report.build import SECTIONS

    heads = {s.heading: s.id for s in SECTIONS}
    for fixture in ("sim_plots", "sim_empty"):
        page_, _ = sim_report(fixture)
        listed = {heads[n] for n, _ in not_applicable_entries(page_) if n in heads}
        rendered = {s.attrs.get("id") for s in page_.sections()}
        assert {s.id for s in SECTIONS} - rendered == listed, fixture
