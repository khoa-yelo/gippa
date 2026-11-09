"""M12-05a AC 10 and AC 11: the report shows the V13 locus views in display order, and the V13
captions say that every candidate is drawn (docs/specs/M12-05a-novel-orf-plot.md § 5, § 6.3).

Interface pinned here (PR A test author, beyond the spec text):

- ``gippa.report.content.display_order(stems, index_rows) -> list[str]`` is pure. ``stems`` are
  the V13 instance stems; ``index_rows`` are the rows of ``plots/data/V13_novel_orfs_locus.tsv``
  as mappings from column name to cell text (``read_tsv`` cells, NA as ``None``; booleans
  ``"true"`` / ``"false"``), in TSV order. It returns panel (a) (``V13_novel_orfs``), the index
  (``V13_novel_orfs_locus``), then ``V13_novel_orfs_locus_<orf_id>`` for each ``rendered = true``
  row in ``order``, then any other stem in alphabetical order;
- ``build_report`` applies it in the V13 slot: the ``<figcaption data-stem>`` sequence of the
  slot follows it, inlined or linked (PR 13).
"""

from __future__ import annotations

import random
from pathlib import Path

import pytest
from helpers import m12_05a as h
from helpers.m12_04_pr2 import flag_tokens, internal_refs
from helpers.report_fixtures import cfg, parse, run_out
from helpers.viz_out_tree import write_schema_table

from gippa.io.tsv import write_tsv

MAIN = "V13_novel_orfs"
INDEX = "V13_novel_orfs_locus"
PREFIX = "V13_novel_orfs_locus_"


def _index_rows(ids: list[str]) -> list[dict[str, str | None]]:
    rows = {str(r["orf_id"]): r for r in h.route_e_rows()}
    out = []
    for k, oid in enumerate(ids, start=1):
        r = rows[oid]
        out.append({"order": str(k), "orf_id": oid, "class": str(r["class"]), "called": "false",
                    "class_fdr": str(r["class_fdr"]), "best_peptide_q": str(r["best_peptide_q"]),
                    "n_peptides": "1", "rank": None, "rendered": "true"})  # fmt: skip
    return out


def _e1_stems() -> list[str]:
    return sorted([MAIN, INDEX, *(PREFIX + i for i in h.E1_ORDER)])


def _want() -> list[str]:
    return [MAIN, INDEX, *(PREFIX + i for i in h.E1_ORDER)]


# --- AC 10: content.display_order -----------------------------------------------------------------


def test_e1_alphabetical_stems_are_not_display_order() -> None:
    """Guard on the fixture: alphabetical order would put ``…130563_m`` first, ``…20324_m``
    before ``…20948_p``; only the index order gives E1."""
    assert _e1_stems() != _want()


def test_display_order_e1() -> None:
    from gippa.report.content import display_order

    assert display_order(_e1_stems(), _index_rows(h.E1_ORDER)) == _want()


def test_display_order_ignores_the_input_stem_order() -> None:
    from gippa.report.content import display_order

    stems = _e1_stems()
    random.Random(5).shuffle(stems)
    assert display_order(stems, _index_rows(h.E1_ORDER)) == _want()


def test_display_order_puts_stale_stems_last_alphabetically() -> None:
    from gippa.report.content import display_order

    stale = [PREFIX + "NOVEL_zzz_1_p", PREFIX + "NOVEL_aaa_1_p"]
    stems = sorted([*_e1_stems(), *stale])
    assert display_order(stems, _index_rows(h.E1_ORDER)) == [*_want(), *sorted(stale)]


def test_display_order_skips_rows_not_rendered() -> None:
    from gippa.report.content import display_order

    rows = _index_rows(h.E1_ORDER)
    rows[0]["rendered"] = "false"  # its file is stale: it goes with the unlisted stems
    first = PREFIX + h.E1_ORDER[0]
    want = [MAIN, INDEX, *(PREFIX + i for i in h.E1_ORDER[1:]), first]
    assert display_order(_e1_stems(), rows) == want


def test_display_order_without_index_rows_keeps_the_input_order() -> None:
    from gippa.report.content import display_order

    stems = _e1_stems()
    assert display_order(stems, []) == stems
    assert display_order(list(reversed(stems)), []) == list(reversed(stems))


def test_display_order_is_pure() -> None:
    from gippa.report.content import display_order

    stems, rows = _e1_stems(), _index_rows(h.E1_ORDER)
    before = (list(stems), [dict(r) for r in rows])
    assert display_order(stems, rows) == display_order(stems, rows)
    assert (stems, rows) == before


# --- AC 10: the V13 slot of build_report ----------------------------------------------------------


def _svg(stem: str, pad: int) -> str:
    return ('<svg xmlns="http://www.w3.org/2000/svg" width="10" height="10">'
            f"<title>{stem}: a finding</title>"
            f'<desc id="pad">{"q" * pad}</desc>'
            '<rect id="r" width="1" height="1"/></svg>\n')  # fmt: skip


def _out(tmp: Path, *, index: bool = True, pad: int = 10) -> Path:
    """An E1-like ``OUT/``: panel (a), the index and the 8 route E views as small SVGs, plus the
    index data TSV in E1 order (``plots/data/V13_novel_orfs_locus.tsv``)."""
    out = run_out(tmp)
    write_schema_table(out / "tables" / "novel_orfs.tsv", "novel_orfs", h.route_e_rows())
    plots = out / "plots"
    (plots / "data").mkdir(parents=True)
    (out / "report").mkdir(exist_ok=True)
    for stem in _e1_stems():
        size = pad if stem.startswith(PREFIX) else 10
        (plots / f"{stem}.svg").write_text(_svg(stem, size), encoding="utf-8")
    if index:
        rows = _index_rows(h.E1_ORDER)
        write_tsv(plots / "data" / f"{INDEX}.tsv", h.INDEX_COLUMNS,
                  [[r[c] for c in h.INDEX_COLUMNS] for r in rows])  # fmt: skip
    return out


def _slot_stems(out: Path) -> list[str]:
    return _slot(out)[0]


def _slot(out: Path) -> tuple[list[str], str]:
    from gippa.io.plot_files import read_plot_files
    from gippa.io.run_view import RunView
    from gippa.report.build import build_report

    html = build_report(RunView.from_outdir(out), read_plot_files(out / "plots"), {}, cfg())
    slots = parse(html).root.by_attr("data-plot", "V13")
    assert slots, "no V13 slot"
    caps = slots[0].by_tag("figcaption")
    return [n.attrs["data-stem"] or "" for n in caps if "data-stem" in n.attrs], html


def test_report_v13_slot_follows_the_index_order(tmp_path: Path) -> None:
    assert _slot_stems(_out(tmp_path)) == _want()


def test_report_v13_slot_without_index_tsv_is_alphabetical(tmp_path: Path) -> None:
    assert _slot_stems(_out(tmp_path, index=False)) == _e1_stems()


def test_report_v13_linked_views_follow_the_index_order(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from gippa.report import build

    out = _out(tmp_path, pad=50_000)
    monkeypatch.setattr(build, "REPORT_LINK_THRESHOLD", 200_000)
    stems, html = _slot(out)
    for oid in h.E1_ORDER:  # linked, not inlined (embed_svg prefixes the ids with the stem)
        assert f'id="{PREFIX}{oid}-pad"' not in html, oid
    assert stems == _want()


# --- AC 11: captions ------------------------------------------------------------------------------


@pytest.mark.parametrize("stem", [MAIN, PREFIX + h.E1_ORDER[0]])
def test_v13_caption(stem: str) -> None:
    from gippa.report.content import CAPTIONS, caption

    assert CAPTIONS["V13"] == h.CAPTION_V13
    assert caption("V13", stem) == h.CAPTION_V13


def test_v13_index_caption() -> None:
    from gippa.report.content import INDEX_CAPTIONS, caption

    assert INDEX_CAPTIONS["V13"] == h.INDEX_CAPTION_V13
    assert caption("V13", INDEX) == f"{h.QUESTION_V13} {h.INDEX_CAPTION_V13}"


@pytest.mark.parametrize("text", [h.CAPTION_V13, h.INDEX_CAPTION_V13])
def test_v13_captions_pass_the_prose_checks(text: str) -> None:
    """Guard on the pinned strings (M12-04 rule 2): no ADR, §, flag token, path or raw class."""
    assert not internal_refs(text), internal_refs(text)
    assert not flag_tokens(text), flag_tokens(text)
    assert "tables/" not in text and "_" not in text
    assert not h.raw_class_tokens(text)
