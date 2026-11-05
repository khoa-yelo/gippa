"""M12-04 PR 13, AC 13.3: V13 locus views are linked, not inlined, when the report would reach
the link threshold (docs/specs/M12-04-outputs-revision.md § PR 13 "Report size").

The threshold is a named constant: 20 MB, 80 % of the 25 MB §10.3 target
(``gippa.steps.plots.REPORT_BUDGET``). When the page with every SVG inlined would be at or over
it, each ``V13_novel_orfs_locus_*`` instance is an ``<img>`` (and an ``<a>``) whose relative
``src`` names its file under ``plots/``; small SVGs stay inlined. Synthetic padded SVGs make the
sizes; no real run is used (the K-12 report on main is 19.8 MB and stays inlined).

Interface pinned here (PR 13 test author, beyond the spec text):
- the constant is ``gippa.report.build.REPORT_LINK_THRESHOLD`` (bytes), read at call time, so a
  test may monkeypatch it down;
- the rule lives in ``build_report``: the size compared is the UTF-8 byte length of the page
  ``build_report`` returns with every SVG inlined (the twins, when given, included); at or over
  the threshold (``>=``) the locus views are linked;
- a linked view: an ``<img>`` in the V13 slot (``data-plot="V13"``) whose ``src`` is relative
  (no scheme, no leading ``/``), ends in ``plots/<stem>.svg`` and resolves from ``OUT/report/``
  to the file; an ``<a href>`` to the same target exists in the slot; the SVG text is not in the
  page. Other plots (here V17) stay inlined.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest
from helpers.report_fixtures import Node, cfg, parse, run_out
from helpers.viz_out_tree import write_schema_table

MB = 1_000_000
LOCUS = ("V13_novel_orfs_locus_NOVEL_chr1_44225_m", "V13_novel_orfs_locus_NOVEL_chr1_45143_m")
MAIN_V13 = "V13_novel_orfs"
V17 = "V17_qc"
LOCUS_SRC = re.compile(r"^(?:\.\./)?plots/V13_novel_orfs_locus_[^/\"]+\.svg$")


def _svg(stem: str, pad: int) -> str:
    return ('<svg xmlns="http://www.w3.org/2000/svg" width="10" height="10">'
            f"<title>{stem}: a finding</title>"
            f'<desc id="pad">{"q" * pad}</desc>'
            '<rect id="r" width="1" height="1"/></svg>\n')  # fmt: skip


def _out(tmp: Path, locus_bytes: int) -> Path:
    """``OUT/`` with a header-only ``novel_orfs.tsv`` (the section applies), two V13 locus SVGs
    padded to ``locus_bytes``, a small V13 main plot and a small V17."""
    out = run_out(tmp)
    write_schema_table(out / "tables" / "novel_orfs.tsv", "novel_orfs", [])
    plots = out / "plots"
    plots.mkdir()
    for stem in LOCUS:
        (plots / f"{stem}.svg").write_text(_svg(stem, locus_bytes), encoding="utf-8")
    for stem in (MAIN_V13, V17):
        (plots / f"{stem}.svg").write_text(_svg(stem, 10), encoding="utf-8")
    return out


def _html(out: Path) -> str:
    from gippa.io.plot_files import read_plot_files
    from gippa.io.run_view import RunView
    from gippa.report.build import build_report

    return build_report(RunView.from_outdir(out), read_plot_files(out / "plots"), {}, cfg())


def inlined(html: str, stem: str) -> bool:
    """The instance's SVG is in the page (``embed_svg`` prefixes its ids with the stem)."""
    return f'id="{stem}-pad"' in html


def _slot(html: str) -> Node:
    found = parse(html).root.by_attr("data-plot", "V13")
    assert found, "no V13 slot"
    return found[0]


def assert_linked(html: str, out: Path) -> None:
    slot = _slot(html)
    srcs = [n.attrs.get("src", "") for n in slot.by_tag("img")]
    hrefs = {n.attrs.get("href", "") for n in slot.by_tag("a")}
    for stem in LOCUS:
        assert not inlined(html, stem), f"{stem} is still inlined"
        mine = [s for s in srcs if s.endswith(f"plots/{stem}.svg")]
        assert len(mine) == 1, f"{stem}: img src {srcs}"
        src = mine[0]
        assert LOCUS_SRC.match(src), src
        assert (out / "report" / src).resolve() == (out / "plots" / f"{stem}.svg").resolve()
        assert (out / "report" / src).is_file(), src
        assert src in hrefs, f"no <a href={src!r}> in the V13 slot"
    assert "q" * 100_000 not in html


# ---------------------------------------------------------------- the constant


def test_threshold_is_a_named_constant_of_20_mb() -> None:
    from gippa.report import build
    from gippa.steps.plots import REPORT_BUDGET

    assert build.REPORT_LINK_THRESHOLD == 20_000_000
    assert int(0.8 * REPORT_BUDGET) == build.REPORT_LINK_THRESHOLD


# ---------------------------------------------------------------- AC 13.3


def test_ac_13_3_large_svgs_are_linked(tmp_path: Path) -> None:
    """Two 10.5 MB locus SVGs push the inlined page over 20 MB: both are linked, V17 stays
    inlined, and the page is well under the threshold."""
    from gippa.report.build import REPORT_LINK_THRESHOLD

    out = _out(tmp_path, int(10.5 * MB))
    html = _html(out)
    assert_linked(html, out)
    assert inlined(html, V17)
    assert len(html.encode("utf-8")) < REPORT_LINK_THRESHOLD


def test_ac_13_3_small_svgs_stay_inlined(tmp_path: Path) -> None:
    out = _out(tmp_path, 1_000)
    html = _html(out)
    for stem in (*LOCUS, V17):
        assert inlined(html, stem), f"{stem} is not inlined"
    assert not [n for n in _slot(html).by_tag("img")
                if "V13_novel_orfs_locus_" in n.attrs.get("src", "")]  # fmt: skip


def test_threshold_boundary_is_at_or_over(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """With the constant patched: the inlined page's own size links (``>=``); one byte more
    keeps it inlined."""
    from gippa.report import build

    out = _out(tmp_path, 50_000)
    monkeypatch.setattr(build, "REPORT_LINK_THRESHOLD", 10**12)
    full = _html(out)
    assert all(inlined(full, s) for s in LOCUS)
    size = len(full.encode("utf-8"))
    monkeypatch.setattr(build, "REPORT_LINK_THRESHOLD", size + 1)
    assert _html(out) == full
    monkeypatch.setattr(build, "REPORT_LINK_THRESHOLD", size)
    assert_linked(_html(out), out)


def test_patched_threshold_links_moderate_svgs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from gippa.report import build

    out = _out(tmp_path, 200_000)
    monkeypatch.setattr(build, "REPORT_LINK_THRESHOLD", 300_000)
    html = _html(out)
    assert_linked(html, out)
    assert inlined(html, V17)


def test_linked_report_on_disk(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Through ``gippa.steps.plots.write_report``: the written page links its locus views and
    every relative link it holds names a file that exists."""
    from gippa.io.run_view import RunView
    from gippa.report import build
    from gippa.steps.plots import write_report

    out = _out(tmp_path, 200_000)
    monkeypatch.setattr(build, "REPORT_LINK_THRESHOLD", 300_000)
    page = write_report(RunView.from_outdir(out), out / "plots", out / "report", cfg())
    html = page.read_text(encoding="utf-8")
    assert_linked(html, out)
    for src in [n.attrs["src"] for n in parse(html).root.by_tag("img") if "src" in n.attrs]:
        assert (out / "report" / src).is_file(), src
