"""M12-04 PR 2 behaviour 5, AC 2.8: an empty-state instance (``Instance.empty``, persisted as
``data-empty="true"`` on its SVG root and read back as ``InstanceFiles.empty``) is shown as a
``<p class="notice">`` holding its title, not as an embedded image
(docs/specs/M12-04-outputs-revision.md § PR 2; lead decision 5).

Interface assumed (PR 2 test author): the notice sits inside the plot's slot
(``div[data-plot="Vnn"]``); a slot whose every instance is empty holds no ``<svg>``, no
``<img>`` and no Plotly twin (``.twin``); a drawn instance in the same slot keeps its SVG. The
image files stay in ``plots/`` (the plots step is unchanged).
"""

from __future__ import annotations

from pathlib import Path

from helpers.m12_04_pr2 import dress, sim_report, words
from helpers.report_fixtures import Node
from unit.report.test_m12_04_pr1_fixtures import TITLE_V02, build

TITLE_V15 = "Too few genes for the determinants model (99 rows, 34 TUs)"


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


def _notices(slot: Node) -> list[str]:
    return [words(n.text()) for n in slot.by_class("notice") if n.tag == "p"]


def test_empty_instance_is_a_notice(tmp_path: Path) -> None:
    run = dress(build("P", tmp_path), tmp_path, empty={"V02": TITLE_V02, "V15": TITLE_V15})
    page = run.page()
    for pid, title in (("V02", TITLE_V02), ("V15", TITLE_V15)):
        slot = page.slot(pid)
        assert any(title in n for n in _notices(slot)), (pid, _notices(slot))
        assert not _outer_svgs(slot), f"{pid}: an <svg> is embedded for an empty state"
        assert not slot.by_tag("img"), f"{pid}: an <img> for an empty state"
        assert not slot.by_class("twin"), f"{pid}: a Plotly twin for an empty state"
    # the control: a drawn plot keeps its SVG and gets no notice
    v01 = page.slot("V01")
    assert _outer_svgs(v01) and not _notices(v01)


def test_empty_files_are_still_written(tmp_path: Path) -> None:
    run = dress(build("P", tmp_path), tmp_path, empty={"V02": TITLE_V02})
    files = run.plots["V02"].instances
    assert files and all(i.empty for i in files)
    assert (tmp_path / "plots-P" / files[0].images[0]).is_file()


def test_sim_empty_empty_states_are_notices() -> None:
    """sim_empty: every instance whose SVG carries ``data-empty="true"``."""
    from gippa.io.plot_files import read_plot_files, svg_title

    page, tmp = sim_report("sim_empty")
    plots = read_plot_files(tmp / "plots")
    seen = 0
    for pid, pf in plots.items():
        empty = [i for i in pf.instances if i.empty]
        if not empty:
            continue
        slots = page.root.by_attr("data-plot", pid)
        if not slots:
            continue  # a plot without a report slot (or a section moved to not-applicable)
        slot = slots[0]
        notices = _notices(slot)
        for inst in empty:
            title = words(svg_title(inst.svg or "") or "")
            assert title and any(title in n for n in notices), (pid, title, notices)
            seen += 1
        drawn = len(pf.instances) - len(empty)
        assert len(_outer_svgs(slot)) == drawn, (pid, drawn, len(_outer_svgs(slot)))
        if not drawn:
            assert not slot.by_tag("img") and not slot.by_class("twin"), pid
    assert seen >= 5, seen
