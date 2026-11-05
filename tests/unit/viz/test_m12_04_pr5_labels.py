"""M12-04 PR 5, P-13: every abundance axis names the run's measure (docs/specs/
M12-04-outputs-revision.md § PR 5, AC 5.4; ``gippa.viz.style.abundance_label`` from PR 0).

Fixtures: ``sim_plots`` (iBAQ, molar) and a copy of it re-labelled MaxLFQ (relative), so a
hard-coded "log2 iBAQ (molar)" fails the second case. The V03/V13 colour bar is tested in
``test_m12_04_pr5_locus_style.py``.

Interface assumed (test-author, PR 5):

- the measure comes from ``genes.tsv`` ``abundance_type`` / ``abundance_scale`` (as the report's
  run box does, PR 1), so V11 and V15 read ``genes_tsv`` too (optional input);
- an axis that shows plain abundance equals the label exactly: V02 (every replicon panel's y
  axis) and V08 (the y axis, abundance mode);
- an axis that qualifies abundance keeps its qualifier and contains the label: V11 panel (c)
  (upstream / downstream partner) and V15 panel (d) (predicted / observed); the bare
  "log2 abundance" phrase is gone from those axes.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest
from helpers import m12_04_pr5 as h

from gippa.viz import style

PIDS = ("V02", "V08", "V11", "V15")
CASES = {"ibaq": ("ibaq", "molar"), "maxlfq": ("maxlfq", "relative")}


@pytest.fixture(scope="module", params=sorted(CASES))
def case(request: pytest.FixtureRequest,
         tmp_path_factory: pytest.TempPathFactory) -> tuple[str, dict[str, Any]]:  # fmt: skip
    atype, scale = CASES[request.param]
    if request.param == "ibaq":
        out = h.SIM_OUT
    else:
        out = h.sim_copy(tmp_path_factory.mktemp("maxlfq"))
        h.set_abundance(out, atype, scale)
    label = style.abundance_label(atype, scale)
    assert label is not None
    return label, h.render(out, PIDS)


def test_v02_panels_use_the_abundance_label(case: tuple[str, dict[str, Any]]) -> None:
    label, results = case
    inst = h.only(results, "V02")
    ylabels = [ax.get_ylabel() for ax in inst.figure.axes if ax.get_ylabel()]
    assert ylabels, "V02 has no labelled panel"
    assert set(ylabels) == {label}


def test_v08_axis_uses_the_abundance_label(case: tuple[str, dict[str, Any]]) -> None:
    label, results = case
    inst = h.only(results, "V08")
    (ax,) = [a for a in inst.figure.axes if a.get_xlabel() == "COG category"]
    assert ax.get_ylabel() == label


@pytest.mark.parametrize(("pid", "panel"), [("V11", "c"), ("V15", "d")])
def test_qualified_abundance_axes_contain_the_label(
    case: tuple[str, dict[str, Any]], pid: str, panel: str
) -> None:
    label, results = case
    inst = h.only(results, pid)
    ax = h.axes(inst.figure, panel)
    for text in (ax.get_xlabel(), ax.get_ylabel()):
        assert label in text, (pid, panel, text)
        assert "log2 abundance" not in text, (pid, panel, text)


@pytest.mark.parametrize(("name", "want"), [("P", "log2 iBAQ (molar)"),
                                            ("R", "log2 MaxLFQ (relative)")])  # fmt: skip
def test_spec_examples_a_and_b_on_the_pr1_fixtures(tmp_path: Path, name: str, want: str) -> None:
    """Spec examples: A "log2 iBAQ (molar)", B "log2 MaxLFQ (relative)" (fixtures P, R)."""
    from unit.report.test_m12_04_pr1_fixtures import build

    results = h.render(build(name, tmp_path).out, ["V02", "V08"])
    v02 = h.only(results, "V02")
    assert {ax.get_ylabel() for ax in v02.figure.axes if ax.get_ylabel()} == {want}
    (ax,) = [a for a in h.only(results, "V08").figure.axes if a.get_xlabel() == "COG category"]
    assert ax.get_ylabel() == want
