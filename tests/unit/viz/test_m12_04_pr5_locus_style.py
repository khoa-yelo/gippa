# ruff: noqa: RUF002  (figure text is Unicode)
"""M12-04 PR 5, P-01: undetected genes in the V03 and V13 locus views are white with a
``#BBBBBB`` outline and a ``////`` hatch (docs/specs/M12-04-outputs-revision.md § PR 5,
AC 5.1; ADR-0075 item 4; Khoa's decision 4 of 2026-11-04). Also the locus colour bar half of
P-13 (AC 5.4).

Tested on the figure artists and the data TSVs, never on pixels.

Interface assumed (test-author, PR 5):

- ``gippa.viz.style.LOCUS_UNDETECTED`` is a mapping with exactly the keys ``facecolor``,
  ``edgecolor`` and ``hatch``, valued ``"white"``, ``"#BBBBBB"`` and ``"////"`` (the three
  matplotlib keyword names the patch takes).
- The locus data frame (``plots/data/V03_locus_<id>.tsv``, also V13's ``locus_<orf>``) gains a
  ``fill`` column: for every ``gene`` row the colour the arrow is filled with, as a matplotlib
  colour string (``"white"`` for undetected genes; a hex string otherwise, e.g. ``"#0072B2"`` for
  a detected gene in a run without abundance). Rows of other tracks may leave it NA. The column
  is written by ``write_instance`` like the others.
- An undetected gene is ``detection_status ∈ {not_detected, partial_excluded}``; the gene arrow
  is a closed polygon on the axes labelled ``"genes"`` (a patch or a patch-collection member)
  spanning exactly the gene's window-frame ``[start, end]`` (origin-crossing genes keep the M9
  unwrapped span, drawn as before).
- The locus key has an entry labelled exactly ``"not detected"`` whose patch is white, edged
  ``#BBBBBB`` and hatched ``////`` (it replaces the grey "undetected" swatch).
- A run without abundance (every ``abundance_log2`` NA, ``abundance_type = none``): no colour
  bar (no ``<colorbar>`` axes), and detected genes are filled ``#0072B2``.
"""

from __future__ import annotations

import tempfile
from pathlib import Path
from typing import Any

import matplotlib
import matplotlib.colors as mcolors
import pandas as pd
import pytest
from helpers import m12_04_pr5 as h
from helpers import viz_locus as vl
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from gippa.config import GippaConfig, resolve_config
from gippa.io.plot_files import write_instance
from gippa.io.run_view import RunView
from gippa.io.tsv import read_tsv
from gippa.viz import style
from gippa.viz.base import Instance
from gippa.viz.locus import render_locus, select_loci

WHITE, GREY, ORANGE, BLUE = "#FFFFFF", "#BBBBBB", "#E69F00", "#0072B2"
HATCH = "////"
UNDETECTED = ("not_detected", "partial_excluded")
STATUSES = ("detected", "not_detected", "partial_excluded", "single_peptide_unestimated")
SCALE = (18.0, 24.0)
CIVIDIS = matplotlib.colormaps["cividis"]
LUT = {h.hexc(CIVIDIS(i)) for i in range(CIVIDIS.N)}


@pytest.fixture(autouse=True)
def _style() -> None:
    style.apply()


@pytest.fixture
def cfg() -> GippaConfig:
    return resolve_config()


def _render(out: Path, cfg: GippaConfig, locus_id: str = "prophage_1",
            scale: tuple[float, float] = SCALE) -> Instance:  # fmt: skip
    view = RunView.from_outdir(out)
    (lo,) = [x for x in select_loci(view, cfg) if x.locus_id == locus_id]
    return render_locus(lo, view, cfg, scale)


def _gene_row(inst: Instance, tag: str) -> Any:
    g = inst.data[(inst.data["track"] == "gene") & (inst.data["id"] == tag)]
    assert len(g) == 1, inst.data
    return g.iloc[0]


def _patch(inst: Instance, tag: str) -> h.GenePatch:
    row = _gene_row(inst, tag)
    found = h.patches_of(h.axes(inst.figure, "genes"), float(row["start"]), float(row["end"]))
    assert len(found) == 1, (tag, found)
    return found[0]


def _assert_undetected(p: h.GenePatch) -> None:
    assert p.face == WHITE, p
    assert p.edge == GREY, p
    assert p.hatch == HATCH, p
    if p.hatch_colour is not None:  # matplotlib >= 3.10 keeps a separate hatch colour
        assert p.hatch_colour == GREY, p


def _region_out(tmp: Path, genes: list[dict[str, object]]) -> Path:
    return vl.build_out(
        tmp,
        genes=genes,
        regions_stats=[vl.region_stat("prophage_1", "prophage", "chr1", 10_001, 20_000, 0.01)],
    )


# --- the style constant -------------------------------------------------------------------------


def test_locus_undetected_holds_the_three_values() -> None:
    assert dict(style.LOCUS_UNDETECTED) == {
        "facecolor": "white",
        "edgecolor": "#BBBBBB",
        "hatch": "////",
    }


# --- AC 5.1 on V03 --------------------------------------------------------------------------------


@pytest.mark.parametrize("status", UNDETECTED)
@pytest.mark.parametrize("strand", ["+", "-"])
def test_undetected_gene_is_white_grey_edged_and_hatched(
    tmp_path: Path, cfg: GippaConfig, status: str, strand: str
) -> None:
    """AC 5.1, and the minus-strand edge case (same rule, arrow to the left)."""
    out = _region_out(tmp_path, [vl.gene("U", "chr1", 12_001, 13_000, strand=strand,
                                         status=status)])  # fmt: skip
    inst = _render(out, cfg)
    _assert_undetected(_patch(inst, "U"))
    assert _gene_row(inst, "U")["fill"] == "white"


def test_detected_genes_keep_the_cividis_fill(tmp_path: Path, cfg: GippaConfig) -> None:
    """Example A: genes at 13–14 log2 keep their cividis colour next to undetected ones."""
    genes = [
        vl.gene("D1", "chr1", 11_001, 12_000, log2=13.2, abundance_scale="molar"),
        vl.gene("D2", "chr1", 12_101, 13_000, log2=14.0, strand="-", abundance_scale="molar"),
        vl.gene("U1", "chr1", 13_101, 14_000, status="not_detected"),
    ]
    inst = _render(_region_out(tmp_path, genes), cfg, scale=(10.0, 20.0))
    for tag, v in (("D1", 13.2), ("D2", 14.0)):
        p = _patch(inst, tag)
        assert p.face == h.hexc(CIVIDIS((v - 10.0) / 10.0)), (tag, p)
        assert p.face in LUT
        assert not p.hatch
        assert h.hexc(_gene_row(inst, tag)["fill"]) == p.face
    _assert_undetected(_patch(inst, "U1"))


def test_single_peptide_unestimated_keeps_the_orange_outline(
    tmp_path: Path, cfg: GippaConfig
) -> None:
    genes = [
        vl.gene("S", "chr1", 11_001, 12_000, status="single_peptide_unestimated"),
        vl.gene("U", "chr1", 12_101, 13_000, status="not_detected"),
    ]
    inst = _render(_region_out(tmp_path, genes), cfg)
    assert _patch(inst, "S").edge == ORANGE
    _assert_undetected(_patch(inst, "U"))


def test_the_legend_not_detected_entry_shows_the_new_style(
    tmp_path: Path, cfg: GippaConfig
) -> None:
    genes = [vl.gene("D", "chr1", 11_001, 12_000, log2=20.0),
             vl.gene("U", "chr1", 12_101, 13_000, status="not_detected")]  # fmt: skip
    inst = _render(_region_out(tmp_path, genes), cfg)
    entries = dict(h.legend_entries(inst.figure))
    assert "not detected" in entries, list(entries)
    assert h.patch_style(entries["not detected"]) == (WHITE, GREY, HATCH)
    assert "undetected" not in entries


def test_the_data_tsv_writes_the_fill_column(tmp_path: Path, cfg: GippaConfig) -> None:
    genes = [vl.gene("D", "chr1", 11_001, 12_000, log2=20.0),
             vl.gene("U", "chr1", 12_101, 13_000, strand="-", status="not_detected"),
             vl.gene("P", "chr1", 13_101, 14_000, status="partial_excluded")]  # fmt: skip
    inst = _render(_region_out(tmp_path / "run", genes), cfg)
    write_instance(inst, "V03", "locus", tmp_path / "plots", ["svg"], 100)
    # the project's own TSV reader: plot data TSVs carry no comment lines and hex fills
    # start with "#", so a comment-stripping reader would cut them
    header, rows = read_tsv(tmp_path / "plots" / "data" / "V03_locus_prophage_1.tsv")
    tsv = pd.DataFrame(rows, columns=header).fillna("")
    assert "fill" in tsv.columns, list(tsv.columns)
    by_tag = dict(zip(tsv.loc[tsv["track"] == "gene", "id"],
                      tsv.loc[tsv["track"] == "gene", "fill"], strict=True))  # fmt: skip
    assert by_tag["U"] == "white"
    assert by_tag["P"] == "white"
    assert by_tag["D"] not in ("white", "NA", "")
    assert h.hexc(by_tag["D"]) in LUT


# --- edge cases ---------------------------------------------------------------------------------


def test_no_abundance_detected_is_blue_and_there_is_no_colour_bar(
    tmp_path: Path, cfg: GippaConfig
) -> None:
    """Edge case C: no abundance; detected genes #0072B2, undetected white and hatched."""
    from helpers.viz_landscape import plot_spec

    genes = [vl.gene("D", "chr1", 11_001, 12_000),
             vl.gene("U", "chr1", 12_101, 13_000, status="not_detected")]  # fmt: skip
    out = _region_out(tmp_path, genes)
    result = plot_spec("V03").fn(RunView.from_outdir(out), cfg)
    (inst,) = [i for i in result if i.suffix == "prophage_1"]
    assert h.colorbar_labels(inst.figure) == []
    d = _patch(inst, "D")
    assert d.face == BLUE and not d.hatch
    assert h.hexc(_gene_row(inst, "D")["fill"]) == BLUE
    _assert_undetected(_patch(inst, "U"))


def test_origin_crossing_undetected_gene_keeps_its_split(tmp_path: Path, cfg: GippaConfig) -> None:
    genes = [
        vl.gene("GX", "chr1", 99_801, 100_200, status="not_detected", crosses_origin=True),
        vl.gene("G2", "chr1", 501, 900, log2=20.0),
    ]
    out = vl.build_out(tmp_path, genes=genes,
                       islands=[vl.island("ISL_1", "chr1", 99_001, 101_500, 0.01)])  # fmt: skip
    inst = _render(out, cfg, "ISL_1")
    row = _gene_row(inst, "GX")
    assert (int(row["start"]), int(row["end"])) == (99_800, 100_200)
    assert row["fill"] == "white"
    drawn = h.patches_of(h.axes(inst.figure, "genes"), 99_800, 100_200)
    assert drawn, h.gene_patches(h.axes(inst.figure, "genes"))
    for p in drawn:
        _assert_undetected(p)


# --- property: any detection_status mix -----------------------------------------------------------


@st.composite
def _genes(draw: Any) -> list[tuple[str, str, float | None]]:
    n = draw(st.integers(1, 8))
    return [(draw(st.sampled_from(STATUSES)), draw(st.sampled_from(["+", "-"])),
             draw(st.one_of(st.none(), st.floats(15.0, 27.0))))
            for _ in range(n)]  # fmt: skip


@settings(max_examples=25, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(spec=_genes())
def test_any_status_mix_gives_the_right_fill_per_gene(
    spec: list[tuple[str, str, float | None]],
) -> None:
    cfg = resolve_config()
    style.apply()
    genes = []
    for i, (status, strand, log2) in enumerate(spec):
        start = 11_001 + 1_000 * i
        value = log2 if status == "detected" else None
        genes.append(vl.gene(f"G{i}", "chr1", start, start + 799, strand=strand, status=status,
                             log2=value, abundance_type="ibaq",
                             abundance_scale="molar"))  # fmt: skip
    with tempfile.TemporaryDirectory() as tmp:
        inst = _render(_region_out(Path(tmp), genes), cfg)
        for i, (status, _strand, log2) in enumerate(spec):
            tag = f"G{i}"
            p, fill = _patch(inst, tag), _gene_row(inst, tag)["fill"]
            assert mcolors.is_color_like(fill), (tag, fill)
            assert h.hexc(fill) == p.face, (tag, fill, p)
            if status in UNDETECTED:
                _assert_undetected(p)
                assert fill == "white"
            elif status == "detected" and log2 is not None:
                assert p.face == h.hexc(CIVIDIS((min(max(log2, 18.0), 24.0) - 18.0) / 6.0))
                assert p.hatch != HATCH
            elif status == "single_peptide_unestimated":
                assert p.edge == ORANGE
                assert p.hatch != HATCH
        inst.figure.clear()


# --- V13 shares the locus renderer, and P-13 on the colour bar ------------------------------------


@pytest.fixture(scope="module")
def sim_locus(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Any]:
    """sim_plots with every third gene made undetected (sim_plots alone puts no undetected
    gene into a V13 window)."""
    out = h.sim_copy(tmp_path_factory.mktemp("sim_locus"))
    count = iter(range(10**6))

    def thin(r: dict[str, str]) -> dict[str, str]:
        if next(count) % 3 == 0:
            r.update(detection_status="not_detected", detected="false", abundance_log2="NA",
                     abundance_raw="NA")  # fmt: skip
        return r

    h.rewrite_table(out, "tables/genes.tsv", thin)
    return h.render(out, ["V03", "V13"])


def test_v13_locus_views_draw_undetected_genes_white_and_hatched(
    sim_locus: dict[str, Any],
) -> None:
    """AC 5.1 on V03 and V13 of the thinned sim_plots copy."""
    seen = 0
    for pid, inst in h.instances(sim_locus):
        if pid not in ("V03", "V13") or "track" not in inst.data.columns:
            continue
        assert "fill" in inst.data.columns, (pid, inst.suffix, list(inst.data.columns))
        ax = h.axes(inst.figure, "genes")
        genes = inst.data[inst.data["track"] == "gene"]
        for r in genes[genes["category"].isin(UNDETECTED)].itertuples(index=False):
            assert r.fill == "white", (pid, inst.suffix, r.id)
            found = h.patches_of(ax, float(r.start), float(r.end))
            assert found, (pid, inst.suffix, r.id)
            for p in found:
                _assert_undetected(p)
            seen += pid == "V13"
    assert seen > 0, "sim_plots should put an undetected gene into a V13 locus view"


def test_locus_colour_bars_use_the_abundance_label(sim_locus: dict[str, Any]) -> None:
    """AC 5.4 (V03 colour bar; V13 draws the same key): sim_plots is iBAQ (molar)."""
    want = style.abundance_label("ibaq", "molar")
    labels = [lab for pid, inst in h.instances(sim_locus) if pid in ("V03", "V13")
              for lab in h.colorbar_labels(inst.figure)]  # fmt: skip
    assert labels, "no locus colour bar drawn on sim_plots"
    assert set(labels) == {want}


def test_locus_colour_bar_label_follows_the_run(tmp_path: Path, cfg: GippaConfig) -> None:
    from helpers.viz_landscape import plot_spec

    genes = [vl.gene("D", "chr1", 11_001, 12_000, log2=20.0, abundance_type="maxlfq",
                     abundance_scale="relative"),
             vl.gene("U", "chr1", 12_101, 13_000, status="not_detected",
                     abundance_type="maxlfq", abundance_scale="relative")]  # fmt: skip
    result = plot_spec("V03").fn(RunView.from_outdir(_region_out(tmp_path, genes)), cfg)
    (inst,) = [i for i in result if i.suffix == "prophage_1"]
    assert h.colorbar_labels(inst.figure) == ["log2 MaxLFQ (relative)"]
