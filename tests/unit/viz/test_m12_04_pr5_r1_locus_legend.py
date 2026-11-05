"""M12-04 PR 5, review round 1 (H1, M2, L3): every gene face drawn in a V03/V13 locus view has
a legend swatch of the same style, and the data TSV's ``fill`` matches the arrow.

Lead decisions (binding, round 1):

- ``single_peptide_unestimated`` genes (not detected) are white with a ``#E69F00`` edge and no
  hatch, with and without abundance; the "single peptide, unestimated" swatch is the same.
- ``#0072B2`` (``locus.DETECTED_PLAIN``) fills detected genes only: every detected gene in a run
  without abundance; in an abundance run, the detected genes whose ``abundance_log2`` is NA.
  These get a legend swatch labelled exactly ``"detected, no abundance value"``.
- Imputed genes keep the grey ``#BBBBBB`` face and ``xxx`` hatch, and their "imputed" swatch.
- ``fill`` is ``matplotlib.colors.to_hex`` of the face (lower case ``#rrggbb``) or ``"white"``.

Interface assumed (test-author, round 1):

- A gene arrow's ``(facecolor, edgecolor, hatch)`` (hex, upper-cased here; hatch ``None`` and
  ``""`` the same) equals the ``(facecolor, edgecolor, hatch)`` of at least one ``Patch`` handle
  of the figure's legends, except a measured gene in an abundance run, which is a cividis colour
  of the colour bar. Swatches may be shown only when such a gene is drawn.
- An abundance run is one where ``style.run_abundance_label`` of ``genes.tsv`` is not ``None``
  (some ``abundance_type`` other than ``none``); a measured gene is detected, not imputed, with a
  finite ``abundance_log2``.
"""

from __future__ import annotations

import re
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
UNDETECTED_STYLE = (WHITE, GREY, "////")
SINGLE_STYLE = (WHITE, ORANGE, None)
IMPUTED_HATCH = "xxx"
NOT_DETECTED = ("not_detected", "partial_excluded")
STATUSES = ("detected", "not_detected", "partial_excluded", "single_peptide_unestimated")
NO_VALUE_KEY = "detected, no abundance value"
SINGLE_KEY = "single peptide, unestimated"
SCALE = (18.0, 24.0)
CIVIDIS = matplotlib.colormaps["cividis"]
LUT = {h.hexc(CIVIDIS(i)) for i in range(CIVIDIS.N)}
HEX = re.compile(r"#[0-9a-f]{6}")

Style = tuple[str, str, str | None]


@pytest.fixture(autouse=True)
def _style() -> None:
    style.apply()


@pytest.fixture
def cfg() -> GippaConfig:
    return resolve_config()


def _genes_kw(abundance: bool) -> dict[str, object]:
    if abundance:
        return {"abundance_type": "ibaq", "abundance_scale": "molar"}
    return {"abundance_type": "none", "abundance_scale": None}


def _render(out: Path, cfg: GippaConfig) -> Instance:
    view = RunView.from_outdir(out)
    (lo,) = [x for x in select_loci(view, cfg) if x.locus_id == "prophage_1"]
    return render_locus(lo, view, cfg, SCALE)


def _out(tmp: Path, genes: list[dict[str, object]]) -> Path:
    return vl.build_out(
        tmp,
        genes=genes,
        regions_stats=[vl.region_stat("prophage_1", "prophage", "chr1", 10_001, 20_000, 0.01)],
    )


def _hatch(value: object) -> str | None:
    return str(value) if value else None


def _gene_style(inst: Instance, tag: str) -> tuple[Style, Any]:
    g = inst.data[(inst.data["track"] == "gene") & (inst.data["id"] == tag)]
    assert len(g) == 1, inst.data
    row = g.iloc[0]
    found = h.patches_of(h.axes(inst.figure, "genes"), float(row["start"]), float(row["end"]))
    assert len(found) == 1, (tag, found)
    p = found[0]
    return (p.face, p.edge, _hatch(p.hatch)), row


def _swatches(inst: Instance) -> dict[str, Style]:
    found = {}
    for label, handle in h.legend_entries(inst.figure):
        if isinstance(handle, matplotlib.patches.Patch):
            face, edge, hatch = h.patch_style(handle)
            found[label] = (face, edge, _hatch(hatch))
    return found


def _assert_fill_format(fill: object, face: str, tag: str) -> None:
    """L3: lower-case ``#rrggbb`` exactly as ``to_hex`` writes it, or the literal ``white``."""
    assert isinstance(fill, str), (tag, fill)
    if fill != "white":
        assert HEX.fullmatch(fill), (tag, fill)
        assert mcolors.to_hex(fill) == fill, (tag, fill)
    assert h.hexc(fill) == face, (tag, fill, face)


# --- H1: the single-peptide face ------------------------------------------------------------------


@pytest.mark.parametrize("abundance", [False, True], ids=["no_abundance", "ibaq"])
def test_single_peptide_gene_is_white_orange_edged_unhatched_and_has_its_swatch(
    tmp_path: Path, cfg: GippaConfig, abundance: bool
) -> None:
    """The reviewer's repro: S was #0072B2 without abundance and #BBBBBB with iBAQ."""
    kw = _genes_kw(abundance)
    genes = [vl.gene("D", "chr1", 11_001, 12_000, log2=20.0 if abundance else None, **kw),
             vl.gene("S", "chr1", 12_101, 13_000, status="single_peptide_unestimated",
                     **kw)]  # fmt: skip
    inst = _render(_out(tmp_path, genes), cfg)
    got, row = _gene_style(inst, "S")
    assert got == SINGLE_STYLE, got
    _assert_fill_format(row["fill"], WHITE, "S")
    keys = _swatches(inst)
    assert keys.get(SINGLE_KEY) == SINGLE_STYLE, keys


# --- M2: #0072B2 only for detected genes, with a swatch -------------------------------------------


def test_abundance_run_detected_gene_without_a_value_is_blue_with_its_own_swatch(
    tmp_path: Path, cfg: GippaConfig
) -> None:
    kw = _genes_kw(True)
    genes = [vl.gene("M", "chr1", 11_001, 12_000, log2=20.0, **kw),
             vl.gene("N", "chr1", 12_101, 13_000, log2=None, **kw),
             vl.gene("U", "chr1", 13_101, 14_000, status="not_detected", **kw)]  # fmt: skip
    inst = _render(_out(tmp_path, genes), cfg)
    (face, edge, hatch), row = _gene_style(inst, "N")
    assert face == BLUE and hatch is None, (face, edge, hatch)
    _assert_fill_format(row["fill"], BLUE, "N")
    keys = _swatches(inst)
    assert keys.get(NO_VALUE_KEY) == (face, edge, hatch), keys
    m_face = _gene_style(inst, "M")[0][0]
    assert m_face in LUT and m_face != BLUE


def test_no_abundance_run_blue_is_only_on_detected_genes(tmp_path: Path, cfg: GippaConfig) -> None:
    kw = _genes_kw(False)
    genes = [vl.gene("D", "chr1", 11_001, 12_000, **kw),
             vl.gene("S", "chr1", 12_101, 13_000, status="single_peptide_unestimated", **kw),
             vl.gene("U", "chr1", 13_101, 14_000, status="not_detected", **kw),
             vl.gene("P", "chr1", 14_101, 15_000, status="partial_excluded", **kw)]  # fmt: skip
    inst = _render(_out(tmp_path, genes), cfg)
    assert _gene_style(inst, "D")[0][0] == BLUE
    for tag in ("S", "U", "P"):
        assert _gene_style(inst, tag)[0][0] != BLUE, tag


@pytest.mark.parametrize("abundance", [False, True], ids=["no_abundance", "ibaq"])
def test_imputed_gene_keeps_grey_cross_hatch_and_its_swatch(
    tmp_path: Path, cfg: GippaConfig, abundance: bool
) -> None:
    """The reviewer's repro: without abundance the imputed arrow had no swatch."""
    kw = _genes_kw(abundance)
    genes = [vl.gene("D", "chr1", 11_001, 12_000, log2=20.0 if abundance else None, **kw),
             vl.gene("I", "chr1", 13_101, 14_000, status="not_detected", imputed=True,
                     log2=15.0 if abundance else None, **kw)]  # fmt: skip
    inst = _render(_out(tmp_path, genes), cfg)
    got, row = _gene_style(inst, "I")
    assert got[0] == GREY and got[2] == IMPUTED_HATCH, got
    _assert_fill_format(row["fill"], GREY, "I")
    assert _swatches(inst).get("imputed") == got, _swatches(inst)


# --- L3: the fill column in the written TSV -------------------------------------------------------


@pytest.mark.parametrize("abundance", [False, True], ids=["no_abundance", "ibaq"])
def test_the_written_fill_column_is_lower_case_hex_or_white(
    tmp_path: Path, cfg: GippaConfig, abundance: bool
) -> None:
    kw = _genes_kw(abundance)
    genes = [vl.gene("D", "chr1", 11_001, 12_000, log2=20.0 if abundance else None, **kw),
             vl.gene("N", "chr1", 12_101, 13_000, **kw),
             vl.gene("S", "chr1", 13_101, 14_000, status="single_peptide_unestimated", **kw),
             vl.gene("U", "chr1", 14_101, 15_000, status="not_detected", **kw),
             vl.gene("I", "chr1", 15_101, 16_000, status="not_detected", imputed=True,
                     log2=15.0 if abundance else None, **kw)]  # fmt: skip
    inst = _render(_out(tmp_path / "run", genes), cfg)
    write_instance(inst, "V03", "locus", tmp_path / "plots", ["svg"], 100)
    header, rows = read_tsv(tmp_path / "plots" / "data" / "V03_locus_prophage_1.tsv")
    tsv = pd.DataFrame(rows, columns=header)
    written = tsv[tsv["track"] == "gene"].set_index("id")["fill"].to_dict()
    for tag in ("D", "N", "S", "U", "I"):
        face = _gene_style(inst, tag)[0][0]
        _assert_fill_format(written[tag], face, tag)
    if not abundance:
        assert written["D"] == "#0072b2", written


# --- property: every face has a swatch ------------------------------------------------------------

Spec = tuple[str, bool, bool, float | None, str]


@st.composite
def _mix(draw: Any) -> tuple[bool, list[Spec]]:
    abundance = draw(st.booleans())
    n = draw(st.integers(1, 7))
    genes = []
    for _ in range(n):
        log2 = draw(st.one_of(st.none(), st.floats(15.0, 27.0))) if abundance else None
        genes.append((draw(st.sampled_from(STATUSES)), draw(st.booleans()),
                      draw(st.booleans()), log2, draw(st.sampled_from(["+", "-"]))))  # fmt: skip
    return abundance, genes


@settings(max_examples=30, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(mix=_mix())
def test_any_mix_every_gene_face_is_a_legend_swatch(mix: tuple[bool, list[Spec]]) -> None:
    """For any mix of status, detected flag, imputed flag, abundance and NA log2: every arrow's
    (face, edge, hatch) is a legend swatch (or a cividis colour-bar colour for a measured gene in
    an abundance run), the per-rule styles hold, and ``fill`` matches the face."""
    abundance, spec = mix
    cfg = resolve_config()
    style.apply()
    kw = _genes_kw(abundance)
    genes = []
    for i, (status, detected, imputed, log2, strand) in enumerate(spec):
        start = 11_001 + 1_000 * i
        row = vl.gene(f"G{i}", "chr1", start, start + 799, strand=strand, status=status,
                      log2=log2, imputed=imputed, **kw)  # fmt: skip
        row["detected"] = detected
        genes.append(row)
    with tempfile.TemporaryDirectory() as tmp:
        inst = _render(_out(Path(tmp), genes), cfg)
        swatches = set(_swatches(inst).values())
        for i, (status, detected, imputed, log2, _strand) in enumerate(spec):
            tag = f"G{i}"
            got, row = _gene_style(inst, tag)
            _assert_fill_format(row["fill"], got[0], tag)
            measured = abundance and detected and not imputed and log2 is not None
            if got[0] in LUT and got[2] is None and abundance:  # the colour bar's colours
                has_value = not imputed and log2 is not None
                assert has_value and (detected or status == "detected"), (tag, spec[i], got)
            else:
                assert got in swatches, (tag, spec[i], got, swatches)
            if got[0] == BLUE:
                assert status == "detected" or detected, (tag, spec[i], got)
            if imputed and status != "single_peptide_unestimated":
                assert got[0] == GREY and got[2] == IMPUTED_HATCH, (tag, spec[i], got)
            elif imputed:
                continue
            elif status in NOT_DETECTED:
                assert got == UNDETECTED_STYLE, (tag, spec[i], got)
            elif status == "single_peptide_unestimated":
                assert got == SINGLE_STYLE, (tag, spec[i], got)
            elif detected and measured:
                want = h.hexc(CIVIDIS((min(max(float(log2), 18.0), 24.0) - 18.0) / 6.0))  # type: ignore[arg-type]
                assert got[0] == want and got[2] is None, (tag, spec[i], got)
            elif detected:
                assert got[0] == BLUE and got[2] is None, (tag, spec[i], got)
        inst.figure.clear()
