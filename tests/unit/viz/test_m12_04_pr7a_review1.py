"""M12-04 PR 7a (#563), review round 1: V01 legend entries for binned ticks, the end-label
rounding rule, transcription-start and midpoint binning, no empty region band, TU starts and the
origin gap (docs/specs/M12-04-outputs-revision.md § PR 7a; lead decisions of the review).

Builders and module fixtures come from ``test_m12_04_pr7a_v01.py`` (A-, C-, full and K-like
runs) and ``test_v01_circos*.py``. Interface assumed (test-author, review round 1):

- **Legend.** Promoter and terminator lines stay in the figure legend when their track is binned;
  their label then ends with "per bin" ("promoters per bin" style); unbinned they read
  "promoter" / "terminator". TU starts are black and have no legend line ("TUs" is in the ring
  key).
- **End label.** On a single circular replicon the last tick label is dropped only when it equals
  the replicon length rounded to the tick unit (4,020,000 bp -> "4 Mb" dropped; 4,641,652 bp ->
  "4 Mb" kept, since the length rounds to 5 Mb).
- **Binning positions.** A binned promoter / terminator counts in the bin of its transcription
  start (BED ``end - 1`` on ``-``, modulo the length); a binned novel ORF counts in the bin of its
  midpoint ``floor((lo - 1 + hi) / 2) mod length``, with ``lo``/``hi`` the 1-based low/high ends
  unwrapped with the strand (``hi + length`` when the ORF crosses the origin).
- **Region band.** Without ``region`` rows no artist of the polar axes lies on the region radii
  (``circos.RADII["region"]``), so nothing unnumbered sits outside ring 1.
- **Origin gap (lead decision 6).** Pinned as a sector gap: a single circular replicon spans
  less than a full turn (gap between 0.5 and 20 degrees), and every ring number's anchor lies in
  that gap, just before the origin.
"""

from __future__ import annotations

import math
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from matplotlib.collections import Collection
from matplotlib.figure import Figure
from matplotlib.patches import Patch
from unit.viz.test_m12_04_pr7a_v01 import (  # noqa: F401  (module fixtures)
    A_REP,
    K_BIN,
    K_LEN,
    K_REP,
    a_inst,
    bed1,
    c_inst,
    drawn_rings,
    full_inst,
    k_inst,
    polar,
    ring_numbers,
    track,
)
from unit.viz.test_v01_circos import (
    FULL_NOVEL,
    FULL_VARIANTS,
    make_run,
    novel,
    render,
    texts,
    tsv,
)
from unit.viz.test_v01_circos_review import calibrate

from gippa.config import GippaConfig, resolve_config
from gippa.viz import circos, style
from gippa.viz.base import Instance

BLACK = "#000000"


@pytest.fixture(autouse=True)
def _style() -> None:
    style.apply()


@pytest.fixture
def cfg() -> GippaConfig:
    return resolve_config()


def legend_labels(fig: Figure) -> list[str]:
    return [t.get_text() for leg in fig.legends for t in leg.get_texts() if t.get_visible()]


def entry(labels: list[str], word: str) -> list[str]:
    return [lab for lab in labels if word in lab.lower()]


# --- 1. binned promoters and terminators keep their legend entries -------------------------


def _per_bin(labels: list[str], word: str) -> None:
    found = entry(labels, word)
    assert found, f"no legend entry for {word}: {labels}"
    assert all(lab.endswith("per bin") for lab in found), found


def test_k_binned_promoters_and_terminators_are_in_the_legend(
    k_inst: tuple[Instance, pd.DataFrame],  # noqa: F811
) -> None:
    labels = legend_labels(k_inst[0].figure)
    _per_bin(labels, "promoter")
    _per_bin(labels, "terminator")


@pytest.fixture(scope="module")
def ticks_binned(tmp_path_factory: pytest.TempPathFactory) -> tuple[Instance, pd.DataFrame]:
    """Only promoters and terminators binned (600 each); 2 ALT sites and 2 novel ORFs not."""
    style.apply()
    tmp = tmp_path_factory.mktemp("ticks_binned")
    out = make_run(
        tmp,
        replicons=A_REP,
        promoters=[bed1("chr1", 101 + 300 * i, i, "P") for i in range(600)],
        terminators=[bed1("chr1", 151 + 300 * i, i, "T") for i in range(600)],
        variants=FULL_VARIANTS,
        novel_orfs=FULL_NOVEL,
    )
    inst = render(out, resolve_config())
    return inst, tsv(inst, tmp)


def test_only_ticks_binned_keeps_every_legend_entry(
    ticks_binned: tuple[Instance, pd.DataFrame],
) -> None:
    inst, data = ticks_binned
    for name in ("promoter", "terminator"):
        assert set(track(data, name)["binned"]) == {"true"}, name
    for name in ("alt_site", "novel_orf"):
        assert set(track(data, name)["binned"]) == {"false"}, name
    labels = legend_labels(inst.figure)
    _per_bin(labels, "promoter")
    _per_bin(labels, "terminator")
    assert "ALT site (high)" in labels and "novel ORF" in labels, labels


def test_unbinned_promoter_and_terminator_entries_unchanged(
    full_inst: tuple[Instance, pd.DataFrame],  # noqa: F811
) -> None:
    labels = legend_labels(full_inst[0].figure)
    assert entry(labels, "promoter") == ["promoter"], labels
    assert entry(labels, "terminator") == ["terminator"], labels


# --- 2. the end label: dropped only when it equals the rounded length ----------------------


@pytest.mark.parametrize(
    ("length", "label", "kept"),
    [
        (4_641_652, "4 Mb", True),  # rounds to 5 Mb: "4 Mb" is a plain tick
        (4_020_000, "4 Mb", False),  # rounds to 4 Mb: the end label sits near the 0 tick
        (200_000, "200 kb", False),  # AC 7a.4
        (5_000_000, "5 Mb", False),
    ],
)
def test_end_label_dropped_only_when_it_equals_the_rounded_length(
    tmp_path: Path, cfg: GippaConfig, length: int, label: str, kept: bool
) -> None:
    inst = render(make_run(tmp_path, replicons=[("chr1", length, "circular")]), cfg)
    shown = texts(inst.figure)
    assert (label in shown) is kept, (label, sorted(t for t in shown if t[:1].isdigit()))


# --- 3. binning positions: transcription start and ORF midpoint ----------------------------


def _straddling_ticks() -> list[tuple[str, int, int, str, int, str]]:
    """600 '-' BED records of 40 bp: ``start`` in bin k, ``end - 1`` (the TSS) in bin k + 1."""
    out = []
    for i in range(600):
        k = (7 * i) % (K_LEN // K_BIN - 1)
        s = K_BIN * k + 9_990 - i % 5
        out.append(("chr1", s, s + 40, f"M{i}", 1, "-"))
    return out


def _orfs() -> list[dict[str, object]]:
    """520 ORFs inside one bin, plus: a '+' ORF and a '-' ORF whose midpoint lies in the bin
    after their low end, and a '+' ORF across the origin (wrapped midpoint in bin 0)."""
    rows = []
    for i in range(520):
        lo = K_BIN * ((13 * i) % 500) + 1_001
        rows.append(novel(f"n{i}", "chr1", "+", lo, lo + 299, True))
    lo = K_BIN * 41 + 9_801
    rows.append(novel("cross_plus", "chr1", "+", lo, lo + 599, True))
    lo = K_BIN * 77 + 9_801
    rows.append(novel("cross_minus", "chr1", "-", lo + 599, lo, True))  # start > end on '-'
    rows.append(novel("origin", "chr1", "+", K_LEN - 99, 700, True))  # 4,999,901 .. 700
    return rows


def _orf_mid(row: dict[str, Any]) -> int:
    start, end = int(row["start_longest"]), int(row["end"])  # type: ignore[call-overload]
    lo, hi = (start, end) if row["strand"] == "+" else (end, start)
    hi += K_LEN if hi < lo else 0
    return math.floor((lo - 1 + hi) / 2) % K_LEN


@pytest.fixture(scope="module")
def positions(tmp_path_factory: pytest.TempPathFactory) -> pd.DataFrame:
    style.apply()
    tmp = tmp_path_factory.mktemp("positions")
    ticks = _straddling_ticks()
    out = make_run(
        tmp,
        replicons=K_REP,
        promoters=ticks,
        terminators=[(c, s, e, "T" + n[1:], v, st) for c, s, e, n, v, st in ticks],
        novel_orfs=_orfs(),
    )
    return tsv(render(out, resolve_config()), tmp)


def _bins(d: pd.DataFrame) -> dict[int, float]:
    """bin start (1-based) -> count, zero bins left out."""
    assert set(d["binned"]) == {"true"}, set(d["binned"])
    got = {int(s): float(v) for s, v in zip(d["start"], d["value"], strict=True)}
    return {s: v for s, v in got.items() if v}


def _expect(pos0: list[int]) -> dict[int, float]:
    return {k * K_BIN + 1: float(n) for k, n in Counter(p // K_BIN for p in pos0).items()}


@pytest.mark.parametrize("name", ["promoter", "terminator"])
def test_binned_minus_ticks_count_at_their_transcription_start(
    positions: pd.DataFrame, name: str
) -> None:
    tss = [(e - 1) % K_LEN for _, _, e, *_ in _straddling_ticks()]
    starts = [s for _, s, *_ in _straddling_ticks()]
    assert all(t // K_BIN != s // K_BIN for t, s in zip(tss, starts, strict=True))  # the setup
    assert _bins(track(positions, name)) == _expect(tss)


def test_binned_novel_orfs_count_at_their_midpoint(positions: pd.DataFrame) -> None:
    rows = _orfs()
    got = _bins(track(positions, "novel_orf"))
    assert got == _expect([_orf_mid(r) for r in rows])
    # the three edge cases, by name: each lands in the bin of its (wrapped) midpoint
    for row, k in ((rows[-3], 42), (rows[-2], 78), (rows[-1], 0)):
        assert _orf_mid(row) // K_BIN == k, row["orf_id"]
        assert got.get(k * K_BIN + 1, 0) >= 1, row["orf_id"]


# --- 4. no region rows, no region band ------------------------------------------------------


def _band_artists(fig: Figure, ring: str) -> list[str]:
    """Artists of the polar axes whose every vertex lies on the radii of ``ring``."""
    r0, r1 = circos.RADII[ring]
    found = []
    for a in polar(fig).get_children():
        if not getattr(a, "get_visible", lambda: False)():
            continue
        if isinstance(a, Collection):
            paths = [np.asarray(p.vertices, dtype=float) for p in a.get_paths()]
        elif isinstance(a, Patch):
            path = a.get_patch_transform().transform_path(a.get_path())
            paths = [np.asarray(path.vertices, dtype=float)]
        else:
            continue
        for v in paths:
            if len(v) and v[:, 1].min() >= r0 - 1e-6 and v[:, 1].max() <= r1 + 1e-6:
                found.append(f"{type(a).__name__} gid={a.get_gid()}")
    return found


@pytest.mark.parametrize("fixture", ["a_inst", "c_inst"])
def test_no_region_rows_draw_no_region_band(fixture: str, request: pytest.FixtureRequest) -> None:
    inst, data = request.getfixturevalue(fixture)
    assert track(data, "region").empty
    assert drawn_rings(data)[0] == "gene_plus"  # ring 1
    assert _band_artists(inst.figure, "region") == []


def test_region_rows_draw_on_the_region_band(
    full_inst: tuple[Instance, pd.DataFrame],  # noqa: F811
) -> None:
    assert _band_artists(full_inst[0].figure, "region")


# --- 5. TU starts: black, no legend line ----------------------------------------------------


def test_tu_starts_are_black_and_have_no_legend_line(
    full_inst: tuple[Instance, pd.DataFrame],  # noqa: F811
) -> None:
    from helpers.viz_tracks import elements

    inst, _ = full_inst
    els = elements(inst.figure, "tu_start")
    assert els and set().union(*(e.colours() for e in els)) == {BLACK}
    assert "TU start" not in legend_labels(inst.figure)


# --- 6. a gap at the origin of a single circular replicon -----------------------------------


def test_single_circular_replicon_leaves_an_origin_gap_holding_the_ring_numbers(
    a_inst: tuple[Instance, pd.DataFrame],  # noqa: F811
    tmp_path: Path,
    cfg: GippaConfig,
) -> None:
    ax = calibrate(tmp_path, cfg, A_REP, "chr1")
    span = ax.k * A_REP[0][1]
    gap = 2 * math.pi - span
    assert math.radians(0.5) <= gap <= math.radians(20), math.degrees(gap)
    inst, data = a_inst
    n = len(drawn_rings(data))
    for k, found in ring_numbers(inst.figure, n).items():
        ((theta, _),) = found
        before = (ax.theta(0) - theta) % (2 * math.pi)  # angle back from the origin
        assert 0 < before < gap, (k, math.degrees(before), math.degrees(gap))
