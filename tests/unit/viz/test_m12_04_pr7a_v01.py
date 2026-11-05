# ruff: noqa: RUF001, RUF002, RUF003  (figure text is Unicode: −)
"""M12-04 PR 7a, P-09: V01 ring key, abundance colour bar, marker binning, terminator colour and
the end tick (docs/specs/M12-04-outputs-revision.md § PR 7a, AC 7a.1-7a.4).

Tested on the figure artists and the data TSV as written (``plots/data/V01_circos.tsv``, 1-based
closed); no pixels, no baselines. Builders come from ``tests/unit/viz/test_v01_circos.py`` and
``test_v01_circos_review.py``. Fixtures, built by hand (R11):

- "A-like": one circular 200 kb replicon with genes on both strands (iBAQ, molar), peptides, TUs,
  promoters, terminators and GC skew, and nothing else: the spec's example key of 7 rings;
- "C-like": the same without abundance (``abundance_type`` ``none``);
- "full": the M9 two-replicon fixture with GC skew support (every ring V01 has);
- "K-like": one circular 5 Mb replicon with 781 novel ORFs (K's count), 600 ALT sites, 2,000
  promoters and 1,500 terminators, every marker inside one 10 kb bin of the 500-bin grid.

Interface assumed (test-author, PR 7a; the spec leaves these open):

- **Rings.** A ring is a group of data-TSV tracks; outer to inner: ``region``, ``gene_plus``,
  ``gene_minus``, ``abundance``, ``peptide_density``, ``tu_start``, ``promoter`` + ``terminator``
  (one ring, "promoters and terminators"), ``island``, ``gc_skew``, ``alt_site``, ``novel_orf``.
  ``sector``, ``oric`` and ``ter`` are not rings. A ring is drawn when one of its tracks has a row
  in the data TSV, and only then. With the spec's example set the order gives exactly the
  example's numbers (genes (+) = 1 … GC skew = 7).
- **Key.** One figure-level artist (a ``Legend`` or a ``Text``) with ``gid == "ring_key"``; its
  visible text lines that read ``"<n> <label>"`` are the entries, numbered 1..N in ring order.
  The spec's labels: "genes (+)", "genes (−)" (U+2212), "abundance", "peptide density", "TUs",
  "promoters and terminators", "GC skew". No entry holds an underscore.
- **Ring numbers.** The bare text ``"<n>"`` is drawn once per ring as a ``Text`` of the polar
  axes (``ax.texts``, polar data coordinates ``(theta, r)``), next to the 0 tick of the **first**
  sector (within 0.2 rad of its origin), at the ring's radius (so ``r`` falls as ``n`` rises).
- **Colour bar.** The abundance ring has a matplotlib colour bar (an axes with ``_colorbar``)
  whose label is ``style.abundance_label(type, scale)`` of the run; without abundance there is
  no colour bar and no ``abundance`` row in the data TSV.
- **Binning.** A marker track (``promoter``, ``terminator``, ``alt_site``, ``novel_orf``) with
  more than 500 rows in total is written as bins: ``binned = true``, ``category = "bin"`` (the
  AC's ``kind = bin``) and ``value`` = the count of markers whose position lies in
  ``[start, end]``. Bins with count 0 may be omitted. Each replicon has at most 500 bins of
  near-equal width (``base.binned(…, max_elements=500)``); on one replicon the edges are
  ``base.bin_edges(length, 500)``. The position of a promoter / terminator is its transcription
  start (as in M9 decision 19); the fixtures place ALT codons and novel ORFs wholly inside one
  bin so any position rule gives the same bin. At most 500 rows: one row per marker, as in 1.0.
- **End tick.** On a single circular replicon, the 0 tick reads ``"0"`` and the end label that
  equals the length (``"200 kb"`` on 200 kb) is not drawn; a linear replicon keeps its end label.
"""

from __future__ import annotations

import math
import re
from collections import Counter
from pathlib import Path
from typing import Any

import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import pytest
from helpers.viz_out_tree import gene_row
from helpers.viz_tracks import elements
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from matplotlib.text import Text
from unit.viz.test_v01_circos import (
    TWO,
    full_run,
    make_run,
    novel,
    pep,
    render,
    texts,
    tsv,
    variant,
)
from unit.viz.test_v01_circos_review import _support, calibrate

from gippa.config import GippaConfig, resolve_config
from gippa.viz import base, style
from gippa.viz.base import Instance

#: data-TSV track -> ring (outer to inner in ``RINGS``)
RING_OF = {
    "region": "region",
    "gene_plus": "gene_plus",
    "gene_minus": "gene_minus",
    "abundance": "abundance",
    "peptide_density": "peptide_density",
    "tu_start": "tu_start",
    "promoter": "promoter_terminator",
    "terminator": "promoter_terminator",
    "island": "island",
    "gc_skew": "gc_skew",
    "alt_site": "alt_site",
    "novel_orf": "novel_orf",
}
RINGS = [
    "region", "gene_plus", "gene_minus", "abundance", "peptide_density", "tu_start",
    "promoter_terminator", "island", "gc_skew", "alt_site", "novel_orf",
]  # fmt: skip
EXAMPLE_KEY = [
    (1, "genes (+)"),
    (2, "genes (−)"),
    (3, "abundance"),
    (4, "peptide density"),
    (5, "TUs"),
    (6, "promoters and terminators"),
    (7, "GC skew"),
]
MARKERS = ("promoter", "terminator", "alt_site", "novel_orf")
CAP = 500
A_LEN = 200_000
A_REP = [("chr1", A_LEN, "circular")]
K_LEN = 5_000_000
K_REP = [("chr1", K_LEN, "circular")]
K_BIN = K_LEN // CAP  # 10 kb: the 500-bin grid of 5 Mb is exact
UNIT = re.compile(r"^\d+(?:\.\d+)? ?(?:bp|kb|Mb)$")


@pytest.fixture(autouse=True)
def _style() -> None:
    style.apply()


@pytest.fixture
def cfg() -> GippaConfig:
    return resolve_config()


# --- builders -------------------------------------------------------------------------------


def g(
    tag: str,
    start: int,
    strand: str,
    log2: float | None,
    *,
    contig: str = "chr1",
    abundance: bool = True,
) -> dict[str, object]:
    """One ``genes.tsv`` row; iBAQ (molar) when ``abundance``, else ``none`` with NA log2."""
    return gene_row(
        locus_tag=tag,
        contig=contig,
        start=start,
        end=start + 899,
        strand=strand,
        feature_type="cds",
        detection_status="detected" if log2 is not None or not abundance else "not_detected",
        detected=log2 is not None or not abundance,
        partial=False,
        crosses_origin=False,
        abundance_log2=log2 if abundance else None,
        abundance_imputed=False,
        abundance_type="ibaq" if abundance else "none",
        abundance_scale="molar" if abundance else None,
    )


def a_like(tmp: Path, *, abundance: bool = True, name: str = "A") -> Path:
    """The spec's example rings: genes (+/−), abundance, peptides, TUs, promoters,
    terminators and GC skew on one circular 200 kb replicon; no regions, islands or markers."""
    genes = [
        g("g1", 1_001, "+", 20.0, abundance=abundance),
        g("g2", 5_001, "-", 22.5, abundance=abundance),
        g("g3", 60_001, "+", 18.0, abundance=abundance),
        g("g4", 120_001, "-", None, abundance=abundance),
    ]
    out = make_run(
        tmp,
        replicons=A_REP,
        genes=genes,
        peptides=[pep("chr1", s) for s in (1_101, 5_201, 60_301)],
        tus=[
            ["TU_1", "chr1", 1_001, 1_900, "+", ["g1"], 801, 1_000],
            ["TU_2", "chr1", 5_001, 5_900, "-", ["g2"], 5_901, 6_100],
        ],
        promoters=[("chr1", 900, 950, "promoter_1", 500, "+")],
        terminators=[("chr1", 1_950, 1_980, "T_1", 800, "+")],
        name=name,
    )
    _support(out, A_REP)
    return out


def bed1(contig: str, pos0: int, i: int, prefix: str) -> tuple[str, int, int, str, int, str]:
    """A 1-bp BED record at 0-based ``pos0``: its transcription start on either strand."""
    return (contig, pos0, pos0 + 1, f"{prefix}{i}", 1, "+" if i % 2 == 0 else "-")


def k_positions() -> dict[str, list[int]]:
    """0-based marker positions of the K-like fixture; ALT codons (3 bp) and novel ORFs (300 bp)
    lie wholly inside one 10 kb bin."""
    return {
        "promoter": [1_234 + 2_497 * i for i in range(2_000)],
        "terminator": [777 + 3_331 * i for i in range(1_500)],
        "alt_site": [K_BIN * ((7 * i) % CAP) + 1_000 + (37 * i) % 8_000 for i in range(600)],
        "novel_orf": [K_BIN * ((13 * i) % CAP) + 1_000 + (53 * i) % 8_000 for i in range(781)],
    }


def k_like(tmp: Path) -> Path:
    pos = k_positions()
    genes = [g(f"g{i}", 1 + 400_000 * i, "+-"[i % 2], 15.0 + i) for i in range(12)]
    orfs = [
        novel(f"n{i}", "chr1", "+", p + 1, p + 300, True) for i, p in enumerate(pos["novel_orf"])
    ]
    return make_run(
        tmp,
        replicons=K_REP,
        genes=genes,
        promoters=[bed1("chr1", p, i, "P") for i, p in enumerate(pos["promoter"])],
        terminators=[bed1("chr1", p, i, "T") for i, p in enumerate(pos["terminator"])],
        variants=[
            variant(f"v{i}", "chr1", p + 1, "ALT", "high") for i, p in enumerate(pos["alt_site"])
        ],
        novel_orfs=orfs,
        name="K",
    )


# --- readers --------------------------------------------------------------------------------


def _draw(fig: Figure) -> None:
    canvas = fig.canvas if isinstance(fig.canvas, FigureCanvasAgg) else FigureCanvasAgg(fig)
    canvas.draw()


def key_entries(fig: Figure) -> list[tuple[int, str]]:
    """The numbered entries ``(n, label)`` of the ``ring_key`` artist(s), in drawing order."""
    _draw(fig)
    keys = fig.findobj(lambda a: a.get_gid() == "ring_key")
    assert keys, "no artist with gid 'ring_key' (the V01 ring key)"
    seen: set[int] = set()
    out = []
    for k in keys:
        for t in k.findobj(Text):
            if id(t) in seen or not t.get_visible():
                continue
            seen.add(id(t))
            for line in t.get_text().splitlines():
                m = re.fullmatch(r"\s*(\d+)\s+(\S.*?)\s*", line)
                if m:
                    out.append((int(m.group(1)), m.group(2)))
    return out


def polar(fig: Figure) -> Any:
    (ax,) = [a for a in fig.axes if a.name == "polar"]
    return ax


def ring_numbers(fig: Figure, n: int) -> dict[int, list[tuple[float, float]]]:
    """``k -> [(theta, r), …]`` of every bare ``"k"`` text of the polar axes, k in 1..n."""
    _draw(fig)
    found: dict[int, list[tuple[float, float]]] = {k: [] for k in range(1, n + 1)}
    for t in polar(fig).texts:
        s = t.get_text().strip()
        if t.get_visible() and s.isdigit() and 1 <= int(s) <= n:
            theta, r = t.get_position()
            found[int(s)].append((float(theta), float(r)))
    return found


def drawn_rings(data: pd.DataFrame) -> list[str]:
    present = {RING_OF[t] for t in data["track"] if t in RING_OF}
    return [r for r in RINGS if r in present]


def colour_bars(fig: Figure) -> list[Any]:
    return [a for a in fig.axes if getattr(a, "_colorbar", None) is not None]


def bar_label(ax: Any) -> str:
    return str(ax.get_ylabel() or ax.get_xlabel())


def angle(a: float, b: float) -> float:
    d = (a - b) % (2 * math.pi)
    return min(d, 2 * math.pi - d)


def track(data: pd.DataFrame, name: str) -> pd.DataFrame:
    return data[data["track"] == name]


# --- module fixtures ------------------------------------------------------------------------


@pytest.fixture(scope="module")
def a_inst(tmp_path_factory: pytest.TempPathFactory) -> tuple[Instance, pd.DataFrame]:
    style.apply()
    tmp = tmp_path_factory.mktemp("a_like")
    inst = render(a_like(tmp), resolve_config())
    return inst, tsv(inst, tmp)


@pytest.fixture(scope="module")
def c_inst(tmp_path_factory: pytest.TempPathFactory) -> tuple[Instance, pd.DataFrame]:
    style.apply()
    tmp = tmp_path_factory.mktemp("c_like")
    inst = render(a_like(tmp, abundance=False, name="C"), resolve_config())
    return inst, tsv(inst, tmp)


@pytest.fixture(scope="module")
def full_inst(tmp_path_factory: pytest.TempPathFactory) -> tuple[Instance, pd.DataFrame]:
    style.apply()
    tmp = tmp_path_factory.mktemp("full")
    out = full_run(tmp)
    _support(out, TWO)
    inst = render(out, resolve_config())
    return inst, tsv(inst, tmp)


@pytest.fixture(scope="module")
def k_inst(tmp_path_factory: pytest.TempPathFactory) -> tuple[Instance, pd.DataFrame]:
    style.apply()
    tmp = tmp_path_factory.mktemp("k_like")
    inst = render(k_like(tmp), resolve_config())
    return inst, tsv(inst, tmp)


FIXTURES = ("a_inst", "c_inst", "full_inst", "k_inst")


# --- AC 7a.1: the ring key ------------------------------------------------------------------


@pytest.mark.parametrize("fixture", FIXTURES)
def test_ac7a1_one_numbered_key_entry_per_drawn_ring(
    fixture: str, request: pytest.FixtureRequest
) -> None:
    inst, data = request.getfixturevalue(fixture)
    rings = drawn_rings(data)
    entries = key_entries(inst.figure)
    assert [n for n, _ in entries] == list(range(1, len(rings) + 1)), (entries, rings)


def test_example_key_reads_the_spec_example(a_inst: tuple[Instance, pd.DataFrame]) -> None:
    inst, data = a_inst
    assert drawn_rings(data) == [
        "gene_plus", "gene_minus", "abundance", "peptide_density", "tu_start",
        "promoter_terminator", "gc_skew",
    ]  # fmt: skip
    assert key_entries(inst.figure) == EXAMPLE_KEY


def test_gc_skew_line_is_labelled_7(a_inst: tuple[Instance, pd.DataFrame]) -> None:
    assert (7, "GC skew") in key_entries(a_inst[0].figure)


def test_no_abundance_drops_the_ring_and_renumbers(c_inst: tuple[Instance, pd.DataFrame]) -> None:
    inst, data = c_inst
    assert track(data, "abundance").empty
    assert elements(inst.figure, "abundance") == []
    assert key_entries(inst.figure) == [
        (1, "genes (+)"),
        (2, "genes (−)"),
        (3, "peptide density"),
        (4, "TUs"),
        (5, "promoters and terminators"),
        (6, "GC skew"),
    ]


@pytest.mark.parametrize("fixture", FIXTURES)
def test_every_ring_number_is_drawn_once_at_the_first_origin_outer_to_inner(
    fixture: str, request: pytest.FixtureRequest, tmp_path: Path, cfg: GippaConfig
) -> None:
    inst, data = request.getfixturevalue(fixture)
    n = len(drawn_rings(data))
    found = ring_numbers(inst.figure, n)
    assert {k: len(v) for k, v in found.items()} == dict.fromkeys(range(1, n + 1), 1), found
    reps = TWO if fixture == "full_inst" else (K_REP if fixture == "k_inst" else A_REP)
    origin = calibrate(tmp_path, cfg, reps, reps[0][0]).theta(0)
    for k, ((theta, _),) in found.items():
        assert angle(theta, origin) <= 0.2, (k, theta, origin)
    radii = [found[k][0][1] for k in range(1, n + 1)]
    assert radii == sorted(radii, reverse=True) and len(set(radii)) == n, radii


def test_several_contigs_share_one_key(full_inst: tuple[Instance, pd.DataFrame]) -> None:
    inst, data = full_inst
    entries = key_entries(inst.figure)
    assert len(entries) == len(set(entries)) == len(drawn_rings(data))


@pytest.mark.parametrize("fixture", FIXTURES)
def test_key_labels_are_plain_words(fixture: str, request: pytest.FixtureRequest) -> None:
    inst, _ = request.getfixturevalue(fixture)
    for _, label in key_entries(inst.figure):
        assert "_" not in label, label
        assert label not in RING_OF, label


# --- colour bar -----------------------------------------------------------------------------


def test_abundance_ring_has_a_colour_bar_with_the_abundance_label(
    a_inst: tuple[Instance, pd.DataFrame],
) -> None:
    inst, _ = a_inst
    bars = colour_bars(inst.figure)
    assert len(bars) == 1, f"{len(bars)} colour bars"
    assert bar_label(bars[0]) == style.abundance_label("ibaq", "molar") == "log2 iBAQ (molar)"


def test_colour_bar_label_follows_the_measure(tmp_path: Path, cfg: GippaConfig) -> None:
    genes = [
        {**g("g1", 1_001, "+", 20.0), "abundance_type": "maxlfq", "abundance_scale": "relative"},
        {**g("g2", 9_001, "-", 21.0), "abundance_type": "maxlfq", "abundance_scale": "relative"},
    ]
    inst = render(make_run(tmp_path, replicons=A_REP, genes=genes), cfg)
    bars = colour_bars(inst.figure)
    assert len(bars) == 1, f"{len(bars)} colour bars"
    assert bar_label(bars[0]) == "log2 MaxLFQ (relative)"


def test_no_abundance_no_colour_bar(c_inst: tuple[Instance, pd.DataFrame]) -> None:
    assert colour_bars(c_inst[0].figure) == []


# --- AC 7a.2: marker binning ----------------------------------------------------------------


@pytest.mark.parametrize("name", MARKERS)
def test_ac7a2_k_marker_tracks_are_bins_on_the_500_grid(
    k_inst: tuple[Instance, pd.DataFrame], name: str
) -> None:
    _, data = k_inst
    d = track(data, name)
    assert len(d) > 0
    assert set(d["binned"]) == {"true"}, set(d["binned"])
    assert set(d["category"]) == {"bin"}, set(d["category"])
    assert len(d) <= CAP
    edges = set(base.bin_edges(K_LEN, CAP).tolist())
    for s, e in zip(d["start"].astype(int), d["end"].astype(int), strict=True):
        assert s - 1 in edges and e in edges and e - (s - 1) == K_BIN, (s, e)


@pytest.mark.parametrize("name", MARKERS)
def test_ac7a2_k_bin_counts_are_the_markers_in_each_bin(
    k_inst: tuple[Instance, pd.DataFrame], name: str
) -> None:
    _, data = k_inst
    d = track(data, name)
    assert set(d["binned"]) == {"true"}, f"{len(d)} element rows, not bins"
    got = {
        (int(s), int(e)): float(v) for s, e, v in zip(d["start"], d["end"], d["value"], strict=True)
    }
    pos = k_positions()[name]
    want = Counter((p // K_BIN) * K_BIN + 1 for p in pos)
    assert {k: v for k, v in got.items() if v != 0} == {
        (s, s + K_BIN - 1): float(n) for s, n in want.items()
    }
    assert sum(got.values()) == len(pos)


@pytest.mark.parametrize("name", MARKERS)
def test_k_binned_marker_tracks_draw_a_density_not_elements(
    k_inst: tuple[Instance, pd.DataFrame], name: str
) -> None:
    inst, data = k_inst
    n = len(elements(inst.figure, name))
    assert 0 < n <= len(track(data, name)) <= CAP, n


def test_k_title_still_counts_elements(k_inst: tuple[Instance, pd.DataFrame]) -> None:
    assert "600 ALT sites, 781 novel ORFs" in k_inst[0].title


def test_exactly_500_markers_are_not_binned(tmp_path: Path, cfg: GippaConfig) -> None:
    proms = [bed1("chr1", 101 + 300 * i, i, "P") for i in range(CAP)]
    alts = [variant(f"v{i}", "chr1", 151 + 300 * i, "ALT", "high") for i in range(CAP)]
    data = tsv(render(make_run(tmp_path, replicons=A_REP, promoters=proms, variants=alts), cfg),
               tmp_path)  # fmt: skip
    for name, cat in (("promoter", {f"P{i}" for i in range(CAP)}), ("alt_site", {"high"})):
        d = track(data, name)
        assert len(d) == CAP, name
        assert set(d["binned"]) == {"false"}, name
        assert set(d["category"]) <= cat, name


def test_501_markers_are_binned(tmp_path: Path, cfg: GippaConfig) -> None:
    proms = [bed1("chr1", 101 + 300 * i, i, "P") for i in range(CAP + 1)]
    data = tsv(render(make_run(tmp_path, replicons=A_REP, promoters=proms), cfg), tmp_path)
    d = track(data, "promoter")
    assert set(d["binned"]) == {"true"}
    assert set(d["category"]) == {"bin"}
    assert d["value"].astype(float).sum() == CAP + 1


# --- property: every marker in exactly one bin ---------------------------------------------


def _marker_run(tmp: Path, name: str, chr_pos: list[int], pl_pos: list[int]) -> Path:
    reps = [("chr1", 100_000, "circular"), ("plasmid_1", 8_000, "circular")]
    placed = [("chr1", p) for p in chr_pos] + [("plasmid_1", p) for p in pl_pos]
    kw: dict[str, Any] = {}
    if name in ("promoter", "terminator"):
        kw[name + "s"] = [bed1(c, p, i, name[0].upper()) for i, (c, p) in enumerate(placed)]
    elif name == "alt_site":
        kw["variants"] = [
            variant(f"v{i}", c, p + 1, "ALT", "high") for i, (c, p) in enumerate(placed)
        ]
    else:
        kw["novel_orfs"] = [
            novel(f"n{i}", c, "+", p + 1, p + 3, True) for i, (c, p) in enumerate(placed)
        ]
    return make_run(tmp, replicons=reps, **kw)


@settings(
    max_examples=10, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture]
)
@given(
    name=st.sampled_from(MARKERS),
    n=st.integers(CAP - 4, CAP + 8),
    share=st.floats(0.0, 0.5),
    seed=st.integers(0, 2**32 - 1),
)
def test_property_markers_fall_in_exactly_one_bin_and_counts_add_up(
    tmp_path_factory: pytest.TempPathFactory, name: str, n: int, share: float, seed: int
) -> None:
    rng = np.random.default_rng(seed)
    n_pl = int(n * share)
    chr_pos = sorted(rng.integers(0, 100_000 - 3, n - n_pl).tolist())
    pl_pos = sorted(rng.integers(0, 8_000 - 3, n_pl).tolist())
    tmp = tmp_path_factory.mktemp("bins")
    data = tsv(render(_marker_run(tmp, name, chr_pos, pl_pos), resolve_config()), tmp)
    d = track(data, name)
    if n <= CAP:  # strictly more than 500 bins
        assert len(d) == n and set(d["binned"]) == {"false"}
        return
    assert set(d["binned"]) == {"true"} and set(d["category"]) == {"bin"}
    v = d["value"].astype(float)
    assert (v >= 0).all() and (v == v.round()).all()
    assert v.sum() == n
    lengths = {"chr1": 100_000, "plasmid_1": 8_000}
    for contig, pos in (("chr1", chr_pos), ("plasmid_1", pl_pos)):
        c = d[d["contig"] == contig]
        s, e = c["start"].astype(int).to_numpy(), c["end"].astype(int).to_numpy()
        assert len(c) <= CAP, contig
        assert (s <= e).all() and (s[1:] > e[:-1]).all(), contig  # sorted, disjoint
        if len(c):
            assert s.min() >= 1 and e.max() <= lengths[contig], contig
            w = e - s + 1
            assert w.max() - w.min() <= 1, contig  # near-equal widths
        assert c["value"].astype(float).sum() == len(pos), contig
        if name in ("promoter", "terminator"):  # 1 bp: the position is unambiguous
            x = np.asarray(pos) + 1
            inside = (x[:, None] >= s[None, :]) & (x[:, None] <= e[None, :])
            assert (inside.sum(axis=1) == 1).all(), contig  # every marker in exactly one bin
            assert inside.sum(axis=0).tolist() == c["value"].astype(float).astype(int).tolist()


# --- AC 7a.3: marker colours ----------------------------------------------------------------


@pytest.mark.parametrize("fixture", ["full_inst", "k_inst"])
def test_ac7a3_terminator_ticks_are_black(fixture: str, request: pytest.FixtureRequest) -> None:
    inst, _ = request.getfixturevalue(fixture)
    els = elements(inst.figure, "terminator")
    assert els
    assert all(e.colours() == {"#000000"} for e in els), {frozenset(e.colours()) for e in els}


@pytest.mark.parametrize("fixture", ["full_inst", "k_inst"])
def test_ac7a3_no_marker_track_uses_a_region_colour(
    fixture: str, request: pytest.FixtureRequest
) -> None:
    inst, _ = request.getfixturevalue(fixture)
    region = {mcolors.to_hex(c).upper() for c in style.REGION.values()}
    for name in MARKERS:
        els = elements(inst.figure, name)
        assert els, name
        used = set().union(*(e.colours() for e in els))
        assert not used & region, (name, sorted(used & region))


def test_terminator_key_entry_is_black(full_inst: tuple[Instance, pd.DataFrame]) -> None:
    """The figure's legend handle for terminators (when one is drawn) shows the new colour."""
    fig = full_inst[0].figure
    _draw(fig)
    for leg in fig.legends:
        for handle, text in zip(leg.legend_handles, leg.get_texts(), strict=False):
            if "terminator" in text.get_text().lower() and hasattr(handle, "get_color"):
                assert mcolors.to_hex(handle.get_color()).upper() == "#000000", text.get_text()


# --- AC 7a.4: the end tick ------------------------------------------------------------------


def test_ac7a4_round_circular_length_is_not_labelled_at_angle_0(
    a_inst: tuple[Instance, pd.DataFrame],
) -> None:
    inst, _ = a_inst
    shown = texts(inst.figure)
    assert "200 kb" not in shown, "the end label equals the length and sits on the 0 tick"
    ax = polar(inst.figure)
    at0 = [
        t.get_text()
        for t in ax.texts
        if t.get_visible()
        and t.get_text().strip()
        and angle(float(t.get_position()[0]), 0.0) < 1e-3
    ]
    lengths = [t for t in at0 if UNIT.match(t.strip())]
    assert lengths == [], at0
    assert "0" in [t.strip() for t in at0], at0


def test_linear_replicon_keeps_its_end_label(tmp_path: Path, cfg: GippaConfig) -> None:
    inst = render(make_run(tmp_path, replicons=[("chr1", 100_000, "linear")]), cfg)
    assert "100 kb" in texts(inst.figure)


def test_inner_ticks_keep_their_labels(a_inst: tuple[Instance, pd.DataFrame]) -> None:
    shown = texts(a_inst[0].figure)
    for label in ("50 kb", "100 kb", "150 kb"):
        assert label in shown, label
