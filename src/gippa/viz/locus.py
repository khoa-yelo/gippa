"""V3 ``locus_<id>``: what happens at this prophage, BGC, island or flagged operon?
(GIPPA_SPEC.md §10.2 V3, §10.1, §6.5; M9-10 and its lead decisions 1-9).

``select_loci`` builds the whole ordered pool (prophage/BGC ``o_minus_e`` rows, islands, operon
flags at ``q <= stats.alpha.operon_flags``; novel ORFs are V13's). The plot renders the first
``report.max_locus_plots`` loci, one ``Instance`` each; the index of the whole pool is an
image-free data table (``base.data_table``: ``data/V03_locus.tsv``) that the report shows as
"Locus views, in display order" (M12-04 PR 4, P-12). Each locus is drawn in a window of its
span +/- ``FLANK``: gene arrows coloured by cividis on one run-wide scale, a peptide pileup,
promoter and terminator marks, region and island shading. Circular windows stay unwrapped
(features are shifted by +/- the replicon length into the window frame, the origin is a dashed
line); linear windows are clipped.
Genes are drawn with plain matplotlib patches on our own ``Figure`` (decision 9: no pyplot).
"""

from __future__ import annotations

import math
import textwrap
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Literal

import matplotlib
import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize, to_hex, to_rgba
from matplotlib.figure import Figure
from matplotlib.patches import Patch, Polygon
from matplotlib.ticker import FuncFormatter, MaxNLocator

from gippa.viz import base, style

if TYPE_CHECKING:
    from gippa.config import GippaConfig
    from gippa.io.run_view import RunView

Kind = Literal["prophage", "bgc", "island", "operon_flag", "novel_orf"]
KIND_ORDER: tuple[str, ...] = ("prophage", "bgc", "island", "operon_flag", "novel_orf")
#: the §8.1 test families in their row order; the V3 order rotates over these (V03-order spec)
FAMILY_ORDER: tuple[str, ...] = ("prophage", "bgc", "island", "operon_discordance",
                                 "missed_operon", "novel_orf")  # fmt: skip
FLAG_FAMILIES: tuple[str, ...] = ("operon_discordance", "missed_operon")
#: the name of each family's own estimate, written in the index (``estimate``); BH q otherwise
ESTIMATE = {"island": "FWER p", "missed_operon": "decoy-FDR q"}
FLANK = 2_000  # bp each side of the span (index decision 20)
SINGLE_MAX_BP = 20_000
MAX_TICKS = 5  # labelled position ticks at most (K-12 survey: 7-digit bp labels collided)
HEIGHT_IN = 2.2
MAX_LABELS = 12
GENE_LABEL = "gene_label"  # gid of the gene name texts, thinned by _thin_labels
UNDETECTED = "#BBBBBB"  # the imputed face (no measured value)
DETECTED_PLAIN = style.DETECTION["detected"]  # a detected gene without a value to colour
NOT_DETECTED = ("not_detected", "partial_excluded")
UNESTIMATED = style.DETECTION["unestimated"]
PROMOTER = style.OKABE_ITO["bluish_green"]
TERMINATOR = style.OKABE_ITO["vermillion"]
NOVEL_ORF = style.OKABE_ITO["reddish_purple"]
ORF_Y = 1.7  # the novel ORF arrow's centre line: its own band above gene labels and marks (V13)
ORF_TOP = 2.2  # the genes axes' top with a novel ORF band
ORF_EXTRA_IN = 0.7  # extra figure height for the band and the key
ORF_KEY = "novel ORF; tick: OSTIR start; purple pileup: its peptides"
IMPUTED_HATCH = "xxx"
#: ``fill`` is the gene arrow's face colour (NA on the other tracks; M12-04 P-01)
COLUMNS = ("track", "id", "contig", "start", "end", "strand", "value", "category", "fill")
INDEX_COLUMNS = ("rank", "locus_id", "kind", "family", "contig", "start", "end", "estimate",
                 "error", "rendered")  # fmt: skip
EMPTY = "no loci to show"
_OPTIONAL = ("regions_stats", "islands", "operon_flags", "tus_tsv", "regions_tsv")
_OPTIONAL_TRACKS = ("promoters_bed", "terminators_bed", "peptide_hits")


@dataclass(frozen=True)
class Locus:
    """One candidate locus; ``start0``/``end0`` are its span, 0-based half-open (``end0`` may
    exceed the length on circular contigs). ``error`` is the locus family's own estimate
    (``None`` = NA): the BH q for prophage and BGC regions, the FWER p for islands, the
    discordance BH q or the missed-operon decoy-FDR q for operon flags. Estimates of different
    families are on different scales and are never compared with each other. ``family`` is the
    §8.1 test family (empty = the same as ``kind``); operon flags carry ``operon_discordance``
    (``pair_class = OP``) or ``missed_operon`` (``NO``)."""

    locus_id: str
    kind: Kind
    contig: str
    start0: int
    end0: int
    error: float | None
    label: str
    family: str = ""

    @property
    def test_family(self) -> str:
        """The §8.1 family that ``error`` belongs to."""
        return self.family or self.kind


@dataclass(frozen=True)
class OrfMark:
    """V13's novel ORF on top of a ``Locus(kind="novel_orf")``: its strand, the
    ``start_best_ostir`` tick in the window frame (``None`` = no tick) and its peptides."""

    strand: str
    tick: int | None = None
    peptides: frozenset[str] = frozenset()


def _kb(bp: float) -> str:
    out = f"{bp / 1_000:.3f}".rstrip("0").rstrip(".")
    return "0" if out == "-0" else out


def wrapped_formatter(length: int | None) -> FuncFormatter:
    """Tick labels in kb; on an unwrapped circular window (``length``) ``x mod length``, never
    negative."""
    if length is None:
        return FuncFormatter(lambda x, _pos: _kb(round(x)))
    return FuncFormatter(lambda x, _pos: _kb(round(x) % length))


def fold_il(peptide: str) -> str:
    """The house I/L-folded peptide key (``gippa.mapping``: I is read as L)."""
    return peptide.upper().replace("I", "L")


def _num(v: object) -> float | None:
    """A finite float, or ``None`` for NA / NaN."""
    if v is None or v is pd.NA:
        return None
    f = float(str(v)) if not isinstance(v, int | float | np.number) else float(v)
    return f if math.isfinite(f) else None


def _pct(x: float | None, n: float | None) -> str:
    if x is None or not n:
        return "NA"
    return f"{int(100 * x / n + 0.5)} %"


def _frame(view: RunView, key: str) -> pd.DataFrame:
    t = view.table(key)
    return pd.DataFrame() if t is None else t


def _region_loci(view: RunView) -> list[Locus]:
    rs = _frame(view, "regions_stats")
    if rs.empty:
        return []
    keep = rs[
        (rs["stat"] == "o_minus_e")
        & (rs["subset"] == "all")
        & rs["region_type"].isin(["prophage", "bgc"])
    ]
    out = []
    for r in keep.itertuples(index=False):
        start, end = _num(r.start), _num(r.end)
        if start is None or end is None:  # no window, so not a V3 locus
            continue
        q = _num(r.q_value)
        obs, exp, n = _num(r.observed), _num(r.expected), _num(r.n_eligible)
        label = (f"{r.region_id}: {_pct(obs, n)} detected vs {_pct(exp, n)} expected "
                 f"({style.q_equals(q)})")  # fmt: skip
        out.append(
            Locus(str(r.region_id), r.region_type, str(r.contig), int(start), int(end), q, label)
        )
    return out


def _island_loci(view: RunView) -> list[Locus]:
    out = []
    for r in _frame(view, "islands").itertuples(index=False):
        p, mean = _num(r.p_value), _num(r.mean_abundance_log2)
        label = (
            f"{r.island_id} ({r.direction}): mean log2 {style.fmt_num(mean, 3)}, "
            f"FWER p {style.fmt_q(p)}"
        )
        out.append(Locus(str(r.island_id), "island", str(r.contig), int(r.start), int(r.end), p,
                         label))  # fmt: skip
    return out


def _flag_loci(view: RunView, alpha: float) -> list[Locus]:
    flags = _frame(view, "operon_flags")
    if flags.empty:
        return []
    genes = _frame(view, "genes_tsv")
    by_tag = {str(g.locus_tag): g for g in genes.itertuples(index=False)} if len(genes) else {}
    tus = _frame(view, "tus_tsv")
    tu_span = {str(t.tu_id): (int(t.start), int(t.end)) for t in tus.itertuples(index=False)}
    kept = []
    for f in flags.itertuples(index=False):
        q = _num(f.q_value)
        a, b = by_tag.get(str(f.gene_a)), by_tag.get(str(f.gene_b))
        if q is None or q > alpha or a is None or b is None:
            continue
        tu_a, tu_b = a.tu_id, b.tu_id
        if isinstance(tu_a, str) and tu_a == tu_b and tu_a in tu_span:
            span = tu_span[tu_a]
        else:
            span = _pair_span(view, a, b)
        clause = f"{f.flags} between {f.gene_a} and {f.gene_b} ({style.q_equals(q)})"
        family = "missed_operon" if str(f.pair_class) == "NO" else "operon_discordance"
        kept.append((q, str(f.pair_id), family, str(a.contig), span, clause))
    # merge within a family only, so a merged locus never takes a minimum across families
    merged: dict[tuple[str, str, int, int], list[tuple[float, str, str]]] = {}
    for q, pid, family, contig, span, clause in sorted(kept, key=lambda k: (k[0], k[1])):
        merged.setdefault((family, contig, *span), []).append((q, pid, clause))
    return [
        Locus(group[0][1], "operon_flag", contig, s, e, group[0][0],
              "; ".join(c for _, _, c in group), family)
        for (family, contig, s, e), group in merged.items()
    ]  # fmt: skip


def _pair_span(view: RunView, a: Any, b: Any) -> tuple[int, int]:
    """The span of genes ``a`` and ``b``; on a circular contig the shorter arc (``b`` unwrapped
    by +/- the length relative to ``a``), with the start shifted into ``[0, length)``."""
    length, circular = _replicon(view, str(a.contig))
    shifts = (0, -length, length) if circular and length else (0,)
    sa, ea = int(a.start), int(a.end)
    spans = [(min(sa, int(b.start) + k), max(ea, int(b.end) + k)) for k in shifts]
    s, e = min(spans, key=lambda t: t[1] - t[0])
    return (s + length, e + length) if s < 0 and length else (s, e)


def family_alpha(family: str, cfg: GippaConfig) -> float:
    """The level of ``family``: ``stats.alpha.operon_flags`` for both flag families,
    ``stats.alpha.context`` otherwise (regions, islands)."""
    alpha = cfg.stats.alpha
    return alpha.operon_flags if family in FLAG_FAMILIES else alpha.context


def tier(lo: Locus, cfg: GippaConfig) -> int:
    """0 when ``error <= alpha`` of the locus family (plain float comparison), 1 above, 2 NA."""
    if lo.error is None:
        return 2
    return 0 if lo.error <= family_alpha(lo.test_family, cfg) else 1


def _order(pool: Sequence[Locus], cfg: GippaConfig) -> list[Locus]:
    """Sort by ``(tier, family_rank, family_index, locus_id)``. ``family_rank`` is the 1-based
    position within the same family and tier by ``(error, locus_id)`` (0 in tier 2), so the
    families take turns within a tier and no error is compared with another family's."""
    groups: dict[tuple[int, str], list[Locus]] = {}
    for lo in pool:
        groups.setdefault((tier(lo, cfg), lo.test_family), []).append(lo)
    keyed: list[tuple[tuple[int, int, int, str], Locus]] = []
    for (t, family), members in groups.items():
        members.sort(key=lambda lo: (lo.error or 0.0, lo.locus_id))
        for rank, lo in enumerate(members, start=1):
            key = (t, 0 if t == 2 else rank, FAMILY_ORDER.index(family), lo.locus_id)
            keyed.append((key, lo))
    return [lo for _, lo in sorted(keyed, key=lambda kl: kl[0])]


def select_loci(view: RunView, cfg: GippaConfig) -> list[Locus]:
    """The whole V3 pool (decision 5) in display order (docs/specs/V03-order-within-kind.md).

    Loci significant at their own family's level come first, then the non-significant ones, then
    those with no estimate. Within the first two tiers the §8.1 families (prophage, BGC, island,
    operon discordance, missed operon) take turns: the best of each family, then the second
    best, and so on; within a family the order is by its own estimate (§10.2 "ordered by q"
    applies within each family). The resulting position (the index ``rank``) is a display order,
    not a ranking of evidence across families: their estimates are on different scales.
    """
    pool = _region_loci(view) + _island_loci(view)
    pool += _flag_loci(view, cfg.stats.alpha.operon_flags)
    return _order(pool, cfg)


# --- rendering ------------------------------------------------------------------------------


def _replicon(view: RunView, contig: str) -> tuple[int | None, bool]:
    pf = view.preflight() or {}
    for r in pf.get("replicons", []):
        if r.get("id") == contig:
            return int(r["length"]), r.get("topology") == "circular"
    return None, False


def _window(lo: Locus, length: int | None, circular: bool) -> tuple[int, int]:
    ws, we = lo.start0 - FLANK, lo.end0 + FLANK
    if circular:
        return ws, we
    return max(ws, 0), we if length is None else min(we, length)


def _in_window(
    frame: pd.DataFrame,
    contig_col: str,
    contig: str,
    window: tuple[int, int],
    length: int | None,
    circular: bool,
) -> list[tuple[Any, int, int]]:
    """``(row, start, end)`` of the rows of ``contig`` overlapping ``window``, shifted by
    +/- ``length`` into the unwrapped window frame on circular contigs. Rows with NA ``start``
    or ``end`` (``regions.tsv`` AMR gene sets, §7.1.7) overlap nothing and are skipped."""
    if frame.empty:
        return []
    shifts = (0, -length, length) if circular and length else (0,)
    out = []
    for r in frame[frame[contig_col] == contig].itertuples(index=False):
        start, end = _num(r.start), _num(r.end)
        if start is None or end is None:
            continue
        for k in shifts:
            s, e = int(start) + k, int(end) + k
            if s < window[1] and e > window[0]:
                out.append((r, s, e))
                break
    return out


def run_scale(genes: pd.DataFrame | None) -> tuple[float, float]:
    """The 2nd-98th percentile of measured (detected, not imputed, finite) ``abundance_log2``."""
    if genes is None or genes.empty:
        return (0.0, 1.0)
    det = genes["detected"].fillna(False).astype(bool)
    imp = genes["abundance_imputed"].fillna(False).astype(bool)
    vals = pd.to_numeric(genes["abundance_log2"], errors="coerce")[det & ~imp].to_numpy(float)
    vals = vals[np.isfinite(vals)]
    if len(vals) == 0:
        return (0.0, 1.0)
    lo, hi = (float(x) for x in np.percentile(vals, [2, 98]))
    return (lo, hi) if hi > lo else (lo - 0.5, lo + 0.5)


def _arrow(s: int, e: int, strand: str, head: float, y0: float = 0.0) -> np.ndarray:
    h, y = min(head, e - s), 0.25
    if strand == "-":
        pts = [(s, 0), (s + h, y), (e, y), (e, -y), (s + h, -y)]
    else:
        pts = [(s, y), (e - h, y), (e, 0), (e - h, -y), (s, -y)]
    return np.asarray([(x, v + y0) for x, v in pts], dtype=float)


def _draw_orf(ax: Axes, pile: Axes, locus: Locus, orf: OrfMark, hits: list[tuple[Any, int, int]],
              head: float) -> list[object]:  # fmt: skip
    """The novel ORF arrow in its own band above the gene arrows, their labels and the ``+``
    marks (decision 11), its OSTIR start tick and its peptides' pileup."""
    s, e = locus.start0, locus.end0
    ax.add_patch(Polygon(_arrow(s, e, orf.strand, head, ORF_Y), closed=True, edgecolor=NOVEL_ORF,
                         facecolor=to_rgba(NOVEL_ORF, 0.35), linewidth=1.0,
                         zorder=3))  # fmt: skip
    if orf.tick is not None:
        ax.plot([orf.tick, orf.tick], [ORF_Y - 0.3, ORF_Y + 0.3], color=NOVEL_ORF, lw=1.0)
    keys = {fold_il(p) for p in orf.peptides}
    mine = [h for h in hits if fold_il(str(h[0].peptide)) in keys]
    for ps, pe, d in _pileup_runs(mine):
        pile.fill_between([ps, pe], 0, d, color=NOVEL_ORF, linewidth=0, zorder=3)
    return ["novel_orf", locus.locus_id, locus.contig, s, e, orf.strand, None, None]


def _gene_style(g: Any, cmap: Any, norm: Normalize, measured: bool) -> tuple[str, dict[str, Any]]:
    """The legend key and patch keywords of one gene arrow; every face equals its key's swatch.

    Detection status comes first: single peptide is white outlined ``UNESTIMATED`` (also when
    imputed); not detected is ``style.LOCUS_UNDETECTED``. Imputed genes are the one exception
    to "a face is its status": grey and cross-hatched whatever their status, since the value
    shown is not measured. A detected gene with a value in an abundance run (``measured``) is
    cividis; every other detected gene is ``DETECTED_PLAIN``."""
    log2, status = _num(g.abundance_log2), str(g.detection_status)
    imputed = bool(g.abundance_imputed) if g.abundance_imputed is not pd.NA else False
    detected = bool(g.detected) if g.detected is not pd.NA else False
    kw: dict[str, Any] = {"edgecolor": "#333333", "linewidth": 0.3, "hatch": None}
    if status == "single_peptide_unestimated":
        return "unestimated", {**kw, "facecolor": "white", "edgecolor": UNESTIMATED,
                               "linewidth": 0.8}  # fmt: skip
    if imputed:
        return "imputed", {**kw, "facecolor": UNDETECTED, "hatch": IMPUTED_HATCH}
    if status in NOT_DETECTED:
        return "not_detected", {**kw, **style.LOCUS_UNDETECTED, "linewidth": 0.5}
    if measured and detected and log2 is not None:
        return "measured", {**kw, "facecolor": cmap(norm(log2))}
    return "detected", {**kw, "facecolor": DETECTED_PLAIN}


def _draw_genes(ax: Axes, rows: list[tuple[Any, int, int]], scale: tuple[float, float],
                head: float, measured: bool = True,
                ) -> tuple[list[list[object]], set[str]]:  # fmt: skip
    """The gene arrows and their data rows, plus the legend keys of the styles drawn."""
    cmap, norm = matplotlib.colormaps[style.LOCUS_ABUNDANCE], Normalize(*scale, clip=True)
    data: list[list[object]] = []
    drawn: set[str] = set()
    step = max(1, math.ceil(len(rows) / MAX_LABELS))
    for i, (g, s, e) in enumerate(sorted(rows, key=lambda t: (t[1], str(t[0].locus_tag)))):
        log2, status = _num(g.abundance_log2), str(g.detection_status)
        imputed = bool(g.abundance_imputed) if g.abundance_imputed is not pd.NA else False
        key, kw = _gene_style(g, cmap, norm, measured)
        drawn.add(key)
        strand = str(g.strand)
        ax.add_patch(Polygon(_arrow(s, e, strand, head), closed=True, **kw))
        if i % step == 0:  # ceil(n / step) <= MAX_LABELS
            name = g.gene if isinstance(g.gene, str) and g.gene else str(g.locus_tag)
            ax.text((s + e) / 2, 0.32, name, ha="center", va="bottom", fontsize=6, rotation=30,
                    gid=GENE_LABEL)  # fmt: skip
        value, category = (None, "imputed") if imputed else (log2, status)
        face = kw["facecolor"]
        fill = "white" if face == "white" else to_hex(face)
        data.append(["gene", str(g.locus_tag), str(g.contig), s, e, strand, value, category,
                     fill])  # fmt: skip
    return data, drawn


def _thin_labels(fig: Figure, ax: Axes) -> None:
    """Hide each gene label whose laid-out box meets the last one kept, left to right (sim100k:
    adjacent short genes collided). Layout is repeated until no label is hidden."""
    for _ in range(3):
        fig.draw_without_rendering()
        renderer = fig._get_renderer()  # type: ignore[attr-defined]
        labels = [t for t in ax.texts if t.get_gid() == GENE_LABEL and t.get_visible()]
        boxes = sorted(((t.get_window_extent(renderer), t) for t in labels),
                       key=lambda bt: (bt[0].x0, bt[1].get_text()))  # fmt: skip
        hidden, right = 0, -math.inf
        for box, t in boxes:
            if box.x0 < right:
                t.set_visible(False)
                hidden += 1
            else:
                right = box.x1
        if not hidden:
            return


def _pileup_runs(rows: list[tuple[Any, int, int]]) -> list[tuple[int, int, int]]:
    """Maximal runs ``(start, end, depth)`` of constant non-zero depth."""
    delta: dict[int, int] = {}
    for _, s, e in rows:
        delta[s] = delta.get(s, 0) + 1
        delta[e] = delta.get(e, 0) - 1
    runs: list[tuple[int, int, int]] = []
    depth, prev = 0, 0
    for pos in sorted(p for p, d in delta.items() if d):
        if depth:
            if runs and runs[-1][1] == prev and runs[-1][2] == depth:
                runs[-1] = (runs[-1][0], pos, depth)
            else:
                runs.append((prev, pos, depth))
        depth += delta[pos]
        prev = pos
    return runs


def _draw_marks(ax: Axes, rows: list[tuple[Any, int, int]], track: str,
                width: float) -> list[list[object]]:  # fmt: skip
    data: list[list[object]] = []
    for r, s, e in rows:
        strand = str(r.strand)
        sign = -1.0 if strand == "-" else 1.0
        y0, y1 = 0.3 * sign, 0.7 * sign
        if track == "promoter":
            x = float(s if strand == "-" else e)
            dx = 0.02 * width * sign
            ax.plot([x, x, x + dx], [y0, y1, y1], color=PROMOTER, lw=0.8)
            ax.plot([x + dx], [y1], marker=">" if sign > 0 else "<", color=PROMOTER, ms=3)
        else:
            x = (s + e) / 2
            ax.plot([x, x], [y0, y1], color=TERMINATOR, lw=0.8)
            ax.plot([x], [y1], marker="o", color=TERMINATOR, ms=3)
        data.append([track, str(r.name), str(r.chrom), s, e, strand, _num(r.score), None])
    return data


def _bed(view: RunView, key: str) -> pd.DataFrame:
    bed = _frame(view, key)
    return bed if bed.empty or "strand" in bed.columns else bed.assign(strand=".")


def render_locus(
    locus: Locus,
    view: RunView,
    cfg: GippaConfig,
    scale: tuple[float, float],
    orf: OrfMark | None = None,
) -> base.Instance:
    """One V3 instance: the locus window drawn on a fresh ``Figure``; data rows 0-based.
    ``orf`` (V13) adds the novel ORF track."""
    length, circular = _replicon(view, locus.contig)
    window = _window(locus, length, circular)
    width = float(window[1] - window[0])

    def near(key: str, contig_col: str = "contig") -> list[tuple[Any, int, int]]:
        return _in_window(_frame(view, key), contig_col, locus.contig, window, length, circular)

    size: Literal["single", "double"] = "single" if width <= SINGLE_MAX_BP else "double"
    chars = int((style.SINGLE_IN if size == "single" else style.DOUBLE_IN) * 72 / (0.6 * 7))
    title = "\n".join(textwrap.fill(c, chars) for c in locus.label.split("; "))
    extra = 0.0 if orf is None else ORF_EXTRA_IN
    fig = style.new_figure(size, HEIGHT_IN + extra + 0.12 * (title.count("\n") + 1) + 0.45,
                           layout="constrained")  # fmt: skip
    gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.8 if orf is None else 2.6])
    pile = fig.add_subplot(gs[0], label="pileup")
    ax = fig.add_subplot(gs[1], label="genes", sharex=pile)
    rows: list[list[object]] = []
    for r, s, e in near("regions_tsv"):
        colour = style.REGION.get(str(r.region_type), style.OKABE_ITO["grey"])
        ax.axvspan(s, e, facecolor=colour, alpha=0.2, edgecolor="none", zorder=0)
        rows.append(["region", str(r.region_id), str(r.contig), s, e, None, _num(r.score),
                     str(r.region_type)])  # fmt: skip
    for r, s, e in near("islands"):
        colour, hatch = style.ISLAND.get(str(r.direction), (style.OKABE_ITO["grey"], "..."))
        ax.axvspan(s, e, fill=False, hatch=hatch, edgecolor=colour, linewidth=0, alpha=0.5,
                   zorder=0)  # fmt: skip
        rows.append(["island", str(r.island_id), str(r.contig), s, e, None,
                     _num(r.mean_abundance_log2), str(r.direction)])  # fmt: skip
    label = style.run_abundance_label(view.table("genes_tsv"))
    gene_rows, drawn = _draw_genes(ax, near("genes_tsv"), scale, 0.015 * width, label is not None)
    rows += gene_rows
    hits = near("peptide_hits")
    if orf is not None:
        rows.append(_draw_orf(ax, pile, locus, orf, hits, 0.015 * width))
    pile_runs = _pileup_runs(hits)
    rows += [["pileup", None, locus.contig, s, e, None, float(d), None] for s, e, d in pile_runs]
    if pile_runs:
        edges = sorted({p for s, e, _ in pile_runs for p in (s, e)})
        depth = {s: d for s, _, d in pile_runs}
        pile.stairs([depth.get(p, 0) for p in edges[:-1]], edges, fill=True,
                    color=style.OKABE_ITO["blue"])  # fmt: skip
    else:
        note = "no peptide placements" if view.has("peptide_hits") else "peptide hits absent"
        pile.text(0.5, 0.5, note, ha="center", va="center", transform=pile.transAxes, fontsize=6)
    for key, track in (("promoters_bed", "promoter"), ("terminators_bed", "terminator")):
        marks = _in_window(_bed(view, key), "chrom", locus.contig, window, length, circular)
        rows += _draw_marks(ax, marks, track, width)
    if circular and length:
        for k in range(math.floor(window[0] / length), math.floor(window[1] / length) + 1):
            if window[0] < k * length < window[1]:
                ax.axvline(k * length, color=style.OKABE_ITO["black"], ls="--", lw=0.6)
    ax.xaxis.set_major_locator(MaxNLocator(nbins=MAX_TICKS - 1, steps=[1, 2, 2.5, 5, 10]))
    ax.xaxis.set_major_formatter(wrapped_formatter(length if circular and length else None))
    ax.set_xlim(*window)
    ax.set_ylim(-1.0, 1.0 if orf is None else ORF_TOP)
    ax.set_yticks([])
    ax.set_xlabel(f"{locus.contig} position (kb)")
    pile.set_ylabel("peptides")
    pile.tick_params(labelbottom=False)
    fig.suptitle(title, fontsize=7)
    _colour_key(fig, [pile, ax], scale, label, drawn)
    if orf is not None:
        _orf_key(pile)
    _thin_labels(fig, ax)
    rows = [r + [None] * (len(COLUMNS) - len(r)) for r in rows]
    data = pd.DataFrame(rows, columns=list(COLUMNS)) if rows else _empty(COLUMNS)
    data = data.astype({"start": "int64", "end": "int64", "value": float})
    return base.Instance(locus.locus_id, fig, data, locus.label, ("start",))


def _colour_key(fig: Figure, axes: list[Axes], scale: tuple[float, float], label: str | None,
                drawn: set[str]) -> None:  # fmt: skip
    """A horizontal cividis colorbar over the run-wide ``scale`` titled with the run's abundance
    ``label`` (none without abundance) plus the no-colour swatches, each the style
    ``_gene_style`` draws: not detected, imputed and single peptide always, the plain detected
    swatch only when such a gene is ``drawn``."""
    u = style.LOCUS_UNDETECTED
    swatches = [Patch(facecolor=u["facecolor"], edgecolor=u["edgecolor"], hatch=u["hatch"],
                      lw=0.5, label="not detected"),
                Patch(facecolor=UNDETECTED, edgecolor="#333333", lw=0.3, hatch=IMPUTED_HATCH,
                      label="imputed"),
                Patch(facecolor="white", edgecolor=UNESTIMATED, lw=0.8,
                      label="single peptide, unestimated")]  # fmt: skip
    if "detected" in drawn:
        plain = "detected" if label is None else "detected, no abundance value"
        swatches.insert(0, Patch(facecolor=DETECTED_PLAIN, edgecolor="#333333", lw=0.3,
                                 label=plain))  # fmt: skip
    if label is not None:
        mappable = ScalarMappable(Normalize(*scale), matplotlib.colormaps[style.LOCUS_ABUNDANCE])
        cbar = fig.colorbar(mappable, ax=axes, location="bottom", shrink=0.5, aspect=30, pad=0.02)
        cbar.set_label(label, fontsize=6)
        cbar.ax.tick_params(labelsize=6)
        cbar.dividers.set_cmap(mappable.cmap)  # the (empty) edge lines share the key's scale
        cbar.dividers.set_norm(mappable.norm)
    fig.legend(handles=swatches, loc="outside lower center", ncols=3, fontsize=6, frameon=False)


def _orf_key(pile: Axes) -> None:
    """The one-entry novel ORF key (V13 decision 13; the BGC shading shares its colour), just
    above the pileup axes so the constrained layout makes room for it under the title."""
    handle = Patch(facecolor=to_rgba(NOVEL_ORF, 0.35), edgecolor=NOVEL_ORF, lw=1.0, label=ORF_KEY)
    pile.legend(handles=[handle], loc="lower center", bbox_to_anchor=(0.5, 1.0), fontsize=6,
                frameon=False, borderaxespad=0.2)  # fmt: skip


def _empty(columns: Iterable[str]) -> pd.DataFrame:
    return pd.DataFrame({c: pd.Series([], dtype=object) for c in columns})


def index_frame(loci: Sequence[Locus], n_rendered: int) -> pd.DataFrame:
    """The locus index in display order: ``rank`` 1..n, the family's ``estimate`` name and value
    (``error``), and whether a locus view was drawn (0-based ``start``)."""
    rows = [
        [i + 1, lo.locus_id, lo.kind, lo.test_family, lo.contig, lo.start0, lo.end0,
         ESTIMATE.get(lo.test_family, "BH q"), math.nan if lo.error is None else lo.error,
         i < n_rendered]
        for i, lo in enumerate(loci)
    ]  # fmt: skip
    if not rows:
        return _empty(INDEX_COLUMNS)
    return pd.DataFrame(rows, columns=list(INDEX_COLUMNS))


def _n_rendered(loci: Sequence[Locus], cfg: GippaConfig) -> int:
    return max(0, min(cfg.report.max_locus_plots, len(loci)))


@base.data_table("V03")
def locus_index(view: RunView, cfg: GippaConfig) -> base.DataTable:
    """``data/V03_locus.tsv``: every selected locus in display order (no image, P-12)."""
    loci = select_loci(view, cfg)
    return base.DataTable(None, index_frame(loci, _n_rendered(loci, cfg)), ("start",))


@base.plot("V03", "locus", needs=["genes_tsv"], optional=[*_OPTIONAL, *_OPTIONAL_TRACKS])
def locus_plots(view: RunView, cfg: GippaConfig) -> list[base.Instance] | base.Skip:
    """Every rendered locus instance; none when no locus is selected (the index table, written
    from ``locus_index``, says so)."""
    loci = select_loci(view, cfg)
    scale = run_scale(view.table("genes_tsv"))
    return [render_locus(lo, view, cfg, scale) for lo in loci[: _n_rendered(loci, cfg)]]
