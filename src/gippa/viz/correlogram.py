"""V10 ``correlogram``: is abundance (and detection) autocorrelated along the replicon?
(GIPPA_SPEC.md §10.2 V10, §7.3.5(c); ADR-0020; M9-09 and its lead decisions 2 and 4-7).

Draws ``tables/autocorrelation.tsv`` and recomputes nothing: per replicon (FASTA order) panel (a)
Moran's I and panel (b) the BB join statistic S_k against the lag, each with the gene-label (light
grey) and TU-block (``#56B4E9``) 2.5-97.5 % envelopes (``null_lo``-``null_hi``) and dashed null
means. A marker is filled when the TU-block row's ``q_value ≤ stats.alpha.context`` (the
conservative null ADR-0020 gates); NA rows, and rows flagged ``not_estimable`` or ``few_values``
whatever their value (lead decision 21), draw nothing, and a panel whose rows are all NA says
"not estimable (<flags>)". Contigs wrap into rows of at most 3 columns, each (b) panel under its
(a) panel (lead decision 23). ``oriter_slope`` rows are context-table rows and are not drawn.
A ``join_count_bb`` row's ``value`` is the raw count BB_k while its ``null_*`` and ``effect`` are on
the Pearson-residual statistic S_k (ADR-0020 Amendment 1), so (b) plots S_k = ``effect +
null_mean`` against the envelopes (docs/specs/stat-m5-joincount-scale.md).

When no ``morans_i`` row is drawable (every one NA or flagged; M12-04 PR 4, P-06) the plot is
one empty state naming the detected genes ("Too few detected genes (4 of 158) to estimate
autocorrelation."). Flags in partial panels are written in plain words (``words.plain``).
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from matplotlib.artist import Artist
from matplotlib.axes import Axes

from gippa.viz import base, style
from gippa.viz.words import plain

if TYPE_CHECKING:
    from gippa.config import GippaConfig
    from gippa.io.run_view import RunView

COLUMNS = ("contig", "statistic", "null_model", "lag", "value", "null_mean", "null_lo",
           "null_hi", "effect", "q_value", "significant", "flags")  # fmt: skip
STATISTICS = ("morans_i", "join_count_bb")  # schema enum order
NULLS = ("gene_label", "tu_block")
_NUM = ("value", "null_mean", "null_lo", "null_hi", "effect", "q_value")
#: null -> (band colour, band alpha, mean-line colour)
NULL_STYLE = {"gene_label": ("#CCCCCC", 0.6, "#777777"), "tu_block": ("#56B4E9", 0.3, "#0072B2")}
NULL_LABEL = {"gene_label": "gene-label", "tu_block": "TU-block"}
YLABEL = {"morans_i": "Moran's I", "join_count_bb": "BB join statistic (Pearson residual)"}
NOT_ESTIMABLE = "Autocorrelation not estimable on this genome"
#: the whole-plot empty state when no Moran's I row is drawable (P-06)
TOO_FEW = "Too few detected genes{count} to estimate autocorrelation."
GREY = "#555555"
NA_FLAGS = frozenset({"not_estimable", "few_values"})
_COLS, _HEAD_IN, _FOOT_IN, _ROW_IN = 3, 0.45, 0.62, 3.4


def data_frame(ac: pd.DataFrame, alpha: float) -> pd.DataFrame:
    """Every ``morans_i`` / ``join_count_bb`` row (NA rows included), sorted by contig,
    statistic (enum order), null (enum order) and lag; ``significant`` = own ``q ≤ alpha``."""
    sel = ac[ac["statistic"].astype(str).isin(STATISTICS)].copy()
    sel["_s"] = [STATISTICS.index(str(s)) for s in sel["statistic"]]
    sel["_n"] = [NULLS.index(str(n)) if str(n) in NULLS else len(NULLS) for n in sel["null_model"]]
    sel["contig"] = sel["contig"].astype(str)
    sel = sel.sort_values(["contig", "_s", "_n", "lag"], kind="stable")
    q = np.asarray(pd.to_numeric(sel["q_value"]), dtype=float)
    out = pd.DataFrame(
        {
            "contig": sel["contig"].to_numpy(dtype=object),
            "statistic": sel["statistic"].astype(str).to_numpy(dtype=object),
            "null_model": sel["null_model"].astype(str).to_numpy(dtype=object),
            "lag": pd.array(sel["lag"], dtype="Int64"),
            **{c: np.asarray(pd.to_numeric(sel[c]), dtype=float) for c in _NUM},
            "significant": np.isfinite(q) & (q <= alpha),
            "flags": [None if pd.isna(f) else str(f) for f in sel["flags"]],
        },
        columns=list(COLUMNS),
    )
    return out.reset_index(drop=True)


def plotted(
    statistic: Sequence[object], value: np.ndarray, effect: np.ndarray, null_mean: np.ndarray
) -> np.ndarray:
    """The drawn statistic: Moran's I ``value``; for ``join_count_bb`` rows S_k = ``effect +
    null_mean`` (the envelopes' scale; ``value`` is the raw BB_k)."""
    bb = np.asarray([str(s) == "join_count_bb" for s in statistic], dtype=bool)
    return np.where(bb, effect + null_mean, value)


def drawn(data: pd.DataFrame) -> pd.DataFrame:
    """``data`` with the flagged rows' numbers set to NA: what the figure and title read."""
    na = [not NA_FLAGS.isdisjoint(_split(f)) for f in data["flags"]]
    out = data.copy()
    out.loc[na, list(_NUM)] = math.nan
    out.loc[na, "significant"] = False
    return out


def _split(flags: object) -> list[str]:
    """The tokens of a ``flags`` cell; NA (``None``, NaN, ``pd.NA``) is no flag."""
    return [t for t in flags.split(";") if t] if isinstance(flags, str) else []


def too_few(genes: pd.DataFrame | None) -> str:
    """The P-06 empty-state title, with "(<detected> of <genes>)" when ``genes.tsv`` is read."""
    if genes is None or "detected" not in genes.columns:
        return TOO_FEW.format(count="")
    detected = int(genes["detected"].astype("boolean").fillna(False).sum())
    return TOO_FEW.format(count=f" ({style.fmt_count(detected)} of {style.fmt_count(len(genes))})")


def contig_order(data: pd.DataFrame, fasta: Sequence[str]) -> list[str]:
    present = set(data["contig"])
    return [c for c in fasta if c in present] + sorted(present - set(fasta))


def _lag1(data: pd.DataFrame, contig: str, null: str) -> pd.Series | None:
    hit = data[
        (data["contig"] == contig)
        & (data["statistic"] == "morans_i")
        & (data["null_model"] == null)
        & (data["lag"] == 1)
    ]
    if not len(hit):
        return None
    row = hit.iloc[0]
    finite = all(math.isfinite(float(row[c])) for c in ("value", "null_lo", "null_hi", "q_value"))
    return row if finite else None


def lag1_moran(data: pd.DataFrame, order: Sequence[str]) -> tuple[str, pd.Series] | None:
    """The lag-1 ``morans_i`` row the V02 and V10 titles both name, with its null: the first
    usable row (``value``, ``q_value``, ``null_lo`` and ``null_hi`` all finite) trying the
    TU-block null first, then the gene-label null (old runs), each over the contigs in
    ``order`` (FASTA order). ``data`` is ``drawn(data_frame(...))``; ``None`` if no row is usable
    (M12-04 PR 5 review 1, M1)."""
    for null in NULLS[::-1]:
        for contig in order:
            row = _lag1(data, contig, null)
            if row is not None:
                return null, row
    return None


def title(data: pd.DataFrame, order: Sequence[str]) -> str:
    """The lag-1 Moran's I of ``lag1_moran``'s row against that row's envelope and q, naming
    the null (lead decision 4; M12-04 P-04)."""
    found = lag1_moran(data, order)
    if found is None:
        return NOT_ESTIMABLE
    null, row = found
    v, lo, hi, q = (float(row[c]) for c in ("value", "null_lo", "null_hi", "q_value"))
    word = "above" if v > hi else "below" if v < lo else "within"
    name = NULL_LABEL[null]
    return (f"Lag-1 Moran's I {style.fmt_num(v)}: {word} the {name} envelope "
            f"({style.q_equals(q)}, {name} null)")  # fmt: skip


def _flags(rows: pd.DataFrame) -> str:
    """The flags of ``rows`` in plain words: those that make a row NA when any, else all."""
    seen: dict[str, None] = {}
    for f in rows["flags"]:
        seen.update(dict.fromkeys(_split(f)))
    tokens = [t for t in seen if t in NA_FLAGS] or list(seen)
    return "; ".join(plain(t) for t in tokens) or "no estimable value"


def _panel(ax: Axes, rows: pd.DataFrame, statistic: str, lags: tuple[int, int]) -> None:
    num = {c: rows[c].to_numpy(dtype=float) for c in ("value", "effect", "null_mean")}
    rows = rows.assign(y=plotted(list(rows["statistic"]), num["value"], num["effect"],
                                 num["null_mean"]))  # fmt: skip
    if not np.isfinite(rows["y"].to_numpy(dtype=float)).any():
        ax.set_axis_off()
        ax.text(0.5, 0.5, f"Not drawn: {_flags(rows)}", ha="center", va="center",
                transform=ax.transAxes, color=GREY, fontsize=6, wrap=True)  # fmt: skip
        return
    finite = np.isfinite(rows["y"].to_numpy(dtype=float))
    na_lags = sorted({int(v) for v in rows["lag"].dropna()} - {int(v) for v in
                     rows.loc[finite, "lag"].dropna()})  # fmt: skip
    if na_lags:
        na_rows = rows[rows["lag"].isin(na_lags)]
        noun = "lag" if len(na_lags) == 1 else "lags"
        note = f"{noun} {', '.join(str(v) for v in na_lags)} not drawn: {_flags(na_rows)}"
        ax.text(0.02, 0.02, note, ha="left", va="bottom", transform=ax.transAxes, color=GREY,
                fontsize=5, wrap=True)  # fmt: skip
    by_null = {n: rows[rows["null_model"] == n] for n in NULLS}
    for null, sub in by_null.items():
        band, a, line = NULL_STYLE[null]
        x = sub["lag"].to_numpy(dtype=float, na_value=np.nan)
        lo, hi, mean = (sub[c].to_numpy(dtype=float) for c in ("null_lo", "null_hi", "null_mean"))
        ok = np.isfinite(x) & np.isfinite(lo) & np.isfinite(hi)
        if ok.any():
            ax.fill_between(x[ok], lo[ok], hi[ok], color=band, alpha=a, linewidth=0,
                            label=f"{NULL_LABEL[null]} 95 % envelope")  # fmt: skip
        ok = np.isfinite(x) & np.isfinite(mean)
        if ok.any():
            ax.plot(x[ok], mean[ok], linestyle="--", color=line, linewidth=0.7)
    # the observed value and the deciding q: the TU-block row, else the gene-label row's value
    tu = by_null["tu_block"].set_index("lag")
    gl = by_null["gene_label"].set_index("lag")
    xs, ys, filled = [], [], []
    for lag in sorted({int(v) for v in rows["lag"].dropna()}):
        v = tu["y"].get(lag, math.nan)
        if not math.isfinite(v):
            v = gl["y"].get(lag, math.nan)
        if math.isfinite(v):
            xs.append(float(lag))
            ys.append(float(v))
            filled.append(bool(tu["significant"].get(lag, False)))
    ax.plot(xs, ys, color="#000000", linewidth=0.6)
    f = np.array(filled, dtype=bool)
    x_arr, y_arr = np.array(xs), np.array(ys)
    ax.plot(x_arr[f], y_arr[f], linestyle="none", marker="o", markersize=3.5, color="#000000",
            label="TU-block q ≤ alpha")  # fmt: skip
    ax.plot(x_arr[~f], y_arr[~f], linestyle="none", marker="o", markersize=3.5,
            markerfacecolor="none", markeredgecolor="#000000", markeredgewidth=0.8,
            label="TU-block q > alpha")  # fmt: skip
    lo, hi = lags
    ax.set_xlim(lo - 0.6, hi + 0.6)
    ax.set_xticks(sorted({lo, hi, *range(5 * math.ceil(lo / 5), hi + 1, 5)} - {0}))
    ax.set_xlabel("lag (genes)")
    ax.set_ylabel(YLABEL[statistic])


@base.plot("V10", "correlogram", needs=("autocorrelation",), optional=("genes_tsv",))
def correlogram(view: RunView, cfg: GippaConfig) -> list[base.Instance] | base.Skip:
    """V10: Moran's I and BB join-count correlograms with both null envelopes."""
    ac = view.table("autocorrelation")
    assert ac is not None  # guaranteed by ``needs``
    data = data_frame(ac, float(cfg.stats.alpha.context))
    if not len(data):
        return [base.empty_state("V10", "correlogram", NOT_ESTIMABLE, COLUMNS, "double")]
    fasta = [str(r["id"]) for r in (view.preflight() or {}).get("replicons", [])]
    order = contig_order(data, fasta)
    lags = (int(cfg.context.moran_lags[0]), int(cfg.context.moran_lags[-1]))
    shown = drawn(data)
    moran = shown[shown["statistic"] == "morans_i"]
    if not np.isfinite(moran["value"].to_numpy(dtype=float)).any():  # P-06
        msg = too_few(view.table("genes_tsv"))
        return [base.empty_state("V10", "correlogram", msg, COLUMNS, "double")]
    text = title(shown, order)
    n_rows, n_cols = -(-len(order) // _COLS), min(_COLS, len(order))
    height = _HEAD_IN + _FOOT_IN + _ROW_IN * n_rows
    fig = style.new_figure("double", height)
    fig.set_label("V10_correlogram")
    grid = fig.add_gridspec(2 * n_rows, n_cols, left=0.09, right=0.98, bottom=_FOOT_IN / height,
                            top=1 - _HEAD_IN / height, hspace=0.6, wspace=0.42)  # fmt: skip
    for k, contig in enumerate(order):
        for i, (panel, statistic) in enumerate(zip("ab", STATISTICS, strict=True)):
            ax = fig.add_subplot(grid[2 * (k // _COLS) + i, k % _COLS], label=f"{panel}:{contig}")
            rows = shown[(shown["contig"] == contig) & (shown["statistic"] == statistic)]
            _panel(ax, rows, statistic, lags)
            ax.set_title(f"({panel}) {contig}", loc="left")
    entries: dict[str, Artist] = {}  # every panel's entries, so an NA panel loses none
    for ax in fig.axes:
        for handle, label in zip(*ax.get_legend_handles_labels(), strict=True):
            entries.setdefault(label, handle)
    if entries:
        fig.legend(list(entries.values()), list(entries), frameon=False, loc="lower center",
                   ncol=4, fontsize=6, bbox_to_anchor=(0.5, 0.0))  # fmt: skip
    fig.suptitle(text, y=1 - 0.05 / height, va="top")
    return [base.Instance(None, fig, data, text)]
