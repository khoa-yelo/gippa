"""V6 ``mapping_fdr``: can I trust tolerant hits and de novo peptides? (GIPPA_SPEC.md §10.2 V6;
M9-05). Draws the M3 ``mapping_fdr.tsv`` estimates; computes none.

(a) targets vs. decoys per ``mapping.length_classes`` class x ``n_diff`` (log scale, 0 at the
floor); (b) the per-m ``fdr_stratum`` (``length_class = all``) with a dashed line at
``stats.alpha.placement`` and ``n_accepted`` above each point, NA strata as open markers at the
floor; (c) de novo precision per score bin with its CI, only with ``denovo_precision`` rows.
Panels (a) and (b) use the ``reverse`` rows, or the ``shuffle`` rows when only those have targets
(M9-05 lead decision 3). The title is built from the panel (b) rows (index decision 5).

Degenerate states (M12-04 PR 4, P-02/P-03): no tolerant targets and no de novo rows is one empty
state; no tolerant targets with de novo rows collapses (a) and (b) to one text line; when every
de novo bin is ``c_assumed``, (c) draws no points (the ``drawn`` data column says so) and states
the assumed c; NA strata of (b) share one note.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from matplotlib.axes import Axes

from gippa.mapping.placement_fdr import length_class_labels
from gippa.viz import base, style

if TYPE_CHECKING:
    from gippa.config import GippaConfig
    from gippa.io.run_view import RunView

COLUMNS = (
    "panel",
    "length_class",
    "n_diff",
    "score_bin",
    "metric",
    "value",
    "ci_low",
    "ci_high",
    "flags",
    "drawn",
)
TARGET = style.OKABE_ITO["blue"]
DECOY = style.OKABE_ITO["vermillion"]
FOOTNOTE = "c assumed 0.99 (no reference peptides)"
EMPTY = "no peptide placements"
#: the whole-figure empty state (``base.empty_state``) when ``mapping_fdr`` has no rows at all
EMPTY_TITLE = "No peptide placements or de novo precision in this run"
#: the tolerant route had no targets (P-02): the empty state and the collapsed (a)/(b) line
NOT_NEEDED = "the tolerant-route FDR is not needed"
#: panel (c) when no bin has a precision value (and the bins are not ``c_assumed``)
NO_PRECISION = "no de novo precision values"
#: panel (c) when every bin is ``c_assumed`` (P-03); ``{c}`` is the assumed value
ASSUMED = "De novo precision not measured (no reference peptides); assumed c = {c}"
_FLOOR_A = 0.5  # panel (a) log-axis floor: 0 counts are drawn here


def _py(value: object) -> object:
    """A plain Python cell: ``None`` for NA, ``int`` / ``float`` for numbers."""
    if value is None or value is pd.NA or (isinstance(value, float) and math.isnan(value)):
        return None
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value)
    return value


def _num(value: object) -> float:
    v = _py(value)
    return math.nan if v is None else float(str(v))


def placement_rows(fdr: pd.DataFrame) -> pd.DataFrame:
    """The ``placement`` rows of one decoy type: ``reverse`` unless only ``shuffle`` has targets."""
    rows = fdr[fdr["section"] == "placement"]
    by_type = [rows[rows["decoy_type"] == t] for t in ("reverse", "shuffle")]
    with_targets = [r for r in by_type if r["targets"].fillna(0).sum() > 0]
    present = [r for r in by_type if len(r)]
    return (with_targets or present or [rows.iloc[0:0]])[0]


def _row(
    panel: str,
    lc: object,
    m: object,
    score_bin: object,
    metric: str,
    value: object,
    lo: object = None,
    hi: object = None,
    flags: object = None,
) -> list[object]:
    return [panel, lc, m, score_bin, metric, *(_py(v) for v in (value, lo, hi, flags)), True]


def _assumed(flags: object) -> bool:
    return "c_assumed" in str(_py(flags) or "").split(";")


def all_assumed(c: pd.DataFrame) -> bool:
    """Every panel (c) row is ``c_assumed`` (P-03: precision not measured)."""
    return bool(len(c)) and all(_assumed(f) for f in c["flags"])


def exact_only(fdr: pd.DataFrame) -> bool:
    """Placement rows exist but none has tolerant-route targets (P-02)."""
    rows = fdr[fdr["section"] == "placement"]
    return (
        bool(len(rows))
        and not pd.to_numeric(rows["targets"], errors="coerce").fillna(0).gt(0).any()
    )


def exact_label(hits: pd.DataFrame | None) -> str:
    """``"All 774 peptides placed by exact match"`` from ``peptide_hits`` (distinct peptides)."""
    if hits is None or "peptide" not in hits.columns or not len(hits):
        return "All peptides placed by exact match"
    n = int(hits["peptide"].dropna().astype(str).nunique())
    noun = "peptide" if n == 1 else "peptides"
    return f"All {style.fmt_count(n)} {noun} placed by exact match"


def data_frame(fdr: pd.DataFrame, classes: Sequence[str]) -> pd.DataFrame:
    """The V06 data TSV: panel (a), (b), (c) rows sorted as in M9-05."""
    rows: list[list[object]] = []
    place = placement_rows(fdr)
    per_class = place[place["length_class"] != "all"]
    extra = [c for c in dict.fromkeys(per_class["length_class"].astype(str)) if c not in classes]
    strata = sorted({int(m) for m in place["n_diff"].dropna()})
    for lc in [*classes, *extra]:
        for m in strata:
            hit = per_class[(per_class["length_class"] == lc) & (per_class["n_diff"] == m)]
            for metric in ("decoys", "targets"):
                value = hit[metric].iloc[0] if len(hit) else None
                rows.append(_row("a", lc, m, None, metric, value))
    totals = place[place["length_class"] == "all"]
    for m in strata:
        hit = totals[totals["n_diff"] == m]
        for metric in ("fdr_stratum", "n_accepted"):
            rows.append(
                _row("b", "all", m, None, metric, hit[metric].iloc[0] if len(hit) else None)
            )
    denovo = fdr[fdr["section"] == "denovo_precision"].sort_values("score_bin", kind="mergesort")
    for r in denovo.itertuples(index=False):
        rows.append(
            _row(
                "c",
                None,
                None,
                _py(r.score_bin),
                "precision",
                r.precision,
                r.ci_low,
                r.ci_high,
                r.flags,
            )
        )
    data = pd.DataFrame(rows, columns=list(COLUMNS), dtype=object)
    c = data["panel"] == "c"
    if all_assumed(data.loc[c]):
        data.loc[c, "drawn"] = False  # P-03: nothing measured, the values are kept
    elif c.any():
        data.loc[c, "drawn"] = [_py(v) is not None for v in data.loc[c, "value"]]
    return data


def _strata(b: pd.DataFrame) -> list[tuple[int, float, int]]:
    """``(n_diff, fdr_stratum or NaN, n_accepted)`` per stratum of the panel (b) rows."""
    out = []
    for m in dict.fromkeys(b["n_diff"]):
        sub = b[b["n_diff"] == m].set_index("metric")["value"]
        acc = _num(sub.get("n_accepted"))
        out.append((int(m), _num(sub.get("fdr_stratum")), 0 if math.isnan(acc) else int(acc)))
    return out


def _join(strata: Sequence[int]) -> str:
    return ", ".join(str(m) for m in strata)


def title(data: pd.DataFrame, alpha: float, exact: str | None = None) -> str:
    """From the panel (b) rows: the strata with FDR ≤ alpha, the accepted count, the NA strata.

    ``exact`` (no tolerant targets, P-02) replaces the (b) summary with that line."""
    strata = _strata(data[data["panel"] == "b"])
    if exact is not None:
        c = data[data["panel"] == "c"]
        if len(c) and not c["drawn"].astype(bool).any():
            return f"{exact}; de novo precision not measured"
        bins = len(c)
        return f"{exact}; de novo precision in {base.plural(bins, 'score bin')}"
    if not strata:
        bins = int((data["panel"] == "c").sum())
        tail = f"; de novo precision in {base.plural(bins, 'score bin')}" if bins else ""
        return f"Placement FDR unestimated: no peptide placements{tail}"
    ok = [m for m, f, _ in strata if not math.isnan(f) and f <= alpha]
    na = [m for m, f, _ in strata if math.isnan(f)]
    accepted = f"{sum(n for _, _, n in strata)} FDR-tested placements accepted"
    if ok:
        parts = [f"Placement FDR ≤ {style.fmt_alpha(alpha)} in strata {_join(ok)}", accepted]
    elif len(na) < len(strata):
        parts = [f"No stratum with placement FDR ≤ {style.fmt_alpha(alpha)}", accepted]
    else:
        parts = [accepted]
    if na:
        noun = "stratum" if len(na) == 1 else "strata"
        parts.append(f"{noun} {_join(na)} unestimated (no targets)")
    return "; ".join(parts)


def _empty(ax: Axes, text: str = EMPTY) -> None:
    ax.set_axis_off()
    ax.text(0.5, 0.5, text, ha="center", va="center", transform=ax.transAxes, color="#555555")


def _panel_a(ax: Axes, a: pd.DataFrame) -> None:
    classes = list(dict.fromkeys(a["length_class"]))
    strata = list(dict.fromkeys(a["n_diff"]))
    width = 0.8 / (2 * len(strata))
    ax.set_yscale("log")
    top = _FLOOR_A
    minor_x, minor_lab = [], []
    for i, lc in enumerate(classes):
        for j, m in enumerate(strata):
            x0 = i - 0.4 + 2 * j * width
            minor_x.append(x0 + width)
            minor_lab.append(str(m))
            for k, (metric, colour) in enumerate((("targets", TARGET), ("decoys", DECOY))):
                hit = a[(a["length_class"] == lc) & (a["n_diff"] == m) & (a["metric"] == metric)]
                v = _num(hit["value"].iloc[0]) if len(hit) else math.nan
                h = max(v, _FLOOR_A) - _FLOOR_A if not math.isnan(v) else 0.0
                top = max(top, _FLOOR_A + h)
                label = metric if (i, j) == (0, 0) else None
                ax.bar(
                    x0 + (k + 0.5) * width,
                    h,
                    width,
                    bottom=_FLOOR_A,
                    color=colour,
                    label=label,
                    linewidth=0,
                )
    ax.set_ylim(_FLOOR_A, max(10.0, top * 2))
    ax.set_xticks(minor_x, minor_lab)
    ax.tick_params(axis="x", length=2)
    for i, lc in enumerate(classes):
        ax.text(i, -0.13, lc, transform=ax.get_xaxis_transform(), ha="center", va="top")
    ax.set_xlabel("differences per length class (aa)", labelpad=12)
    ax.set_ylabel("placements (0 at the floor)")
    ax.legend(frameon=False, loc="upper left")
    ax.set_title("(a) targets vs. decoys", loc="left")


def _panel_b(ax: Axes, b: pd.DataFrame, alpha: float) -> None:
    strata = _strata(b)
    finite = [f for _, f, _ in strata if not math.isnan(f)]
    ax.set_xlim(-0.6, len(strata) - 0.4 if strata else 0.6)
    ax.set_ylim(0.0, 1.3 * max([alpha, *finite]))
    ax.axhline(alpha, linestyle="--", color="#555555", linewidth=0.8)
    est = [(m, f) for m, f, _ in strata if not math.isnan(f)]
    if est:
        ax.plot(
            [m for m, _ in est],
            [f for _, f in est],
            marker="o",
            linestyle="none",
            color=TARGET,
            markersize=4,
        )
    na = [m for m, f, _ in strata if math.isnan(f)]
    if na:
        ax.plot(
            na,
            [0.0] * len(na),
            marker="o",
            linestyle="none",
            markerfacecolor="none",
            markeredgecolor="#999999",  # Okabe-Ito grey (§10.1)
            markersize=4,
            clip_on=False,
        )
    for m, f, n in strata:
        y = 0.0 if math.isnan(f) else f
        ax.annotate(
            str(n), (m, y), xytext=(0, 4), textcoords="offset points", ha="center", fontsize=6
        )
    if na:  # one stacked note, not one overlapping label per stratum (P-02)
        noun = "stratum" if len(na) == 1 else "strata"
        ax.text(
            0.02,
            0.97,
            f"NA (no targets):\n{noun} {_join(na)}",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=6,
            color="#555555",
        )
    ax.set_xticks([m for m, _, _ in strata])
    ax.set_xlabel("differences (n_diff)")
    ax.set_ylabel("stratum FDR")
    ax.set_title(f"(b) stratum FDR (dashed: {style.fmt_alpha(alpha)})", loc="left")


def _c_value(fdr: pd.DataFrame) -> str:
    """The assumed c of the de novo rows (``c_value``), 0.99 when not recorded."""
    rows = fdr[fdr["section"] == "denovo_precision"]
    values = pd.to_numeric(rows.get("c_value", pd.Series([], dtype=float)), errors="coerce")
    finite = values.dropna()
    return base.sig(float(finite.iloc[0]) if len(finite) else 0.99, 3)


def _panel_c_assumed(ax: Axes, c_value: str) -> None:
    ax.set_axis_off()
    ax.set_title("(c) de novo precision", loc="left")
    text = ASSUMED.format(c=c_value).replace("; ", ";\n").replace(" (no", "\n(no")
    ax.text(0.0, 0.6, text, ha="left", va="center", transform=ax.transAxes, color="#555555")


def _panel_c(ax: Axes, c: pd.DataFrame) -> bool:
    """Draw panel (c); ``True`` when a row is ``c_assumed`` (the footnote is needed)."""
    x = np.array([_num(v) for v in c["score_bin"]])
    y = np.array([_num(v) for v in c["value"]])
    lo = np.array([_num(v) for v in c["ci_low"]])
    hi = np.array([_num(v) for v in c["ci_high"]])
    assumed = np.array(["c_assumed" in str(f).split(";") for f in c["flags"]], dtype=bool)
    ci = np.isfinite(y) & np.isfinite(lo) & np.isfinite(hi) & (hi > lo)
    if ci.any():
        ax.errorbar(
            x[ci],
            y[ci],
            yerr=[y[ci] - lo[ci], hi[ci] - y[ci]],
            fmt="none",
            ecolor="#000000",
            elinewidth=0.8,
            capsize=0,
        )
    for hollow in (False, True):
        sel = np.isfinite(y) & (assumed == hollow)
        if sel.any():
            ax.plot(
                x[sel],
                y[sel],
                marker="o",
                linestyle="none",
                markersize=4,
                markeredgecolor=TARGET,
                markerfacecolor="none" if hollow else TARGET,
            )
    vals = np.concatenate([y[np.isfinite(y)], lo[np.isfinite(lo)]])
    ax.set_ylim(max(0.0, float(vals.min()) - 0.05) if len(vals) else 0.0, 1.01)
    ax.set_xticks([int(v) for v in x if math.isfinite(v)])
    ax.set_xlabel("de novo score bin")
    ax.set_ylabel("precision (CI)" if ci.any() else "precision")
    ax.set_title("(c) de novo precision", loc="left")
    return bool(assumed.any())


@base.plot("V06", "mapping_fdr", needs=("mapping_fdr",), optional=("peptide_hits",))
def mapping_fdr(view: RunView, cfg: GippaConfig) -> list[base.Instance] | base.Skip:
    """V6: placement targets/decoys, stratum FDR against alpha, de novo precision."""
    fdr = view.table("mapping_fdr")
    assert fdr is not None  # guaranteed by ``needs``
    alpha = float(cfg.stats.alpha.placement)
    classes = length_class_labels([list(c) for c in cfg.mapping.length_classes])
    data = data_frame(fdr, classes)
    if not len(data):
        return [base.empty_state("V06", "mapping_fdr", EMPTY_TITLE, COLUMNS, "double")]
    exact = exact_label(view.table("peptide_hits")) if exact_only(fdr) else None
    if exact is not None and not (data["panel"] == "c").any():
        msg = f"{exact}; {NOT_NEEDED}."
        return [base.empty_state("V06", "mapping_fdr", msg, COLUMNS, "double")]
    if exact is not None:  # (a) and (b) are one text line: their rows are kept, not drawn
        data.loc[data["panel"].isin(["a", "b"]), "drawn"] = False
    sub = {p: data[data["panel"] == p] for p in ("a", "b", "c")}
    text_only = exact is not None and not sub["c"]["drawn"].astype(bool).any()
    fig = style.new_figure("double", 1.6 if text_only else 2.8)
    fig.set_label("V06_mapping_fdr")
    if exact is not None:  # (a) and (b) collapse to one text line (P-02)
        panels, ratios = ["ab", "c"], [2.6, 1.6]
    else:
        panels = [p for p in ("a", "b", "c") if len(sub[p]) or p != "c"]
        ratios = [2.6, 1.5, 1.3][: len(panels)]
    grid = fig.add_gridspec(
        1,
        len(panels),
        width_ratios=ratios,
        left=0.08,
        right=0.98,
        bottom=0.24,
        top=0.8,
        wspace=0.35,
    )
    axes = {p: fig.add_subplot(grid[0, i], label=p) for i, p in enumerate(panels)}
    if exact is not None:
        axes["ab"].set_axis_off()
        line = f"{exact}:\n{NOT_NEEDED}."
        axes["ab"].text(0.0, 0.6, line, ha="left", va="center", transform=axes["ab"].transAxes)
    elif len(sub["a"]) and len(sub["b"]):
        _panel_a(axes["a"], sub["a"])
        _panel_b(axes["b"], sub["b"], alpha)
        if (placement_rows(fdr)["decoy_type"] == "shuffle").any():  # lead decision 3 fallback
            for p in ("a", "b"):
                # on a second line: one line runs into the next panel's title (§14.4)
                axes[p].set_title(axes[p].get_title(loc="left") + "\n(shuffle decoys)", loc="left")
    else:
        _empty(axes["a"])
        _empty(axes["b"])
    if "c" in axes:
        if all_assumed(sub["c"]):
            _panel_c_assumed(axes["c"], _c_value(fdr))
        elif not sub["c"]["drawn"].astype(bool).any():
            _empty(axes["c"], NO_PRECISION)
        elif _panel_c(axes["c"], sub["c"]):
            fig.text(0.98, 0.02, FOOTNOTE, ha="right", va="bottom", fontsize=6, color="#555555")
    text = title(data, alpha, exact)
    fig.suptitle(text, y=0.97)
    return [base.Instance(None, fig, data, text)]
