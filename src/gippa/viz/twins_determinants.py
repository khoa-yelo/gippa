"""Plotly twin of V15 ``determinants`` (M9-20; GIPPA_SPEC.md §10.2 V15): the four panels from the
data TSV, (b) one sub-plot per PD feature; (a) dodged numeric y (circle elastic net, open diamond
LightGBM), reversed, "CI unestimated" marked; the static footnotes as an annotation (decision 5,
S3). Pure: no file I/O."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from gippa.determinants import CAVEATS
from gippa.viz.determinants import (
    A_HEADER,
    BANNER,
    BLOCKS,
    CAT_COLOUR,
    CAT_LABEL,
    COLUMNS,
    D_HEADER,
    DELTA_FOOT,
    GATED_SUFFIX,
    GREY,
    LEGEND,
    MODEL_COLOUR,
    MODELS,
    PD_FEATURES,
    PD_XLABEL,
)
from gippa.viz.interactive import interactive
from gippa.viz.style import OKABE_ITO
from gippa.viz.twins import DOT, Meta, abundance, ci, finish, floats, fmt, points, rgba, titles

_SYMBOL = {"elastic_net": "circle", "lightgbm": "diamond-open"}

_SPECS = [[{"colspan": 2}, None, {"colspan": 2}, None], [{"colspan": 4}, None, None, None],
          [{}, {}, {}, {}]]  # fmt: skip
_AT = {"a": (1, 1), "d": (1, 3), "c": (2, 1)}
_AT_GATED = {"a": (1, 1), "d": (1, 2)}  # no predictive power: (a) and (d) only


def _b_and_c(fig: Any, data: pd.DataFrame, panel: pd.Series) -> None:
    """(b) one sub-plot per PD feature, the mean and its spread; (c) the SHAP points."""
    b = data[(panel == "b") & data["category"].astype(str).eq("pd_mean")]
    for k, f in enumerate(PD_FEATURES):
        g = b[b["feature"].astype(str) == f]
        x, lo, hi = floats(g["x"]), floats(g["ci_low"]), floats(g["ci_high"])
        fig.add_trace(go.Scatter(x=[*x, *x[::-1]], y=[*lo, *hi[::-1]], fill="toself",
                                 fillcolor=rgba(OKABE_ITO["sky_blue"], 0.5), mode="lines",
                                 line={"width": 0}, showlegend=False, hoverinfo="skip"),
                      3, k + 1)  # fmt: skip
        fig.add_trace(go.Scatter(x=x, y=floats(g["y"]), mode="lines", showlegend=False,
                                 line={"color": OKABE_ITO["blue"]}), 3, k + 1)  # fmt: skip
        ylab = "PD (log2 abundance)" if k == 0 else None
        titles(fig, {"row": 3, "col": k + 1}, PD_XLABEL[f], ylab)
    c = data[panel == "c"]
    hover = [f"{t}{DOT}{f}{DOT}SHAP {fmt(s)}" for t, f, s in zip(c["locus_tag"], c["feature"],
                                                                 c["x"], strict=True)]  # fmt: skip
    dots = {"color": floats(c["ci_low"]), "size": 3, "colorscale": "Viridis"}
    fig.add_trace(points(len(c), x=floats(c["x"]), y=floats(c["y"]), mode="markers", text=hover,
                         showlegend=False, marker=dots), *_AT["c"])  # fmt: skip
    feats = list(dict.fromkeys(c["feature"].astype(str)))
    fig.update_yaxes(tickvals=list(np.arange(len(feats))[::-1]), ticktext=feats, row=2, col=1)
    titles(fig, {"row": 2, "col": 1}, "SHAP (log2 abundance)", None)


@interactive("V15", COLUMNS)
def v15(data: pd.DataFrame, title: str, meta: Meta = None) -> Any:
    """(a) ΔR² with CI bars per model; (b) PD mean and spread; (c) SHAP points; (d) observed vs.
    predicted by category, the locus tag in the hover text. A gated title (no predictive power,
    as V15 writes it) gives (a) and (d) only, under the banner."""
    gated = title.endswith(GATED_SUFFIX)
    if gated:
        at = _AT_GATED
        fig = make_subplots(rows=1, cols=2, subplot_titles=[A_HEADER, D_HEADER])
    else:
        at = _AT
        names = [A_HEADER, D_HEADER, "(c) SHAP (LightGBM, secondary)",
                 *(f"(b) {PD_XLABEL[f]}" for f in PD_FEATURES)]  # fmt: skip
        fig = make_subplots(rows=3, cols=4, specs=_SPECS, subplot_titles=names,
                            vertical_spacing=0.1)  # fmt: skip
    panel = data["panel"].astype(str)
    a = data[panel == "a"]
    blocks = [b for b in BLOCKS if b in set(a["block"].astype(str))]
    for k, m in enumerate(MODELS):
        d = a[a["model"].astype(str) == m]
        x, lo, hi = floats(d["x"]), floats(d["ci_low"]), floats(d["ci_high"])
        y = np.array([blocks.index(str(b)) for b in d["block"]], dtype=float) + (k - 0.5) * 0.3
        hover = [f"{b}{DOT}ΔR² {fmt(v)}{DOT}{ci(lo_, hi_)}"
                 for b, v, lo_, hi_ in zip(d["block"], d["x"], d["ci_low"], d["ci_high"],
                                           strict=True)]  # fmt: skip
        fig.add_trace(points(len(d), x=x, y=y, mode="markers", name=LEGEND[m], text=hover,
                             marker={"color": MODEL_COLOUR[m], "symbol": _SYMBOL[m]},
                             error_x={"type": "data", "symmetric": False, "array": hi - x,
                                      "arrayminus": x - lo}), *at["a"])  # fmt: skip
        na = ~(np.isfinite(lo) & np.isfinite(hi))
        if na.any():
            fig.add_trace(go.Scatter(x=x[na], y=y[na], mode="text", text=["CI unestimated"]
                                     * int(na.sum()), textposition="middle right",
                                     textfont={"color": GREY, "size": 9}, showlegend=False,
                                     hoverinfo="skip"), *at["a"])  # fmt: skip
    fig.update_yaxes(tickvals=list(range(len(blocks))), ticktext=blocks, autorange="reversed",
                     row=1, col=1)  # fmt: skip
    titles(fig, {"row": 1, "col": 1}, "ΔR² (drop-one-block, CV)", None)
    if not gated:
        _b_and_c(fig, data, panel)
    d = data[panel == "d"]
    for cat, colour in CAT_COLOUR.items():
        s = d[d["category"].astype(str) == cat]
        hover = [f"{t}{DOT}predicted {fmt(x)}{DOT}observed {fmt(y)}"
                 for t, x, y in zip(s["locus_tag"], s["x"], s["y"], strict=True)]  # fmt: skip
        fig.add_trace(points(len(s), x=floats(s["x"]), y=floats(s["y"]), mode="markers",
                             name=CAT_LABEL[cat], text=hover, marker={"color": colour, "size": 4}),
                      *at["d"])  # fmt: skip
    xd = abundance(meta, "predicted ") + " (out-of-fold)"
    titles(fig, {"row": at["d"][0], "col": at["d"][1]}, xd, abundance(meta, "observed "))
    notes = [DELTA_FOOT] if len(a) else []
    fig.add_annotation(text="<br>".join([*notes, CAVEATS[0]]), x=0.0, y=-0.05, xref="paper",
                       yref="paper", xanchor="left", yanchor="top", align="left",
                       showarrow=False, font={"color": GREY, "size": 9})  # fmt: skip
    if gated:
        fig.add_annotation(text=BANNER, x=0.5, y=1.12, xref="paper", yref="paper",
                           showarrow=False, font={"color": GREY, "size": 11})  # fmt: skip
    finish(fig, title, 520 if gated else 1000)
    fig.update_layout(margin={"t": 70, "l": 60, "r": 20, "b": 110})
    return fig
