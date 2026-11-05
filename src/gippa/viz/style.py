"""The one style module for every GIPPA plot (GIPPA_SPEC.md §10.1; M9-01, M9-index decisions 6-10).

Pure constants and functions, no file I/O. ``apply()`` sets the rcParams once per process (the
``plots`` step and every test fixture call it); ``new_figure`` makes a figure at one of the two
§10.1 widths. Figures are created without pyplot so no global figure manager holds them.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, MutableMapping
from decimal import Decimal
from types import MappingProxyType
from typing import Any, Literal, cast

import matplotlib
import numpy as np
import pandas as pd
from cycler import cycler
from matplotlib.figure import Figure

SINGLE_IN: float = 3.5
DOUBLE_IN: float = 7.2

OKABE_ITO: Mapping[str, str] = MappingProxyType(
    {
        "black": "#000000",
        "orange": "#E69F00",
        "sky_blue": "#56B4E9",
        "bluish_green": "#009E73",
        "yellow": "#F0E442",
        "blue": "#0072B2",
        "vermillion": "#D55E00",
        "reddish_purple": "#CC79A7",
        "grey": "#999999",
    }
)
DETECTION: Mapping[str, str] = MappingProxyType(
    {"detected": "#0072B2", "not_detected": "#BBBBBB", "unestimated": "#E69F00"}
)
#: ``genes.tsv`` ``detection_status`` -> colour (decision 9; ``partial_excluded`` is hatched)
STATUS_COLOUR: Mapping[str, str] = MappingProxyType(
    {
        "detected": DETECTION["detected"],
        "not_detected": DETECTION["not_detected"],
        "partial_excluded": DETECTION["not_detected"],
        "single_peptide_unestimated": DETECTION["unestimated"],
    }
)
REGION: Mapping[str, str] = MappingProxyType(
    {
        "prophage": "#D55E00",
        "plasmid": "#009E73",
        "bgc": "#CC79A7",
        "crispr": "#56B4E9",
        "amr": "#F0E442",
    }
)
#: island direction -> (colour, hatch) (decision 10)
ISLAND: Mapping[str, tuple[str, str]] = MappingProxyType(
    {"high": ("#0072B2", "///"), "low": ("#D55E00", "\\\\\\")}
)

_GROUPS = {
    "information": "JAKLB",
    "cellular": "DYVTMNZWUOX",
    "metabolism": "CGEFHIPQ",
    "poorly": "RS",
}
#: COG letter -> super-category (the four NCBI groups; ``data/cog_categories.tsv``, decision 8)
COG_SUPER: Mapping[str, str] = MappingProxyType(
    {letter: group for group, letters in _GROUPS.items() for letter in letters}
)
COG_SUPER_COLOUR: Mapping[str, str] = MappingProxyType(
    {"information": "#0072B2", "cellular": "#009E73", "metabolism": "#E69F00", "poorly": "#999999"}
)
#: neutral greys drawn as chart furniture (gene-arrow edges, label envelopes, ring backgrounds)
NEUTRAL_GREY: Mapping[str, str] = MappingProxyType(
    {"edge": "#333333", "envelope": "#CCCCCC", "sector": "#DDDDDD", "ring": "#F4F4F4"}
)
CONTINUOUS: str = "viridis"
LOCUS_ABUNDANCE: str = "cividis"

#: the category cycle: Okabe-Ito order without black and grey
_CYCLE = [v for k, v in OKABE_ITO.items() if k not in ("black", "grey")]

RC: Mapping[str, Any] = MappingProxyType(
    {
        "font.family": ["DejaVu Sans"],
        "font.size": 7,
        "axes.titlesize": 8,
        "axes.labelsize": 7,
        "xtick.labelsize": 6,
        "ytick.labelsize": 6,
        "legend.fontsize": 6,
        "figure.titlesize": 8,
        "svg.fonttype": "none",
        "svg.hashsalt": "gippa",  # ids are renumbered on write (ADR-0064 item 7, amended)
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "path.simplify": True,
        # pinned: importing pyCirclize sets a tight bbox and 0.5-in padding (M9-08 decision 16)
        "savefig.bbox": None,
        "savefig.pad_inches": 0.1,
    }
)


def apply() -> None:
    """Select the Agg backend and set the GIPPA rcParams (behaviour 1)."""
    matplotlib.use("Agg", force=True)
    # RcParams is typed with literal keys; RC is checked by the style tests instead
    cast("MutableMapping[str, Any]", matplotlib.rcParams).update(RC)
    matplotlib.rcParams["axes.prop_cycle"] = cycler(color=_CYCLE)


def new_figure(width: Literal["single", "double"], height_in: float, **kw: Any) -> Figure:
    """A figure ``SINGLE_IN`` or ``DOUBLE_IN`` wide and ``height_in`` tall."""
    if width not in ("single", "double"):
        raise ValueError(f"width must be 'single' or 'double', not {width!r}")
    w = SINGLE_IN if width == "single" else DOUBLE_IN
    return Figure(figsize=(w, height_in), **kw)


def cog_colour(letter: str) -> str:
    """The super-category colour of one COG letter; unknown letters are poorly characterised."""
    return COG_SUPER_COLOUR[COG_SUPER.get(letter, "poorly")]


def plotly_template() -> Any:
    """The §10.1 look for the Plotly twins: DejaVu Sans (sans-serif in the browser), the rcParams
    font sizes, white background, the category cycle as the colorway (M9-15 behaviour 1)."""
    import plotly.graph_objects as go

    axis = {"tickfont": {"size": RC["xtick.labelsize"]}, "linecolor": "#000000",
            "gridcolor": "#EEEEEE", "zeroline": False, "showline": True, "ticks": "outside",
            "title": {"font": {"size": RC["axes.labelsize"]}}}  # fmt: skip
    font = {"family": "DejaVu Sans, sans-serif", "size": RC["font.size"], "color": "#000000"}
    return go.layout.Template(
        layout={
            "font": font,
            "title": {"font": {"size": RC["axes.titlesize"]}, "x": 0.01},
            "xaxis": axis,
            "yaxis": axis,
            "legend": {"font": {"size": RC["legend.fontsize"]}},
            "paper_bgcolor": "#FFFFFF",
            "plot_bgcolor": "#FFFFFF",
            "colorway": list(_CYCLE),
            "hoverlabel": {"font": {"size": RC["font.size"] + 2}},
        }
    )


# --- numbers and labels in figure text and the report (M12-04 PR 0) ----------------------------

MINUS = "\u2212"
_ABUNDANCE_NAME: Mapping[str, str] = MappingProxyType(
    {"ibaq": "iBAQ", "maxlfq": "MaxLFQ", "top3": "Top3", "intensity": "intensity"}
)


def _is_na(x: object) -> bool:
    """``None``, ``pd.NA`` and NaN (float or numpy) are NA."""
    if x is None or x is pd.NA:
        return True
    try:
        return bool(math.isnan(cast(float, x)))
    except (TypeError, ValueError):
        return False


def _signed(text: str, negative: bool) -> str:
    return f"{MINUS}{text}" if negative else text


def _sig_text(x: float, sig: int) -> str:
    """``|x|`` to ``sig`` significant digits, trailing zeros kept: plain inside
    ``1e-4 <= |x| < 1e6`` after rounding, ``1.2e-12`` style outside it."""
    rounded = f"{abs(x):.{sig - 1}e}"
    mant, exp = rounded.split("e")
    if -4 <= int(exp) < 6:
        return format(Decimal(rounded), "f")
    return f"{mant}e{int(exp):+d}" if int(exp) > 0 else f"{mant}e{int(exp)}"


def fmt_num(x: object, sig: int = 2) -> str:
    """A number for prose: integers exact (comma-grouped from 10,000), other values to ``sig``
    significant digits with trailing zeros kept; a true minus (U+2212), never on zero; NA and
    infinity -> ``"NA"``."""
    if _is_na(x):
        return "NA"
    if isinstance(x, (int, np.integer)):
        n = int(x)
        body = f"{abs(n):,}" if abs(n) >= 10_000 else str(abs(n))
        return _signed(body, n < 0)
    v = float(cast(float, x))
    if not math.isfinite(v):
        return "NA"
    if v == 0.0:
        return _sig_text(0.0, sig)
    return _signed(_sig_text(v, sig), v < 0)


def fmt_q(q: object) -> str:
    """A q- or p-value: ``"< 0.001"`` below 0.001, otherwise 2 significant digits with trailing
    zeros kept (``0.05`` -> ``"0.050"``); NA -> ``"NA"``."""
    if _is_na(q):
        return "NA"
    v = float(cast(float, q))
    if not math.isfinite(v):
        return "NA"
    return "< 0.001" if v < 0.001 else _sig_text(v, 2)


def fmt_pct(n: int, d: int) -> str:
    """``n / d`` as a whole percent (``"12 %"``); ``d == 0`` -> ``"NA"``; a partial share never
    rounds to 0 % or 100 %."""
    if d == 0:
        return "NA"
    pct = int(100 * n / d + 0.5)
    if 0 < n < d:
        pct = min(max(pct, 1), 99)
    return f"{pct} %"


def abundance_label(atype: object, scale: object) -> str | None:
    """The abundance axis label from ``abundance_type`` and ``abundance_scale``; ``None`` when
    the run has no abundance (``none`` or NA), so the caller draws no abundance axis."""
    if _is_na(atype) or atype == "none":
        return None
    name = _ABUNDANCE_NAME.get(str(atype))
    if name is None or _is_na(scale):
        return "log2 abundance (measure unknown)" if name is None else f"log2 {name}"
    return f"log2 {name} ({scale})"
