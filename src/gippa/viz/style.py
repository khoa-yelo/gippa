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
    {
        "edge": "#333333",
        "envelope": "#CCCCCC",
        "sector": "#DDDDDD",
        "ring": "#F4F4F4",
        "no_category": "#E0E0E0",  # V08 genes without a COG letter, hatched (M12-04 P-19)
    }
)
CONTINUOUS: str = "viridis"
LOCUS_ABUNDANCE: str = "cividis"
#: an undetected gene arrow in the V03/V13 locus views (M12-04 P-01): the patch keywords
LOCUS_UNDETECTED: Mapping[str, str] = MappingProxyType(
    {"facecolor": "white", "edgecolor": "#BBBBBB", "hatch": "////"}
)

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
    """``None``, ``pd.NA``, ``pd.NaT`` and NaN (float or numpy) are NA."""
    if x is None or x is pd.NA or x is pd.NaT:
        return True
    try:
        return bool(math.isnan(cast(float, x)))
    except (TypeError, ValueError):
        return False


def _as_number(x: object) -> int | float | None:
    """``x`` as an ``int`` (Python, numpy and bool integers) or a ``float``; ``None`` for NA and
    anything that does not parse as a number."""
    if _is_na(x):
        return None
    if isinstance(x, (bool, np.bool_)):
        return int(bool(x))
    if isinstance(x, (int, np.integer)):
        return int(x)
    try:
        v = float(cast(float, x))
    except (TypeError, ValueError, OverflowError):
        return None
    return None if math.isnan(v) else v


def _signed(text: str, negative: bool) -> str:
    return f"{MINUS}{text}" if negative else text


def _grouped(n: int) -> str:
    """An integer exactly, comma-grouped from 10,000, with a true minus."""
    body = f"{abs(n):,}" if abs(n) >= 10_000 else str(abs(n))
    return _signed(body, n < 0)


def _sig_text(x: float, sig: int) -> str:
    """``|x|`` to ``sig`` significant digits, trailing zeros kept: plain inside
    ``1e-4 <= |x| < 1e6`` after rounding, ``1.2e<U+2212>12`` style outside it (a true minus, no
    ``+`` and no leading zero in the exponent)."""
    rounded = f"{abs(x):.{sig - 1}e}"
    mant, exp_text = rounded.split("e")
    exp = int(exp_text)
    if -4 <= exp < 6:
        return format(Decimal(rounded), "f")
    return f"{mant}e{_signed(str(abs(exp)), exp < 0)}"


def fmt_num(x: object, sig: int = 2) -> str:
    """A number for prose: integers exact (comma-grouped from 10,000), other values to ``sig``
    significant digits with trailing zeros kept; a true minus (U+2212), never on zero; NA,
    infinity and non-numbers -> ``"NA"``. ``sig < 1`` raises ``ValueError``."""
    if sig < 1:
        raise ValueError(f"sig must be >= 1, not {sig}")
    v = _as_number(x)
    if v is None:
        return "NA"
    if isinstance(v, int):
        return _grouped(v)
    if not math.isfinite(v):
        return "NA"
    if v == 0.0:
        return _sig_text(0.0, sig)
    return _signed(_sig_text(v, sig), v < 0)


def fmt_count(x: object) -> str:
    """A count: an integer or integer-valued float printed exactly, comma-grouped from 10,000;
    NA -> ``"NA"``. A non-integer value raises ``ValueError``."""
    v = _as_number(x)
    if v is None:
        if _is_na(x):
            return "NA"
        raise ValueError(f"not a count: {x!r}")
    if isinstance(v, float):
        if not math.isfinite(v):
            return "NA"
        if not v.is_integer():
            raise ValueError(f"not a whole count: {x!r}")
        v = int(v)
    return _grouped(v)


def fmt_alpha(alpha: object) -> str:
    """An FDR level as its shortest exact decimal (``0.05`` -> ``"0.05"``, ``0.025`` ->
    ``"0.025"``), never in scientific notation; NA -> ``"NA"``."""
    v = _as_number(alpha)
    if v is None or not math.isfinite(v):
        return "NA"
    text = format(Decimal(repr(float(v))), "f")
    return _signed(text.lstrip("-"), v < 0)


def fmt_q(q: object) -> str:
    """A q- or p-value: ``"< 0.001"`` below 0.001, otherwise 2 significant digits with trailing
    zeros kept (``0.05`` -> ``"0.050"``); NA and non-numbers -> ``"NA"``."""
    v = _as_number(q)
    if v is None or not math.isfinite(v):
        return "NA"
    return "< 0.001" if v < 0.001 else _sig_text(float(v), 2)


def q_equals(q: object, name: str = "q") -> str:
    """``name`` and ``fmt_q(q)`` joined for prose: ``"q = 0.43"``, ``"q < 0.001"``."""
    text = fmt_q(q)
    return f"{name} {text}" if text.startswith("<") else f"{name} = {text}"


def fmt_pct(n: int, d: int) -> str:
    """``n / d`` as a whole percent (``"12 %"``); ``d == 0`` -> ``"NA"``; a partial share never
    rounds to 0 % or 100 %."""
    if d == 0:
        return "NA"
    pct = int(100 * n / d + 0.5)
    if 0 < n < d:
        pct = min(max(pct, 1), 99)
    return f"{pct} %"


def shorten_middle(text: str, limit: int = 24) -> str:
    """``text`` if it has at most ``limit`` characters, else its head and tail joined by one
    middle ellipsis, ``limit`` characters in all (long contig names in titles and labels)."""
    if len(text) <= limit:
        return text
    head = (limit - 1) // 2
    return f"{text[:head]}\u2026{text[len(text) - (limit - 1 - head) :]}"


def abundance_label(atype: object, scale: object) -> str | None:
    """The abundance axis label from ``abundance_type`` (``ibaq``, ``maxlfq``, ``top3``,
    ``intensity``, any case) and ``abundance_scale`` (``molar`` or ``relative``), e.g.
    ``"log2 iBAQ (molar)"``; an NA or empty scale drops the parentheses. ``None`` when the run
    has no abundance (``none`` or NA), so the caller draws no abundance axis."""
    if _is_na(atype) or str(atype).strip().lower() == "none":
        return None
    name = _ABUNDANCE_NAME.get(str(atype).strip().lower())
    if name is None:
        return "log2 abundance (measure unknown)"
    if _is_na(scale) or not str(scale).strip():
        return f"log2 {name}"
    return f"log2 {name} ({str(scale).strip()})"


#: the plain-abundance axis text when the run's measure is unknown
ABUNDANCE_FALLBACK = "log2 abundance"


def run_abundance_label(genes: pd.DataFrame | None) -> str | None:
    """``abundance_label`` of a run: the first ``abundance_type`` of ``genes.tsv`` that is not
    ``none``/NA, with the scale of the first gene of that type that has one; ``None`` when no
    gene carries a measure (or the columns are absent)."""
    if genes is None or "abundance_type" not in genes:
        return None
    types = genes["abundance_type"]
    keep = types.notna() & (types.astype(str).str.strip().str.lower() != "none")
    if not keep.any():
        return None
    atype = types[keep].iloc[0]
    scale: object = None
    if "abundance_scale" in genes:
        found = genes.loc[keep & (types == atype), "abundance_scale"].dropna()
        found = found[found.astype(str).str.strip() != ""]
        scale = None if found.empty else found.iloc[0]
    return abundance_label(atype, scale)
