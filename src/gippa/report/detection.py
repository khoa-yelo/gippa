"""The Gene detection views of ``tables/detection_summary.tsv`` (M12-04 PR 13, T-07).

Pure: the long file as read by ``RunView.table`` in, header and rows of two pivoted tables out.
The long file is unchanged and no table file is written (Q8).

- **Tier view:** one row per ``q_threshold`` group of the ``tier`` section, ascending (a group
  that is no finite number after the finite ones), with ``n_detected``, ``n_genes``, ``fdr``,
  ``fdr_2pep`` and ``fp``; a column NA on every row is dropped (``decoy_picked``:
  ``q | detected | FDR``). A tier without such rows shows its counts only (``group_type = tier``
  rows): ``protein_list`` its detected count, ``two_peptide`` also the 2-peptide count.
- **Detectability view:** observed over expected by group (the ``o_over_e`` rows of the
  ``basic`` section, ``gippa.context.summary.basic_rows``), in file order.

Counts are exact (``fmt_num(int(v))``); every other value is ``fmt_num(v)``; NA reads "NA".
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import Any

import pandas as pd

from gippa.viz.style import fmt_num

#: a pivoted table: header cells and body rows
Grid = tuple[list[str], list[list[str]]]

TIER_METRICS = ("n_detected", "n_genes", "fdr", "fdr_2pep", "fp")
TIER_HEADER = ["q", "detected", "genes", "FDR", "FDR (2-pep)", "FP"]
COUNT_METRICS = ("n_detected", "n_detected_2pep")
COUNT_HEADER = ["tier", "detected", "detected with ≥ 2 peptides"]
#: the count-only tiers that report a 2-peptide count
TWO_PEP_TIERS = frozenset({"two_peptide"})
OE_METRICS = ("n_genes", "n_detected", "expected", "o_over_e")
OE_HEADER = ["grouping", "group", "genes", "detected (O)", "expected (E)", "O/E"]
_COUNTS = frozenset({"n_detected", "n_detected_2pep", "n_genes"})
#: ``group_type`` -> its label in the detectability view
_GROUPINGS = {"replicon": "Replicon", "feature_class": "Feature class",
              "cog_category": "COG category"}  # fmt: skip


def cell(metric: str, value: Any) -> str:
    """``value`` as shown: counts exact, other numbers to 2 significant digits, NA as "NA"."""
    v = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
    if v is None or pd.isna(v) or not math.isfinite(float(v)):
        return "NA"
    v = float(v)
    return fmt_num(int(v)) if metric in _COUNTS and v.is_integer() else fmt_num(v)


def _col(frame: pd.DataFrame, name: str) -> pd.Series:
    return frame[name].astype(str)


def _pivot(rows: pd.DataFrame, metrics: Sequence[str]) -> list[str]:
    """One group's ``metrics`` cells (a later duplicate row wins; a missing metric is NA)."""
    by = dict(zip(_col(rows, "metric"), rows["value"], strict=True))
    return [cell(m, by.get(m)) for m in metrics]


def _float_key(group: str) -> tuple[bool, float, str]:
    """Finite thresholds ascending, then any other group by its text (review r1 L5)."""
    try:
        q = float(group)
    except ValueError:
        q = math.nan
    finite = math.isfinite(q)
    return not finite, q if finite else 0.0, group


def _drop_na(header: list[str], rows: list[list[str]]) -> Grid:
    """``header`` and ``rows`` without the columns (but the first) NA on every row."""
    keep = [0] + [k for k in range(1, len(header)) if any(r[k] != "NA" for r in rows)]
    return [header[k] for k in keep], [[r[k] for k in keep] for r in rows]


def tier_view(frame: pd.DataFrame) -> Grid | None:
    """The tier view; ``None`` when the ``tier`` section has no row."""
    tier = frame[_col(frame, "section") == "tier"]
    curve = tier[_col(tier, "group_type") == "q_threshold"]
    if len(curve):
        groups = sorted(set(_col(curve, "group")), key=_float_key)
        rows = [[g, *_pivot(curve[_col(curve, "group") == g], TIER_METRICS)] for g in groups]
        return _drop_na(TIER_HEADER, rows)
    counts = tier[_col(tier, "group_type") == "tier"]
    if not len(counts):
        return None
    groups = list(dict.fromkeys(_col(counts, "group")))
    width = 3 if TWO_PEP_TIERS & set(groups) else 2
    rows = [[g, *_pivot(counts[_col(counts, "group") == g], COUNT_METRICS)] for g in groups]
    return COUNT_HEADER[:width], [r[:width] for r in rows]


def oe_view(frame: pd.DataFrame) -> Grid | None:
    """The detectability view; ``None`` when no ``o_over_e`` row exists."""
    basic = frame[_col(frame, "section") == "basic"]
    oe = basic[_col(basic, "metric") == "o_over_e"]
    if not len(oe):
        return None
    keys = list(dict.fromkeys(zip(_col(oe, "group_type"), _col(oe, "group"), strict=True)))
    rows = []
    for gtype, group in keys:
        sub = basic[(_col(basic, "group_type") == gtype) & (_col(basic, "group") == group)]
        label = _GROUPINGS.get(gtype, gtype.replace("_", " ").capitalize())
        rows.append([label, group, *_pivot(sub, OE_METRICS)])
    return OE_HEADER, rows


def views(frame: pd.DataFrame) -> dict[str, Grid]:
    """``{"tier": ..., "oe": ...}``, each view present only when it has a row."""
    found = {"tier": tier_view(frame), "oe": oe_view(frame)}
    return {k: v for k, v in found.items() if v is not None}
