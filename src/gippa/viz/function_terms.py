"""Function-term coverage of a run, read from its ``function.collections`` event (M12-04 PR 8,
F-01).

The ``function`` step logs one ``function.collections`` event: collection ->
``{"n_terms", "coverage", "n_unknown_ids"}``. The last event with at least one collection is the
run's coverage; a run without one (an old run, ``function`` off) has none, and its report and
plots behave as before. Pure apart from the event read a ``RunView`` does.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Any

from gippa.io.run_view import RunView
from gippa.viz.style import fmt_pct

__all__ = [
    "NO_TERMS",
    "collections",
    "coverage_note",
    "has_terms",
    "kegg_terms",
    "names",
    "no_terms",
    "no_terms_sentence",
]

EVENT = "function.collections"
#: the skip reason of V08 when no collection has a term (``SKIP_TEXT``)
NO_TERMS = "no_terms:function"
GO = ("go_bp", "go_mf", "go_cc", "go_slim")
NOT_INSTALLED = "not installed"

Collections = Mapping[str, Mapping[str, Any]]


def collections(view: RunView) -> Collections | None:
    """The last ``function.collections`` event's collections, or ``None`` when the run has no
    such event with at least one collection (or its events cannot be read)."""
    try:
        found = view.events(EVENT)
    except Exception:
        return None
    if not found:
        return None
    data = {str(k): v for k, v in found[-1].items() if isinstance(v, Mapping)}
    return data or None


def _coverage(entry: Mapping[str, Any]) -> float | None:
    cov = entry.get("coverage")
    if isinstance(cov, bool) or not isinstance(cov, int | float) or not math.isfinite(cov):
        return None
    return float(cov)


def _n_terms(entry: Mapping[str, Any]) -> int | None:
    n = entry.get("n_terms")
    return int(n) if isinstance(n, int | float) and not isinstance(n, bool) else None


def has_terms(entry: Mapping[str, Any]) -> bool:
    """Whether a collection has a term: neither its coverage nor its ``n_terms`` is 0 (a null
    coverage with 0 terms has none)."""
    return _coverage(entry) != 0.0 and _n_terms(entry) != 0


def no_terms(colls: Collections | None) -> bool:
    """Whether an event with at least one collection exists and no collection has a term."""
    return bool(colls) and not any(has_terms(e) for e in (colls or {}).values())


def kegg_terms(colls: Collections) -> bool:
    """Whether the event's ``kegg`` collection has at least one term."""
    n = _n_terms(colls.get("kegg") or {})
    return n is not None and n > 0


def names(colls: Collections) -> str:
    """The collections named in the no-terms sentences: COG, GO and KEGG, or COG and GO when the
    event has no ``kegg`` key."""
    return "COG, GO or KEGG" if "kegg" in colls else "COG or GO"


def no_terms_sentence(colls: Collections) -> str:
    """The first sentence of the Function notice."""
    return f"No {names(colls)} terms were assigned to any gene."


def _pct(cov: float | None) -> str:
    return "NA" if cov is None else fmt_pct(round(cov * 10_000), 10_000)


def coverage_note(colls: Collections) -> str:
    """Coverage per collection as whole percents: ``"COG 72 %, GO 41 %, KEGG 66 %"``; GO is its
    best-covered collection; "KEGG not installed" when the event has no ``kegg`` key."""
    parts: list[str] = []
    if "cog" in colls:
        parts.append(f"COG {_pct(_coverage(colls['cog']))}")
    go = [c for g in GO if g in colls for c in (_coverage(colls[g]),) if c is not None]
    if any(g in colls for g in GO):
        parts.append(f"GO {_pct(max(go) if go else None)}")
    parts.append(f"KEGG {_pct(_coverage(colls['kegg']))}" if "kegg" in colls
                 else f"KEGG {NOT_INSTALLED}")  # fmt: skip
    known = {"cog", "kegg", *GO}
    parts += [f"{c} {_pct(_coverage(e))}" for c, e in colls.items() if c not in known]
    return ", ".join(parts)
