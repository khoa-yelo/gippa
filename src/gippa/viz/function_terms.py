"""Function-term coverage of a run, read from its ``function.collections`` events (M12-04 PR 8,
F-01).

The ``function`` step logs a ``function.collections`` event: collection ->
``{"n_terms", "coverage", "n_unknown_ids"}``. A collection is well-formed when it is a mapping
holding a numeric ``coverage`` or ``n_terms``; the run's coverage is the last event with at least
one well-formed collection (its malformed collections are dropped). When no event qualifies (an
old run, ``function`` off, malformed events) coverage is unknown and the report and plots behave
as before. Pure apart from the event read a ``RunView`` does.
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
    "well_formed",
]

EVENT = "function.collections"
#: the skip reason of V08 when no collection has a term (``SKIP_TEXT``)
NO_TERMS = "no_terms:function"
GO = ("go_bp", "go_mf", "go_cc", "go_slim")
NOT_INSTALLED = "not installed"

Collections = Mapping[str, Mapping[str, Any]]


def collections(view: RunView) -> Collections | None:
    """The well-formed collections of the last ``function.collections`` event that has at least
    one, or ``None`` when no event qualifies (or the events cannot be read)."""
    try:
        found = view.events(EVENT)
    except Exception:
        return None
    for data in reversed(found):
        good = {str(k): v for k, v in data.items() if well_formed(v)}
        if good:
            return good
    return None


def _coverage(entry: Mapping[str, Any]) -> float | None:
    cov = entry.get("coverage")
    if isinstance(cov, bool) or not isinstance(cov, int | float) or not math.isfinite(cov):
        return None
    return float(cov)


def _n_terms(entry: Mapping[str, Any]) -> int | None:
    n = entry.get("n_terms")
    return int(n) if isinstance(n, int | float) and not isinstance(n, bool) else None


def well_formed(entry: object) -> bool:
    """Whether ``entry`` is a mapping holding a numeric ``coverage`` or ``n_terms``."""
    return isinstance(entry, Mapping) and (
        _coverage(entry) is not None or _n_terms(entry) is not None
    )


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
    """The collections whose keys the event holds, as named in the no-terms sentences: COG
    (``cog``), GO (any ``go_*``), KEGG (``kegg``), joined "A or B" / "A, B or C"."""
    found = [name for name, hit in (("COG", "cog" in colls),
                                    ("GO", any(str(k).startswith("go_") for k in colls)),
                                    ("KEGG", "kegg" in colls)) if hit]  # fmt: skip
    if not found:
        return "function"
    return found[0] if len(found) == 1 else f"{', '.join(found[:-1])} or {found[-1]}"


def no_terms_sentence(colls: Collections) -> str:
    """The first sentence of the Function notice."""
    return f"No {names(colls)} terms were assigned to any gene."


_PCT_D = 10_000


def _pct(cov: float | None) -> str:
    """A whole percent; a share strictly between 0 and 1 never reads 0 % or 100 %."""
    if cov is None:
        return "NA"
    n = round(cov * _PCT_D)
    if 0 < cov < 1:
        n = min(max(n, 1), _PCT_D - 1)
    return fmt_pct(n, _PCT_D)


def coverage_note(colls: Collections) -> str:
    """Coverage per collection as whole percents: ``"COG 72 %, GO 41 %, KEGG 66 %"``; "KEGG not
    installed" when the event has no ``kegg`` key. GO is shown as its best-covered collection:
    the event has no union over the GO collections, so this is a lower bound on the share of
    genes with any GO term."""
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
