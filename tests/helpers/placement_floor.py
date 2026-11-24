"""Costs tables for the ADR-0079 placement floor tests (tests/unit/test_placement_floor_*.py)."""

from __future__ import annotations

import dataclasses
from collections.abc import Iterable, Mapping
from typing import Any

import pandas as pd

from gippa.config import GippaConfig
from gippa.mapping.placement_fdr import PeptideCosts

DATA_KEYS = {"decoy_type", "n_diff", "n_targets", "n_decoys", "q_floor", "alpha"}
ALPHABET19 = "ACDEFGHKLMNPQRSTVWY"
PREFIX = "MDEFGHKLNPQRSTVWYACDEFGHKLMNPQRSTVWY"


def pep_of(i: int, length: int = 9) -> str:
    """Distinct, ordinary peptide ``i`` of ``length`` aa (as in tests/unit/test_placement_fdr)."""
    assert length >= 7 and i < 19**3
    digits = "".join(ALPHABET19[(i // 19**k) % 19] for k in (2, 1, 0))
    return PREFIX[: length - 3] + digits


def costs_frame(rows: Iterable[dict[str, Any]]) -> pd.DataFrame:
    out = []
    for r in rows:
        r = dict(r)
        r.setdefault("length", len(r["peptide"]))
        r.setdefault("source", "database")
        r.setdefault("search_db", "other")
        for k in ("m_target", "m_decoy", "e_target", "e_decoy", "pep", "q_value", "pi_denovo"):
            r.setdefault(k, None)
        out.append(dataclasses.asdict(PeptideCosts(**r)))
    return pd.DataFrame(out, columns=[f.name for f in dataclasses.fields(PeptideCosts)])


class Ids:
    """Hands out distinct peptide numbers across the lists of one costs table."""

    def __init__(self) -> None:
        self.i = 0

    def take(self, n: int) -> list[int]:
        out = list(range(self.i, self.i + n))
        self.i += n
        return out


def stratum(
    ids: Ids, *, m: int, t: int, d: int, source: str = "database", pi: float | None = None
) -> list[dict[str, Any]]:
    """One list: ``t`` target-only peptides (E = 1e-4 .. ) and ``d`` decoy-only peptides with E
    below every target (so every decoy counts at every threshold)."""
    kw: dict[str, Any] = {"source": source}
    if pi is not None:
        kw["pi_denovo"] = pi
    rows = [
        {"peptide": pep_of(i), "m_target": m, "e_target": 1e-4 * (k + 1), **kw}
        for k, i in enumerate(ids.take(t))
    ]
    rows += [
        {"peptide": pep_of(i), "m_decoy": m, "e_decoy": 1e-6 * (k + 1), **kw}
        for k, i in enumerate(ids.take(d))
    ]
    return rows


def cfg(alpha: float = 0.01) -> GippaConfig:
    return GippaConfig.model_validate({"stats": {"alpha": {"placement": alpha}}})


def row9_costs() -> pd.DataFrame:
    """Both decoy types, all three strata: reverse m=0 large (no breach), reverse m=1 T=13
    (breach), reverse m=2 T=0, shuffle m=0 T=150 (no breach), shuffle m=1 T=200 (no breach),
    shuffle m=2 T=68 (breach). Listed out of order on purpose."""
    ids = Ids()
    rows = stratum(ids, m=2, t=68, d=0, source="denovo", pi=0.6)
    rows += stratum(ids, m=1, t=13, d=7)
    rows += stratum(ids, m=1, t=200, d=1, source="denovo", pi=0.6)
    rows += stratum(ids, m=0, t=150, d=0, source="denovo", pi=0.6)
    rows += stratum(ids, m=0, t=400, d=2)
    rows += stratum(ids, m=2, t=0, d=4)
    return costs_frame(rows)


def as_dict(item: Any) -> dict[str, Any]:
    """A floor breach (mapping, dataclass or named tuple) as a plain dict."""
    if dataclasses.is_dataclass(item) and not isinstance(item, type):
        return dataclasses.asdict(item)
    if isinstance(item, Mapping):
        return dict(item)
    if hasattr(item, "_asdict"):
        return dict(item._asdict())
    raise AssertionError(f"a floor breach must be a mapping or a dataclass, got {type(item)}")
