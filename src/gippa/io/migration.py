"""Reading 1.0 output trees: the 1.0 -> 1.1 column renames (ADR-0075; M12-04 PR 11b).

GIPPA 1.1 renamed ten columns in five tables (one name per quantity). No writer emits a 1.0
name and no 1.0 name is kept as an alias; this module only lets ``RunView`` (and so ``gippa
report`` and ``gippa compare``) read a tree written by 1.0. It is the one place in ``src/`` that
spells the 1.0 names.

``COLUMN_RENAMES`` maps a table (its schema stem) to ``{1.0 name: 1.1 name}``. ``VALUE_RENAMES``
maps the 1.0 feature names held as *values* (``determinants_model.feature_or_block``, the
``feature`` column of the ``determinants_pd`` and ``determinants_shap`` support files) to their
1.1 names, the same renames as the ``determinants_features`` columns. Pure functions only.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import pandas as pd

__all__ = [
    "COLUMN_RENAMES",
    "VALUE_COLUMNS",
    "VALUE_RENAMES",
    "migrate_columns",
    "migrate_values",
]

#: table (schema file stem) -> {1.0 column name: 1.1 column name}
COLUMN_RENAMES: Mapping[str, Mapping[str, str]] = {
    "novel_orfs": {"ostir_tir": "tir_ostir"},
    "determinants_features": {
        "mw": "mw_da",
        "promotech_max": "promoter_score",
        "promotech_tu": "promoter_score_tu",
    },
    "peptide_hits": {"peptide_len": "peptide_len_aa"},
    "pairs": {
        "class": "pair_class",
        "promoter_between": "promoter_score_between",
        "terminator_between": "terminator_score_between",
    },
    "operon_flags": {"uniop_class": "pair_class", "flag": "flags"},
}

#: 1.0 feature name -> 1.1 feature name, for the columns of ``VALUE_COLUMNS``
VALUE_RENAMES: Mapping[str, str] = dict(COLUMN_RENAMES["determinants_features"])

#: table or support file stem -> the column holding feature names
VALUE_COLUMNS: Mapping[str, str] = {
    "determinants_model": "feature_or_block",
    "determinants_pd": "feature",
    "determinants_shap": "feature",
}


def migrate_columns(table: str, columns: Sequence[str]) -> list[str]:
    """``columns`` with each 1.0 name of ``table`` replaced by its 1.1 name, in place; other
    columns and other tables unchanged. A header holding both a 1.0 name and its 1.1 name is
    ambiguous and raises ``ValueError``."""
    renames = COLUMN_RENAMES.get(table, {})
    present = set(columns)
    for old, new in renames.items():
        if old in present and new in present:
            raise ValueError(
                f"table {table}: the header has both the 1.0 column {old!r} and its 1.1 name "
                f"{new!r}; cannot tell which to read"
            )
    return [renames.get(c, c) for c in columns]


def migrate_values(table: str, frame: pd.DataFrame) -> pd.DataFrame:
    """``frame`` with the 1.0 feature names in its ``VALUE_COLUMNS[table]`` column replaced by
    their 1.1 names (a copy when anything changes); other tables are returned unchanged."""
    col = VALUE_COLUMNS.get(table)
    if col is None or col not in frame.columns:
        return frame
    values = frame[col]
    hit = values.isin(list(VALUE_RENAMES))
    if not bool(hit.any()):
        return frame
    out = frame.copy()
    out.loc[hit, col] = values[hit].map(dict(VALUE_RENAMES))
    return out
