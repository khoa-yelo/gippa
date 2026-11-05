"""M12-04 PR 11b: the 1.0 -> 1.1 column renames (docs/specs/M12-04-outputs-revision.md § 5.1).

The single test-side copy of the rename table, shared by the ``test_m12_04_pr11b_*`` modules.
Lead defaults (binding): Q1 ``pair_class`` in both tables, Q2 ``operon_flags.flag`` -> ``flags``,
Q15 ``terminator_between`` -> ``terminator_score_between`` and ``promotech_*`` ->
``promoter_score*`` included.
"""

from __future__ import annotations

from collections.abc import Mapping

#: table (schema file stem) -> {1.0 name: 1.1 name}
RENAMES: Mapping[str, Mapping[str, str]] = {
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

#: the published path of each renamed table under ``OUT/``, and its ``RunView`` key
PATHS: Mapping[str, tuple[str, str]] = {
    "novel_orfs": ("tables/novel_orfs.tsv", "novel_orfs"),
    "determinants_features": ("tables/determinants_features.tsv", "determinants_features"),
    "peptide_hits": ("mapping/peptide_hits.tsv", "peptide_hits"),
    "pairs": ("tables/pairs.tsv", "pairs"),
    "operon_flags": ("tables/operon_flags.tsv", "operon_flags"),
}

#: feature values renamed the same way (determinants_model.feature_or_block, the
#: determinants_pd / determinants_shap ``feature`` column)
VALUE_RENAMES: Mapping[str, str] = {
    "promotech_max": "promoter_score",
    "promotech_tu": "promoter_score_tu",
    "mw": "mw_da",
}

#: every (table, old, new) triple
TRIPLES: tuple[tuple[str, str, str], ...] = tuple(
    (t, old, new) for t, m in RENAMES.items() for old, new in m.items()
)
#: the 1.0 names that are distinctive enough for a source grep (``class``, ``flag`` and ``mw``
#: are generic words; ``peptide_len`` is also the config key ``input.peptide_len``)
GREP_OLD = (
    "ostir_tir",
    "promotech_max",
    "promotech_tu",
    "uniop_class",
    "peptide_len",
    "promoter_between",
    "terminator_between",
)


def to_1_0(table: str, columns: list[str]) -> list[str]:
    """``columns`` (either version) with the 1.1 names of ``table`` put back to 1.0."""
    back = {new: old for old, new in RENAMES.get(table, {}).items()}
    return [back.get(c, c) for c in columns]
