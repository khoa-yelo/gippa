# ruff: noqa: RUF001  (the report's wording is Unicode: en dash, ≤)
"""M12-04 PR 11c readers: V12 lines and the report cards read the number, and NA shows in plain
words via ``gippa.viz.words.plain`` (AC 11c.4; the V12 footnote for 1 unchanged).

docs/specs/M12-04-outputs-revision.md > PR 11c behaviour 2 and AC 11c.4; § 1 rule 2 (no raw
token in prose).

Interface pinned here (test-author reading):

- "Plain words via ``plain()``": for an NA cell the reader's text holds ``plain(token)`` for
  ``token`` in ``{"unestimated", "<column flag>"}`` (``class_error_unestimated`` for variants,
  ``class_fdr_unestimated`` for novel ORFs). The tests swap ``FLAG_TEXT`` for sentinels, so a
  hard-coded word does not pass.
- V12: ``gippa.viz.variants.title(variants)`` and ``_header(row)`` (the site line) with a float
  NaN ``class_error_rate``; a rate of exactly 1 keeps its ``*`` mark and the footnote trigger.
- Report: the ``alt_calls`` and ``novel_orfs`` cards of ``build_report`` (``_alt_estimate`` and
  ``_novel_fdr``), read through ``RunView`` from a 1.1 tree (NA cells) and from a 1.0 tree (the
  word ``unestimated``, mapped on read).
- Grep (AST) over ``src/gippa``: in a module that names either column, nothing is compared with
  the string ``"unestimated"`` (``==``, ``!=``, ``in``, ``.eq``, ``.isin``, ``.replace``), except
  in ``gippa/io/migration.py`` (the one place that reads 1.0 cells) and regime comparisons of the
  class-FDR side table (an expression naming ``regime``).
"""

from __future__ import annotations

import ast
import math
from collections.abc import Iterator
from pathlib import Path
from types import MappingProxyType

import pandas as pd
import pytest
from helpers.report_content import page, words
from helpers.report_fixtures import run_out
from helpers.viz_out_tree import write_schema_table

import gippa.viz.words as wordsmod
from gippa.viz.variants import DEGENERATE_MARK, _has_one, _header, title

REPO = Path(__file__).resolve().parents[2]
SRC = REPO / "src" / "gippa"
CE, CE_FLAG = "class_error_rate", "class_error_unestimated"
FDR, FDR_FLAG = "class_fdr", "class_fdr_unestimated"
S_WORD, S_CE, S_FDR = "zzplainword", "zzplaince", "zzplainfdr"


@pytest.fixture
def sentinels(monkeypatch: pytest.MonkeyPatch) -> Iterator[None]:
    """``plain`` returns a sentinel for the NA tokens, so the text must come through it."""
    text = dict(wordsmod.FLAG_TEXT)
    text.update({"unestimated": S_WORD, CE_FLAG: S_CE, FDR_FLAG: S_FDR})
    monkeypatch.setattr(wordsmod, "FLAG_TEXT", MappingProxyType(text))
    yield


def _plain_ce(text: str) -> bool:
    return S_WORD in text or S_CE in text


def _plain_fdr(text: str) -> bool:
    return S_WORD in text or S_FDR in text


# ---------------------------------------------------------------- V12


def _site(rate: object, n: object = 7) -> dict[str, object]:
    return {"gene": "g1", "aa_pos": 81, "contig": "chr1", "codon_start": 100, "strand": "+",
            "ref_aa": "Y", "alt_aa": "D", "alt_type": "substitution", "confidence": "high",
            "n_alt": 8, "n_ref": 0, CE: rate, "class_error_n": n}  # fmt: skip


def _alt(rates: list[float]) -> pd.DataFrame:
    n = len(rates)
    return pd.DataFrame({"call": ["ALT"] * n, "confidence": ["high"] * n,
                         CE: pd.Series(rates, dtype=float)})  # fmt: skip


def test_v12_site_line_shows_na_in_plain_words(sentinels: None) -> None:
    for na in (math.nan, None, pd.NA):
        line = _header(_site(na, None))
        assert _plain_ce(line), line
        assert "class error NA" not in line and "nan" not in line.lower(), line
        assert "unestimated" not in line, line  # the raw word is not hard-coded


def test_v12_title_shows_all_na_in_plain_words(sentinels: None) -> None:
    t = title(_alt([math.nan, math.nan]))
    assert _plain_ce(t), t
    assert "nan" not in t.lower() and "unestimated" not in t, t


def test_v12_title_range_ignores_na_rows() -> None:
    t = title(_alt([0.2, math.nan, 0.4]))
    assert "0.2" in t and "0.4" in t and "nan" not in t.lower(), t


def test_v12_rate_one_keeps_the_mark_and_footnote_trigger() -> None:
    """Edge case: exactly 1 stays a number; the V12 footnote logic is unchanged."""
    assert _has_one(_alt([0.2, 1.0]))
    from gippa.viz.style import fmt_num

    # PR 5 (P-18): 1 prints through fmt_num ("1.0")
    assert f"{fmt_num(1.0)}{DEGENERATE_MARK}" in title(_alt([0.2, 1.0]))
    assert not _has_one(_alt([math.nan, 0.5]))
    assert DEGENERATE_MARK in _header(_site(1.0))


def test_v12_zero_is_a_number() -> None:
    assert "class error 0" in _header(_site(0.0))
    assert "class error 0" in title(_alt([0.0]))


# ---------------------------------------------------------------- report cards


def _variants(out: Path, rates: list[object], flags: list[object]) -> None:
    rows = [{"variant_id": f"VAR_{i + 1}", "contig": "chr1", "codon_start": 100 + 3 * i,
             "strand": "+", "gene": "g00001", "aa_pos": 10 + i, "ref_codon": "GCT",
             "ref_aa": "A", "alt_aa": "V", "alt_type": "substitution", "call": "ALT",
             "confidence": "high", "n_ref": 0, "n_alt": 3, "nt_changes": 1,
             "alt_codon_vcf": "GTT", CE: r, "class_error_n": None if r is None else 50,
             "source": "denovo", "flags": f}
            for i, (r, f) in enumerate(zip(rates, flags, strict=True))]  # fmt: skip
    write_schema_table(out / "tables" / "variants.tsv", "variants", rows)


def _novel(out: Path, rows: list[tuple[str, bool, object, object]]) -> None:
    body = [{"orf_id": i, "class": "intergenic", "called": k, FDR: f, "flags": fl}
            for i, k, f, fl in rows]  # fmt: skip
    write_schema_table(out / "tables" / "novel_orfs.tsv", "novel_orfs", body)


def _estimate(out: Path, card: str) -> str:
    found = page(out).cards()[card].by_class("estimate")
    assert found, f"card {card}: no .estimate"
    return words(found[0].text())


def test_alt_card_all_na_in_plain_words(tmp_path: Path, sentinels: None) -> None:
    out = run_out(tmp_path)
    _variants(out, [None, None], [CE_FLAG, f"snv;{CE_FLAG}"])
    est = _estimate(out, "alt_calls")
    assert _plain_ce(est), est


def test_alt_card_partial_na_counts_them_in_plain_words(tmp_path: Path, sentinels: None) -> None:
    out = run_out(tmp_path)
    _variants(out, [0.012, None, 0.034], ["snv", CE_FLAG, None])
    est = _estimate(out, "alt_calls")
    assert "0.012–0.034" in est, est
    assert "1" in est and _plain_ce(est), est


def test_alt_card_reads_a_1_0_tree_with_the_word(tmp_path: Path, sentinels: None) -> None:
    out = run_out(tmp_path)
    _variants(out, ["unestimated"], ["snv"])
    est = _estimate(out, "alt_calls")
    assert _plain_ce(est), est


def test_novel_card_called_na_in_plain_words(tmp_path: Path, sentinels: None) -> None:
    out = run_out(tmp_path)
    _novel(out, [("N1", True, None, FDR_FLAG), ("N2", True, 0.02, None)])
    est = _estimate(out, "novel_orfs")
    assert _plain_fdr(est), est
    assert "0.02" not in est, est  # one unestimated called row: no bound (M9-16 B1)


def test_novel_card_numeric_zero_is_a_bound(tmp_path: Path) -> None:
    out = run_out(tmp_path)
    _novel(out, [("N1", True, 0.0, None), ("N2", False, None, FDR_FLAG)])
    est = _estimate(out, "novel_orfs")
    assert "class FDR" in est and "unestimated" not in est, est


def test_novel_card_reads_a_1_0_tree_with_the_word(tmp_path: Path, sentinels: None) -> None:
    out = run_out(tmp_path)
    _novel(out, [("N1", True, "unestimated", FDR_FLAG)])
    est = _estimate(out, "novel_orfs")
    assert _plain_fdr(est), est


# ---------------------------------------------------------------- grep: no string comparison

_COLUMNS = frozenset({CE, FDR})
_ALLOWED = {"io/migration.py"}
_METHODS = frozenset({"eq", "ne", "isin", "replace", "contains", "startswith", "endswith"})


def _unestimated_names(tree: ast.AST) -> set[str]:
    """Module names bound to the string ``"unestimated"`` (assigned or imported)."""
    names = {t.id for n in ast.walk(tree) if isinstance(n, ast.Assign)
             and isinstance(n.value, ast.Constant) and n.value.value == "unestimated"
             for t in n.targets if isinstance(t, ast.Name)}  # fmt: skip
    for n in ast.walk(tree):
        if isinstance(n, ast.ImportFrom):
            names |= {a.asname or a.name for a in n.names if a.name == "UNESTIMATED"}
    return names


def _is_word(node: ast.AST, names: set[str], seq: bool = False) -> bool:
    if isinstance(node, ast.Constant):
        return node.value == "unestimated"
    if isinstance(node, ast.Name):
        return node.id in names
    if isinstance(node, ast.Attribute):
        return node.attr == "UNESTIMATED"
    if seq and isinstance(node, ast.List | ast.Tuple | ast.Set):  # ``.isin([...])`` arguments
        return any(_is_word(e, names) for e in node.elts)
    return False


def _names_a_column(tree: ast.AST) -> bool:
    return any((isinstance(n, ast.Constant) and n.value in _COLUMNS)
               or (isinstance(n, ast.Attribute) and n.attr in _COLUMNS)
               for n in ast.walk(tree))  # fmt: skip


def _offences(path: Path) -> list[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    if not _names_a_column(tree):
        return []
    names = _unestimated_names(tree)
    bad = []
    for n in ast.walk(tree):
        if isinstance(n, ast.Compare):
            hit = any(_is_word(o, names) for o in (n.left, *n.comparators))
        elif isinstance(n, ast.Call) and isinstance(n.func, ast.Attribute):
            hit = n.func.attr in _METHODS and any(_is_word(a, names, seq=True) for a in n.args)
        else:
            continue
        src = ast.unparse(n)
        if hit and "regime" not in src:
            bad.append(f"{path.relative_to(SRC)}:{n.lineno}: {src}")
    return bad


def test_no_reader_compares_the_columns_with_the_word_unestimated() -> None:
    bad = [o for p in sorted(SRC.rglob("*.py"))
           if str(p.relative_to(SRC)) not in _ALLOWED for o in _offences(p)]  # fmt: skip
    assert not bad, "compare the number (NA), not the word:\n" + "\n".join(bad)
