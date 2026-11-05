"""M12-04 PR 9, T-01: every published TSV writer formats floats ``%.6g`` (docs/specs/
M12-04-outputs-revision.md > PR 9; GIPPA_SPEC.md §6.6 "floats ``%.6g``"; ADR-0075).

The review's correction (spec > Corrections, T-01): counted as *more than 6 significant digits*,
only ``tables/detection_summary.tsv`` (``gippa.io.summary_tsv.write_detection_summary``) and
``tables/regions_stats.tsv`` (``gippa.io.regions_tsv.write_regions_stats``) break the rule; both
write the shortest ``repr`` of a float today. PR 9 routes them through ``gippa.io.tsv``.

Pinned here (test author):

- "Through ``gippa.io.tsv`` with ``%.6g``" means the written cell of a ``number`` column equals
  ``gippa.io.tsv.format_value(float(v))`` for every float ``v``: ``%.6g``, ``0`` for zero, ``NA``
  for NaN and the ADR-0006 ``INF``/``-INF`` for infinities (the spec's edge case says "as
  ``gippa.io.tsv`` does today"). An integral float therefore follows ``%.6g`` too
  (``4.000000000000001`` → ``4``; ``1e20`` → ``1e+20``).
- ``integer`` columns are exempt and keep every digit (``end = 1114548``).
- The significant-digit count of a cell strips the sign, the exponent, the decimal point, and
  leading and trailing zeros: ``0.00012345`` has 5 (spec: it passes), ``1234567`` has 7.
- The property over writers covers the writers whose frames can be built from the schema alone:
  ``detection_summary``, ``regions_stats``, ``autocorrelation`` and ``islands`` (the last two
  already comply and guard against a regression). Every other published TSV is covered by the
  committed goldens scan below and by the AC 9.1 integration run
  (tests/integration/test_m12_04_pr9_numbers_run.py).
- The goldens scan reads every ``tests/data/golden/**/*.tsv`` that has a shipped table schema.
  Before PR 9, ``sim50k/m4`` and ``sim100k_real/m4`` ``detection_summary.tsv`` and
  ``regions_stats.tsv`` fail it (the spec's § 2 list of goldens PR 9 changes).

No statistical claim is made (a formatting rule), so there is no calibration test (§8.4).
"""

from __future__ import annotations

import json
import math
from collections.abc import Callable
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

REPO = Path(__file__).resolve().parents[2]
TABLE_DIR = REPO / "src" / "gippa" / "schemas" / "tables"
GOLDEN = REPO / "tests" / "data" / "golden"
EXEMPT = frozenset({"NA", "INF", "-INF", "inf", "-inf", ""})
SETTINGS = settings(max_examples=80, deadline=None, suppress_health_check=[HealthCheck.too_slow])


def sig_digits(cell: str) -> int:
    """Significant digits of a numeric cell (sign, exponent, point, leading/trailing 0 removed)."""
    mantissa = cell.strip().lstrip("+-").lower().split("e")[0].replace(".", "")
    return len(mantissa.lstrip("0").rstrip("0"))


def schema_fields(table: str) -> list[dict[str, Any]]:
    doc = json.loads((TABLE_DIR / f"{table}.schema.json").read_text(encoding="utf-8"))
    return list(doc["fields"])


def number_columns(table: str) -> frozenset[str]:
    return frozenset(f["name"] for f in schema_fields(table) if f["type"] == "number")


def read_cells(path: Path) -> tuple[list[str], list[list[str]]]:
    lines = path.read_text(encoding="utf-8").split("\n")
    if lines and lines[-1] == "":
        lines.pop()
    return lines[0].split("\t"), [ln.split("\t") for ln in lines[1:]]


def violations(path: Path, table: str) -> list[tuple[int, str, str]]:
    """(row, column, cell) of every ``number`` cell with more than 6 significant digits."""
    nums = number_columns(table)
    header, rows = read_cells(path)
    out = []
    for i, row in enumerate(rows, start=1):
        for col, cell in zip(header, row, strict=True):
            if col in nums and cell not in EXEMPT and sig_digits(cell) > 6:
                out.append((i, col, cell))
    return out


def test_sig_digits_helper() -> None:
    assert sig_digits("0.00012345") == 5
    assert sig_digits("0.000648689") == 6
    assert sig_digits("5.8708e-07") == 5
    assert sig_digits("-12.383906517164503") == 17
    assert sig_digits("1234567") == 7
    assert sig_digits("4") == 1
    assert sig_digits("1e+20") == 1


# ---------------------------------------------------------------- the spec's before → after


def _summary_frame(values: list[tuple[float, float, float]]) -> pd.DataFrame:
    from gippa.context.summary import SUMMARY_COLUMNS

    rows = [("basic", "tier", f"g{i}", "fdr", v, lo, hi, 10, None)
            for i, (v, lo, hi) in enumerate(values)]  # fmt: skip
    return pd.DataFrame(rows, columns=list(SUMMARY_COLUMNS))


def _write_summary(tmp_path: Path, df: pd.DataFrame) -> Path:
    from gippa.io.summary_tsv import write_detection_summary

    path = tmp_path / "detection_summary.tsv"
    write_detection_summary(path, df)
    return path


@pytest.mark.parametrize(
    "value, want",
    [
        (0.0006486894965844569, "0.000648689"),  # A, tier fdr at q 0.001
        (5.870803769695086e-07, "5.8708e-07"),  # A, fdr_2pep
        (4.000000000000001, "4"),  # C
        (0.00012345, "0.00012345"),  # 8 decimals but 5 significant digits: unchanged
    ],
)
def test_detection_summary_spec_examples(tmp_path: Path, value: float, want: str) -> None:
    path = _write_summary(tmp_path, _summary_frame([(value, value, value)]))
    header, rows = read_cells(path)
    got = dict(zip(header, rows[0], strict=True))
    assert (got["value"], got["ci_low"], got["ci_high"]) == (want, want, want)


def test_detection_summary_na_and_infinities(tmp_path: Path) -> None:
    from gippa.io.tsv import format_value

    path = _write_summary(tmp_path, _summary_frame([(math.nan, -math.inf, math.inf)]))
    header, rows = read_cells(path)
    got = dict(zip(header, rows[0], strict=True))
    assert got["value"] == "NA"
    assert (got["ci_low"], got["ci_high"]) == (format_value(-math.inf), format_value(math.inf))
    assert got["n"] == "10"


def _regions_row(**over: Any) -> dict[str, Any]:
    from gippa.context.regions_table import REGIONS_STATS_COLUMNS

    row: dict[str, Any] = {
        "region_id": "prophage_1", "region_type": "prophage", "subset": "all", "contig": "chr1",
        "start": 1114547, "end": 1158000, "n_genes": 40, "n_eligible": 38,
        "observed": 30.0, "expected": 17.6, "o_minus_e": 12.4, "o_over_e": 1.7,
        "hl_shift_log2": 0.25, "stat": "o_minus_e", "effect": 12.4, "ci_low": 8.0,
        "ci_high": 16.0, "p_value": 0.01, "q_value": 0.02, "method": "circular_shift",
        "null_model": "circular_shift", "n_perm": 999, "label": "x", "flags": None,
    }  # fmt: skip
    row.update(over)
    assert set(row) == set(REGIONS_STATS_COLUMNS)
    return row


def _write_regions(tmp_path: Path, rows: list[dict[str, Any]]) -> Path:
    from gippa.context.regions_table import REGIONS_STATS_COLUMNS
    from gippa.io.regions_tsv import write_regions_stats

    path = tmp_path / "regions_stats.tsv"
    write_regions_stats(path, pd.DataFrame(rows, columns=list(REGIONS_STATS_COLUMNS)))
    return path


def test_regions_stats_spec_examples(tmp_path: Path) -> None:
    """A ``regions_stats.tsv``: ``42.38390651716451`` → ``42.3839``, ``-12.383906517164503`` →
    ``-12.3839``; the 7-digit ``end`` (an integer column) keeps every digit."""
    path = _write_regions(tmp_path, [_regions_row(observed=42.38390651716451,
                                                  o_minus_e=-12.383906517164503,
                                                  effect=-12.383906517164503)])  # fmt: skip
    header, rows = read_cells(path)
    got = dict(zip(header, rows[0], strict=True))
    assert got["observed"] == "42.3839"
    assert got["o_minus_e"] == "-12.3839"
    assert got["effect"] == "-12.3839"
    assert got["start"] == "1114548"  # 1-based, integer column, not subject to the rule
    assert got["end"] == "1158000"
    assert got["n_perm"] == "999"


# ---------------------------------------------------------------- property over writers


def _detection_summary(path: Path, df: pd.DataFrame) -> None:
    from gippa.io.summary_tsv import write_detection_summary

    write_detection_summary(path, df)


def _regions_stats(path: Path, df: pd.DataFrame) -> None:
    from gippa.io.regions_tsv import write_regions_stats

    write_regions_stats(path, df)


def _autocorrelation(path: Path, df: pd.DataFrame) -> None:
    from gippa.io.context_tsv import write_autocorrelation

    write_autocorrelation(path, df)


def _islands(path: Path, df: pd.DataFrame) -> None:
    from gippa.io.context_tsv import write_islands

    write_islands(path, df)


WRITERS: dict[str, Callable[[Path, pd.DataFrame], None]] = {
    "detection_summary": _detection_summary,
    "regions_stats": _regions_stats,
    "autocorrelation": _autocorrelation,
    "islands": _islands,
}


def _frame_from_schema(table: str, numbers: list[float]) -> pd.DataFrame:
    """One row per entry of ``numbers``: every ``number`` field holds it, every other field a
    fixed valid value (the first enum value, ``1`` for integers, no flags)."""
    fields = schema_fields(table)
    rows = []
    for x in numbers:
        row: dict[str, Any] = {}
        for f in fields:
            name, kind = f["name"], f["type"]
            enum = f.get("constraints", {}).get("enum")
            if name == "flags":
                row[name] = None
            elif name == "genes":
                row[name] = ("G1",)
            elif kind == "number":
                row[name] = x
            elif kind == "integer":
                row[name] = 1
            elif enum:
                row[name] = enum[0]
            else:
                row[name] = "x"
        rows.append(row)
    return pd.DataFrame(rows, columns=[f["name"] for f in fields])


finite_or_not = st.one_of(
    st.floats(allow_nan=False, allow_infinity=False),
    st.floats(min_value=-1e3, max_value=1e3, allow_nan=False),
    st.floats(min_value=1e-12, max_value=1e-2),
    st.sampled_from([math.nan, math.inf, -math.inf, 0.0, -0.0, 4.000000000000001, 1e20]),
)


@pytest.mark.parametrize("table", sorted(WRITERS))
@SETTINGS
@given(values=st.lists(finite_or_not, min_size=1, max_size=6))
def test_number_cells_are_format_value(table: str, values: list[float], tmp_path_factory:
                                       pytest.TempPathFactory) -> None:  # fmt: skip
    """Property (T-01): every ``number`` cell a published writer emits is
    ``gippa.io.tsv.format_value`` of the float, so it never has more than 6 significant
    digits."""
    from gippa.io.tsv import format_value

    path = tmp_path_factory.mktemp(table) / f"{table}.tsv"
    WRITERS[table](path, _frame_from_schema(table, values))
    nums = number_columns(table)
    header, rows = read_cells(path)
    assert len(rows) == len(values)
    for x, row in zip(values, rows, strict=True):
        want = format_value(float(x))
        for col, cell in zip(header, row, strict=True):
            if col in nums:
                assert cell == want, (table, col, x, cell, want)
                assert cell in EXEMPT or sig_digits(cell) <= 6, (table, col, cell)


@pytest.mark.parametrize("table", ["detection_summary", "regions_stats"])
@SETTINGS
@given(n=st.integers(min_value=0, max_value=10**12))
def test_integer_columns_keep_every_digit(table: str, n: int, tmp_path_factory:
                                          pytest.TempPathFactory) -> None:  # fmt: skip
    """Edge case: integer columns (``n``, ``end``, ``n_perm``, …) are not subject to the rule."""
    path = tmp_path_factory.mktemp(table) / f"{table}.tsv"
    df = _frame_from_schema(table, [0.5])
    ints = [f["name"] for f in schema_fields(table) if f["type"] == "integer"]
    for c in ints:
        df[c] = n
    WRITERS[table](path, df)
    header, rows = read_cells(path)
    got = dict(zip(header, rows[0], strict=True))
    for c in ints:
        want = n + 1 if c == "start" else n  # regions_stats writes start 1-based
        assert got[c] == str(want), (table, c, got[c])


# ---------------------------------------------------------------- committed goldens


def _goldens() -> list[Path]:
    return sorted(p for p in GOLDEN.rglob("*.tsv")
                  if (TABLE_DIR / f"{p.stem}.schema.json").is_file())  # fmt: skip


@pytest.mark.parametrize("path", _goldens(), ids=lambda p: str(p.relative_to(GOLDEN)))
def test_golden_tables_follow_6g(path: Path) -> None:
    """The committed goldens are writer output: none holds a ``number`` cell with more than 6
    significant digits (PR 9 regenerates ``sim50k/m4`` and ``sim100k_real/m4``
    ``detection_summary.tsv`` and ``regions_stats.tsv``, values only)."""
    bad = violations(path, path.stem)
    assert not bad, f"{len(bad)} cells over %.6g in {path.relative_to(REPO)}, e.g. {bad[:3]}"


def test_golden_scan_covers_the_two_offending_tables() -> None:
    names = {p.name for p in _goldens()}
    assert {"detection_summary.tsv", "regions_stats.tsv"} <= names
