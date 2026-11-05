"""M12-04 PR 0: shared number helpers, the abundance label and ``Instance.empty``
(docs/specs/M12-04-outputs-revision.md § PR 0, AC 0.1, 0.4, 0.5; ADR-0075).

Interface assumed beyond the spec text:
- ``fmt_q``, ``fmt_num``, ``fmt_pct`` and ``abundance_label`` live in ``gippa.viz.style`` (the
  spec names ``gippa.viz.style.fmt_q`` and lists the others in the same block).
- Integers (``int`` and numpy integers) are printed exactly, not rounded to ``sig`` digits: the
  spec's examples are counts (``4,373`` genes, ``200,000`` bp, ``158`` genes).
- Grouping of integers >= 10,000 is a comma, a thin space (U+2009) or a narrow no-break space
  (U+202F); the spec allows "a thin-space or comma grouping". Integers below 10,000 may or may
  not be grouped (the spec's example ``4,373`` sits below its own threshold).
- NA is ``None``, ``float('nan')``, ``numpy.nan`` or ``pd.NA``.
- Review round 1 (PR #549) pins the comma: ``10000`` -> ``"10,000"``, ``4373`` -> ``"4373"``.
"""

from __future__ import annotations

import math
import re

import numpy as np
import pandas as pd
import pytest
from hypothesis import given
from hypothesis import strategies as st
from matplotlib.figure import Figure

from gippa.viz import base, style

MINUS = "\u2212"
GROUP_SEPARATORS = (",", "\u2009", "\u202f")
NA_VALUES = [None, float("nan"), np.nan, pd.NA]


def _parse(text: str) -> float:
    """The number a formatted string stands for: true minus to ASCII, grouping removed."""
    out = text.replace(MINUS, "-")
    for sep in GROUP_SEPARATORS:
        out = out.replace(sep, "")
    return float(out)


def _half_unit(x: float, sig: int) -> float:
    """Half a unit in the ``sig``-th significant digit of ``x`` (the rounding bound)."""
    return 0.5 * 10 ** (math.floor(math.log10(abs(x))) - sig + 1)


def _sig_digits(text: str) -> int:
    """Significant digits of a plain decimal string, trailing zeros counted."""
    digits = re.sub(r"[^0-9]", "", text.replace(MINUS, ""))
    return len(digits.lstrip("0"))


# --- AC 0.1: fmt_q -------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("q", "expected"),
    [
        (0.432, "0.43"),  # spec example
        (1, "1.0"),  # spec example: trailing zero kept
        (1.0, "1.0"),
        (0.05, "0.050"),  # spec example: trailing zero kept
        (0.001, "0.0010"),  # edge case: exactly 0.001 is not "< 0.001"
        (0.0009999, "< 0.001"),
        (0.0001, "< 0.001"),
        (1e-300, "< 0.001"),
        (0.0, "< 0.001"),  # edge case: q = 0
        (0.00123, "0.0012"),
        (0.0124093, "0.012"),  # the run-A Mayu FDR of the PR 1 example
        (0.414585, "0.41"),  # the run-B null p of the PR 1 example
        (0.5, "0.50"),
        (0.999, "1.0"),
    ],
)
def test_fmt_q_examples(q: float, expected: str) -> None:
    assert style.fmt_q(q) == expected


@pytest.mark.parametrize("na", NA_VALUES)
def test_fmt_q_na(na: object) -> None:
    assert style.fmt_q(na) == "NA"


def test_fmt_q_accepts_numpy_floats() -> None:
    assert style.fmt_q(np.float64(0.432)) == "0.43"
    assert style.fmt_q(np.float32(0.05)) == "0.050"


@given(st.floats(min_value=0.001, max_value=1.0, allow_nan=False))
def test_fmt_q_property_two_significant_digits_and_round_trip(q: float) -> None:
    out = style.fmt_q(q)
    assert not out.startswith("<")
    assert "-" not in out and MINUS not in out
    assert _sig_digits(out) == 2, out
    assert abs(_parse(out) - q) <= _half_unit(q, 2) * (1 + 1e-9), out


@given(st.floats(min_value=0.0, max_value=0.001, exclude_max=True, allow_nan=False))
def test_fmt_q_property_below_threshold(q: float) -> None:
    assert style.fmt_q(q) == "< 0.001"


# --- AC 0.1: fmt_num -------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("x", "sig", "expected"),
    [
        (-0.147806, 2, f"{MINUS}0.15"),  # spec example, true minus U+2212
        (14.2137, 3, "14.2"),  # spec example: sig=3 for log2 abundances
        (-0.0174037, 2, f"{MINUS}0.017"),  # the run-B model R² mean of the PR 1 example
        (0.0656274, 2, "0.066"),  # its SD
        (0.5, 2, "0.50"),  # trailing zeros kept (as for fmt_q)
        (-0.0, 2, "0.0"),  # edge case: negative zero has no minus
        (0.0, 2, "0.0"),
    ],
)
def test_fmt_num_examples(x: float, sig: int, expected: str) -> None:
    assert style.fmt_num(x, sig=sig) == expected


def test_fmt_num_default_sig_is_two() -> None:
    assert style.fmt_num(-0.147806) == f"{MINUS}0.15"
    assert style.fmt_num(0.432) == "0.43"


@pytest.mark.parametrize("na", NA_VALUES)
def test_fmt_num_na(na: object) -> None:
    assert style.fmt_num(na) == "NA"
    assert style.fmt_num(na, sig=3) == "NA"


@pytest.mark.parametrize(("n", "digits"), [(10_000, "10000"), (200_000, "200000"),
                                           (1_234_567, "1234567")])  # fmt: skip
def test_fmt_num_groups_integers_from_ten_thousand(n: int, digits: str) -> None:
    out = style.fmt_num(n)
    assert any(sep in out for sep in GROUP_SEPARATORS), out
    assert re.sub(r"[^0-9]", "", out) == digits
    # groups of three from the right
    assert re.fullmatch(r"\d{1,3}([,\u2009\u202f]\d{3})+", out), out


def test_fmt_num_run_box_examples() -> None:
    assert style.fmt_num(200_000) in ("200,000", "200\u2009000", "200\u202f000")
    assert style.fmt_num(4_373) in ("4,373", "4373", "4\u2009373", "4\u202f373")


@pytest.mark.parametrize("n", [0, 4, 21, 99, 158, 774, 9_999])
def test_fmt_num_small_integers_are_exact(n: int) -> None:
    out = style.fmt_num(n)
    assert _parse(out) == n
    assert re.fullmatch(r"\d[\d,\u2009\u202f]*", out), out


def test_fmt_num_negative_integer_uses_true_minus() -> None:
    out = style.fmt_num(-4_373)
    assert out.startswith(MINUS)
    assert _parse(out) == -4_373


def test_fmt_num_numpy_integers_match_python_integers() -> None:
    assert style.fmt_num(np.int64(200_000)) == style.fmt_num(200_000)
    assert style.fmt_num(np.int64(158)) == style.fmt_num(158)


@given(st.floats(allow_nan=False, allow_infinity=False), st.integers(min_value=1, max_value=6))
def test_fmt_num_property_never_ascii_minus_sign(x: float, sig: int) -> None:
    """Review M4: no ASCII hyphen anywhere, the exponent included."""
    out = style.fmt_num(x, sig=sig)
    assert "-" not in out, out


@given(st.floats())
def test_fmt_q_property_never_ascii_minus_sign(q: float) -> None:
    assert "-" not in style.fmt_q(q)


@given(st.integers(min_value=-(10**12), max_value=10**12))
def test_fmt_num_property_integers_round_trip_exactly(n: int) -> None:
    out = style.fmt_num(n)
    assert "-" not in out, out
    assert _parse(out) == n, out


@given(
    st.floats(min_value=1e-3, max_value=1e4, exclude_max=True, allow_nan=False),
    st.booleans(),
    st.integers(min_value=1, max_value=5),
)
def test_fmt_num_property_round_trip_within_sig_rounding(
    magnitude: float, negative: bool, sig: int
) -> None:
    x = -magnitude if negative else magnitude
    out = style.fmt_num(x, sig=sig)
    assert "-" not in out, out
    assert out.startswith(MINUS) == negative, out
    assert abs(_parse(out) - x) <= _half_unit(x, sig) * (1 + 1e-9), out


# --- fmt_pct and the title_pct alias -----------------------------------------------------------


@pytest.mark.parametrize(
    ("n", "d", "expected"),
    [
        (99, 158, "63 %"),  # run A: 99 of 158 genes detected
        (0, 0, "NA"),  # d = 0
        (5, 0, "NA"),
        (0, 10, "0 %"),
        (10, 10, "100 %"),
        (1, 1000, "1 %"),  # a partial share never rounds to 0 %
        (999, 1000, "99 %"),  # nor to 100 %
    ],
)
def test_fmt_pct(n: int, d: int, expected: str) -> None:
    assert style.fmt_pct(n, d) == expected


@given(st.integers(min_value=0, max_value=10_000), st.integers(min_value=0, max_value=10_000))
def test_title_pct_stays_an_alias_of_fmt_pct(n: int, d: int) -> None:
    assert base.title_pct(n, d) == style.fmt_pct(n, d)


# --- AC 0.5: abundance_label -------------------------------------------------------------------


@pytest.mark.parametrize(
    ("atype", "scale", "expected"),
    [
        ("ibaq", "molar", "log2 iBAQ (molar)"),
        ("maxlfq", "relative", "log2 MaxLFQ (relative)"),
        ("top3", "relative", "log2 Top3 (relative)"),
        ("top3", "molar", "log2 Top3 (molar)"),
        ("intensity", "relative", "log2 intensity (relative)"),
        ("intensity", "molar", "log2 intensity (molar)"),
        ("unknown", "relative", "log2 abundance (measure unknown)"),
        ("unknown", "molar", "log2 abundance (measure unknown)"),
        ("unknown", None, "log2 abundance (measure unknown)"),
    ],
)
def test_abundance_label_table(atype: str, scale: str | None, expected: str) -> None:
    assert style.abundance_label(atype, scale) == expected


@pytest.mark.parametrize("scale", [None, float("nan"), pd.NA, "relative", "molar"])
def test_abundance_label_none_type_gives_none(scale: object) -> None:
    assert style.abundance_label("none", scale) is None


@pytest.mark.parametrize("atype", [None, float("nan"), pd.NA])
@pytest.mark.parametrize("scale", [None, pd.NA, "relative"])
def test_abundance_label_na_type_gives_none(atype: object, scale: object) -> None:
    assert style.abundance_label(atype, scale) is None


# --- AC 0.4: Instance.empty --------------------------------------------------------------------


def test_empty_state_instance_is_empty() -> None:
    inst = base.empty_state("V05", "sankey", "No peptides in this run.", ["a", "b"], "single")
    assert inst.empty is True


def test_ordinary_instance_is_not_empty() -> None:
    inst = base.Instance(suffix=None, figure=Figure(), data=pd.DataFrame(), title="t")
    assert inst.empty is False


def test_instance_empty_can_be_set_and_is_frozen() -> None:
    inst = base.Instance(suffix=None, figure=Figure(), data=pd.DataFrame(), title="t", empty=True)
    assert inst.empty is True
    with pytest.raises(AttributeError):
        inst.empty = False  # type: ignore[misc]


# --- review round 1 (PR #549) ------------------------------------------------------------------
# Lead decisions: M4 (true minus in the exponent), M5 (fmt_count), M6 (fmt_alpha), L items and
# the mutation gaps. ``fmt_count`` and ``fmt_alpha`` are assumed to live in gippa.viz.style next
# to fmt_num.


@pytest.mark.parametrize(
    ("x", "sig", "expected"),
    [
        (1.2e-12, 2, f"1.2e{MINUS}12"),  # M4: U+2212 in the exponent
        (-1.2e-12, 2, f"{MINUS}1.2e{MINUS}12"),
        (1.2e-5, 2, f"1.2e{MINUS}5"),  # no leading zero: never "e<minus>05"
        (2e6, 2, "2.0e6"),  # no "+" and no leading zero: never "e+06"
        (1e6, 2, "1.0e6"),  # 1e6 exactly (a float) is outside the plain range
        (999999.5, 2, "1.0e6"),  # rounds up to 1e6, so scientific
        (999999.5, 7, "999999.5"),  # below 1e6 after rounding: plain
        (1e-4, 2, "0.00010"),  # 1e-4 exactly is inside the plain range
        (9.99e-5, 3, f"9.99e{MINUS}5"),  # just below 1e-4: scientific
        (9.9e-5, 2, f"9.9e{MINUS}5"),
        (9.99e-5, 2, "0.00010"),  # rounds up to 1.0e-4, so plain (the range is after rounding)
        (9.995, 2, "10"),  # rounding boundaries at 2 significant digits
        (99.95, 2, "100"),
        (0.0995, 2, "0.10"),
    ],
)
def test_fmt_num_scientific_bounds_and_rounding(x: float, sig: int, expected: str) -> None:
    assert style.fmt_num(x, sig=sig) == expected


@given(
    st.one_of(
        st.floats(min_value=1e-300, max_value=1e-5),
        st.floats(min_value=2e6, max_value=1e300),
    ),
    st.booleans(),
)
def test_fmt_num_property_exponent_has_no_plus_or_leading_zero(x: float, negative: bool) -> None:
    out = style.fmt_num(-x if negative else x)
    assert "e" in out, out
    exp = out.split("e", 1)[1]
    assert re.fullmatch(rf"{MINUS}?[1-9]\d*", exp), out


@pytest.mark.parametrize("x", [math.inf, -math.inf, np.inf, -np.inf])
def test_fmt_num_and_fmt_q_infinity_is_na(x: float) -> None:
    assert style.fmt_num(x) == "NA"
    assert style.fmt_num(x, sig=3) == "NA"
    assert style.fmt_q(x) == "NA"


@pytest.mark.parametrize(("n", "expected"), [(4_373, "4373"), (9_999, "9999"),
                                             (10_000, "10,000"), (10_001, "10,001")])  # fmt: skip
def test_fmt_num_grouping_threshold_is_ten_thousand(n: int, expected: str) -> None:
    assert style.fmt_num(n) == expected


@pytest.mark.parametrize("bad", ["NA", "abc", "", pd.NaT, None])
def test_fmt_num_and_fmt_q_non_numbers_are_na(bad: object) -> None:
    assert style.fmt_num(bad) == "NA"
    assert style.fmt_q(bad) == "NA"


@pytest.mark.parametrize("b", [True, False])
def test_fmt_num_numpy_bool_is_a_bool(b: bool) -> None:
    assert style.fmt_num(np.bool_(b)) == style.fmt_num(b)


@pytest.mark.parametrize("x", [1.5, 158, np.int64(3)])
@pytest.mark.parametrize("sig", [0, -1])
def test_fmt_num_sig_below_one_raises(x: object, sig: int) -> None:
    with pytest.raises(ValueError):
        style.fmt_num(x, sig=sig)


# --- M5: fmt_count -----------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("x", "expected"),
    [
        (0, "0"),
        (158, "158"),
        (np.int64(158), "158"),
        (np.int32(99), "99"),
        (158.0, "158"),  # integer-valued float
        (np.float64(774.0), "774"),
        (4_373, "4373"),
        (9_999, "9999"),
        (10_000, "10,000"),  # grouped from 10,000
        (200_000.0, "200,000"),
        (1_234_567, "1,234,567"),
    ],
)
def test_fmt_count_examples(x: object, expected: str) -> None:
    assert style.fmt_count(x) == expected


@pytest.mark.parametrize("na", NA_VALUES)
def test_fmt_count_na(na: object) -> None:
    assert style.fmt_count(na) == "NA"


@pytest.mark.parametrize("x", [158.5, 0.1, np.float64(2.25), -3.5])
def test_fmt_count_non_integer_float_raises(x: float) -> None:
    with pytest.raises(ValueError):
        style.fmt_count(x)


@given(st.integers(min_value=0, max_value=10**12))
def test_fmt_count_property_matches_fmt_num_on_integers(n: int) -> None:
    assert style.fmt_count(n) == style.fmt_num(n)
    assert style.fmt_count(float(n)) == style.fmt_num(n)  # 10**12 < 2**53: exact as a float


# --- M6: fmt_alpha -----------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("alpha", "expected"),
    [(0.01, "0.01"), (0.05, "0.05"), (0.025, "0.025"), (0.1, "0.1"), (0.001, "0.001"),
     (np.float64(0.05), "0.05")],
)  # fmt: skip
def test_fmt_alpha_shortest_exact_decimal(alpha: float, expected: str) -> None:
    assert style.fmt_alpha(alpha) == expected


@given(st.floats(min_value=0.001, max_value=1.0))
def test_fmt_alpha_property_round_trips_exactly(alpha: float) -> None:
    out = style.fmt_alpha(alpha)
    assert "e" not in out.lower() and "-" not in out, out
    assert float(out) == alpha, out


# --- L items: abundance_label ------------------------------------------------------------------


def test_abundance_label_empty_scale_has_no_empty_parentheses() -> None:
    out = style.abundance_label("ibaq", "")
    assert out is not None and "()" not in out
    assert out == "log2 iBAQ"


@pytest.mark.parametrize("scale", [None, float("nan"), np.nan, pd.NA])
@pytest.mark.parametrize(
    ("atype", "name"),
    [("ibaq", "iBAQ"), ("maxlfq", "MaxLFQ"), ("top3", "Top3"), ("intensity", "intensity")],
)
def test_abundance_label_known_type_na_scale(atype: str, name: str, scale: object) -> None:
    assert style.abundance_label(atype, scale) == f"log2 {name}"


@pytest.mark.parametrize(
    ("atype", "scale", "expected"),
    [
        ("IBAQ", "molar", "log2 iBAQ (molar)"),
        ("iBAQ", "molar", "log2 iBAQ (molar)"),
        ("MaxLFQ", "relative", "log2 MaxLFQ (relative)"),
        ("Top3", "relative", "log2 Top3 (relative)"),
        ("Intensity", "molar", "log2 intensity (molar)"),
        ("Unknown", "relative", "log2 abundance (measure unknown)"),
        ("NONE", "relative", None),
        ("None", None, None),
    ],
)
def test_abundance_label_type_is_case_insensitive(
    atype: str, scale: str | None, expected: str | None
) -> None:
    assert style.abundance_label(atype, scale) == expected


def test_abundance_label_docstring_states_the_scale_vocabulary() -> None:
    doc = style.abundance_label.__doc__ or ""
    assert "molar" in doc and "relative" in doc
    assert style.abundance_label("ibaq", "molar") == "log2 iBAQ (molar)"
    assert style.abundance_label("ibaq", "relative") == "log2 iBAQ (relative)"
