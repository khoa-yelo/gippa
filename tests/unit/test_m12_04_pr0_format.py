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


_HYPHEN_AS_SIGN = re.compile(r"(?<![eE])-")


@given(st.floats(allow_nan=False, allow_infinity=False), st.integers(min_value=1, max_value=6))
def test_fmt_num_property_never_ascii_minus_sign(x: float, sig: int) -> None:
    out = style.fmt_num(x, sig=sig)
    assert _HYPHEN_AS_SIGN.search(out) is None, out


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
    assert _HYPHEN_AS_SIGN.search(out) is None, out
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
