"""V12 pooled class-error labels (M9-18 amendment; statistician note local/stat/v12-class-error.md):
the site line names the class n, the title says the rate is pooled only when a rate is numeric,
and a rate of exactly 1 on an ALT site gets the ``*`` mark and the bias-1 footnote."""

from __future__ import annotations

import math

import pandas as pd

from gippa.viz import style
from gippa.viz.variants import DEGENERATE_NOTE, POOLED, _header, title
from gippa.viz.words import plain


def _site(rate: object, n: object) -> dict[str, object]:
    return {"gene": "g1", "aa_pos": 81, "contig": "chr1", "codon_start": 100, "strand": "+",
            "ref_aa": "Y", "alt_aa": "D", "alt_type": "substitution", "confidence": "high",
            "n_alt": 8, "n_ref": 0, "class_error_rate": rate, "class_error_n": n}  # fmt: skip


def _alt(rates: list[object], call: str = "ALT") -> pd.DataFrame:
    n = len(rates)
    return pd.DataFrame({"call": [call] * n, "confidence": ["high"] * n,
                         "class_error_rate": rates})  # fmt: skip


def test_site_line_names_the_class_n() -> None:
    assert _header(_site(0.25, 7)).endswith("class error 0.25 (pooled, class of 7 paired reads)")


def test_site_line_with_na_n() -> None:
    assert _header(_site(0.25, None)).endswith("(pooled, class of NA paired reads)")


def test_unestimated_site_has_no_pooled_suffix() -> None:
    """M12-04 PR 11c: an NA rate reads as the flag's plain words."""
    text = _header(_site(math.nan, 7))
    assert text.endswith(plain("class_error_unestimated")) and "pooled" not in text


def test_pooled_clause_only_with_a_numeric_rate() -> None:
    assert POOLED in title(_alt([0.2, 0.4]))
    assert POOLED not in title(_alt([math.nan]))


def test_rate_one_marks_the_title_and_below_one_does_not() -> None:
    assert f"{style.fmt_num(1.0)}*" in title(_alt([0.2, 1.0]))  # "1.0*" (M12-04 P-18)
    assert "*" not in title(_alt([0.2, 0.9]))


def test_footnote_text_names_bias_1() -> None:
    assert DEGENERATE_NOTE.startswith("* ") and "bias 1" in DEGENERATE_NOTE
