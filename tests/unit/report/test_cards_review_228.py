"""PR #228 review B1 (lead decision 15), revised by M12-04 PR 1 (AC 1.4).

#228 had the peptide card state the placement claim per family and list every ``mapping_fdr``
stratum row as a descriptive value with its full key. M12-04 PR 1 moves the per-stratum FDR string
off the cards: it stays in V06 and ``mapping/mapping_fdr.tsv``, and the report text no longer
states the per-family claim anywhere. These tests pin that a run with many strata (two decoy
types, three difference counts, two length classes) still gives a short card with only the
placement threshold, and that none of the stratum rows reaches it.
"""

from __future__ import annotations

import re
from pathlib import Path

from helpers.report_content import page, words, write_mapping
from helpers.report_fixtures import Page, run_out
from helpers.viz_out_tree import write_schema_table


def _rows() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for d, base in (("reverse", 0.001), ("shuffle", 0.002)):
        for n in (0, 1, 2):
            for k, lc in enumerate(("7-8", "9-12")):
                rows.append(
                    {
                        "section": "placement",
                        "length_class": lc,
                        "n_diff": n,
                        "targets": 100,
                        "decoys": 1,
                        "n_accepted": 90,
                        "decoy_type": d,
                        "fdr_stratum": round(base + 0.0001 * (10 * n + k), 6),
                    }
                )
    return rows  # fmt: skip


def _page(tmp: Path) -> Page:
    out = run_out(tmp)
    write_mapping(out, placed_q=[0.001], unmapped_peptides=0, proteins=["P1"],
                  unmapped_proteins=0, strata={"7-8": 0.0072})  # fmt: skip
    write_schema_table(out / "mapping" / "mapping_fdr.tsv", "mapping_fdr", _rows())
    return page(out)


def _card(p: Page) -> str:
    return words(p.cards()["peptides_placed"].text())


def test_peptide_card_states_only_the_placement_threshold(tmp_path: Path) -> None:
    """AC 1.4: with 12 stratum rows the card is still short and names no stratum."""
    text = _card(_page(tmp_path))
    assert re.search(r"placement q ≤ 0\.010*(?![\d])", text), text
    assert "stratum" not in text and "not q-values" not in text, text
    assert "per family" not in text and "decoy type" not in text, text
    assert len(text) <= 160, text


def test_no_stratum_row_reaches_the_card(tmp_path: Path) -> None:
    """No stratum key or value is on the card, and the card no longer reads ``mapping_fdr``;
    the stratum detail is V06's (the mapping section keeps its slot)."""
    p = _page(tmp_path)
    text = _card(p)
    assert "diff " not in text, text
    for r in _rows():
        key = f"{r['decoy_type']} {r['n_diff']} diff {r['length_class']}"
        assert key not in text, (key, text)
        assert r["decoy_type"] not in text, (r, text)
    src = p.cards()["peptides_placed"].attrs.get("data-source", "")
    assert "mapping_fdr.tsv" not in src, src
    assert p.section("mapping").by_attr("data-plot", "V06"), "V06 holds the stratum FDR"
