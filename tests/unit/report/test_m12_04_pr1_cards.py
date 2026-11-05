# ruff: noqa: RUF001, RUF002, RUF003  (the report's wording is Unicode: ≤, ≥, ², −, ·)
"""M12-04 PR 1, behaviour 3: summary cards that follow the inputs (AC 1.1 to 1.5;
docs/specs/M12-04-outputs-revision.md § PR 1, R-02).

Fixtures P (peptides, mayu), R (proteins only, protein_list) and D (de novo, two_peptide, no
abundance) are in ``test_m12_04_pr1_fixtures``.

Interface assumed beyond the spec text (PR 1 test author):
- card names stay the M9-16 ones (``genes_detected``, ``peptides_placed``, ``proteins_placed``,
  ``model_r2``); "the card whose name is peptides" of AC 1.1 is ``peptides_placed``;
- "the run had peptides / proteins" means the manifest inputs list that role; a role that is
  present always gets its card unless its input count is 0 (no "0 of 0");
- "3 significant digits on cards" against the examples' 2 ("0.012", "−0.017 ± 0.066"): a card
  number is accepted when it is ``fmt_num(x, 2)``, ``fmt_num(x, 3)`` or ``fmt_q(x)`` of the
  table value; AC 1.3's "≤ 3 significant digits" is checked on every number of the card;
- the "this run: 99 genes, 34 TUs" counts are the genes with a finite ``abundance_log2`` and their
  distinct ``tu_id`` (the fixture's univariate ``n`` is the same 99).
"""

from __future__ import annotations

import re
import tempfile
from pathlib import Path
from typing import Any

import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st
from unit.report.test_m12_04_pr1_fixtures import (
    FIXTURES,
    Run,
    build,
    fixture_d,
    fixture_p,
    fixture_r,
    text,
    write_fitted,
    write_peptides,
    write_proteins,
    write_run_manifest,
    write_unmapped,
)

from gippa.viz import style

ZERO_OF_ZERO = re.compile(r"\b0 (peptides|proteins) placed of 0\b")
NUMBER = re.compile(r"(?<![\w.])[−-]?\d[\d,  ]*(?:\.\d+)?")
INTERNAL = re.compile(r"ADR-\d{4}|§\s?\d|GIPPA_SPEC|docs/methods")
MINUS = "−"


def _cards(run: Run) -> dict[str, str]:
    return {k: text(v) for k, v in run.page().cards().items()}


def _card(run: Run, name: str) -> str:
    found = _cards(run)
    assert name in found, f"{run.name}: no {name} card; cards: {sorted(found)}"
    return found[name]


def _shown(x: float) -> set[str]:
    return {style.fmt_num(x, sig=2), style.fmt_num(x, sig=3), style.fmt_q(x)}


def sig_digits(token: str) -> int:
    digits = re.sub(r"[^0-9.]", "", token)
    if "." in digits:
        return len(digits.replace(".", "").lstrip("0")) or 1
    return len(digits.rstrip("0")) or 1


# ---------------------------------------------------------------- AC 1.1, 1.2: input-aware


def test_ac_1_1_proteins_only_has_no_peptide_card(tmp_path: Path) -> None:
    run = fixture_r(tmp_path)
    cards = _cards(run)
    assert "peptides_placed" not in cards, cards
    assert not any("two-peptide" in t for t in cards.values()), cards
    assert "peptides_placed" not in {c["name"] for c in run.cards()}


def test_ac_1_1_peptides_only_has_no_protein_card(tmp_path: Path) -> None:
    run = fixture_p(tmp_path)
    cards = _cards(run)
    assert "proteins_placed" not in cards, cards
    assert "proteins_placed" not in {c["name"] for c in run.cards()}
    assert "peptides_placed" in cards, cards


def test_de_novo_run_has_no_protein_card(tmp_path: Path) -> None:
    cards = _cards(fixture_d(tmp_path))
    assert "proteins_placed" not in cards and "peptides_placed" in cards, cards


@pytest.mark.parametrize("name", FIXTURES)
def test_ac_1_2_no_zero_of_zero_card(name: str, tmp_path: Path) -> None:
    cards = _cards(build(name, tmp_path))
    bad = {k: t for k, t in cards.items() if ZERO_OF_ZERO.search(t) or re.search(r"\b0 of 0\b", t)}
    assert not bad, bad


def test_both_inputs_show_both_cards(tmp_path: Path) -> None:
    run = build("both", tmp_path)
    cards = _cards(run)
    assert "peptides_placed" in cards and "proteins_placed" in cards, cards
    assert re.search(r"\b12 of 15 proteins placed\b", cards["proteins_placed"]), cards


def test_proteins_card_example_b(tmp_path: Path) -> None:
    """Spec example B: "106 of 106 proteins placed"."""
    card = _card(fixture_r(tmp_path), "proteins_placed")
    assert "106 of 106 proteins placed" in card, card


def test_peptides_card_example_c(tmp_path: Path) -> None:
    """Spec example C: "21 of 21 peptides placed · placement q ≤ 0.01"."""
    card = _card(fixture_d(tmp_path), "peptides_placed")
    assert "21 of 21 peptides placed" in card, card
    assert "placement q ≤ 0.01" in card, card


@settings(max_examples=25, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(
    pep=st.booleans(), n_pep=st.integers(0, 4), un_pep=st.integers(0, 2),
    prot=st.booleans(), n_prot=st.integers(0, 4), un_prot=st.integers(0, 2),
)  # fmt: skip
def test_property_cards_follow_the_input_roles(
    pep: bool, n_pep: int, un_pep: int, prot: bool, n_prot: int, un_prot: int
) -> None:
    """A card per input role given, none for a role not given, and never "0 of 0"."""
    n_pep, un_pep = (n_pep, un_pep) if pep else (0, 0)
    n_prot, un_prot = (n_prot, un_prot) if prot else (0, 0)
    roles = ["genome"] + (["peptides"] if pep else []) + (["proteins"] if prot else [])
    with tempfile.TemporaryDirectory() as tmp:
        run = fixture_p(Path(tmp), roles=tuple(roles))
        write_peptides(run.out, n_pep, source="database")
        write_proteins(run.out, n_prot)
        write_unmapped(run.out, n_peptides=un_pep, n_proteins=un_prot)
        cards = _cards(run)
    assert not any(ZERO_OF_ZERO.search(t) for t in cards.values()), cards
    assert ("peptides_placed" in cards) == (pep and n_pep + un_pep > 0), (roles, cards)
    assert ("proteins_placed" in cards) == (prot and n_prot + un_prot > 0), (roles, cards)


# ---------------------------------------------------------------- the genes card per tier


def test_two_peptide_line_only_in_mayu_and_decoy_picked(tmp_path: Path) -> None:
    for name, shown in (("P", True), ("decoy", True), ("R", False), ("D", False)):
        card = _card(build(name, tmp_path / name), "genes_detected")
        assert ("also pass" in card or "≥ 2 distinct" in card) == shown, (name, card)
        if shown:
            assert re.search(r"(?<![\d.])91(?![\d.])", card), (name, card)


def test_genes_card_example_a_n_of_total(tmp_path: Path) -> None:
    card = _card(fixture_p(tmp_path), "genes_detected")
    assert "99 of 158 genes detected" in card, card


def test_genes_card_example_b_protein_list_says_no_fdr(tmp_path: Path) -> None:
    """Spec example B: "106 of 158 genes detected · protein list: no FDR estimated"."""
    card = _card(fixture_r(tmp_path), "genes_detected")
    assert "106 of 158 genes detected" in card, card
    assert "no FDR estimated" in card, card
    assert "two-peptide" not in card, card


# ---------------------------------------------------------------- AC 1.3: the mayu FDR, honestly


def test_ac_1_3_mayu_card_names_the_gene_level_fdr(tmp_path: Path) -> None:
    card = _card(fixture_p(tmp_path), "genes_detected")
    assert "gene-level FDR" in card, card
    assert "q ≤ 0.01" in card, card
    assert any(v in card for v in _shown(0.0124093)), card
    assert "0.0124093" not in card and "0.0124 " not in card, card


def test_ac_1_3_mayu_fdr_is_not_shown_as_measured_at_the_peptide_threshold(tmp_path: Path) -> None:
    """Run A's card read "Mayu FDR 0.0124 at q ≤ 0.01": a gene FDR above the threshold it seemed
    to be controlled at. The q belongs to the peptides, and the card must say so."""
    card = _card(fixture_p(tmp_path), "genes_detected")
    assert not re.search(r"FDR\s+[\d.]+\s+at q ≤", card), card
    fdr = card.index("gene-level FDR")
    q = card.index("q ≤ 0.01", fdr)
    assert "peptide" in card[fdr:q], card


def test_ac_1_3_every_number_has_at_most_three_significant_digits(tmp_path: Path) -> None:
    card = _card(fixture_p(tmp_path), "genes_detected")
    tokens = NUMBER.findall(card)
    assert tokens, card
    assert all(sig_digits(t) <= 3 for t in tokens), (tokens, card)


def test_decoy_picked_fdr_is_formatted(tmp_path: Path) -> None:
    card = _card(build("decoy", tmp_path), "genes_detected")
    assert any(v in card for v in _shown(0.0098)), card
    assert all(sig_digits(t) <= 3 for t in NUMBER.findall(card)), card


# ---------------------------------------------------------------- AC 1.4: short, no strata


@pytest.mark.parametrize("name", FIXTURES)
def test_ac_1_4_no_card_longer_than_160_characters(name: str, tmp_path: Path) -> None:
    long = {k: (len(t), t) for k, t in _cards(build(name, tmp_path)).items() if len(t) > 160}
    assert not long, long


def test_ac_1_4_de_novo_card_has_no_stratum_string(tmp_path: Path) -> None:
    cards = _cards(fixture_d(tmp_path))
    assert not any("stratum" in t for t in cards.values()), cards
    assert not any(re.search(r"\bdiff\b", t) for t in cards.values()), cards


@pytest.mark.parametrize("name", ["P", "D", "both"])
def test_peptide_card_sources_drop_the_stratum_columns(name: str, tmp_path: Path) -> None:
    """The per-stratum FDR moves off the cards (it stays in V06 and mapping_fdr.tsv)."""
    found = {c["name"]: c for c in build(name, tmp_path).cards()}
    assert not any("fdr_stratum" in s for s in found["peptides_placed"]["sources"]), found


@pytest.mark.parametrize("name", FIXTURES)
def test_cards_hold_no_internal_reference(name: str, tmp_path: Path) -> None:
    cards = _cards(build(name, tmp_path))
    assert not {k: t for k, t in cards.items() if INTERNAL.search(t)}, cards


# ---------------------------------------------------------------- AC 1.5: the R² card


def test_ac_1_5_r2_card_label_first_example_b(tmp_path: Path) -> None:
    """Spec example B: "Model R²: −0.017 ± 0.066 (CV), no predictive power; null p 0.41
    (uncalibrated)"."""
    card = _card(fixture_r(tmp_path), "model_r2")
    assert card.startswith("Model R²:"), card
    assert any(f"{style.fmt_num(-0.0174037, sig=s)} ± {style.fmt_num(0.0656274, sig=s)} (CV)"
               in card for s in (2, 3)), card  # fmt: skip
    assert "no predictive power" in card, card
    assert any(f"null p {v}" in card for v in _shown(0.414585)), card
    assert "(uncalibrated)" in card, card
    assert "-0.0174037" not in card and "0.0656274" not in card, card
    assert not re.search(r"(?<![\w])-\d", card), "a true minus, not a hyphen"


def test_ac_1_5_value_comes_after_the_label(tmp_path: Path) -> None:
    card = _card(fixture_r(tmp_path), "model_r2")
    assert "Model R²:" in card and MINUS in card, card
    assert card.index("Model R²:") < card.index(MINUS), card


def test_ac_1_5_positive_r2_has_no_no_power_note(tmp_path: Path) -> None:
    run = fixture_r(tmp_path)
    write_fitted(run.out, 0.21, 0.03, 0.001)
    card = _card(run, "model_r2")
    assert card.startswith("Model R²:"), card
    assert any(f"{style.fmt_num(0.21, sig=s)} ± {style.fmt_num(0.03, sig=s)}" in card
               for s in (2, 3)), card  # fmt: skip
    assert "no predictive power" not in card, card


def test_ac_1_5_mean_minus_sd_at_or_below_zero_is_no_power(tmp_path: Path) -> None:
    run = fixture_r(tmp_path)
    write_fitted(run.out, 0.05, 0.06, 0.2)
    assert "no predictive power" in _card(run, "model_r2")


@settings(max_examples=20, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(r2=st.integers(-500, 900), sd=st.integers(0, 500))
def test_property_no_power_iff_mean_or_mean_minus_sd_at_most_zero(r2: int, sd: int) -> None:
    """The PR 7b rule on the card: m ≤ 0 or m − s ≤ 0 (values in thousandths, exact in %.6g)."""
    m, s = r2 / 1000, sd / 1000
    with tempfile.TemporaryDirectory() as tmp:
        run = fixture_r(Path(tmp))
        write_fitted(run.out, m, s, 0.3)
        card = _card(run, "model_r2")
    assert card.startswith("Model R²:"), card
    assert ("no predictive power" in card) == (r2 <= 0 or r2 - sd <= 0), (m, s, card)


def test_ac_1_5_not_fitted_states_both_minimums(tmp_path: Path) -> None:
    from gippa.determinants.models import MIN_ROWS, MIN_TUS

    card = _card(fixture_p(tmp_path), "model_r2")
    assert card.startswith("Model R²:"), card
    assert "not fitted" in card, card
    assert re.search(rf"≥ {MIN_ROWS} genes", card), card
    assert re.search(rf"≥ {MIN_TUS} TUs", card), card
    assert "unestimated" not in card, card


def test_ac_1_5_not_fitted_states_this_runs_counts(tmp_path: Path) -> None:
    """Spec example A: "(this run: 99 genes, 34 TUs)"."""
    card = _card(fixture_p(tmp_path), "model_r2")
    assert re.search(r"this run: 99 genes, 34 TUs", card), card


def test_ac_1_5_no_abundance_says_so(tmp_path: Path) -> None:
    """Spec example C: "Model R²: not fitted: no abundance in this run"."""
    card = _card(fixture_d(tmp_path), "model_r2")
    assert card.startswith("Model R²:"), card
    assert "not fitted" in card and "no abundance" in card, card
    assert "NA" not in card.split(), card


# ---------------------------------------------------------------- manifest without inputs


def test_old_manifest_without_inputs_still_builds_cards(tmp_path: Path) -> None:
    """An old run (no manifest ``inputs``): the cards fall back to the tables, and an empty
    table still gives no "0 of 0" card."""
    run = fixture_r(tmp_path)
    write_run_manifest(run.out, roles=None, tier="protein_list", modules=run.modules)
    cards: dict[str, Any] = _cards(run)
    assert "proteins_placed" in cards, cards
    assert not any(ZERO_OF_ZERO.search(t) for t in cards.values()), cards
