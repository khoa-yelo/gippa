"""M12-04 PR 10: out-of-gene start artefacts get ``putative_start`` and are capped at ``candidate``
(docs/specs/M12-04-outputs-revision.md "PR 10" Rule 1-5, Edge cases, AC 10.1; statistician
sign-off ``local/handoff/pr10_statistician.md`` C1, C2, C4, C5).

Interface assumed (test author): the rule runs inside the ``variants`` chain
(``gippa.steps.variants.run_variants`` / ``run_variants_qc``, signatures unchanged), so every test
here drives the chain on a hand-built world (``helpers.putative_start_toys``). A flagged row has
``putative_start`` in ``flags`` (§9.4 order of ``gippa.variants.flags.FLAG_ORDER``, after
``database_derived``) and ``confidence = candidate``; ``n_ref``, ``n_alt``, ``call`` and the class
error columns are those of the uncapped chain. Every world also carries a positive control
(:data:`CONTROL`, an out-of-gene GTG read as M at peptide position 0 by two peptides), so a
"not flagged" case cannot pass on an implementation that flags nothing.

Not covered: a ``block`` whose ALT contains M cannot occur (no §7.2.4 block or near-block holds
M), so only the swap form of "a block or swap containing M" is tested.
"""

from __future__ import annotations

from typing import Any

import pandas as pd
import pytest
from helpers.putative_start_toys import (
    Site,
    World,
    cfg_of,
    flags,
    peptide_of,
    row_at,
    run_world,
)

FLAG = "putative_start"
CONTROL = Site(100)  # out-of-gene GTG after G, M at peptide position 0 in both peptides
GENE = ("G1", 20, 90)  # in-gene sites sit at codon 60


def with_control(**kw: Any) -> World:
    sites = (*kw.pop("sites"), CONTROL)
    return World(sites=sites, **kw)


def assert_flagged(row: dict[str, Any], *, n_alt: int = 2, call: str = "ALT") -> None:
    assert FLAG in flags(row), row
    assert row["confidence"] == "candidate", row
    assert (int(row["n_alt"]), row["call"]) == (n_alt, call), row


def assert_not_flagged(row: dict[str, Any], confidence: str = "high") -> None:
    assert FLAG not in flags(row), row
    assert row["confidence"] == confidence, row


def run_case(w: World, s: int, **cfg: Any) -> tuple[dict[str, Any], dict[str, Any]]:
    """(case row at codon ``s``, control row)."""
    table, _, _ = run_world(w, cfg_of(table=w.table, **cfg))
    return row_at(table, w, s), row_at(table, w, CONTROL.s)


# ---------------------------------------------------------------- AC 10.1, minus strand


def test_ac_10_1_minus_out_of_gene_gtg_two_n_terminal_peptides() -> None:
    w = World(strand="-", sites=(Site(60),))
    table, _, _ = run_world(w)
    assert_flagged(row_at(table, w, 60))


def test_ac_10_1_same_site_with_m_at_peptide_position_3_is_not_flagged() -> None:
    w = with_control(strand="-", sites=(Site(60, peptides=((3, 9), (3, 12))),))
    case, ctrl = run_case(w, 60)
    assert_not_flagged(case)
    assert_flagged(ctrl)


def test_var_107_minus_ttg_flags_in_section_9_4_order() -> None:
    """Spec example VAR_107 (no gene, ``TTG L→M``, high before): flags
    ``snv;strain_difference;database_derived;putative_start`` after, then M12-04 PR 11c's
    ``class_error_unestimated`` (this run has no spectra), which ends the §9.4 order."""
    w = World(strand="-", sites=(Site(60, ref_codon="TTG"),))
    table, _, _ = run_world(w)
    row = row_at(table, w, 60)
    assert flags(row) == (
        "snv", "strain_difference", "database_derived", FLAG, "class_error_unestimated"
    )  # fmt: skip
    assert row["confidence"] == "candidate"


def test_minus_strand_m_at_the_last_peptide_position_is_not_flagged() -> None:
    """Peptide order is transcription order: the lowest coordinate is the last residue."""
    w = with_control(strand="-", sites=(Site(60, peptides=((7, 8), (10, 11))),))
    case, ctrl = run_case(w, 60)
    assert_not_flagged(case)
    assert_flagged(ctrl)


def test_plus_strand_m_at_the_last_peptide_position_is_not_flagged() -> None:
    w = with_control(sites=(Site(60, peptides=((7, 8), (10, 11))),))
    case, ctrl = run_case(w, 60)
    assert_not_flagged(case)
    assert_flagged(ctrl)


# ---------------------------------------------------------------- condition 5: trypsin guard


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize(("prev", "flagged"), [("K", False), ("R", False), ("G", True)])
def test_in_gene_guard_on_the_preceding_residue(strand: str, prev: str, flagged: bool) -> None:
    w = with_control(strand=strand, genes=(GENE,), sites=(Site(60, prev=prev),))
    case, ctrl = run_case(w, 60)
    assert case["gene"] == "G1"
    if flagged:
        assert_flagged(case)
    else:
        assert_not_flagged(case)
    assert_flagged(ctrl)


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("prev", ["K", "R"])
def test_out_of_gene_after_k_or_r_is_flagged(strand: str, prev: str) -> None:
    w = World(strand=strand, sites=(Site(60, prev=prev),))
    table, _, _ = run_world(w)
    row = row_at(table, w, 60)
    assert pd.isna(row["gene"])
    assert_flagged(row)


def test_in_gene_with_no_preceding_codon_on_a_linear_contig_counts_as_not_k_r() -> None:
    """Rule 5: no preceding codon on a linear contig is "not K/R". A pseudogene frames the codon
    (``gene`` set) but has no annotated start (``flags.start_codons`` skips pseudogenes)."""
    w = with_control(genes=(("P1", 0, 40, "pseudogene"),), sites=(Site(0, prev=None),))
    case, ctrl = run_case(w, 0)
    assert case["gene"] == "P1"
    assert_flagged(case)
    assert_flagged(ctrl)


# ---------------------------------------------------------------- condition 4: every peptide first


def test_mixed_peptide_positions_are_not_flagged() -> None:
    w = with_control(sites=(Site(60, peptides=((0, 8), (2, 11))),))
    case, ctrl = run_case(w, 60)
    assert_not_flagged(case)
    assert_flagged(ctrl)


def test_reference_m_before_the_site_is_not_an_alt_first_residue() -> None:
    """A matched reference M at peptide position 0 with the ALT M at 1: not flagged (Edge
    cases: "M at position 1 but the site at position 2 or later")."""
    w = with_control(sites=(Site(60, prev="M", peptides=((1, 8), (1, 11))),))
    case, ctrl = run_case(w, 60)
    assert_not_flagged(case)
    assert_flagged(ctrl)


def test_one_peptide_with_cap_from_high_min_peptides_1() -> None:
    """``high_min_peptides = 1``: one N-terminal peptide would be ``high``; capped."""
    w = World(sites=(Site(60, peptides=((0, 9),)),))
    table, _, _ = run_world(w, cfg_of(high_min_peptides=1))
    assert_flagged(row_at(table, w, 60), n_alt=1)


def test_the_rule_applies_whatever_the_call() -> None:
    """A REF peptide at the codon makes the call ``mixed``; the row is still flagged and capped,
    with ``n_ref``, ``n_alt`` and ``call`` unchanged."""
    w = World(sites=(Site(60, ref_peptides=((3, 10),)),))
    table, _, _ = run_world(w)
    row = row_at(table, w, 60)
    assert_flagged(row, call="mixed")
    assert int(row["n_ref"]) == 1


# ---------------------------------------------------------------- conditions 1-2: the start set


def test_ctg_with_the_default_start_codons_is_not_flagged() -> None:
    w = with_control(sites=(Site(60, ref_codon="CTG"),))
    case, ctrl = run_case(w, 60)
    assert_not_flagged(case)
    assert_flagged(ctrl)


def test_ctg_is_flagged_when_configured() -> None:
    w = World(sites=(Site(60, ref_codon="CTG"),))
    table, _, _ = run_world(w, cfg_of(start_codons=["ATG", "GTG", "TTG", "CTG"]))
    assert_flagged(row_at(table, w, 60))


def test_gtg_is_not_flagged_when_dropped_from_start_codons() -> None:
    """The start set is the run's ``annotation.start_codons``, never a constant."""
    w = with_control(sites=(Site(60, ref_codon="TTG"),))
    table, _, _ = run_world(w, cfg_of(start_codons=["ATG", "TTG"]))
    assert_flagged(row_at(table, w, 60))
    assert_not_flagged(row_at(table, w, CONTROL.s))  # the control is GTG


def test_table_4_tta_is_flagged_when_configured() -> None:
    w = World(table=4, sites=(Site(60, ref_codon="TTA"),))
    table, _, _ = run_world(w, cfg_of(table=4, start_codons=["ATG", "GTG", "TTG", "TTA"]))
    assert_flagged(row_at(table, w, 60))


def test_table_4_tta_with_the_default_start_codons_is_not_flagged() -> None:
    w = with_control(table=4, sites=(Site(60, ref_codon="TTA"),))
    case, ctrl = run_case(w, 60)
    assert_not_flagged(case)
    assert_flagged(ctrl)


def test_table_25_gtg_is_flagged() -> None:
    w = World(table=25, sites=(Site(60),))
    table, _, _ = run_world(w, cfg_of(table=25))
    assert_flagged(row_at(table, w, 60))


def test_k_to_m_at_aag_is_not_flagged() -> None:
    w = with_control(sites=(Site(60, ref_codon="AAG"),))
    case, ctrl = run_case(w, 60)
    assert case["ref_aa"] == "K"
    assert_not_flagged(case)
    assert_flagged(ctrl)


def test_swap_containing_m_is_not_flagged() -> None:
    """A swap over ``VM`` reads the site's ATG as V; the M side of the swap is ATG. Neither the
    swap row nor any row at the ATG codon is flagged."""
    w = with_control(sites=(Site(60, kind="swap"),))
    table, _, _ = run_world(w, cfg_of())
    rows = table[table["codon_start"].astype(int).isin([w.codon_start(60), w.codon_start(61)])]
    assert len(rows) == 1 and rows.iloc[0]["alt_type"] == "swap", rows
    assert FLAG not in flags(rows.iloc[0].to_dict())
    assert rows.iloc[0]["confidence"] == "high"
    assert_flagged(row_at(table, w, CONTROL.s))


# ---------------------------------------------------------------- condition 3 and the origin


@pytest.mark.parametrize("at", [1, 2])
def test_plus_codon_across_the_origin_is_flagged(at: int) -> None:
    """``+`` codon at L - 1 or L - 2 of a circular contig (spans the origin)."""
    w0 = World(topology="circular", sites=(Site(40),))
    w = World(topology="circular", off=w0.length - at - 3 * 40, sites=(Site(40), CONTROL))
    assert w.codon_start(40) == w.length - at
    table, _, _ = run_world(w)
    assert_flagged(row_at(table, w, 40))
    assert_flagged(row_at(table, w, CONTROL.s))


@pytest.mark.parametrize(("prev", "flagged"), [("K", False), ("G", True)])
def test_preceding_residue_wraps_the_origin(prev: str, flagged: bool) -> None:
    """In-gene ``+`` codon at 0: the preceding codon is at L - 3."""
    n = 120
    w = World(topology="circular", off=3 * n - 3 * 40, genes=(("G1", 20, 80),),
              sites=(Site(40, prev=prev), CONTROL))  # fmt: skip
    assert w.codon_start(40) == 0 and w.codon_start(39) == w.length - 3
    table, _, _ = run_world(w)
    case = row_at(table, w, 40)
    assert case["gene"] == "G1"
    if flagged:
        assert_flagged(case)
    else:
        assert_not_flagged(case)
    assert_flagged(row_at(table, w, CONTROL.s))


def test_minus_peptides_unwrapped_past_the_origin_are_matched_mod_l() -> None:
    """``-`` site at codon_start 3: its peptides run through 0, so their placements are unwrapped
    (``end`` > L) and the alignment codon is ``L + 3``; the site is still flagged."""
    n = 120
    w = World(strand="-", topology="circular", off=6 + 3 * 50, sites=(Site(50), CONTROL),
              n_codons=n)  # fmt: skip
    assert w.codon_start(50) == 3
    table, _, _ = run_world(w)
    assert_flagged(row_at(table, w, 50))


def test_minus_peptide_with_a_later_m_unwrapped_past_the_origin_vetoes() -> None:
    """One N-terminal peptide inside [0, L) and one with the M at position 4 whose placement
    crosses the origin (codon ``L + c``): not every peptide starts at the M, so not flagged."""
    n = 120
    w = World(strand="-", topology="circular", off=30 + 3 + 3 * 50,
              sites=(Site(50, peptides=((0, 8), (4, 20))), CONTROL), n_codons=n)  # fmt: skip
    assert w.codon_start(50) == 30
    lo, hi = w.span(46, 66)
    assert hi > w.length, (lo, hi)  # the q = 4 placement is unwrapped past the origin
    table, _, _ = run_world(w)
    assert_not_flagged(row_at(table, w, 50))
    assert_flagged(row_at(table, w, CONTROL.s))


def test_annotated_start_across_the_origin_is_start_codon_not_putative_start() -> None:
    n = 120
    w = World(topology="circular", off=3 * n - 2 - 3 * 40, genes=(("G1", 40, 80),),
              sites=(Site(40), CONTROL))  # fmt: skip
    assert w.codon_start(40) == w.length - 2
    table, _, _ = run_world(w)
    case = row_at(table, w, 40)
    assert (case["gene"], int(case["aa_pos"])) == ("G1", 1)
    assert "start_codon" in flags(case)
    assert_not_flagged(case)
    assert_flagged(row_at(table, w, CONTROL.s))


@pytest.mark.parametrize("strand", ["+", "-"])
def test_annotated_non_atg_start_is_start_codon_not_putative_start(strand: str) -> None:
    w = with_control(strand=strand, genes=(("G1", 60, 90),), sites=(Site(60),))
    case, ctrl = run_case(w, 60)
    assert "start_codon" in flags(case)
    assert_not_flagged(case)
    assert_flagged(ctrl)


# ---------------------------------------------------------------- C1: id_mixups


def test_c1_same_db_peptide_on_a_putative_start_site_is_an_id_mixup() -> None:
    from gippa.steps.variants import qc_frame
    from gippa.variants.qc import EXEMPT_FLAGS

    assert frozenset({"start_codon", "stop_codon"}) == EXEMPT_FLAGS
    w = World(sites=(Site(60, sources=(("database", "same"), ("denovo", "same"))),))
    table, _, support = run_world(w, cfg_of("same"))
    row = row_at(table, w, 60)
    assert_flagged(row)
    qc = qc_frame(table, support)
    got = qc.set_index("metric").loc["id_mixups"]
    assert (int(got["value"]), int(got["n"])) == (1, 1)


# ---------------------------------------------------------------- C4: class error rates


def _spectra(peptides: list[str], errors: set[str]) -> pd.DataFrame:
    """Two de novo spectra per peptide; the reference reads the same peptide, except the second
    spectrum of each peptide in ``errors`` (one error each)."""
    rows = []
    for i, p in enumerate(peptides):
        for j in range(2):
            sid = f"s{i}_{j}"
            ref = "WWWWWWW" if (p in errors and j == 1) else p
            rows.append({"run": "r1", "spectrum_id": sid, "peptide": p, "source": "denovo",
                         "q_value": 0.001, "is_decoy": False})  # fmt: skip
            rows.append({"run": "r1", "spectrum_id": sid, "peptide": ref, "source": "database",
                         "q_value": 0.001, "is_decoy": False})  # fmt: skip
    return pd.DataFrame(rows)


def test_c4_class_error_is_pooled_over_capped_and_uncapped_rows() -> None:
    """The cap changes no class rate: a ``putative_start`` row and an unflagged row of the same
    class ``(substitution, mass_equivalent)`` carry the same rate and n, counted over the paired
    reads of both rows' peptides (8 reads, 3 errors)."""
    capped = Site(30, sources=(("denovo", "same"),) * 2)
    plain = Site(60, prev="K", sources=(("denovo", "same"),) * 2)  # in-gene after K: no flag
    w = World(genes=(("G1", 45, 90),), sites=(capped, plain))
    peps = [peptide_of(w, st, q, n) for st in (capped, plain) for q, n in st.peptides]
    spectra = _spectra(peps, {peps[0], peps[2], peps[3]})
    table, _, _ = run_world(w, cfg_of("same"), spectra)
    a, b = row_at(table, w, 30), row_at(table, w, 60)
    assert_flagged(a)
    assert_not_flagged(b)
    assert "mass_equivalent" in flags(a) and "mass_equivalent" in flags(b)
    for r in (a, b):
        assert int(r["class_error_n"]) == 8, r
        assert float(r["class_error_rate"]) == pytest.approx(3 / 8), r


# ---------------------------------------------------------------- C5: VCF CONF = TSV confidence


def test_c5_vcf_conf_equals_tsv_confidence_on_every_row() -> None:
    w = World(genes=(GENE,), sites=(
        Site(10), Site(30, peptides=((3, 9), (3, 12))), Site(60, prev="K"), Site(75),
        Site(100, ref_peptides=((2, 10),)),
    ))  # fmt: skip
    table, records, _ = run_world(w)
    by_pos = {(r.contig, r.pos0): r for r in records}
    n_capped = 0
    for row in table.to_dict("records"):
        if "vcf_unrepresentable" in flags(row):
            continue
        rec = by_pos[(row["contig"], int(row["codon_start"]))]
        i = list(rec.info["AA_ALT"]).index(row["alt_aa"])
        assert rec.info["CONF"][i] == row["confidence"], (row, rec.info)
        n_capped += FLAG in flags(row)
        if FLAG in flags(row):
            assert FLAG in rec.info["FLAGS"], rec.info
    assert n_capped == 3, table[["codon_start", "confidence", "flags"]]
