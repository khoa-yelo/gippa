"""ADR-0078: every q-passing decoy with protein IDs competes in the picked gene FDR.

Spec: docs/adr/ADR-0078-decoy-pairing-independent-of-acceptance.md ("Spec for tests", examples
1-11), docs/methods/detection-tiers.md; GIPPA_SPEC.md §7.3.2 (decision 4), §8.1, §8.5. Keeps
ADR-0047, M4-12 decision 12 and M7-12 decisions 23/28/29; retires ADR-0009's 1/f.

Notation (ADR "Spec for tests"): C2 genes G1-G4, non-partial; protein_hits empty; default
prefixes (``rev_``, ``CON__``). Key S = {q}: targets and decoys carry q only (PEP and score NA).
Accepted unique targets, one per gene, sit on proteins P1-P4 (voted to G1-G4) with q 0.001,
0.002, 0.003, 0.004. Every decoy passes q <= 0.01.

Rule pinned here:
- ``pair_decoys`` gains a ``unit`` column: the gene when every counted stripped ID maps to one and
  the same gene; otherwise the decoy's own unit, keyed by its sorted, de-duplicated counted
  stripped IDs joined by ``;``. Decoys with one key share one unit (its best decoy counts). Own
  units never take a target and contribute to D only.
- Picked q = (D(t) + 1)/T(t), made monotone, over genes + own units (no 1/f).
- Summary: ``n_decoys_paired``, ``n_decoys_unpaired`` and ``decoy_paired_share`` as in 1.1.0;
  ``decoy_win_weight`` 1.0 whenever there is a pairing row.
- S = ``shared_key_columns(accepted, <all pairing rows>)`` in ``detect`` and ``choose_tier``.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import detection_tiers_fixtures as fx
from detection_tiers_fixtures import Builder, cfg

GENES = ["G1", "G2", "G3", "G4"]
TARGET_Q = [0.001, 0.002, 0.003, 0.004]
EXT = "GIPPAORF_nterm_extension_chr1_50_190_p"
INTERGENIC = "GIPPAORF_intergenic_c1_100_200_p"
TWO_THIRDS = 2.0 / 3.0


def _guard(share: float = 0.0) -> Any:
    """Config whose decoy-pairing guard is ``share`` (max of min_identity / min_coverage)."""
    return cfg(mapping={"proteins": {"min_identity": share, "min_coverage": share}})


def _builder(tags: list[str] | None = None, partial: set[str] | None = None) -> Builder:
    """G1-G4 (or ``tags``) with the ADR's four accepted targets on P1-P4 (q only)."""
    b = Builder(tags=list(tags or GENES), partial=set(partial or ()))
    for i, (g, q) in enumerate(zip(GENES, TARGET_Q, strict=True)):
        b.target(g, q, None, None, protein_ids=[f"P{i + 1}"])
    return b


def _pairing(r: fx.Run, c: Any = None) -> pd.DataFrame:
    """``pair_decoys`` as ``detect`` calls it (vote map over the accepted rows, q filter on)."""
    from gippa.genetable.detection import pair_decoys, protein_gene_map

    c = c if c is not None else cfg()
    pmap = protein_gene_map(r.genes, r.peptides, r.peptide_hits, r.protein_hits,
                            accepted=r.accepted,
                            contaminant_prefixes=c.input.contaminant_prefixes)  # fmt: skip
    return pair_decoys(r.peptides, pmap, c.input.decoy_prefixes, q_value_max=c.input.q_value_max,
                       contaminant_prefixes=c.input.contaminant_prefixes)  # fmt: skip


def _gene_unit(pairing: pd.DataFrame) -> list[tuple[str | None, str]]:
    """(gene or None, unit) per row, sorted, so the row order does not matter."""
    assert "unit" in pairing.columns, f"pair_decoys lacks the ADR-0078 unit column: {list(pairing)}"
    rows = [
        (None if fx.is_na(g) else str(g), str(u))
        for g, u in zip(pairing["gene"].tolist(), pairing["unit"].tolist(), strict=True)
    ]
    return sorted(rows, key=lambda x: (x[0] or "", x[1]))


def _qs(r: fx.Run, tags: list[str] | None = None) -> list[float]:
    return [r.q(g) for g in (tags or GENES)]


def _assert_q(r: fx.Run, expected: list[float], tags: list[str] | None = None) -> None:
    got = _qs(r, tags)
    for g, q, e in zip(tags or GENES, got, expected, strict=True):
        if math.isnan(e):
            assert math.isnan(q), (g, q)
        else:
            assert q == pytest.approx(e, rel=1e-9), (g, got, expected)


# ---------------------------------------------------------------- example 1 (the defect)


def _ex1() -> Builder:
    b = _builder()
    b.decoy(["rev_P2"], 0.0025, None, None)
    b.decoy(["rev_P9"], 0.0015, None, None)  # P9 carries no target
    return b


def test_ex1_pairing_rows_gene_and_own_unit() -> None:
    r = _ex1().run(_guard(0.5))
    assert _gene_unit(_pairing(r)) == [(None, "P9"), ("G2", "G2")]


def test_ex1_targetless_decoy_counts_as_a_decoy_win() -> None:
    """Units G1 T .001, P9 D .0015, G2 T .002, G3 T .003, G4 T .004 -> q 0.5 for all
    (1.1.0: f 0.5, D 0 -> 0.25)."""
    r = _ex1().run(_guard(0.5))
    assert r.det.tier == "decoy_picked"
    _assert_q(r, [0.5] * 4)


def test_ex1_summary_rows() -> None:
    r = _ex1().run(_guard(0.5))
    assert fx.tier_metric(r.det, "n_decoys_paired") == 1
    assert fx.tier_metric(r.det, "n_decoys_unpaired") == 1
    assert fx.tier_metric(r.det, "decoy_paired_share") == pytest.approx(0.5)
    assert fx.tier_metric(r.det, "decoy_win_weight") == 1.0, "1/f is retired (ADR-0078 rule 5)"


def test_ex1_default_guard_still_falls_to_mayu() -> None:
    """Rule 4: the 80 % guard is unchanged (share 0.5 -> mayu, ``decoy_pairing_low``)."""
    r = _ex1().run()
    assert r.det.tier == "mayu"
    assert "decoy_pairing_low" in r.det.warnings


# ---------------------------------------------------------------- example 2 (multi-gene, partial)


def _ex2() -> Builder:
    b = _builder()
    b.decoy(["rev_P2"], 0.0025, None, None)
    b.decoy(["rev_P1;rev_P2"], 0.0015, None, None)
    b.decoy(["rev_P7", "rev_P1"], 0.0035, None, None)  # P7 unmapped
    return b


def test_ex2_unit_keys_are_sorted_counted_ids() -> None:
    r = _ex2().run(_guard())
    assert _gene_unit(_pairing(r)) == [(None, "P1;P2"), (None, "P1;P7"), ("G2", "G2")]


def test_ex2_multi_gene_and_partly_mapped_decoys_compete() -> None:
    """G1 T, P1;P2 D, G2 T, G3 T, P1;P7 D, G4 T -> G1-G3 2/3, G4 0.75 (1.1.0: 0.25)."""
    r = _ex2().run(_guard())
    assert r.det.tier == "decoy_picked"
    _assert_q(r, [TWO_THIRDS, TWO_THIRDS, TWO_THIRDS, 0.75])


# ---------------------------------------------------------------- example 3 (same key twice)


def _ex3() -> Builder:
    b = _builder()
    b.decoy(["rev_P9"], 0.0015, None, None)
    b.decoy(["rev_P9"], 0.0030, None, None)
    return b


def test_ex3_same_key_is_one_unit() -> None:
    r = _ex3().run(_guard())
    assert _gene_unit(_pairing(r)) == [(None, "P9"), (None, "P9")]


def test_ex3_only_the_best_decoy_of_a_unit_counts() -> None:
    r = _ex3().run(_guard())
    assert r.det.tier == "decoy_picked"
    _assert_q(r, [0.5] * 4)


# ---------------------------------------------------------------- example 4 (contaminants)


def _ex4() -> Builder:
    b = _builder()
    b.decoy(["rev_CON__X"], 0.0005, None, None)
    b.decoy(["rev_P9", "rev_CON__X"], 0.0015, None, None)
    b.decoy(["rev_P2", "rev_CON__X"], 0.0025, None, None)
    return b


def test_ex4_contaminant_decoys_unitised_on_non_contaminant_ids() -> None:
    r = _ex4().run(_guard())
    assert _gene_unit(_pairing(r)) == [(None, "P9"), ("G2", "G2")]


def test_ex4_mixed_contaminant_decoy_competes_in_its_own_unit() -> None:
    r = _ex4().run(_guard())
    assert fx.tier_metric(r.det, "n_decoys_paired") == 1
    assert fx.tier_metric(r.det, "n_decoys_unpaired") == 1
    _assert_q(r, [0.5] * 4)


# ---------------------------------------------------------------- example 5 (ORF entries)


def test_ex5_orf_entries_unchanged() -> None:
    b = _builder()
    b.target("G1", 0.001, None, None, protein_ids=[EXT])  # votes the extension entry to G1
    b.decoy([f"rev_{INTERGENIC}"], 0.0005, None, None)
    b.decoy([f"rev_{EXT}"], 0.0025, None, None)
    r = b.run(_guard())
    assert _gene_unit(_pairing(r)) == [("G1", "G1")]
    assert fx.tier_metric(r.det, "n_decoys_paired") == 1
    assert fx.tier_metric(r.det, "n_decoys_unpaired") == 0


# ---------------------------------------------------------------- example 6 (locus tags)


def test_ex6_locus_tag_run_unchanged() -> None:
    """Regression guard (passes on 1.1.0): G3 decoy win, G1 G2 G4 at 2/3."""
    b = Builder(tags=list(GENES))
    for g, q in zip(GENES, TARGET_Q, strict=True):
        b.target(g, q, None, None, protein_ids=[g])
    b.decoy(["rev_G2"], 0.0025, None, None)
    b.decoy(["rev_G3"], 0.0005, None, None)
    r = b.run()
    assert r.det.tier == "decoy_picked"
    _assert_q(r, [TWO_THIRDS, TWO_THIRDS, math.nan, TWO_THIRDS])
    assert fx.tier_metric(r.det, "n_decoys_paired") == 2
    assert fx.tier_metric(r.det, "decoy_win_weight") == 1.0


# ---------------------------------------------------------------- example 7 (namespace)


def _ex7(decoy_q: float = 0.0035) -> Builder:
    tags = ["G1", "G2", "P1;P2", "G4"]
    b = Builder(tags=tags)
    for i, (g, q) in enumerate(zip(tags, TARGET_Q, strict=True)):
        b.target(g, q, None, None, protein_ids=[f"P{i + 1}"])
    b.decoy(["rev_P1;rev_P2"], decoy_q, None, None)  # default: worse than gene "P1;P2"'s target
    return b


def test_ex7_own_unit_never_takes_a_target() -> None:
    """Own unit ``P1;P2`` (D .0035) stays apart from gene ``P1;P2`` (T .003): FDR 1, 1/2, 1/3, 2/4
    -> q 1/3, 1/3, 1/3, 0.5. Merged with the gene it would lose: q 0.25 for all."""
    r = _ex7().run(_guard())
    assert r.det.tier == "decoy_picked"
    tags = ["G1", "G2", "P1;P2", "G4"]
    assert _gene_unit(_pairing(r)) == [(None, "P1;P2")]
    _assert_q(r, [1 / 3, 1 / 3, 1 / 3, 0.5], tags)


def test_ex7_own_unit_beating_the_namesake_target_stays_apart() -> None:
    """Own unit ``P1;P2`` (D .002) beats gene ``P1;P2``'s target (.003) but is not its decoy: the
    gene keeps a finite q and the own unit adds one to D. FDR 1, 2/2, 2/3, 2/4 -> q 0.5 for all.
    Code review of PR #607: a mutant grouping own units with genes by key
    survived the test above.
    """
    r = _ex7(0.002).run(_guard())
    assert r.det.tier == "decoy_picked"
    _assert_q(r, [0.5, 0.5, 0.5, 0.5], ["G1", "G2", "P1;P2", "G4"])


# ---------------------------------------------------------------- example 8 (ADR-0047)


def test_ex8_decoy_above_q_value_max_has_no_row_when_targets_carry_q() -> None:
    b = _builder()
    b.decoy(["rev_P1"], 0.009, None, None)
    b.decoy(["rev_P9"], 0.02, None, None)
    r = b.run()
    assert fx.tier_metric(r.det, "n_decoys_paired") == 1
    assert fx.tier_metric(r.det, "n_decoys_unpaired") == 0


def test_ex8_decoy_above_q_value_max_has_a_row_when_no_target_carries_q() -> None:
    b = Builder(tags=list(GENES))
    for i, g in enumerate(GENES):
        b.target(g, None, 0.001 * (i + 1), None, protein_ids=[f"P{i + 1}"])
    b.decoy(["rev_P1"], 0.009, 0.5, None)
    b.decoy(["rev_P9"], 0.02, 0.5, None)
    r = b.run()
    assert fx.tier_metric(r.det, "n_decoys_paired") == 1
    assert fx.tier_metric(r.det, "n_decoys_unpaired") == 1


# ---------------------------------------------------------------- example 9 (random loss)


def _ex9() -> Builder:
    """Four losing gene decoys keep the share at 0.8 (default guard) next to ``rev_LOST_0``."""
    b = _builder()
    for i in range(4):
        b.decoy([f"rev_P{i + 1}"], 0.009, None, None)
    b.decoy(["rev_LOST_0"], 0.0015, None, None)
    return b


def test_ex9_lost_decoy_gets_its_own_unit() -> None:
    r = _ex9().run()
    assert (None, "LOST_0") in _gene_unit(_pairing(r))


def test_ex9_lost_decoy_counts_as_a_decoy_win() -> None:
    """G1 T, LOST_0 D, G2-G4 T -> q 0.5 for all (1.1.0: f 0.8, D 0 -> 0.25)."""
    r = _ex9().run()
    assert r.det.tier == "decoy_picked"
    _assert_q(r, [0.5] * 4)
    assert fx.tier_metric(r.det, "decoy_paired_share") == pytest.approx(0.8)
    assert fx.tier_metric(r.det, "n_decoys_unpaired") == 1
    assert fx.tier_metric(r.det, "decoy_win_weight") == 1.0


# ---------------------------------------------------------------- example 10 (ranking key)


def _ex10() -> Builder:
    b = Builder(tags=list(GENES))
    for i, (g, q) in enumerate(zip(GENES, TARGET_Q, strict=True)):
        b.target(g, q, None, 50.0 - 10 * i, protein_ids=[f"P{i + 1}"])
    b.decoy(["rev_P1"], 0.001, None, 10.0)  # ties G1 on q, loses on score
    for i in range(1, 4):
        b.decoy([f"rev_P{i + 1}"], 0.009, None, 1.0)
    b.decoy(["rev_P9"], 0.0015, None, None)  # own unit, score NA
    return b


def test_ex10_own_unit_na_score_drops_score_from_s() -> None:
    from gippa.genetable.detection import choose_tier

    r = _ex10().run()
    tier, warns = choose_tier(r.peptides, r.peptide_hits, r.protein_hits, _pairing(r),
                              accepted=r.accepted)  # fmt: skip
    assert tier == "decoy_picked"
    assert "decoy_unscored" in warns
    assert r.det.tier == "decoy_picked"
    assert "decoy_unscored" in r.det.warnings


def test_ex10_ranking_on_q_only() -> None:
    """S = {q}: the G1 tie goes to the decoy; with P9 two decoy wins lie above G2-G4 -> q 1."""
    r = _ex10().run()
    _assert_q(r, [math.nan, 1.0, 1.0, 1.0])


# ---------------------------------------------------------------- example 11 (partial gene)


def _ex11() -> Builder:
    b = _builder(tags=[*GENES, "G5"], partial={"G5"})
    b.decoy(["rev_G5"], 0.0005, None, None)  # pairs to the partial gene by locus tag
    return b


def test_ex11_decoy_of_a_partial_gene_leaves_the_competition() -> None:
    """Regression guard (passes on 1.1.0): D 0 -> q 0.25 for G1-G4; G5 partial_excluded."""
    r = _ex11().run()
    assert r.det.tier == "decoy_picked"
    _assert_q(r, [0.25] * 4)
    assert r.status("G5") == "partial_excluded"
    assert fx.tier_metric(r.det, "n_decoys_paired") == 1


def test_ex11_partial_gene_decoy_unit_is_the_gene() -> None:
    r = _ex11().run()
    assert _gene_unit(_pairing(r)) == [("G5", "G5")]


# ---------------------------------------------------------------- properties (oracle)

POOL = ["P1", "P2", "P3", "P4", "P7", "P8", "P9", "CON__X"]
Q_GRID = [0.0005, 0.001, 0.0015, 0.002, 0.0025, 0.003, 0.0035, 0.004, 0.006, 0.009]
VOTED = {f"P{i + 1}": g for i, g in enumerate(GENES)}


def _oracle(
    decoys: list[tuple[list[str], float]],
) -> tuple[list[tuple[str | None, str]], list[float]]:
    """ADR-0078 rules 1-2 from scratch: rows (gene, unit) and gene q over genes + own units."""
    from gippa.stats.tda import picked_fdr

    rows: list[tuple[str | None, str]] = []
    best: dict[tuple[str, str], float] = {}
    for ids, q in decoys:
        counted = sorted({x.removeprefix("rev_") for x in ids} - {"CON__X"})
        if not counted:
            continue
        hit = {VOTED.get(x) for x in counted}
        gene = next(iter(hit)) if len(hit) == 1 and None not in hit else None
        unit = ("gene", gene) if gene is not None else ("own", ";".join(counted))
        rows.append((gene, gene if gene is not None else unit[1]))
        best[unit] = max(best.get(unit, -np.inf), -q)
    own = sorted(u for u in best if u[0] == "own")
    t = [-q for q in TARGET_Q] + [np.nan] * len(own)
    d = [best.get(("gene", g), np.nan) for g in GENES] + [best[u] for u in own]
    q = picked_fdr(np.asarray(t), np.asarray(d)).q[: len(GENES)]
    return sorted(rows, key=lambda x: (x[0] or "", x[1])), [float(v) for v in q]


@settings(max_examples=20, deadline=None)
@given(
    decoys=st.lists(
        st.tuples(
            st.lists(st.sampled_from(POOL), min_size=1, max_size=3, unique=True),
            st.sampled_from(Q_GRID),
        ),
        min_size=1,
        max_size=6,
    )
)
def test_property_units_and_q_match_the_adr0078_oracle(
    decoys: list[tuple[list[str], float]],
) -> None:
    """Any mix of gene-paired, multi-gene, unmapped and contaminant decoys: the unit column and
    the gene q follow rules 1-2 (best decoy per unit, ties to the decoy, no 1/f)."""
    b = _builder()
    for ids, q in decoys:
        b.decoy([f"rev_{x}" for x in ids], q, None, None)
    rows, q_exp = _oracle(decoys)
    r = b.run(_guard())
    assert _gene_unit(_pairing(r)) == rows
    if not rows:
        return
    assert r.det.tier == "decoy_picked"
    _assert_q(r, q_exp)


@settings(max_examples=15, deadline=None)
@given(
    base=st.lists(st.sampled_from(Q_GRID), min_size=0, max_size=3),
    extra_q=st.sampled_from(Q_GRID),
    name=st.sampled_from(["P7", "P8", "LOST_3", "P1;P9"]),
)
def test_property_an_extra_unmappable_decoy_never_lowers_gene_q(
    base: list[float], extra_q: float, name: str
) -> None:
    """Rule 3: an unmappable decoy can only add to D, so no gene q falls (counting it as a win
    is conservative)."""
    b0, b1 = _builder(), _builder()
    for i, q in enumerate(base):
        for b in (b0, b1):
            b.decoy([f"rev_P{i + 1}"], q, None, None)
    b1.decoy([f"rev_{x}" for x in name.split(";")], extra_q, None, None)
    q0 = np.array(_qs(b0.run(_guard())))
    r1 = b1.run(_guard())
    q1 = np.array(_qs(r1))
    assert r1.det.tier == "decoy_picked"
    both = np.isfinite(q0) & np.isfinite(q1)
    assert (np.isfinite(q0) == np.isfinite(q1)).all(), "an own unit never takes a gene"
    assert (q1[both] >= q0[both] - 1e-12).all(), (q0, q1)
