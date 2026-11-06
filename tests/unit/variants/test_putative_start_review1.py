"""M12-04 PR 10 review round 1 (#566): test gaps in the ``putative_start`` rule
(docs/specs/M12-04-outputs-revision.md "PR 10" Rule 1 and 4; sign-off C2).

- Rule 4 reads **every** alignment of every ALT peptide at the codon: a peptide placed twice over
  the site codon (a tandem repeat), once with the M at peptide position 0 and once at position 3,
  vetoes the flag even when its q = 0 placement comes first; the same peptide's other placement
  at a different locus does not.
- Rule 1: the ALT must be M. A GTG read as A, or as I (``J``), by N-terminal peptides is not
  flagged.

Worlds are ``helpers.putative_start_toys``; each "not flagged" world carries the positive
:data:`CONTROL` site so a build that flags nothing cannot pass.
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
    gene_model_of,
    genome_of,
    hits_of,
    row_at,
    run_world,
)
from helpers.readthrough_toys import CODON

FLAG = "putative_start"
CONTROL = Site(100)
S = 60  # the site codon (transcription index, + strand)
#: ``MST`` repeats around the site: codon S - 3 is ATG, S is the GTG site, S + 3 and S + 6 ATG
REPEAT = {S - 3: "ATG", S - 2: CODON["S"], S - 1: CODON["T"], S + 1: CODON["S"],
          S + 2: CODON["T"], S + 3: "ATG", S + 4: CODON["S"], S + 5: CODON["T"],
          S + 6: "ATG", S + 7: CODON["S"], S + 8: CODON["T"]}  # fmt: skip
P = "MSTMSTMST"  # placed at [S, S + 9) (M at q 0) and at [S - 3, S + 6) (M at q 3)
Q = "MSTMSTMS"  # a second ALT peptide, placed once, M at q 0


def _hit(w: World, pep: str, a: int, n_diff: int, uniqueness: str = "multi_identical"
         ) -> dict[str, Any]:  # fmt: skip
    start, end = w.span(a, a + len(pep))
    return {
        "peptide": pep, "peptide_len_aa": len(pep), "source": "database", "contig": "c1",
        "start": start, "end": end, "strand": w.strand, "n_diff": n_diff,
        "diff_types": "substitution" if n_diff else pd.NA, "placement_q": 0.0,
        "uniqueness": uniqueness,
    }  # fmt: skip


def _run(w: World, rows: list[dict[str, Any]]) -> pd.DataFrame:
    from gippa.steps.variants import run_variants_qc

    ctrl_hits, ctrl_peps = hits_of(w)
    hits = pd.concat([pd.DataFrame(rows), ctrl_hits], ignore_index=True)
    own = pd.DataFrame({"peptide": sorted({r["peptide"] for r in rows}), "search_db": "other"})
    peps = pd.concat([own, ctrl_peps], ignore_index=True)
    _, table, _, _ = run_variants_qc(hits, None, peps, gene_model_of(w), {"c1": genome_of(w)},
                                     cfg_of())  # fmt: skip
    return table


def _world(extra: dict[int, str]) -> World:
    codons = dict(REPEAT)
    codons[S] = "GTG"
    codons.update(extra)
    return World(sites=(CONTROL,), extra_codons=codons)


def _assert_control(table: pd.DataFrame, w: World) -> None:
    ctrl = row_at(table, w, CONTROL.s)
    assert FLAG in flags(ctrl) and ctrl["confidence"] == "candidate", ctrl


@pytest.mark.parametrize("q0_first", [True, False])
def test_tandem_repeat_placement_with_a_later_m_vetoes(q0_first: bool) -> None:
    """``P`` is placed over the site codon twice: at q 0 and at q 3. Not every alignment at the
    codon has the M first, so no flag, whichever placement comes first in ``peptide_hits``."""
    w = _world({})
    rows = [_hit(w, P, S, 1), _hit(w, P, S - 3, 1)]
    if not q0_first:
        rows.reverse()
    rows.append(_hit(w, Q, S, 1, "unique"))
    table = _run(w, rows)
    row = row_at(table, w, S)
    assert int(row["n_alt"]) == 2, row
    assert FLAG not in flags(row), row
    assert row["confidence"] == "high", row
    _assert_control(table, w)


def test_the_same_peptide_placed_at_another_locus_does_not_veto() -> None:
    """``P`` at the site (q 0) and, exactly, at codons [20, 29): only alignments at this codon
    count, so the site is flagged."""
    other = {20 + i: ("ATG" if a == "M" else CODON[a]) for i, a in enumerate(P)}
    w = _world(other)
    rows = [_hit(w, P, 20, 0), _hit(w, P, S, 1), _hit(w, Q, S, 1, "unique")]
    table = _run(w, rows)
    row = row_at(table, w, S)
    assert int(row["n_alt"]) == 2, row
    assert FLAG in flags(row), row
    assert row["confidence"] == "candidate", row


@pytest.mark.parametrize(("alt", "row_alt"), [("A", "A"), ("I", "J")])
def test_alt_other_than_m_is_never_flagged(alt: str, row_alt: str) -> None:
    """Rule 1: an out-of-gene GTG read as A (or I, ``J`` in the table) by two N-terminal
    peptides is an ordinary substitution."""
    w = World(sites=(Site(S, alt=alt), CONTROL))
    table, _, _ = run_world(w)
    row = row_at(table, w, S, row_alt)
    assert (row["alt_type"], row["ref_codon"]) == ("substitution", "GTG"), row
    assert FLAG not in flags(row), row
    assert row["confidence"] == "high", row
    _assert_control(table, w)
