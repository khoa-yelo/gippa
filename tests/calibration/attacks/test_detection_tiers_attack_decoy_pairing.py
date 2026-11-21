"""Attack null for ADR-0078: picked gene FDR when decoys pair only through accepted-peptide votes
(M4-12; GIPPA_SPEC.md §7.3.2, §8.4, §8.5;
docs/adr/ADR-0078-decoy-pairing-independent-of-acceptance.md "Calibration attack").

Generator: the ``gene_fdr`` replicate of tests/helpers/detection_tiers_calsim.py (M4-07 -> M4-08 ->
M4-09, pi_wrong = 0.05, decoys emitted) at 4,500 genes in both modes (2,000 genes lack power in
the fast mode). Every target protein ID ``GSIM_<k>`` and its decoy ``rev_GSIM_<k>`` is re-labelled
to ``ACC<perm(k):06d>`` / ``rev_ACC<...>`` with a permutation drawn from the replicate rng; C11
placements stay, so targets still place by sequence. protein_hits is empty, so the protein -> gene
map is vote-only (as in campaign run N2): a decoy whose protein has no accepted target cannot pair.
Under 1.1.0 those decoys (exactly the ones that would win their gene) leave D and 1/f does not
repair the loss (about 0.0128 +- 0.0005 FDP at q = 0.01 against a bound of about 0.011).

Claims (criterion ``fdr``, ADR-0015 two-block confirmation, suite ``gene_fdr``; stream
``gene_fdr.attack_decoy_pairing.reid100``; sim stream ``gene_fdr.attack_decoy_pairing.sim.reid100``
(+ ``.full``)), each at q in {0.01, 0.05, 0.10} with detection rerun at
``stats.alpha.identification`` = q: mean FDP <= q + k SE (k = ``se_multiplier``).

- The tier must be ``decoy_picked`` in every replicate (gene-paired share >= 0.86 in 50/50
  simulated). A replicate in another tier fails the test through ``pytest.fail``, which the
  ``xfail(raises=AssertionError)`` marks below never absorb.
- q = 0.01 and q = 0.05 are confirmed fails on 1.1.0: ``xfail(strict=True)`` until the ADR-0078
  fix lands (then the marks go). q = 0.10 is saturated by the peptide q cap and passes either way.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import _criteria as C
import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import detection_tiers_calsim as cs

MODES = [
    pytest.param(
        False, marks=[pytest.mark.calibration, pytest.mark.calibration_nightly], id="fast"
    ),
    pytest.param(True, marks=pytest.mark.calibration_full, id="full"),
]
SUITE = "gene_fdr"
PI = 0.05
N_GENES = 4500  # both modes (ADR-0078)
STREAM = "gene_fdr.attack_decoy_pairing.reid100"
SIM_STREAM = "gene_fdr.attack_decoy_pairing.sim.reid100"
XFAIL = pytest.mark.xfail(
    strict=True,
    raises=AssertionError,
    reason="ADR-0078: 1.1.0 drops decoys whose protein has no accepted target (vote-only "
    "pairing), so picked gene FDR is anti-conservative; fixed in 1.1.1",
)
LEVELS = [
    pytest.param(0.01, marks=XFAIL, id="q0.01"),
    pytest.param(0.05, marks=XFAIL, id="q0.05"),
    pytest.param(0.10, id="q0.10"),
]
_MEMO: dict[tuple[bool, int], tuple[dict[float, float], set[str]]] = {}


def relabel(peptides: pd.DataFrame, rng: np.random.Generator) -> pd.DataFrame:
    """Every ``GSIM_<k>`` / ``rev_GSIM_<k>`` protein ID -> ``ACC<perm>`` / ``rev_ACC<perm>``, one
    permutation over the distinct stripped IDs (sorted) drawn from ``rng``."""
    p = peptides.copy()
    ids = p["protein_ids"].tolist()
    base = sorted({str(x).removeprefix("rev_") for c in ids for x in c})
    perm = rng.permutation(len(base))
    new = {b: f"ACC{int(k):06d}" for b, k in zip(base, perm, strict=True)}

    def one(x: str) -> str:
        return f"rev_{new[x[4:]]}" if x.startswith("rev_") else new[x]

    p["protein_ids"] = pd.Series([[one(str(x)) for x in c] for c in ids], index=p.index,
                                 dtype=object)  # fmt: skip
    return p


def _replicate(full: bool, rep: int) -> tuple[dict[float, float], set[str]]:
    key = (full, rep)
    if key not in _MEMO:
        rng = C.replicate_rng(SIM_STREAM + (".full" if full else ""), rep)
        lay, e = cs.simulate(rng, PI, N_GENES)
        assert not len(cs.protein_hits(e)), "the attack needs an empty protein_hits (vote-only)"
        p, h = cs.masked(e, "picked")
        p = relabel(p, rng)
        _MEMO[key] = cs.run_variant(lay, e, p, h)
    return _MEMO[key]


def test_relabel_leaves_no_locus_tag_id() -> None:
    """Unmarked guard: after relabelling no protein ID is a C2 locus tag, decoys keep ``rev_``
    and each target ID keeps one image (a bijection)."""
    rng = np.random.default_rng(78)
    p = pd.DataFrame(
        {"protein_ids": [["GSIM_00005"], ["GSIM_00010", "GSIM_00005"], ["rev_GSIM_00010"]]}
    )
    out = relabel(p, rng)["protein_ids"].tolist()
    flat = [x for c in out for x in c]
    assert not any("GSIM_" in x for x in flat)
    assert out[2][0] == f"rev_{out[1][0]}"
    assert out[0][0] == out[1][1] and out[0][0] != out[1][0]


@pytest.mark.calsuite(SUITE)
@pytest.mark.parametrize("q", LEVELS)
@pytest.mark.parametrize("full", MODES)
def test_picked_fdr_with_vote_only_decoy_pairing(full: bool, q: float) -> None:
    R = C.replicates(SUITE, full)

    def fn(rng: np.random.Generator) -> Any:
        rep = cs.rep_index(STREAM, rng, R)
        fdps, tiers = _replicate(full, rep)
        if tiers != {"decoy_picked"}:
            pytest.fail(f"replicate {rep}: tier {sorted(tiers)}, ADR-0078 needs decoy_picked")
        return fdps[q]

    rec = C.calibrate(SUITE, "fdr", q, fn, R, stream=STREAM)
    assert rec["status"] != "fail", rec
