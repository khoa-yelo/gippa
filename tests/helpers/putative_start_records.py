"""Replicate records of the gippasim scenario ``alt_start_initiation`` (M12-04 PR 10; statistician
sign-off ``local/handoff/pr10_statistician.md`` C6) for the calibration test
(tests/calibration/test_putative_start_calibration.py) and its attack
(tests/calibration/attacks/test_attack_putative_start.py).

One replicate = one ``gippasim.alt_start.alt_start_initiation`` world through the real ``variants``
chain (``gippa.steps.variants.run_variants_qc`` and ``qc_frame``). :func:`record` evaluates the
sign-off's thresholds (``putative_start_oracle.checks``) on the chain's output, or on the output a
mutant of the rule would give:

- rule mutants 1-5 (``putative_start_oracle.RULE_MUTANTS``): the mutant's flagged set from the
  oracle on the chain's alignments; a row it flags is ``candidate``, every other row has its
  uncapped confidence, and ``id_mixups`` is recomputed on the mutant's flags;
- ``no_cap`` (6): the chain's flags with the uncapped confidence;
- ``exempt_in_id_mixups`` (7): ``id_mixups`` with ``putative_start`` in ``EXEMPT_FLAGS``.

Seeds: ``_criteria.suite_seed(STREAM, i)``; the arm and tier pick the scenario arguments.
"""

from __future__ import annotations

from functools import lru_cache
from typing import Any

import pandas as pd

from helpers import putative_start_oracle as O

STREAM = "variants_putative_start"
#: scenario arms: the two genome relations (sign-off) and the table-4 TTA arm (mutant 3)
ARMS: dict[str, dict[str, Any]] = {
    "related": {"genome_relation": "related"},
    "same": {"genome_relation": "same"},
    "table4": {"genome_relation": "related", "translation_table": 4,
               "start_codons": ("ATG", "GTG", "TTG", "TTA")},
}  # fmt: skip
SIZES = {False: "100k", True: "1M"}


def seed_of(i: int) -> int:
    import _criteria as C

    return C.suite_seed(STREAM, i)


def scenario(i: int, arm: str, full: bool) -> Any:
    from gippasim.alt_start import alt_start_initiation

    return alt_start_initiation(seed_of(i), size=SIZES[full], **ARMS[arm])


@lru_cache(maxsize=4)
def replicate(i: int, arm: str, full: bool) -> dict[str, Any]:
    """The chain's output on replicate ``i``: scenario, table, support, alignments, id_mixups."""
    from gippa.steps.variants import qc_frame, run_variants_qc

    scn = scenario(i, arm, full)
    _, table, records, support = run_variants_qc(
        scn.peptide_hits, None, scn.peptides, scn.gm, scn.genome, scn.cfg
    )
    qc = qc_frame(table, support).set_index("metric").loc["id_mixups"]
    aligns = O.chain_alignments(scn.peptide_hits, scn.genome, scn.gm, scn.cfg)
    return {"scn": scn, "table": table, "records": records, "support": support,
            "aligns": aligns, "id_mixups": (int(qc["value"]), int(qc["n"]))}  # fmt: skip


def _flagged(table: pd.DataFrame) -> set[str]:
    return {str(v) for v, f in zip(table["variant_id"], table["flags"], strict=True)
            if O.FLAG in O.row_flags(f)}  # fmt: skip


def _with_flags(table: pd.DataFrame, flagged: set[str]) -> pd.DataFrame:
    out = table.copy()
    out["flags"] = [
        tuple(x for x in O.row_flags(f) if x != O.FLAG) + ((O.FLAG,) if str(v) in flagged else ())
        for v, f in zip(out["variant_id"], out["flags"], strict=True)
    ]
    return out


def outcome(rep: dict[str, Any], mutant: str | None = None) -> dict[str, Any]:
    """``(flagged, confidence, id_mixups)`` of the chain or of ``mutant``."""
    scn, table, support = rep["scn"], rep["table"], rep["support"]
    uncapped = dict(zip(table["variant_id"].astype(str), O.uncapped_confidence(table, scn.cfg),
                        strict=True))  # fmt: skip
    own = dict(zip(table["variant_id"].astype(str), table["confidence"].astype(str), strict=True))
    if mutant is None:
        return {"flagged": _flagged(table), "confidence": own, "id_mixups": rep["id_mixups"]}
    if mutant == "no_cap":
        return {"flagged": _flagged(table), "confidence": uncapped, "id_mixups": rep["id_mixups"]}
    if mutant == "exempt_in_id_mixups":
        return {"flagged": _flagged(table), "confidence": own,
                "id_mixups": O.id_mixups_exempting(table, support, {O.FLAG})}  # fmt: skip
    from gippa.variants.qc import id_mixups

    flagged = O.oracle_flagged(table, support, rep["aligns"], scn.genome, scn.gm, scn.cfg, mutant)
    conf = {v: ("candidate" if v in flagged else c) for v, c in uncapped.items()}
    return {"flagged": flagged, "confidence": conf,
            "id_mixups": id_mixups(_with_flags(table, flagged), support)}  # fmt: skip


def record(i: int, arm: str, full: bool, mutant: str | None = None) -> dict[str, Any]:
    """The sign-off's per-replicate counts and verdicts (``putative_start_oracle.checks``), plus
    ``oracle_agrees`` (the oracle's flagged set equals the chain's) and the unguarded B2 count."""
    rep = replicate(i, arm, full)
    scn, table = rep["scn"], rep["table"]
    out = outcome(rep, mutant)
    unguarded = O.oracle_flagged(table, rep["support"], rep["aligns"], scn.genome, scn.gm,
                                 scn.cfg, "no_kr_guard")  # fmt: skip
    rec = O.checks(table, scn.truth, out["flagged"], out["confidence"], out["id_mixups"],
                   unguarded)  # fmt: skip
    oracle = O.oracle_flagged(table, rep["support"], rep["aligns"], scn.genome, scn.gm, scn.cfg)
    rec |= {"replicate": i, "arm": arm, "full": full, "mutant": mutant,
            "oracle_agrees": oracle == _flagged(table), "n_rows": len(table)}  # fmt: skip
    return rec


def failed(rec: dict[str, Any]) -> list[str]:
    return [k for k, ok in rec["ok"].items() if not ok]
