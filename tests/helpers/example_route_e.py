"""AC 22 of M12-01 route E (spec docs/specs/M12-01-example-dataset.md, amendment 2026-11-06): the
route-E-only assertions of the slow route test, as a pure function of the measured keys
(``tests/helpers/example_routes.measure``) so it is unit-tested without a run."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

TIER_E = "mayu"  # Q-L10


def route_e_misses(measured: Mapping[str, Any]) -> list[str]:
    """Empty when route E shows the six-frame workflow as AC 22 requires: tier ``mayu``,
    ``n_novel_called = 0``, and, when there are candidates, fewer unpaired candidates than
    candidates and at least one ``picked_decoy`` class (the rule-9 rewrite makes them pair)."""
    misses = []
    if measured.get("detection_tier") != TIER_E:
        misses.append(f"detection_tier: got {measured.get('detection_tier')!r}, want {TIER_E!r}")
    if measured.get("n_novel_called") != 0:
        misses.append(f"n_novel_called: got {measured.get('n_novel_called')!r}, want 0")
    n = measured.get("n_novel_candidates")
    if not isinstance(n, int):
        return [*misses, f"n_novel_candidates: not measured ({n!r})"]
    if n > 0:
        unpaired = measured.get("n_novel_unpaired")
        if not (isinstance(unpaired, int) and unpaired < n):
            misses.append(f"n_novel_unpaired: got {unpaired!r}, want < {n} (rule 9 pairs them)")
        classes = measured.get("novel_classes") or []
        if not any(isinstance(c, Mapping) and c.get("regime") == "picked_decoy" for c in classes):
            misses.append("novel_classes: no picked_decoy entry")
    return misses
