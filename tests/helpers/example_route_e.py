"""AC 22 of M12-01 route E (spec docs/specs/M12-01-example-dataset.md, amendments 2026-11-06): the
route-E-only assertions of the slow route test, as a pure function of the measured keys
(``tests/helpers/example_routes.measure``) so it is unit-tested without a run."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

TIER_E = "decoy_picked"  # Q-L15, measured in trial job 1943965


def _is_int(v: Any) -> bool:
    return isinstance(v, int) and not isinstance(v, bool)


def route_e_misses(measured: Mapping[str, Any]) -> list[str]:
    """Empty when route E shows the six-frame workflow as AC 22 requires (amended 2026-11-06 (b),
    route E's own region): tier ``decoy_picked``; ``n_novel_called = 0``; at least one candidate; an
    integer ``n_novel_unpaired`` (0 when the FDR table is empty or absent, never null) below the
    number of candidates; and at least one ``picked_decoy`` class (the rule-9 rewrite makes them
    pair)."""
    misses = []
    if measured.get("detection_tier") != TIER_E:
        misses.append(f"detection_tier: got {measured.get('detection_tier')!r}, want {TIER_E!r}")
    if measured.get("n_novel_called") != 0:
        misses.append(f"n_novel_called: got {measured.get('n_novel_called')!r}, want 0")
    n = measured.get("n_novel_candidates")
    if not _is_int(n):
        return [*misses, f"n_novel_candidates: not measured ({n!r})"]
    unpaired = measured.get("n_novel_unpaired")
    if not _is_int(unpaired):
        misses.append(f"n_novel_unpaired: got {unpaired!r}, want an integer (0, never null)")
    if n < 1:
        return [*misses, f"n_novel_candidates: got {n}, want >= 1 (route E's region, AC 22)"]
    if _is_int(unpaired) and not unpaired < n:
        misses.append(f"n_novel_unpaired: got {unpaired!r}, want < {n} (rule 9 pairs them)")
    classes = measured.get("novel_classes") or []
    if not any(isinstance(c, Mapping) and c.get("regime") == "picked_decoy" for c in classes):
        misses.append("novel_classes: no picked_decoy entry")
    return misses
