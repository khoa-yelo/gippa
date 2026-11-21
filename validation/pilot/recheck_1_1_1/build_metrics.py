"""Build validation/pilot/metrics/metrics_1_1_1.json from the recheck outputs in this directory.

The 1.1.1 addendum (ADR-0078, PR #607): pilot metrics whose value changes when only
gippa/genetable/detection.py is swapped from 1.1.0 to the fix. metrics.json stays the 1.0.0
freeze record; this file lists the superseding values under the same key paths. Aggregates only.
Usage: build_metrics.py [OUT]   (default ../metrics/metrics_1_1_1.json)
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PILOT = HERE.parent
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else PILOT / "metrics" / "metrics_1_1_1.json"
ARM_KEYS = ("genes_q01", "maxq", "gene_paired", "unpaired_own", "own_units", "decoy_wins")


def _jsonl(name: str) -> dict[str, dict]:
    rows = (json.loads(ln) for ln in (HERE / name).read_text().splitlines() if ln.strip())
    return {r["run"]: r for r in rows}


def _share(row: list) -> dict:
    k, n, v, lo, hi = row
    d = {"k": k, "n": n, "value": v, "ci_low": lo, "ci_high": hi}
    if n == 0:
        d["reason"] = "empty_denominator"
    return d


def _changes(old: object, new: object, path: str) -> list[dict]:
    if isinstance(old, dict) and isinstance(new, dict):
        out: list[dict] = []
        for k in sorted(set(old) | set(new)):
            out += _changes(old.get(k), new.get(k), f"{path}.{k}" if path else k)
        return out
    return [] if old == new else [{"key": path, "old": old, "new": new}]


def main() -> None:
    before, after = _jsonl("out_v110.jsonl"), _jsonl("out_fix.jsonl")
    ctl = json.loads((HERE / "derived_ctl.json").read_text())
    fix = json.loads((HERE / "derived_fix.json").read_text())
    rec = json.loads((PILOT / "metrics" / "metrics.json").read_text())
    arms = [
        {"run": run, "commit": b["commit"], "tier": b["tier"], "eligible": b["eligible"]}
        | {f"{k}_before": b[k] for k in ARM_KEYS}
        | {f"{k}_after": after[run][k] for k in ARM_KEYS}
        for run, b in before.items()
    ]
    m3 = {a: {"subsets": {s: _share(r) for s, r in fix["m3"][a].items()}} for a in ("A1", "A2")}
    new = {"m3": m3, "m4": fix["m4"], "m5": fix["m5"], "h1": fix["h1"], "h3": fix["h3"]}
    old = {
        "m3": {a: {"subsets": {s: _share(r) for s, r in ctl["m3"][a].items()}} for a in m3},
        "m4": ctl["m4"],
        "m5": ctl["m5"],
        "h1": ctl["h1"],
        "h3": ctl["h3"],
    }
    # the control reproduces the recorded 1.0.0 values (shares: k, n and value)
    for a in m3:
        for s, d in old["m3"][a]["subsets"].items():
            r = rec["m3"][a]["subsets"][s]["share"]
            assert (d["k"], d["n"], d["value"]) == (r["k"], r["n"], r["value"]), (a, s)
    assert old["m4"] == rec["m4"] and old["h1"] == rec["h1"], "control differs from metrics.json"
    assert old["m5"]["enriched_categories"] == rec["m5"]["enriched_categories"]
    doc = {
        "schema": "gippa-pilot-addendum/1",
        "addendum": "1.1.1 detection recheck (ADR-0078)",
        "supersedes": "validation/pilot/metrics/metrics.json",
        "adr": "ADR-0078",
        "pr": 607,
        "code_before": "a2a95c76 (1.1.0)",
        "code_after": "7d4c4bfd (ADR-0078 fix, squash-merged as 92fbf3af)",
        "freeze_commit": rec["gippa_commit"],
        "method": (
            "freeze code fa49b8f with gippa/genetable/detection.py swapped; gene_table, "
            "detectability, context_stats, function and compare rerun on the freeze outputs"
        ),
        "slurm": {"control": 2001992, "fix": 2002006},
        "control_matches_recorded": True,
        "detection_recheck": arms,
        **new,
        "changed": _changes(old, new, ""),
        "verdicts_unchanged": ["h1", "h3", "te3", "te4", "m3 plasmid", "m3 prophage"],
    }
    OUT.write_text(json.dumps(doc, indent=1, sort_keys=False) + "\n", encoding="utf-8")
    print(f"{OUT}: {len(doc['changed'])} changed keys")


if __name__ == "__main__":
    main()
