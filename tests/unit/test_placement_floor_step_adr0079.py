"""ADR-0079 Option 1 through the registered ``placement`` step: one WARNING
``placement.stratum_q_floor_above_alpha`` per list whose q floor c/T is strictly above
``stats.alpha.placement``, with data ``{decoy_type, n_diff, n_targets, n_decoys, q_floor, alpha}``,
ordered by decoy type then m; and the manifest route (``placement:<event>`` in ``warnings``, not
a run flag).

Inputs: one 50-kb gippasim database chain (``helpers.placement_sims.run_chain``, as
tests/unit/test_placement_step_adr0019.py). Two checks on it:

* the natural run: the logged warnings are exactly the breaches implied by the step's own outputs
  (``mapping_fdr`` ``all`` rows; c = 1 for m >= 1, the largest ``stratum0_classes`` pi-hat for
  database stratum 0);
* a rerun of the step with ``best_costs`` replaced by the ADR row-9 costs table
  (``helpers.placement_floor.row9_costs``): exactly two warnings, reverse m = 1 (T = 13, D = 7)
  then shuffle m = 2 (T = 68, D = 0).

The message is checked for its fixed parts only (number formats are open).
"""

from __future__ import annotations

import logging
import math
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from helpers.manifest_data import build, event, validate
from helpers.mapping_steps import fresh_steps, key_like, mapping_module
from helpers.placement_floor import DATA_KEYS, row9_costs
from helpers.placement_sims import Chain, run_chain

SEED = 20261124
EVENT = "placement.stratum_q_floor_above_alpha"


class _Sink(logging.Handler):
    def __init__(self) -> None:
        super().__init__(logging.INFO)
        self.records: list[logging.LogRecord] = []

    def emit(self, record: logging.LogRecord) -> None:
        self.records.append(record)


class _capture:
    def __enter__(self) -> list[logging.LogRecord]:
        self.sink = _Sink()
        self.log = logging.getLogger("gippa.steps")
        self.old = self.log.level
        self.log.addHandler(self.sink)
        self.log.setLevel(logging.INFO)
        return self.sink.records

    def __exit__(self, *exc: object) -> None:
        self.log.removeHandler(self.sink)
        self.log.setLevel(self.old)


@pytest.fixture(scope="module")
def chain(tmp_path_factory: pytest.TempPathFactory) -> Iterator[tuple[Chain, list]]:
    from gippasim.proteomics import ProteomicsConfig

    with _capture() as records:
        ch = run_chain(
            "database", SEED, tmp_path_factory.mktemp("floor"), length=50_000,
            proteomics=ProteomicsConfig(pi_wrong=0.05, search_db="other", n_runs=1),
        )  # fmt: skip
    yield ch, list(records)


def floor_events(records: list[logging.LogRecord]) -> list[logging.LogRecord]:
    return [r for r in records if getattr(r, "event", None) == EVENT]


def key(data: dict[str, Any]) -> tuple[str, int, int, int]:
    return (str(data["decoy_type"]), int(data["n_diff"]), int(data["n_targets"]),
            int(data["n_decoys"]))  # fmt: skip


def check_record(r: logging.LogRecord, alpha: float) -> dict[str, Any]:
    assert r.levelno == logging.WARNING, r.levelname
    assert r.name == "gippa.steps.placement", r.name
    data = dict(r.data)  # type: ignore[attr-defined]
    assert set(data) == DATA_KEYS, data
    assert float(data["alpha"]) == alpha
    msg = r.getMessage()
    dt, m, t = data["decoy_type"], int(data["n_diff"]), int(data["n_targets"])
    assert msg.startswith("placement: "), msg
    assert f"stratum {dt} m={m} has {t} targets" in msg, msg
    assert "q floor" in msg and "alpha" in msg and "none can be accepted" in msg, msg
    return data


def expected_from_outputs(ch: Chain) -> list[tuple[tuple[str, int, int, int], float]]:
    alpha = float(ch.cfg.stats.alpha.placement)
    fdr = pd.read_parquet(ch.avail["mapping_fdr"])
    fdr = fdr[(fdr["section"] == "placement") & (fdr["length_class"] == "all")]
    classes = pd.read_parquet(ch.avail[key_like(list(ch.avail), "stratum0_classes")])
    present = classes.loc[classes["n"] > 0, "pi_wrong"]
    c_db0 = float(present.max()) if len(present) else 1.0
    out = []
    for r in fdr.sort_values(["decoy_type", "n_diff"], kind="stable").to_dict("records"):
        t = int(r["targets"])
        if t == 0:
            continue
        if int(r["n_diff"]) >= 1:
            c = 1.0
        elif r["decoy_type"] == "reverse":
            c = c_db0
        else:
            pytest.skip("a database chain has no de novo stratum-0 targets")
        if c / t > alpha:
            out.append(((str(r["decoy_type"]), int(r["n_diff"]), t, int(r["decoys"])), c / t))
    return out


def test_natural_run_warns_exactly_for_the_lists_above_alpha(chain: tuple[Chain, list]) -> None:
    ch, records = chain
    alpha = float(ch.cfg.stats.alpha.placement)
    want = expected_from_outputs(ch)
    got = [check_record(r, alpha) for r in floor_events(records)]
    assert [key(d) for d in got] == [w[0] for w in want], (got, want)
    for d, (_, floor) in zip(got, want, strict=True):
        assert float(d["q_floor"]) == pytest.approx(floor, rel=1e-9)
    fdr = pd.read_parquet(ch.avail["mapping_fdr"])
    for d in got:
        sel = (fdr["section"] == "placement") & (fdr["length_class"] == "all")
        row = fdr[sel & (fdr["decoy_type"] == d["decoy_type"]) & (fdr["n_diff"] == d["n_diff"])]
        assert int(row["n_accepted"].iloc[0]) == 0


def test_row9_costs_through_the_step_give_two_ordered_warnings(
    chain: tuple[Chain, list], tmp_path: Path
) -> None:
    ch, _ = chain
    alpha = float(ch.cfg.stats.alpha.placement)
    assert alpha == 0.01
    costs = row9_costs()
    with pytest.MonkeyPatch.context() as mp:
        reg = fresh_steps(mp)
        mp.setattr(mapping_module(), "best_costs", lambda *a, **k: costs.copy())
        spec = reg.get("placement")
        ins = {k: ch.avail[k] for k in spec.inputs}
        ins |= {k: ch.avail[k] for k in spec.optional_inputs if k in ch.avail}
        outs = {k: tmp_path / k for k in spec.outputs}
        with _capture() as records:
            reg.run_step("placement", ch.cfg, ins, outs)
        got = [check_record(r, alpha) for r in floor_events(records)]
    assert [key(d) for d in got] == [("reverse", 1, 13, 7), ("shuffle", 2, 68, 0)], got
    assert math.isclose(float(got[0]["q_floor"]), 1 / 13, rel_tol=1e-12)
    assert math.isclose(float(got[1]["q_floor"]), 1 / 68, rel_tol=1e-12)
    fdr = pd.read_parquet(outs["mapping_fdr"])
    allr = fdr[(fdr["section"] == "placement") & (fdr["length_class"] == "all")]
    acc = {(r["decoy_type"], int(r["n_diff"])): int(r["n_accepted"]) for _, r in allr.iterrows()}
    assert acc[("reverse", 1)] == 0 and acc[("shuffle", 2)] == 0
    assert acc[("reverse", 0)] == 400 and acc[("shuffle", 1)] == 200


# ---------------------------------------------------------------------------- manifest route


def _floor_event(**over: Any) -> dict[str, Any]:
    data = {"decoy_type": "reverse", "n_diff": 1, "n_targets": 13, "n_decoys": 7,
            "q_floor": 1 / 13, "alpha": 0.01}  # fmt: skip
    data.update(over)
    return event("placement", EVENT, data, level="WARNING")


def test_the_warning_reaches_the_manifest_but_not_the_flags() -> None:
    m = build(events=[_floor_event(), _floor_event(decoy_type="shuffle", n_diff=2, n_targets=68,
                                                   n_decoys=0, q_floor=1 / 68)])  # fmt: skip
    assert f"placement:{EVENT}" in m["warnings"]
    assert m["warnings"].count(f"placement:{EVENT}") == 1, m["warnings"]
    assert not any("q_floor" in f for f in m["flags"]), m["flags"]
    validate(m)
