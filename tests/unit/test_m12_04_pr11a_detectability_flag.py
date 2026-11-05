"""M12-04 PR 11a, T-06 (second half): ``detectability_constant`` moves from the Gene rows to the run
(docs/specs/M12-04-outputs-revision.md § PR 11a, § 5.4).

AC 11a.2: no row of any ``genes.tsv`` contains ``detectability_constant``; on a
constant-detectability fixture the manifest ``flags`` contains it once. The ``detection_summary``
row ``detectability / overall / all / detectability_constant`` stays.

Interface assumed (test-author reading): the token reaches the manifest the way every other run
flag does (tests/unit/steps/test_run_flag_events.py): the ``detectability`` step logs one record
with ``extra={"event": <any name>, "data": {"flag": "detectability_constant"}}`` when the model is
constant, and ``gippa.manifest.build_manifest`` unions ``data.flag``. The token is registered in
``gippa.manifest.RUN_FLAGS`` and leaves ``GENE_FLAGS`` and the genes ``flags`` pattern (§9.4 moves
it from Gene to Run; ``docs/reference/flags.md`` follows the schema).

Fixture: ``gene_table_fixture`` (20 genes, fewer than 50 per class, so M4-13 fits the constant
model), gene 5 partial.
"""

from __future__ import annotations

import logging
import re
import sys
from collections.abc import Iterator, Mapping
from contextlib import contextmanager
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from helpers.manifest_data import build, event
from hypothesis import given, settings
from hypothesis import strategies as st

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "tests" / "helpers"))
sys.path.insert(0, str(REPO / "tests" / "unit" / "genetable"))
import gene_table_fixture as gt  # noqa: E402
import m4_contract_fixture as mf  # noqa: E402
from gene_table_fixture import as_list, tag  # noqa: E402

TOKEN = "detectability_constant"
FLAGS_MD = REPO / "docs" / "reference" / "flags.md"
GENES_SCHEMA = REPO / "src" / "gippa" / "schemas" / "tables" / "genes.schema.json"


class _Recorder(logging.Handler):
    def __init__(self) -> None:
        super().__init__(level=logging.DEBUG)
        self.records: list[logging.LogRecord] = []

    def emit(self, record: logging.LogRecord) -> None:
        self.records.append(record)


@contextmanager
def _recording() -> Iterator[list[logging.LogRecord]]:
    rec = _Recorder()
    log = logging.getLogger("gippa.steps")
    old = log.level
    log.addHandler(rec)
    log.setLevel(logging.DEBUG)
    try:
        yield rec.records
    finally:
        log.removeHandler(rec)
        log.setLevel(old)


def _flag_records(records: list[logging.LogRecord]) -> list[logging.LogRecord]:
    out = []
    for r in records:
        ev, data = getattr(r, "event", None), getattr(r, "data", None)
        if isinstance(ev, str) and isinstance(data, Mapping) and data.get("flag") == TOKEN:
            out.append(r)
    return out


def _as_event(r: logging.LogRecord) -> dict[str, Any]:
    return event(r.name.rsplit(".", 1)[-1], str(r.event), dict(r.data), level=r.levelname)  # type: ignore[attr-defined]


def _frames() -> dict[str, Any]:
    fr = gt.frames_with()
    fr["genes"].loc[5, "partial"] = True
    fr["layers"].at[5, "flags"] = ["partial"]
    return fr


def _write_work(root: Path, frames: dict[str, Any]) -> None:
    from gippa.io.work_tables import default_paths

    mf.write_work(frames, default_paths(root) | {"preflight": root / "preflight.json"})


@pytest.fixture
def constant_run(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> tuple[pd.DataFrame, list]:
    """gene_table + detectability on the 20-gene fixture: (genes.parquet frame, step records)."""
    from gippa.config import GippaConfig
    from gippa.steps.detectability import run_detectability
    from gippa.steps.gene_table import run_gene_table

    (tmp_path / "db").mkdir()
    monkeypatch.setenv("GIPPA_DB_DIR", str(tmp_path / "db"))
    work = tmp_path / "work"
    work.mkdir()
    _write_work(work, _frames())
    cfg = GippaConfig(seed=3)
    run_gene_table(work, cfg)
    with _recording() as records:
        out = run_detectability(work, cfg)
    parts = work / "genetable" / "detection_summary_parts" / "detectability.parquet"
    (tmp_path / "summary_part").write_bytes(parts.read_bytes())
    return pd.read_parquet(out), list(records)


# ----------------------------------------------------------------------------- vocabulary


def test_token_is_no_gene_flag() -> None:
    from gippa.genetable.assemble import GENE_FLAGS

    assert TOKEN not in GENE_FLAGS


def test_token_is_a_registered_run_flag() -> None:
    from gippa.manifest import RUN_FLAGS

    assert TOKEN in RUN_FLAGS


def test_genes_flags_pattern_rejects_the_token() -> None:
    import json

    fields = {f["name"]: f for f in json.loads(GENES_SCHEMA.read_text("utf-8"))["fields"]}
    pattern = fields["flags"]["constraints"]["pattern"]
    assert re.fullmatch(pattern, "partial")
    assert not re.fullmatch(pattern, TOKEN)
    assert not re.fullmatch(pattern, f"partial;{TOKEN}")


def test_flags_md_genes_section_lacks_the_token() -> None:
    text = FLAGS_MD.read_text(encoding="utf-8")
    m = re.search(r"^## `genes`\n(.*?)(?=^## |\Z)", text, re.S | re.M)
    assert m, "no genes section in flags.md"
    assert f"`{TOKEN}`" not in m.group(1)


# ----------------------------------------------------------------------------- rows


def test_assemble_with_constant_model_sets_no_row_flag() -> None:
    t = gt.build_world(frames=_frames(), constant=True).assemble()
    assert not any(TOKEN in as_list(f) for f in t["flags"])
    # the other flags stay (M4-15 example: gene 7 is a pseudogene with top3_lt3)
    assert as_list(t.set_index("locus_tag").loc[tag(7), "flags"]) == ["pseudogene", "top3_lt3"]


def test_detectability_step_writes_no_row_flag(constant_run: tuple[pd.DataFrame, list]) -> None:
    df, _ = constant_run
    assert len(df) == 20
    assert not any(TOKEN in as_list(f) for f in df["flags"])
    assert as_list(df.set_index("locus_tag").loc[tag(5), "flags"]) == ["partial"]


def test_constant_model_still_fills_p_detect(constant_run: tuple[pd.DataFrame, list]) -> None:
    """Only the flag moves: p_detect is still the overall detected share (M4-13)."""
    df, _ = constant_run
    by = df.set_index("locus_tag")
    eligible = ~by["partial"].astype(bool)
    rate = float(by.loc[eligible, "detected"].astype(bool).mean())
    assert (by.loc[eligible, "p_detect"].astype(float) - rate).abs().max() < 1e-12


def test_genes_tsv_written_from_the_constant_table_has_no_token(
    tmp_path: Path, constant_run: tuple[pd.DataFrame, list]
) -> None:
    from gippa.io.genes_tsv import read_genes_parquet, write_genes_tsv

    df, _ = constant_run
    pq = tmp_path / "g.parquet"
    df.to_parquet(pq)
    path = tmp_path / "genes.tsv"
    write_genes_tsv(path, read_genes_parquet(pq))
    assert TOKEN not in path.read_text(encoding="utf-8")


# ----------------------------------------------------------------------------- run flag


def test_detectability_step_emits_the_run_flag_once(
    constant_run: tuple[pd.DataFrame, list],
) -> None:
    _, records = constant_run
    hits = _flag_records(records)
    assert len(hits) == 1, [(r.name, r.getMessage()) for r in records]
    assert hits[0].name == "gippa.steps.detectability"


def test_manifest_flags_contain_the_token_once(constant_run: tuple[pd.DataFrame, list]) -> None:
    _, records = constant_run
    events = [_as_event(r) for r in records if isinstance(getattr(r, "event", None), str)]
    m = build(events=events)
    assert m["flags"].count(TOKEN) == 1, m["flags"]


@settings(max_examples=30, deadline=None)
@given(
    n=st.integers(min_value=1, max_value=6),
    others=st.lists(
        st.sampled_from(["c_assumed", "bakta_db_light", "anti_sd_default"]), max_size=4
    ),
)
def test_property_manifest_lists_the_run_flag_once(n: int, others: list[str]) -> None:
    """However many times a resumed run logs it, the manifest lists the token once."""
    evs = [event("detectability", "detectability.run_flag", {"flag": TOKEN})] * n
    evs += [event("placement", "placement.run_flag", {"flag": o}) for o in others]
    m = build(events=evs)
    assert m["flags"].count(TOKEN) == 1
    assert m["flags"] == sorted(set(m["flags"]))


def test_detection_summary_overall_row_stays(
    tmp_path: Path, constant_run: tuple[pd.DataFrame, list]
) -> None:
    diag = pd.read_parquet(tmp_path / "summary_part")
    row = diag[(diag["section"] == "detectability") & (diag["metric"] == TOKEN)]
    assert len(row) == 1 and float(row["value"].iloc[0]) == 1.0, diag


# ----------------------------------------------------------------------------- fitted model


@pytest.fixture
def fitted_run(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> tuple[pd.DataFrame, list]:
    """gene_table on a 200-gene fixture, then ``genes_pre`` gets 100 detected / 99 undetected
    eligible genes with varied covariates (gene 5 partial, ~100 TUs), so M4-13 fits the model."""
    import numpy as np

    from gippa.config import GippaConfig
    from gippa.genetable.covariates import COVARIATE_COLUMNS
    from gippa.io.genes_tsv import read_genes_parquet, write_genes_parquet
    from gippa.steps.detectability import run_detectability
    from gippa.steps.gene_table import run_gene_table

    (tmp_path / "db").mkdir()
    monkeypatch.setenv("GIPPA_DB_DIR", str(tmp_path / "db"))
    work = tmp_path / "work"
    work.mkdir()
    fr = gt.frames_with(200)
    fr["genes"].loc[5, "partial"] = True
    fr["layers"].at[5, "flags"] = ["partial"]
    _write_work(work, fr)
    cfg = GippaConfig(seed=3)
    pre_path = run_gene_table(work, cfg)
    pre = read_genes_parquet(pre_path)
    n = len(pre)
    r = np.random.default_rng(11)
    for c in COVARIATE_COLUMNS:
        pre[c] = r.uniform(1.0, 500.0, n)
    pre["detected"] = np.arange(n) % 2 == 0
    pre["tu_id"] = [f"TU_{i // 2 + 1}" for i in range(n)]
    write_genes_parquet(pre_path, pre)
    with _recording() as records:
        out = run_detectability(work, cfg)
    return pd.read_parquet(out), list(records)


def test_fitted_model_emits_no_run_flag(fitted_run: tuple[pd.DataFrame, list]) -> None:
    """Kills ``if res.constant:`` → ``if True:``: only a constant model logs the run flag."""
    df, records = fitted_run
    by = df.set_index("locus_tag")
    eligible = ~by["partial"].astype(bool)
    p = by.loc[eligible, "p_detect"].astype(float)
    assert p.notna().all() and p.nunique() > 1, "fixture did not fit a non-constant model"
    assert _flag_records(records) == [], [(r.name, r.getMessage()) for r in records]
    events = [_as_event(r) for r in records if isinstance(getattr(r, "event", None), str)]
    m = build(events=events)
    assert TOKEN not in m["flags"], m["flags"]
    assert not any(TOKEN in str(w) for w in m.get("warnings", [])), m.get("warnings")
