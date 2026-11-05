"""M12-04 PR 12 review round 1, M1: ``scripts/example/summarize_runs.py`` reads the regions table
in both layouts: ``OUT/intermediate/regions.parquet`` (1.1, docs/specs/M12-04-outputs-revision.md
> PR 12 "Publish rule") and ``OUT/regions.parquet`` (1.0, § 5.5)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import ModuleType

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[2]
REGIONS = pd.DataFrame({"region_id": ["r1", "r2", "r3"],
                        "region_type": ["intergenic", "intergenic", "antisense"]})  # fmt: skip
WANT = {"intergenic": 2, "antisense": 1}


def _module() -> ModuleType:
    path = ROOT / "scripts" / "example" / "summarize_runs.py"
    spec = importlib.util.spec_from_file_location("gippa_example_summarize_runs", path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.mark.parametrize("rel", ["intermediate/regions.parquet", "regions.parquet"],
                         ids=["1.1", "1.0"])  # fmt: skip
def test_regions_counted_in_either_layout(tmp_path: Path, rel: str) -> None:
    out = tmp_path / "out"
    (out / rel).parent.mkdir(parents=True, exist_ok=True)
    REGIONS.to_parquet(out / rel)
    res = _module().results(out)
    assert res.get("regions") == WANT, f"{rel}: {res}"


def test_no_regions_table_no_regions_key(tmp_path: Path) -> None:
    out = tmp_path / "out"
    out.mkdir()
    assert "regions" not in _module().results(out)
