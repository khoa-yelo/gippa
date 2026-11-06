"""Fast checks of the M12-01 slow route test's own logic (tests/e2e/test_m12_01_example_routes.py):
reading ``expected/<route>.json`` tolerances, measuring a run directory and splitting the printed
command block. No tool runs; the slow test itself is Slurm-only (spec section 5, AC 17)."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path
from types import ModuleType

import pytest
from hypothesis import given
from hypothesis import strategies as st

E2E = Path(__file__).resolve().parents[1] / "e2e" / "test_m12_01_example_routes.py"


def _load() -> ModuleType:
    spec = importlib.util.spec_from_file_location("m12_01_example_routes_e2e", E2E)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod  # dataclasses look their module up
    spec.loader.exec_module(mod)
    return mod


R = _load()


def test_exact_set_and_numeric_tolerances_nested_layout() -> None:
    doc = {
        "status": {"value": "completed", "tolerance": "exact"},
        "layers_missing": {"value": ["tx_rates", "rbs_novel"], "tolerance": "set"},
        "n_genes": {"value": 158, "tolerance": "max(1, 1 %)"},
        "n_detected": {"value": 99, "tolerance": {"abs": 1, "rel": 0.02}},
        "warnings": {"value": ["x"], "tolerance": "informational"},
    }
    ok = {"status": "completed", "layers_missing": ["rbs_novel", "tx_rates"], "n_genes": 159,
          "n_detected": 98, "warnings": []}  # fmt: skip
    assert R.compare(ok, doc) == []
    bad = dict(ok, status="failed", n_genes=160, n_detected=96, layers_missing=["tx_rates"])
    misses = R.compare(bad, doc)
    assert {m.split(":")[0] for m in misses} == {"status", "n_genes", "n_detected",
                                                 "layers_missing"}  # fmt: skip


def test_flat_layout_and_missing_tolerance_or_key_fail() -> None:
    tol = {"detection_tier": "exact", "n_genes": {"abs": 1, "rel": 0.01}}
    doc = {"detection_tier": "mayu", "n_genes": 100, "convert_rows": 5, "tolerance": tol}
    misses = R.compare({"detection_tier": "mayu", "n_genes": 101}, doc)
    assert misses == ["convert_rows: not measured (expected 5)"]
    misses = R.compare({"detection_tier": "mayu", "n_genes": 101, "convert_rows": 5}, doc)
    assert misses and misses[0].startswith("convert_rows: no tolerance")


@given(exp=st.integers(0, 10_000), delta=st.integers(-500, 500))
def test_max_abs_rel_bound(exp: int, delta: int) -> None:
    bound = max(1, 0.02 * exp)
    assert R._within(exp + delta, exp, "max(1, 2 %)") == (abs(delta) <= bound)


def test_measure_reads_manifest_genes_hits_and_peptides(tmp_path: Path) -> None:
    (tmp_path / "tables").mkdir()
    (tmp_path / "mapping").mkdir()
    (tmp_path / "run_manifest.json").write_text(json.dumps({
        "status": "completed", "detection_tier": "mayu", "layers_missing": ["b", "a"],
        "warnings": ["s:w"]}))  # fmt: skip
    (tmp_path / "tables" / "genes.tsv").write_text(
        "locus_tag\tdetected\tdetected_2pep\tabundance_raw\tabundance_imputed\n"
        "g1\ttrue\ttrue\t10.5\tfalse\n"
        "g2\ttrue\tfalse\t3\ttrue\n"
        "g3\tfalse\tfalse\tNA\tfalse\n"
    )
    (tmp_path / "mapping" / "peptide_hits.tsv").write_text(
        "peptide\tgene\nAAK\tg1\nAAK\tg2\nCCR\tg1\n"
    )
    (tmp_path / "peptides.tsv").write_text("peptide\trun\nAAK\tWT_1\nCCR\tWT_1\n")
    got = R.measure("A", tmp_path)
    assert got == {"status": "completed", "detection_tier": "mayu", "layers_missing": ["a", "b"],
                   "warnings": ["s:w"], "convert_rows": 2, "n_genes": 3, "n_detected": 2,
                   "n_detected_2pep": 1, "n_abundance": 1, "n_peptides_placed": 2}  # fmt: skip
    assert "convert_rows" not in R.measure("B", tmp_path)
    assert "n_abundance" not in R.measure("C", tmp_path)


def test_printed_routes_split_by_yaml_titles() -> None:
    from gippa import example

    routes = example._routes()
    lines = example.route_lines("1.1", Path("/x/gippa-example"), routes)
    got = R.printed_routes("\n".join(lines), routes)
    assert sorted(got) == ["A", "B", "C", "D"]
    for r in routes:
        assert len(got[str(r["id"])]) == len(r["commands"])
        assert all(c.startswith("gippa ") and "{" not in c for c in got[str(r["id"])])


@pytest.mark.parametrize("extra", [[], ["--profile", "farmshare_local,apptainer"]])
def test_run_args_appended_to_gippa_run_only(extra: list[str]) -> None:
    assert R._argv("gippa convert --from diann a --out b", "/bin/g", extra)[-1] == "b"
    run = R._argv("gippa run --genome g --outdir o", "/bin/g", extra)
    assert run == ["/bin/g", "run", "--genome", "g", "--outdir", "o", *extra]
