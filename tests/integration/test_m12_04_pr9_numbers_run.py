"""M12-04 PR 9, AC 9.1: after a real ``gippa run`` on the sim, no ``number`` value of a published
TSV has more than 6 significant digits (docs/specs/M12-04-outputs-revision.md > PR 9; §6.6).

Run as tests/integration/test_plots_report_run.py (M9-13): gippasim ``SIM100K`` at
``SIM100K_SEED`` with M3-19 proteomics and every precomputed layer, peptides and abundance,
``--modules core``, ``--seed 7``, the test DB of tests/helpers/mapping_steps.py, local executor.

Pinned here (test author):

- Scanned: every ``OUT/tables/*.tsv`` and ``OUT/mapping/*.tsv`` with a shipped table schema
  (``src/gippa/schemas/tables/<stem>.schema.json``): its ``number`` fields are checked; its
  ``integer`` and ``string`` fields are exempt (spec). ``OUT/plots/data/*.tsv`` have no schema:
  a cell counts as a number there when it parses as a float and is not a plain integer
  (``-?[0-9]+``), so integer columns and labels stay exempt.
- A published table under ``tables/`` or ``mapping/`` without a schema fails the test: the AC
  speaks of "its schema", so every published TSV must have one.
- ``NA``, ``INF``/``-INF`` (and ``inf``/``-inf``) are not numbers with digits.
- Significant digits as in tests/unit/test_m12_04_pr9_numbers.py (sign, exponent, point, leading
  and trailing zeros removed).
"""

from __future__ import annotations

import dataclasses
import json
import re
from pathlib import Path
from typing import Any

import pytest
from click.testing import CliRunner

pytestmark = [pytest.mark.slow, pytest.mark.integration]

REPO = Path(__file__).resolve().parents[2]
TABLE_DIR = REPO / "src" / "gippa" / "schemas" / "tables"
SEED = 7
EXEMPT = frozenset({"NA", "INF", "-INF", "inf", "-inf", "nan", ""})
INTEGER = re.compile(r"^-?[0-9]+$")


def sig_digits(cell: str) -> int:
    mantissa = cell.strip().lstrip("+-").lower().split("e")[0].replace(".", "")
    return len(mantissa.lstrip("0").rstrip("0"))


def _is_float(cell: str) -> bool:
    try:
        float(cell)
    except ValueError:
        return False
    return True


def _cells(path: Path) -> tuple[list[str], list[list[str]]]:
    lines = path.read_text(encoding="utf-8").split("\n")
    if lines and lines[-1] == "":
        lines.pop()
    return lines[0].split("\t"), [ln.split("\t") for ln in lines[1:]]


@pytest.fixture(scope="module")
def out(tmp_path_factory: pytest.TempPathFactory) -> Path:
    from helpers.mapping_steps import write_test_db

    import gippasim
    from gippa.cli import cli
    from gippasim.proteomics import ProteomicsConfig

    cfg = dataclasses.replace(
        gippasim.SIM100K, proteomics=ProteomicsConfig(), write_precomputed=True
    )
    sim = gippasim.simulate(cfg, gippasim.SIM100K_SEED, tmp_path_factory.mktemp("sim") / "sim")
    db = Path(write_test_db(tmp_path_factory.mktemp("test_db_pr9")))
    outdir = tmp_path_factory.mktemp("LOCAL") / "OUT"
    args = ["run", "--genome", str(sim.genome_fasta), "--peptides", str(sim.peptides),
            "--abundance", str(sim.abundance), "--annotation", str(sim.bakta_json),
            "--config", str(sim.config_yaml), "--outdir", str(outdir), "--modules", "core",
            "--seed", str(SEED), "--executor", "local"]  # fmt: skip
    with pytest.MonkeyPatch.context() as mp:
        mp.setenv("GIPPA_DB_DIR", str(db))
        res = CliRunner().invoke(cli, args)
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        raise res.exception
    assert res.exit_code == 0, res.output
    return outdir


def _published(out: Path) -> list[Path]:
    return sorted([*(out / "tables").glob("*.tsv"), *(out / "mapping").glob("*.tsv")])


def test_published_tables_have_schemas(out: Path) -> None:
    paths = _published(out)
    assert any(p.name == "detection_summary.tsv" for p in paths)
    assert any(p.name == "regions_stats.tsv" for p in paths)
    missing = [p.name for p in paths if not (TABLE_DIR / f"{p.stem}.schema.json").is_file()]
    assert not missing, f"published TSVs without a table schema: {missing}"


def test_ac_9_1_number_columns_follow_6g(out: Path) -> None:
    bad: dict[str, list[tuple[int, str, str]]] = {}
    for path in _published(out):
        schema = TABLE_DIR / f"{path.stem}.schema.json"
        if not schema.is_file():
            continue
        fields = json.loads(schema.read_text(encoding="utf-8"))["fields"]
        nums = {f["name"] for f in fields if f["type"] == "number"}
        header, rows = _cells(path)
        for i, row in enumerate(rows, start=1):
            for col, cell in zip(header, row, strict=True):
                if col in nums and cell not in EXEMPT and sig_digits(cell) > 6:
                    bad.setdefault(f"{path.parent.name}/{path.name}", []).append((i, col, cell))
    summary = {k: (len(v), v[:2]) for k, v in bad.items()}
    assert not bad, f"number cells over %.6g: {summary}"


def test_ac_9_1_plot_data_follow_6g(out: Path) -> None:
    data = sorted((out / "plots" / "data").glob("*.tsv"))
    assert data, "the run wrote no plots/data/*.tsv"
    bad: dict[str, list[tuple[int, str, str]]] = {}
    for path in data:
        header, rows = _cells(path)
        for i, row in enumerate(rows, start=1):
            for col, cell in zip(header, row, strict=True):
                if cell in EXEMPT or INTEGER.match(cell) or not _is_float(cell):
                    continue
                if sig_digits(cell) > 6:
                    bad.setdefault(path.name, []).append((i, col, cell))
    summary: dict[str, Any] = {k: (len(v), v[:2]) for k, v in bad.items()}
    assert not bad, f"plot data cells over %.6g: {summary}"
