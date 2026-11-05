"""M12-04 PR 11a: a ``genes.tsv`` in the 1.0 column order stays readable (Q3 default "yes";
AC 11a.4; docs/specs/M12-04-outputs-revision.md § PR 11a edge cases, § 5.5 "``gippa report`` and
``gippa compare`` still read 1.0 output folders").

- ``RunView.table("genes_tsv")`` reads columns by name: a 1.0-ordered file gives the same frame,
  in the 1.1 (schema) column order, as the same rows written in the 1.1 order.
- AC 11a.4: ``gippa compare`` of a 1.0-ordered run against a 1.1-ordered one succeeds and writes
  the same ``compare_detection.tsv`` as two 1.1-ordered runs.
- ``gippa report`` on an ``OUT/`` whose ``genes.tsv`` is 1.0-ordered succeeds and writes the same
  report as on the 1.1-ordered file (fake plot registry, tests/helpers/report_fixtures.py).

The 1.1 order is ``ORDER_1_1`` (§ 5.2) from test_m12_04_pr11a_genes_order.py; the files are
written by hand under an explicit header, so these tests do not depend on the writer.
"""

from __future__ import annotations

import re
from collections.abc import Mapping, Sequence
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from click.testing import CliRunner, Result
from helpers.compare_data import SHA_A, gene, write_preflight
from helpers.report_fixtures import swap_registry
from unit.test_m12_04_pr11a_genes_order import ORDER_1_0, ORDER_1_1

from gippa.exitcodes import ExitCode

N = 12
FAST_CMP = ["--set", "stats.bootstrap=200"]
FAST_REPORT = ["--set", "report.formats=[svg]", "--set", "report.dpi=50"]


@pytest.fixture(autouse=True)
def _no_db(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    (tmp_path / "db_empty").mkdir()
    monkeypatch.setenv("GIPPA_DB_DIR", str(tmp_path / "db_empty"))
    monkeypatch.setenv("GIPPA_ROOT", str(tmp_path / "db_empty"))


def _genes(shift: float, detected_upto: int = N) -> list[dict[str, object]]:
    vals = np.sin(np.arange(N) * 1.7) * 3 + np.arange(N) * 0.4
    return [
        gene(f"g{i:02d}", start=1000 * i, detected=i < detected_upto,
             log2=float(v + shift * (i % 3)), tu=f"TU_{i // 2}", detection_tier="mayu",
             annotation_source="bakta_light", flags="partial" if i == 3 else None)
        for i, v in enumerate(vals)
    ]  # fmt: skip


def _cell(v: object) -> str:
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return "NA"
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, float):
        return f"{v:.6g}"
    return str(v)


def _write_genes(out: Path, rows: Sequence[Mapping[str, object]], order: Sequence[str]) -> Path:
    """``OUT/tables/genes.tsv`` under the header ``order`` (1-based ``start`` on disk)."""
    path = out / "tables" / "genes.tsv"
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["\t".join(order)]
    for r in rows:
        row = dict(r)
        row["start"] = int(row["start"]) + 1  # type: ignore[call-overload]
        lines.append("\t".join(_cell(row.get(c)) for c in order))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _run(out: Path, rows: Sequence[Mapping[str, object]], order: Sequence[str]) -> Path:
    _write_genes(out, rows, order)
    write_preflight(out, SHA_A)
    return out


def test_the_two_orders_differ() -> None:
    assert ORDER_1_0 != ORDER_1_1 and set(ORDER_1_0) == set(ORDER_1_1)


# ----------------------------------------------------------------------------- RunView


def test_run_view_reads_a_1_0_ordered_file_in_the_1_1_order(tmp_path: Path) -> None:
    from gippa.io.run_view import RunView

    old = _run(tmp_path / "old" / "OUT", _genes(0.0), ORDER_1_0)
    new = _run(tmp_path / "new" / "OUT", _genes(0.0), ORDER_1_1)
    a = RunView.from_outdir(old).table("genes_tsv")
    b = RunView.from_outdir(new).table("genes_tsv")
    assert a is not None and b is not None
    assert tuple(a.columns) == ORDER_1_1
    pd.testing.assert_frame_equal(a, b)


# ----------------------------------------------------------------------------- AC 11a.4


def _compare(*args: str) -> Result:
    from gippa.cli import cli

    res = CliRunner().invoke(cli, ["compare", *args])
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        raise res.exception
    return res


def _detection(tmp: Path, order_a: Sequence[str], order_b: Sequence[str]) -> bytes:
    a = _run(tmp / "runA", _genes(0.0), order_a)
    b = _run(tmp / "runB", _genes(0.5, detected_upto=10), order_b)
    cmp = tmp / "CMP"
    res = _compare(str(a), str(b), "--outdir", str(cmp), "--labels", "A1,A2", *FAST_CMP)
    assert res.exit_code == ExitCode.SUCCESS, res.output
    return (cmp / "tables" / "compare_detection.tsv").read_bytes()


def test_compare_1_0_against_1_1_equals_two_1_1_runs(tmp_path: Path) -> None:
    mixed = _detection(tmp_path / "mixed", ORDER_1_0, ORDER_1_1)
    both_new = _detection(tmp_path / "new", ORDER_1_1, ORDER_1_1)
    assert mixed == both_new


def test_compare_two_1_0_runs_equals_two_1_1_runs(tmp_path: Path) -> None:
    assert _detection(tmp_path / "old", ORDER_1_0, ORDER_1_0) == _detection(
        tmp_path / "new", ORDER_1_1, ORDER_1_1
    )


# ----------------------------------------------------------------------------- gippa report


def _report(out: Path) -> Result:
    from gippa.cli import cli

    res = CliRunner().invoke(cli, ["report", str(out), *FAST_REPORT])
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        raise res.exception
    return res


def _report_html(tmp: Path, order: Sequence[str]) -> str:
    out = _run(tmp / "OUT", _genes(0.0), order)
    res = _report(out)
    assert res.exit_code == ExitCode.SUCCESS, res.output
    html = (out / "report" / "report.html").read_text(encoding="utf-8")
    return html.replace(str(tmp), "<TMP>")


def _without_tables(html: str) -> str:
    return re.sub(r"<table\b.*?</table>", "<table/>", html, flags=re.S)


def test_report_on_a_1_0_ordered_genes_tsv_matches_the_1_1_report(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Cards, notes and sections agree; only the raw file previews (which show each file's own
    header, a by-design copy of the file head) may differ in column order."""
    swap_registry(monkeypatch)
    old = _report_html(tmp_path / "old", ORDER_1_0)
    new = _report_html(tmp_path / "new", ORDER_1_1)
    assert "genes detected" in new
    assert _without_tables(old) == _without_tables(new)
