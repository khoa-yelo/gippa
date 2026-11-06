"""``gippa convert`` creates missing parent folders of ``--out`` and ``--abundance-out``.

M12-01 dry run (jobs 1939608, 1939670): the example's printed commands write to
``gippa-example/runs/<route>/peptides.tsv`` before any run folder exists; convert crashed with
exit 5 ("internal error") instead of writing the table.
"""

from __future__ import annotations

from pathlib import Path

from click.testing import CliRunner
from test_convert_review import ROW, _tsv  # type: ignore[import-not-found]

from gippa.cli import cli
from gippa.exitcodes import ExitCode


def test_out_in_a_missing_folder_is_created(tmp_path: Path) -> None:
    src = _tsv(tmp_path, [ROW])
    out = tmp_path / "runs" / "A" / "peptides.tsv"
    result = CliRunner().invoke(cli, ["convert", "--from", "diann", str(src), "--out", str(out)])
    assert result.exit_code == ExitCode.SUCCESS, result.output
    assert out.read_text().startswith("peptide\t")


def test_abundance_out_in_a_missing_folder_is_created(tmp_path: Path) -> None:
    src = _tsv(tmp_path, [ROW])
    out = tmp_path / "peptides.tsv"
    ab = tmp_path / "abundance" / "D" / "abundance.tsv"
    argv = ["convert", "--from", "diann", str(src), "--out", str(out), "--abundance-out", str(ab)]
    result = CliRunner().invoke(cli, argv)
    assert result.exit_code == ExitCode.SUCCESS, result.output
    assert ab.is_file()
