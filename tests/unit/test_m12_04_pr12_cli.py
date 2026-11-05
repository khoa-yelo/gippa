"""M12-04 PR 12 (O-01) through the CLI: ``gippa run`` writes ``OUT/README.txt`` (also after a failed
run), its manifest reads the moved intermediates, ``--resume`` into a 1.0 directory warns once,
and ``gippa report`` / ``gippa compare`` read both layouts (docs/specs/M12-04-outputs-revision.md
> PR 12 "Reads", "README", edge cases; AC 12.1, 12.3, 12.4; § 5.5 "gippa report and gippa
compare still read 1.0 output folders").

The fake ``nextflow`` (as tests/unit/manifest/test_manifest_cli_nextflow.py) exits
``$FAKE_NF_RC``, copies ``$FAKE_NF_TRACE`` to the ``-with-trace`` path and copies the tree
``$FAKE_NF_FILES`` into ``OUT`` (what GIPPA_STEP's publishDir would have published). Nothing real
is launched; ``$GIPPA_ROOT`` and ``$GIPPA_DB_DIR`` are fakes from ``helpers.fake_env``.

Interface pinned here (test author):
- the legacy warning is one ``OUT/logs/events.jsonl`` line at level ``WARNING`` whose
  ``event`` (or ``<step>.<event>``) is ``outputs.legacy_intermediates`` and whose ``data`` holds the
  count as an integer value (the fixture has 4 legacy top-level entries holding 4 files, so a
  count of entries or of files both read 4);
- ``README.txt`` equals ``gippa.io.out_readme.render_readme(OUT)`` of the finished tree.
"""

from __future__ import annotations

import json
import math
import os
import random
import stat
from pathlib import Path
from typing import Any

import pyarrow as pa  # type: ignore[import-untyped]
import pyarrow.parquet as pq  # type: ignore[import-untyped]
import pytest
from click.testing import CliRunner, Result
from helpers.fake_env import fake_gippa_root, write_setup_manifest

from gippa.exitcodes import ExitCode

TOP_LEVEL = frozenset({
    "README.txt", "run_manifest.json", "resolved_config.yaml", "preflight", "annotation",
    "mapping", "tables", "plots", "report", "schemas", "logs", "intermediate", "work",
})  # fmt: skip
FAKE_NF = """#!/bin/sh
printf '%s\\n' "$@" > "$FAKE_NF_ARGS"
if [ -n "$FAKE_NF_TRACE" ]; then
  prev=""
  for a in "$@"; do
    if [ "$prev" = "-with-trace" ]; then mkdir -p "$(dirname "$a")"; cp "$FAKE_NF_TRACE" "$a"; fi
    prev="$a"
  done
fi
if [ -n "$FAKE_NF_FILES" ]; then cp -R "$FAKE_NF_FILES"/. .; fi
exit "${FAKE_NF_RC:-0}"
"""
COLS = ["task_id", "hash", "native_id", "name", "tag", "status", "exit", "attempt", "realtime",
        "%cpu", "peak_rss", "workdir", "container"]  # fmt: skip
DIAG = [("chr1", True, "auto", 2, 5, 0.33, []),
        ("pA", False, "not_scanned", None, None, math.nan, ["few_values"])]  # fmt: skip
#: 1.0 top-level intermediates (4 entries, 4 files)
LEGACY = ("codon.parquet", "genes.coords", "genetable/genes.parquet",
          "context/gene_context.parquet")  # fmt: skip
LEGACY_EVENT = "outputs.legacy_intermediates"


def _write_diag(path: Path) -> None:
    schema = pa.schema([("contig", pa.string()), ("scanned", pa.bool_()),
                        ("block_rule", pa.string()), ("l_hat", pa.int64()),
                        ("block_len", pa.int64()), ("morans_i_lag1", pa.float64()),
                        ("flags", pa.list_(pa.string()))])  # fmt: skip
    cols = list(zip(*DIAG, strict=True))
    path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(pa.Table.from_arrays([pa.array(list(c), type=f.type) for c, f in
                                         zip(cols, schema, strict=True)], schema=schema),
                   path)  # fmt: skip


@pytest.fixture
def nf(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> dict[str, Path]:
    monkeypatch.setenv("GIPPA_ROOT", str(fake_gippa_root(tmp_path / "gippa_root")))
    db = tmp_path / "db"
    db.mkdir()
    write_setup_manifest(db)
    monkeypatch.setenv("GIPPA_DB_DIR", str(db))
    for var in ("SLURM_CPUS_ON_NODE", "SLURM_CPUS_PER_TASK", "SLURM_MEM_PER_NODE"):
        monkeypatch.delenv(var, raising=False)
    bindir = tmp_path / "fakebin"
    bindir.mkdir()
    exe = bindir / "nextflow"
    exe.write_text(FAKE_NF)
    exe.chmod(exe.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    monkeypatch.setenv("PATH", f"{bindir}{os.pathsep}{os.environ.get('PATH', '')}")
    monkeypatch.setenv("FAKE_NF_ARGS", str(tmp_path / "nf_args.txt"))
    for var in ("FAKE_NF_RC", "FAKE_NF_TRACE", "FAKE_NF_FILES"):
        monkeypatch.delenv(var, raising=False)
    rng = random.Random(3)
    genome = tmp_path / "in" / "genome.fa"
    genome.parent.mkdir()
    seq = "".join(rng.choice("ACGT") for _ in range(6000))
    genome.write_text(">chr1 [topology=circular]\n" + "\n".join(seq[i:i + 60]
                      for i in range(0, 6000, 60)) + "\n")  # fmt: skip
    aa = "ACDEFGHIKLMNPQRSTVWY"
    proteins = tmp_path / "in" / "proteins.faa"
    proteins.write_text("".join(f">prot{i}\nM{''.join(rng.choice(aa) for _ in range(150))}\n"
                                for i in range(3)))  # fmt: skip
    return {"genome": genome, "proteins": proteins}


def _published(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, diag_rel: str) -> None:
    """What the GIPPA_STEP publishDir writes: the island diagnostics at ``diag_rel`` and one
    table."""
    files = tmp_path / "published"
    _write_diag(files / diag_rel)
    (files / "tables").mkdir(parents=True)
    (files / "tables" / "sensitivity.tsv").write_text(
        "family\tvariant\tn_primary\tn_variant\tn_shared\tjaccard\tagreement\n"
    )
    monkeypatch.setenv("FAKE_NF_FILES", str(files))


def _run(nf: dict[str, Path], out: Path, *extra: str) -> tuple[int, str]:
    from gippa.cli import cli

    argv = ["run", "--genome", str(nf["genome"]), "--proteins", str(nf["proteins"]),
            "--outdir", str(out), *extra]  # fmt: skip
    res = CliRunner().invoke(cli, argv)
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        raise res.exception
    return res.exit_code, res.output + (res.stderr if res.stderr_bytes is not None else "")


def _top(out: Path) -> set[str]:
    return {p.name for p in out.iterdir() if not p.name.startswith(".")}


def _events(out: Path) -> list[dict[str, Any]]:
    path = out / "logs" / "events.jsonl"
    return [json.loads(ln) for ln in path.read_text(encoding="utf-8").splitlines() if ln.strip()]


def _legacy_events(out: Path) -> list[dict[str, Any]]:
    return [e for e in _events(out)
            if LEGACY_EVENT in (e.get("event"), f"{e.get('step')}.{e.get('event')}")]  # fmt: skip


def _replicons_want() -> list[dict[str, Any]]:
    names = ["contig", "scanned", "block_rule", "l_hat", "block_len", "morans_i_lag1", "flags"]
    return [{k: (None if isinstance(v, float) and math.isnan(v) else v)
             for k, v in zip(names, r, strict=True)} for r in DIAG]  # fmt: skip


def _legacy_tree(out: Path) -> None:
    for rel in LEGACY:
        (out / rel).parent.mkdir(parents=True, exist_ok=True)
        (out / rel).write_bytes(b"legacy\n")


# ---------------------------------------------------------------- README after gippa run


def test_run_writes_the_readme_and_only_published_top_level_names(
    nf: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from gippa.io.out_readme import render_readme

    _published(tmp_path, monkeypatch, "intermediate/context/island_diagnostics.parquet")
    out = tmp_path / "OUT"
    code, log = _run(nf, out)
    assert code == ExitCode.SUCCESS, log
    readme = out / "README.txt"
    assert readme.is_file(), f"no README.txt; OUT has {sorted(_top(out))}"
    assert readme.read_text(encoding="ascii") == render_readme(out)
    assert _top(out) <= TOP_LEVEL, sorted(_top(out) - TOP_LEVEL)


def test_failed_run_still_writes_the_readme(
    nf: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from gippa.io.out_readme import render_readme

    monkeypatch.setenv("FAKE_NF_RC", "1")
    out = tmp_path / "OUT"
    code, log = _run(nf, out)
    assert code == ExitCode.STEP_FAILED, log
    readme = out / "README.txt"
    assert readme.is_file(), f"a failed run wrote no README.txt; OUT has {sorted(_top(out))}"
    assert readme.read_text(encoding="ascii") == render_readme(out)


# ---------------------------------------------------------------- the manifest reads intermediate/


def test_manifest_context_reads_intermediate_island_diagnostics(
    nf: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _published(tmp_path, monkeypatch, "intermediate/context/island_diagnostics.parquet")
    out = tmp_path / "OUT"
    code, log = _run(nf, out)
    assert code == ExitCode.SUCCESS, log
    m = json.loads((out / "run_manifest.json").read_text(encoding="utf-8"))
    assert m["context"]["replicons"] == _replicons_want()


def test_manifest_context_still_reads_a_1_0_location(
    nf: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _published(tmp_path, monkeypatch, "context/island_diagnostics.parquet")
    out = tmp_path / "OUT"
    code, log = _run(nf, out, "--resume")
    assert code == ExitCode.SUCCESS, log
    m = json.loads((out / "run_manifest.json").read_text(encoding="utf-8"))
    assert m["context"]["replicons"] == _replicons_want()


# ---------------------------------------------------------------- --resume into a 1.0 directory


def test_resume_into_a_1_0_directory_keeps_the_old_files_and_warns_once(
    nf: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    out = tmp_path / "OUT"
    _legacy_tree(out)
    _published(tmp_path, monkeypatch, "intermediate/context/island_diagnostics.parquet")
    code, log = _run(nf, out, "--resume")
    assert code == ExitCode.SUCCESS, log
    for rel in LEGACY:
        assert (out / rel).read_bytes() == b"legacy\n", f"{rel} was removed or rewritten"
    found = _legacy_events(out)
    assert len(found) == 1, f"{LEGACY_EVENT} events: {found}"
    (ev,) = found
    assert ev["level"] == "WARNING", ev
    counts = [v for v in ev["data"].values() if isinstance(v, int) and not isinstance(v, bool)]
    assert len(LEGACY) in counts, f"no count {len(LEGACY)} in {ev['data']}"


def test_resume_into_a_1_1_directory_does_not_warn(
    nf: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    out = tmp_path / "OUT"
    (out / "intermediate" / "genetable").mkdir(parents=True)
    (out / "intermediate" / "codon.parquet").write_bytes(b"new\n")
    _published(tmp_path, monkeypatch, "intermediate/context/island_diagnostics.parquet")
    code, log = _run(nf, out, "--resume")
    assert code == ExitCode.SUCCESS, log
    assert _legacy_events(out) == []


def test_fresh_run_does_not_warn(
    nf: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _published(tmp_path, monkeypatch, "intermediate/context/island_diagnostics.parquet")
    out = tmp_path / "OUT"
    code, log = _run(nf, out)
    assert code == ExitCode.SUCCESS, log
    assert _legacy_events(out) == []


# ---------------------------------------------------------------- gippa report on both layouts


FAST = ["--set", "report.formats=[svg]", "--set", "report.dpi=50"]


def _report_fixture(tmp_path: Path, layout: str, monkeypatch: pytest.MonkeyPatch) -> Path:
    from helpers.report_content import write_diagnostics
    from helpers.report_fixtures import register_fake, run_out, swap_registry

    (tmp_path / "db_empty").mkdir()
    monkeypatch.setenv("GIPPA_DB_DIR", str(tmp_path / "db_empty"))
    monkeypatch.setenv("GIPPA_ROOT", str(tmp_path / "db_empty"))
    swap_registry(monkeypatch)
    register_fake("V07")
    out = run_out(tmp_path, islands=3)
    rows = [{"contig": "chr1", "scanned": True, "block_rule": "auto", "l_hat": 4,
             "block_len": 9, "morans_i_lag1": 0.4213}]  # fmt: skip
    path = write_diagnostics(out, rows)  # the 1.0 location, OUT/context/
    if layout == "1.1":
        dst = out / "intermediate" / "context" / path.name
        dst.parent.mkdir(parents=True)
        path.rename(dst)
        path.parent.rmdir()
    return out


def _report(out: Path) -> Result:
    from gippa.cli import cli

    res = CliRunner().invoke(cli, ["report", str(out), *FAST])
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        raise res.exception
    return res


@pytest.mark.parametrize("layout", ["1.0", "1.1"])
def test_report_reads_island_diagnostics_in_either_layout(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, layout: str
) -> None:
    """AC 12.4: ``gippa report`` succeeds on a 1.0 tree and on a 1.1 tree, and finds the
    island diagnostics (an intermediate) in both."""
    out = _report_fixture(tmp_path, layout, monkeypatch)
    res = _report(out)
    assert res.exit_code == 0, res.output
    html = (out / "report" / "report.html").read_text(encoding="utf-8")
    assert 'data-table="island_diagnostics"' in html, f"{layout}: island diagnostics not read"


# ---------------------------------------------------------------- gippa compare on both layouts


def test_compare_reads_a_1_0_and_a_1_1_run_together(tmp_path: Path) -> None:
    """§ 5.5: ``gippa compare`` still reads 1.0 output folders, also next to a 1.1 one."""
    import numpy as np
    from helpers.compare_data import gene, write_run

    from gippa.cli import cli

    def genes(shift: float) -> list[dict[str, object]]:
        vals = np.sin(np.arange(8) * 1.7) * 3 + np.arange(8) * 0.4
        return [gene(f"g{i:02d}", start=1000 * i, log2=float(v + shift * (i % 3)),
                     tu=f"TU_{i // 2}") for i, v in enumerate(vals)]  # fmt: skip

    old, new = write_run(tmp_path / "run10", genes(0.0)), write_run(tmp_path / "run11", genes(0.5))
    _legacy_tree(old)
    for rel in LEGACY:
        (new / "intermediate" / rel).parent.mkdir(parents=True, exist_ok=True)
        (new / "intermediate" / rel).write_bytes(b"new\n")
    args = [str(old), str(new), "--outdir", str(tmp_path / "CMP"), "--labels", "A,B",
            "--set", "stats.bootstrap=50"]  # fmt: skip
    res = CliRunner().invoke(cli, ["compare", *args])
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        raise res.exception
    assert res.exit_code == 0, res.output
    for name in ("compare_detection", "compare_abundance", "compare_stability"):
        assert (tmp_path / "CMP" / "tables" / f"{name}.tsv").is_file(), name


# ---------------------------------------------------------------- gippa demo (local executor)


@pytest.mark.slow
def test_demo_layout_and_readme(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Edge cases "The local developer executor" and "gippa demo output": same layout."""
    from gippa.cli import cli

    monkeypatch.delenv("GIPPA_DB_DIR", raising=False)
    out = tmp_path / "OUT"
    res = CliRunner().invoke(cli, ["demo", "--outdir", str(out)])
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        raise res.exception
    assert res.exit_code == 0, res.output
    assert _top(out) <= TOP_LEVEL, sorted(_top(out) - TOP_LEVEL)
    assert (out / "intermediate" / "genetable" / "genes.parquet").is_file()
    assert (out / "tables" / "genes.tsv").is_file()
    assert (out / "report" / "report.html").is_file()
    from gippa.io.out_readme import render_readme

    assert (out / "README.txt").read_text(encoding="ascii") == render_readme(out)
