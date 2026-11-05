# ruff: noqa: E501  (EXAMPLE_A is the spec README verbatim; one of its lines is 101 characters)
"""M12-04 PR 12 (O-01): ``OUT/README.txt`` (docs/specs/M12-04-outputs-revision.md > PR 12 "README",
"Example: A's top level", edge case "A failed run", AC 12.3; § 1 rules 2 and 3; D-01 docs half).

Interface pinned here (test author; the spec says only "written … by ``gippa.io``"):
- ``gippa.io.out_readme.render_readme(out: Path) -> str``: the README text for the entries that
  exist in ``out`` (hidden entries and ``README.txt`` itself are ignored);
- ``gippa.io.out_readme.write_readme(out: Path) -> Path``: writes ``out / "README.txt"`` (that
  text, UTF-8 = ASCII) and returns its path.
An entry line is two spaces, the entry (``name/`` for a folder, the file name for a file, or a
path below one for "Start here"), two or more spaces and its sentence, as in the spec example.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

import gippa
from gippa.io.run_view import PATHS
from gippa.runner.local import PUBLISHED

ROOT = Path(__file__).resolve().parents[2]
FOLDERS = ("tables", "mapping", "annotation", "plots", "report", "preflight", "schemas", "logs",
           "intermediate", "work")  # fmt: skip
FILES = ("run_manifest.json", "resolved_config.yaml")
ENTRY = re.compile(r"^  (\S+)\s{2,}\S")

#: the spec's README for run A (every entry present); the version line follows gippa.__version__
EXAMPLE_A = """GIPPA results

Start here:
  report/report.html    the report: open it in a web browser
  tables/genes.tsv      one row per gene: detection, abundance, function, context, features

Folders:
  tables/        result tables (TSV); every column is described in schemas/datapackage.json
  mapping/       where each peptide or protein was placed in the genome, and what was not placed
  annotation/    the gene model and annotation layers this run used; tools/ holds native tool outputs
  plots/         figures (SVG, PNG, PDF); plots/data/ holds the numbers behind each figure
  report/        report.html, self-contained
  preflight/     the input checks made before the analysis
  schemas/       a machine-readable description of every table (Frictionless Data Package)
  logs/          step logs, events and timings
  intermediate/  files passed between steps; not results; their layout may change between versions
  work/          the workflow engine's working files; needed only to resume this run

Files:
  run_manifest.json     versions, databases, parameters, input checksums and seeds
  resolved_config.yaml  the configuration this run used

Coordinates: TSV, GFF3 and VCF are 1-based closed; BED and proBed are 0-based half-open.
Table reference: https://khoa-yelo.github.io/gippa/reference/tables/
Written by GIPPA {version}.
"""


def _render(out: Path) -> str:
    from gippa.io.out_readme import render_readme

    return render_readme(out)


def _tree(
    out: Path, folders: tuple[str, ...] | set[str], files: tuple[str, ...] | set[str]
) -> Path:
    out.mkdir(parents=True, exist_ok=True)
    for d in folders:
        (out / d).mkdir(exist_ok=True)
    for f in files:
        (out / f).write_text("x\n")
    if "report" in folders:
        (out / "report" / "report.html").write_text("<html></html>\n")
    if "tables" in folders:
        (out / "tables" / "genes.tsv").write_text("locus_tag\n")
    return out


def _full(out: Path) -> Path:
    return _tree(out, FOLDERS, FILES)


def _named_heads(text: str) -> set[str]:
    """Top-level names the README's entry lines name (``tables/genes.tsv`` names ``tables``)."""
    heads = set()
    for line in text.splitlines():
        m = ENTRY.match(line)
        if m:
            heads.add(m.group(1).split("/")[0])
    return heads


# ---------------------------------------------------------------- the spec example


def test_example_a_text(tmp_path: Path) -> None:
    out = _full(tmp_path / "OUT")
    assert _render(out) == EXAMPLE_A.replace("{version}", gippa.__version__)


def test_write_readme_writes_the_rendered_text(tmp_path: Path) -> None:
    from gippa.io.out_readme import write_readme

    out = _full(tmp_path / "OUT")
    path = write_readme(out)
    assert path == out / "README.txt"
    assert path.read_bytes() == _render(out).encode("ascii")
    assert _render(out) == path.read_text(encoding="ascii"), "README.txt names itself"


def test_folders_and_files_in_the_fixed_order(tmp_path: Path) -> None:
    text = _render(_full(tmp_path / "OUT"))
    folders = text.split("Folders:\n", 1)[1].split("\n\n", 1)[0]
    files = text.split("Files:\n", 1)[1].split("\n\n", 1)[0]
    assert [ENTRY.match(ln).group(1) for ln in folders.splitlines()] == [f"{d}/" for d in FOLDERS]  # type: ignore[union-attr]
    assert [ENTRY.match(ln).group(1) for ln in files.splitlines()] == list(FILES)  # type: ignore[union-attr]


def test_ends_with_the_docs_link_and_the_version(tmp_path: Path) -> None:
    lines = _render(_full(tmp_path / "OUT")).rstrip("\n").splitlines()
    assert lines[-2] == "Table reference: https://khoa-yelo.github.io/gippa/reference/tables/"
    assert lines[-1] == f"Written by GIPPA {gippa.__version__}."


# ---------------------------------------------------------------- AC 12.3


#: top-level names GIPPA 1.0 published at the top of ``OUT/`` that 1.1 moves under
#: ``intermediate/`` (review round 1, M2): the first path parts of the 1.0 Nextflow ``key=path``
#: outputs (``git show v1.0.0:workflow/**/*.nf``, so also of 1.0 ``PUBLISHED_FILES``) outside the
#: published folders, ``work/`` and the plot support directories, plus the spec's "Example: A's
#: top level" list. ``legacy_entries`` must treat each of these, and no user file, as legacy.
LEGACY_1_0 = (
    "annotation_meta.json", "codon.parquet", "context", "features.parquet", "genes.coords",
    "genes.parquet", "genetable", "hits_exact.parquet", "hits_tolerant.parquet",
    "layers.parquet", "mapping_fdr.parquet", "novel", "novel_candidates.parquet",
    "novel_promoter_windows.parquet", "novel_windows", "novel_windows_map.parquet",
    "operon_meta.json", "ostir_input.csv", "ostir_input_novel.csv", "pairs_raw.parquet",
    "placement.parquet", "placement_costs.parquet", "promoter_windows.parquet",
    "proteins.leftover.faa", "rbs.parquet", "rbs_novel.parquet", "region_genes.parquet",
    "regions.parquet", "replicons.parquet", "rrna16s.fa", "sixframe.npz", "smoke.done",
    "stratum0_classes.parquet", "terminators.parquet", "terminators_crosscheck.tsv",
    "trna_counts.tsv", "tu_upstream", "tus.parquet", "tx_rates.parquet", "uniop_cds.gff",
    "uniop_map.tsv", "windows", "windows_map.parquet", "xrefs.parquet",
)  # fmt: skip
#: 1.0 names that were folders (the rest were files)
LEGACY_DIRS = frozenset({"context", "genetable", "novel", "novel_windows", "tu_upstream",
                         "windows"})  # fmt: skip
#: a user's own files next to the results: named in the README, never "can be removed"
FOREIGN = ("samplesheet.csv", "notes.txt", "slurm-123.out")
#: non-ASCII top-level names (a user's file); the README stays ASCII and names them escaped
NON_ASCII = ("r\u00e9sum\u00e9.txt", "\u00fcbersicht", "\u6570\u636e.csv", "lambda_\u03bb.tsv",
             "smile_\U0001f600.txt")  # fmt: skip
REMOVABLE = "can be removed"


def _ascii(name: str) -> str:
    """How the README writes a top-level name: backslash-escaped outside ASCII."""
    return name.encode("ascii", "backslashreplace").decode("ascii")


def _names(text: str, name: str, is_dir: bool) -> bool:
    """``text`` names the entry ``name`` as a whole token (``name/`` for a folder), so
    ``genes.parquet`` is not found in ``genes.parquet.bak`` nor ``nes.tsv`` in ``genes.tsv``."""
    token = re.escape(_ascii(name) + ("/" if is_dir else ""))
    after = r"(?=$|[\s,;:)]|\.(?:\s|$))" if is_dir else r"(?=$|[\s,;:)/]|\.(?:\s|$))"
    return re.search(r"(?<![\w./\\-])" + token + after, text, re.MULTILINE) is not None


def _paragraphs_naming(text: str, name: str, is_dir: bool) -> list[str]:
    return [p for p in text.split("\n\n") if _names(p, name, is_dir)]


def _extras(out: Path, dirs: set[str], files: set[str]) -> None:
    for d in dirs:
        (out / d).mkdir()
        (out / d / "part.parquet").write_bytes(b"x")
    for f in files:
        (out / f).write_bytes(b"x")


_FOREIGN_NAME = st.one_of(
    st.sampled_from(FOREIGN + NON_ASCII),
    st.from_regex(r"\Amy_[a-z0-9]{2,8}(\.(txt|csv|out|log|tsv))?\Z", fullmatch=True),
)


@settings(
    max_examples=60, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture]
)
@given(
    st.sets(st.sampled_from(FOLDERS)),
    st.sets(st.sampled_from(FILES)),
    st.sets(st.sampled_from(LEGACY_1_0)),
    st.sets(_FOREIGN_NAME, max_size=6),
    st.data(),
)
def test_names_exactly_the_entries_that_exist(
    tmp_path_factory: pytest.TempPathFactory,
    folders: set[str],
    files: set[str],
    legacy: set[str],
    foreign: set[str],
    data: st.DataObject,
) -> None:
    """Property (AC 12.3): the README names every top-level entry present (except itself) and
    nothing absent, whatever a run (or a failed run, or a ``--resume`` into a 1.0 folder holding a
    user's own files) left: published entries on their own lines, 1.0 step intermediates and
    foreign entries by name (non-ASCII names backslash-escaped, the text stays ASCII); only the
    1.0 intermediates are called removable; ≤ 60 lines."""
    from gippa.io.out_readme import legacy_entries

    out = _tree(tmp_path_factory.mktemp("rd") / "OUT", folders, files)
    (out / "README.txt").write_text("stale\n")
    (out / ".nextflow.log").write_text("hidden\n")
    foreign_dirs = {f for f in foreign if data.draw(st.booleans(), label=f"{f} is a folder")}
    _extras(out, (legacy & LEGACY_DIRS) | foreign_dirs,
            (legacy - LEGACY_DIRS) | (foreign - foreign_dirs))  # fmt: skip
    text = _render(out)
    assert text.isascii()
    assert len(text.splitlines()) <= 60, f"{len(text.splitlines())} lines"
    published = set(FOLDERS) | set(FILES)
    assert _named_heads(text) & published == folders | files
    assert sorted(legacy_entries(out)) == sorted(legacy), "legacy_entries: 1.0 names only"
    for name in legacy:
        is_dir = name in LEGACY_DIRS
        assert _names(text, name, is_dir), f"1.0 entry {name!r} not named:\n{text}"
        assert any(REMOVABLE in p for p in _paragraphs_naming(text, name, is_dir)), (
            f"1.0 entry {name!r} not called removable:\n{text}"
        )
    for name in foreign:
        is_dir = name in foreign_dirs
        assert _names(text, name, is_dir), f"entry {_ascii(name)!r} not named:\n{text}"
        assert all("removed" not in p for p in _paragraphs_naming(text, name, is_dir)), (
            f"a user's entry {_ascii(name)!r} is called removable:\n{text}"
        )
    for name in set(LEGACY_1_0) - legacy:
        is_dir = name in LEGACY_DIRS
        assert not _names(text, name, is_dir), f"absent entry {name!r} named:\n{text}"


def test_fifty_legacy_entries_all_named_within_60_lines(tmp_path: Path) -> None:
    """H1: a ``--resume`` into a 1.0 folder with ~50 top-level intermediates still gives a README
    that names each of them, in at most 60 lines."""
    out = _full(tmp_path / "OUT")
    extra = tuple(f"step{i:02d}.parquet" for i in range(50 - len(LEGACY_1_0)))
    _extras(out, set(LEGACY_DIRS), (set(LEGACY_1_0) - LEGACY_DIRS) | set(extra))
    assert len(LEGACY_1_0) + len(extra) == 50
    text = _render(out)
    assert text.isascii()
    assert len(text.splitlines()) <= 60, f"{len(text.splitlines())} lines:\n{text}"
    missing = [n for n in LEGACY_1_0 if not _names(text, n, n in LEGACY_DIRS)]
    missing += [n for n in extra if not _names(text, n, False)]
    assert not missing, f"not named: {missing}\n{text}"


@pytest.mark.parametrize("name", ["genetable", "codon.parquet", "context", "novel"])
def test_a_known_1_0_name_is_legacy_and_removable(tmp_path: Path, name: str) -> None:
    from gippa.io.out_readme import legacy_entries

    out = _full(tmp_path / "OUT")
    is_dir = name in LEGACY_DIRS
    _extras(out, {name} if is_dir else set(), set() if is_dir else {name})
    assert legacy_entries(out) == [name]
    text = _render(out)
    assert _names(text, name, is_dir), text
    assert any(REMOVABLE in p for p in _paragraphs_naming(text, name, is_dir)), text


@pytest.mark.parametrize("name", [*FOREIGN, *NON_ASCII, "my_results"])
def test_a_users_own_entry_is_named_not_legacy_not_removable(tmp_path: Path, name: str) -> None:
    """M2: only 1.0 step intermediates are legacy; a user's file next to the results is named
    but never called removable (and triggers no ``legacy_intermediates`` warning)."""
    from gippa.io.out_readme import legacy_entries, write_readme

    out = _full(tmp_path / "OUT")
    is_dir = name == "my_results"
    _extras(out, {name} if is_dir else set(), set() if is_dir else {name})
    assert legacy_entries(out) == []
    path = write_readme(out)  # must not raise on a non-ASCII name
    text = path.read_bytes().decode("ascii")
    assert text == _render(out)
    assert _names(text, name, is_dir), f"{_ascii(name)!r} not named:\n{text}"
    assert "removed" not in text, text


def test_users_entry_and_legacy_entry_together(tmp_path: Path) -> None:
    out = _full(tmp_path / "OUT")
    _extras(out, {"genetable"}, {"samplesheet.csv"})
    text = _render(out)
    assert any(REMOVABLE in p for p in _paragraphs_naming(text, "genetable", True)), text
    assert _names(text, "samplesheet.csv", False), text
    assert all("removed" not in p for p in _paragraphs_naming(text, "samplesheet.csv", False)), text


def test_failed_run_without_report_names_no_report(tmp_path: Path) -> None:
    out = _tree(tmp_path / "OUT", ("tables", "preflight", "logs", "work"), FILES)
    text = _render(out)
    assert "report" not in _named_heads(text)
    assert "report.html" not in text
    assert "  tables/genes.tsv" in text


@pytest.mark.parametrize("parent", ["home/u/farmshare/GIPPA_ROOT/local", "plain"])
def test_plain_ascii_no_timestamp_no_path(tmp_path: Path, parent: str) -> None:
    out = _full(tmp_path / parent / "OUT")
    text = _render(out)
    assert text.isascii()
    assert len(text.splitlines()) <= 60
    assert str(out) not in text and str(tmp_path) not in text
    assert "/home" not in text and "local/" not in text
    assert "farmshare" not in text.lower() and "$GIPPA_ROOT" not in text
    assert not re.search(r"\d{4}-\d{2}-\d{2}|\d{1,2}:\d{2}", text), "timestamp in README"


def test_no_internal_references(tmp_path: Path) -> None:
    """§ 1 rule 2: no ADR, § or spec references and no flag tokens in user-facing prose."""
    text = _render(_full(tmp_path / "OUT"))
    assert not re.search(r"ADR-\d{4}|§\d|GIPPA_SPEC\.md|docs/methods", text)


def test_render_is_deterministic(tmp_path: Path) -> None:
    out = _full(tmp_path / "OUT")
    first = _render(out)
    (out / "README.txt").write_text(first)
    assert _render(out) == first


# ---------------------------------------------------------------- in sync with the publish map


def test_every_published_top_level_name_has_a_sentence(tmp_path: Path) -> None:
    """Every top-level name the publish map can produce (``PATHS``, ``PUBLISHED``) has its own
    README line, so a new published folder cannot be added without one."""
    heads = {rel.split("/")[0] for rel in PATHS.values()} | set(PUBLISHED.values())
    folders = {h for h in heads if "." not in h} | {"schemas", "intermediate", "work"}
    files = {h for h in heads if "." in h}
    out = _tree(tmp_path / "OUT", tuple(sorted(folders)), tuple(sorted(files)))
    missing = (folders | files) - _named_heads(_render(out))
    assert not missing, f"README has no line for {sorted(missing)}"


# ---------------------------------------------------------------- docs (D-01, outputs.md)


def test_outputs_guide_shows_intermediate_and_readme() -> None:
    text = (ROOT / "docs" / "user-guide" / "outputs.md").read_text(encoding="utf-8")
    assert "intermediate/" in text
    assert "README.txt" in text
    assert "Some files sit outside this list" not in text
