"""``OUT/README.txt``: a plain-text guide to the output folder, written at the end of every
``gippa run`` and ``gippa demo`` (also after a failed run).

``render_readme`` names only the top-level entries present in ``OUT/``, in a fixed order, each
with a fixed sentence; it ends with the coordinate note, the table reference link and the GIPPA
version. The text is ASCII, with no timestamp and no path outside ``OUT/``, so two runs that leave
the same entries get the same README. ``legacy_entries`` lists the top-level step intermediates
an earlier GIPPA version left at the top of ``OUT/`` (they now live under ``intermediate/``).
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import gippa
from gippa.runner.local import INTERMEDIATE, PUBLISHED_TOP

README = "README.txt"
DOCS_TABLES = "https://khoa-yelo.github.io/gippa/reference/tables/"

#: (path below OUT/, sentence): the files to open first, when present
START_HERE: tuple[tuple[str, str], ...] = (
    ("report/report.html", "the report: open it in a web browser"),
    ("tables/genes.tsv", "one row per gene: detection, abundance, function, context, features"),
)
#: (folder, sentence), in README order
FOLDERS: tuple[tuple[str, str], ...] = (
    ("tables", "result tables (TSV); every column is described in schemas/datapackage.json"),
    ("mapping", "where each peptide or protein was placed in the genome, and what was not placed"),
    ("annotation", "the gene model and annotation layers this run used; tools/ holds native tool "
                   "outputs"),
    ("plots", "figures (SVG, PNG, PDF); plots/data/ holds the numbers behind each figure"),
    ("report", "report.html, self-contained"),
    ("preflight", "the input checks made before the analysis"),
    ("schemas", "a machine-readable description of every table (Frictionless Data Package)"),
    ("logs", "step logs, events and timings"),
    (INTERMEDIATE, "files passed between steps; not results; their layout may change between "
                   "versions"),
    ("work", "the workflow engine's working files; needed only to resume this run"),
)  # fmt: skip
#: (file, sentence), in README order
FILES: tuple[tuple[str, str], ...] = (
    ("run_manifest.json", "versions, databases, parameters, input checksums and seeds"),
    ("resolved_config.yaml", "the configuration this run used"),
)
COORDINATES = (
    "Coordinates: TSV, GFF3 and VCF are 1-based closed; BED and proBed are 0-based half-open."
)
#: name column widths: "Start here" and "Files" share one, "Folders" has its own
_WIDE = max(len(n) for n, _ in (*START_HERE, *FILES))
_FOLDER_WIDE = max(len(n) + 1 for n, _ in FOLDERS)
#: top-level names GIPPA 1.1 writes or Nextflow keeps (everything else is legacy or foreign)
_KNOWN = PUBLISHED_TOP | {INTERMEDIATE, "work"}


def _entries(out: Path) -> set[str]:
    """Non-hidden top-level names of ``out``, without ``README.txt``."""
    if not out.is_dir():
        return set()
    return {p.name for p in out.iterdir() if not p.name.startswith(".") and p.name != README}


def _section(title: str, rows: Sequence[tuple[str, str]], width: int) -> list[str]:
    if not rows:
        return []
    return [f"{title}:", *(f"  {name.ljust(width)}  {text}" for name, text in rows), ""]


def legacy_entries(out: Path) -> list[str]:
    """Top-level entries of ``out`` that are not published names, ``intermediate/`` or ``work/``
    (hidden entries ignored), sorted: the step intermediates an earlier version published there."""
    return sorted(_entries(out) - _KNOWN)


def count_files(out: Path, names: Sequence[str]) -> int:
    """Files in the top-level entries ``names`` of ``out`` (a file counts 1, a folder its files)."""
    total = 0
    for name in names:
        path = out / name
        total += sum(1 for p in path.rglob("*") if p.is_file()) if path.is_dir() else 1
    return total


def render_readme(out: Path) -> str:
    """The README text for the entries present in ``out`` (hidden entries and ``README.txt``
    itself are ignored)."""
    present = _entries(out)
    start = [(rel, text) for rel, text in START_HERE if (out / rel).is_file()]
    folders = [(f"{n}/", t) for n, t in FOLDERS if n in present and (out / n).is_dir()]
    files = [(n, t) for n, t in FILES if n in present and (out / n).is_file()]
    lines = ["GIPPA results", ""]
    lines += _section("Start here", start, _WIDE)
    lines += _section("Folders", folders, _FOLDER_WIDE)
    lines += _section("Files", files, _WIDE)
    legacy = len(legacy_entries(out))
    if legacy:
        lines += [f"Also here: {legacy} other top-level entries, not written by this version "
                  "(an earlier version's step intermediates can be removed).", ""]  # fmt: skip
    lines += [COORDINATES, f"Table reference: {DOCS_TABLES}",
              f"Written by GIPPA {gippa.__version__}."]  # fmt: skip
    return "\n".join(lines) + "\n"


def write_readme(out: Path) -> Path:
    """Write ``out / "README.txt"`` (ASCII) and return its path."""
    path = out / README
    out.mkdir(parents=True, exist_ok=True)
    path.write_bytes(render_readme(out).encode("ascii"))
    return path
