"""``OUT/README.txt``: a plain-text guide to the output folder, written at the end of every
``gippa run`` and ``gippa demo`` (also after a failed run).

``render_readme`` names only the top-level entries present in ``OUT/``, in a fixed order, each
with a fixed sentence; it ends with the coordinate note, the table reference link and the GIPPA
version. The text is ASCII, with no timestamp and no path outside ``OUT/``, so two runs that leave
the same entries get the same README. Every other non-hidden top-level entry is named too
(non-ASCII names backslash-escaped): the step intermediates GIPPA 1.0 left at the top of
``OUT/`` (``LEGACY_TOP``, now under ``intermediate/``) in a paragraph saying they can be
removed, anything else (a user's own files) in a separate neutral paragraph.
``legacy_entries`` lists the ``LEGACY_TOP`` names present.
"""

from __future__ import annotations

import textwrap
from collections.abc import Iterable, Sequence
from pathlib import Path

import gippa
from gippa.runner.local import INTERMEDIATE

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
#: top-level step intermediates GIPPA 1.0 published at the top of ``OUT/`` (1.1 writes them under
#: ``intermediate/``): the first path parts of the 1.0 workflow's published outputs outside the
#: published folders and ``work/``
LEGACY_TOP: frozenset[str] = frozenset({
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
})  # fmt: skip
#: width of the wrapped name lists (two-space indent included)
_WRAP = 96


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
    """The ``LEGACY_TOP`` names present at the top of ``out``, sorted: the step intermediates
    GIPPA 1.0 published there (a user's own files are never legacy)."""
    return sorted(_entries(out) & LEGACY_TOP)


def _token(out: Path, name: str) -> str:
    """How the README writes a top-level entry: ``name/`` for a folder, the file name for a file,
    backslash-escaped outside ASCII."""
    text = name.encode("ascii", "backslashreplace").decode("ascii")
    return f"{text}/" if (out / name).is_dir() else text


def _name_list(title: str, out: Path, names: Iterable[str]) -> list[str]:
    """``title`` then the entries ``names``, comma-separated and wrapped, then a blank line."""
    tokens = ", ".join(_token(out, n) for n in names)
    wrapped = textwrap.wrap(tokens, _WRAP, initial_indent="  ", subsequent_indent="  ",
                            break_on_hyphens=False, break_long_words=False)  # fmt: skip
    return [title, *wrapped, ""]


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
    legacy = legacy_entries(out)
    if legacy:
        lines += _name_list(
            "From GIPPA 1.0 (now under intermediate/); they can be removed:", out, legacy
        )
    shown = {n for n, _ in FOLDERS if n in present and (out / n).is_dir()}
    shown |= {n for n, _ in FILES if n in present and (out / n).is_file()}
    other = sorted(present - shown - set(legacy))
    if other:
        lines += _name_list("Also here, not written by GIPPA:", out, other)
    lines += [COORDINATES, f"Table reference: {DOCS_TABLES}",
              f"Written by GIPPA {gippa.__version__}."]  # fmt: skip
    return "\n".join(lines) + "\n"


def write_readme(out: Path) -> Path:
    """Write ``out / "README.txt"`` (ASCII) and return its path."""
    path = out / README
    out.mkdir(parents=True, exist_ok=True)
    path.write_bytes(render_readme(out).encode("ascii"))
    return path
