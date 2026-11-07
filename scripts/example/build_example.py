"""Build the quick-start example bundle: one Nissle 1917 chromosome region, five input routes.

The example is a cut-down subset of the M10 pilot (validation/pilot/README.md, pass 3): E. coli
Nissle 1917 (GCF_021559835.1, chromosome NZ_CP058217.1), PXD078459 runs WT_1-WT_3. Nothing is
re-searched. The script subsets the pilot's existing outputs to one contiguous region and writes
the release tree OUT/<name>/gippa-example-<X.Y>/ (spec M12-01 section 2 > Layout):

    genome/<name>.fna.gz          the shared region (routes A-D) as one linear contig (no tag)
    genome/<name>_e.fna.gz        route E's region, written by the same rules
    A_diann/report.parquet        DIA-NN 2.3.1 report rows whose proteins all lie in the region
    B_proteins/proteins.faa       the region proteins DIA-NN quantified (A2 strain database)
    B_proteins/abundance.tsv      their PG.MaxLFQ (median over runs), as `gippa convert` computes it
    C_cascadia/WT_{1,2,3}.ssl     Cascadia rows whose peptide occurs in the region's six frames
    D_generic/peptides.tsv        a generic peptide table, peptide level, built from route A's rows
    E_sixframe/report.parquet     DIA-NN rows of the six-frame search A4 in the region (rule 3E, 9)
    SELECTION.json                the rule, the counts and the sha256 of every source and output
    README.md, PROVENANCE.md, LICENSE-DATA.md, bundle.json, CHECKSUMS.sha256, expected/ (optional)

then packs it as BUNDLE_OUT/gippa-example-<X.Y>.tar.gz (deterministic: sorted members, mtime,
uid and gid 0, empty owner names, modes 0644/0755, gzip without name and mtime) and prints the
archive's sha256 and size, its db_sources.yaml row and its validation/example/bundles.tsv row.
Sizing tars per route stay in OUT/<name>/bundles/, outside the release tree.

The tree is built in a staging dir next to it and replaces the old tree only after the archive is
released; a failed or refused build leaves the old tree, the old archive and no partial file.

Selection rule (recorded verbatim in SELECTION.json):

* Region: NZ_CP058217.1 [start, end) in 0-based half-open coordinates, shared by routes A-D.
* Route E's region: NZ_CP058217.1 [start_e, end_e), written by the same rules to its own genome
  file; the route E rules below use it in place of the shared region.
* Region proteins: records of the pass-3 A2 strain database (`db_id_map.tsv`, 1-based closed)
  on NZ_CP058217.1 with start - 1 >= region start and end <= region end (CDS and sORF).
* Route A: a DIA-NN row is kept iff every ID in `Protein.Ids` (split on ';') is a region protein.
  This applies to target and decoy rows alike (DIA-NN decoys carry the target IDs with Decoy = 1),
  so the decoys keep their `protein_ids`. Rows with a contaminant (`Cont_`) or out-of-region ID are
  dropped. All 71 columns are kept. q-values, PEPs and MaxLFQ are DIA-NN's, computed on the whole
  proteome search; they are not recomputed for the subset.
* Route B: the protein groups of route A's target rows with a non-missing, non-zero PG.MaxLFQ;
  value = median over runs of PG.MaxLFQ (the `gippa convert --abundance-out` rule). A group "A;B" is
  written as one row per member ID; a protein in several groups is written once, with the value of
  the first of its groups in sorted order. The FASTA holds each written protein's sequence from the
  strain database (a protein missing from it fails the build).
* Route C: a Cascadia row is kept iff its stripped sequence (I and L equated) occurs as an exact
  substring of a stop-free stretch of the region's six-frame translation (genetic code 11). The
  `file` column is reduced to the mzML base name (the run name is its stem).
* Route D: route A's rows grouped by (stripped peptide, run): q_value = min Global.Q.Value,
  pep = min PEP, intensity = sum of Precursor.Normalised (targets), protein_ids = the union,
  is_decoy = all rows decoy; modified_peptide = the first Modified.Sequence in ProForma
  (`(UniMod:n)` -> `[UNIMOD:n]`). Columns: peptide, modified_peptide, protein_ids, q_value, pep,
  is_decoy, intensity, run. Numbers keep full precision in every table.

* Route E: the DIA-NN rows of the six-frame search A4 (`gippa makedb --six-frame-orfs`
  database). Region entries are the A4 `db_id_map.tsv` records on NZ_CP058217.1 with
  start - 1 >= region start and end <= region end, of kind cds, sorf or orf, on either strand (an
  ORF is placed by its map coordinates). A row is kept iff every ID in `Protein.Ids` (split on
  ';') is a region entry, targets and decoys alike; rows with a contaminant (`Cont_`) ID, an
  out-of-region protein or an out-of-region ORF (including one crossing a region end) are
  dropped. All columns are kept.
* Route E ORF IDs in region coordinates: every ORF ID
  `GIPPAORF_{class}_NZ_CP058217.1_{s1}_{e1}_{p|m}` of a kept row, in any text column, becomes
  `GIPPAORF_{class}_{C}_{s1 - S}_{e1 - S}_{p|m}` (C = route E's genome contig, S = route E's 0-based
  region start). Class, strand and the ID order within a cell are unchanged; CDS and sORF IDs are
  unchanged. The build fails if an ID's coordinates differ from its `db_id_map.tsv` row, if the
  rewritten ID does not parse back (gippa.novel.fdr.parse_orf_id) to its region interval, if the
  region sequence there (code 11, on its strand) differs from the ID's `strain_db.faa` record beyond
  the first residue (not checked for nterm_extension), or if an ORF ID naming NZ_CP058217.1 is left.
* Scrubbing (R12): the `.ssl` `file` column and every DIA-NN text column holding an absolute
  path are reduced to base names; SELECTION.json lists the rewritten columns. No other value is
  changed, except by the route E ORF-ID rewrite. The build fails if any bundle file (a cell,
  parquet metadata or an `expected/` file) still holds a personal or build path.
* No empty route: the build fails if any route keeps zero rows or zero proteins, if route A
  keeps no decoy, or if route E keeps no target row, no decoy row or no row with a region ORF ID.

Usage (venv python, PYTHONPATH=<checkout>/src; small, runs in seconds). Defaults: shared region
NZ_CP058217.1:1965001-2165000 (--start/--end move it), route E's region
NZ_CP058217.1:1400001-1600000 (START_E/END_E, not an option), --name nissle_r200k, --out and
--bundle-out $GIPPA_ROOT/example, bundle version = the installed GIPPA MAJOR.MINOR. A checkout
with uncommitted changes fails unless --allow-dirty (the commit is then recorded as
<sha>-dirty). --check rebuilds in a temporary dir and compares the archive with its
validation/example/bundles.tsv row; it writes nothing to --out or --bundle-out:

    python scripts/example/build_example.py [--expected DIR] [--bundle-out DIR] [--allow-dirty]
    python scripts/example/build_example.py --check [--bundles-tsv TSV] [--expected DIR]
"""

from __future__ import annotations

import argparse
import contextlib
import csv
import getpass
import gzip
import hashlib
import importlib.util
import io
import json
import os
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
from collections.abc import Iterable, Mapping
from importlib.resources import files as package_files
from pathlib import Path
from types import ModuleType
from typing import Any

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

CHROM = "NZ_CP058217.1"
RUNS = ("WT_1", "WT_2", "WT_3")
START, END, NAME = 1_965_000, 2_165_000, "nissle_r200k"
# Route E's own region (rule 1E); read at build time, so a test may set them on the module
START_E, END_E = 1_400_000, 1_600_000
HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
BUNDLES_TSV = REPO / "validation" / "example" / "bundles.tsv"
_ABS = re.compile(r"^(/|[A-Za-z]:[\\/]|\\\\)")  # POSIX, Windows drive or UNC
# a personal or cluster path anywhere in a value (R12); URLs do not match (a word char precedes)
_PERSONAL = re.compile(
    r"(?<![\w.~-])/(?:home|Users|scratch|farmshare|oak|gpfs|mnt|media|private|tmp|var/folders)/"
    r"|[A-Za-z]:\\Users\\"
)
_GIT_ENV = ("GIT_DIR", "GIT_WORK_TREE", "GIT_INDEX_FILE")

# Static digests for PROVENANCE.md (M3): (label, path under $GIPPA_ROOT or the checkout, sha256).
# The file's own sha256 is used where it exists; the pinned value otherwise.
STATIC_DIGESTS: dict[str, tuple[str, str, str]] = {
    "diann_image": (
        "diann_2.3.1.sif",
        "",  # not under $GIPPA_ROOT/pilot; validation/pilot/diann_image.sha256
        "1946452d56d2fb0f66e05ced473df88e53484e28040f1ecbd83c4a43405a6381",
    ),
    "diann_cfg": (
        "validation/pilot/diann_A2.cfg",
        "repo:validation/pilot/diann_A2.cfg",
        "c0d145e87e582e450295517ac393d2c86093587605dd3da66c6f6a923c1f98e4",
    ),
    "a2_abundance": (
        "pilot/pass3/inputs/a2_abundance.tsv.gz",
        "root:pilot/pass3/inputs/a2_abundance.tsv.gz",
        "3970344b9c0edc0fcf756b267063a1d1c30c05f6b6bdb6654d973eae1e6f888e",
    ),
    "cascadia_weights": (
        "cascadia_0.0.7.ckpt",
        "",  # $GIPPA_DB_DIR/cascadia; validation/pilot/a3_assets.tsv
        "0f5e91413c716b22019699c8cb25655216f9c96d0f2aeb6ec42084bb65652192",
    ),
}


class BuildError(Exception):
    """The inputs cannot give a valid bundle (rule 7, a leak, an existing different archive)."""


def _load_text() -> ModuleType:
    """``bundle_text.py`` next to this script (loaded by path: scripts/ is not a package)."""
    spec = importlib.util.spec_from_file_location(
        "gippa_example_bundle_text", HERE / "bundle_text.py"
    )
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


T = _load_text()


def _root() -> Path:
    return Path(os.environ.get("GIPPA_ROOT", "/home/classes/bios/270/khoa/giap/local"))


A4_DIR = Path("runs") / "pilot" / "sixframe-fa49b8f"  # the pilot's six-frame arm (route E)
# Route E's PROVENANCE extras under $GIPPA_ROOT/A4_DIR: sha256 recorded when present (not sources)
A4_EXTRAS = ("jobs/diann_A4.cfg", "db/makedb_manifest.json", "A4/diann_version.txt")
A4_COMMANDS = ("A4/command_step1.txt", "A4/command_step2.txt")


def sources(root: Path) -> dict[str, Path]:
    p3 = root / "pilot" / "pass3"
    a3 = root / "pilot" / "A3_cascadia_cpu"
    a4 = root / A4_DIR
    return {
        "genome": p3 / "genome" / "nissle1917_wt.fna.gz",
        "db_id_map": p3 / "db_a2" / "db_id_map.tsv",
        "strain_db": p3 / "db_a2" / "strain_db.faa.gz",
        "diann_report": p3 / "A2" / "A2_report.parquet",
        **{f"ssl_{r}": a3 / r / f"{r}.ssl" for r in RUNS},
        "a4_report": a4 / "A4" / "A4_report.parquet",
        "a4_db_id_map": a4 / "db" / "db_id_map.tsv",
        "a4_strain_db": a4 / "db" / "strain_db.faa",
    }


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_fasta(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    name: str | None = None
    seq: list[str] = []
    with (
        gzip.open(path, "rt", encoding="utf-8")
        if path.suffix == ".gz"
        else path.open(encoding="utf-8")
    ) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    out[name] = "".join(seq)
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line.strip())
    if name is not None:
        out[name] = "".join(seq)
    return out


def write_fasta(path: Path, records: Iterable[tuple[str, str]], width: int = 60) -> None:
    buf = io.StringIO()
    for name, seq in records:
        buf.write(f">{name}\n")
        for i in range(0, len(seq), width):
            buf.write(seq[i : i + width] + "\n")
    data = buf.getvalue().encode("utf-8")
    if path.suffix == ".gz":
        # mtime=0 and no file name: byte-identical rebuilds
        with (
            path.open("wb") as raw,
            gzip.GzipFile(fileobj=raw, mode="wb", mtime=0, filename="") as gz,
        ):
            gz.write(data)
    else:
        path.write_bytes(data)


_CODON = {}
_BASES = "TCAG"
_AA11 = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for _i, _a in enumerate(_BASES):
    for _j, _b in enumerate(_BASES):
        for _k, _c in enumerate(_BASES):
            _CODON[_a + _b + _c] = _AA11[16 * _i + 4 * _j + _k]


def translate(nt: str) -> str:
    return "".join(_CODON.get(nt[i : i + 3], "X") for i in range(0, len(nt) - 2, 3))


def revcomp(nt: str) -> str:
    return nt.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def six_frames(nt: str) -> list[str]:
    rc = revcomp(nt)
    return [translate(s[f:]) for s in (nt, rc) for f in range(3)]


def strip_mods(seq: str) -> str:
    seq = re.sub(r"\[[^\]]*\]|\([^)]*\)", "", seq)
    return re.sub(r"[^A-Z]", "", seq)


def region_proteins(db_map: Path, start: int, end: int) -> pd.DataFrame:
    m = pd.read_csv(db_map, sep="\t", dtype={"db_id": str}, encoding="utf-8")
    inside = (m["contig"] == CHROM) & (m["start"] - 1 >= start) & (m["end"] <= end)
    return m[inside].reset_index(drop=True)


def subset_diann(report: Path, keep_ids: set[str]) -> tuple[pa.Table, dict[str, int]]:
    t = pq.read_table(report)
    for i, f in enumerate(t.schema):  # Arrow cannot filter string_view columns: read as string
        if pa.types.is_string_view(f.type):
            t = t.set_column(i, pa.field(f.name, pa.string()), t.column(i).cast(pa.string()))
    ids = t.column("Protein.Ids").to_pylist()
    mask = [bool(x) and all(p in keep_ids for p in str(x).split(";")) for x in ids]
    sub = t.filter(pa.array(mask))
    decoy = sub.column("Decoy").to_pylist()
    counts = {
        "rows_total": t.num_rows,
        "rows_kept": sub.num_rows,
        "decoy_rows_kept": int(sum(1 for d in decoy if d)),
        "target_rows_kept": int(sum(1 for d in decoy if not d)),
        "rows_mixed_dropped": int(
            sum(
                1
                for x in ids
                if x
                and any(p in keep_ids for p in str(x).split(";"))
                and not all(p in keep_ids for p in str(x).split(";"))
            )
        ),
    }
    return sub, counts


_ORF_CHROM = re.compile(rf"GIPPAORF_(\w+?)_{re.escape(CHROM)}_(\d+)_(\d+)_([pm])")
_ORF_ID = re.compile(rf"^{_ORF_CHROM.pattern}$")


def region_entries(db_map: Path, start: int, end: int) -> pd.DataFrame:
    """Rule 3E: the A4 map records (cds, sorf, orf) contained in the region, either strand."""
    m = region_proteins(db_map, start, end)
    return m[m["kind"].isin(["cds", "sorf", "orf"])].reset_index(drop=True)


def _is_text(typ: pa.DataType) -> bool:
    """A string column: plain, large, view or dictionary-encoded."""
    if pa.types.is_dictionary(typ):
        return _is_text(typ.value_type)
    return bool(pa.types.is_string(typ) or pa.types.is_large_string(typ)
                or pa.types.is_string_view(typ))  # fmt: skip


def _text_columns(t: pa.Table) -> list[int]:
    return [i for i, f in enumerate(t.schema) if _is_text(f.type)]


def _set_text(t: pa.Table, i: int, values: list[Any]) -> pa.Table:
    """Column ``i`` of ``t`` replaced by ``values`` in the column's own type (BuildError when the
    type cannot hold them)."""
    field = t.schema.field(i)
    try:
        if pa.types.is_dictionary(field.type):
            arr = pa.array(values, type=field.type.value_type).dictionary_encode()
            arr = arr.cast(field.type)
        else:
            arr = pa.array(values, type=field.type)
    except (pa.ArrowInvalid, pa.ArrowNotImplementedError, pa.ArrowTypeError) as exc:
        raise BuildError(f"column {field.name} ({field.type}) cannot be rewritten: {exc}") from exc
    return t.set_column(i, field, arr)


def orf_id_rewrite(
    t: pa.Table, orfs: pd.DataFrame, region_nt: str, start: int, contig: str,
    records: Mapping[str, str],
) -> tuple[pa.Table, list[str], int]:  # fmt: skip
    """Rule 9: the region ORF IDs of ``t`` (``orfs``: their map rows) in region coordinates;
    (table, rewritten columns, distinct IDs rewritten). A failed check raises BuildError."""
    from gippa.mapping.sixframe import translate_text
    from gippa.novel.fdr import parse_orf_id

    by_id = {str(r.db_id): r for r in orfs.itertuples(index=False)}
    new_id: dict[str, str] = {}

    def renamed(old: str) -> str:
        if old in new_id:
            return new_id[old]
        m = _ORF_ID.match(old)
        row = by_id[old]
        if m is None:
            raise BuildError(f"route E: region ORF ID {old} is not "
                             f"GIPPAORF_{{class}}_{CHROM}_{{s1}}_{{e1}}_{{p|m}}")  # fmt: skip
        cls, s1, e1, sfx = m.group(1), int(m.group(2)), int(m.group(3)), m.group(4)
        strand = "+" if sfx == "p" else "-"
        new = f"GIPPAORF_{cls}_{contig}_{s1 - start}_{e1 - start}_{sfx}"
        if (int(row.start), int(row.end), str(row.strand)) != (s1, e1, strand):
            raise BuildError(f"route E: ORF ID {old} (rewritten {new}) disagrees with its "
                             f"db_id_map.tsv row {row.start}-{row.end} {row.strand}")  # fmt: skip
        want = (cls, contig, s1 - start - 1, e1 - start, strand)
        if parse_orf_id(new) != want:
            raise BuildError(f"route E: rewritten ORF ID {new} (from {old}) does not parse back "
                             f"to {want}")  # fmt: skip
        if cls != "nterm_extension":
            nt = region_nt[s1 - start - 1 : e1 - start]
            aa = translate_text(revcomp(nt) if strand == "-" else nt, 11)
            rec = records.get(old)
            if rec is None or len(rec) != len(aa) or rec[1:] != aa[1:]:
                raise BuildError(f"route E: the region translation of {new} (from {old}) differs "
                                 "from its strain_db.faa record")  # fmt: skip
        new_id[old] = new
        return new

    rewritten = []
    for i in _text_columns(t):
        vals = t.column(i).to_pylist()
        out = [";".join(renamed(p) if p in by_id else p for p in v.split(";"))
               if isinstance(v, str) else v for v in vals]  # fmt: skip
        if out != vals:
            t = _set_text(t, i, out)
            rewritten.append(t.schema.field(i).name)
    left = sorted({m.group(0) for i in range(t.num_columns) for v in t.column(i).to_pylist()
                   if isinstance(v, str) for m in _ORF_CHROM.finditer(v)})  # fmt: skip
    if left:
        raise BuildError(f"route E: ORF IDs naming {CHROM} left after the rewrite (not "
                         f"region entries, or in a column not rewritten): {left[:5]}")  # fmt: skip
    return t, rewritten, len(new_id)


# a path argument: absolute, ~/ or ../; alone, after `=` or `:`, or glued to a flag (`-B/x:/y`)
_PATH_ARG = re.compile(
    r"(?:^|(?<=[\s=:'\"]))(?P<flag>-{1,2}[A-Za-z][\w-]*)?"
    r"(?P<path>(?:~|\.\.)?/(?!/)[^\s:'\"=,;]*)",
    re.M,
)
_EXTENSION = re.compile(r"[^.]\.[A-Za-z0-9]+$")


def shorten_paths(text: str) -> str:
    """R12: each path argument of ``text`` becomes its base name when its last component has a
    file extension, and ``<dir>`` otherwise (a trailing ``/`` or no extension)."""

    def short(m: re.Match[str]) -> str:
        path = m.group("path")
        last = path.rsplit("/", 1)[-1]
        name = last if _EXTENSION.search(last) else "<dir>"
        return (m.group("flag") or "") + name

    return _PATH_ARG.sub(short, text)


def login_names() -> list[str]:
    """The building user's login names (``getpass.getuser()``, the home directory's name)."""
    names = []
    with contextlib.suppress(Exception):  # getuser raises when no name is set
        names.append(getpass.getuser())
    names.append(Path.home().name)
    return sorted({n for n in names if n and n not in (".", "/")})


def a4_extras(root: Path) -> dict[str, Any]:
    """PROVENANCE's route E extras (present files only): sha256s and the shortened commands."""
    a4 = root / A4_DIR
    digests = {rel: sha256(a4 / rel) for rel in A4_EXTRAS if (a4 / rel).is_file()}
    commands = {Path(rel).name: shorten_paths((a4 / rel).read_text(encoding="utf-8")).strip()
                for rel in A4_COMMANDS if (a4 / rel).is_file()}  # fmt: skip
    for rel in (*A4_EXTRAS, *A4_COMMANDS):
        if not (a4 / rel).is_file():
            print(f"build_example: warning: $GIPPA_ROOT/{(A4_DIR / rel).as_posix()} not found; "
                  "PROVENANCE.md omits it", file=sys.stderr)  # fmt: skip
    return {"digests": digests, "commands": commands}


def abundance(sub: pd.DataFrame) -> tuple[pd.DataFrame, int]:
    """Route B: one row per protein (first group in sorted order wins); (table, number of
    proteins that sit in more than one group)."""
    t = sub[(sub["Decoy"] == 0)]
    t = t[t["PG.MaxLFQ"].notna() & (t["PG.MaxLFQ"] > 0)]
    per_run = t.groupby(["Protein.Group", "Run"])["PG.MaxLFQ"].first().reset_index()
    med = per_run.groupby("Protein.Group")["PG.MaxLFQ"].median()
    value: dict[str, float] = {}
    groups: dict[str, int] = {}
    for group, v in sorted(med.items()):
        for pid in str(group).split(";"):
            value.setdefault(pid, float(v))
            groups[pid] = groups.get(pid, 0) + 1
    frame = pd.DataFrame(list(value.items()), columns=["protein_id", "maxlfq"])
    return frame, sum(1 for n in groups.values() if n > 1)


def generic_table(sub: pd.DataFrame) -> pd.DataFrame:
    d = sub.copy()
    d["peptide"] = d["Stripped.Sequence"].astype(str)
    d["proforma"] = (
        d["Modified.Sequence"]
        .astype(str)
        .str.replace(r"\(UniMod:(\d+)\)", r"[UNIMOD:\1]", regex=True)
    )
    d["is_decoy"] = d["Decoy"].astype(int) == 1
    d["inten"] = d["Precursor.Normalised"].where(~d["is_decoy"])
    d["inten"] = d["inten"].where(d["inten"] > 0)
    rows = []
    for (pep, run), g in d.groupby(["peptide", "Run"], sort=True):
        pids = sorted({p for x in g["Protein.Ids"].astype(str) for p in x.split(";")})
        inten = g["inten"].dropna()
        rows.append(
            {
                "peptide": pep,
                "modified_peptide": g["proforma"].iloc[0],
                "protein_ids": ";".join(pids),
                "q_value": float(g["Global.Q.Value"].min()),
                "pep": float(g["PEP"].min()),
                "is_decoy": "true" if bool(g["is_decoy"].all()) else "false",
                "intensity": float(inten.sum()) if len(inten) else "",
                "run": run,
            }
        )
    return pd.DataFrame(rows)


def subset_ssl(path: Path, frames_il: list[str]) -> tuple[list[dict[str, str]], list[str], int]:
    with path.open(encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        header = list(reader.fieldnames or [])
        rows = list(reader)
    kept = []
    for r in rows:
        pep = strip_mods(r["sequence"]).replace("I", "L")
        if len(pep) >= 1 and any(pep in f for f in frames_il):
            r = dict(r)
            r["file"] = Path(r["file"]).name
            kept.append(r)
    return kept, header, len(rows)


def write_tsv_rows(path: Path, header: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=header, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def write_table(frame: pd.DataFrame, path: Path) -> None:
    """A TSV with full float precision (``repr``: every value round-trips)."""
    frame.to_csv(path, sep="\t", index=False, lineterminator="\n", encoding="utf-8")


def scrub_paths(t: pa.Table) -> tuple[pa.Table, list[str]]:
    """Rule 8: every text column holding an absolute path, with each such value reduced to its
    base name (``/`` or backslash separated); the names of the rewritten columns."""
    rewritten = []
    for i in _text_columns(t):
        vals = t.column(i).to_pylist()
        if not any(isinstance(v, str) and _ABS.match(v) for v in vals):
            continue
        new = [re.split(r"[\\/]", v)[-1] if isinstance(v, str) and _ABS.match(v) else v
               for v in vals]  # fmt: skip
        t = _set_text(t, i, new)
        rewritten.append(t.schema.field(i).name)
    return t, rewritten


def _leaks(text: str, literals: Iterable[str]) -> bool:
    return bool(_PERSONAL.search(text)) or any(lit and lit in text for lit in literals)


def personal_paths(
    top: Path, literals: Iterable[str | Path], words: Iterable[str] = ()
) -> list[str]:
    """R12 leak guard: ``<file> (<where>)`` for every bundle file holding a personal or build
    path (a parquet cell or metadata entry, a gzipped or plain text file), or, in a text file,
    one of ``words`` (the login names) as a whole word. Values are not reported."""
    lits = [str(x) for x in literals if str(x) not in ("", "/", ".")]
    word = re.compile("|".join(rf"(?<!\w){re.escape(w)}(?!\w)" for w in words if w) or r"(?!)")
    hits: list[str] = []
    for p in sorted(top.rglob("*")):
        if not p.is_file():
            continue
        rel = p.relative_to(top).as_posix()
        if p.suffix == ".parquet":
            t = pq.read_table(p)
            for name in t.column_names:
                vals = t.column(name).to_pylist()
                if any(isinstance(v, str) and _leaks(v, lits) for v in vals):
                    hits.append(f"{rel} (column {name})")
            meta = pq.ParquetFile(p).metadata.metadata or {}
            for k, v in sorted(meta.items()):
                if k == b"ARROW:schema":
                    continue
                if any(_leaks(x.decode("utf-8", "replace"), lits) for x in (k, v)):
                    hits.append(f"{rel} (parquet metadata {k.decode('utf-8', 'replace')})")
            continue
        raw = p.read_bytes()
        if p.suffix == ".gz":
            raw = gzip.decompress(raw)
        text = raw.decode("utf-8", "replace")
        if _leaks(text, lits):
            hits.append(rel)
        elif word.search(text):
            hits.append(f"{rel} (the login name)")
    return hits


def bundle_version(installed: str) -> str:
    """The bundle version: the ``MAJOR.MINOR`` of the installed GIPPA."""
    m = re.match(r"(\d+)\.(\d+)", installed)
    if m is None:
        raise BuildError(f"cannot read the GIPPA version {installed!r}")
    return f"{m.group(1)}.{m.group(2)}"


def _git(*args: str) -> subprocess.CompletedProcess[str] | None:
    env = {k: v for k, v in os.environ.items() if k not in _GIT_ENV}
    try:
        return subprocess.run(["git", "-C", str(HERE), *args], capture_output=True, text=True,
                              check=False, timeout=60, env=env, encoding="utf-8")  # fmt: skip
    except (OSError, subprocess.SubprocessError):
        return None


def git_state() -> tuple[str, bool | None]:
    """(HEAD sha, dirty) of the checkout holding this script; ("", None) outside a checkout.
    Dirty = a tracked change or an untracked file other than Python byte code."""
    head = _git("rev-parse", "HEAD")
    if head is None or head.returncode != 0:
        return "", None
    status = _git("status", "--porcelain")
    if status is None or status.returncode != 0:
        return head.stdout.strip(), None
    changed = [ln for ln in status.stdout.splitlines()
               if ln.strip() and "__pycache__" not in ln and not ln.endswith(".pyc")]  # fmt: skip
    return head.stdout.strip(), bool(changed)


def build_commit(allow_dirty: bool) -> str:
    """The commit recorded in the bundle: ``<sha>``, or ``<sha>-dirty`` with ``--allow-dirty``."""
    sha, dirty = git_state()
    if dirty is None:
        if not allow_dirty:
            raise BuildError(f"{HERE} is not in a git checkout, so it cannot be proven clean "
                             "(not dirty); pass --allow-dirty to build anyway")  # fmt: skip
        return f"{sha}-dirty" if sha else "unknown-dirty"
    if dirty:
        if not allow_dirty:
            raise BuildError(
                "the checkout is dirty (uncommitted changes); commit them or pass "
                "--allow-dirty (the commit is then recorded as <sha>-dirty)"
            )
        return f"{sha}-dirty"
    return sha


def routes() -> list[dict[str, Any]]:
    """``gippa/data/example_routes.yaml`` (the single source of the route commands)."""
    import yaml

    text = package_files("gippa.data").joinpath("example_routes.yaml").read_text(encoding="utf-8")
    out: list[dict[str, Any]] = yaml.safe_load(text)["routes"]
    return out


def route_paths(rs: Iterable[Mapping[str, Any]]) -> set[str]:
    return {p for r in rs for c in r["commands"] for p in re.findall(r"\{dir\}/(\S+)", str(c))}


def _tar_info(name: str, path: Path) -> tarfile.TarInfo:
    info = tarfile.TarInfo(name)
    info.mtime, info.uid, info.gid, info.uname, info.gname = 0, 0, 0, "", ""
    if path.is_dir():
        info.type, info.mode = tarfile.DIRTYPE, 0o755
    else:
        info.size, info.mode = path.stat().st_size, 0o644
    return info


def deterministic_tar(out: Path, base: Path, members: Iterable[Path]) -> None:
    """``members`` (files or directories, no symlinks) as ``out``, named relative to
    ``base.parent``: sorted, mtime/uid/gid 0, empty owners, modes 0644/0755, gzip with no file
    name and mtime 0. Written to ``out.partial`` first, then renamed."""
    named = sorted((m.relative_to(base.parent).as_posix(), m) for m in members)
    bad = [n for n, m in named if m.is_symlink() or not (m.is_file() or m.is_dir())]
    if bad:
        raise BuildError(f"not a regular file or directory: {bad}")
    tmp = out.with_name(out.name + ".partial")
    out.parent.mkdir(parents=True, exist_ok=True)
    try:
        with (
            tmp.open("wb") as raw,
            gzip.GzipFile(fileobj=raw, mode="wb", mtime=0, filename="") as gz,
            tarfile.open(fileobj=gz, mode="w", format=tarfile.PAX_FORMAT) as tar,
        ):
            for name, m in named:
                info = _tar_info(name, m)
                if m.is_dir():
                    tar.addfile(info)
                else:
                    with m.open("rb") as fh:
                        tar.addfile(info, fh)
        os.replace(tmp, out)
    finally:
        tmp.unlink(missing_ok=True)


def checksums(top: Path) -> str:
    """``sha256sum`` lines of every file under ``top`` but ``CHECKSUMS.sha256``, sorted by path."""
    rels = sorted(p.relative_to(top).as_posix() for p in top.rglob("*") if p.is_file())
    return "".join(f"{sha256(top / r)}  {r}\n" for r in rels if r != "CHECKSUMS.sha256")


def _check_routes(a: Mapping[str, Any], c: Mapping[str, Mapping[str, int]], n_b: int, n_d: int,
                  n_prot: int, e: Mapping[str, Any]) -> None:  # fmt: skip
    """Rule 7: no empty route, route A keeps a decoy, route E a target, a decoy and an ORF row."""
    empty = [k for k, n in (("region proteins", n_prot), ("A rows", a["target_rows_kept"]),
                            ("A decoys", a["decoy_rows_kept"]), ("B proteins", n_b),
                            ("D rows", n_d), ("E target rows", e["target_rows_kept"]),
                            ("E decoys", e["decoy_rows_kept"]),
                            ("E rows with an ORF ID", e["rows_with_orf_ids"]))
             if n == 0]  # fmt: skip
    empty += [f"C rows of {r}" for r, v in c.items() if v["rows_kept"] == 0]
    if empty:
        raise BuildError("empty route (rule 7): " + ", ".join(empty))


def _rel(path: Path, root: Path) -> str:
    """A source path relative to ``$GIPPA_ROOT``, lexically (R12: no absolute path and no
    ``..`` in the bundle). A source that resolves outside ``$GIPPA_ROOT`` (a symlink) fails."""
    rel = path.relative_to(root).as_posix()
    if not path.exists():
        raise BuildError(f"source $GIPPA_ROOT/{rel} does not exist")
    if not path.resolve().is_relative_to(root.resolve()):
        raise BuildError(f"source $GIPPA_ROOT/{rel} resolves outside $GIPPA_ROOT (a symlink?); "
                         "copy the pilot under $GIPPA_ROOT")  # fmt: skip
    return rel


def static_digests(root: Path) -> dict[str, dict[str, str]]:
    """PROVENANCE.md's tool and table digests: the file's sha256 where it exists under
    ``$GIPPA_ROOT`` (``root:``) or the checkout (``repo:``), the pinned value otherwise."""
    out: dict[str, dict[str, str]] = {}
    for key, (label, where, pinned) in STATIC_DIGESTS.items():
        digest = pinned
        if where:
            kind, rel = where.split(":", 1)
            p = (root if kind == "root" else REPO) / rel
            if p.is_file():
                digest = sha256(p)
                if digest != pinned:
                    print(f"build_example: warning: {label} sha256 differs from the pinned value",
                          file=sys.stderr)  # fmt: skip
        out[key] = {"label": label, "sha256": digest}
    return out


def region_info(start: int, end: int) -> dict[str, Any]:
    """``SELECTION.json``'s description of the region ``[start, end)`` of ``CHROM``."""
    return {
        "contig": CHROM,
        "start_0based": start,
        "end_exclusive": end,
        "length_bp": end - start,
        "one_based_closed": f"{CHROM}:{start + 1}-{end}",
        "output_contig_id": f"nissle1917_{CHROM.replace('.', '_')}_{start}_{end}",
        "topology": "linear (no tag)",
    }


def build(
    start: int, end: int, name: str, stage: Path, root: Path, *, version: str, commit: str,
    expected: Path | None = None, expected_commit: str = "",
) -> dict[str, Any]:  # fmt: skip
    """The release tree ``stage/gippa-example-<version>/`` and the per-route sizing tars in
    ``stage/bundles/`` (``stage`` must not exist); returns the counts and the tree. Route E's
    region is the module's ``START_E``/``END_E`` (read here, at call time)."""
    import gippa

    start_e, end_e = START_E, END_E
    if not 0 <= start_e < end_e:
        raise BuildError(f"route E's region [{start_e}, {end_e}) is empty or negative")

    src = sources(root)
    src_rel = {k: _rel(v, root) for k, v in src.items()}  # before anything is read (L2)
    out = stage / f"gippa-example-{version}"
    for sub in ("genome", "A_diann", "B_proteins", "C_cascadia", "D_generic", "E_sixframe"):
        (out / sub).mkdir(parents=True)
    (stage / "bundles").mkdir()

    genome = read_fasta(src["genome"])
    if CHROM not in genome:
        raise BuildError(f"{CHROM} is not in $GIPPA_ROOT/{src_rel['genome']}")
    chrom = genome[CHROM].upper()
    if end > len(chrom):
        raise BuildError(f"--end {end} is beyond {CHROM} ({len(chrom)} bp)")
    if end_e > len(chrom):
        raise BuildError(f"route E's region end {end_e} is beyond {CHROM} ({len(chrom)} bp)")
    region, region_e = region_info(start, end), region_info(start_e, end_e)
    region_nt = chrom[start:end]
    g_path = out / "genome" / f"{name}.fna.gz"
    write_fasta(g_path, [(region["output_contig_id"], region_nt)])
    region_nt_e, contig_e = chrom[start_e:end_e], region_e["output_contig_id"]
    g_path_e = out / "genome" / f"{name}_e.fna.gz"
    write_fasta(g_path_e, [(contig_e, region_nt_e)])

    prot = region_proteins(src["db_id_map"], start, end)
    keep_ids = set(prot["db_id"])

    # Route A
    sub_t, a_counts = subset_diann(src["diann_report"], keep_ids)
    sub_t, a_scrubbed = scrub_paths(sub_t)
    a_path = out / "A_diann" / "report.parquet"
    pq.write_table(sub_t, a_path, compression="zstd")
    sub = sub_t.to_pandas()
    db = read_fasta(src["strain_db"])

    # Route B
    ab, n_shared = abundance(sub)
    missing = [pid for pid in ab["protein_id"] if pid not in db]
    if missing:
        raise BuildError(f"route B proteins missing from $GIPPA_ROOT/{src_rel['strain_db']}: "
                         + ", ".join(missing))  # fmt: skip
    ab_path = out / "B_proteins" / "abundance.tsv"
    write_table(ab, ab_path)
    f_path = out / "B_proteins" / "proteins.faa"
    write_fasta(f_path, [(pid, db[pid]) for pid in ab["protein_id"]])

    # Route C
    frames_il = [f.replace("I", "L") for f in six_frames(region_nt)]
    c_counts: dict[str, dict[str, int]] = {}
    for r in RUNS:
        kept, header, n = subset_ssl(src[f"ssl_{r}"], frames_il)
        write_tsv_rows(out / "C_cascadia" / f"{r}.ssl", header, kept)
        c_counts[r] = {"rows_total": n, "rows_kept": len(kept)}

    # Route D
    gen = generic_table(sub)
    d_path = out / "D_generic" / "peptides.tsv"
    write_table(gen, d_path)

    # Route E (rules 3E, 8, 9)
    ent = region_entries(src["a4_db_id_map"], start_e, end_e)
    orfs = ent[ent["kind"] == "orf"]
    orf_ids = set(orfs["db_id"])
    e_t, e_kept = subset_diann(src["a4_report"], set(ent["db_id"]))
    e_counts: dict[str, Any] = dict(e_kept)
    e_counts["rows_with_orf_ids"] = sum(
        1 for x in e_t.column("Protein.Ids").to_pylist() if set(str(x).split(";")) & orf_ids
    )
    e_counts["rows_per_run"] = {
        str(k): int(v) for k, v in e_t.column("Run").to_pandas().value_counts().sort_index().items()
    }
    _check_routes(a_counts, c_counts, len(ab), len(gen), len(prot), e_counts)
    e_t, e_scrubbed = scrub_paths(e_t)
    e_t, e_rw_cols, n_rw = orf_id_rewrite(
        e_t, orfs, region_nt_e, start_e, contig_e, read_fasta(src["a4_strain_db"])
    )
    e_path = out / "E_sixframe" / "report.parquet"
    pq.write_table(e_t, e_path, compression="zstd")

    targets = sub[sub["Decoy"] == 0]
    counts: dict[str, Any] = {
        "region_proteins": len(prot),
        "region_proteins_cds": int((prot["kind"] == "cds").sum()),
        "region_proteins_sorf": int((prot["kind"] == "sorf").sum()),
        "A": {
            **a_counts,
            "distinct_target_peptides": int(targets["Stripped.Sequence"].nunique()),
            "distinct_target_peptides_q01": int(
                targets.loc[targets["Global.Q.Value"] <= 0.01, "Stripped.Sequence"].nunique()
            ),
            "proteins_with_rows": len(
                {p for x in targets["Protein.Ids"] for p in str(x).split(";")}
            ),
            "rows_per_run": {k: int(v) for k, v in sub["Run"].value_counts().sort_index().items()},
        },
        "B": {
            "abundance_rows": len(ab),
            "fasta_records": len(ab),
            "proteins_in_several_groups": n_shared,
        },
        "C": c_counts,
        "D": {"rows": len(gen), "decoy_rows": int((gen["is_decoy"] == "true").sum())},
        "E": {
            **e_counts,
            "region_entries": len(ent),
            "region_orf_entries": len(orfs),
            "orf_entries_named": n_rw,
        },
    }
    rewritten = {"A_diann/report.parquet": a_scrubbed, "C_cascadia/*.ssl": ["file"],
                 "E_sixframe/report.parquet": e_scrubbed}  # fmt: skip
    rewritten = {k: v for k, v in rewritten.items() if v}
    src_meta = {k: {"path": src_rel[k], "sha256": sha256(v)} for k, v in src.items()}
    sel: dict[str, Any] = {
        "name": name,
        "region": region,
        "region_e": region_e,
        "rule": __doc__.split("Selection rule (recorded verbatim in SELECTION.json):")[1]
        .split("Usage")[0]
        .strip(),
        "fdr_note": (
            "q-values, PEPs and MaxLFQ in routes A, B and D are DIA-NN's, computed on the "
            "whole-proteome A2 search (4,685 strain proteins + contaminants, three runs); route "
            "E's are DIA-NN's on the whole six-frame search A4 (147,753 strain records + 381 "
            "contaminants, the same runs); Cascadia scores are per spectrum. Nothing was "
            "re-searched or re-scored for the subset."
        ),
        "counts": counts,
        "rewritten_columns": rewritten,
        "orf_id_rewrite": {
            "rule": "route E ORF IDs in region coordinates (see rule)",
            "file": "E_sixframe/report.parquet",
            "contig": contig_e,
            "offset": start_e,
            "columns": e_rw_cols,
            "distinct_ids": n_rw,
        },
        "sources": src_meta,
        "sources_relative_to": "$GIPPA_ROOT",
    }
    route_files = {
        "A_diann": [a_path],
        "B_proteins": [f_path, ab_path],
        "C_cascadia": [out / "C_cascadia" / f"{r}.ssl" for r in RUNS],
        "D_generic": [d_path],
        "E_sixframe": [g_path_e, e_path],
    }
    data_files = [g_path, *(p for ps in route_files.values() for p in ps)]
    sel["outputs"] = {
        p.relative_to(out).as_posix(): {"bytes": p.stat().st_size, "sha256": sha256(p)}
        for p in data_files
    }
    sel_path = out / "SELECTION.json"
    sel_path.write_text(json.dumps(sel, indent=2) + "\n", encoding="utf-8")

    rs = routes()
    absent = sorted(r for r in route_paths(rs) if not (out / r).is_file())
    if absent:
        print(f"build_example: warning: example_routes.yaml names files this build did not write "
              f"(--name {name}?): {absent}", file=sys.stderr)  # fmt: skip
    if expected is not None:
        for p in sorted(expected.rglob("*")):
            if p.is_file() and not p.is_symlink():
                dest = out / "expected" / p.relative_to(expected)
                dest.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(p, dest)
    installed = str(gippa.__version__)
    exp_note = (f"from {expected_commit or 'an unrecorded commit'}" if expected is not None
                else "not included in this build")  # fmt: skip
    (out / "README.md").write_text(T.readme(version, region["one_based_closed"],
                                            region["length_bp"], rs, expected is not None,
                                            region_e, counts["E"], f"{name}.fna.gz",
                                            f"{name}_e.fna.gz"), encoding="utf-8")  # fmt: skip
    (out / "PROVENANCE.md").write_text(T.provenance(
        version, region, src_meta, counts,
        {"commit": commit, "gippa_version": installed, "expected": exp_note}, rewritten,
        static_digests(root), a4_extras(root), n_rw, region_e, sha256(g_path_e),
    ), encoding="utf-8")  # fmt: skip
    (out / "LICENSE-DATA.md").write_text(T.LICENSE_DATA, encoding="utf-8")
    meta = {
        "bundle_version": version,
        "gippa_version": installed,
        "build_commit": commit,
        "expected_commit": expected_commit if expected is not None else None,
        "region": region["one_based_closed"],
        "region_0based_half_open": [CHROM, start, end],
        "region_e": region_e["one_based_closed"],
        "region_e_0based_half_open": [CHROM, start_e, end_e],
        "name": name,
        "counts": counts,
        "tools": T.TOOLS,
        "routes": [r["id"] for r in rs],
    }
    (out / "bundle.json").write_text(json.dumps(meta, indent=2, sort_keys=True) + "\n",
                                     encoding="utf-8")  # fmt: skip
    (out / "CHECKSUMS.sha256").write_text(checksums(out), encoding="utf-8")

    leaks = personal_paths(out, [root, root.resolve(), Path.home(), stage.resolve()],
                           login_names())  # fmt: skip
    if leaks:
        raise BuildError("personal or build path in the bundle (R12): " + "; ".join(leaks))

    sizes: dict[str, int] = {}
    for route, ps in route_files.items():
        tar = stage / "bundles" / f"{name}_{route}.tar.gz"
        shared = [] if route == "E_sixframe" else [g_path]  # route E ships its own genome
        deterministic_tar(tar, out, [*shared, sel_path, *ps])
        sizes[route] = tar.stat().st_size
    (stage / "bundles" / "sizes.json").write_text(json.dumps(sizes, indent=2) + "\n",
                                                  encoding="utf-8")  # fmt: skip
    return {"selection": counts, "route_bytes": sizes, "tree": out}


def pack(top: Path, out: Path) -> tuple[str, int, int]:
    """``top`` as the archive ``out``; (sha256, size, file count)."""
    members = [top, *top.rglob("*")]
    deterministic_tar(out, top, members)
    return sha256(out), out.stat().st_size, sum(1 for m in members if m.is_file())


def release(top: Path, archive: Path, force: bool) -> tuple[str, int, int]:
    """Pack ``top`` as ``archive`` (left read-only); (sha256, size, file count). An existing
    archive with other content is kept unless ``force``."""
    tmp = archive.with_name(archive.name + ".new")
    try:
        digest, size, n_files = pack(top, tmp)
        if archive.exists() and not force and sha256(archive) != digest:
            raise BuildError(f"{archive} exists with another sha256; remove it or pass --force")
        os.replace(tmp, archive)
    finally:
        tmp.unlink(missing_ok=True)
    archive.chmod(0o444)
    return digest, size, n_files


def install(stage: Path, work: Path, version: str) -> None:
    """Move the staged tree and sizing tars to ``work/`` (after the release succeeded)."""
    work.mkdir(parents=True, exist_ok=True)
    for part in (f"gippa-example-{version}", "bundles"):
        dest = work / part
        if dest.exists():
            old = work / f"{part}.partial-old-{os.getpid()}"
            os.replace(dest, old)
            shutil.rmtree(old)
        os.replace(stage / part, dest)


def recorded_row(tsv: Path, version: str) -> dict[str, str]:
    """The ``bundles.tsv`` row of ``version``."""
    if not tsv.is_file():
        raise BuildError(f"--bundles-tsv {tsv}: no such file")
    with tsv.open(encoding="utf-8", newline="") as fh:
        rows = [r for r in csv.DictReader(fh, delimiter="\t") if r.get("version") == version]
    if len(rows) != 1:
        raise BuildError(f"{tsv}: {len(rows)} rows for version {version}, need exactly 1")
    if not rows[0].get("sha256"):
        raise BuildError(f"{tsv}: the {version} row has no sha256")
    return rows[0]


def check(a: argparse.Namespace, root: Path, version: str, commit: str) -> int:
    """``--check``: rebuild in a temporary dir and compare with the ``bundles.tsv`` row."""
    row = recorded_row(a.bundles_tsv or BUNDLES_TSV, version)
    want = row["sha256"].removeprefix("sha256:")
    with tempfile.TemporaryDirectory(prefix="gippa-example-check-") as tmp:
        stage = Path(tmp) / "stage"
        res = build(a.start, a.end, a.name, stage, root, version=version, commit=commit,
                    expected=a.expected, expected_commit=a.expected_commit)  # fmt: skip
        digest, size, _ = pack(res["tree"], Path(tmp) / T.archive_name(version))
    if row.get("build_commit") and row["build_commit"] != commit:
        print(f"build_example: note: bundles.tsv records build_commit {row['build_commit']}, "
              f"this checkout is {commit}", file=sys.stderr)  # fmt: skip
    problems = []
    if digest != want:
        problems.append(f"sha256 mismatch: rebuilt {digest}, bundles.tsv {want}")
    if row.get("size_bytes") and row["size_bytes"] != str(size):
        problems.append(f"size mismatch: rebuilt {size}, bundles.tsv {row['size_bytes']}")
    if problems:
        print("build_example: --check failed: " + "; ".join(problems), file=sys.stderr)
        return 1
    print(f"--check ok: gippa-example-{version} sha256 {digest}, {size} bytes")
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--start", type=int, default=START,
                    help=f"0-based start of the shared region (routes A-D) on {CHROM} "
                         f"(default {START}); route E's region is START_E")  # fmt: skip
    ap.add_argument("--end", type=int, default=END,
                    help=f"exclusive end of the shared region (default {END})")  # fmt: skip
    ap.add_argument("--name", default=NAME, help=f"default {NAME}")
    ap.add_argument("--out", type=Path, default=None, help="work dir; default $GIPPA_ROOT/example")
    ap.add_argument("--bundle-out", type=Path, default=None,
                    help="where the archive goes; default $GIPPA_ROOT/example")  # fmt: skip
    ap.add_argument("--expected", type=Path, default=None,
                    help="copy DIR/* into the bundle's expected/ (none without it)")  # fmt: skip
    ap.add_argument("--expected-commit", default="", help="the commit that produced --expected")
    ap.add_argument("--force", action="store_true", help="replace an archive with another sha256")
    ap.add_argument("--allow-dirty", action="store_true",
                    help="build from a dirty checkout (records <sha>-dirty)")  # fmt: skip
    ap.add_argument("--check", action="store_true",
                    help="rebuild in a temporary dir and compare with the bundles.tsv row; "
                         "writes nothing to --out or --bundle-out")  # fmt: skip
    ap.add_argument("--bundles-tsv", type=Path, default=None,
                    help="--check reads it (default validation/example/bundles.tsv)")  # fmt: skip
    a = ap.parse_args(argv)
    if not 0 <= a.start < a.end:
        ap.error("need 0 <= start < end")
    if a.expected is not None and not a.expected.is_dir():
        ap.error(f"--expected {a.expected}: not a directory")
    import gippa

    root = _root()
    stage: Path | None = None
    try:
        version = bundle_version(str(gippa.__version__))
        commit = build_commit(a.allow_dirty)
        if a.check:
            return check(a, root, version, commit)
        out_root = a.out or root / "example"
        work = out_root / a.name
        out_root.mkdir(parents=True, exist_ok=True)
        stage = Path(tempfile.mkdtemp(prefix=f".{a.name}.partial-", dir=out_root))
        shutil.rmtree(stage)  # build() creates it
        res = build(a.start, a.end, a.name, stage, root, version=version, commit=commit,
                    expected=a.expected, expected_commit=a.expected_commit)  # fmt: skip
        archive = (a.bundle_out or root / "example") / T.archive_name(version)
        digest, size, n_files = release(res["tree"], archive, a.force)
        install(stage, work, version)
    except BuildError as exc:
        print(f"build_example: {exc}", file=sys.stderr)
        return 1
    finally:
        if stage is not None and stage.exists():
            shutil.rmtree(stage)
    job = os.environ.get("SLURM_JOB_ID", "")
    print(json.dumps({"selection": res["selection"], "route_bytes": res["route_bytes"]}, indent=2))
    print(f"\narchive: {archive}\nsha256: {digest}\nsize_bytes: {size}\nfiles: {n_files}")
    print("\n# db_sources.yaml row\n" + T.db_sources_row(version, digest, size, job))
    print("\n# validation/example/bundles.tsv row\n" + "\t".join(T.TSV_COLUMNS))
    print(T.bundles_tsv_row(version, digest, size, commit, job,
                            a.expected_commit if a.expected is not None else ""))  # fmt: skip
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
