"""Build the quick-start example bundle: one Nissle 1917 chromosome region, four input routes.

The example is a cut-down subset of the M10 pilot (validation/pilot/README.md, pass 3): E. coli
Nissle 1917 (GCF_021559835.1, chromosome NZ_CP058217.1), PXD078459 runs WT_1-WT_3. Nothing is
re-searched. The script subsets the pilot's existing outputs to one contiguous region and writes
the release tree OUT/<name>/gippa-example-<X.Y>/ (spec M12-01 section 2 > Layout):

    genome/<name>.fna.gz          the region as one linear contig (no topology tag)
    A_diann/report.parquet        DIA-NN 2.3.1 report rows whose proteins all lie in the region
    B_proteins/proteins.faa       the region proteins DIA-NN quantified (A2 strain database)
    B_proteins/abundance.tsv      their PG.MaxLFQ (median over runs), as `gippa convert` computes it
    C_cascadia/WT_{1,2,3}.ssl     Cascadia rows whose peptide occurs in the region's six frames
    D_generic/peptides.tsv        a generic peptide table, peptide level, built from route A's rows
    SELECTION.json                the rule, the counts and the sha256 of every source and output
    README.md, PROVENANCE.md, LICENSE-DATA.md, bundle.json, CHECKSUMS.sha256, expected/ (optional)

then packs it as BUNDLE_OUT/gippa-example-<X.Y>.tar.gz (deterministic: sorted members, mtime,
uid and gid 0, empty owner names, modes 0644/0755, gzip without name and mtime) and prints the
archive's sha256 and size, its db_sources.yaml row and its validation/example/bundles.tsv row.
Sizing tars per route stay in OUT/<name>/bundles/, outside the release tree.

Selection rule (recorded verbatim in SELECTION.json):

* Region: NZ_CP058217.1 [start, end) in 0-based half-open coordinates.
* Region proteins: records of the pass-3 A2 strain database (`db_id_map.tsv`, 1-based closed)
  on NZ_CP058217.1 with start - 1 >= region start and end <= region end (CDS and sORF).
* Route A: a DIA-NN row is kept iff every ID in `Protein.Ids` (split on ';') is a region protein.
  This applies to target and decoy rows alike (DIA-NN decoys carry the target IDs with Decoy = 1),
  so the decoys keep their `protein_ids`. Rows with a contaminant (`Cont_`) or out-of-region ID are
  dropped. All 71 columns are kept. q-values, PEPs and MaxLFQ are DIA-NN's, computed on the whole
  proteome search; they are not recomputed for the subset.
* Route B: the protein groups of route A's target rows with a non-missing, non-zero PG.MaxLFQ;
  value = median over runs of PG.MaxLFQ (the `gippa convert --abundance-out` rule). A group "A;B" is
  written as one row per member ID; the FASTA holds each member's sequence from the strain database.
* Route C: a Cascadia row is kept iff its stripped sequence (I and L equated) occurs as an exact
  substring of a stop-free stretch of the region's six-frame translation (genetic code 11). The
  `file` column is reduced to the mzML base name (the run name is its stem).
* Route D: route A's rows grouped by (stripped peptide, run): q_value = min Global.Q.Value,
  pep = min PEP, intensity = sum of Precursor.Normalised (targets), protein_ids = the union,
  is_decoy = all rows decoy; modified_peptide = the first Modified.Sequence in ProForma
  (`(UniMod:n)` -> `[UNIMOD:n]`). Columns: peptide, modified_peptide, protein_ids, q_value, pep,
  is_decoy, intensity, run.

* Scrubbing (R12): the `.ssl` `file` column and every DIA-NN text column holding an absolute
  path are reduced to base names; SELECTION.json lists the rewritten columns. No other value is
  changed.
* No empty route: the build fails if any route keeps zero rows or zero proteins, or if route A
  keeps no decoy.

Usage (venv python, PYTHONPATH=<checkout>/src; small, runs in seconds). Defaults: region
NZ_CP058217.1:1965001-2165000, --name nissle_r200k, --out and --bundle-out $GIPPA_ROOT/example,
bundle version = the installed GIPPA MAJOR.MINOR:

    python scripts/example/build_example.py [--expected DIR] [--bundle-out DIR]
"""

from __future__ import annotations

import argparse
import csv
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
HERE = Path(__file__).resolve().parent
_ABS = re.compile(r"^(/|[A-Za-z]:[\\/]|\\\\)")  # POSIX, Windows drive or UNC


class BuildError(Exception):
    """The inputs cannot give a valid bundle (rule 7, an existing different archive)."""


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


def sources(root: Path) -> dict[str, Path]:
    p3 = root / "pilot" / "pass3"
    a3 = root / "pilot" / "A3_cascadia_cpu"
    return {
        "genome": p3 / "genome" / "nissle1917_wt.fna.gz",
        "db_id_map": p3 / "db_a2" / "db_id_map.tsv",
        "strain_db": p3 / "db_a2" / "strain_db.faa.gz",
        "diann_report": p3 / "A2" / "A2_report.parquet",
        **{f"ssl_{r}": a3 / r / f"{r}.ssl" for r in RUNS},
    }


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_fasta(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    opener = gzip.open if path.suffix == ".gz" else open
    name, seq = None, []
    with opener(path, "rt") as fh:
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
    data = buf.getvalue().encode()
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
    m = pd.read_csv(db_map, sep="\t", dtype={"db_id": str})
    inside = (m["contig"] == CHROM) & (m["start"] - 1 >= start) & (m["end"] <= end)
    return m[inside].reset_index(drop=True)


def subset_diann(report: Path, keep_ids: set[str]) -> tuple[pa.Table, dict[str, int]]:
    t = pq.read_table(report)
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


def abundance(sub: pd.DataFrame) -> pd.DataFrame:
    t = sub[(sub["Decoy"] == 0)]
    t = t[t["PG.MaxLFQ"].notna() & (t["PG.MaxLFQ"] > 0)]
    per_run = t.groupby(["Protein.Group", "Run"])["PG.MaxLFQ"].first().reset_index()
    med = per_run.groupby("Protein.Group")["PG.MaxLFQ"].median()
    rows = []
    for group, v in sorted(med.items()):
        for pid in str(group).split(";"):
            rows.append((pid, float(v)))
    return pd.DataFrame(rows, columns=["protein_id", "maxlfq"])


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
    with path.open() as fh:
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
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=header, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def scrub_paths(t: pa.Table) -> tuple[pa.Table, list[str]]:
    """Rule 8: every text column holding an absolute path, with each such value reduced to its
    base name (``/`` or backslash separated); the names of the rewritten columns."""
    rewritten = []
    for i, field in enumerate(t.schema):
        if not (pa.types.is_string(field.type) or pa.types.is_large_string(field.type)):
            continue
        vals = t.column(i).to_pylist()
        if not any(isinstance(v, str) and _ABS.match(v) for v in vals):
            continue
        new = [re.split(r"[\\/]", v)[-1] if isinstance(v, str) and _ABS.match(v) else v
               for v in vals]  # fmt: skip
        t = t.set_column(i, field, pa.array(new, type=field.type))
        rewritten.append(field.name)
    return t, rewritten


def bundle_version(installed: str) -> str:
    """The bundle version: the ``MAJOR.MINOR`` of the installed GIPPA."""
    m = re.match(r"(\d+)\.(\d+)", installed)
    if m is None:
        raise BuildError(f"cannot read the GIPPA version {installed!r}")
    return f"{m.group(1)}.{m.group(2)}"


def git_commit() -> str:
    """The commit of the checkout holding this script ("" outside a git checkout)."""
    try:
        res = subprocess.run(["git", "-C", str(HERE), "rev-parse", "HEAD"], capture_output=True,
                             text=True, check=False, timeout=30)  # fmt: skip
    except (OSError, subprocess.SubprocessError):
        return ""
    return res.stdout.strip() if res.returncode == 0 else ""


def routes() -> list[dict[str, Any]]:
    """``gippa/data/example_routes.yaml`` (the single source of the route commands)."""
    import yaml

    raw = yaml.safe_load(package_files("gippa.data").joinpath("example_routes.yaml").read_text())
    out: list[dict[str, Any]] = raw["routes"]
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


def checksums(top: Path) -> str:
    """``sha256sum`` lines of every file under ``top`` but ``CHECKSUMS.sha256``, sorted by path."""
    rels = sorted(p.relative_to(top).as_posix() for p in top.rglob("*") if p.is_file())
    return "".join(f"{sha256(top / r)}  {r}\n" for r in rels if r != "CHECKSUMS.sha256")


def _check_routes(a: Mapping[str, Any], c: Mapping[str, Mapping[str, int]], n_b: int, n_d: int,
                  n_prot: int) -> None:  # fmt: skip
    """Rule 7: no empty route, and route A keeps a decoy."""
    empty = [k for k, n in (("region proteins", n_prot), ("A rows", a["target_rows_kept"]),
                            ("A decoys", a["decoy_rows_kept"]), ("B proteins", n_b),
                            ("D rows", n_d)) if n == 0]  # fmt: skip
    empty += [f"C rows of {r}" for r, v in c.items() if v["rows_kept"] == 0]
    if empty:
        raise BuildError("empty route (rule 7): " + ", ".join(empty))


def _rel(path: Path, root: Path) -> str:
    """A source path relative to ``$GIPPA_ROOT`` (R12: no absolute path in the bundle)."""
    try:
        return path.resolve().relative_to(root.resolve()).as_posix()
    except ValueError:
        return Path(os.path.relpath(path.resolve(), root.resolve())).as_posix()


def build(
    start: int, end: int, name: str, out_root: Path, root: Path, *, version: str,
    expected: Path | None = None, expected_commit: str = "",
) -> dict[str, Any]:  # fmt: skip
    """The release tree ``out_root/name/gippa-example-<version>/`` (rebuilt from scratch) and the
    per-route sizing tars in ``out_root/name/bundles/``; returns the counts and the tree."""
    import gippa

    src = sources(root)
    work = out_root / name
    out = work / f"gippa-example-{version}"
    if out.exists():
        shutil.rmtree(out)
    for sub in ("genome", "A_diann", "B_proteins", "C_cascadia", "D_generic"):
        (out / sub).mkdir(parents=True, exist_ok=True)
    (work / "bundles").mkdir(parents=True, exist_ok=True)

    genome = read_fasta(src["genome"])
    chrom = genome[CHROM].upper()
    region_nt = chrom[start:end]
    contig = f"nissle1917_{CHROM.replace('.', '_')}_{start}_{end}"
    g_path = out / "genome" / f"{name}.fna.gz"
    write_fasta(g_path, [(contig, region_nt)])

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
    ab = abundance(sub)
    ab_path = out / "B_proteins" / "abundance.tsv"
    ab.to_csv(ab_path, sep="\t", index=False, lineterminator="\n", float_format="%.6g")
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
    gen.to_csv(d_path, sep="\t", index=False, lineterminator="\n", float_format="%.6g")
    _check_routes(a_counts, c_counts, len(ab), len(gen), len(prot))

    targets = sub[sub["Decoy"] == 0]
    region = {
        "contig": CHROM,
        "start_0based": start,
        "end_exclusive": end,
        "length_bp": end - start,
        "one_based_closed": f"{CHROM}:{start + 1}-{end}",
        "output_contig_id": contig,
        "topology": "linear (no tag)",
    }
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
        "B": {"abundance_rows": len(ab), "fasta_records": len(ab)},
        "C": c_counts,
        "D": {"rows": len(gen), "decoy_rows": int((gen["is_decoy"] == "true").sum())},
    }
    rewritten = {"A_diann/report.parquet": a_scrubbed, "C_cascadia/*.ssl": ["file"]}
    rewritten = {k: v for k, v in rewritten.items() if v}
    src_meta = {k: {"path": _rel(v, root), "sha256": sha256(v)} for k, v in src.items()}
    sel: dict[str, Any] = {
        "name": name,
        "region": region,
        "rule": __doc__.split("Selection rule (recorded verbatim in SELECTION.json):")[1]
        .split("Usage")[0]
        .strip(),
        "fdr_note": (
            "q-values, PEPs and MaxLFQ in routes A, B and D are DIA-NN's, computed on the "
            "whole-proteome A2 search (4,685 strain proteins + contaminants, three runs); Cascadia "
            "scores are per spectrum. Nothing was re-searched or re-scored for the subset."
        ),
        "counts": counts,
        "rewritten_columns": rewritten,
        "sources": src_meta,
        "sources_relative_to": "$GIPPA_ROOT",
    }
    route_files = {
        "A_diann": [a_path],
        "B_proteins": [f_path, ab_path],
        "C_cascadia": [out / "C_cascadia" / f"{r}.ssl" for r in RUNS],
        "D_generic": [d_path],
    }
    data_files = [g_path, *(p for ps in route_files.values() for p in ps)]
    sel["outputs"] = {
        p.relative_to(out).as_posix(): {"bytes": p.stat().st_size, "sha256": sha256(p)}
        for p in data_files
    }
    sel_path = out / "SELECTION.json"
    sel_path.write_text(json.dumps(sel, indent=2) + "\n")

    rs = routes()
    missing = sorted(r for r in route_paths(rs) if not (out / r).is_file())
    if missing:
        print(f"build_example: warning: example_routes.yaml names files this build did not write "
              f"(--name {name}?): {missing}", file=sys.stderr)  # fmt: skip
    if expected is not None:
        for p in sorted(expected.rglob("*")):
            if p.is_file() and not p.is_symlink():
                dest = out / "expected" / p.relative_to(expected)
                dest.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(p, dest)
    commit = git_commit()
    installed = str(gippa.__version__)
    exp_note = (f"from {expected_commit or 'an unrecorded commit'}" if expected is not None
                else "not included in this build")  # fmt: skip
    (out / "README.md").write_text(T.readme(version, region["one_based_closed"], rs,
                                            expected is not None))  # fmt: skip
    (out / "PROVENANCE.md").write_text(T.provenance(
        version, region, src_meta, counts,
        {"commit": commit, "gippa_version": installed, "expected": exp_note}, rewritten,
    ))  # fmt: skip
    (out / "LICENSE-DATA.md").write_text(T.LICENSE_DATA)
    meta = {
        "bundle_version": version,
        "gippa_version": installed,
        "build_commit": commit,
        "expected_commit": expected_commit if expected is not None else None,
        "region": region["one_based_closed"],
        "region_0based_half_open": [CHROM, start, end],
        "name": name,
        "counts": counts,
        "tools": T.TOOLS,
        "routes": [r["id"] for r in rs],
    }
    (out / "bundle.json").write_text(json.dumps(meta, indent=2, sort_keys=True) + "\n")
    (out / "CHECKSUMS.sha256").write_text(checksums(out))

    sizes: dict[str, int] = {}
    for route, ps in route_files.items():
        tar = work / "bundles" / f"{name}_{route}.tar.gz"
        deterministic_tar(tar, out, [g_path, sel_path, *ps])
        sizes[route] = tar.stat().st_size
    (work / "bundles" / "sizes.json").write_text(json.dumps(sizes, indent=2) + "\n")
    return {"selection": counts, "route_bytes": sizes, "tree": out, "commit": commit}


def release(top: Path, archive: Path, force: bool) -> tuple[str, int, int]:
    """Pack ``top`` as ``archive`` (left read-only); (sha256, size, file count). An existing
    archive with other content is kept unless ``force``."""
    members = [top, *top.rglob("*")]
    tmp = archive.with_name(archive.name + ".new")
    deterministic_tar(tmp, top, members)
    digest = sha256(tmp)
    if archive.exists() and not force and sha256(archive) != digest:
        tmp.unlink()
        raise BuildError(f"{archive} exists with another sha256; remove it or pass --force")
    os.replace(tmp, archive)
    archive.chmod(0o444)
    return digest, archive.stat().st_size, sum(1 for m in members if m.is_file())


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--start", type=int, default=START,
                    help=f"0-based region start on {CHROM} (default {START})")  # fmt: skip
    ap.add_argument("--end", type=int, default=END, help=f"exclusive region end (default {END})")
    ap.add_argument("--name", default=NAME, help=f"default {NAME}")
    ap.add_argument("--out", type=Path, default=None, help="work dir; default $GIPPA_ROOT/example")
    ap.add_argument("--bundle-out", type=Path, default=None,
                    help="where the archive goes; default $GIPPA_ROOT/example")  # fmt: skip
    ap.add_argument("--expected", type=Path, default=None,
                    help="copy DIR/* into the bundle's expected/ (none without it)")  # fmt: skip
    ap.add_argument("--expected-commit", default="", help="the commit that produced --expected")
    ap.add_argument("--force", action="store_true", help="replace an archive with another sha256")
    a = ap.parse_args(argv)
    if not 0 <= a.start < a.end:
        ap.error("need 0 <= start < end")
    if a.expected is not None and not a.expected.is_dir():
        ap.error(f"--expected {a.expected}: not a directory")
    import gippa

    root = _root()
    try:
        version = bundle_version(str(gippa.__version__))
        res = build(a.start, a.end, a.name, a.out or root / "example", root, version=version,
                    expected=a.expected, expected_commit=a.expected_commit)  # fmt: skip
        archive = (a.bundle_out or root / "example") / T.archive_name(version)
        digest, size, n_files = release(res["tree"], archive, a.force)
    except BuildError as exc:
        print(f"build_example: {exc}", file=sys.stderr)
        return 1
    job = os.environ.get("SLURM_JOB_ID", "")
    print(json.dumps({"selection": res["selection"], "route_bytes": res["route_bytes"]}, indent=2))
    print(f"\narchive: {archive}\nsha256: {digest}\nsize_bytes: {size}\nfiles: {n_files}")
    print("\n# db_sources.yaml row\n" + T.db_sources_row(version, digest, size, job))
    print("\n# validation/example/bundles.tsv row\n" + "\t".join(T.TSV_COLUMNS))
    print(T.bundles_tsv_row(version, digest, size, res["commit"], job,
                            a.expected_commit if a.expected is not None else ""))  # fmt: skip
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
