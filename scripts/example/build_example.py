"""Build the quick-start example bundle: one Nissle 1917 chromosome region, four input routes.

The example is a cut-down subset of the M10 pilot (validation/pilot/README.md, pass 3): E. coli
Nissle 1917 (GCF_021559835.1, chromosome NZ_CP058217.1), PXD078459 runs WT_1-WT_3. Nothing is
re-searched. The script subsets the pilot's existing outputs to one contiguous region and writes,
under OUT/<name>/:

    genome/<name>.fna.gz          the region as one linear contig (no topology tag)
    A_diann/report.parquet        DIA-NN 2.3.1 report rows whose proteins all lie in the region
    A_diann/search_db.faa         the region's records of the search database (`--proteins`): they
                                  map the report's protein IDs to genes, so the decoys pair
    B_proteins/proteins.faa       the region proteins DIA-NN quantified (A2 strain database)
    B_proteins/abundance.tsv      their PG.MaxLFQ (median over runs), as `gippa convert` computes it
    C_cascadia/WT_{1,2,3}.ssl     Cascadia rows whose peptide occurs in the region's six frames
    D_generic/peptides.tsv        a generic peptide table, peptide level, built from route A's rows
    SELECTION.json                the rule, the counts and the sha256 of every source and output
    bundles/<route>.tar.gz        genome + one route's files (+ SELECTION.json), for sizing

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

Usage (venv python, PYTHONPATH=<checkout>/src; small, runs in seconds):

    python scripts/example/build_example.py --start 1965000 --end 2065000 --name nissle_r100k \
        --out $GIPPA_ROOT/example
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import os
import re
import tarfile
from collections.abc import Iterable
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

CHROM = "NZ_CP058217.1"
RUNS = ("WT_1", "WT_2", "WT_3")


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


def deterministic_tar(out: Path, base: Path, members: list[Path]) -> None:
    with (
        out.open("wb") as raw,
        gzip.GzipFile(fileobj=raw, mode="wb", mtime=0, filename="") as gz,
        tarfile.open(fileobj=gz, mode="w", format=tarfile.PAX_FORMAT) as tar,
    ):
        for m in sorted(members):
            info = tar.gettarinfo(str(m), arcname=str(m.relative_to(base.parent)))
            info.mtime, info.uid, info.gid, info.uname, info.gname = 0, 0, 0, "", ""
            info.mode = 0o644
            with m.open("rb") as fh:
                tar.addfile(info, fh)


def build(start: int, end: int, name: str, out_root: Path, root: Path) -> dict[str, object]:
    src = sources(root)
    out = out_root / name
    for sub in ("genome", "A_diann", "B_proteins", "C_cascadia", "D_generic", "bundles"):
        (out / sub).mkdir(parents=True, exist_ok=True)

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
    a_path = out / "A_diann" / "report.parquet"
    pq.write_table(sub_t, a_path, compression="zstd")
    sub = sub_t.to_pandas()
    db = read_fasta(src["strain_db"])
    s_path = out / "A_diann" / "search_db.faa"
    write_fasta(s_path, [(pid, db[pid]) for pid in prot["db_id"]])

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

    targets = sub[sub["Decoy"] == 0]
    sel: dict[str, object] = {
        "name": name,
        "region": {
            "contig": CHROM,
            "start_0based": start,
            "end_exclusive": end,
            "length_bp": end - start,
            "one_based_closed": f"{CHROM}:{start + 1}-{end}",
            "output_contig_id": contig,
            "topology": "linear (no tag)",
        },
        "rule": __doc__.split("Selection rule (recorded verbatim in SELECTION.json):")[1]
        .split("Usage")[0]
        .strip(),
        "fdr_note": (
            "q-values, PEPs and MaxLFQ in routes A, B and D are DIA-NN's, computed on the "
            "whole-proteome A2 search (4,685 strain proteins + contaminants, three runs); Cascadia "
            "scores are per spectrum. Nothing was re-searched or re-scored for the subset."
        ),
        "counts": {
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
                "rows_per_run": {
                    k: int(v) for k, v in sub["Run"].value_counts().sort_index().items()
                },
            },
            "B": {"abundance_rows": len(ab), "fasta_records": len(ab)},
            "C": c_counts,
            "D": {
                "rows": len(gen),
                "decoy_rows": int((gen["is_decoy"] == "true").sum()),
            },
        },
        "sources": {k: {"path": str(v), "sha256": sha256(v)} for k, v in src.items()},
    }
    routes = {
        "A_diann": [a_path, s_path],
        "B_proteins": [f_path, ab_path],
        "C_cascadia": [out / "C_cascadia" / f"{r}.ssl" for r in RUNS],
        "D_generic": [d_path],
    }
    files = [g_path, *(p for ps in routes.values() for p in ps)]
    sel["outputs"] = {
        str(p.relative_to(out)): {"bytes": p.stat().st_size, "sha256": sha256(p)} for p in files
    }
    sel_path = out / "SELECTION.json"
    sel_path.write_text(json.dumps(sel, indent=2) + "\n")
    sizes: dict[str, int] = {}
    for route, ps in routes.items():
        tar = out / "bundles" / f"{name}_{route}.tar.gz"
        deterministic_tar(tar, out, [g_path, sel_path, *ps])
        sizes[route] = tar.stat().st_size
    allp = out / "bundles" / f"{name}_all.tar.gz"
    deterministic_tar(allp, out, [sel_path, *files])
    sizes["all"] = allp.stat().st_size
    (out / "bundles" / "sizes.json").write_text(json.dumps(sizes, indent=2) + "\n")
    return {"selection": sel["counts"], "bundle_bytes": sizes}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument(
        "--start", type=int, required=True, help="0-based region start on NZ_CP058217.1"
    )
    ap.add_argument("--end", type=int, required=True, help="exclusive region end")
    ap.add_argument("--name", required=True)
    ap.add_argument("--out", type=Path, default=None, help="default $GIPPA_ROOT/example")
    a = ap.parse_args(argv)
    if not 0 <= a.start < a.end:
        ap.error("need 0 <= start < end")
    root = _root()
    res = build(a.start, a.end, a.name, a.out or root / "example", root)
    print(json.dumps(res, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
