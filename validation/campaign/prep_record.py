# ruff: noqa: E501  (report-line f-strings; a validation script, not package code)
"""Campaign Phase 2 prep record: aggregates only (no peptide data).

For each dataset in genomes.tsv, read the fetch outputs ($GIPPA_ROOT/campaign/<acc>/raw/SHA256SUMS),
the genome downloads (genome/<label>/genome.sha256, md5checksums.txt) and the Bakta outputs
(bakta/<label>/out/<label>.txt, runtime.txt, versions.txt, command.txt), then
  * add a `sha256` column to validation/campaign/<acc>/files.tsv for the fetched files, and
  * write validation/campaign/<acc>/PREP.md.

    python -I validation/campaign/prep_record.py --root $GIPPA_ROOT/campaign \
        --fetch-job 2001986 --bakta-job 2001987 [--note ACC=text ...]
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
COUNT_KEYS = [
    "CDSs",
    "sORFs",
    "pseudogenes",
    "hypotheticals",
    "tRNAs",
    "rRNAs",
    "ncRNAs",
    "CRISPR arrays",
]


def read_tsv(p: Path) -> list[dict[str, str]]:
    with p.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def sha_map(p: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    if p.exists():
        for line in p.read_text().splitlines():
            if line.strip():
                h, name = line.split(None, 1)
                out[name.strip().lstrip("*")] = h
    return out


def bakta_counts(txt: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    if not txt.exists():
        return out
    for line in txt.read_text().splitlines():
        m = re.match(r"^([A-Za-z][A-Za-z ]*): (.+)$", line.strip())
        if m:
            out.setdefault(m.group(1), m.group(2))
    return out


def refseq_counts(gff: Path) -> tuple[int, int]:
    cds_ids: set[str] = set()
    pseudo = 0
    with gzip.open(gff, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            attrs = f[8]
            if f[2] == "CDS" and "pseudo=true" not in attrs:
                m = re.search(r"(?:^|;)ID=([^;]+)", attrs)
                cds_ids.add(m.group(1) if m else f"{f[0]}:{f[3]}:{f[4]}")
            elif f[2] == "pseudogene":
                pseudo += 1
    return len(cds_ids), pseudo


def kv(p: Path) -> str:
    return p.read_text().strip() if p.exists() else "missing"


def update_files_tsv(acc: str, shas: dict[str, str]) -> None:
    p = HERE / acc / "files.tsv"
    rows = read_tsv(p)
    fields = list(rows[0].keys())
    if "sha256" not in fields:
        fields.append("sha256")
    for r in rows:
        name = Path(r["file_name"]).name
        r["sha256"] = shas.get(name, r.get("sha256") or "NA")
    with p.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, required=True)
    ap.add_argument("--fetch-job", required=True)
    ap.add_argument("--bakta-job", required=True)
    ap.add_argument("--note", action="append", default=[], help="ACC=text, appended under Notes")
    a = ap.parse_args()
    notes: dict[str, list[str]] = {}
    for n in a.note:
        k, v = n.split("=", 1)
        notes.setdefault(k, []).append(v)

    genomes = read_tsv(HERE / "genomes.tsv")
    tasks = read_tsv(HERE / "fetch_tasks.tsv")
    sel = read_tsv(HERE / "run_selection.tsv")
    for acc in dict.fromkeys(g["dataset"] for g in genomes):
        man = {r["file_name"]: r for r in read_tsv(HERE / acc / "files.tsv")}
        shas = sha_map(a.root / acc / "raw" / "SHA256SUMS")
        update_files_tsv(acc, shas)
        mytasks = [t for t in tasks if t["accession"] == acc]
        wanted = [f for t in mytasks for f in t["files"].split()]
        label = next(
            (r["label"] for r in sel if r["accession"] == acc and r["file_name"] in wanted), ""
        )
        L = [f"# {acc}: Phase 2 preparation record", ""]
        L += [
            "Aggregates only. Data under `$GIPPA_ROOT/campaign/"
            + acc
            + "/` (not committed). No GIPPA run yet.",
            "",
        ]
        L += ["## Raw files", ""]
        L += [
            f"Selection: {label} (`run_selection.tsv`). Source: PRIDE over HTTPS (ADR-0077), resumable "
            f"transfer, byte size checked against the PRIDE file list, deposit SHA-1 checked where the "
            f"deposit has `checksum.txt`, sha256 recorded (also in `files.tsv`, column `sha256`). "
            f"Fetch job {a.fetch_job}, array task(s) {', '.join(t['task'] for t in mytasks)}.",
            "",
        ]
        L += ["| File | Bytes | Deposit checksum | Verified | sha256 |", "|---|---|---|---|---|"]
        got_bytes = 0
        for f in wanted:
            m = man[f]
            name = Path(f).name
            ok = name in shas
            raw = a.root / acc / "raw" / name
            if ok:
                got_bytes += raw.stat().st_size if raw.exists() else 0
            ck = m["checksum"] if m["checksum"] != "NA" else "none (no deposit checksum)"
            ver = (
                ("size + SHA-1" if m["checksum"].startswith("sha1:") else "size")
                if ok
                else "NOT FETCHED"
            )
            L.append(f"| {name} | {m['bytes']} | {ck} | {ver} | {shas.get(name, 'NA')} |")
        L += [
            "",
            f"Fetched {sum(Path(f).name in shas for f in wanted)}/{len(wanted)} files, "
            f"{got_bytes / 1e9:.1f} GB.",
            "",
        ]
        L += ["## Genome and Bakta annotation", ""]
        L += [
            "NCBI RefSeq assembly over HTTPS, MD5 checked against `md5checksums.txt`; replicons NCBI "
            "nuccore marks circular get `[topology=circular]`; Bakta 1.12.2, DB 6.0 full, pinned image "
            "and `gippa.annotate.bakta_cmd.bakta_command` argv (`--keep-contig-headers`) at the "
            f"campaign pin, as Phase 1b job 1998570. Bakta job {a.bakta_job} (array, one task per genome).",
            "",
        ]
        for g in [g for g in genomes if g["dataset"] == acc]:
            lab = g["label"]
            gd = a.root / acc / "genome" / lab
            bd = a.root / acc / "bakta" / lab
            asm = Path(g["ftp_path"]).name
            md5 = {}
            if (gd / "md5checksums.txt").exists():
                for line in (gd / "md5checksums.txt").read_text().splitlines():
                    h, n = line.split(None, 1)
                    md5[n.strip().lstrip("./")] = h
            gsha = sha_map(gd / "genome.sha256")
            c = bakta_counts(bd / "out" / f"{lab}.txt")
            L += [
                f"### {lab}: {g['assembly']} ({g['assembly_name']}, {g['assembly_level']}, strain {g['strain']})",
                "",
            ]
            L.append(
                f"- Replicons (accession:bp:NCBI topology): {g['replicons']}; tagged circular: {g['circular']}"
            )
            for suf in ("_genomic.fna.gz", "_genomic.gff.gz"):
                n = asm + suf
                L.append(
                    f"- `{n}`: MD5 {md5.get(n, 'missing')} (= md5checksums.txt, checked in job), sha256 {gsha.get(n, 'missing')}"
                )
            L.append(
                f"- Tagged genome `{lab}.fna.gz`: sha256 {gsha.get(lab + '.fna.gz', 'missing')}"
            )
            L.append(f"- Bakta runtime: {kv(bd / 'runtime.txt')}")
            vers = kv(bd / "versions.txt").splitlines()
            L.append(f"- Versions: {'; '.join(v for v in vers[:2])}")
            if c:
                L.append(
                    f"- Genome: {c.get('Length', '?')} bp, {c.get('Count', '?')} replicon(s), GC {c.get('GC', '?')} %, coding density {c.get('coding density', '?')} %"
                )
                L.append("- Bakta counts: " + ", ".join(f"{k} {c.get(k, '?')}" for k in COUNT_KEYS))
            else:
                L.append("- Bakta counts: NOT AVAILABLE (job not finished or failed)")
            gff = gd / (asm + "_genomic.gff.gz")
            if gff.exists():
                rc, rp = refseq_counts(gff)
                L.append(
                    f"- RefSeq annotation (answer key, not an input): {rc} CDS, {rp} pseudogenes"
                    + (
                        f"; Bakta CDS + sORF / RefSeq CDS = {(int(c['CDSs']) + int(c.get('sORFs', 0))) / rc:.3f}"
                        if c.get("CDSs")
                        else ""
                    )
                )
            L.append("")
        if acc in notes:
            L += ["## Notes", ""] + [f"- {n}" for n in notes[acc]] + [""]
        (HERE / acc / "PREP.md").write_text("\n".join(L))
        print(acc, "written")


if __name__ == "__main__":
    main()
