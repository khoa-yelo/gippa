"""Collect per-dataset provenance (recipe step 10) into one JSON: GIPPA commit and dirty flag,
container digests, input sha256s, job IDs and wall times of every step, DIA-NN per-run
precursor counts and protein groups, database and library sizes. Aggregates only.

    collect_provenance.py D OUT_JSON      (D = $GIPPA_ROOT/campaign/<accession>)
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path
from typing import Any


def kv(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    if path.is_file():
        for tok in path.read_text().split():
            if "=" in tok:
                k, v = tok.split("=", 1)
                out[k] = v
    return out


def diann(dn: Path) -> dict[str, Any]:
    p = dn / "provenance_report.json"
    d = json.loads(p.read_text()) if p.is_file() else {}
    log = (dn / "step2_report.log").read_text() if (dn / "step2_report.log").is_file() else ""
    s1 = (dn / "step1.log").read_text() if (dn / "step1.log").is_file() else ""
    prec = [int(x) for x in re.findall(r"Precursors at 1% peptidoform FDR: (\d+)", log)]
    pg = re.findall(r"Protein groups with global q-value <= 0.01: (\d+)", log)
    lib = re.findall(r"(\d+) precursors generated", s1)
    return {
        "job_id": d.get("job_id"),
        "exit_code": d.get("exit_code"),
        "step1": d.get("step1"),
        "step2": d.get("step2"),
        "library_bytes": d.get("library_bytes"),
        "library_precursors": int(lib[-1]) if lib else None,
        "image_sha256": d.get("image_sha256"),
        "cfg_sha256": d.get("cfg_sha256"),
        "fasta_sha256": d.get("fasta_sha256"),
        "n_runs": d.get("n_runs"),
        "precursors_1pct_per_run": prec,
        "protein_groups_global_q001": int(pg[-1]) if pg else None,
        "version": (dn / "diann_version.txt").read_text().strip()
        if (dn / "diann_version.txt").is_file()
        else None,
    }


def run(out: Path) -> dict[str, Any]:
    m = json.loads((out / "run_manifest.json").read_text())
    return {
        k: m.get(k)
        for k in (
            "git_commit",
            "git_dirty",
            "gippa_version",
            "exit_code",
            "status",
            "started_at",
            "finished_at",
            "profile",
            "modules_off",
            "warnings",
            "layers_missing",
        )
    } | {
        "master_seed": (m.get("seeds") or {}).get("master"),
        "inputs": {i.get("role"): i.get("sha256") for i in m.get("inputs", [])},
        "tools": {
            t["name"]: {"version": t.get("version"), "container_digest": t.get("container_digest")}
            for t in m.get("tools", [])
        },
    }


def main() -> int:
    D, dest = Path(sys.argv[1]), Path(sys.argv[2])
    doc: dict[str, Any] = {"dataset": D.name}
    for r in ("a", "b"):
        rt = D / f"runtime_{r}.txt"
        if not rt.is_file():
            continue
        cp = D / f"in_{r}" / "convert_provenance.json"
        db = D / f"db_{r}" / "strain_db.faa"
        doc[f"route_{r}"] = {
            "runtime": kv(rt),
            "makedb_fasta_bytes": db.stat().st_size if db.is_file() else None,
            "makedb_records": sum(1 for ln in db.open() if ln.startswith(">"))
            if db.is_file()
            else None,
            "diann": diann(D / f"diann_{r}"),
            "convert_preflight": json.loads(cp.read_text()) if cp.is_file() else None,
            "gippa_run": run(D / f"out_{r}")
            if (D / f"out_{r}" / "run_manifest.json").is_file()
            else None,
        }
    pr = D / "post_runtime.txt"
    doc["post"] = pr.read_text().splitlines() if pr.is_file() else None
    sums = D / "raw" / "SHA256SUMS"
    doc["raw_sha256"] = (
        dict(reversed(ln.split()) for ln in sums.read_text().splitlines())
        if sums.is_file()
        else None
    )
    dest.write_text(json.dumps(doc, indent=1, default=str) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
