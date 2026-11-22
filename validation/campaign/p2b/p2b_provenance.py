"""Standard recipe step 10 for one Phase 2 group B dataset: write ``<ACC>/provenance.json``.

Aggregates only: commit, image digests, input sha256s, job IDs, wall times, database and library
sizes, DIA-NN per-run precursor and protein counts (``report.stats.tsv`` columns, no rows of the
report itself). Reads ``lane_steps.tsv`` (lane.sbatch), ``db_*/strain_db.faa``,
``diann_*/{step1.done,step1.log,provenance_report.json,report.stats.tsv}``,
``in_*/convert_provenance.json``, ``out_*/run_manifest.json``, ``raw/SHA256SUMS`` and ``sacct``.

    python3 -I p2b_provenance.py $GIPPA_ROOT/campaign/ACC validation/campaign/ACC/provenance.json
"""

import csv
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path


def sha(p: Path) -> str | None:
    if not p.is_file():
        return None
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def sacct(job: str) -> dict:
    out = subprocess.run(["sacct", "-j", job, "-X", "-n", "-P", "-o",
                          "JobID,JobName,State,Elapsed,ExitCode,NodeList,AllocCPUS,Start,End"],
                         capture_output=True, text=True).stdout.strip().splitlines()
    rss = subprocess.run(["sacct", "-j", job, "-n", "-P", "-o", "MaxRSS"], capture_output=True, text=True).stdout.split()
    if not out:
        return {"job_id": job, "state": "unknown"}
    k = ["job_id", "name", "state", "elapsed", "exit_code", "node", "alloc_cpus", "start", "end"]
    d = dict(zip(k, out[0].split("|")))
    d["max_rss"] = max(rss, key=lambda s: float(re.sub(r"[^0-9.]", "", s) or 0) * {"K": 1, "M": 1e3, "G": 1e6}.get(s[-1:], 1e-3), default=None)
    return d


def db(p: Path) -> dict:
    f = p / "strain_db.faa"
    if not f.is_file():
        return {}
    n = cont = 0
    with open(f) as fh:
        for ln in fh:
            if ln.startswith(">"):
                n += 1
                cont += ln.startswith(">CONT_")
    return {"records": n, "contaminant_records": cont, "bytes": f.stat().st_size, "sha256": sha(f)}


def diann(p: Path) -> dict:
    d: dict = {}
    s1 = p / "step1.log"
    if s1.is_file():
        m = re.search(r"(\d+) precursors generated", s1.read_text())
        d["library_precursors"] = int(m.group(1)) if m else None
    done = p / "step1.done"
    if done.is_file():
        m = re.search(r"job_id=(\S*) wall_s=(\d+)", done.read_text())
        d["step1_job"], d["step1_wall_s"] = (m.group(1), int(m.group(2))) if m else (None, None)
    lib = p / "lib.predicted.speclib"
    d["library_bytes"] = lib.stat().st_size if lib.is_file() else None
    pr = p / "provenance_report.json"
    if pr.is_file():
        j = json.loads(pr.read_text())
        d |= {"step2_job": j.get("job_id"), "step2_wall_s": j["step2"]["wall_s"], "n_runs": j.get("n_runs"),
              "image_sha256": j.get("image_sha256"), "cfg_sha256": j.get("cfg_sha256"),
              "fasta_sha256": j.get("fasta_sha256")}
    st = p / "report.stats.tsv"
    if st.is_file():
        rows = list(csv.DictReader(open(st), delimiter="\t"))
        keep = [c for c in ("Precursors.Identified", "Proteins.Identified") if rows and c in rows[0]]
        d["per_run"] = [{"run": Path(r.get("File.Name", "")).name, **{c: r[c] for c in keep}} for r in rows]
    rp = p / "report.parquet"
    d["report_sha256"] = sha(rp)
    return d


def run(p: Path) -> dict:
    f = p / "run_manifest.json"
    if not f.is_file():
        return {}
    m = json.loads(f.read_text())
    return {"git_commit": m.get("git_commit"), "git_dirty": m.get("git_dirty"), "gippa_version": m.get("gippa_version"),
            "profile": m.get("profile"), "exit_code": m.get("exit_code"), "started_at": m.get("started_at"),
            "finished_at": m.get("finished_at"), "seeds": m.get("seeds"), "modules_off": m.get("modules_off"),
            "layers_missing": m.get("layers_missing"), "warnings": m.get("warnings"),
            "inputs": {i.get("role"): i.get("sha256") for i in m.get("inputs", [])},
            "tools": {t["name"]: {"version": t.get("version"), "digest": t.get("container_digest")} for t in m.get("tools", [])}}


def main(d: str, out: str) -> int:
    D = Path(d)
    steps = [dict(zip(("step", "exit", "wall_s", "job"), ln.rstrip("\n").split("\t")))
             for ln in open(D / "lane_steps.tsv")] if (D / "lane_steps.tsv").is_file() else []
    jobs = sorted({s["job"] for s in steps if s.get("job")})
    raw = [ln.split() for ln in open(D / "raw" / "SHA256SUMS")] if (D / "raw" / "SHA256SUMS").is_file() else []
    doc = {"dataset": D.name, "lane_steps": steps, "jobs": {j: sacct(j) for j in jobs},
           "raw_sha256": {Path(x[-1]).name: x[0] for x in raw if len(x) >= 2},
           "db": {r: db(D / f"db_{r}") for r in ("a", "b")},
           "diann": {r: diann(D / f"diann_{r}") for r in ("a", "b")},
           "convert": {r: json.loads((D / f"in_{r}" / "convert_provenance.json").read_text())
                       for r in ("a", "b") if (D / f"in_{r}" / "convert_provenance.json").is_file()},
           "runs": {r: run(D / f"out_{r}") for r in ("a", "b", "a_te7", "r2") if (D / f"out_{r}").is_dir()}}
    Path(out).write_text(json.dumps(doc, indent=1) + "\n")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:3]))
