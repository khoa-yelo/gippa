"""Write ``$GIPPA_ROOT/campaign/ecoli_k12/<arm>/provenance.json`` for each Phase 1 arm.

Standard recipe step 10: commit, image digests, input sha256s, job IDs and wall times. Reads the
DIA-NN ``provenance_<stem>.json`` files of ``diann_search.sbatch``, ``convert_provenance.json``
of ``convert_preflight.sbatch``, the GIPPA ``run_manifest.json`` (commit, tool images, inputs)
and ``sacct`` (state, elapsed, MaxRSS). No peptide or DIA-NN row is read.

    python3 -I collect_provenance.py SPEC.json      # SPEC: {"arm": {"jobs": {...}, ...}, ...}

SPEC per arm: ``dir`` (provenance.json goes here), ``jobs`` (label -> Slurm job ID),
``files`` (label -> path, sha256 recorded), ``diann`` (paths of provenance_<stem>.json),
``convert`` (path of convert_provenance.json), ``run`` (GIPPA run directory).
"""

import hashlib
import json
import subprocess
import sys
from pathlib import Path


def sha(p: str) -> str | None:
    if not Path(p).is_file():
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
    rss = subprocess.run(["sacct", "-j", job, "-n", "-P", "-o", "MaxRSS"],
                         capture_output=True, text=True).stdout.split()
    if not out:
        return {"job_id": job, "state": "unknown"}
    k = ["job_id", "name", "state", "elapsed", "exit_code", "node", "alloc_cpus", "start", "end"]
    d = dict(zip(k, out[0].split("|")))
    d["max_rss"] = max(rss, key=lambda s: float(s.rstrip("KMG") or 0) * {"K": 1, "M": 1e3, "G": 1e6}.get(s[-1:], 1e-3), default=None)
    return d


def run_info(run: Path) -> dict:
    m = json.loads((run / "run_manifest.json").read_text())
    tools = m.get("tools") or {}
    return {"git_commit": m.get("git_commit"), "git_dirty": m.get("git_dirty"),
            "gippa_version": m.get("gippa_version") or m.get("giap_version"),
            "status": m.get("status"), "exit_code": m.get("exit_code"),
            "started_at": m.get("started_at"), "finished_at": m.get("finished_at"),
            "profile": m.get("profile"), "seeds_master": (m.get("seeds") or {}).get("master"),
            "detection_tier": m.get("detection_tier"), "modules_off": m.get("modules_off"),
            "inputs": [{k: i.get(k) for k in ("role", "path", "sha256")} for i in m.get("inputs", [])],
            "tools": tools, "databases": m.get("databases"), "warnings": m.get("warnings"),
            "layers_missing": m.get("layers_missing")}


def main() -> int:
    spec = json.loads(Path(sys.argv[1]).read_text())
    for arm, s in spec.items():
        d = {"arm": arm, "dataset": "PXD058422 / PXD066241 (PMC12529367)",
             "genome": "GCF_000005845.2 = U00096.3 by sequence",
             "jobs": {lab: sacct(j) for lab, j in s.get("jobs", {}).items()},
             "files_sha256": {lab: {"path": p, "sha256": sha(p)} for lab, p in s.get("files", {}).items()},
             "diann": [json.loads(Path(p).read_text()) for p in s.get("diann", []) if Path(p).is_file()],
             "convert": json.loads(Path(s["convert"]).read_text()) if s.get("convert") and Path(s["convert"]).is_file() else None,
             "gippa_run": run_info(Path(s["run"])) if s.get("run") and (Path(s["run"]) / "run_manifest.json").is_file() else None,
             "notes": s.get("notes")}
        out = Path(s["dir"]) / "provenance.json"
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(d, indent=1) + "\n")
        print(f"wrote {out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
