"""Per-arm provenance (standard recipe step 10) for Phase 2 group A; counts and digests only.

    python3 -I provenance.py ARM_DIR ROUTE OUT.json      # ARM_DIR = $GIPPA_ROOT/campaign/ACC/TAG

Reads diann_R/provenance_report.json, step1.log/step2_report.log (precursor and protein-group
counts only), in_R/convert_provenance.json, in_R/runtime_out_R.txt and out_R/run_manifest.json
(commit, tool images, inputs, warnings). ROUTE a | b | te7 | x_<REF>. No peptide or DIA-NN row.
"""

import json
import re
import sys
from pathlib import Path


def jl(p: Path):
    return json.loads(p.read_text()) if p.is_file() else None


def diann(d: Path) -> dict | None:
    prov = jl(d / "provenance_report.json")
    if prov is None:
        return None
    s1 = (d / "step1.log").read_text() if (d / "step1.log").is_file() else ""
    s2 = (d / "step2_report.log").read_text() if (d / "step2_report.log").is_file() else ""
    m = re.search(r"(\d+) precursors generated", s1)
    pg = re.search(r"Protein groups with global q-value <= 0.01: (\d+)", s2)
    return {k: prov.get(k) for k in ("job_id", "node", "exit_code", "image_sha256", "cfg_sha256",
                                      "campaign_commit", "fasta_sha256", "n_runs", "step1", "step2",
                                      "wall_s", "library_bytes")} | {
        "library_precursors": int(m.group(1)) if m else None,
        "ids_at_1pct_per_run": [int(x) for x in re.findall(r"Number of IDs at 0.01 FDR: (\d+)", s2)],
        "protein_groups_global_q01": int(pg.group(1)) if pg else None,
        "runs": [Path(r).name for r in (prov.get("run_list") and Path(prov["run_list"]).read_text().split() or [])]}


def run_info(run: Path) -> dict | None:
    m = jl(run / "run_manifest.json")
    if m is None:
        return None
    return {k: m.get(k) for k in ("git_commit", "git_dirty", "gippa_version", "status", "exit_code",
                                  "started_at", "finished_at", "profile", "detection_tier",
                                  "modules_off", "warnings", "layers_missing", "tools", "databases")} | {
        "seeds_master": (m.get("seeds") or {}).get("master"),
        "inputs": [{k: i.get(k) for k in ("role", "path", "sha256")} for i in m.get("inputs", [])]}


def main() -> int:
    arm, route, out = Path(sys.argv[1]), sys.argv[2], Path(sys.argv[3])
    base = "a" if route in ("a", "te7") or route.startswith("x_") else route
    in_dir = arm / ("in_" + route if route.startswith("x_") else "in_" + base)
    out_dir = arm / ("out_a_te7" if route == "te7" else "out_" + route)
    rt = in_dir / f"runtime_{out_dir.name}.txt"
    doc = {"arm_dir": str(arm), "route": route, "diann": diann(arm / f"diann_{base}"),
           "convert": jl(in_dir / "convert_provenance.json"),
           "gippa_run_runtime": rt.read_text().strip() if rt.is_file() else None,
           "gippa_run": run_info(out_dir)}
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(doc, indent=1) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
