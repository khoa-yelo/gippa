"""PXD065022 partial aggregates while gippa run is blocked (#611 fixed in #612; #614 open):
DIA-NN counts per strain and route, convert/preflight results, T1 host share (route-A peptide
tables), T6 truth-set counts, and the failed gippa runs' exit/step. Aggregates only.

    PYTHONPATH=<pin>:<pin>/src python tpa_partial.py OUT.json
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import p2_metrics as m
import provenance as pv

C = Path("/home/classes/bios/270/khoa/giap/local/campaign/PXD065022")
TAGS = ["tpa_dal1", "tpa_grady", "tpa_haitib", "tpa_madras", "tpa_mexicoa", "tpa_phil1", "tpa_ss14"]


def failed_run(out: Path) -> dict | None:
    man = out / "run_manifest.json"
    if not man.is_file():
        return None
    d = json.loads(man.read_text())
    return {
        "git_commit": d.get("git_commit"),
        "exit_code": d.get("exit_code"),
        "status": d.get("status"),
    }


def main() -> int:
    doc: dict = {"schema": "gippa-campaign-phase2a-tpa-partial/1", "strains": {}}
    for t in TAGS:
        s: dict = {}
        for r in ("a", "b"):
            s[f"diann_{r}"] = pv.diann(C / t / f"diann_{r}")
            cp = C / t / f"in_{r}" / "convert_provenance.json"
            s[f"convert_{r}"] = json.loads(cp.read_text()) if cp.is_file() else None
            s[f"run_{r}_bfc74ff6"] = failed_run(C / t / f"out_{r}")
            s[f"run_{r}_cff0fd68"] = failed_run(C / t / f"out_{r}_fix")
        pep = C / t / "in_a" / "peptides.tsv.gz"
        s["T1_host"] = m.host_share(pep) if pep.is_file() else None
        tc = C / "truth" / t / "counts.json"
        s["T6_truth_counts"] = json.loads(tc.read_text()) if tc.is_file() else None
        doc["strains"][t] = s
    Path(sys.argv[1]).write_text(json.dumps(doc, indent=1, default=str) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
