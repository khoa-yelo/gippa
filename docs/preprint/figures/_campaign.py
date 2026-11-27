"""Readers for the demo and discovery campaign evidence (validation/campaign/).

Until the campaign branches merge, the evidence lives in the worktrees campaign-p2a (Phase 2
summary, group A, novel_pc/) and campaign-ecoli (Phase 1). Every number comes from a committed
file: PHASE2_SUMMARY.md tables, per-dataset metrics.json, ecoli/metrics_phase1_fix.json and
novel_pc/novel_pc_threshold.json / passers.tsv. Aggregates only.
"""

from __future__ import annotations

import csv
import json
import re
from pathlib import Path

from _common import REPO, WORKTREES


def campaign_dir(worktree: str) -> Path:
    """validation/campaign/ on main once merged, else in the named worktree."""
    main = REPO / "validation" / "campaign"
    if (main / "PHASE2_SUMMARY.md").is_file() and worktree == "campaign-p2a":
        return main
    if (main / "ecoli" / "metrics_phase1_fix.json").is_file() and worktree == "campaign-ecoli":
        return main
    return WORKTREES / worktree / "validation" / "campaign"


P2A = campaign_dir("campaign-p2a")
ECOLI = campaign_dir("campaign-ecoli")
SUMMARY = P2A / "PHASE2_SUMMARY.md"

#: summary row label -> (short label, run id used in novel_pc)
DATASETS = [
    ("E. coli", "E. coli K-12 (Ph. 1)", "ecoli_e2"),
    ("C. crescentus", "C. crescentus", "ccre"),
    ("B. theta", "B. thetaiotaomicron", "bthe"),
    ("Synechococcus", "Synechococcus", "syn"),
    ("V. natriegens", "V. natriegens", "vnat"),
    ("M. xanthus", "M. xanthus", "mxan"),
    ("C. difficile", "C. difficile", "cdif"),
    ("M. tuberculosis", "M. tuberculosis", "mtb"),
    ("H. volcanii", "H. volcanii", "hvo"),
    ("S. islandicus", "S. islandicus", "sisl"),
]
TPA = [("tpa_dal1", "DAL-1"), ("tpa_grady", "Grady"), ("tpa_haitib", "Haiti B"),
       ("tpa_madras", "Madras"), ("tpa_mexicoa", "Mexico A"), ("tpa_phil1", "Phil 1"),
       ("tpa_ss14", "SS14")]  # fmt: skip


def section(md: str, heading: str) -> list[list[str]]:
    """Rows (cells) of the first Markdown table after ``heading``."""
    start = md.index(heading)
    rows: list[list[str]] = []
    for line in md[start:].splitlines()[1:]:
        if line.startswith("|"):
            cells = [c.strip() for c in line.strip().strip("|").split("|")]
            if set("".join(cells)) <= set("-: "):
                continue
            rows.append(cells)
        elif rows:
            break
    return rows[1:]  # drop the header


def _num(s: str) -> int:
    return int(s.replace(",", ""))


def detection() -> list[dict]:
    """Route-A genes detected / coding genes per genome (§2.3, T. pallidum per strain)."""
    md = SUMMARY.read_text()
    out = []
    by_key = {r[0].strip("*"): r for r in section(md, "### 2.3 Detection depth")}
    for key, label, run in DATASETS:
        row = next(v for k, v in by_key.items() if k.startswith(key) or key in k)
        m = re.search(r"([\d,]+) / ([\d,]+) =", row[3])
        if not m:
            raise SystemExit(f"no k / n in §2.3 row {key}")
        out.append({"dataset": label, "run": run, "k": _num(m.group(1)), "n": _num(m.group(2)),
                    "source": f"{SUMMARY.name} §2.3"})  # fmt: skip
    tpa = json.loads((P2A / "PXD065022" / "metrics.json").read_text())["T2"]["per_strain"]
    for run, name in TPA:
        out.append(
            {
                "dataset": f"T. pallidum {name}",
                "run": run,
                "k": tpa[run]["k"],
                "n": tpa[run]["n"],
                "source": "PXD065022/metrics.json › T2.per_strain",
            }
        )
    return out


def retention() -> list[dict]:
    """Share of route-A genes that route B keeps (k / n)."""
    out = []
    e = json.loads((ECOLI / "ecoli" / "metrics_phase1_fix.json").read_text())
    kb = e["row14_route_b"]["E1_kept_by_E2"]
    out.append({"dataset": DATASETS[0][1], "run": "ecoli_e2", "k": kb["k"], "n": kb["n"],
                "source": "ecoli/metrics_phase1_fix.json › row14_route_b"})  # fmt: skip
    md = SUMMARY.read_text()
    rows = section(md, "### 2.4 Route B: gene retention and cost")
    for key, label, run in DATASETS[1:]:
        row = next(r for r in rows if r[0].startswith(key) or key in r[0])
        m = re.match(r"([\d.]+) % \[([\d.]+), ([\d.]+)\]", row[1])
        if not m:
            raise SystemExit(f"no kept share in §2.4 row {key}")
        v, lo, hi = (float(x) / 100 for x in m.groups())
        out.append({"dataset": label, "run": run, "share": v, "ci_low": lo, "ci_high": hi,
                    "source": f"{SUMMARY.name} §2.4"})  # fmt: skip
    tpa = json.loads((P2A / "PXD065022" / "metrics.json").read_text())["strains"]
    for run, name in TPA:
        kb = tpa[run]["metrics"]["S9_route_b"]["A_kept_by_B"]
        out.append(
            {
                "dataset": f"T. pallidum {name}",
                "run": run,
                "k": kb["k"],
                "n": kb["n"],
                "source": "PXD065022/metrics.json › strains.*.metrics.S9_route_b",
            }
        )
    return out


def te4() -> list[dict]:
    """log2(RplL/RplJ) iBAQ and top3, and log2(AtpA/AtpD), per genome (§2.5)."""
    md = SUMMARY.read_text()
    out = []
    for r in section(md, "### 2.5 RplL/RplJ stoichiometry (TE4)"):
        name = r[0].strip("*")
        m = re.match(r"\**(−?-?[\d.]+)\** \((−?-?[\d.]+)\)", r[1])
        if not m:
            continue  # the archaea row (NT)
        f = lambda s: float(s.replace("−", "-"))  # noqa: E731
        atp = r[5].split()[0].replace("−", "-")
        out.append({"genome": name.replace("*", ""), "rpl_ibaq": f(m.group(1)),
                    "rpl_top3": f(m.group(2)), "verdict": r[2].strip("*"),
                    "atp": float(atp) if re.match(r"-?[\d.]+$", atp) else None,
                    "source": f"{SUMMARY.name} §2.5"})  # fmt: skip
    return out


def verdicts() -> list[dict]:
    """Pass / fail / not testable or not run / reported counts per dataset (§1)."""
    md = SUMMARY.read_text()
    out = []
    start = md.index("| Dataset | TE1 |")
    lines = md[start:].splitlines()[2:]
    for line in lines:
        if not line.startswith("|"):
            break
        r = [c.strip() for c in line.strip().strip("|").split("|")]
        name = r[0].strip("*")
        m = re.match(r"\**(\d+) / (\d+) / (\d+)(?: NT, \d+ dropped)? / (\d+)", r[15])
        if not m or name.startswith("Total"):
            continue
        p, f, nt, rep = (int(x) for x in m.groups())
        out.append({"dataset": name.strip("*").replace("*", ""), "pass": p, "fail": f,
                    "nt_nr": nt, "reported": rep, "source": f"{SUMMARY.name} §1"})  # fmt: skip
    return out


def novel_pc() -> dict:
    return json.loads((P2A / "novel_pc" / "novel_pc_threshold.json").read_text())


def candidates() -> list[dict]:
    """Every candidate with a numeric class FDR: the tau_all tier of passers.tsv (q <= 1)."""
    with open(P2A / "novel_pc" / "passers.tsv") as fh:
        return [r for r in csv.DictReader(fh, delimiter="\t") if r["threshold"] == "tau_all"]
