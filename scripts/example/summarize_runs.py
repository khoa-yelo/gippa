"""Summarise the example timing runs: wall time, per-process times and result sanity.

Reads every run directory under $GIPPA_ROOT/example/runs/ (scripts/slurm/example_run.sbatch) and
prints one JSON object per run: timing.txt, the Nextflow trace (realtime per process, summed
over tasks, the top processes and the critical-path-ish span from first submit to last complete),
and from the outputs: genes, detected genes, the detection tier and gene-level FDR rows, the counts
of the tables that small linear contigs make degenerate (operons, operon flags, islands, regions,
autocorrelation, enrichment) and the skipped-plot list of the manifest.

Usage: python scripts/example/summarize_runs.py [RUN_DIR ...]
"""

from __future__ import annotations

import json
import os
import re
import sys
from pathlib import Path

import pandas as pd

_UNITS = {"ms": 0.001, "s": 1.0, "m": 60.0, "h": 3600.0, "d": 86400.0}


def duration_s(text: str) -> float:
    if not isinstance(text, str) or text in ("-", ""):
        return 0.0
    return sum(float(v) * _UNITS[u] for v, u in re.findall(r"([\d.]+)(ms|s|m|h|d)", text))


def _rows(path: Path) -> int | None:
    if not path.exists():
        return None
    try:
        return len(pd.read_csv(path, sep="\t", comment="#", low_memory=False))
    except pd.errors.EmptyDataError:
        return 0


def trace_summary(trace: Path, top: int = 8) -> dict[str, object]:
    t = pd.read_csv(trace, sep="\t")
    t["proc"] = t["name"].str.replace(r" \(.*\)$", "", regex=True)
    t["real_s"] = t["realtime"].map(duration_s)
    t["dur_s"] = t["duration"].map(duration_s)
    t["submit_dt"] = pd.to_datetime(t["submit"], errors="coerce")
    end = (t["submit_dt"] + pd.to_timedelta(t["dur_s"], unit="s")).max()
    by = t.groupby("proc").agg(
        tasks=("real_s", "size"), real_s=("real_s", "sum"), dur_s=("dur_s", "max")
    )
    by = by.sort_values("real_s", ascending=False)
    return {
        "tasks": len(t),
        "failed": int((t["status"] == "FAILED").sum()),
        "nf_span_s": round((end - t["submit_dt"].min()).total_seconds(), 1),
        "sum_realtime_s": round(float(t["real_s"].sum()), 1),
        "top": [
            {
                "process": p,
                "tasks": int(r.tasks),
                "realtime_s": round(r.real_s, 1),
                "max_task_duration_s": round(r.dur_s, 1),
            }
            for p, r in by.head(top).iterrows()
        ],
    }


def results(out: Path) -> dict[str, object]:
    res: dict[str, object] = {}
    genes = out / "tables" / "genes.tsv"
    if genes.exists():
        g = pd.read_csv(genes, sep="\t", comment="#", low_memory=False)
        res["genes"] = len(g)
        res["genes_partial"] = int(g["partial"].astype(str).str.lower().eq("true").sum())
        det = g["detected"].astype(str).str.lower().eq("true")
        res["detected"] = int(det.sum())
        res["detection_tier"] = sorted(g["detection_tier"].dropna().astype(str).unique().tolist())
        res["tu_ge2"] = int(g.loc[g["tu_size"].fillna(0) >= 2, "tu_id"].nunique())
    ds = out / "tables" / "detection_summary.tsv"
    if ds.exists():
        d = pd.read_csv(ds, sep="\t", comment="#")
        tier = d[d["section"] == "tier"]
        res["tier_rows"] = {
            f"{r.group}:{r.metric}": r.value for r in tier.itertuples() if r.group_type == "tier"
        }
        q01 = tier[(tier["group"].astype(str) == "0.01")]
        res["q01"] = {r.metric: r.value for r in q01.itertuples()}
    for name in (
        "operons",
        "operon_flags",
        "islands",
        "regions_stats",
        "autocorrelation",
        "enrichment_detection",
        "enrichment_gsea",
        "variants",
        "novel_orfs",
        "determinants_model",
        "sensitivity",
        "category_summary",
    ):
        res[f"rows_{name}"] = _rows(out / "tables" / f"{name}.tsv")
    rp = out / "regions.parquet"
    if rp.exists():
        r = pd.read_parquet(rp)
        res["regions"] = r["region_type"].value_counts().to_dict()
    man = out / "run_manifest.json"
    if man.exists():
        m = json.loads(man.read_text())
        for k in ("status", "modules_off", "layers_missing", "warnings", "flags"):
            if k in m:
                v = m[k]
                res[f"manifest_{k}"] = v if not isinstance(v, list) or len(v) < 40 else v[:40]
    return res


def summarise(run: Path) -> dict[str, object]:
    s: dict[str, object] = {"run": run.name}
    tt = run / "timing.txt"
    if tt.exists():
        s.update(dict(kv.split("=", 1) for kv in tt.read_text().split() if "=" in kv))
    trace = run / "out" / "logs" / "trace.tsv"
    if trace.exists():
        s["trace"] = trace_summary(trace)
    s["results"] = results(run / "out")
    return s


def main(argv: list[str]) -> int:
    root = Path(os.environ.get("GIPPA_ROOT", "/home/classes/bios/270/khoa/giap/local"))
    runs = [Path(a) for a in argv] or sorted((root / "example" / "runs").glob("*"))
    for r in runs:
        print(json.dumps(summarise(r), default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
