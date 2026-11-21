"""Phase 1b check 1b-5: TE8-style decoy placement FDR of N2's peptides on MG1655 at defaults.

Runs the pinned TE8 chain (validation/proxies/te_proxies.py: preflight + the local executor to
the ``placement`` step, default config, search_db other, genome_relation related) on the N2
evidence peptides (target, q <= input.q_value_max, not contaminant) and applies TE8's pooled
estimator ``pooled_decoy_fdr``. Writes counts only (JSON) to OUT.

    PYTHONPATH=<pin>/src:<pin>/validation/proxies python placement_fdr.py PEPTIDES GENOME JSON OUT
"""

from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path

import pandas as pd
import te_proxies as te
from gippa.config import resolve_config


def main(pep: Path, genome: Path, bakta_json: Path, out: Path) -> None:
    cfg = resolve_config()
    t = pd.read_csv(pep, sep="\t", dtype=str, keep_default_na=False)
    dec = t.get("is_decoy", pd.Series(["false"] * len(t))).str.lower().isin(["true", "1"])
    q = pd.to_numeric(t.get("q_value", pd.Series(["0"] * len(t))), errors="coerce").fillna(0.0)
    pids = t.get("protein_ids", pd.Series([""] * len(t)))
    cont = pids.str.split(";").apply(lambda xs: bool(xs) and all(
        any(x.startswith(p) for p in cfg.input.contaminant_prefixes) for x in xs if x))
    keep = t[(~dec) & (q <= cfg.input.q_value_max) & (~cont)]
    sources: dict[str, set[str]] = {}
    for p, ids in zip(keep["peptide"], keep["protein_ids"]):
        sources.setdefault(p, set()).update(i for i in ids.split(";") if i)
    with tempfile.TemporaryDirectory(prefix="gippa-1b5-") as tmp:
        w = te.run_chain(cfg, genome, bakta_json, sources, "other", Path(tmp), ["placement"])
        per = pd.read_parquet(w / "placement")
        costs = te._costs(w)
    acc = per[per["accepted"].astype(bool)]
    m = acc["m_target"].to_numpy()
    res = {
        "peptides": len(sources),
        "placed_rows": int(len(per)),
        "accepted": int(len(acc)),
        "accepted_exact": int((m == 0).sum()),
        "accepted_m1": int((m == 1).sum()),
        "accepted_m2": int((m == 2).sum()),
        "rejected": int(len(per) - len(acc)),
        "unplaced": int(len(sources) - len(per)),
        "pooled_decoy_fdr": te.pooled_decoy_fdr(costs, per),
        "threshold": te.TE8_FDR_MAX,
    }
    res["passed"] = res["pooled_decoy_fdr"] <= te.TE8_FDR_MAX
    out.write_text(json.dumps(res, indent=2) + "\n")
    print(json.dumps(res))


if __name__ == "__main__":
    main(*map(Path, sys.argv[1:5]))
