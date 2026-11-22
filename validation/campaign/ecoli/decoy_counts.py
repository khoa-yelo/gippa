"""Decoy competition counts of a finished ``gippa run`` (ADR-0078 recheck; aggregates only).

Finds the run's GENE_TABLE Nextflow task under ``OUT/work``, resolves the inputs it staged
(``ln -sfn`` lines of ``.command.run``), and recomputes the gene detection with the code on
``PYTHONPATH`` (the pinned tree): ``accepted_peptides`` -> ``gene_evidence`` -> ``detect`` and
``pair_decoys``. Prints a JSON of counts: decoy rows at q <= input.q_value_max, gene-paired vs
own-unit rows, own units, decoy wins (in genes, in genes without an accepted target, own units),
target wins, recomputed detections / max q vs ``tables/genes.tsv``. No peptide is written.

    PYTHONPATH=<pin>/src python -I decoy_counts.py OUT [OUT ...] > counts.json
"""

from __future__ import annotations

import json
import re
import shlex
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from gippa.config import resolve_config
from gippa.genetable.covariates import covariate_table
from gippa.genetable.detection import (
    detect,
    pair_decoys,
    protein_gene_map,
    shared_key_columns,
    _key_rank,
)
from gippa.genetable.evidence import accepted_peptides, gene_evidence
from gippa.io.work_tables import load_m4_inputs
from gippa.steps.gene_table import m4_paths


def task_inputs(out: Path) -> dict[str, Path]:
    cmds = [p for p in (out / "work").glob("*/*/.command.sh") if "_step gene_table" in p.read_text()]
    if len(cmds) != 1:
        raise SystemExit(f"{out}: {len(cmds)} gene_table tasks")
    sh = shlex.split(cmds[0].read_text().splitlines()[1])
    rel: dict[str, str] = {}
    for i, tok in enumerate(sh):
        if tok == "--in":
            k, v = sh[i + 1].split("=", 1)
            rel[k] = v
    links = dict(
        (m.group(2), m.group(1))
        for m in re.finditer(r"ln -sfn (\S+) (\S+)", (cmds[0].parent / ".command.run").read_text())
    )
    return {k: Path(links[v]) for k, v in rel.items()}


def counts(out: Path) -> dict:
    cfg = resolve_config(out / "resolved_config.yaml")
    inputs = task_inputs(out)
    m4, _ = load_m4_inputs(m4_paths(inputs))
    cov = covariate_table(m4.genes, cfg)
    acc = accepted_peptides(m4.peptide_hits, m4.peptides, cfg)
    ev = gene_evidence(m4.genes, m4.peptide_hits, acc, m4.replicons)
    det = detect(m4.genes, ev, acc, m4.peptides, m4.peptide_hits, m4.protein_hits, cov, cfg)
    pmap = protein_gene_map(m4.genes, m4.peptides, m4.peptide_hits, m4.protein_hits, accepted=acc,
                            contaminant_prefixes=cfg.input.contaminant_prefixes)
    pairing = pair_decoys(m4.peptides, pmap, cfg.input.decoy_prefixes,
                          q_value_max=cfg.input.q_value_max,
                          contaminant_prefixes=cfg.input.contaminant_prefixes)
    paired = pairing["gene"].notna().to_numpy(dtype=bool)
    # Win counting as _picked_q: best key rank per gene; own units are decoy wins by construction
    key = shared_key_columns(acc, pairing)
    n_t = len(acc)
    cols = {c: np.r_[acc[c].astype(float), pairing[c].astype(float)] for c in key}
    rank = _key_rank(cols, n_t + len(pairing))
    tags = [t for t, p in zip(m4.genes["locus_tag"].astype(str), m4.genes["partial"].astype(bool)) if not p]
    t_best = pd.Series(rank[:n_t]).groupby(np.asarray(acc["gene"].astype(str).tolist(), dtype=object)).max().reindex(tags)
    dr = pd.Series(rank[n_t:])
    g_best = dr[paired].groupby(np.asarray(pairing["gene"][paired].astype(str).tolist(), dtype=object)).max().reindex(tags)
    # as gippa.stats.tda.picked_fdr: the decoy wins ties (d >= t) and genes without a target
    t, d = t_best.to_numpy(float), g_best.to_numpy(float)
    dwin = np.isfinite(d) & (~np.isfinite(t) | (d >= np.where(np.isfinite(t), t, -np.inf)))
    twin = np.isfinite(t) & ~dwin
    own_units = pairing.loc[~paired, "unit"].nunique()
    tab = det.table
    g = pd.read_csv(out / "tables" / "genes.tsv", sep="\t", comment="#", low_memory=False)
    gc = g[g["feature_type"].astype(str).isin(["cds", "sorf"])]
    gq = pd.to_numeric(gc.loc[gc["detection_status"] == "detected", "gene_q_value"], errors="coerce")
    sm = det.summary_rows
    srow = {r.metric: r.value for r in sm.itertuples() if r.group_type == "tier"}
    return {
        "decoy_rows": int(len(pairing)), "gene_paired_rows": int(paired.sum()),
        "own_unit_rows": int((~paired).sum()), "own_units": int(own_units),
        "gene_unit_decoy_wins": int(dwin.sum()),
        "gene_unit_decoy_wins_no_accepted_target": int((dwin & ~np.isfinite(t)).sum()),
        "decoy_wins_total": int(dwin.sum() + own_units),
        "target_wins": int(twin.sum()),
        "recomputed_detected": int(tab["detected"].sum()),
        "recomputed_max_q": float(tab.loc[tab["detected"], "gene_q_value"].max()),
        "table_detected": int(len(gq)), "table_max_q": float(gq.max()),
        "decoy_paired_share": float(srow.get("decoy_paired_share", np.nan)),
        "decoy_win_weight": float(srow.get("decoy_win_weight", np.nan)),
        "tier": det.tier,
    }


if __name__ == "__main__":
    print(json.dumps({a: counts(Path(a)) for a in sys.argv[1:]}, indent=1))
