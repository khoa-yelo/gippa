# ruff: noqa: E741, SIM115  (copied as run for the recheck)
"""Recompute decoy_picked gene detections (aggregates only) with whatever gippa is on PYTHONPATH.
Adapted from the statistician's n12/real.py. Usage: recheck.py LABEL RUN_DIR [RUN_DIR...]"""

import inspect
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import gippa
from gippa.config import resolve_config
from gippa.genetable import detection as det
from gippa.genetable.evidence import accepted_peptides
from gippa.io.work_tables import load_m4_inputs
from gippa.stats.tda import picked_fdr

label = sys.argv[1]
NEW = "f" not in inspect.signature(det._picked_q).parameters
res = []
for run in sys.argv[2:]:
    o = Path(run)
    I = o / "intermediate" if (o / "intermediate").is_dir() else o
    man = json.load(open(o / "run_manifest.json"))
    cfg = resolve_config(base=man["config"])
    paths = {
        k: I / f
        for k, f in [
            ("replicons", "replicons.parquet"),
            ("genes", "genes.parquet"),
            ("xrefs", "xrefs.parquet"),
            ("annotation_meta", "annotation_meta.json"),
            ("tus", "tus.parquet"),
            ("layers", "layers.parquet"),
        ]
    }
    paths["peptides"] = o / "work/preflight/peptides.norm.parquet"
    for k in ("peptide_hits", "protein_hits", "unmapped"):
        paths[k] = o / "mapping" / f"{k}.tsv"
    ph = pd.read_csv(paths["peptide_hits"], sep="\t", nrows=0)
    if "peptide_len" in ph.columns and "peptide_len_aa" not in ph.columns:  # 1.1 rename only
        tmp = (
            Path(__file__).parent / "renamed" / str(o.relative_to(o.parents[2])).replace("/", "__")
        )
        tmp.mkdir(parents=True, exist_ok=True)
        full = pd.read_csv(paths["peptide_hits"], sep="\t", dtype=str, keep_default_na=False)
        full.rename(columns={"peptide_len": "peptide_len_aa"}).to_csv(
            tmp / "peptide_hits.tsv", sep="\t", index=False
        )
        paths["peptide_hits"] = tmp / "peptide_hits.tsv"
    m4, issues = load_m4_inputs(paths)
    errs = [i.code for i in issues if i.severity == "error" and "region" not in i.code]
    pep = m4.peptides
    acc = accepted_peptides(m4.peptide_hits, pep, cfg)
    pmap = det.protein_gene_map(
        m4.genes,
        pep,
        m4.peptide_hits,
        m4.protein_hits,
        accepted=acc,
        contaminant_prefixes=cfg.input.contaminant_prefixes,
    )
    tq = bool(pep.loc[~pep.is_decoy.astype(bool), "q_value"].astype(float).notna().any())
    pair = det.pair_decoys(
        pep,
        pmap,
        cfg.input.decoy_prefixes,
        q_value_max=cfg.input.q_value_max if tq else None,
        contaminant_prefixes=cfg.input.contaminant_prefixes,
    )
    mp = max(cfg.mapping.proteins.min_identity, cfg.mapping.proteins.min_coverage)
    tier, w = det.choose_tier(
        pep, m4.peptide_hits, m4.protein_hits, pair, accepted=acc, min_pairing=mp
    )
    tags = m4.genes.locus_tag.astype(str)
    elig = m4.genes[~m4.genes.partial.astype(bool)].locus_tag.astype(str).tolist()
    # route diagnostics: share of target protein IDs that are C2 locus tags
    tset = set(tags)
    tgt = pep[~pep.is_decoy.astype(bool)]
    ids = [
        str(x)
        for v in tgt.protein_ids
        if v is not None
        for x in (v if not isinstance(v, str) else v.split(";"))
        if str(x)
    ]
    lt_share = float(np.mean([x in tset for x in ids])) if ids else float("nan")
    n_ph = len(m4.protein_hits)
    paired_mask = pair.gene.notna()
    if NEW:
        key = det.shared_key_columns(acc, pair)
        q = det._picked_q(elig, acc, pair, key)
    else:
        key = det.shared_key_columns(acc, pair[paired_mask])
        f = paired_mask.sum() / len(pair) if len(pair) else np.nan
        q = det._picked_q(elig, acc, pair, f, key)
    # decoy wins, recomputed independently on the same key ranking
    use = pair if NEW else pair[paired_mask]
    n_t = len(acc)
    cols = {c: np.r_[acc[c].astype(float), use[c].astype(float)] for c in key}
    rank = det._key_rank(cols, n_t + len(use))
    tb = (
        pd.Series(rank[:n_t])
        .groupby(acc.gene.astype(str).to_numpy())
        .max()
        .reindex(elig)
        .to_numpy(float)
    )
    dr = pd.Series(rank[n_t:])
    pm = use.gene.notna().to_numpy()
    gb = dr[pm].groupby(use.gene[pm].astype(str).to_numpy()).max().reindex(elig).to_numpy(float)
    own = (
        dr[~pm].groupby(use["unit"][~pm].astype(str).to_numpy()).max().to_numpy(float)
        if NEW
        else np.array([])
    )
    pr = picked_fdr(np.r_[tb, np.full(own.size, np.nan)], np.r_[gb, own])
    alpha = cfg.stats.alpha.identification
    sel = q <= alpha
    r = dict(
        label=label,
        code=getattr(gippa, "__version__", "?"),
        run=str(o.relative_to(o.parents[2])),
        commit=man["git_commit"][:8],
        recorded_tier=man.get("detection_tier"),
        tier=tier,
        warn=w,
        load_errors=errs,
        protein_hits=n_ph,
        target_id_locus_tag_share=round(lt_share, 4),
        decoy_rows=len(pair),
        gene_paired=int(paired_mask.sum()),
        unpaired_own=int((~paired_mask).sum()),
        own_units=int(
            pair.loc[~paired_mask].shape[0]
            and (pair.loc[~paired_mask, "unit"].nunique() if "unit" in pair else 0)
        ),
        f=round(float(paired_mask.mean()), 4) if len(pair) else None,
        key=list(map(str, key)),
        eligible=len(elig),
        alpha=alpha,
        genes_q01=int(sel.sum()),
        maxq=round(float(np.nanmax(q[sel])), 4) if sel.any() else None,
        target_wins=int(np.isfinite(q).sum()),
        decoy_wins=int(pr.n_decoy_wins),
    )
    dump = Path(__file__).parent / "detected" / label
    dump.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"locus_tag": elig, "q": q}).to_csv(
        dump / (r["run"].replace("/", "__") + ".tsv"), sep="\t", index=False
    )
    print(json.dumps(r), flush=True)
    res.append(r)
