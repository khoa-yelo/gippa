"""Campaign Phase 2 group B metrics, read from finished GIPPA runs (routes A and B).

Writes one JSON of **aggregates only** (counts, shares with Wilson 95 % CIs, q-values, locus
tags, 1-based coordinates); no peptide sequence and no DIA-NN row leaves ``$GIPPA_ROOT``. Row
names follow ``validation/campaign/<ACC>/ANSWER_KEY.md`` (S1-S13 standard, then the dataset's
own rows). TE1, TE3, TE4 and TE5 reuse the pilot's functions at the pinned commit, as Phase 1
did (``validation/campaign/ecoli/phase1_metrics.py``, from which the shared helpers are taken).
Run with the pinned tree on ``PYTHONPATH`` (``<pin>:<pin>/src``)::

    p2b_metrics.py --acc ACC --dir $GIPPA_ROOT/campaign/ACC --bakta-json J --genome G \
        --refseq-gff GFF --contaminants FASTA --out metrics.json
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
import re
import sys
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from validation.pilot import metrics as pm
from validation.pilot import te_real as tr
from validation.proxies.te_proxies import non_atg_nterm, read_coding_genes

from gippa.io.peptides import read_peptide_table
from gippa.stats.cutoffs import wilson_ci

RIBO_ATP = re.compile(r"^(rp[lms][A-Z]|atp[A-I])")
SEED = 20261014
TE4 = {  # (RplL, RplJ, AtpA, AtpD) locus tags, ANSWER_KEY.md loci tables
    "PXD049476": ("KGNFJJ_00023", "KGNFJJ_00022", "KGNFJJ_02722", "KGNFJJ_02724"),
    "PXD065346": ("FLEEGO_00126", "FLEEGO_00125", "FLEEGO_03723", "FLEEGO_03721"),
    "PXD072227": ("AAHINN_02716", "AAHINN_02717", "AAHINN_02307", "AAHINN_02272"),
    "PXD060688": ("BFAOEK_03063", "BFAOEK_03062", "BFAOEK_06952", "BFAOEK_06849"),
}
SYN_MARKERS = {"psaA": "AAHINN_00391", "psaB": "AAHINN_00392", "psbB": "AAHINN_00529",
               "cbbL": "AAHINN_01966", "cbbS": "AAHINN_01965", "mpeA": "AAHINN_00497",
               "mpeB": "AAHINN_00498"}
VNAT_CHR = {"chr1": "NZ_CP009977.1", "chr2": "NZ_CP009978.1"}
CDIF = "NZ_CP016318.1"


# ---------------------------------------------------------------- helpers (as Phase 1)
def share(k: int, n: int) -> dict[str, Any]:
    if n == 0:
        return {"k": 0, "n": 0, "value": None, "ci_low": None, "ci_high": None}
    lo, hi = wilson_ci(k, n)
    return {"k": int(k), "n": int(n), "value": k / n, "ci_low": float(lo), "ci_high": float(hi)}


def newcombe(a: dict[str, Any], b: dict[str, Any]) -> dict[str, Any]:
    """Difference a - b of two Wilson shares, Newcombe hybrid-score 95 % CI."""
    if not a["n"] or not b["n"]:
        return {"diff": None, "ci_low": None, "ci_high": None}
    d = a["value"] - b["value"]
    lo = d - math.sqrt((a["value"] - a["ci_low"]) ** 2 + (b["ci_high"] - b["value"]) ** 2)
    hi = d + math.sqrt((a["ci_high"] - a["value"]) ** 2 + (b["value"] - b["ci_low"]) ** 2)
    return {"diff": d, "ci_low": lo, "ci_high": hi}


def clean(v: Any) -> Any:
    c = pm.clean(v)
    if isinstance(c, float) and math.isnan(c):
        return None
    return c


def tbl(out: Path, stem: str, sub: str = "tables") -> pd.DataFrame | None:
    p = out / sub / f"{stem}.tsv"
    return pm._table(p, stem) if p.is_file() else None


def manifest(out: Path) -> dict[str, Any]:
    return json.loads((out / "run_manifest.json").read_text())


def coding(g: pd.DataFrame) -> pd.DataFrame:
    return g[g["feature_type"].astype(str).isin(["cds", "sorf"])]


def is_det(s: pd.Series) -> pd.Series:
    return s.astype(str) == "detected"


def arm_summary(out: Path) -> dict[str, Any]:
    man = manifest(out)
    cod = coding(tbl(out, "genes"))
    st = Counter(cod["detection_status"].astype(str))
    tier = Counter(cod.loc[is_det(cod["detection_status"]), "detection_tier"].astype(str))
    two = int((is_det(cod["detection_status"]) & cod["detected_2pep"].map(lambda v: str(v).lower() == "true")).sum())
    ds = tbl(out, "detection_summary")
    fdr = None
    if ds is not None:
        r = ds[(ds["section"].astype(str) == "tier") & (ds["group"].astype(str).isin(["0.01", "0.010"]))
               & (ds["metric"].astype(str) == "fdr")]
        fdr = clean(pd.to_numeric(r["value"], errors="coerce").iloc[0]) if len(r) else None
    mf = tbl(out, "mapping_fdr", "mapping")
    mrows = []
    if mf is not None:
        for r in mf[mf["length_class"].astype(str) == "all"].itertuples(index=False):
            mrows.append({k: clean(getattr(r, k)) for k in ("section", "n_diff", "targets", "decoys",
                          "n_accepted", "precision", "ci_low", "ci_high", "flags")})
    return {"commit": man.get("git_commit"), "git_dirty": man.get("git_dirty"),
            "exit_code": man.get("exit_code"), "status": man.get("status"),
            "detection_tier": man.get("detection_tier"), "started_at": man.get("started_at"),
            "finished_at": man.get("finished_at"), "warnings": man.get("warnings"),
            "layers_missing": man.get("layers_missing"), "modules_off": man.get("modules_off"),
            "seeds": man.get("seeds"),
            "inputs": {i.get("role"): i.get("sha256") for i in man.get("inputs", [])},
            "n_coding_genes": int(len(cod)), "detection_status": dict(st), "detected_by_tier": dict(tier),
            "n_detected": int(st.get("detected", 0)), "n_detected_2pep": two,
            "gene_fdr_at_q01": fdr, "mapping_fdr_all": mrows}


def detected(out: Path) -> set[str]:
    g = tbl(out, "genes")
    return set(g.loc[is_det(g["detection_status"]), "locus_tag"].astype(str))


def gene_row(out: Path, tag: str) -> dict[str, Any] | None:
    g = tbl(out, "genes")
    r = g[g["locus_tag"].astype(str) == tag]
    if not len(r):
        return None
    r = r.iloc[0]
    return {k: clean(r[k]) for k in ("locus_tag", "gene", "contig", "start", "end", "strand",
            "length_aa", "detection_status", "detection_tier", "gene_q_value", "n_peptides_unique",
            "n_peptides_shared", "coverage_frac", "abundance_log2", "abundance_scale", "flags")}


def hits_in(out: Path, contig: str, lo: int, hi: int, strand: str) -> pd.DataFrame:
    """peptide_hits rows overlapping [lo, hi] (1-based closed, as the table)."""
    h = tbl(out, "peptide_hits", "mapping")
    s = pd.to_numeric(h["start"], errors="coerce")
    e = pd.to_numeric(h["end"], errors="coerce")
    return h[(h["contig"].astype(str) == contig) & (h["strand"].astype(str) == strand) & (s <= hi) & (e >= lo)]


def refseq(gff: Path) -> pd.DataFrame:
    rows = []
    with gzip.open(gff, "rt") as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("CDS", "pseudogene", "gene"):
                continue
            attr = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            if f[2] == "gene" and attr.get("gene_biotype") not in ("protein_coding", "pseudogene"):
                continue
            rows.append({"type": f[2], "contig": f[0], "start": int(f[3]), "end": int(f[4]), "strand": f[6],
                         "name": attr.get("gene") or attr.get("Name"), "locus_tag": attr.get("locus_tag"),
                         "product": attr.get("product"),
                         "pseudo": attr.get("pseudo") == "true" or f[2] == "pseudogene"
                         or attr.get("gene_biotype") == "pseudogene"})
    return pd.DataFrame(rows)


def _best(ref: pd.DataFrame, contig: str, strand: str, a: int, b: int) -> tuple[int, Any] | None:
    cds = ref[(ref["type"].isin(["CDS", "pseudogene"]) | ref["pseudo"]) & (ref["contig"] == contig)
              & (ref["strand"] == strand) & (ref["start"] <= b) & (ref["end"] >= a)]
    return max(((min(b, c["end"]) - max(a, c["start"]) + 1, c) for _, c in cds.iterrows()),
               key=lambda x: x[0], default=None)


def _orf_span(r: pd.Series) -> tuple[int, int] | None:
    """1-based closed span; without a start codon it is derived from the longest length."""
    st, e = clean(r["start_longest"]), int(r["end"])
    if st is not None:
        a, b = sorted((int(st), e))
        return a, b
    n = clean(r["length_aa_longest"])
    if n is None:
        return None
    return (e, e + 3 * int(n) + 2) if str(r["strand"]) == "-" else (e - 3 * int(n) - 2, e)


def novel(out: Path, ref: pd.DataFrame) -> dict[str, Any]:
    n = tbl(out, "novel_orfs")
    if n is None:
        return {"reason": "missing:novel_orfs.tsv"}
    called = n[n["called"].map(lambda v: str(v).lower() == "true")]
    rows, k = [], 0
    for _, r in called.iterrows():
        sp = _orf_span(r)
        best = _best(ref, str(r["contig"]), str(r["strand"]), *sp) if sp else None
        L = sp[1] - sp[0] + 1 if sp else 0
        match = best is not None and best[0] >= 0.5 * L
        k += match
        rows.append({"orf_id": clean(r["orf_id"]), "class": clean(r["class"]), "contig": clean(r["contig"]),
                     "start": sp[0] if sp else None, "end": sp[1] if sp else int(r["end"]),
                     "strand": r["strand"], "length_aa": clean(r["length_aa_longest"]),
                     "n_peptides": clean(r["n_peptides"]), "best_peptide_q": clean(r["best_peptide_q"]),
                     "class_fdr": clean(r["class_fdr"]), "flags": clean(r["flags"]),
                     "overlapping_gene": clean(r["overlapping_gene"]),
                     "matches_refseq": bool(match), "status": "known_refseq" if match else "new",
                     "refseq_name": None if best is None else best[1]["name"],
                     "refseq_locus": None if best is None else best[1]["locus_tag"],
                     "refseq_pseudo": None if best is None else bool(best[1]["pseudo"]),
                     "overlap_frac": None if best is None or not L else best[0] / L})
    cls = Counter(n["class"].astype(str))
    fdr = {c: clean(pd.to_numeric(n.loc[n["class"].astype(str) == c, "class_fdr"], errors="coerce").min())
           for c in cls}
    # exploratory: every candidate with >= 2 peptides, same matching rule
    multi = []
    for _, r in n[pd.to_numeric(n["n_peptides"], errors="coerce") >= 2].iterrows():
        sp = _orf_span(r)
        best = _best(ref, str(r["contig"]), str(r["strand"]), *sp) if sp else None
        L = sp[1] - sp[0] + 1 if sp else 0
        multi.append({"orf_id": clean(r["orf_id"]), "class": clean(r["class"]), "contig": clean(r["contig"]),
                      "start": sp[0] if sp else None, "end": sp[1] if sp else int(r["end"]), "strand": r["strand"],
                      "length_aa": clean(r["length_aa_longest"]), "n_peptides": clean(r["n_peptides"]),
                      "best_peptide_q": clean(r["best_peptide_q"]), "class_fdr": clean(r["class_fdr"]),
                      "flags": clean(r["flags"]), "overlapping_gene": clean(r["overlapping_gene"]),
                      "refseq_name": None if best is None else best[1]["name"],
                      "refseq_locus": None if best is None else best[1]["locus_tag"],
                      "refseq_product": None if best is None else best[1]["product"],
                      "overlap_frac": None if best is None or not L else round(best[0] / L, 3)})
    return {"n_candidates": int(len(n)), "candidates_by_class": dict(cls), "n_called": int(len(called)),
            "called_by_class": dict(Counter(called["class"].astype(str))), "min_class_fdr_by_class": fdr,
            "share_matching_refseq": share(k, len(called)), "called": rows,
            "candidates_ge2_peptides": multi}


def q_floor(out: Path, nov: dict[str, Any]) -> dict[str, Any]:
    w = [x for x in (manifest(out).get("warnings") or []) if "class_q_floor" in str(x)]
    need = [c for c, n in nov.get("candidates_by_class", {}).items()
            if n < 100 and (nov["min_class_fdr_by_class"].get(c) or 1) > 0.01]
    verdict = "pass" if (not need or w) else "fail"
    return {"warnings": w, "classes_below_100_and_fdr_above_alpha": need, "verdict": verdict}


def prophages(out: Path) -> dict[str, Any]:
    rs = tbl(out, "regions_stats")
    if rs is None:
        return {"reason": "missing:regions_stats.tsv", "verdict": "not testable"}
    regs = []
    for rtype in ("prophage", "bgc", "crispr"):
        p = rs[(rs["region_type"].astype(str) == rtype) & (rs["subset"].astype(str) == "all")]
        for rid, grp in p.groupby("region_id"):
            oe = grp[grp["method"].astype(str) == "circular_shift_o_minus_e"]
            r0 = grp.iloc[0]
            row = {"region_id": rid, "region_type": rtype, "contig": clean(r0["contig"]),
                   "start": clean(r0["start"]), "end": clean(r0["end"]), "n_genes": clean(r0["n_genes"]),
                   "n_eligible": clean(r0["n_eligible"]), "observed": clean(r0["observed"]),
                   "expected": clean(r0["expected"]), "o_over_e": clean(r0["o_over_e"])}
            if len(oe):
                o = oe.iloc[0]
                row |= {"o_minus_e": clean(o["effect"]), "ci_low": clean(o["ci_low"]),
                        "ci_high": clean(o["ci_high"]), "p_value": clean(o["p_value"]), "q_value": clean(o["q_value"])}
            regs.append(row)
    ph = [r for r in regs if r["region_type"] == "prophage"]
    lower = [r for r in ph if r.get("o_minus_e") is not None and r["o_minus_e"] < 0]
    sig = [r for r in lower if r.get("q_value") is not None and r["q_value"] <= 0.05]
    verdict = "not testable" if not ph else "pass" if len(sig) == len(ph) else "fail"
    isl = tbl(out, "islands")
    islands = [] if isl is None else [
        {k: clean(r[k]) for k in ("island_id", "direction", "contig", "start", "end", "n_genes",
                                  "detected_frac", "mean_abundance_log2", "effect", "q_value")
         if k in isl.columns} for _, r in isl.iterrows()]
    sig_all = [r for r in regs if r.get("q_value") is not None and r["q_value"] <= 0.05]
    return {"regions": regs, "n_prophage": len(ph), "n_prophage_lower": len(lower),
            "n_prophage_lower_q05": len(sig), "n_regions_q05_any_type": len(sig_all),
            "verdict": verdict, "islands": islands}


def operons(out: Path) -> dict[str, Any]:
    g = tbl(out, "genes")
    name = dict(zip(g["locus_tag"].astype(str), g["gene"].map(lambda v: clean(v) or "")))
    ra = {t for t, n in name.items() if RIBO_ATP.match(n)}
    pairs, fl = tbl(out, "pairs"), tbl(out, "operon_flags")
    res: dict[str, Any] = {"n_ribo_atp_genes": len(ra)}
    disc: list[dict[str, Any]] = []
    if pairs is not None:
        ins = pairs[pairs["gene_a"].astype(str).isin(ra) & pairs["gene_b"].astype(str).isin(ra)]
        q = pd.to_numeric(ins["q_value"], errors="coerce")
        d = pd.to_numeric(ins["delta_log2"], errors="coerce").abs()
        res |= {"n_pairs": int(len(ins)), "pair_class": dict(Counter(ins["pair_class"].astype(str))),
                "n_pairs_q_le_0.05": int((q <= 0.05).sum()), "median_abs_delta_log2": clean(d.median()),
                "min_q": clean(q.min())}
        for _, r in ins[q <= 0.05].iterrows():
            disc.append({"gene_a": name.get(str(r["gene_a"])), "gene_b": name.get(str(r["gene_b"])),
                         "delta_log2": clean(r["delta_log2"]), "q_value": clean(r["q_value"]), "source": "pairs"})
    if fl is not None:
        f = fl[fl["gene_a"].astype(str).isin(ra) & fl["gene_b"].astype(str).isin(ra)]
        res["operon_flags"] = [{"gene_a": name.get(str(r["gene_a"])), "gene_b": name.get(str(r["gene_b"])),
                                "flags": clean(r["flags"]), "abundance_pattern": clean(r["abundance_pattern"]),
                                "q_value": clean(r["q_value"])} for _, r in f.iterrows()]
        disc += [x | {"source": "operon_flags"} for x in res["operon_flags"]
                 if x["abundance_pattern"] not in ("concordant", "coabundant", "none", None)]
        res["n_flags_all"] = int(len(fl))
        res["flag_counts_all"] = dict(Counter(fl["flags"].astype(str)))
    res["discordant"] = disc
    res["verdict"] = "not testable" if pairs is None else ("pass" if not disc else "fail")
    return res


def te4_tags(out: Path, tags: tuple[str, str, str, str]) -> dict[str, Any]:
    rows = {t: gene_row(out, t) for t in tags}
    res: dict[str, Any] = {"genes": rows, "pairs": {}}
    for (a, b, lo, hi, lab) in ((tags[0], tags[1], 1.0, 3.0, "RplL/RplJ"), (tags[2], tags[3], -1.0, 1.0, "AtpA/AtpD")):
        ra, rb = rows[a], rows[b]
        ok_meas = all(r and r["abundance_log2"] is not None and r["detection_status"] == "detected" for r in (ra, rb))
        v = (ra["abundance_log2"] - rb["abundance_log2"]) if ok_meas else None
        scale = {r["abundance_scale"] for r in (ra, rb) if r}
        res["pairs"][lab] = {"log2_ratio": v, "low": lo, "high": hi, "scale": sorted(map(str, scale)),
                             "status": "not_evaluable" if v is None else ("pass" if lo <= v <= hi else "fail")}
    st = [p["status"] for p in res["pairs"].values()]
    res["verdict"] = "pass" if all(s == "pass" for s in st) else ("not testable" if all(s == "not_evaluable" for s in st) else "fail")
    return res


def variants(out: Path) -> dict[str, Any]:
    v = tbl(out, "variants")
    if v is None or not len(v):
        return {"n_rows": 0, "n_high_alt": 0, "verdict": "pass"}
    fl = v["flags"].astype(str)
    high = v[(v["confidence"].astype(str) == "high") & (v["call"].astype(str).isin(["ALT", "mixed"]))]
    site = high["flags"].astype(str).str.contains("start_codon|stop_codon")
    return {"n_rows": int(len(v)), "n_high_alt": int(len(high)), "n_high_outside_start_stop": int((~site).sum()),
            "n_mass_equivalent": int(fl.str.contains("mass_equivalent").sum()),
            "verdict": "pass" if int((~site).sum()) == 0 else "fail"}


def bootstrap_median_diff(x: np.ndarray, y: np.ndarray, n: int = 2000) -> dict[str, Any]:
    rng = np.random.default_rng(SEED)
    d = float(np.median(x) - np.median(y))
    bs = [np.median(rng.choice(x, len(x))) - np.median(rng.choice(y, len(y))) for _ in range(n)]
    lo, hi = np.percentile(bs, [2.5, 97.5])
    return {"median_x": float(np.median(x)), "median_y": float(np.median(y)), "diff": d,
            "ci_low": float(lo), "ci_high": float(hi), "n_x": int(len(x)), "n_y": int(len(y)), "n_boot": n}


def measured(g: pd.DataFrame) -> pd.Series:
    return pd.to_numeric(g["abundance_log2"], errors="coerce")


def q75(g: pd.DataFrame) -> float:
    return float(measured(coding(g)).dropna().quantile(0.75))


# ---------------------------------------------------------------- dataset rows
def vnat(a: Path, b: Path) -> dict[str, Any]:
    ga, gb = coding(tbl(a, "genes")), coding(tbl(b, "genes"))
    res: dict[str, Any] = {"per_chromosome": {}}
    detb = set(gb.loc[is_det(gb["detection_status"]), "locus_tag"].astype(str))
    for k, c in VNAT_CHR.items():
        s = ga[ga["contig"].astype(str) == c]
        d = is_det(s["detection_status"])
        res["per_chromosome"][k] = {"contig": c, "n_coding": int(len(s)), "detected": share(int(d.sum()), len(s)),
                                    "median_abundance_log2": clean(measured(s[d]).median()),
                                    "B_kept": share(len(set(s.loc[d, "locus_tag"].astype(str)) & detb), int(d.sum()))}
    c1, c2 = res["per_chromosome"]["chr1"], res["per_chromosome"]["chr2"]
    v1 = newcombe(c1["detected"], c2["detected"])
    res["V1"] = v1 | {"verdict": "pass" if v1["diff"] is not None and v1["ci_low"] > 0 else "fail"}
    x = measured(ga[(ga["contig"].astype(str) == VNAT_CHR["chr1"]) & is_det(ga["detection_status"])]).dropna().to_numpy()
    y = measured(ga[(ga["contig"].astype(str) == VNAT_CHR["chr2"]) & is_det(ga["detection_status"])]).dropna().to_numpy()
    v2 = bootstrap_median_diff(x, y)
    res["V2"] = v2 | {"verdict": "pass" if v2["diff"] > 0 and v2["ci_low"] > 0 else "fail"}
    return res


def syn(a: Path) -> dict[str, Any]:
    g = tbl(a, "genes")
    rows = {k: gene_row(a, t) for k, t in SYN_MARKERS.items()}
    det = [k for k, r in rows.items() if r and r["detection_status"] == "detected"]
    ab = [r["abundance_log2"] for r in rows.values() if r and r["abundance_log2"] is not None]
    med = float(np.median(ab)) if ab else None
    thr = q75(g)
    ok = len(det) == len(rows) and med is not None and med >= thr
    pct = float((measured(coding(g)).dropna() <= med).mean()) if med is not None else None
    return {"markers": rows, "n_detected": len(det), "n": len(rows), "median_abundance_log2": med,
            "q75_measured": thr, "median_percentile": pct, "verdict": "pass" if ok else "fail"}


def mxan(a: Path, b: Path, rs: dict[str, Any]) -> dict[str, Any]:
    g = coding(tbl(a, "genes"))
    rid = g["region_ids"].astype(str)
    inb = rid.str.contains("bgc_")
    d = is_det(g["detection_status"])
    si, so = share(int((d & inb).sum()), int(inb.sum())), share(int((d & ~inb).sum()), int((~inb).sum()))
    bgc = [r for r in rs["regions"] if r["region_type"] == "bgc"]
    per = []
    for r in bgc:
        m = rid.map(lambda v, k=r["region_id"]: k in re.split(r"[;,]", str(v)))
        per.append(r | {"n_coding_in_table": int(m.sum()), "n_detected": int((m & d).sum())})
    lm = {k: manifest(x).get("layers_missing") for k, x in (("A", a), ("B", b))}
    bad = any(any("antismash" in str(z).lower() or "bgc" in str(z).lower() for z in (v or [])) for v in lm.values())
    return {"M2": {"inside": si, "outside": so, "diff_inside_minus_outside": newcombe(si, so),
                   "n_regions": len(bgc), "regions_with_detected": sum(1 for r in per if r["n_detected"] > 0),
                   "n_regions_q05": sum(1 for r in per if (r.get("q_value") or 1) <= 0.05), "regions": per},
            "M3": {"layers_missing": lm, "verdict": "fail" if bad or not bgc else "pass"}}


def cdif(a: Path, b: Path, pa: Path, pb: Path, nov_b: dict[str, Any]) -> dict[str, Any]:
    res: dict[str, Any] = {}
    for lab, out in (("A", a), ("B", b)):
        g = coding(tbl(out, "genes"))
        L = pd.to_numeric(g["length_aa"], errors="coerce")
        s = g[L <= 100]
        res.setdefault("C1", {})[lab] = share(int(is_det(s["detection_status"]).sum()), len(s))
    res["C2"] = {"called_le_100aa": [r for r in nov_b.get("called", []) if (r.get("length_aa") or 1e9) <= 100],
                 "n_called": nov_b.get("n_called"), "n_new": sum(1 for r in nov_b.get("called", []) if r["status"] == "new")}
    c3 = {}
    for lab, out, pep in (("A", a, pa), ("B", b, pb)):
        P = read_peptide_table(pep)[0]
        acc = pm._target(P) & (pm._q(P) <= pm.Q_ACCEPT)
        sp: dict[str, set[str]] = {}
        for p, ids, ok in zip(P["peptide"].astype(str), P["protein_ids"], acc, strict=True):
            if ok:
                for i in pm.items(ids):
                    if i.startswith("SPIKE_"):
                        sp.setdefault(i, set()).add(p)
        allp = set().union(*sp.values()) if sp else set()
        h = tbl(out, "peptide_hits", "mapping")
        placed = set(h["peptide"].astype(str)) & allp if h is not None else set()
        c3[lab] = {"per_protein": {k: len(v) for k, v in sorted(sp.items())}, "n_peptides": len(allp),
                   "n_placed_on_genome": len(placed),
                   "verdict": "not testable" if not allp else ("pass" if not placed else "fail")}
    res["C3"] = c3
    # C4 GrdB: both sides of the Sec TGA 2,677,827-2,677,829 (-)
    gr = gene_row(a, "FLEEGO_02514")
    h = hits_in(a, CDIF, 2677566, 2678876, "-")
    ex = h[h["route"].astype(str) == "exact"]
    s = pd.to_numeric(ex["start"], errors="coerce")
    e = pd.to_numeric(ex["end"], errors="coerce")
    up = int(ex.loc[s > 2677829, "peptide"].nunique())      # minus strand: N-terminal side is right
    dn = int(ex.loc[e < 2677827, "peptide"].nunique())
    ac = int(ex.loc[(s <= 2677827) & (e >= 2677829), "peptide"].nunique())
    det = gr is not None and gr["detection_status"] == "detected"
    res["C4"] = {"gene": gr, "peptides_nterm_side": up, "peptides_cterm_side": dn, "peptides_across": ac,
                 "verdict": "not testable" if not det else ("pass" if up and dn and not ac else "fail")}
    # C5 GrdA fragments and gap
    c5 = {"genes": {t: gene_row(a, t) for t in ("FLEEGO_02516", "FLEEGO_02515")}}
    for lab, out in (("A", a), ("B", b)):
        c5[lab] = {w: int(hits_in(out, CDIF, lo, hi, "-")["peptide"].nunique())
                   for w, (lo, hi) in (("nterm_fragment", (2679250, 2679381)), ("gap", (2679232, 2679246)),
                                        ("cterm_fragment", (2678908, 2679231)))}
    c5["B_candidates"] = novel_overlapping(b, CDIF, 2678908, 2679381, "-")
    res["C5"] = c5
    # C6 PrdB C-terminal tail 3,749,919-3,750,188 (-)
    hb = hits_in(b, CDIF, 3749919, 3750188, "-")
    cands = novel_overlapping(b, CDIF, 3749919, 3750188, "-")
    n = int(hb["peptide"].nunique())
    res["C6"] = {"gene_A": gene_row(a, "FLEEGO_03446"), "B_peptides_in_tail": n,
                 "B_locus_class": dict(Counter(hb["locus_class"].astype(str))) if n else {},
                 "B_candidates": cands,
                 "verdict": "not testable" if not n else ("pass" if cands and all(
                     str(x) not in ("cds", "annotated") for x in hb["locus_class"].astype(str)) else "fail")}
    return res


def novel_overlapping(out: Path, contig: str, lo: int, hi: int, strand: str) -> list[dict[str, Any]]:
    n = tbl(out, "novel_orfs")
    if n is None:
        return []
    s = pd.to_numeric(n["start_longest"], errors="coerce")
    e = pd.to_numeric(n["end"], errors="coerce")
    a = pd.concat([s, e], axis=1).min(axis=1)
    b = pd.concat([s, e], axis=1).max(axis=1)
    sel = n[((n["contig"].astype(str) == contig) & (n["strand"].astype(str) == strand) & (a <= hi) & (b >= lo)).fillna(False).astype(bool)]
    keep = ("orf_id", "class", "strand", "start_longest", "end", "length_aa_longest", "n_peptides",
            "best_peptide_q", "class_fdr", "called", "overlapping_gene", "flags")
    return [{k: clean(r[k]) for k in keep} for _, r in sel.iterrows()]


# ---------------------------------------------------------------- main
def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--acc", required=True)
    ap.add_argument("--dir", type=Path, required=True)
    ap.add_argument("--bakta-json", type=Path, required=True)
    ap.add_argument("--genome", type=Path, required=True)
    ap.add_argument("--refseq-gff", type=Path, required=True)
    ap.add_argument("--contaminants", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    x = ap.parse_args(argv)
    D = x.dir
    A, B = D / "out_a", D / "out_b"
    pa, pb = D / "in_a" / "peptides.tsv.gz", D / "in_b" / "peptides.tsv.gz"
    doc: dict[str, Any] = {"schema": "gippa-campaign-p2b/1", "acc": x.acc}
    ref = refseq(x.refseq_gff)
    doc["A"], doc["B"] = arm_summary(A), arm_summary(B)
    man = manifest(A)
    cfg = pm._cfg(man, A / "run_manifest.json")
    gA = tbl(A, "genes")
    P = read_peptide_table(pa)[0]
    hits, unm = tbl(A, "peptide_hits", "mapping"), tbl(A, "unmapped", "mapping")
    cod, circ = read_coding_genes(x.bakta_json, x.genome)
    te1, exc = pm.te1_real(P, hits, unm, cod, circ, non_atg_nterm(cod), cfg)
    doc["S1_te1"] = te1 | {"n_exceptions": len(exc)}
    t3 = tr.te3(gA)
    t3.pop("_rows", None)
    doc["S2_te3"] = t3
    doc["S3_te4"] = te4_tags(A, TE4[x.acc]) | {"te4_name_based": tr.te4(gA)}
    heads = [ln[1:] for ln in x.contaminants.read_text().splitlines() if ln.startswith(">")]
    fam = {f: set(ids) | {"CONT_" + i for i in ids} for f, ids in tr.contaminant_families(heads).items()}
    t5 = tr.te5_arm(P, hits, unm, fam, cfg)
    t5.pop("_violations", None)
    doc["S4_te5"] = t5
    fa, fb = doc["A"]["gene_fdr_at_q01"], doc["B"]["gene_fdr_at_q01"]
    doc["S5_gene_fdr"] = {"A": fa, "B": fb, "verdict": "pass" if fa is not None and fb is not None
                          and fa <= 0.01 and fb <= 0.01 else "fail"}
    doc["S6_regions"] = prophages(A)
    doc["S7_operons"] = operons(A)
    nov = novel(B, ref)
    doc["S8_novel"] = nov
    doc["S9_q_floor"] = q_floor(B, nov)
    d1, d2 = detected(A), detected(B)
    doc["S10_route_b"] = {"A_detected": len(d1), "B_detected": len(d2), "A_kept_by_B": share(len(d1 & d2), len(d1)),
                          "B_only": len(d2 - d1)}
    cmp_dir = D / "cmp_ab" / "tables"
    if cmp_dir.is_dir():
        doc["compare"] = {}
        for f in sorted(cmp_dir.glob("compare_*.tsv")):
            t = pd.read_csv(f, sep="\t")
            if "section" in t.columns:  # drop the per-gene detection matrix
                t = t[t["section"].astype(str) != "matrix"]
            doc["compare"][f.stem] = [{k: clean(v) for k, v in r.items()} for r in t.head(60).to_dict("records")]
    doc["S12_variants"] = {"A": variants(A), "B": variants(B)}
    te7 = D / "te7.diff"
    if te7.exists():
        doc["S13_te7"] = {"diff_bytes": te7.stat().st_size,
                          "n_files": len(list((D / "out_a_te7" / "tables").glob("*"))),
                          "verdict": "pass" if te7.stat().st_size == 0 else "fail"}
    if (D / "out_r2").is_dir():
        doc["S11_round2"] = arm_summary(D / "out_r2")
    if x.acc == "PXD049476":
        doc["V"] = vnat(A, B)
    elif x.acc == "PXD072227":
        doc["Y1"] = syn(A)
    elif x.acc == "PXD060688":
        doc |= mxan(A, B, doc["S6_regions"])
    elif x.acc == "PXD065346":
        doc["C"] = cdif(A, B, pa, pb, nov)
    x.out.write_text(json.dumps(doc, indent=1, default=str) + "\n")
    print(f"wrote {x.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
