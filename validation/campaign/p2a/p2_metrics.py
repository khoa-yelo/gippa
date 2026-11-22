"""Campaign Phase 2 group A metrics, read from finished GIPPA runs (aggregates only).

Generalised from ``validation/campaign/ecoli/phase1_metrics.py`` (branch validation/campaign-ecoli).
Writes one JSON of counts, shares with Wilson 95 % CIs, q-values, locus tags and 1-based
coordinates: no peptide sequence and no DIA-NN row leaves ``$GIPPA_ROOT``. Row names follow the
dataset's ``ANSWER_KEY.md`` (S1-S11 standard; T*, B*, C* dataset checks). TE1/TE3/TE4/TE5 reuse the
pilot's pure functions at the pinned commit. Run with ``PYTHONPATH=<pin>:<pin>/src``::

    p2_metrics.py --a OUT_A --pep-a PEP_A [--b OUT_B --cmp CMP] --bakta-json J --genome G \
        --refseq-gff GFF --contaminants FASTA --loci answer_loci.tsv --tag TAG --out metrics.json

Verdict rules are those of the answer key; S5 as Phase 1 row 6, S7 as Phase 1 row 9.
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

import pandas as pd
from validation.pilot import metrics as pm
from validation.pilot import te_real as tr
from validation.proxies.te_proxies import non_atg_nterm, read_coding_genes

from gippa.io.peptides import read_peptide_table
from gippa.stats.cutoffs import wilson_ci

RIBO_ATP = re.compile(r"^(rp[lms][A-Z]|atp[A-I])")


def share(k: int, n: int) -> dict[str, Any]:
    if n == 0:
        return {"k": 0, "n": 0, "value": None, "ci_low": None, "ci_high": None}
    lo, hi = wilson_ci(k, n)
    return {"k": k, "n": n, "value": k / n, "ci_low": float(lo), "ci_high": float(hi)}


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


def arm_summary(out: Path) -> dict[str, Any]:
    man = manifest(out)
    g = tbl(out, "genes")
    cod = g[g["feature_type"].astype(str).isin(["cds", "sorf"])]
    st = Counter(cod["detection_status"].astype(str))
    tier = Counter(cod.loc[cod["detection_status"].astype(str) == "detected", "detection_tier"].astype(str))
    q = pd.to_numeric(cod["gene_q_value"], errors="coerce")
    fdr = tbl(out, "mapping_fdr", "mapping")
    fdr_rows = []
    if fdr is not None:
        f = fdr[(fdr["length_class"].astype(str) == "all")]
        for r in f.itertuples(index=False):
            fdr_rows.append({k: clean(getattr(r, k)) for k in ("section", "n_diff", "targets",
                             "decoys", "n_accepted", "precision", "ci_low", "ci_high", "flags")})
    return {
        "commit": man.get("git_commit"), "git_dirty": man.get("git_dirty"),
        "exit_code": man.get("exit_code"), "status": man.get("status"),
        "detection_tier": man.get("detection_tier"), "started_at": man.get("started_at"),
        "finished_at": man.get("finished_at"), "warnings": man.get("warnings"),
        "layers_missing": man.get("layers_missing"), "modules_off": man.get("modules_off"),
        "seeds": man.get("seeds"),
        "inputs": {i.get("role"): i.get("sha256") for i in man.get("inputs", [])},
        "n_coding_genes": int(len(cod)), "detection_status": dict(st), "detected_by_tier": dict(tier),
        "n_detected": int(st.get("detected", 0)),
        "n_gene_q_le_0.01": int((q <= 0.01).sum()),
        "mapping_fdr_all": fdr_rows,
    }


def detected(out: Path) -> set[str]:
    g = tbl(out, "genes")
    return set(g.loc[g["detection_status"].astype(str) == "detected", "locus_tag"].astype(str))


def gene_row(out: Path, tag: str) -> dict[str, Any] | None:
    g = tbl(out, "genes")
    r = g[g["locus_tag"].astype(str) == tag]
    if not len(r):
        return None
    r = r.iloc[0]
    return {k: clean(r[k]) for k in ("locus_tag", "gene", "start", "end", "strand", "detection_status",
            "detection_tier", "gene_q_value", "n_peptides_unique", "n_peptides_shared",
            "coverage_frac", "abundance_log2", "flags")}


def hits_in(out: Path, contig: str, lo: int, hi: int, strand: str) -> pd.DataFrame:
    h = tbl(out, "peptide_hits", "mapping")
    s = pd.to_numeric(h["start"], errors="coerce")
    e = pd.to_numeric(h["end"], errors="coerce")
    return h[(h["contig"].astype(str) == contig) & (h["strand"].astype(str) == strand)
             & (s <= hi) & (e >= lo)]


def window_counts(out: Path, contig: str, windows: dict[str, tuple[int, int]], strand: str) -> dict[str, Any]:
    res = {}
    for name, (lo, hi) in windows.items():
        h = hits_in(out, contig, lo, hi, strand)
        ex = h[h["route"].astype(str) == "exact"]
        res[name] = {"window": [lo, hi], "strand": strand,
                     "distinct_peptides_exact": int(ex["peptide"].nunique()),
                     "distinct_peptides_any_route": int(h["peptide"].nunique()),
                     "locus_class": dict(Counter(ex["locus_class"].astype(str))),
                     "best_q": clean(pd.to_numeric(ex["best_q_input"], errors="coerce").min())
                     if len(ex) else None}
    return res


def novel_in(out: Path, lo: int, hi: int, strand: str) -> list[dict[str, Any]]:
    n = tbl(out, "novel_orfs")
    if n is None:
        return []
    s = pd.to_numeric(n["start_longest"], errors="coerce")
    e = pd.to_numeric(n["end"], errors="coerce")
    lo_c = pd.concat([s, e], axis=1).min(axis=1)
    hi_c = pd.concat([s, e], axis=1).max(axis=1)
    sel = n[((n["strand"].astype(str) == strand) & (lo_c <= hi) & (hi_c >= lo)).fillna(False).astype(bool)]
    keep = ("orf_id", "class", "strand", "start_longest", "start_best_ostir", "end", "length_aa_longest",
            "n_peptides", "best_peptide_q", "class_fdr", "called", "overlapping_gene", "flags")
    return [{k: clean(r[k]) for k in keep} for _, r in sel.iterrows()]


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
            rows.append({"type": f[2], "start": int(f[3]), "end": int(f[4]), "strand": f[6],
                         "name": attr.get("gene") or attr.get("Name"), "locus_tag": attr.get("locus_tag"),
                         "pseudo": attr.get("pseudo") == "true" or f[2] == "pseudogene"
                         or attr.get("gene_biotype") == "pseudogene"})
    return pd.DataFrame(rows)


def novel_vs_refseq(out: Path, ref: pd.DataFrame) -> dict[str, Any]:
    n = tbl(out, "novel_orfs")
    if n is None:
        return {"reason": "missing:novel_orfs.tsv"}
    called = n[n["called"].map(lambda v: str(v).lower() == "true")]
    res_all = _match(n, ref) if True else None
    cds = ref[ref["type"].isin(["CDS", "pseudogene"]) | ref["pseudo"]]
    out_rows, k = [], 0
    for _, r in called.iterrows():
        a, b = sorted((int(r["start_longest"]), int(r["end"])))
        L = b - a + 1
        same = cds[(cds["strand"] == r["strand"]) & (cds["start"] <= b) & (cds["end"] >= a)]
        best = None
        for _, c in same.iterrows():
            ov = min(b, c["end"]) - max(a, c["start"]) + 1
            if best is None or ov > best[0]:
                best = (ov, c)
        match = best is not None and best[0] >= 0.5 * L
        k += match
        out_rows.append({"orf_id": clean(r["orf_id"]), "class": clean(r["class"]), "start": a, "end": b,
                         "strand": r["strand"], "n_peptides": clean(r["n_peptides"]),
                         "class_fdr": clean(r["class_fdr"]), "matches_refseq": bool(match),
                         "refseq_name": None if best is None else best[1]["name"],
                         "refseq_locus": None if best is None else best[1]["locus_tag"],
                         "refseq_pseudo": None if best is None else bool(best[1]["pseudo"]),
                         "overlap_frac": None if best is None else best[0] / L})
    cls = Counter(n["class"].astype(str))
    ccls = Counter(called["class"].astype(str))
    fdr = {c: clean(pd.to_numeric(n.loc[n["class"].astype(str) == c, "class_fdr"], errors="coerce").min())
           for c in cls}
    return {"n_candidates": int(len(n)), "candidates_by_class": dict(cls), "n_called": int(len(called)),
            "called_by_class": dict(ccls), "min_class_fdr_by_class": fdr,
            "share_matching_refseq": share(k, len(called)), "called": out_rows,
            "exploratory_all_candidates": res_all}


def _match(n: pd.DataFrame, ref: pd.DataFrame) -> dict[str, Any]:
    """Exploratory (not the pre-registered row 9): every candidate, called or not, against the
    RefSeq CDS/pseudogenes with the same >= 50 % same-strand overlap rule, by class and by
    n_peptides >= 2."""
    cds = ref[ref["type"].isin(["CDS", "pseudogene"]) | ref["pseudo"]]
    rows = []
    for _, r in n.iterrows():
        e = int(r["end"])
        st = clean(r["start_longest"])
        a, b = sorted((int(st), e)) if st is not None else (e, e)
        L = b - a + 1
        same = cds[(cds["strand"] == r["strand"]) & (cds["start"] <= b) & (cds["end"] >= a)]
        best = max(((min(b, c["end"]) - max(a, c["start"]) + 1, c) for _, c in same.iterrows()),
                   key=lambda x: x[0], default=None)
        rows.append({"class": str(r["class"]), "n_peptides": int(clean(r["n_peptides"]) or 0),
                     "match": best is not None and st is not None and best[0] >= 0.5 * L,
                     "pseudo": bool(best[1]["pseudo"]) if best is not None else False})
    out: dict[str, Any] = {}
    for key, sub in (("all", rows), ("n_peptides_ge_2", [x for x in rows if x["n_peptides"] >= 2])):
        out[key] = share(sum(x["match"] for x in sub), len(sub))
        out[key + "_by_class"] = {c: share(sum(x["match"] for x in sub if x["class"] == c),
                                           sum(x["class"] == c for x in sub))
                                  for c in sorted({x["class"] for x in sub})}
        out[key + "_pseudogene_matches"] = sum(x["match"] and x["pseudo"] for x in sub)
    return out


def operons(out: Path) -> dict[str, Any]:
    """Pairs of adjacent genes that are both ribosomal-protein or ATP-synthase genes: the
    within-operon abundance test of ``pairs.tsv`` (q <= 0.05 = discordant) and any
    ``operon_flags.tsv`` row whose abundance pattern is not concordant/coabundant."""
    g = tbl(out, "genes")
    name = dict(zip(g["locus_tag"].astype(str), g["gene"].map(lambda v: clean(v) or "")))
    ra = {t for t, n in name.items() if RIBO_ATP.match(n)}
    pairs = tbl(out, "pairs")
    fl = tbl(out, "operon_flags")
    res: dict[str, Any] = {"n_ribo_atp_genes": len(ra)}
    disc: list[dict[str, Any]] = []
    if pairs is not None:
        inside = pairs[pairs["gene_a"].astype(str).isin(ra) & pairs["gene_b"].astype(str).isin(ra)]
        q = pd.to_numeric(inside["q_value"], errors="coerce")
        d = pd.to_numeric(inside["delta_log2"], errors="coerce").abs()
        both = inside["detected_a"].map(lambda v: str(v).lower() == "true") & inside["detected_b"].map(lambda v: str(v).lower() == "true")
        res |= {"n_pairs": int(len(inside)), "n_pairs_both_detected": int(both.sum()),
                "pair_class": dict(Counter(inside["pair_class"].astype(str))),
                "n_pairs_q_le_0.05": int((q <= 0.05).sum()),
                "median_abs_delta_log2": clean(d.median()), "max_abs_delta_log2": clean(d.max()),
                "min_q": clean(q.min())}
        for _, r in inside[q <= 0.05].iterrows():
            disc.append({"gene_a": name.get(str(r["gene_a"])), "gene_b": name.get(str(r["gene_b"])),
                         "delta_log2": clean(r["delta_log2"]), "q_value": clean(r["q_value"]), "source": "pairs"})
    if fl is not None:
        f = fl[fl["gene_a"].astype(str).isin(ra) & fl["gene_b"].astype(str).isin(ra)]
        res["operon_flags"] = [{"gene_a": name.get(str(r["gene_a"])), "gene_b": name.get(str(r["gene_b"])),
                                "flags": clean(r["flags"]), "abundance_pattern": clean(r["abundance_pattern"]),
                                "q_value": clean(r["q_value"])} for _, r in f.iterrows()]
        for x in res["operon_flags"]:
            if x["abundance_pattern"] not in ("concordant", "coabundant", "none", None):
                disc.append(x | {"source": "operon_flags"})
        res["n_flags_all"] = int(len(fl))
        res["flag_counts_all"] = dict(Counter(fl["flags"].astype(str)))
    res["discordant"] = disc
    res["verdict"] = "not testable" if pairs is None else ("pass" if not disc else "fail")
    return res




def coding_genes(out: Path) -> pd.DataFrame:
    g = tbl(out, "genes")
    return g[g["feature_type"].astype(str).isin(["cds", "sorf"])].copy()


def is_det(g: pd.DataFrame) -> pd.Series:
    return g["detection_status"].astype(str) == "detected"


def host_share(pep: Path, prefixes: tuple[str, ...] = ("CONT_",)) -> dict[str, Any]:
    """T1: distinct accepted (q <= 0.01) target peptides whose proteins are all CONT_ records,
    split into rabbit host (CONT_RABIT_) and the contaminant library."""
    P = read_peptide_table(pep)[0]
    acc = pm._target(P) & (pm._q(P) <= pm.Q_ACCEPT)
    kinds: dict[str, set[str]] = {"all": set(), "cont_only": set(), "rabbit_only": set(), "library_only": set()}
    for p, ids, ok in zip(P["peptide"].astype(str), P["protein_ids"], acc, strict=True):
        if not ok:
            continue
        it = [i for i in pm.items(ids) if i]
        kinds["all"].add(p)
        if it and all(i.startswith(prefixes) for i in it):
            kinds["cont_only"].add(p)
            if all(i.startswith("CONT_RABIT_") for i in it):
                kinds["rabbit_only"].add(p)
            elif not any(i.startswith("CONT_RABIT_") for i in it):
                kinds["library_only"].add(p)
    n = len(kinds["all"])
    return {"n_accepted_peptides": n, "cont_only": share(len(kinds["cont_only"]), n),
            "rabbit_only": share(len(kinds["rabbit_only"]), n),
            "library_only": share(len(kinds["library_only"]), n)}


def loci_status(outs: dict[str, Path], loci: list[dict[str, str]]) -> list[dict[str, Any]]:
    """Answer-key loci: the same-stop Bakta CDS's gene row in each given run (route label -> out),
    plus exact peptides inside the locus span on its strand and frame (route B when present)."""
    res = []
    for r in loci:
        row: dict[str, Any] = {k: r[k] for k in ("key", "label", "genome", "contig", "start", "end", "strand",
                                                 "bakta_locus", "match")}
        for lab, out in outs.items():
            gr = gene_row(out, r["bakta_locus"]) if r["bakta_locus"] and r["match"] == "same_stop" else None
            row[f"{lab}_gene"] = None if gr is None else {k: gr[k] for k in (
                "detection_status", "gene_q_value", "n_peptides_unique", "n_peptides_shared", "abundance_log2", "flags")}
            a, b = int(r["start"]), int(r["end"])
            h = hits_in(out, r["contig"], a, b, r["strand"])
            h = h[h["route"].astype(str) == "exact"]
            s = pd.to_numeric(h["start"], errors="coerce")
            e = pd.to_numeric(h["end"], errors="coerce")
            inframe = ((s - a) % 3 == 0) if r["strand"] == "+" else ((b - e) % 3 == 0)
            # peptide_hits coordinates are 1-based closed (writer convention), as the loci
            row[f"{lab}_peptides_in_span_inframe"] = int(h.loc[inframe.fillna(False).astype(bool), "peptide"].nunique())
            row[f"{lab}_peptides_in_span"] = int(h["peptide"].nunique())
        res.append(row)
    return res


def detected_q(gr: dict[str, Any] | None) -> bool:
    if not gr or gr.get("detection_status") != "detected":
        return False
    q = gr.get("gene_q_value")
    return q is None or float(q) <= 0.01


def per_contig(out: Path) -> dict[str, Any]:
    g = coding_genes(out)
    d = is_det(g)
    return {str(c): share(int(d[g["contig"].astype(str) == c].sum()), int((g["contig"].astype(str) == c).sum()))
            for c in sorted(g["contig"].astype(str).unique())}


def abundance_percentile(out: Path, tags: list[str]) -> dict[str, Any]:
    g = coding_genes(out)
    g = g[is_det(g)].copy()
    ab = pd.to_numeric(g["abundance_log2"], errors="coerce")
    g["pct"] = ab.rank(pct=True)
    sub = g[g["locus_tag"].astype(str).isin(tags)]
    pct = sub["pct"].dropna()
    return {"n_tags": len(tags), "n_detected": int(len(sub)),
            "median_percentile": clean(pct.median()) if len(pct) else None,
            "percentiles": {str(t): clean(p) for t, p in zip(sub["locus_tag"], sub["pct"])}}


def regions_summary(out: Path) -> dict[str, Any]:
    rs = tbl(out, "regions_stats")
    isl = tbl(out, "islands")
    res: dict[str, Any] = {}
    if rs is not None:
        q = pd.to_numeric(rs["q_value"], errors="coerce")
        sig = rs[q <= 0.05]
        res["n_region_tests"] = int(len(rs))
        res["n_q_le_0.05"] = int(len(sig))
        res["by_type_tested"] = dict(Counter(rs["region_type"].astype(str)))
        keep = [c for c in ("region_id", "region_type", "subset", "method", "start", "end", "n_genes",
                            "observed", "expected", "o_over_e", "effect", "q_value") if c in rs.columns]
        res["significant"] = [{k: clean(r[k]) for k in keep} for _, r in sig.iterrows()]
    if isl is not None:
        res["islands"] = [{k: clean(r[k]) for k in ("island_id", "direction", "start", "end", "n_genes",
                           "detected_frac", "mean_abundance_log2", "effect", "q_value") if k in isl.columns}
                          for _, r in isl.iterrows()]
    fl = tbl(out, "operon_flags")
    if fl is not None:
        res["operon_flags"] = {"n": int(len(fl)), "by_flag": dict(Counter(fl["flags"].astype(str)))}
    return res


def gene_fdr(out: Path) -> dict[str, Any]:
    ds = tbl(out, "detection_summary")
    return {"detection_summary": None if ds is None else ds.to_dict("records")}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--a", type=Path, required=True)
    ap.add_argument("--pep-a", type=Path, required=True)
    ap.add_argument("--b", type=Path)
    ap.add_argument("--cmp", type=Path)
    ap.add_argument("--patched", type=Path)
    ap.add_argument("--bakta-json", type=Path, required=True)
    ap.add_argument("--genome", type=Path, required=True)
    ap.add_argument("--refseq-gff", type=Path, required=True)
    ap.add_argument("--contaminants", type=Path, required=True)
    ap.add_argument("--loci", type=Path)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--host", action="store_true", help="T1 host share (rabbit)")
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    doc: dict[str, Any] = {"schema": "gippa-campaign-phase2a/1", "tag": a.tag}
    ref = refseq(a.refseq_gff)
    A = a.a
    man = manifest(A)
    cfg = pm._cfg(man, A / "run_manifest.json")
    genes = tbl(A, "genes")
    P = read_peptide_table(a.pep_a)[0]
    hits, unm = tbl(A, "peptide_hits", "mapping"), tbl(A, "unmapped", "mapping")
    doc["A"] = arm_summary(A)
    doc["A_per_contig"] = per_contig(A)
    doc["S10_A"] = gene_fdr(A)
    coding, circular = read_coding_genes(a.bakta_json, a.genome)
    te1, exc = pm.te1_real(P, hits, unm, coding, circular, non_atg_nterm(coding), cfg)
    doc["S1_te1"] = te1 | {"n_exceptions": len(exc)}
    t3 = tr.te3(genes)
    t3.pop("_rows", None)
    doc["S2_te3"] = t3
    doc["S3_te4"] = tr.te4(genes)
    heads = [ln[1:] for ln in a.contaminants.read_text().splitlines() if ln.startswith(">")]
    fam = {f: set(ids) | {"CONT_" + i for i in ids} for f, ids in tr.contaminant_families(heads).items()}
    t5 = tr.te5_arm(P, hits, unm, fam, cfg)
    t5.pop("_violations", None)
    doc["S4_te5"] = t5
    doc["S5_operons"] = operons(A)
    doc["S6_context_A"] = regions_summary(A)
    if a.host:
        doc["T1_host"] = host_share(a.pep_a)
    outs = {"A": A}
    if a.b:
        B = a.b
        outs["B"] = B
        doc["B"] = arm_summary(B)
        doc["S10_B"] = gene_fdr(B)
        doc["S6_context_B"] = regions_summary(B)
        doc["S7_novel"] = novel_vs_refseq(B, ref)
        doc["S8_q_floor"] = {"warnings": [x for x in (doc["B"]["warnings"] or []) if "class_q_floor" in str(x)]}
        d1, d2 = detected(A), detected(B)
        doc["S9_route_b"] = {"A_detected": len(d1), "B_detected": len(d2),
                             "A_kept_by_B": share(len(d1 & d2), len(d1)), "B_only": len(d2 - d1)}
    if a.cmp:
        ab = a.cmp / "tables" / "compare_abundance.tsv"
        st = a.cmp / "tables" / "compare_stability.tsv"
        doc["compare"] = {"abundance": pd.read_csv(ab, sep="\t").to_dict("records") if ab.is_file() else None,
                          "stability": pd.read_csv(st, sep="\t").to_dict("records") if st.is_file() else None}
    if a.loci:
        import csv
        rows = [r for r in csv.DictReader(a.loci.open(), delimiter="\t") if r["genome"] == a.tag]
        doc["loci"] = loci_status(outs, rows)
    if a.patched:
        doc["S11_round2"] = arm_summary(a.patched)
    a.out.write_text(json.dumps(doc, indent=1, default=str) + "\n")
    print(f"wrote {a.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
