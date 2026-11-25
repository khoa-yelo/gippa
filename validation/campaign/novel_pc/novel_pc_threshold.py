#!/usr/bin/env python3
"""Positive-control q threshold for route-B novel-ORF candidates (campaign analysis, 2026-11-24).

Khoa's request: use the class-FDR q of the prfB N-terminal extension (and of other candidates
that a known true gene independently confirms) as an extra, clearly labelled reporting tier for
novel ORFs in the other route-B runs.  This changes no GIPPA call and no alpha (R3): GIPPA's calls
at alpha 0.01 stay as they are.  Reads ``tables/novel_orfs.tsv`` and ``mapping/peptide_hits.tsv``
of each route-B run and the RefSeq GFF of its genome.

The q used is the ``class_fdr`` column: ``gippa.novel.evidence._call`` calls a candidate with a
numeric class FDR when ``class_fdr <= alpha`` (and the class rule holds).

Output (aggregates only; never peptide sequences or DIA-NN rows):
    novel_pc_threshold.json   controls, thresholds, per-run tallies, passer lists
    passers.tsv               one row per (threshold, passing candidate)

Usage:
    python novel_pc_threshold.py --root $GIPPA_ROOT/campaign --outdir DIR
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import pandas as pd

# ------------------------------------------------------------------ runs (route B, six-frame)
# (dataset label, accession, run label, run dir rel. to root, genome dir rel. to root,
#  contig alias novel_orfs contig -> GFF seqid)
RUNS: list[tuple[str, str, str, str, str, dict[str, str]]] = [
    ("E. coli K-12 (Phase 1 E2)", "ecoli", "ecoli_e2", "ecoli_k12/runs/out_e2_fix", "ecoli_k12/ref",
     {"U00096.3": "NC_000913.3"}),
    *[("T. pallidum", "PXD065022", s, f"PXD065022/{s}/out_b_fix2", f"PXD065022/genome/{s}", {})
      for s in ("tpa_dal1", "tpa_grady", "tpa_haitib", "tpa_madras", "tpa_mexicoa", "tpa_phil1",
                "tpa_ss14")],
    ("C. crescentus", "PXD069004", "ccre", "PXD069004/ccre_na1000/out_b", "PXD069004/genome/ccre_na1000", {}),
    ("B. thetaiotaomicron", "PXD055395", "bthe", "PXD055395/bthe_vpi5482/out_b",
     "PXD055395/genome/bthe_vpi5482", {}),
    ("Synechococcus CC9311", "PXD072227", "syn", "PXD072227/out_b", "PXD072227/genome/syn_cc9311", {}),
    ("V. natriegens", "PXD049476", "vnat", "PXD049476/out_b", "PXD049476/genome/vnat_atcc14048", {}),
    ("M. xanthus", "PXD060688", "mxan", "PXD060688/out_b", "PXD060688/genome/mxan_dk1622", {}),
    ("C. difficile", "PXD065346", "cdif", "PXD065346/out_b", "PXD065346/genome/cdif_630derm", {}),
    ("M. tuberculosis", "PXD057335", "mtb", "PXD057335/out_b", "PXD057335/genome/mtb_h37rv_circ", {}),
    ("H. volcanii", "PXD066008", "hvo", "PXD066008/out_b", "PXD066008/genome/hvol_ds2", {}),
    ("S. islandicus", "PXD079213", "sisl", "PXD079213/out_b", "PXD079213/genome/sisl_rey15a", {}),
]

# ------------------------------------------------------------------ controls (pre-registered)
# Written down before any passer list was looked at.  Each control is a RefSeq gene whose protein
# is independently established; a route-B candidate is that control when it meets the rule.
#   kind "prfB":      RefSeq prfB (programmed +1 frameshift).  Candidate on prfB's strand whose ORF
#                     span overlaps the prfB gene.  Class nterm_extension = primary (the in-frame
#                     part upstream of Bakta's start, ending at the frameshift-site stop); any
#                     other class over the 25-codon leader = secondary ("prfB leader").
#   kind "infC":      Bakta's late infC start: nterm_extension candidate overlapping RefSeq infC.
#   kind "prdB":      PrdB selenoprotein read-through past the Sec (UGA) codon: candidate
#                     overlapping RefSeq prdB, any class.
#   kind "missed":    RefSeq gene Bakta did not call, with peptides: candidate span covers >= 50 %
#                     of the RefSeq CDS, same strand.
# Datasets are listed where the control is expected; prfB is searched in every run.
CONTROLS: list[dict[str, Any]] = [
    {"name": "prfB N-terminal extension", "kind": "prfB", "match": {"gene": "prfB"}, "runs": "all",
     "reason": "PrfB (RF2) is essential; its N-terminus is made by a programmed +1 frameshift that "
               "Bakta cannot model, so peptides upstream of Bakta's start are true by biology"},
    {"name": "infC late start", "kind": "infC", "match": {"gene": "infC"}, "runs": ["ecoli_e2", "cdif"],
     "reason": "IF3 starts at the AUU codon upstream of Bakta's later ATG/GTG start; peptides there "
               "are true N-terminal sequence"},
    {"name": "PrdB read-through (Sec)", "kind": "prdB", "match": {"gene": "prdB"}, "runs": ["cdif"],
     "reason": "PrdB is a selenoprotein; translation continues through the UGA (Sec) codon"},
    {"name": "rpmB (L28)", "kind": "missed", "match": {"product": "50S ribosomal protein L28"},
     "runs": ["tpa_dal1", "tpa_grady", "tpa_haitib", "tpa_madras", "tpa_mexicoa", "tpa_phil1", "tpa_ss14"],
     "reason": "RefSeq ribosomal protein L28 not called by Bakta"},
    {"name": "L36 Rv3461c", "kind": "missed", "match": {"locus_tag": "Rv3461c"}, "runs": ["mtb"],
     "reason": "RefSeq ribosomal protein L36 (MS-supported in RefSeq) not called by Bakta"},
    {"name": "L34 Rv3924c", "kind": "missed", "match": {"locus_tag": "Rv3924c"}, "runs": ["mtb"],
     "reason": "RefSeq ribosomal protein L34 not called by Bakta"},
    {"name": "L37ae", "kind": "missed", "match": {"product": "50S ribosomal protein L37ae"}, "runs": ["sisl"],
     "reason": "RefSeq ribosomal protein L37ae not called by Bakta"},
    {"name": "L40e", "kind": "missed", "match": {"product": "50S ribosomal protein L40e"}, "runs": ["sisl"],
     "reason": "RefSeq ribosomal protein L40e not called by Bakta"},
    {"name": "cdzC", "kind": "missed", "match": {"locus_tag": "CCNA_03933"}, "runs": ["ccre"],
     "reason": "RefSeq bacteriocin cdzC (LC-MS / Ribo-seq evidence in RefSeq) not called by Bakta"},
    {"name": "cdzD", "kind": "missed", "match": {"locus_tag": "CCNA_03932"}, "runs": ["ccre"],
     "reason": "RefSeq bacteriocin cdzD not called by Bakta"},
    {"name": "cdzI", "kind": "missed", "match": {"locus_tag": "CCNA_03931"}, "runs": ["ccre"],
     "reason": "RefSeq bacteriocin immunity cdzI not called by Bakta"},
]

TWO_PEPTIDE_CLASSES = {"cds_outframe", "antisense"}  # gippa.novel.evidence (read-only reference)


def wilson(k: int, n: int, z: float = 1.959963984540054) -> list[float | None]:
    if n == 0:
        return [None, None, None]
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return [p, max(0.0, c - h), min(1.0, c + h)]


def read_gff(gdir: Path) -> pd.DataFrame:
    gff = next(gdir.glob("*_genomic.gff.gz"))
    rows = []
    genes: dict[str, dict[str, str]] = {}
    with gzip.open(gff, "rt") as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            attr = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            if f[2] in ("gene", "pseudogene"):
                genes[attr.get("ID", "")] = attr
            if f[2] not in ("CDS", "pseudogene"):
                continue
            g = genes.get(attr.get("Parent", ""), {})
            olt = (attr.get("old_locus_tag") or g.get("old_locus_tag") or "").replace("%2C", ",")
            rows.append({"seqid": f[0], "type": f[2], "start": int(f[3]), "end": int(f[4]), "strand": f[6],
                         "id": attr.get("ID"), "gene": attr.get("gene") or g.get("gene") or g.get("Name"),
                         "locus_tag": attr.get("locus_tag") or g.get("locus_tag"),
                         "old_locus_tag": olt, "product": attr.get("product", ""),
                         "pseudo": attr.get("pseudo") == "true" or f[2] == "pseudogene"})
    return pd.DataFrame(rows)


def num(v: Any) -> float | None:
    try:
        x = float(v)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(x) else x


def load_novel(run: Path) -> pd.DataFrame:
    n = pd.read_csv(run / "tables" / "novel_orfs.tsv", sep="\t", dtype=str, keep_default_na=False)
    n = n.drop(columns=["peptides"])  # never carry peptide sequences
    n["q"] = n["class_fdr"].map(num)
    n["npep"] = n["n_peptides"].astype(int)
    n["has_start"] = n["start_longest"].map(lambda v: num(v) is not None)
    lo, hi = [], []
    for _, r in n.iterrows():
        e = int(r["end"])
        s = num(r["start_longest"])
        if s is not None:
            a, b = sorted((int(s), e))
        else:  # no start codon: ORF span from the stop and the ORF length (approximate)
            span = 3 * (int(r["length_aa_longest"]) + 1)
            a, b = (e - span + 1, e) if r["strand"] == "+" else (e, e + span - 1)
        lo.append(a)
        hi.append(b)
    n["lo"], n["hi"] = lo, hi
    return n


def peptide_spans(run: Path) -> dict[str, tuple[int, int]]:
    h = pd.read_csv(run / "mapping" / "peptide_hits.tsv", sep="\t", dtype=str, keep_default_na=False,
                    usecols=["start", "end", "novel_candidate_id"])
    h = h[h["novel_candidate_id"].astype(str).str.startswith("NOVEL_")]
    out: dict[str, tuple[int, int]] = {}
    for cid, g in h.groupby("novel_candidate_id"):
        out[str(cid)] = (int(pd.to_numeric(g["start"]).min()), int(pd.to_numeric(g["end"]).max()))
    return out


def best_same_strand(r: pd.Series, ref: pd.DataFrame, seqid: str, lo: int, hi: int) -> tuple[int, Any]:
    same = ref[(ref["seqid"] == seqid) & (ref["strand"] == r["strand"]) & (ref["start"] <= hi) & (ref["end"] >= lo)]
    best = (0, None)
    for _, c in same.iterrows():
        ov = min(hi, c["end"]) - max(lo, c["start"]) + 1
        if ov > best[0]:
            best = (ov, c)
    return best


def refseq_match(r: pd.Series, ref: pd.DataFrame, seqid: str) -> dict[str, Any]:
    """Report rule: >= 50 % of the candidate (start_longest..end) on one same-strand RefSeq CDS or
    pseudogene; no start codon = no match.  Also the span variant (no-start ORFs by their span)."""
    lo, hi = int(r["lo"]), int(r["hi"])
    ov, c = best_same_strand(r, ref, seqid, lo, hi)
    frac = ov / (hi - lo + 1)
    return {"strict": bool(r["has_start"] and c is not None and frac >= 0.5),
            "span": bool(c is not None and frac >= 0.5),
            "refseq_gene": None if c is None else (c["gene"] or c["locus_tag"]),
            "refseq_locus": None if c is None else c["locus_tag"], "overlap_frac": round(frac, 3)}


def find_controls(label: str, n: pd.DataFrame, ref: pd.DataFrame, alias: dict[str, str]) -> list[dict[str, Any]]:
    out = []
    for ctl in CONTROLS:
        if ctl["runs"] != "all" and label not in ctl["runs"]:
            continue
        key, val = next(iter(ctl["match"].items()))
        feats = ref[(ref[key] == val) | ((key == "locus_tag") & ref["old_locus_tag"].str.split(",").map(lambda x: val in x))]
        if feats.empty:
            continue
        for _, r in n.iterrows():
            seqid = alias.get(r["contig"], r["contig"])
            f = feats[(feats["seqid"] == seqid) & (feats["strand"] == r["strand"])]
            if f.empty:
                continue
            glo, ghi = int(f["start"].min()), int(f["end"].max())
            lo, hi = int(r["lo"]), int(r["hi"])
            if hi < glo or lo > ghi:
                continue
            cls = r["class"]
            if ctl["kind"] == "prfB":
                role = "primary" if cls == "nterm_extension" else "secondary (prfB leader)"
            elif ctl["kind"] == "infC":
                if cls != "nterm_extension":
                    continue
                role = "secondary"
            elif ctl["kind"] == "prdB":  # stop = prdB stop: past the Sec; else the part up to the Sec
                past = int(r["end"]) == (glo if r["strand"] == "-" else ghi)
                role = "secondary" if past else "secondary (prdB N-terminal part, ends at the Sec codon)"
            else:  # missed: candidate covers >= 50 % of the RefSeq CDS
                cov = max(0, min(hi, ghi) - max(lo, glo) + 1) / (ghi - glo + 1)
                if cov < 0.5:
                    continue
                role = "secondary"
            out.append({"control": ctl["name"], "kind": ctl["kind"], "role": role, "run": label,
                        "orf_id": r["orf_id"], "class": cls, "q": r["q"], "n_peptides": int(r["npep"]),
                        "has_start": bool(r["has_start"]), "flags": r["flags"],
                        "refseq_feature": f"{val} {glo}-{ghi} ({r['strand']})"})
    return out


def class_rule(r: pd.Series) -> str:
    """The visible part of GIPPA's class rule (rule 4); start_upstream_nt is not in the table."""
    c = r["class"]
    if c in TWO_PEPTIDE_CLASSES:
        return "met" if r["npep"] >= 2 else "not met (needs >= 2 peptides)"
    if c in ("nterm_extension", "intergenic"):
        ok = r["has_start"] and str(r["rbs_supported"]).lower() == "true"
        return "met" if ok else "not met (needs start codon + RBS)"
    return "n/a"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    a = ap.parse_args()
    a.outdir.mkdir(parents=True, exist_ok=True)

    data: dict[str, tuple[pd.DataFrame, pd.DataFrame, dict[str, str], tuple[str, str]]] = {}
    controls: list[dict[str, Any]] = []
    for ds, acc, label, rd, gd, alias in RUNS:
        n = load_novel(a.root / rd)
        spans = peptide_spans(a.root / rd)
        n["pep_lo"] = n["orf_id"].map(lambda i: spans.get(i, (None, None))[0])
        n["pep_hi"] = n["orf_id"].map(lambda i: spans.get(i, (None, None))[1])
        ref = read_gff(a.root / gd)
        data[label] = (n, ref, alias, (ds, acc))
        controls += find_controls(label, n, ref, alias)

    ctl_ids = {(c["run"], c["orf_id"]) for c in controls}
    prfb_q = [c["q"] for c in controls if c["kind"] == "prfB" and c["role"] == "primary" and c["q"] is not None]
    all_q = [c["q"] for c in controls if c["q"] is not None]
    tau = {"tau_prfB": max(prfb_q), "tau_all": max(all_q)}
    per_class: dict[str, float] = {}
    for c in controls:
        if c["q"] is not None:
            per_class[c["class"]] = max(per_class.get(c["class"], 0.0), c["q"])
    tau.update({f"tau_class:{k}": v for k, v in sorted(per_class.items())})
    # post hoc sensitivity (not pre-registered): tau_all without the T. pallidum rpmB controls,
    # whose classes have 17-41 candidates and q in {0.5, 1.0}
    tau["tau_all_excl_tpa_rpmB"] = max(c["q"] for c in controls if c["q"] is not None
                                       and not (c["control"] == "rpmB (L28)"))
    tiers: dict[str, dict[str, float]] = {
        "tau_prfB": defaultdict(lambda: tau["tau_prfB"]),
        "tau_all": defaultdict(lambda: tau["tau_all"]),
        "tau_class": defaultdict(lambda: -1.0, per_class),  # classes without a control: none pass
        "tau_all_excl_tpa_rpmB": defaultdict(lambda: tau["tau_all_excl_tpa_rpmB"]),
    }

    runs_out: dict[str, Any] = {}
    passer_rows = []
    for tname, tmap in tiers.items():
        for label, (n, ref, alias, (ds, acc)) in data.items():
            p = n[[r["q"] is not None and r["q"] <= tmap[r["class"]] + 1e-12 for _, r in n.iterrows()]]
            rec: dict[str, Any] = {"dataset": ds, "accession": acc, "n_candidates": int(len(n)),
                                   "n_pass": int(len(p)), "pass_by_class": dict(Counter(p["class"])),
                                   "class_q_max_of_passers": {c: float(g["q"].max()) for c, g in p.groupby("class")},
                                   "candidates_by_class": dict(Counter(n["class"]))}
            exp_false = sum(len(g) * float(g["q"].max()) for _, g in p.groupby("class"))
            rec["gippa_class_fdr_expected_false"] = exp_false
            ks = ksp = nc = 0
            lst = []
            for _, r in p.iterrows():
                is_ctl = (label, r["orf_id"]) in ctl_ids
                m = refseq_match(r, ref, alias.get(r["contig"], r["contig"]))
                row = {"threshold": tname, "tau": tmap[r["class"]], "dataset": ds, "accession": acc, "run": label,
                       "orf_id": r["orf_id"], "class": r["class"], "contig": r["contig"], "strand": r["strand"],
                       "lo": int(r["lo"]), "hi": int(r["hi"]), "has_start": bool(r["has_start"]),
                       "aa": int(r["length_aa_longest"]), "n_peptides": int(r["npep"]), "q": float(r["q"]),
                       "control": is_ctl, "refseq_match": m["strict"], "refseq_match_span": m["span"],
                       "refseq_gene": m["refseq_gene"], "refseq_locus": m["refseq_locus"],
                       "overlap_frac": m["overlap_frac"],
                       "needs_2_peptides": r["class"] in TWO_PEPTIDE_CLASSES, "class_rule": class_rule(r),
                       "flags": r["flags"]}
                passer_rows.append(row)
                if not is_ctl:
                    nc += 1
                    ks += m["strict"]
                    ksp += m["span"]
                    lst.append(row)
            rec.update({"n_controls_passing": int(len(p)) - nc, "n_noncontrol_passing": nc,
                        "noncontrol_refseq_match": [ks, nc, *wilson(ks, nc)],
                        "noncontrol_refseq_match_span": [ksp, nc, *wilson(ksp, nc)]})
            runs_out.setdefault(tname, {})[label] = rec
        tot = runs_out[tname]
        K = sum(r["noncontrol_refseq_match"][0] for r in tot.values())
        Ks = sum(r["noncontrol_refseq_match_span"][0] for r in tot.values())
        N = sum(r["n_noncontrol_passing"] for r in tot.values())
        byc: Counter[str] = Counter()
        for r in tot.values():
            byc.update(r["pass_by_class"])
        runs_out[tname + "_total"] = {"n_pass": sum(r["n_pass"] for r in tot.values()), "pass_by_class": dict(byc),
                                      "n_controls_passing": sum(r["n_controls_passing"] for r in tot.values()),
                                      "n_noncontrol_passing": N,
                                      "noncontrol_refseq_match": [K, N, *wilson(K, N)],
                                      "noncontrol_refseq_match_span": [Ks, N, *wilson(Ks, N)],
                                      "gippa_class_fdr_expected_false": sum(r["gippa_class_fdr_expected_false"] for r in tot.values()),
                                      "n_candidates": sum(r["n_candidates"] for r in tot.values())}

    doc = {"q_column": "class_fdr", "controls_spec": CONTROLS, "controls_found": controls, "tau": tau,
           "results": runs_out}
    (a.outdir / "novel_pc_threshold.json").write_text(json.dumps(doc, indent=1, default=str))
    pd.DataFrame(passer_rows).to_csv(a.outdir / "passers.tsv", sep="\t", index=False)
    print(json.dumps({"tau": tau, "n_controls": len(controls)}, indent=1))


if __name__ == "__main__":
    main()
