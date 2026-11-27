"""Campaign Phase 2 group C metrics (PXD057335, PXD066008, PXD079213), read from finished runs.

Writes one JSON of **aggregates only** (counts, shares with Wilson 95 % CIs, q-values, locus tags,
1-based coordinates); no peptide sequence and no DIA-NN row. Row numbers follow each dataset's
``ANSWER_KEY.md``. The standard checks reuse the Phase 1 helpers (``phase1_lib.py`` = a copy of
``validation/campaign/ecoli/phase1_metrics.py`` from branch validation/campaign-ecoli at 96f1a8be)
and the pilot's pure functions at the pinned commit. Run with ``PYTHONPATH=<pin>:<pin>/src:<this dir>``::

    p2_metrics.py --dataset mtb|hvol|sisl --a OUT_A --pep-a PEP_A --b OUT_B --pep-b PEP_B \
        --cmp CMP --bakta-json J --genome G --refseq-gff GFF --contaminants FASTA --out JSON

Rules fixed before reading the data (also in the answer keys):
- RefSeq names are carried to a Bakta CDS/sORF only when both share the stop codon (contig,
  strand, 3' end). For the archaea, TE3 and the operon check use these RefSeq names.
- Operon check: adjacent pairs (``pairs.tsv``) of genes both in the ribosomal-protein /
  ATP-synthase set; discordant = within-operon q <= 0.05 or an ``operon_flags`` row whose
  abundance pattern is not concordant/coabundant/none.
- Prophage/provirus regions: every geNomad region O < E and q <= 0.05 -> pass.
- Novel ORF match: >= 50 % same-strand overlap with a RefSeq CDS or pseudogene.
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
import sys
from collections import Counter
from pathlib import Path
from typing import Any

import pandas as pd
import phase1_lib as p1
from scipy.stats import fisher_exact
from validation.pilot import metrics as pm
from validation.pilot import te_real as tr
from validation.proxies.te_proxies import non_atg_nterm, read_coding_genes

from gippa.io.peptides import read_peptide_table

RIBO = re.compile(r"^(30S|50S) ribosomal protein")
ATP = re.compile(r"ATP synthase", re.IGNORECASE)
PEPPE = re.compile(r"\bPPE\b|\bPE\b|PE-PGRS|PE family|PPE family")
share, clean, tbl = p1.share, p1.clean, p1.tbl


def refseq_cds(gff: Path) -> pd.DataFrame:
    """Non-pseudo RefSeq CDS, one row per ID, with the stop key and names."""
    parts: dict[str, list[Any]] = {}
    with gzip.open(gff, "rt") as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "CDS":
                continue
            at = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            if at.get("pseudo") == "true":
                continue
            parts.setdefault(at.get("ID", f"{f[0]}:{f[3]}"), []).append((f, at))
    rows = []
    for segs in parts.values():
        f, at = segs[0]
        s = min(int(x[0][3]) for x in segs)
        e = max(int(x[0][4]) for x in segs)
        rows.append(
            {
                "contig": f[0],
                "strand": f[6],
                "start": s,
                "end": e,
                "stop": e if f[6] == "+" else s,
                "product": at.get("product", ""),
                "gene": at.get("gene", ""),
                "locus_tag": at.get("locus_tag", ""),
                "length_aa": (e - s + 1) // 3 - 1,
            }
        )
    return pd.DataFrame(rows)


def carry_names(genes: pd.DataFrame, ref: pd.DataFrame) -> pd.DataFrame:
    """genes.tsv rows with RefSeq product/gene/locus where the stop codon is shared."""
    g = genes.copy()
    s = pd.to_numeric(g["start"], errors="coerce")
    e = pd.to_numeric(g["end"], errors="coerce")
    g["_stop"] = e.where(g["strand"].astype(str) == "+", s)
    key = dict(
        zip(zip(ref["contig"], ref["strand"], ref["stop"], strict=False), ref.index, strict=False)
    )
    idx = [
        key.get((str(c), str(st), int(sp)) if pd.notna(sp) else None)
        for c, st, sp in zip(g["contig"], g["strand"], g["_stop"], strict=False)
    ]
    g["refseq_idx"] = idx
    g["refseq_product"] = [ref.at[i, "product"] if i is not None else "" for i in idx]
    g["refseq_gene"] = [ref.at[i, "gene"] if i is not None else "" for i in idx]
    g["refseq_locus"] = [ref.at[i, "locus_tag"] if i is not None else "" for i in idx]
    return g


def coding(g: pd.DataFrame) -> pd.DataFrame:
    return g[g["feature_type"].astype(str).isin(["cds", "sorf"])]


def det(g: pd.DataFrame) -> pd.Series:
    return g["detection_status"].astype(str) == "detected"


def gene_fdr(out: Path) -> dict[str, Any]:
    ds = tbl(out, "detection_summary")
    r = ds[
        (ds["section"].astype(str) == "tier")
        & (ds["group_type"].astype(str) == "q_threshold")
        & (ds["group"].astype(str) == "0.01")
    ]
    m = {str(x["metric"]): clean(x["value"]) for _, x in r.iterrows()}
    t = ds[(ds["section"].astype(str) == "tier") & (ds["group_type"].astype(str) == "tier")]
    tier = {f"{x['group']}:{x['metric']}": clean(x["value"]) for _, x in t.iterrows()}
    fdr = m.get("fdr")
    return {
        "fdr_at_q0.01": fdr,
        "n_detected_at_q0.01": m.get("n_detected"),
        "tier": tier,
        "verdict": "not testable" if fdr is None else ("pass" if fdr <= 0.01 else "fail"),
    }


def by_replicon(g: pd.DataFrame) -> dict[str, Any]:
    c = coding(g)
    return {str(k): share(int(det(x).sum()), len(x)) for k, x in c.groupby("contig")}


def operons_set(out: Path, tags: set[str], name: dict[str, str]) -> dict[str, Any]:
    pairs, fl = tbl(out, "pairs"), tbl(out, "operon_flags")
    res: dict[str, Any] = {"n_set_genes": len(tags)}
    disc: list[dict[str, Any]] = []
    if pairs is None:
        return res | {"verdict": "not testable", "reason": "missing:pairs.tsv"}
    inside = pairs[pairs["gene_a"].astype(str).isin(tags) & pairs["gene_b"].astype(str).isin(tags)]
    q = pd.to_numeric(inside["q_value"], errors="coerce")
    d = pd.to_numeric(inside["delta_log2"], errors="coerce").abs()
    res |= {
        "n_pairs": len(inside),
        "pair_class": dict(Counter(inside["pair_class"].astype(str))),
        "n_pairs_q_le_0.05": int((q <= 0.05).sum()),
        "min_q": clean(q.min()),
        "median_abs_delta_log2": clean(d.median()),
    }
    for _, r in inside[q <= 0.05].iterrows():
        disc.append(
            {
                "gene_a": name.get(str(r["gene_a"])),
                "gene_b": name.get(str(r["gene_b"])),
                "delta_log2": clean(r["delta_log2"]),
                "q_value": clean(r["q_value"]),
            }
        )
    if fl is not None:
        f = fl[fl["gene_a"].astype(str).isin(tags) & fl["gene_b"].astype(str).isin(tags)]
        res["operon_flags"] = [
            {
                "gene_a": name.get(str(r["gene_a"])),
                "gene_b": name.get(str(r["gene_b"])),
                "flags": clean(r["flags"]),
                "abundance_pattern": clean(r["abundance_pattern"]),
                "q_value": clean(r["q_value"]),
            }
            for _, r in f.iterrows()
        ]
        disc += [
            x
            for x in res["operon_flags"]
            if x["abundance_pattern"] not in ("concordant", "coabundant", "none", None)
        ]
    res["discordant"] = disc
    res["verdict"] = "not testable" if not len(inside) else ("pass" if not disc else "fail")
    return res


def regions(out: Path) -> dict[str, Any]:
    rs = tbl(out, "regions_stats")
    if rs is None:
        return {"verdict": "not testable", "reason": "missing:regions_stats.tsv"}
    types = dict(Counter(rs.loc[rs["subset"].astype(str) == "all", "region_type"].astype(str)))
    p = rs[
        rs["region_type"].astype(str).isin(["prophage", "provirus"])
        & (rs["subset"].astype(str) == "all")
    ]
    regs = []
    for rid, grp in p.groupby("region_id"):
        r0 = grp.iloc[0]
        oe = grp[grp["method"].astype(str) == "circular_shift_o_minus_e"]
        row = {
            "region_id": rid,
            "contig": clean(r0.get("contig")),
            "start": clean(r0["start"]),
            "end": clean(r0["end"]),
            "n_genes": clean(r0["n_genes"]),
            "observed": clean(r0["observed"]),
            "expected": clean(r0["expected"]),
            "o_over_e": clean(r0["o_over_e"]),
        }
        if len(oe):
            o = oe.iloc[0]
            row |= {
                "o_minus_e": clean(o["effect"]),
                "ci_low": clean(o["ci_low"]),
                "ci_high": clean(o["ci_high"]),
                "q_value": clean(o["q_value"]),
            }
        regs.append(row)
    lower = [r for r in regs if (r.get("o_minus_e") or 0) < 0]
    sig = [r for r in lower if r.get("q_value") is not None and r["q_value"] <= 0.05]
    isl = tbl(out, "islands")
    islands = (
        []
        if isl is None
        else [
            {
                k: clean(r[k])
                for k in (
                    "island_id",
                    "direction",
                    "contig",
                    "start",
                    "end",
                    "n_genes",
                    "detected_frac",
                    "effect",
                    "q_value",
                )
                if k in r
            }
            for _, r in isl.iterrows()
        ]
    )
    sigreg = rs[
        (rs["subset"].astype(str) == "all")
        & (pd.to_numeric(rs["q_value"], errors="coerce") <= 0.05)
    ]
    return {
        "region_types": types,
        "prophage_regions": regs,
        "n_lower": len(lower),
        "n_lower_q05": len(sig),
        "verdict": "not testable" if not regs else ("pass" if len(sig) == len(regs) else "fail"),
        "islands": islands,
        "significant_rows_q05": [
            {k: clean(r[k]) for k in ("region_type", "region_id", "method", "effect", "q_value")}
            for _, r in sigreg.iterrows()
        ],
    }


def named_genes(out: Path, g: pd.DataFrame, tags: dict[str, str]) -> dict[str, Any]:
    res = {}
    for name, tag in tags.items():
        res[name] = p1.gene_row(out, tag) if tag else None
    return res


def arm(out: Path) -> dict[str, Any]:
    s = p1.arm_summary(out)
    s["gene_fdr"] = gene_fdr(out)
    v = tbl(out, "variants")
    s["n_variant_rows"] = 0 if v is None else len(v)
    s["variant_confidence"] = {} if v is None else dict(Counter(v["confidence"].astype(str)))
    return s


def locus_by_stop(g: pd.DataFrame, contig: str, strand: str, stop: int) -> str | None:
    r = g[
        (g["contig"].astype(str) == contig)
        & (g["strand"].astype(str) == strand)
        & (g["_stop"] == stop)
    ]
    return str(r.iloc[0]["locus_tag"]) if len(r) else None


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    for k in ("a", "pep-a", "bakta-json", "genome", "refseq-gff", "contaminants", "out"):
        ap.add_argument(f"--{k}", type=Path, required=True)
    for k in ("b", "pep-b", "cmp", "patched"):
        ap.add_argument(f"--{k}", type=Path)
    ap.add_argument("--dataset", choices=["mtb", "hvol", "sisl"], required=True)
    a = ap.parse_args(argv)
    doc: dict[str, Any] = {"schema": "gippa-campaign-p2c/1", "dataset": a.dataset}
    ref_all = p1.refseq(a.refseq_gff)
    rc = refseq_cds(a.refseq_gff)
    A = a.a
    man = p1.manifest(A)
    cfg = pm._cfg(man, A / "run_manifest.json")
    gA = carry_names(tbl(A, "genes"), rc)
    PA = read_peptide_table(a.pep_a)[0]
    hitsA, unmA = tbl(A, "peptide_hits", "mapping"), tbl(A, "unmapped", "mapping")
    doc["A"] = arm(A)
    cod, circ = read_coding_genes(a.bakta_json, a.genome)
    te1, exc = pm.te1_real(PA, hitsA, unmA, cod, circ, non_atg_nterm(cod), cfg, q_crossings=0)
    doc["row_te1"] = te1 | {"n_exceptions": len(exc)}
    t3 = tr.te3(tbl(A, "genes"))
    t3.pop("_rows", None)
    doc["row_te3_gippa_set"] = t3
    gR = gA.copy()
    gR["product"] = gR["refseq_product"]
    gR["gene"] = gR["refseq_gene"]
    t3r = tr.te3(gR.drop(columns=["_stop", "refseq_idx"]))
    t3r.pop("_rows", None)
    doc["row_te3_refseq_set_pilot_rule"] = t3r
    # the answer keys' set: every RefSeq "(30S|50S) ribosomal protein" CDS carried by stop (the
    # pilot's strict product rule drops archaeal names such as "S19e"; reported above)
    Gc = tr._coding(gR.drop(columns=["_stop", "refseq_idx"]))
    mask = Gc["product"].astype(str).str.match(RIBO).to_numpy(dtype=bool)
    full = tr._te3_clauses(Gc, mask)
    full["status"] = (
        "not_evaluable"
        if not full["n_ribosomal"]
        else "pass"
        if full["detection"]["holds"] and full["abundance"]["holds"]
        else "fail"
    )
    full["not_detected"] = [
        {"locus_tag": r.locus_tag, "refseq_locus": r.refseq_locus, "product": r.product}
        for r in Gc[
            mask & (Gc["detection_status"].astype(str) != "detected").to_numpy()
        ].itertuples()
    ]
    doc["row_te3_refseq_set"] = full
    doc["row_te4"] = tr.te4(tbl(A, "genes"))
    heads = [ln[1:] for ln in a.contaminants.read_text().splitlines() if ln.startswith(">")]
    fam = {
        f: set(ids) | {"CONT_" + i for i in ids}
        for f, ids in tr.contaminant_families(heads).items()
    }
    t5 = tr.te5_arm(PA, hitsA, unmA, fam, cfg)
    t5.pop("_violations", None)
    doc["row_te5"] = t5
    doc["row_gene_fdr"] = doc["A"]["gene_fdr"]
    # operon set: Bakta names (mtb) or RefSeq names carried by stop (archaea)
    prod = gA["product"].astype(str) if a.dataset == "mtb" else gA["refseq_product"].astype(str)
    gname = gA["gene"].map(lambda v: clean(v) or "")
    use_name = (
        (gname.str.match(r"^(rp[lms][A-Z]|atp[A-I])$"))
        if a.dataset == "mtb"
        else pd.Series(False, index=gA.index)
    )
    sel = coding(gA)[
        (
            prod.loc[coding(gA).index].str.match(RIBO)
            | prod.loc[coding(gA).index].str.contains(ATP)
            | use_name.loc[coding(gA).index]
        )
    ]
    name = {
        str(t): (str(n) or str(p))
        for t, n, p in zip(gA["locus_tag"], gA["refseq_gene"], gA["refseq_product"], strict=False)
    }
    doc["row_operons"] = operons_set(A, set(sel["locus_tag"].astype(str)), name)
    doc["row_regions_A"] = regions(A)
    doc["by_replicon_A"] = by_replicon(gA)
    doc["refseq_mapped_by_stop"] = int(gA["refseq_idx"].notna().sum())
    if a.dataset == "mtb":
        tags = {
            n: locus_by_stop(gA, "NC_000962.3", "+", s)
            for n, s in (("esxA", 4352896), ("esxB", 4352576), ("dnaA", 1524))
        }
        doc["row_esx"] = named_genes(A, gA, tags)
        doc["row_esx"]["verdict"] = (
            "pass"
            if all(
                (doc["row_esx"][n] or {}).get("detection_status") == "detected"
                and (doc["row_esx"][n]["gene_q_value"] or 1) <= 0.01
                for n in ("esxA", "esxB")
            )
            else "fail"
        )
        c = coding(gA)
        pe = c["refseq_product"].astype(str).str.contains(PEPPE)
        d = det(c)
        k1, n1, k0, n0 = int((d & pe).sum()), int(pe.sum()), int((d & ~pe).sum()), int((~pe).sum())
        odds, p = fisher_exact([[k1, n1 - k1], [k0, n0 - k0]])
        doc["row_peppe"] = {
            "pe_ppe": share(k1, n1),
            "other": share(k0, n0),
            "fisher_p": float(p),
            "odds_ratio": float(odds),
            "verdict": "pass" if (k1 / max(n1, 1) < k0 / max(n0, 1) and p <= 0.05) else "fail",
        }
        rep = None
        for cand in sorted(A.glob("work/**/replicons.parquet")):
            rep = pd.read_parquet(cand).to_dict("records")
            break
        doc["row_topology"] = {
            "replicons_parquet": rep,
            "genes_crosses_origin": int(
                gA["crosses_origin"].map(lambda v: str(v).lower() == "true").sum()
            ),
        }
    if a.dataset == "hvol":
        c = coding(gA)
        doc["row_phv2_A"] = share(
            int(det(c[c["contig"].astype(str) == "NC_013965.1"]).sum()),
            int((c["contig"].astype(str) == "NC_013965.1").sum()),
        )
        doc["row_pyrE2_A"] = p1.gene_row(A, locus_by_stop(gA, "NC_013967.1", "+", 302281) or "")
    if a.dataset == "sisl":
        doc["row_negctl_A"] = {
            n: p1.gene_row(A, locus_by_stop(gA, "NC_017276.1", "-", s) or "")
            for n, s in (("pyrF", 1307768), ("pyrE", 1308423), ("lacS", 2067746))
        }
    if a.b:
        B = a.b
        gB = carry_names(tbl(B, "genes"), rc)
        doc["B"] = arm(B)
        doc["row_novel"] = p1.novel_vs_refseq(B, ref_all)
        doc["q_floor_warnings"] = [
            x for x in (doc["B"]["warnings"] or []) if "class_q_floor" in str(x)
        ]
        dA, dB = p1.detected(A), p1.detected(B)
        doc["row_route_b"] = {
            "A_detected": len(dA),
            "B_detected": len(dB),
            "A_kept_by_B": share(len(dA & dB), len(dA)),
            "B_only": len(dB - dA),
        }
        doc["by_replicon_B"] = by_replicon(gB)
        doc["row_regions_B"] = regions(B)
        n = tbl(B, "novel_orfs")
        # RefSeq CDS without a Bakta stop: is there a route-B candidate overlapping them same-strand?
        bstop = set(
            zip(gA["contig"].astype(str), gA["strand"].astype(str), gA["_stop"], strict=False)
        )
        miss = rc[
            [
                (c_, s_, sp) not in bstop
                for c_, s_, sp in zip(rc["contig"], rc["strand"], rc["stop"], strict=False)
            ]
        ]
        rec = []
        hB = tbl(B, "peptide_hits", "mapping")
        hs, he = (
            pd.to_numeric(hB["start"], errors="coerce"),
            pd.to_numeric(hB["end"], errors="coerce"),
        )
        for _, r in miss.iterrows():
            cand = []
            if n is not None and len(n):
                lo = pd.concat(
                    [
                        pd.to_numeric(n["start_longest"], errors="coerce"),
                        pd.to_numeric(n["end"], errors="coerce"),
                    ],
                    axis=1,
                )
                ov = n[
                    (n["strand"].astype(str) == r["strand"])
                    & (n["contig"].astype(str) == r["contig"] if "contig" in n else True)
                    & (lo.min(axis=1) <= r["end"])
                    & (lo.max(axis=1) >= r["start"])
                ]
                cand = ov["orf_id"].astype(str).tolist()
            hh = hB[
                (hB["contig"].astype(str) == r["contig"])
                & (hB["strand"].astype(str) == r["strand"])
                & (hs >= r["start"])
                & (he <= r["end"])
                & (hB["route"].astype(str) == "exact")
            ]
            rec.append(
                {
                    "refseq_locus": r["locus_tag"],
                    "contig": r["contig"],
                    "start": int(r["start"]),
                    "end": int(r["end"]),
                    "strand": r["strand"],
                    "length_aa": int(r["length_aa"]),
                    "product": r["product"],
                    "route_b_candidates": cand,
                    "route_b_exact_peptides": int(hh["peptide"].nunique()),
                }
            )
        doc["refseq_cds_missed_by_bakta"] = {
            "n": len(rec),
            "n_le_70aa": sum(x["length_aa"] <= 70 for x in rec),
            "with_route_b_peptides": sum(x["route_b_exact_peptides"] > 0 for x in rec),
            "le_70aa_with_route_b_peptides": sum(
                x["route_b_exact_peptides"] > 0 and x["length_aa"] <= 70 for x in rec
            ),
            "with_candidate": sum(bool(x["route_b_candidates"]) for x in rec),
            "rows_with_evidence": [
                x for x in rec if x["route_b_exact_peptides"] > 0 or x["route_b_candidates"]
            ],
        }
        if n is not None and len(n):
            keep = (
                "orf_id",
                "class",
                "contig",
                "strand",
                "start_longest",
                "end",
                "length_aa_longest",
                "n_peptides",
                "best_peptide_q",
                "class_fdr",
                "called",
                "overlapping_gene",
                "flags",
            )
            top = n[pd.to_numeric(n["n_peptides"], errors="coerce") >= 2]
            doc["candidates_ge2_peptides"] = [
                {k: clean(r[k]) for k in keep if k in r} for _, r in top.iterrows()
            ]
        if a.dataset == "hvol":
            c = coding(gB)
            doc["row_phv2_B"] = share(
                int(det(c[c["contig"].astype(str) == "NC_013965.1"]).sum()),
                int((c["contig"].astype(str) == "NC_013965.1").sum()),
            )
            doc["row_pyrE2_B"] = p1.gene_row(B, locus_by_stop(gB, "NC_013967.1", "+", 302281) or "")
        if a.dataset == "sisl":
            doc["row_negctl_B"] = {
                n_: p1.gene_row(B, locus_by_stop(gB, "NC_017276.1", "-", s) or "")
                for n_, s in (("pyrF", 1307768), ("pyrE", 1308423), ("lacS", 2067746))
            }
    if a.cmp:
        ab = a.cmp / "tables" / "compare_abundance.tsv"
        st = a.cmp / "tables" / "compare_stability.tsv"
        doc["compare"] = {
            "abundance": pd.read_csv(ab, sep="\t").to_dict("records") if ab.is_file() else None,
            "stability": pd.read_csv(st, sep="\t").to_dict("records") if st.is_file() else None,
        }
    if a.patched:
        doc["round2"] = arm(a.patched)
    a.out.write_text(json.dumps(doc, indent=1, default=str) + "\n")
    print(f"wrote {a.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
