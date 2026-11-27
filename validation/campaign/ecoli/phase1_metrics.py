"""Phase 1 (E. coli K-12 MG1655) acceptance metrics, read from finished GIPPA runs.

Reads the run directories of arms E1-E4 and writes one JSON of **aggregates only** (counts,
shares with Wilson 95 % CIs, q-values, locus tags, 1-based coordinates): no peptide sequence and
no DIA-NN row leaves ``$GIPPA_ROOT``. The row numbers follow ``ANSWER_KEY.md``. TE1, TE3, TE4,
TE5 and TE2 reuse the pilot's pure functions (``validation/pilot/metrics.py``, ``te_real.py``) at
the pinned commit, so the campaign judges exactly as the pilot did. Run with the pinned tree on
``PYTHONPATH`` (``<pin>:<pin>/src``)::

    phase1_metrics.py --e1 OUT --pep-e1 PEP [--e2 OUT --pep-e2 PEP --cmp CMP] \
        [--e3 LABEL=OUT=PEP ...] [--e4 OUT] [--te7 OUT_REP] --bakta-json J --genome G \
        --refseq-gff GFF --contaminants FASTA --out metrics.json

Verdict rules fixed before the data were read (stated again in REPORT.md):
- row 5: every geNomad prophage region of E1 has observed < expected detection
  (``circular_shift_o_minus_e``) and q <= 0.05 -> pass; any region with O >= E -> fail; all
  lower but some q > 0.05 -> fail (direction holds, significance does not).
- row 6: no gene pair inside a ribosomal-protein or ATP-synthase transcription unit carries an
  abundance-discordance flag in ``operon_flags.tsv`` -> pass.
- row 9: an ORF "matches" when a RefSeq CDS or pseudogene on the same strand overlaps >= 50 % of
  its length (share and Wilson CI over the called ORFs).
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

LOCI = {  # Bakta 1.12.2 locus tags on U00096.3 (ANSWER_KEY.md loci table)
    "crl": "KOIAIB_00263",
    "prfB": "KOIAIB_03083",
    "fdnG": "KOIAIB_01574",
    "fdoG": "KOIAIB_04147",
    "fdhF": "KOIAIB_04360",
    "mntS": "KOIAIB_00882",
    "lon": "KOIAIB_00480",
    "clpP": "KOIAIB_00477",
}
SEC = {"fdnG": (1547986, 1547988), "fdoG": (4085235, 4085237), "fdhF": (4298947, 4298949)}
CRL_WINDOWS = {  # 1-based closed, + strand (ANSWER_KEY.md)
    "fragment1": (257829, 257899),
    "nterm_extension": (258676, 258733),
    "bakta_cds": (258734, 259006),
}
PRFB_WINDOWS = {  # - strand
    "leader_other_frame": (3036208, 3036282),
    "nterm_extension": (3036066, 3036206),
    "bakta_cds": (3035184, 3036065),
}
PROPHAGE_INT = {
    "intD": "DLP12",
    "intE": "e14",
    "intR": "Rac",
    "intQ": "Qin",
    "intA": "CP4-57",
    "intF": "CP4-6",
}
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
    tier = Counter(
        cod.loc[cod["detection_status"].astype(str) == "detected", "detection_tier"].astype(str)
    )
    q = pd.to_numeric(cod["gene_q_value"], errors="coerce")
    fdr = tbl(out, "mapping_fdr", "mapping")
    fdr_rows = []
    if fdr is not None:
        f = fdr[(fdr["length_class"].astype(str) == "all")]
        for r in f.itertuples(index=False):
            fdr_rows.append(
                {
                    k: clean(getattr(r, k))
                    for k in (
                        "section",
                        "n_diff",
                        "targets",
                        "decoys",
                        "n_accepted",
                        "precision",
                        "ci_low",
                        "ci_high",
                        "flags",
                    )
                }
            )
    return {
        "commit": man.get("git_commit"),
        "git_dirty": man.get("git_dirty"),
        "exit_code": man.get("exit_code"),
        "status": man.get("status"),
        "detection_tier": man.get("detection_tier"),
        "started_at": man.get("started_at"),
        "finished_at": man.get("finished_at"),
        "warnings": man.get("warnings"),
        "layers_missing": man.get("layers_missing"),
        "modules_off": man.get("modules_off"),
        "seeds": man.get("seeds"),
        "inputs": {i.get("role"): i.get("sha256") for i in man.get("inputs", [])},
        "n_coding_genes": len(cod),
        "detection_status": dict(st),
        "detected_by_tier": dict(tier),
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
    return {
        k: clean(r[k])
        for k in (
            "locus_tag",
            "gene",
            "start",
            "end",
            "strand",
            "detection_status",
            "detection_tier",
            "gene_q_value",
            "n_peptides_unique",
            "n_peptides_shared",
            "coverage_frac",
            "abundance_log2",
            "flags",
        )
    }


def hits_in(out: Path, contig: str, lo: int, hi: int, strand: str) -> pd.DataFrame:
    h = tbl(out, "peptide_hits", "mapping")
    s = pd.to_numeric(h["start"], errors="coerce")
    e = pd.to_numeric(h["end"], errors="coerce")
    return h[
        (h["contig"].astype(str) == contig)
        & (h["strand"].astype(str) == strand)
        & (s <= hi)
        & (e >= lo)
    ]


def window_counts(out: Path, windows: dict[str, tuple[int, int]], strand: str) -> dict[str, Any]:
    res = {}
    for name, (lo, hi) in windows.items():
        h = hits_in(out, "U00096.3", lo, hi, strand)
        ex = h[h["route"].astype(str) == "exact"]
        res[name] = {
            "window": [lo, hi],
            "strand": strand,
            "distinct_peptides_exact": int(ex["peptide"].nunique()),
            "distinct_peptides_any_route": int(h["peptide"].nunique()),
            "locus_class": dict(Counter(ex["locus_class"].astype(str))),
            "best_q": clean(pd.to_numeric(ex["best_q_input"], errors="coerce").min())
            if len(ex)
            else None,
        }
    return res


def novel_in(out: Path, lo: int, hi: int, strand: str) -> list[dict[str, Any]]:
    n = tbl(out, "novel_orfs")
    if n is None:
        return []
    s = pd.to_numeric(n["start_longest"], errors="coerce")
    e = pd.to_numeric(n["end"], errors="coerce")
    lo_c = pd.concat([s, e], axis=1).min(axis=1)
    hi_c = pd.concat([s, e], axis=1).max(axis=1)
    sel = n[
        ((n["strand"].astype(str) == strand) & (lo_c <= hi) & (hi_c >= lo))
        .fillna(False)
        .astype(bool)
    ]
    keep = (
        "orf_id",
        "class",
        "strand",
        "start_longest",
        "start_best_ostir",
        "end",
        "length_aa_longest",
        "n_peptides",
        "best_peptide_q",
        "class_fdr",
        "called",
        "overlapping_gene",
        "flags",
    )
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
            rows.append(
                {
                    "type": f[2],
                    "start": int(f[3]),
                    "end": int(f[4]),
                    "strand": f[6],
                    "name": attr.get("gene") or attr.get("Name"),
                    "locus_tag": attr.get("locus_tag"),
                    "pseudo": attr.get("pseudo") == "true"
                    or f[2] == "pseudogene"
                    or attr.get("gene_biotype") == "pseudogene",
                }
            )
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
        out_rows.append(
            {
                "orf_id": clean(r["orf_id"]),
                "class": clean(r["class"]),
                "start": a,
                "end": b,
                "strand": r["strand"],
                "n_peptides": clean(r["n_peptides"]),
                "class_fdr": clean(r["class_fdr"]),
                "matches_refseq": bool(match),
                "refseq_name": None if best is None else best[1]["name"],
                "refseq_locus": None if best is None else best[1]["locus_tag"],
                "refseq_pseudo": None if best is None else bool(best[1]["pseudo"]),
                "overlap_frac": None if best is None else best[0] / L,
            }
        )
    cls = Counter(n["class"].astype(str))
    ccls = Counter(called["class"].astype(str))
    fdr = {
        c: clean(
            pd.to_numeric(n.loc[n["class"].astype(str) == c, "class_fdr"], errors="coerce").min()
        )
        for c in cls
    }
    return {
        "n_candidates": len(n),
        "candidates_by_class": dict(cls),
        "n_called": len(called),
        "called_by_class": dict(ccls),
        "min_class_fdr_by_class": fdr,
        "share_matching_refseq": share(k, len(called)),
        "called": out_rows,
        "exploratory_all_candidates": res_all,
    }


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
        best = max(
            ((min(b, c["end"]) - max(a, c["start"]) + 1, c) for _, c in same.iterrows()),
            key=lambda x: x[0],
            default=None,
        )
        rows.append(
            {
                "class": str(r["class"]),
                "n_peptides": int(clean(r["n_peptides"]) or 0),
                "match": best is not None and st is not None and best[0] >= 0.5 * L,
                "pseudo": bool(best[1]["pseudo"]) if best is not None else False,
            }
        )
    out: dict[str, Any] = {}
    for key, sub in (("all", rows), ("n_peptides_ge_2", [x for x in rows if x["n_peptides"] >= 2])):
        out[key] = share(sum(x["match"] for x in sub), len(sub))
        out[key + "_by_class"] = {
            c: share(
                sum(x["match"] for x in sub if x["class"] == c), sum(x["class"] == c for x in sub)
            )
            for c in sorted({x["class"] for x in sub})
        }
        out[key + "_pseudogene_matches"] = sum(x["match"] and x["pseudo"] for x in sub)
    return out


def prophages(out: Path, ref: pd.DataFrame) -> dict[str, Any]:
    rs = tbl(out, "regions_stats")
    if rs is None:
        return {"reason": "missing:regions_stats.tsv"}
    p = rs[(rs["region_type"].astype(str) == "prophage") & (rs["subset"].astype(str) == "all")]
    ints = ref[(ref["type"] == "gene") & ref["name"].isin(list(PROPHAGE_INT))]
    regs = []
    for rid, grp in p.groupby("region_id"):
        oe = grp[grp["method"].astype(str) == "circular_shift_o_minus_e"]
        hl = grp[grp["method"].astype(str) == "circular_shift_hodges_lehmann"]
        r0 = grp.iloc[0]
        a, b = int(r0["start"]), int(r0["end"])
        names = sorted(
            {PROPHAGE_INT[x] for x in ints.loc[(ints["start"] <= b) & (ints["end"] >= a), "name"]}
        )
        row = {
            "region_id": rid,
            "start": a,
            "end": b,
            "named_by_integrase": names,
            "n_genes": clean(r0["n_genes"]),
            "n_eligible": clean(r0["n_eligible"]),
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
                "p_value": clean(o["p_value"]),
                "q_value": clean(o["q_value"]),
            }
        if len(hl):
            h = hl.iloc[0]
            row |= {"hl_shift_log2": clean(h["effect"]), "hl_q_value": clean(h["q_value"])}
        regs.append(row)
    lower = [r for r in regs if r.get("o_minus_e") is not None and r["o_minus_e"] < 0]
    sig = [r for r in lower if r.get("q_value") is not None and r["q_value"] <= 0.05]
    verdict = "not testable" if not regs else "pass" if len(sig) == len(regs) else "fail"
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
                    "start",
                    "end",
                    "n_genes",
                    "detected_frac",
                    "mean_abundance_log2",
                    "effect",
                    "q_value",
                )
            }
            for _, r in isl.iterrows()
        ]
    )
    return {
        "regions": regs,
        "n_regions": len(regs),
        "n_lower": len(lower),
        "n_lower_q05": len(sig),
        "verdict": verdict,
        "islands": islands,
    }


def operons(out: Path) -> dict[str, Any]:
    """Pairs of adjacent genes that are both ribosomal-protein or ATP-synthase genes: the
    within-operon abundance test of ``pairs.tsv`` (q <= 0.05 = discordant) and any
    ``operon_flags.tsv`` row whose abundance pattern is not concordant/coabundant."""
    g = tbl(out, "genes")
    name = dict(
        zip(g["locus_tag"].astype(str), g["gene"].map(lambda v: clean(v) or ""), strict=False)
    )
    ra = {t for t, n in name.items() if RIBO_ATP.match(n)}
    pairs = tbl(out, "pairs")
    fl = tbl(out, "operon_flags")
    res: dict[str, Any] = {"n_ribo_atp_genes": len(ra)}
    disc: list[dict[str, Any]] = []
    if pairs is not None:
        inside = pairs[pairs["gene_a"].astype(str).isin(ra) & pairs["gene_b"].astype(str).isin(ra)]
        q = pd.to_numeric(inside["q_value"], errors="coerce")
        d = pd.to_numeric(inside["delta_log2"], errors="coerce").abs()
        both = inside["detected_a"].map(lambda v: str(v).lower() == "true") & inside[
            "detected_b"
        ].map(lambda v: str(v).lower() == "true")
        res |= {
            "n_pairs": len(inside),
            "n_pairs_both_detected": int(both.sum()),
            "pair_class": dict(Counter(inside["pair_class"].astype(str))),
            "n_pairs_q_le_0.05": int((q <= 0.05).sum()),
            "median_abs_delta_log2": clean(d.median()),
            "max_abs_delta_log2": clean(d.max()),
            "min_q": clean(q.min()),
        }
        for _, r in inside[q <= 0.05].iterrows():
            disc.append(
                {
                    "gene_a": name.get(str(r["gene_a"])),
                    "gene_b": name.get(str(r["gene_b"])),
                    "delta_log2": clean(r["delta_log2"]),
                    "q_value": clean(r["q_value"]),
                    "source": "pairs",
                }
            )
    if fl is not None:
        f = fl[fl["gene_a"].astype(str).isin(ra) & fl["gene_b"].astype(str).isin(ra)]
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
        for x in res["operon_flags"]:
            if x["abundance_pattern"] not in ("concordant", "coabundant", "none", None):
                disc.append(x | {"source": "operon_flags"})
        res["n_flags_all"] = len(fl)
        res["flag_counts_all"] = dict(Counter(fl["flags"].astype(str)))
    res["discordant"] = disc
    res["verdict"] = "not testable" if pairs is None else ("pass" if not disc else "fail")
    return res


def seleno(out: Path) -> dict[str, Any]:
    res = {}
    for gname, (a, b) in SEC.items():
        gr = gene_row(out, LOCI[gname])
        strand = gr["strand"] if gr else ("+" if gname == "fdnG" else "-")
        h = (
            hits_in(out, "U00096.3", int(gr["start"]), int(gr["end"]), strand)
            if gr
            else pd.DataFrame()
        )
        ex = h[h["route"].astype(str) == "exact"] if len(h) else h
        s = pd.to_numeric(ex["start"], errors="coerce") if len(ex) else pd.Series(dtype=float)
        e = pd.to_numeric(ex["end"], errors="coerce") if len(ex) else pd.Series(dtype=float)
        left = ex[e < a] if len(ex) else ex
        right = ex[s > b] if len(ex) else ex
        across = ex[(s <= a) & (e >= b)] if len(ex) else ex
        res[gname] = {
            "gene": gr,
            "sec_codon": [a, b],
            "peptides_5prime_genomic_left": int(left["peptide"].nunique()) if len(ex) else 0,
            "peptides_genomic_right": int(right["peptide"].nunique()) if len(ex) else 0,
            "peptides_across_sec": int(across["peptide"].nunique()) if len(ex) else 0,
        }
        r = res[gname]
        r["both_sides"] = r["peptides_5prime_genomic_left"] > 0 and r["peptides_genomic_right"] > 0
    det = [g for g, r in res.items() if r["gene"] and r["gene"]["detection_status"] == "detected"]
    verdict = (
        "not testable"
        if not det
        else ("pass" if all(res[g]["both_sides"] for g in det) else "fail")
    )
    return {"genes": res, "detected": det, "verdict": verdict}


def e3_set(out: Path, pep: Path, ko: list[str]) -> dict[str, Any]:
    P = read_peptide_table(pep)[0]
    acc = pm._target(P) & (pm._q(P) <= pm.Q_ACCEPT)
    cas_peps: dict[str, set[str]] = {}
    for p, ids, ok in zip(P["peptide"].astype(str), P["protein_ids"], acc, strict=True):
        if not ok:
            continue
        for i in pm.items(ids):
            if i.startswith("CASSETTE_"):
                cas_peps.setdefault(i, set()).add(p)
    um = tbl(out, "unmapped", "mapping")
    hits = tbl(out, "peptide_hits", "mapping")
    lab = (
        {}
        if um is None
        else dict(zip(um["sequence_id"].astype(str), um["label"].astype(str), strict=False))
    )
    routes = pm._routes(hits) if hits is not None else {}
    cassette = {}
    for prot, peps in sorted(cas_peps.items()):
        cassette[prot] = {
            "distinct_peptides": len(peps),
            "labels": dict(Counter(lab.get(p, "absent") for p in peps)),
            "exact_on_genome": sum("exact" in routes.get(p, set()) for p in peps),
        }
    genes = {g: gene_row(out, LOCI[g]) for g in ("lon", "clpP", "mntS")}
    return {"knockouts": ko, "genes": genes, "cassette": cassette, "summary": arm_summary(out)}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--e1", type=Path, required=True)
    ap.add_argument("--pep-e1", type=Path, required=True)
    ap.add_argument("--e2", type=Path)
    ap.add_argument("--pep-e2", type=Path)
    ap.add_argument("--cmp", type=Path)
    ap.add_argument("--e3", action="append", default=[], help="LABEL=OUT=PEP=KO1,KO2")
    ap.add_argument("--e4", type=Path)
    ap.add_argument("--patched", type=Path, help="round-2 run (makedb --patch)")
    ap.add_argument("--bakta-json", type=Path, required=True)
    ap.add_argument("--genome", type=Path, required=True)
    ap.add_argument("--refseq-gff", type=Path, required=True)
    ap.add_argument("--contaminants", type=Path, required=True)
    ap.add_argument("--q-crossings", type=int, default=0)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    doc: dict[str, Any] = {"schema": "gippa-campaign-phase1/1"}
    ref = refseq(a.refseq_gff)
    e1 = a.e1
    man1 = manifest(e1)
    cfg1 = pm._cfg(man1, e1 / "run_manifest.json")
    genes1 = tbl(e1, "genes")
    P1 = read_peptide_table(a.pep_e1)[0]
    hits1, unm1 = tbl(e1, "peptide_hits", "mapping"), tbl(e1, "unmapped", "mapping")
    doc["E1"] = arm_summary(e1)
    coding, circular = read_coding_genes(a.bakta_json, a.genome)
    te1, exc = pm.te1_real(
        P1, hits1, unm1, coding, circular, non_atg_nterm(coding), cfg1, q_crossings=a.q_crossings
    )
    doc["row1_te1"] = te1 | {"n_exceptions": len(exc)}
    t3 = tr.te3(genes1)
    t3.pop("_rows", None)
    doc["row2_te3"] = t3
    doc["row3_te4"] = tr.te4(genes1)
    heads = [ln[1:] for ln in a.contaminants.read_text().splitlines() if ln.startswith(">")]
    # gippa makedb writes the contaminants as CONT_<id> without their descriptions, so the
    # pilot's family rule (on the description) is applied to the source FASTA and the IDs are
    # given makedb's prefix; both spellings are kept.
    fam = {
        f: set(ids) | {"CONT_" + i for i in ids}
        for f, ids in tr.contaminant_families(heads).items()
    }
    t5 = tr.te5_arm(P1, hits1, unm1, fam, cfg1)
    t5.pop("_violations", None)
    doc["row4_te5"] = t5
    doc["row5_prophages"] = prophages(e1, ref)
    doc["row6_operons"] = operons(e1)
    doc["row10_seleno"] = seleno(e1)
    doc["row7_crl"] = {
        "E1_gene": gene_row(e1, LOCI["crl"]),
        "E1_windows": window_counts(e1, CRL_WINDOWS, "+"),
    }
    doc["row8_prfb"] = {
        "E1_gene": gene_row(e1, LOCI["prfB"]),
        "E1_windows": window_counts(e1, PRFB_WINDOWS, "-"),
    }
    if a.e2:
        doc["E2"] = arm_summary(a.e2)
        doc["row7_crl"] |= {
            "E2_gene": gene_row(a.e2, LOCI["crl"]),
            "E2_windows": window_counts(a.e2, CRL_WINDOWS, "+"),
            "E2_novel": novel_in(a.e2, 257829, 259006, "+"),
        }
        doc["row8_prfb"] |= {
            "E2_gene": gene_row(a.e2, LOCI["prfB"]),
            "E2_windows": window_counts(a.e2, PRFB_WINDOWS, "-"),
            "E2_novel": novel_in(a.e2, 3035184, 3036282, "-"),
        }
        doc["row9_novel"] = novel_vs_refseq(a.e2, ref)
        w = [x for x in (doc["E2"]["warnings"] or []) if "class_q_floor" in x]
        doc["row15_q_floor"] = {"warnings": w}
        d1, d2 = detected(e1), detected(a.e2)
        doc["row14_route_b"] = {
            "E1_detected": len(d1),
            "E2_detected": len(d2),
            "E1_kept_by_E2": share(len(d1 & d2), len(d1)),
            "E2_only": len(d2 - d1),
        }
        doc["row10_seleno_E2"] = seleno(a.e2)
    if a.cmp:
        ab = a.cmp / "tables" / "compare_abundance.tsv"
        st = a.cmp / "tables" / "compare_stability.tsv"
        doc["compare"] = {
            "abundance": pd.read_csv(ab, sep="\t").to_dict("records") if ab.is_file() else None,
            "stability": pd.read_csv(st, sep="\t").to_dict("records") if st.is_file() else None,
        }
    if a.e3:
        doc["row11_e3"] = {}
        for spec in a.e3:
            label, out, pep, ko = spec.split("=")
            doc["row11_e3"][label] = e3_set(Path(out), Path(pep), [k for k in ko.split(",") if k])
    if a.e4:
        v = tbl(a.e4, "variants")
        g4 = tbl(a.e4, "genes")
        ds = tbl(a.e4, "detection_summary")
        t2 = tr.te2(v, g4, ds)
        t2.pop("_rows", None)
        doc["E4"] = arm_summary(a.e4)
        doc["row12_te2"] = t2
        if v is not None:
            fl = v["flags"].astype(str)
            me = v[fl.str.contains("mass_equivalent")]
            doc["row13_mass_equivalent"] = {
                "n_alt_sites": len(v),
                "n_mass_equivalent": len(me),
                "by_alt_type": dict(Counter(me["alt_type"].astype(str))),
                "by_ref_alt": dict(
                    Counter(me["ref_aa"].astype(str) + ">" + me["alt_aa"].astype(str))
                ),
                "by_confidence": dict(Counter(me["confidence"].astype(str))),
                "high_not_mass_equivalent": int(
                    (
                        (v["confidence"].astype(str) == "high")
                        & ~fl.str.contains("mass_equivalent")
                    ).sum()
                ),
            }
    if a.patched:
        doc["round2"] = arm_summary(a.patched)
    a.out.write_text(json.dumps(doc, indent=1, default=str, sort_keys=False) + "\n")
    print(f"wrote {a.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
