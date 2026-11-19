"""Phase 1b scoring (checks 1b-1 to 1b-4 and the ortholog half of 1b-6; ANSWER_KEY_1b.md).

    python -I score.py D SEED      # D = $GIPPA_ROOT/campaign/ncm3722

Reads $D/truth/{truth_sites,orthologs}.tsv, $D/out_n1 and $D/out_n2 tables; writes
$D/score/score.json (aggregates) and $D/score/candidates.tsv (coordinates, flags, counts; no
peptide sequences). Nothing here is committed except the aggregates copied into REPORT_1b.md.
"""

from __future__ import annotations

import csv
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

MAXD = {10: 1, 11: 1}  # mapping.tolerant.max_diff: 1 at 10-11 aa, 2 at 12-40 aa


def rd(p: Path) -> list[dict[str, str]]:
    with p.open() as fh:
        return list(csv.DictReader((ln for ln in fh if not ln.startswith("#")), delimiter="\t"))


def wilson(k: int, n: int, z: float = 1.96) -> dict[str, float | int | None]:
    if n == 0:
        return {"k": k, "n": n, "p": None, "lo": None, "hi": None}
    p = k / n
    c = (p + z * z / (2 * n)) / (1 + z * z / n)
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return {"k": k, "n": n, "p": round(p, 4), "lo": round(max(0, c - h), 4), "hi": round(min(1, c + h), 4)}


def fold(s: str) -> str:
    return s.replace("I", "L").replace("J", "L")


def fnum(x: str) -> float:
    try:
        return float(x)
    except ValueError:
        return math.nan


def covers(h: dict[str, str], contig: str, strand: str, c: int) -> bool:
    if h["contig"] != contig or h["strand"] != strand or h["route"] != "exact":
        return False
    s, e = int(h["start"]), int(h["end"])
    if not (s <= c and c + 2 <= e):
        return False
    return (c - s) % 3 == 0 if strand == "+" else (e - (c + 2)) % 3 == 0


def main(d: Path, seed: int) -> None:
    out = d / "score"
    out.mkdir(exist_ok=True)
    truth = rd(d / "truth" / "truth_sites.tsv")
    orth = {r["ncm_locus"]: r["mg_locus"] for r in rd(d / "truth" / "orthologs.tsv")}
    h1 = [h for h in rd(d / "out_n1" / "mapping" / "peptide_hits.tsv") if fnum(h["best_q_input"]) <= 0.01
          or h["best_q_input"] in ("", "NA")]
    h2 = rd(d / "out_n2" / "mapping" / "peptide_hits.tsv")
    var = rd(d / "out_n2" / "tables" / "variants.tsv")
    res: dict[str, object] = {}

    # ---- coverage of truth sites by N1 peptides (exact, in the ortholog's frame)
    by_contig: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for h in h1:
        by_contig[(h["contig"], h["strand"])].append(h)
    sites_by_ncm: dict[str, list[dict[str, str]]] = defaultdict(list)
    for t in truth:
        sites_by_ncm[t["ncm_locus"]].append(t)
    covered, eligible = [], []
    for t in truth:
        c = int(t["ncm_codon_start"])
        hs = [h for h in by_contig[(t["ncm_contig"], t["ncm_strand"])] if covers(h, t["ncm_contig"], t["ncm_strand"], c)]
        if not hs:
            continue
        covered.append(t)
        ok = False
        for h in hs:
            ln = int(h["peptide_len_aa"])
            if ln < 10:
                continue
            nd = sum(covers(h, u["ncm_contig"], u["ncm_strand"], int(u["ncm_codon_start"])) for u in sites_by_ncm[t["ncm_locus"]])
            ok |= nd <= MAXD.get(ln, 2)
        if ok:
            eligible.append(t)

    # ---- N2 calls matched to truth
    tkey = {(t["mg_contig"], t["mg_strand"], int(t["mg_codon_start"])): t for t in truth}

    def match(v: dict[str, str]) -> list[tuple[str, str, int]]:
        if v["alt_type"] == "frameshift":
            return []
        c0, st, alt = int(v["codon_start"]), v["strand"], v["alt_aa"]
        step = 3 if st == "+" else -3
        hits = []
        for k, a in enumerate(alt):
            key = (v["contig"], st, c0 + k * step)
            t = tkey.get(key)
            if t is not None:
                if fold(a) != fold(t["alt_aa"]):
                    return []
                hits.append(key)
        if v["alt_type"] != "substitution" and len(alt) * 3 != len(v["ref_codon"]):
            return []
        return hits

    called_any, called_high = set(), set()
    for v in var:
        if v["call"] not in ("ALT", "mixed"):
            continue
        for key in match(v):
            called_any.add(key)
            if v["confidence"] == "high":
                called_high.add(key)

    def k_of(ts: list[dict[str, str]], s: set) -> int:
        return sum((t["mg_contig"], t["mg_strand"], int(t["mg_codon_start"])) in s for t in ts)

    res["truth_sites"] = len(truth)
    res["covered_sites"] = len(covered)
    res["eligible_sites"] = len(eligible)
    res["recall_any_covered"] = wilson(k_of(covered, called_any), len(covered))
    res["recall_high_covered"] = wilson(k_of(covered, called_high), len(covered))
    res["recall_any_eligible"] = wilson(k_of(eligible, called_any), len(eligible))
    res["recall_high_eligible"] = wilson(k_of(eligible, called_high), len(eligible))

    # ---- precision of high ALT calls; candidates = high calls not in the truth set
    high_alt = [v for v in var if v["confidence"] == "high" and v["call"] == "ALT"]
    high_mixed = [v for v in var if v["confidence"] == "high" and v["call"] == "mixed"]
    tp = [v for v in high_alt if match(v)]
    res["precision_high_alt"] = wilson(len(tp), len(high_alt))
    res["precision_high_mixed"] = wilson(sum(bool(match(v)) for v in high_mixed), len(high_mixed))
    res["variants_rows"] = len(var)
    res["variants_by_call_conf"] = {f"{c}_{q}": sum(v["call"] == c and v["confidence"] == q for v in var)
                                    for c in ("ALT", "mixed") for q in ("high", "candidate")}
    cands = []
    for v in high_alt + high_mixed:
        if match(v):
            continue
        fl = set(v["flags"].split(";")) - {""}
        if "mass_equivalent" in fl:
            cls = "mass_equivalent"
        elif fl & {"start_codon", "stop_codon", "putative_start"}:
            cls = "start_or_stop_codon"
        else:
            cls = "unexplained"
        cands.append({"variant_id": v["variant_id"], "contig": v["contig"], "codon_start": v["codon_start"],
                      "strand": v["strand"], "gene": v["gene"], "aa_pos": v["aa_pos"], "ref_aa": v["ref_aa"],
                      "alt_aa": v["alt_aa"], "alt_type": v["alt_type"], "call": v["call"], "n_ref": v["n_ref"],
                      "n_alt": v["n_alt"], "flags": v["flags"], "class": cls})
    res["candidates_by_class"] = {c: sum(x["class"] == c for x in cands)
                                  for c in ("mass_equivalent", "start_or_stop_codon", "unexplained")}
    if cands:
        with (out / "candidates.tsv").open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(cands[0]), delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(cands)

    # ---- unplaced share of N1 peptides on MG1655, by N1 replicon
    origin: dict[str, set[str]] = defaultdict(set)
    for h in h1:
        origin[fold(h["peptide"])].add(h["contig"])
    placed2 = {fold(h["peptide"]) for h in h2}
    contigs = sorted({c for s in origin.values() for c in s})
    groups: dict[str, list[str]] = defaultdict(list)
    for p, cs in origin.items():
        groups[next(iter(cs)) if len(cs) == 1 else "both"].append(p)
    res["unplaced_by_origin"] = {g: wilson(sum(p not in placed2 for p in ps), len(ps)) for g, ps in sorted(groups.items())}
    res["unplaced_all"] = wilson(sum(p not in placed2 for p in origin), len(origin))
    res["n1_contigs"] = contigs
    res["n2_route"] = {r: sum(h["route"] == r for h in h2) for r in ("exact", "tolerant")}
    um = rd(d / "out_n2" / "mapping" / "unmapped.tsv")
    res["n2_unmapped_reasons"] = {}
    for u in um:
        k = f"{u['label']}:{u['reason']}"
        res["n2_unmapped_reasons"][k] = res["n2_unmapped_reasons"].get(k, 0) + 1

    # ---- run-level mapping_fdr placement totals
    res["n2_mapping_fdr_all"] = [
        {k: r[k] for k in ("decoy_type", "n_diff", "targets", "decoys", "pi_wrong", "fdr_stratum", "n_accepted")}
        for r in rd(d / "out_n2" / "mapping" / "mapping_fdr.tsv")
        if r["section"] == "placement" and r["length_class"] == "all"]

    # ---- ortholog-joined comparison (1b-6)
    g1 = {r["locus_tag"]: r for r in rd(d / "out_n1" / "tables" / "genes.tsv")}
    g2 = {r["locus_tag"]: r for r in rd(d / "out_n2" / "tables" / "genes.tsv")}
    for tag, g in (("n1", g1), ("n2", g2)):
        res[f"{tag}_genes"] = len(g)
        res[f"{tag}_detected"] = sum(r["detected"] == "true" for r in g.values())
        res[f"{tag}_tier"] = next(iter(g.values()))["detection_tier"]
    pairs = [(n, m) for n, m in orth.items() if n in g1 and m in g2]
    res["ortholog_pairs_in_tables"] = len(pairs)
    det_both = [(n, m) for n, m in pairs if g1[n]["detected"] == "true" and g2[m]["detected"] == "true"]
    det1 = {n for n, m in pairs if g1[n]["detected"] == "true"}
    det2 = {n for n, m in pairs if g2[m]["detected"] == "true"}
    res["detected_jaccard_orthologs"] = round(len(det1 & det2) / max(1, len(det1 | det2)), 4)
    res["n1_detected_orthologs_kept_in_n2"] = wilson(len(det1 & det2), len(det1))
    ab = [(fnum(g1[n]["abundance_log2"]), fnum(g2[m]["abundance_log2"])) for n, m in det_both
          if g1[n]["abundance_imputed"] != "true" and g2[m]["abundance_imputed"] != "true"]
    ab = [(a, b) for a, b in ab if not (math.isnan(a) or math.isnan(b))]
    x = np.array(ab)
    rho = float(spearmanr(x[:, 0], x[:, 1]).statistic)
    rng = np.random.default_rng(seed)
    boots = []
    for _ in range(1000):
        i = rng.integers(0, len(x), len(x))
        boots.append(spearmanr(x[i, 0], x[i, 1]).statistic)
    res["spearman_orthologs"] = {"n": len(x), "rho": round(rho, 4), "lo": round(float(np.quantile(boots, 0.025)), 4),
                                 "hi": round(float(np.quantile(boots, 0.975)), 4), "draws": 1000, "seed": seed}
    m2n = {m: n for n, m in orth.items()}
    for tab in ("islands",):
        s1, s2 = set(), set()
        for r in rd(d / "out_n1" / "tables" / f"{tab}.tsv"):
            s1 |= {(r["direction"], gname) for gname in r["genes"].split(";") if gname in orth}
        for r in rd(d / "out_n2" / "tables" / f"{tab}.tsv"):
            s2 |= {(r["direction"], m2n[gname]) for gname in r["genes"].split(";") if gname in m2n}
        res["islands_n1"] = len(rd(d / "out_n1" / "tables" / "islands.tsv"))
        res["islands_n2"] = len(rd(d / "out_n2" / "tables" / "islands.tsv"))
        res["island_gene_jaccard_orthologs"] = {"n1": len(s1), "n2": len(s2), "shared": len(s1 & s2),
                                                "jaccard": round(len(s1 & s2) / max(1, len(s1 | s2)), 4)}
    f1 = {(frozenset((r["gene_a"], r["gene_b"])), r["flags"]) for r in rd(d / "out_n1" / "tables" / "operon_flags.tsv")
          if r["gene_a"] in orth and r["gene_b"] in orth}
    f2 = {(frozenset((m2n[r["gene_a"]], m2n[r["gene_b"]])), r["flags"]) for r in rd(d / "out_n2" / "tables" / "operon_flags.tsv")
          if r["gene_a"] in m2n and r["gene_b"] in m2n}
    res["operon_flags_n1"] = len(rd(d / "out_n1" / "tables" / "operon_flags.tsv"))
    res["operon_flags_n2"] = len(rd(d / "out_n2" / "tables" / "operon_flags.tsv"))
    res["operon_flag_jaccard_orthologs"] = {"n1": len(f1), "n2": len(f2), "shared": len(f1 & f2),
                                            "jaccard": round(len(f1 & f2) / max(1, len(f1 | f2)), 4)}
    (out / "score.json").write_text(json.dumps(res, indent=2) + "\n")
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main(Path(sys.argv[1]), int(sys.argv[2]) if len(sys.argv) > 2 else 20261014)
