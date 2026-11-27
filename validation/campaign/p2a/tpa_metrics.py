"""PXD065022 metrics.json from the _fix2 runs (pin rel111c-pin, main 38e1cdba); aggregates only.

Reads the per-strain outputs of ``tpa_score.sbatch`` ($C/score_fix2/TAG/{p2,prov_*,placement_fdr}.json,
t6/score.json, te7_diff.txt) and the route-A tables it needs for the dataset checks T2-T5, and the
Nextflow task logs for the gap-fix warnings (#611/#614/#615). Verdict rules: ANSWER_KEY.md.
No peptide sequence and no DIA-NN row is written.

    PYTHONPATH=<pin>/src python tpa_metrics.py OUT.json
"""

from __future__ import annotations

import csv
import json
import math
import re
import statistics
import sys
from pathlib import Path
from typing import Any

C = Path("/home/classes/bios/270/khoa/giap/local/campaign/PXD065022")
KEY = Path(__file__).resolve().parent.parent / "PXD065022"
SUF = "_fix2"
TAGS = ["tpa_dal1", "tpa_grady", "tpa_haitib", "tpa_madras", "tpa_mexicoa", "tpa_phil1", "tpa_ss14"]
NICHOLS = ["tpa_dal1", "tpa_haitib", "tpa_madras"]
SS14 = ["tpa_ss14", "tpa_mexicoa", "tpa_grady", "tpa_phil1"]
# T5: DAL-1 locus of each paper protein and the paper's direction (+1 = higher in Nichols-like)
T5 = {"TP0477_pgl": ("CLKMHN_00468", 1), "TP0487": ("CLKMHN_00477", 1), "TP0597": ("CLKMHN_00595", 1),
      "TP0594": ("CLKMHN_00592", -1), "TP0732_folD": ("CLKMHN_00732", -1)}
WARN = {
    "tu_upstream_skipped": re.compile(r"operon_pairs: (\d+) TU upstream region\(s\) contain non-ACGT"),
    "promotech_windows_skipped": re.compile(r"promoter_layer: (\d+) window position\(s\) contain non-ACGT"),
    "terminatornet": re.compile(r"TerminatorNet: the genome has (\d+) non-ACGT base\(s\), removed before the "
                                r"scan; .* and (\d+) call\(s\) whose hairpin"),
}


def wilson(k: int, n: int, z: float = 1.96) -> dict[str, Any]:
    if n == 0:
        return {"k": k, "n": n, "value": None, "ci_low": None, "ci_high": None}
    p = k / n
    c = (p + z * z / (2 * n)) / (1 + z * z / n)
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return {"k": k, "n": n, "value": p, "ci_low": max(0.0, c - h), "ci_high": min(1.0, c + h)}


def jl(p: Path) -> Any:
    return json.loads(p.read_text()) if p.is_file() else None


def rd(p: Path) -> list[dict[str, str]]:
    with p.open() as fh:
        return list(csv.DictReader((ln for ln in fh if not ln.startswith("#")), delimiter="\t"))


def fnum(x: Any) -> float:
    try:
        return float(x)
    except (TypeError, ValueError):
        return math.nan


def gapfix(out: Path) -> dict[str, Any] | None:
    """The three gap-fix warnings, read from the run's Nextflow task logs (.command.err)."""
    if not (out / "work").is_dir():
        return None
    res: dict[str, Any] = {"tu_upstream_skipped": None, "promotech_windows_skipped": None,
                           "terminatornet_removed_bases": None, "terminatornet_dropped_calls": None,
                           "terminatornet_shards_reporting": 0}
    for err in sorted(out.glob("work/*/*/.command.err")):
        txt = err.read_text(errors="replace")
        for k, rx in WARN.items():
            for m in rx.finditer(txt):
                if k == "terminatornet":
                    res["terminatornet_shards_reporting"] += 1
                    res["terminatornet_removed_bases"] = int(m.group(1))
                    res["terminatornet_dropped_calls"] = (res["terminatornet_dropped_calls"] or 0) + int(m.group(2))
                else:
                    res[k] = (res[k] or 0) + int(m.group(1))
    return res


def genes(tag: str, route: str = "a") -> dict[str, dict[str, str]]:
    return {r["locus_tag"]: r for r in rd(C / tag / f"out_{route}{SUF}" / "tables" / "genes.tsv")}


def det(r: dict[str, str] | None) -> bool:
    return bool(r) and r["detection_status"] == "detected" and not (fnum(r["gene_q_value"]) > 0.01)


def t3_tp0134a() -> dict[str, Any]:
    """Accepted route-B peptides (exact, q <= 0.01) lying inside 156,669-157,301 (+) in that frame."""
    lo, hi = 156669, 157301
    hits = rd(C / "tpa_dal1" / f"out_b{SUF}" / "mapping" / "peptide_hits.tsv")
    inside, inframe = set(), set()
    for h in hits:
        if h["contig"] != "NC_016844.1" or h["strand"] != "+" or h["route"] != "exact":
            continue
        if fnum(h["best_q_input"]) > 0.01:
            continue
        s, e = int(h["start"]), int(h["end"])
        if s >= lo and e <= hi:
            inside.add(h["peptide"])
            if (s - lo) % 3 == 0:
                inframe.add(h["peptide"])
    return {"window": [lo, hi], "strand": "+", "accepted_peptides_inside": len(inside),
            "accepted_peptides_inside_inframe": len(inframe), "recovered": len(inframe) >= 1}


def main(out: Path) -> int:
    doc: dict[str, Any] = {"schema": "gippa-campaign-phase2a-tpa/1", "dataset": "PXD065022",
                           "runs_suffix": SUF, "strains": {}}
    for t in TAGS:
        S = C / f"score{SUF}" / t
        s: dict[str, Any] = {"metrics": jl(S / "p2.json")}
        for r in ("a", "b", "x", "te7"):
            p = jl(S / f"prov_{r}.json")
            if p is not None:
                g = p.get("gippa_run") or {}
                s[f"prov_{r}"] = {"diann": p.get("diann"), "convert": p.get("convert"),
                                  "gippa_run_runtime": p.get("gippa_run_runtime"),
                                  "gippa_run": {k: v for k, v in g.items() if k != "inputs"} | {
                                      "inputs": [{k: i.get(k) for k in ("role", "sha256")}
                                                 for i in g.get("inputs") or []]}}
        s["gapfix"] = {r: gapfix(C / t / d) for r, d in (("a", f"out_a{SUF}"), ("b", f"out_b{SUF}"),
                                                          ("x_tpa_dal1", f"out_x_tpa_dal1{SUF}"),
                                                          ("te7", f"out_a{SUF}_te7"))
                       if (C / t / d).is_dir()}
        if t != "tpa_dal1":
            s["T6_score"] = jl(S / "t6" / "score.json")
            s["T6_placement_fdr"] = jl(S / "placement_fdr.json")
            s["T6_truth_counts"] = jl(C / "truth" / t / "counts.json")
        doc["strains"][t] = s

    # T2 detection per strain
    t2 = {}
    for t in TAGS:
        a = (doc["strains"][t]["metrics"] or {}).get("A") or {}
        t2[t] = wilson(a.get("n_detected", 0), a.get("n_coding_genes", 0))
    ok_a = all(v["value"] is not None and 0.50 <= v["value"] <= 0.95 for v in t2.values())
    fewest = min(t2, key=lambda k: t2[k]["k"])
    doc["T2"] = {"per_strain": t2, "a_all_in_range": ok_a, "fewest": fewest, "b_grady_fewest": fewest == "tpa_grady",
                 "verdict": "pass" if ok_a and fewest == "tpa_grady" else "fail"}

    # T3 the 12 newly annotated genes
    loci = [r for r in csv.DictReader((KEY / "answer_loci.tsv").open(), delimiter="\t") if r["key"] == "new12"]
    rows, k = [], 0
    for r in loci:
        if r["label"] == "TP0134a":
            x = t3_tp0134a()
            row = {"label": r["label"], "genome": r["genome"], "rule": "route-B peptide in window, in frame"} | x
        else:
            ga, gb = genes(r["genome"], "a").get(r["bakta_locus"]), None
            if (C / r["genome"] / f"out_b{SUF}" / "tables" / "genes.tsv").is_file():
                gb = genes(r["genome"], "b").get(r["bakta_locus"])
            row = {"label": r["label"], "genome": r["genome"], "bakta_locus": r["bakta_locus"],
                   "A": None if ga is None else {k2: ga[k2] for k2 in ("detection_status", "gene_q_value",
                                                                       "n_peptides_unique", "n_peptides_shared")},
                   "B": None if gb is None else {k2: gb[k2] for k2 in ("detection_status", "gene_q_value",
                                                                       "n_peptides_unique", "n_peptides_shared")},
                   "recovered": det(ga) or det(gb)}
        k += row["recovered"]
        rows.append(row)
    doc["T3"] = {"rows": rows, "recovered": wilson(k, len(rows)), "threshold": 9,
                 "verdict": "pass" if k >= 9 else "fail"}

    # T4 translated pseudogene distal parts (route A)
    t4 = {}
    for lab, tag, loc in (("TP0865a", "tpa_haitib", "FOKJEF_00855"), ("TP0924a", "tpa_madras", "BIIDIP_00910")):
        g = genes(tag, "a").get(loc)
        t4[lab] = {"strain": tag, "bakta_locus": loc, "detected_q01": det(g),
                   "row": None if g is None else {k2: g[k2] for k2 in (
                       "detection_status", "gene_q_value", "n_peptides_unique", "n_peptides_shared",
                       "coverage_frac", "abundance_log2")}}
    n4 = sum(v["detected_q01"] for v in t4.values())
    doc["T4"] = {"rows": t4, "n_detected": n4, "verdict": "pass" if n4 == 2 else "fail"}

    # T5 Nichols-like vs SS14-like (route A, median-centred on genes detected in all seven)
    orth = {t: {r["mg_locus"]: r["ncm_locus"] for r in rd(C / "truth" / t / "orthologs.tsv")} for t in TAGS[1:]}
    orth["tpa_dal1"] = {g: g for g in genes("tpa_dal1")}
    G = {t: genes(t) for t in TAGS}

    def strain_row(t: str, dal: str) -> dict[str, str] | None:
        loc = orth[t].get(dal)
        return None if loc is None else G[t].get(loc)

    core = [d for d in orth["tpa_dal1"]
            if all((r := strain_row(t, d)) is not None and r["detection_status"] == "detected"
                   and not math.isnan(fnum(r["abundance_log2"])) for t in TAGS)]
    centre = {t: statistics.median(fnum(strain_row(t, d)["abundance_log2"]) for d in core) for t in TAGS}
    t5rows, agree, testable = {}, 0, 0
    for lab, (dal, sign) in T5.items():
        vals: dict[str, float | None] = {}
        for t in TAGS:
            r = strain_row(t, dal)
            vals[t] = (fnum(r["abundance_log2"]) - centre[t]) if r is not None and r["detection_status"] == "detected" \
                and not math.isnan(fnum(r["abundance_log2"])) else None
        n_ = [vals[t] for t in NICHOLS if vals[t] is not None]
        s_ = [vals[t] for t in SS14 if vals[t] is not None]
        ok = len(n_) >= 2 and len(s_) >= 2
        diff = statistics.mean(n_) - statistics.mean(s_) if ok else None
        ag = ok and diff is not None and diff * sign > 0
        testable += ok
        agree += ag
        t5rows[lab] = {"dal1_locus": dal, "paper_sign": sign, "centred_log2": vals, "n_nichols": len(n_),
                       "n_ss14": len(s_), "testable": ok, "mean_diff_nichols_minus_ss14": diff, "agrees": ag}
    doc["T5"] = {"n_core_genes": len(core), "centre_log2": centre, "rows": t5rows, "testable": testable,
                 "agree": agree, "verdict": "not testable" if testable < 3 else ("pass" if agree >= 4 else "fail")}

    # T6 per strain: pooled TE8 placement FDR <= 1 %
    t6 = {}
    for t in TAGS[1:]:
        pf = doc["strains"][t].get("T6_placement_fdr")
        t6[t] = None if pf is None else {"pooled_decoy_fdr": pf["pooled_decoy_fdr"], "passed": pf["passed"]}
    doc["T6"] = {"per_strain": t6, "verdict": "pass" if all(v and v["passed"] for v in t6.values()) else
                 ("fail" if any(v and not v["passed"] for v in t6.values()) else "open")}

    # T7 determinism
    dfile = C / f"score{SUF}" / "tpa_dal1" / "te7_diff.txt"
    if dfile.is_file():
        n = len(dfile.read_text().splitlines())
        sa = (C / f"score{SUF}" / "tpa_dal1" / "te7_tables_a.sha256").read_text().splitlines()
        doc["T7"] = {"diff_lines": n, "n_files": len(sa), "verdict": "pass" if n == 0 else "fail"}
    out.write_text(json.dumps(doc, separators=(",", ":"), default=str) + "\n")
    print(f"wrote {out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main(Path(sys.argv[1])))
