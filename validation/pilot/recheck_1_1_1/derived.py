# ruff: noqa: B023, E402, E501, E731, SIM115  (copied as run for the recheck)
"""Derived pilot metrics (m3 coverage, h1, m4, m5 enriched categories, h3) from regenerated runs,
with the freeze validation/pilot/metrics.py functions. Aggregates only. Usage: derived.py VARIANT..."""

import json
import sys
from pathlib import Path

WT = Path("/home/classes/bios/270/khoa/giap/local/worktrees/pilot-fa49b8f145ea")
sys.path.insert(0, str(WT))
from validation.pilot import metrics as M
from validation.proxies.te_proxies import read_coding_genes

from gippa.io.protein_fasta import read_protein_fasta

S = Path(__file__).parent
P3 = Path("/home/classes/bios/270/khoa/giap/local/pilot/pass3")
genes, _ = read_coding_genes(
    P3 / "inputs/nissle1917_wt.bakta.json", P3 / "genome/nissle1917_wt.fna.gz"
)
k12 = {
    r.id: r.sequence
    for r in read_protein_fasta(
        Path("/home/classes/bios/270/khoa/giap/local/data/proteomes/UP000000625.fasta")
    )
}
classes = M.strain_specific_set({g.locus_tag: g.aa for g in genes}, k12)
for v in sys.argv[1:]:
    R = S / "runs" / v
    tab = lambda a, n: M._table(R / a / "out" / "tables" / f"{n}.tsv", n)
    g1, g2 = tab("a1", "genes"), tab("a2", "genes")
    cmp = {s: M._table(R / "cmp" / "tables" / f"{s}.tsv", s) for s in M.COMPARE}
    doc = {"variant": v}
    doc["m3"] = {
        a: {
            k: [
                x["detected"],
                x["genes"],
                x["share"]["value"],
                x["share"]["ci_low"],
                x["share"]["ci_high"],
            ]
            for k, x in M.coverage(g, classes)["subsets"].items()
        }
        for a, g in (("A1", g1), ("A2", g2))
    }
    h = M.hypothesis1(g1, g2, cmp["compare_detection"], M.gene_subsets(g2, classes))
    doc["h1"] = {k: h[k] for k in h}
    doc["m4"] = M.abundance_agreement(cmp["compare_abundance"], M.ARMS)
    doc["m5"] = M.stability(cmp["compare_stability"], M.ARMS)
    doc["h3"] = M.hypothesis3(doc["m4"], doc["m5"])
    json.dump(doc, open(S / f"derived_{v}.json", "w"), indent=1, default=str)
    print(v, "written")
