"""Figure 6: genome context at one locus: E. coli K-12 prfB, whose N-terminal part is made by a
programmed +1 frameshift that Bakta does not model.

Reads the Phase 1 route-B run on GIPPA 1.1.1 code ($GIPPA_ROOT/campaign/ecoli_k12/runs/out_e2_fix:
mapping/peptide_hits.tsv coordinates and locus classes only, tables/genes.tsv and
tables/novel_orfs.tsv) and the RefSeq GFF of NC_000913.3 in the same directory tree. No peptide
sequence is read into the figure or its data table: the peptide track is a coverage depth.
(a) RefSeq and Bakta genes, GIPPA's two novel-ORF candidates and peptide coverage by locus
class; (b) the gene-table rows of the three genes in the window.
"""

from __future__ import annotations

import csv
import gzip
import os
from pathlib import Path

import numpy as np

from _common import OKABE, new_figure, out_dir, panel_label, parse_args, save

LO, HI = 3034950, 3036600
CLASSES = {"cds_inframe": ("in-frame CDS", OKABE["blue"]),
           "nterm_extension": ("N-terminal extension", OKABE["vermillion"]),
           "intergenic": ("intergenic", OKABE["bluish_green"])}  # fmt: skip
GENES = ["KOIAIB_03082", "KOIAIB_03083", "KOIAIB_03084"]
CANDS = ["NOVEL_U00096.3_3035184_m", "NOVEL_U00096.3_3036205_m"]


def _root() -> Path:
    root = os.environ.get("GIPPA_ROOT")
    if not root:
        raise SystemExit("GIPPA_ROOT is not set (source scripts/farmshare_env.sh)")
    return Path(root) / "campaign" / "ecoli_k12"


def _tsv(p: Path) -> list[dict]:
    with open(p) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _gene(ax, lo: int, hi: int, y: float, strand: str, col: str, label: str, h: float = 0.32,
          hatch: str | None = None) -> None:  # fmt: skip
    """A gene arrow on [lo, hi] (1-based closed), clipped to the window."""
    a, b = max(lo, LO), min(hi, HI)
    head = 45
    if strand == "-":
        xs = [a, a + head, b, b, a + head] if lo >= LO else [a, a, b, b, a]
        ys = [y, y + h / 2, y + h / 2, y - h / 2, y - h / 2] if lo >= LO else [
            y - h / 2, y + h / 2, y + h / 2, y - h / 2, y - h / 2]  # fmt: skip
    else:
        xs = [a, b - head, b, b - head, a]
        ys = [y + h / 2, y + h / 2, y, y - h / 2, y - h / 2]
    ax.fill(xs, ys, color=col if hatch is None else "white", ec=col, lw=0.8, hatch=hatch)
    ax.text((a + b) / 2, y + h / 2 + 0.06, label, ha="center", va="bottom", fontsize=6)


def main() -> None:
    args = parse_args(__doc__)
    d = out_dir(args)
    run = _root() / "runs" / "out_e2_fix"
    rows: list[dict] = []

    ref = []
    with gzip.open(_root() / "ref" / "GCF_000005845.2_ASM584v2_genomic.gff.gz", "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) == 9 and f[2] == "CDS" and int(f[4]) >= LO and int(f[3]) <= HI:
                gene = next((x[5:] for x in f[8].split(";") if x.startswith("gene=")), "?")
                ref.append((int(f[3]), int(f[4]), f[6], gene))
    genes = {g["locus_tag"]: g for g in _tsv(run / "tables" / "genes.tsv")
             if g["locus_tag"] in GENES}  # fmt: skip
    cands = {
        c["orf_id"]: c for c in _tsv(run / "tables" / "novel_orfs.tsv") if c["orf_id"] in CANDS
    }
    depth = {k: np.zeros(HI - LO + 1, dtype=int) for k in CLASSES}
    ends: dict[str, int] = {}
    with open(run / "mapping" / "peptide_hits.tsv") as fh:
        for h in csv.DictReader(fh, delimiter="\t"):
            if h["contig"] != "U00096.3":
                continue
            s, e = int(h["start"]), int(h["end"])
            if e < LO or s > HI or h["locus_class"] not in depth:
                continue
            depth[h["locus_class"]][max(s, LO) - LO : min(e, HI) - LO + 1] += 1
            if h["novel_candidate_id"] in CANDS:
                cid = h["novel_candidate_id"]
                ends[cid] = max(ends.get(cid, 0), e)

    fig = new_figure("double", 3.3)
    ax = fig.add_axes((0.09, 0.25, 0.56, 0.6))
    x = np.arange(LO, HI + 1)
    # peptide coverage, stacked by class (scaled into y 0..1.1)
    base = np.zeros_like(x, dtype=float)
    dmax = max(int(sum(v for v in depth.values()).max()), 1)
    scale = 1.1 / dmax
    for k, (lab, col) in CLASSES.items():
        top = base + depth[k] * scale
        ax.fill_between(x, base, top, step="mid", color=col, lw=0, label=lab)
        base = top
        rows.append({"panel": "a", "track": f"coverage_{k}", "max_depth": int(depth[k].max()),
                     "nt_covered": int((depth[k] > 0).sum())})  # fmt: skip
    ax.text(LO + 10, 1.15, f"peptide coverage (distinct placements; max {dmax})", fontsize=6,
            va="bottom", color="#333333")  # fmt: skip
    # tracks
    y_ref, y_bak, y_cand = 3.45, 2.6, 1.8
    for lo, hi, strand, gene in ref:
        lab = gene if gene != "prfB" else ""
        _gene(ax, lo, hi, y_ref, strand, "#777777", lab)
        rows.append({"panel": "a", "track": "refseq", "gene": gene, "start": lo, "end": hi,
                     "strand": strand})  # fmt: skip
    prfb = sorted((lo, hi) for lo, hi, _s, g in ref if g == "prfB")
    ax.text((prfb[0][0] + prfb[-1][1]) / 2, y_ref + 0.22, "prfB (2 segments)",
            ha="center", va="bottom", fontsize=6)  # fmt: skip
    slip = prfb[0][1]  # 3' end of the first segment in genome order (the frameshift site)
    ax.annotate("+1 frameshift\n(codon 26)", xy=(slip + 1, y_ref - 0.17), xytext=(slip - 260, 3.95),
                fontsize=6, ha="center", arrowprops={"arrowstyle": "-", "lw": 0.6})  # fmt: skip
    for lt, g in genes.items():
        _gene(ax, int(g["start"]), int(g["end"]), y_bak, g["strand"], OKABE["blue"],
              g["gene"] if lt != "KOIAIB_03083" else f"prfB ({lt})")  # fmt: skip
        rows.append({"panel": "a", "track": "bakta", "locus_tag": lt, "gene": g["gene"],
                     "start": g["start"], "end": g["end"], "strand": g["strand"]})  # fmt: skip
    c1 = cands[CANDS[0]]
    lo1 = int(genes["KOIAIB_03083"]["end"]) + 1
    _gene(ax, lo1, ends[CANDS[0]], y_cand, "-", OKABE["vermillion"], "", hatch="////")
    ax.text(lo1 - 40, y_cand,
            f"N-terminal extension, {c1['n_peptides']} peptides,\nclass FDR "
            f"{float(c1['class_fdr']):.3f}, no start codon", fontsize=6, ha="right", va="center",
            color=OKABE["vermillion"])  # fmt: skip
    c2 = cands[CANDS[1]]
    lo2, hi2 = int(c2["end"]), int(c2["start_longest"])
    _gene(ax, lo2, hi2, y_cand, "-", OKABE["bluish_green"], "", hatch="\\\\\\\\")
    ax.text(hi2 + 25, y_cand, f"intergenic,\n{c2['n_peptides']} peptide,\nFDR "
            f"{float(c2['class_fdr']):.3f}", fontsize=6, ha="left", va="center",
            color="#007A5E")  # fmt: skip
    for cid, c in cands.items():
        rows.append({"panel": "a", "track": "candidate", "orf_id": cid, "class": c["class"],
                     "n_peptides": c["n_peptides"], "class_fdr": c["class_fdr"],
                     "called": c["called"], "flags": c["flags"]})  # fmt: skip
    ax.text(LO - 30, y_ref, "RefSeq", ha="right", va="center", fontsize=6)
    ax.text(LO - 30, y_bak, "Bakta", ha="right", va="center", fontsize=6)
    ax.text(LO - 30, y_cand, "GIPPA\ncandidates", ha="right", va="center", fontsize=6)
    ax.set_xlim(LO, HI)
    ax.set_ylim(0, 4.4)
    ax.set_yticks([])
    for s in ("left", "right", "top"):
        ax.spines[s].set_visible(False)
    ticks = np.arange(3035000, HI, 500)
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{t:,}" for t in ticks])
    ax.set_xlabel("U00096.3 position (nt; genes on the minus strand point left)")
    ax.legend(loc="upper left", bbox_to_anchor=(0.0, -0.17), frameon=False, ncol=3, fontsize=6,
              handlelength=1.0, columnspacing=0.8)  # fmt: skip
    ax.set_title("prfB: peptides upstream of the Bakta start and in the leader", loc="left")
    panel_label(ax, "a")

    # (b) gene-table rows
    axb = fig.add_axes((0.69, 0.25, 0.3, 0.6))
    axb.axis("off")
    cols = ["gene", "detected", "gene q", "unique\npeptides", "log2\niBAQ", "TU"]
    cells = []
    for lt in GENES:
        g = genes[lt]
        cells.append([g["gene"], g["detected"],
                      f"{float(g['gene_q_value']):.1e}", g["n_peptides_unique"],
                      f"{float(g['abundance_log2']):.1f}", g["tu_id"]])  # fmt: skip
        rows.append({"panel": "b", **{k: g[k] for k in ("locus_tag", "gene", "detected",
                     "detection_tier", "gene_q_value", "n_peptides_unique", "abundance_log2",
                     "abundance_type", "tu_id", "uniop_p_prev", "uniop_p_next")}})  # fmt: skip
    tab = axb.table(
        cellText=cells,
        colLabels=cols,
        loc="upper center",
        cellLoc="center",
        colWidths=[0.13, 0.17, 0.19, 0.17, 0.13, 0.21],
    )
    tab.auto_set_font_size(False)
    tab.set_fontsize(6)
    tab.scale(1.0, 1.9)
    for (r, _c), cell in tab.get_celld().items():
        cell.set_linewidth(0.4)
        if r == 0:
            cell.set_facecolor("#EEEEEE")
    axb.text(0.0, 0.25, "genes.tsv rows, route B (gene q from the\ndecoy_picked tier; "
             "locus tags KOIAIB_03082-4).\nprfB's 17 unique "
             "peptides include the 4\nextension peptides (flag nterm_extension_\nevidence). "
             "Neither candidate is called\nat α = 0.01.", fontsize=6, va="top",
             transform=axb.transAxes)  # fmt: skip
    axb.set_title("Gene table", loc="left")
    panel_label(axb, "b")
    save(fig, d, "fig6_prfb_locus", rows)


if __name__ == "__main__":
    main()
