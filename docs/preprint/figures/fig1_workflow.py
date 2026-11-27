"""Figure 1: the GIPPA workflow and where each error estimate is made (draft schematic).

Drawn programmatically as a draft for Khoa; no data are read. Stage contents follow
docs/preprint/methods.md Sections 1-19 and the significance levels of Table M1
(src/gippa/config/defaults.yaml › stats.alpha).
"""

from __future__ import annotations

from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

from _common import OKABE, new_figure, out_dir, parse_args, save

#: (column, title, lines, error-control tag or None)
STAGES = [
    (
        0,
        "Inputs",
        [
            "genome FASTA",
            "Bakta JSON, GFF3",
            "  or GenBank",
            "peptide table",
            "  (gippa convert:",
            "  DIA-NN, Cascadia)",
            "or proteins +",
            "  abundance",
        ],
        "search-engine q\n≤ 0.01 (inherited)",
    ),
    (
        1,
        "Annotation",
        [
            "preflight checks",
            "Bakta genes, oriC",
            "geNOMAD prophages",
            "antiSMASH BGCs",
            "UniOP operons, TUs",
            "promoters,",
            "  terminators, RBS",
            "codon usage",
        ],
        None,
    ),
    (
        2,
        "Mapping",
        [
            "six-frame exact",
            "tolerant (≤ 2",
            "  differences)",
            "proteins: miniprot,",
            "  tblastn",
            "reversed or",
            "  shuffled genome",
            "  decoys",
        ],
        "placement FDR\nq ≤ 0.01",
    ),
    (
        3,
        "Gene table",
        [
            "one row per",
            "  Bakta locus tag",
            "detection tier",
            "  (picked, MAYU,",
            "  2-peptide, list)",
            "iBAQ abundance",
            "detectability",
        ],
        "gene-level FDR\nq ≤ 0.01",
    ),
    (
        4,
        "Context, calls",
        [
            "region tests",
            "autocorrelation",
            "islands",
            "operons",
            "variants",
            "novel ORFs",
            "ORA, GSEA",
            "determinants",
        ],
        "BH, FWER, class\nFDR; or flagged",
    ),
    (
        5,
        "Outputs",
        [
            "22 result tables",
            "  (TSV + schema)",
            "15 plots",
            "HTML report",
            "run manifest",
            "gippa compare",
        ],
        None,
    ),
]
TAG = OKABE["vermillion"]


def main() -> None:
    args = parse_args(__doc__)
    d = out_dir(args)
    fig = new_figure("double", 3.0)
    ax = fig.add_axes((0.0, 0.0, 1.0, 1.0))
    ax.set_xlim(0, 6)
    ax.set_ylim(0, 3.0)
    ax.axis("off")
    w, x0 = 0.86, 0.07
    rows = []
    for col, title, lines, tag in STAGES:
        x = col + x0
        box = FancyBboxPatch((x, 0.55), w, 1.95, boxstyle="round,pad=0.02,rounding_size=0.05",
                             fc="#F4F4F4", ec="#333333", lw=0.6)  # fmt: skip
        ax.add_patch(box)
        ax.text(x + w / 2, 2.42, title, ha="center", va="top", fontsize=7, fontweight="bold")
        for i, line in enumerate(lines):
            ax.text(x + 0.05, 2.22 - 0.148 * i, line, ha="left", va="top", fontsize=6)
        if tag:
            ax.add_patch(FancyBboxPatch((x + 0.04, 0.62), w - 0.08, 0.36,
                                        boxstyle="round,pad=0.01,rounding_size=0.03",
                                        fc="white", ec=TAG, lw=0.8))  # fmt: skip
            ax.text(x + w / 2, 0.80, tag, ha="center", va="center", fontsize=5.5, color=TAG)
        if col < 5:
            ax.add_patch(FancyArrowPatch((x + w + 0.01, 1.55), (x + 1 - 0.01, 1.55),
                                         arrowstyle="-|>", mutation_scale=7, lw=0.8,
                                         color="#333333"))  # fmt: skip
        rows.append({"column": col, "stage": title, "contents": "; ".join(lines),
                     "error_control": (tag or "").replace("\n", " ")})  # fmt: skip
    # side route: gippa makedb feeds the user's own spectrum search
    ax.add_patch(FancyBboxPatch((0.07, 0.08), 1.86, 0.32,
                                boxstyle="round,pad=0.02,rounding_size=0.04",
                                fc="white", ec=OKABE["blue"], lw=0.7, ls="--"))  # fmt: skip
    ax.text(1.0, 0.24, "gippa makedb: strain FASTA\n(± six-frame ORFs) for your search",
            ha="center", va="center", fontsize=5.5, color=OKABE["blue"])  # fmt: skip
    ax.add_patch(FancyArrowPatch((0.5, 0.42), (0.5, 0.55), arrowstyle="-|>", mutation_scale=6,
                                 lw=0.7, color=OKABE["blue"]))  # fmt: skip
    ax.text(5.95, 0.24, "red boxes: where GIPPA reports an error estimate",
            ha="right", va="center", fontsize=5.5, color=TAG)  # fmt: skip
    ax.text(
        3.0,
        2.88,
        "GIPPA turns one genome and any search output into one gene table, "
        "with an error estimate at every step",
        ha="center",
        va="top",
        fontsize=8,
    )
    save(fig, d, "fig1_workflow", rows)


if __name__ == "__main__":
    main()
