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
        ("search q ≤ 0.01\n(inherited)", "none"),
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
        ("presence calls only;\nmissing layer flagged", "none"),
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
        ("placement FDR\nq ≤ 0.01", "est"),
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
        ("gene FDR q ≤ 0.01\n(2-peptide: none)", "est"),
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
        ("BH FDR, FWER, class\nFDR; ORA flagged", "est"),
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
        ("every value: q, CI\nor a flag", "info"),
    ),
]
TAG = OKABE["vermillion"]
NONE = "#666666"


def _box(ax, x, y, w, h, ec, ls="-", fc="white", lw=0.7) -> None:
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.04",
                                fc=fc, ec=ec, lw=lw, ls=ls))  # fmt: skip


def _arrow(ax, a, b, col="#333333", **kw) -> None:
    ax.add_patch(FancyArrowPatch(a, b, arrowstyle="-|>", mutation_scale=7, lw=0.8, color=col,
                                 **kw))  # fmt: skip


def main() -> None:
    args = parse_args(__doc__)
    d = out_dir(args)
    fig = new_figure("double", 3.4)
    ax = fig.add_axes((0.0, 0.0, 1.0, 1.0))
    ax.set_xlim(0, 6)
    ax.set_ylim(0, 3.4)
    ax.axis("off")
    w, x0, yb, hb = 0.86, 0.07, 0.98, 2.0
    rows = []
    for col, title, lines, (tag, kind) in STAGES:
        x = col + x0
        _box(ax, x, yb, w, hb, "#333333", fc="#F4F4F4", lw=0.6)
        ax.text(x + w / 2, yb + hb - 0.07, title, ha="center", va="top", fontsize=7,
                fontweight="bold")  # fmt: skip
        for i, line in enumerate(lines):
            ax.text(x + 0.05, yb + hb - 0.27 - 0.15 * i, line, ha="left", va="top", fontsize=6)
        col_t = {"est": TAG, "none": NONE, "info": "#333333"}[kind]
        _box(ax, x + 0.03, yb + 0.06, w - 0.06, 0.36, col_t, ls="--" if kind == "none" else "-",
             lw=0.9 if kind == "est" else 0.6)  # fmt: skip
        ax.text(x + w / 2, yb + 0.24, tag, ha="center", va="center", fontsize=6, color=col_t)
        if col < 5:
            _arrow(ax, (x + w + 0.01, yb + 1.1), (x + 1 - 0.01, yb + 1.1))
        rows.append({"column": col, "stage": title, "contents": "; ".join(lines),
                     "error_control": tag.replace("\n", " "), "kind": kind})  # fmt: skip
    # side loop: gippa makedb -> the user's own search -> peptide table
    _box(ax, 1.07, 0.18, 1.86, 0.45, OKABE["blue"], ls="--")
    ax.text(2.0, 0.405, "gippa makedb (optional): strain FASTA from the\nBakta annotation, "
            "± six-frame ORFs, contaminants", ha="center", va="center", fontsize=6,
            color=OKABE["blue"])  # fmt: skip
    _box(ax, 0.07, 0.18, 0.86, 0.45, "#333333", ls=":")
    ax.text(0.5, 0.405, "your search engine\n(DIA-NN, Cascadia;\noutside GIPPA)", ha="center",
            va="center", fontsize=6)  # fmt: skip
    _arrow(ax, (1.5, yb - 0.01), (1.5, 0.64), col=OKABE["blue"])
    _arrow(ax, (1.06, 0.405), (0.94, 0.405), col=OKABE["blue"])
    _arrow(ax, (0.5, 0.64), (0.5, yb - 0.01))
    # key
    _box(ax, 3.55, 0.36, 0.22, 0.12, TAG, lw=0.9)
    ax.text(3.83, 0.42, "error estimate reported", va="center", fontsize=6, color=TAG)
    _box(ax, 4.85, 0.36, 0.22, 0.12, NONE, ls="--")
    ax.text(5.13, 0.42, "no estimate; flagged", va="center", fontsize=6, color=NONE)
    title = ("GIPPA: one genome and any search output in, one gene table out, with an error "
             "estimate or a flag at every step")  # fmt: skip
    ax.text(3.0, 3.3, title, ha="center", va="top", fontsize=8)
    save(fig, d, "fig1_workflow", rows)


if __name__ == "__main__":
    main()
