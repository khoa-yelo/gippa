"""Figure 3: the E. coli Nissle 1917 pilot (PXD078459), arms A1 (K-12 reference) and A2 (strain).

Reads only committed aggregates: validation/pilot/recheck_1_1_1/metrics_1_1_1.json (detection,
H1; 1.1.1 detection code on the frozen pilot inputs, Slurm 2002006) and
validation/pilot/metrics/metrics.json (free entrapment m8, unchanged by ADR-0078). No peptide
sequence, no per-gene value. (a) detection share by gene subset; (b) the pre-registered H1;
(c) A1 free-entrapment FDP estimates against the 1 % target.
"""

from __future__ import annotations

import json

import numpy as np

from _common import OKABE, REPO, new_figure, out_dir, panel_label, parse_args, save

ADD = REPO / "validation/pilot/recheck_1_1_1/metrics_1_1_1.json"
MET = REPO / "validation/pilot/metrics/metrics.json"
SUBSETS = [
    ("all", "all genes"),
    ("k12_shared", "K-12-shared"),
    ("strain_specific", "strain-specific"),
    ("bgc", "BGC"),
    ("prophage", "prophage"),
    ("plasmid", "plasmid"),
]
ARMS = {"A1": ("A1: K-12 proteome", OKABE["orange"]), "A2": ("A2: strain database",
                                                             OKABE["blue"])}  # fmt: skip


def main() -> None:
    args = parse_args(__doc__)
    d = out_dir(args)
    add = json.loads(ADD.read_text())
    met = json.loads(MET.read_text())
    rows: list[dict] = []
    fig = new_figure("double", 2.6)
    gs = fig.add_gridspec(1, 3, left=0.145, right=0.975, bottom=0.2, top=0.76, wspace=0.8,
                          width_ratios=[1.25, 1.0, 1.0])  # fmt: skip

    # (a) detection share by subset, Wilson 95 % CI
    ax = fig.add_subplot(gs[0, 0])
    y = np.arange(len(SUBSETS))[::-1]
    for j, (arm, (lab, col)) in enumerate(ARMS.items()):
        sub = add["m3"][arm]["subsets"]
        for yi, (key, _name) in zip(y, SUBSETS, strict=True):
            s = sub[key]
            yy = yi + (0.17 if j == 0 else -0.17)
            ax.errorbar(s["value"], yy, xerr=[[s["value"] - s["ci_low"]],
                        [s["ci_high"] - s["value"]]], fmt="o", ms=3, color=col, lw=0.8,
                        capsize=1.5, label=lab if yi == y[0] else None)  # fmt: skip
            rows.append(
                {
                    "panel": "a",
                    "arm": arm,
                    "subset": key,
                    "k": s["k"],
                    "n": s["n"],
                    "share": s["value"],
                    "ci_low": s["ci_low"],
                    "ci_high": s["ci_high"],
                }
            )
    ax.set_yticks(y)
    ax.set_yticklabels([f"{n} ({add['m3']['A2']['subsets'][k]['n']:,})" for k, n in SUBSETS])
    ax.set_xlim(-0.02, 1.0)
    ax.set_xlabel("share of genes detected (gene q ≤ 0.01)")
    ax.set_title("Strain genes need the strain database", loc="left")
    ax.set_ylim(-0.6, len(SUBSETS) + 0.1)
    ax.legend(
        loc="upper right", ncol=2, frameon=False, fontsize=5, handletextpad=0.2, columnspacing=0.8
    )
    ax.text(0.02, y[2] - 0.45, "A1: 0 by construction", fontsize=5, color=OKABE["orange"])
    panel_label(ax, "a", y=1.04)

    # (b) H1: who detects the discordant genes, and are they mostly plasmid/prophage/BGC?
    ax = fig.add_subplot(gs[0, 1])
    h1 = add["h1"]
    items = [("prop_a2_only", f"A2-only of discordant\n({h1['b']} of {h1['b'] + h1['c']})"),
             ("strain_specific_share", "strain-specific\namong A2-only"),
             ("mostly", "plasmid, prophage\nor BGC among A2-only"),
             ("mostly_base_rate", "same, all genes\n(base rate)")]  # fmt: skip
    yb = np.arange(len(items))[::-1]
    for yi, (key, _lab) in zip(yb, items, strict=True):
        s = h1[key]
        col = OKABE["vermillion"] if key == "mostly" else OKABE["blue"]
        if key == "mostly_base_rate":
            col = "#999999"
        ax.errorbar(s["value"], yi, xerr=[[s["value"] - s["ci_low"]], [s["ci_high"] - s["value"]]],
                    fmt="o", ms=3, color=col, lw=0.8, capsize=1.5)  # fmt: skip
        rows.append({"panel": "b", "quantity": key, "k": s["k"], "n": s["n"], "share": s["value"],
                     "ci_low": s["ci_low"], "ci_high": s["ci_high"]})  # fmt: skip
    ax.axvline(0.5, color=OKABE["vermillion"], lw=0.6, ls="--")
    ax.text(0.51, yb[-1] - 0.45, "'mostly':\n> 0.5", fontsize=5, color=OKABE["vermillion"])
    ax.set_yticks(yb)
    ax.set_yticklabels([lab for _, lab in items], fontsize=5.5)
    ax.set_xlim(0, 1.0)
    ax.set_ylim(-0.8, len(items) - 0.4)
    ax.set_xlabel("share (Wilson 95 % CI)")
    ax.set_title("H1: 'mostly' fails", loc="left")
    panel_label(ax, "b", y=1.04)

    # (c) A1 free entrapment
    ax = fig.add_subplot(gs[0, 2])
    m8 = met["m8"]
    lv = [("precursor", "precursor\n(assessed)"), ("peptide", "peptide\n(not assessed)")]
    for yi, (key, _lab) in zip([1, 0], lv, strict=True):
        e = m8[key]
        lo, est, ub = e["fdp_lower"], e["fdp_combined"], e["ci"]["ucb95"]
        ax.plot([lo, ub], [yi, yi], color=OKABE["blue"], lw=1.0)
        ax.plot(lo, yi, "|", color=OKABE["blue"], ms=6)
        ax.plot(ub, yi, "|", color=OKABE["blue"], ms=6)
        ax.plot(est, yi, "o", color=OKABE["blue"], ms=3.5)
        rows.append({"panel": "c", "level": key, "N_T": e["N_T"], "N_E": e["N_E"],
                     "fdp_lower": lo, "fdp_combined": est, "ucb95": ub,
                     "assessment": e["assessment"]})  # fmt: skip
    ax.axvline(0.01, color=OKABE["vermillion"], lw=0.6, ls="--")
    ax.text(0.0105, 1.45, "1 % target", fontsize=5, color=OKABE["vermillion"])
    ax.set_yticks([1, 0])
    ax.set_yticklabels([lab for _, lab in lv], fontsize=5.5)
    ax.set_ylim(-1.0, 1.7)
    ax.set_xlim(0, 0.035)
    ax.set_xlabel("A1 entrapment FDP estimate")
    ax.set_title("A1: inconclusive", loc="left")
    ax.text(
        0.0005,
        -0.85,
        "bar: lower bound to ucb95\ndot: combined estimate",
        fontsize=5,
        color="#666666",
    )
    panel_label(ax, "c", y=1.04)

    fig.suptitle("Nissle 1917 pilot (PXD078459, 3 runs): A2 detects 3,051 of 4,685 genes, "
                 "A1 2,645; GIPPA 1.1.1 detection on the frozen inputs", x=0.01, ha="left",
                 y=0.985)  # fmt: skip
    save(fig, d, "fig3_pilot", rows)


if __name__ == "__main__":
    main()
