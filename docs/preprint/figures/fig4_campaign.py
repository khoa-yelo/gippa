"""Figure 4: the demo and discovery campaign, Phases 1 and 2 (11 datasets, 17 genomes).

Reads the committed campaign evidence through _campaign.py (PHASE2_SUMMARY.md §1, §2.3-2.5;
PXD065022/metrics.json; ecoli/metrics_phase1_fix.json). (a) route-A detection share per genome;
(b) share of route-A genes kept by the six-frame route B; (c) the TE4 RplL/RplJ ratio against its
band; (d) answer-key verdicts per dataset.
"""

from __future__ import annotations

import numpy as np

import _campaign as cam
from _common import OKABE, new_figure, out_dir, panel_label, parse_args, save, wilson


def _order(rows: list[dict]) -> list[dict]:
    return rows  # Phase 1 first, then Phase 2 in summary order, T. pallidum strains last


def main() -> None:
    args = parse_args(__doc__)
    d = out_dir(args)
    det = _order(cam.detection())
    ret = {r["run"]: r for r in cam.retention()}
    rows: list[dict] = []
    fig = new_figure("double", 5.0)
    gl = fig.add_gridspec(1, 2, left=0.15, right=0.6, bottom=0.11, top=0.89, wspace=0.08)
    gr = fig.add_gridspec(2, 1, left=0.725, right=0.985, bottom=0.11, top=0.89, hspace=0.75)
    y = np.arange(len(det))[::-1]
    labels = [r["dataset"] for r in det]

    # (a) detection share
    ax = fig.add_subplot(gl[0, 0])
    for yi, r in zip(y, det, strict=True):
        v, lo, hi = wilson(r["k"], r["n"])
        col = OKABE["orange"] if r["run"] == "ecoli_e2" else OKABE["blue"]
        ax.errorbar(v, yi, xerr=[[v - lo], [hi - v]], fmt="o", ms=3, color=col, lw=0.8,
                    capsize=1.2)  # fmt: skip
        ax.text(0.205, yi, f"{r['k']:,}/{r['n']:,}", fontsize=5, va="center", ha="left",
                color="#555555")  # fmt: skip
        rows.append({"panel": "a", **r, "share": v, "ci_low": lo, "ci_high": hi})
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=6)
    ax.set_xlim(0.2, 1.0)
    ax.set_ylim(-0.7, len(det) - 0.3)
    ax.set_xlabel("genes detected, route A\n(share of coding genes, q ≤ 0.01)")
    ax.set_title("Detection depth", loc="left")
    panel_label(ax, "a")

    # (b) route-B retention
    ax = fig.add_subplot(gl[0, 1], sharey=ax)
    for yi, r0 in zip(y, det, strict=True):
        r = ret[r0["run"]]
        if "k" in r:
            v, lo, hi = wilson(r["k"], r["n"])
        else:
            v, lo, hi = r["share"], r["ci_low"], r["ci_high"]
        col = OKABE["orange"] if r["run"] == "ecoli_e2" else OKABE["blue"]
        ax.errorbar(v, yi, xerr=[[v - lo], [hi - v]], fmt="o", ms=3, color=col, lw=0.8,
                    capsize=1.2)  # fmt: skip
        rows.append({"panel": "b", **r, "share": v, "ci_low": lo, "ci_high": hi})
    ax.axvline(0.95, color=OKABE["vermillion"], lw=0.6, ls="--")
    ax.text(0.951, len(det) - 0.6, "prior ≥ 0.95", fontsize=5, color=OKABE["vermillion"])
    ax.tick_params(labelleft=False)
    ax.set_xlim(0.93, 1.004)
    ax.set_xlabel("route-A genes kept by route B\n(share, Wilson 95 % CI)")
    ax.set_title("Six-frame search keeps genes", loc="left")
    panel_label(ax, "b")

    # (c) TE4
    ax = fig.add_subplot(gr[0, 0])
    te = cam.te4()
    x = np.arange(len(te))
    ax.axhspan(1, 3, color="#E8E8E8", lw=0)
    for xi, r in zip(x, te, strict=True):
        col = OKABE["vermillion"] if r["verdict"] == "fail" else OKABE["blue"]
        ax.plot(xi, r["rpl_ibaq"], "o", color=col, ms=3.5)
        ax.plot(xi, r["rpl_top3"], "o", mfc="white", mec=col, ms=3.2, mew=0.7)
        rows.append({"panel": "c", **r})
    ax.set_xticks(x)
    ax.set_xticklabels([r["genome"].replace("T. pallidum ", "Tpa ").replace(" (Phase 1)", "")
                        for r in te], rotation=60, ha="right", fontsize=5)  # fmt: skip
    ax.set_ylabel("log2(RplL/RplJ)")
    ax.set_ylim(-1.6, 3.8)
    ax.text(len(te) - 0.5, 2.9, "TE4 band [1, 3]", fontsize=5, ha="right", va="top",
            color="#555555")  # fmt: skip
    ax.plot([], [], "o", color="#555555", ms=3, label="iBAQ (verdict)")
    ax.plot([], [], "o", mfc="white", mec="#555555", ms=3, label="Top3")
    ax.legend(loc="lower left", frameon=False, fontsize=5, ncol=2)
    ax.set_title("TE4 fails in 6 of 15 genomes", loc="left")
    panel_label(ax, "c")

    # (d) verdicts
    ax = fig.add_subplot(gr[1, 0])
    ver = cam.verdicts()
    ver = [ver[-1], *ver[:-1]]
    names = [v["dataset"].replace("Phase 1 E. coli, for reference", "E. coli (Ph. 1)")
             .replace(" (7 strains)", " (7)").replace("thetaiotaomicron", "theta")
             .replace(" CC9311", "") for v in ver]  # fmt: skip
    yv = np.arange(len(ver))[::-1]
    left = np.zeros(len(ver))
    for key, lab, col in [("pass", "pass", OKABE["blue"]), ("fail", "fail", OKABE["vermillion"]),
                          ("nt_nr", "not testable / run", "#BBBBBB"),
                          ("reported", "reported only", OKABE["sky_blue"])]:  # fmt: skip
        vals = np.array([v[key] for v in ver], dtype=float)
        ax.barh(yv, vals, left=left, color=col, height=0.7, label=lab)
        left += vals
    for v in ver:
        rows.append({"panel": "d", **v})
    ax.set_yticks(yv)
    ax.set_yticklabels(names, fontsize=5.5)
    ax.set_xlabel("answer-key checks (count)")
    ax.legend(loc="lower right", frameon=False, fontsize=5, bbox_to_anchor=(1.0, 1.0), ncol=4,
              handlelength=0.8, columnspacing=0.6)  # fmt: skip
    ax.set_title("Pre-registered checks", loc="left", pad=14)
    panel_label(ax, "d", y=1.13)

    fig.suptitle("Campaign: GIPPA 1.1.1 code on 11 public DIA datasets (17 genomes); "
                 "every route-A peptide maps exactly in every genome", x=0.01, ha="left",
                 y=0.995)  # fmt: skip
    save(fig, d, "fig4_campaign", rows)


if __name__ == "__main__":
    main()
