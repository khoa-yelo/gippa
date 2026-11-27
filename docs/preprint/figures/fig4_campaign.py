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

ECOLI_COL = OKABE["black"]


def _name(s: str) -> str:
    """One name per genome in every panel."""
    s = s.replace("E. coli K-12 (Ph. 1)", "E. coli K-12").replace(
        "E. coli (Phase 1)", "E. coli K-12"
    )
    s = s.replace("Phase 1 E. coli, for reference", "E. coli K-12")
    if s == "B. theta":
        s = "B. thetaiotaomicron"
    return s.replace(" CC9311", "")


def main() -> None:
    args = parse_args(__doc__)
    d = out_dir(args)
    det = cam.detection()
    ret = {r["run"]: r for r in cam.retention()}
    te = {_name(r["genome"]): r for r in cam.te4()}
    rows: list[dict] = []
    fig = new_figure("double", 6.0)
    gt = fig.add_gridspec(1, 3, left=0.175, right=0.985, bottom=0.38, top=0.92, wspace=0.12,
                          width_ratios=[1.15, 1.0, 1.0])  # fmt: skip
    y = np.arange(len(det))[::-1]
    labels = [_name(r["dataset"]) for r in det]

    # (a) detection share
    ax = fig.add_subplot(gt[0, 0])
    for yi, r in zip(y, det, strict=True):
        v, lo, hi = wilson(r["k"], r["n"])
        col = ECOLI_COL if r["run"] == "ecoli_e2" else OKABE["blue"]
        ax.errorbar(v, yi, xerr=[[v - lo], [hi - v]], fmt="o", ms=3, color=col, lw=0.8,
                    capsize=1.2)  # fmt: skip
        tx = hi + 0.015 if v < 0.75 else lo - 0.015
        ax.text(tx, yi, f"{r['k']:,} / {r['n']:,}", fontsize=6, va="center",
                ha="left" if v < 0.75 else "right",
                color="#555555")  # fmt: skip
        rows.append({"panel": "a", **r, "share": v, "ci_low": lo, "ci_high": hi})
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=6)
    ax.set_xlim(0.2, 1.0)
    ax.set_ylim(-0.7, len(det) - 0.3)
    ax.set_xlabel("route-A genes detected\n(share of coding genes, gene q ≤ 0.01)")
    ax.set_title("Detection depth", loc="left")
    panel_label(ax, "a")

    # (b) route-B retention
    ax = fig.add_subplot(gt[0, 1], sharey=ax)
    for yi, r0 in zip(y, det, strict=True):
        r = ret[r0["run"]]
        if "k" in r:
            v, lo, hi = wilson(r["k"], r["n"])
        else:
            v, lo, hi = r["share"], r["ci_low"], r["ci_high"]
        col = ECOLI_COL if r["run"] == "ecoli_e2" else OKABE["blue"]
        ax.errorbar(v, yi, xerr=[[v - lo], [hi - v]], fmt="o", ms=3, color=col, lw=0.8,
                    capsize=1.2)  # fmt: skip
        rows.append({"panel": "b", **r, "share": v, "ci_low": lo, "ci_high": hi})
    ax.axvline(0.95, color=OKABE["vermillion"], lw=0.6, ls="--")
    ax.text(0.952, -0.55, "prior\n≥ 0.95", fontsize=6, color=OKABE["vermillion"], va="bottom")
    ax.tick_params(labelleft=False)
    ax.set_xlim(0.93, 1.004)
    ax.set_xlabel("route-A genes also detected\nby route B (share)")
    ax.set_title("Six-frame route keeps genes", loc="left")
    panel_label(ax, "b")

    # (c) TE4, same genome order
    ax = fig.add_subplot(gt[0, 2], sharey=ax)
    ax.axvspan(1, 3, color="#E8E8E8", lw=0)
    n_fail = n_test = 0
    for yi, lab in zip(y, labels, strict=True):
        r = te.get(lab)
        if r is None:
            ax.text(-1.45, yi, "not testable (archaeal stalk)", fontsize=6, va="center",
                    color="#777777")  # fmt: skip
            continue
        n_test += 1
        n_fail += r["verdict"] == "fail"
        col = OKABE["vermillion"] if r["verdict"] == "fail" else OKABE["blue"]
        mk = "^" if r["verdict"] == "fail" else "o"
        ax.plot(r["rpl_ibaq"], yi, mk, color=col, ms=3.5)
        ax.plot(r["rpl_top3"], yi, mk, mfc="white", mec=col, ms=3.3, mew=0.7)
        rows.append({"panel": "c", "genome": lab, **{k: v for k, v in r.items() if k != "genome"}})
    ax.tick_params(labelleft=False)
    ax.set_xlim(-1.6, 3.8)
    ax.set_xlabel("TE4: log2(RplL/RplJ)\n(band [1, 3], grey)")
    ax.plot([], [], "o", color="#555555", ms=3, label="iBAQ (verdict)")
    ax.plot([], [], "o", mfc="white", mec="#555555", ms=3, label="Top3")
    ax.plot([], [], "^", color=OKABE["vermillion"], ms=3, label="fail")
    ax.legend(loc="upper left", bbox_to_anchor=(0.0, -0.16), frameon=False, fontsize=6, ncol=3,
              handletextpad=0.2, columnspacing=0.6)  # fmt: skip
    ax.set_title(f"TE4 fails in {n_fail} of {n_test} genomes", loc="left")
    panel_label(ax, "c")

    # (d) verdicts per dataset, same order (T. pallidum strains as one dataset)
    ax = fig.add_axes((0.175, 0.055, 0.5, 0.2))
    ver = {_name(v["dataset"]): v for v in cam.verdicts()}
    order = [n for n in labels if not n.startswith("T. pallidum")] + ["T. pallidum (7 strains)"]
    vv = [ver[n] for n in order]
    yv = np.arange(len(vv))[::-1]
    left = np.zeros(len(vv))
    for key, lab, col in [("pass", "pass", OKABE["blue"]), ("fail", "fail", OKABE["vermillion"]),
                          ("nt_nr", "not testable or not run", "#BBBBBB"),
                          ("reported", "reported, no threshold", OKABE["sky_blue"])]:  # fmt: skip
        vals = np.array([v[key] for v in vv], dtype=float)
        ax.barh(yv, vals, left=left, color=col, height=0.7, label=lab)
        left += vals
    for n, v in zip(order, vv, strict=True):
        rows.append({"panel": "d", **v, "dataset": n})
    ax.set_yticks(yv)
    ax.set_yticklabels(order, fontsize=6)
    ax.set_xlabel("pre-registered answer-key checks (count)")
    ax.set_xticks(np.arange(0, 19, 2))
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, fontsize=6)
    ax.set_title("Answer-key checks per dataset", loc="left")
    panel_label(ax, "d")
    save(fig, d, "fig4_campaign", rows)


if __name__ == "__main__":
    main()
