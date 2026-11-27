"""Figure 5: novel ORFs in the campaign: why GIPPA calls none at α = 0.01, and the post hoc
positive-control q tier (an analysis tier, not a GIPPA call or a change of α).

Reads validation/campaign/novel_pc/novel_pc_threshold.json and passers.tsv (route-B
novel_orfs.tsv aggregates of 17 runs; 1,581 candidates). (a) the smallest class FDR of each
(run, class) against the class size, with a 1/n guide line (not a floor of the estimator);
(b) passers of the τ_prfB tier per run,
split into positive controls and non-controls with or without a RefSeq match.
"""

from __future__ import annotations

from collections import defaultdict

import numpy as np

import _campaign as cam
from _common import OKABE, new_figure, out_dir, panel_label, parse_args, plain_log, save
from fig4_campaign import _name

# One novel-ORF class palette, shared with Fig. 6.
CLASSES = {"nterm_extension": ("N-terminal extension", OKABE["reddish_purple"], "o"),
           "intergenic": ("intergenic", OKABE["bluish_green"], "s"),
           "cds_outframe": ("CDS out of frame", OKABE["black"], "D"),
           "antisense": ("antisense", OKABE["sky_blue"], "^")}  # fmt: skip
THR = "#444444"  # threshold lines
ALPHA = 0.01


def main() -> None:
    args = parse_args(__doc__)
    d = out_dir(args)
    pc = cam.novel_pc()
    tau = pc["tau"]["tau_prfB"]
    res = pc["results"]["tau_prfB"]
    rows: list[dict] = []
    minq: dict[tuple[str, str], float] = defaultdict(lambda: float("inf"))
    for c in cam.candidates():
        k = (c["run"], c["class"])
        minq[k] = min(minq[k], float(c["q"]))

    fig = new_figure("double", 3.5)
    gs = fig.add_gridspec(1, 2, left=0.075, right=0.985, bottom=0.3, top=0.9, wspace=0.5,
                          width_ratios=[1.0, 1.15])  # fmt: skip

    # (a) smallest class FDR against class size
    ax = fig.add_subplot(gs[0, 0])
    nn = np.logspace(0.1, 3, 50)
    ax.plot(nn, 1 / nn, color="#333333", lw=0.6, ls="-")
    ax.text(150, 0.0052, "1/n (guide)", fontsize=6, color="#333333", ha="right")
    ax.axhline(ALPHA, color=THR, lw=0.7, ls="--")
    ax.text(1.7, ALPHA * 0.82, "α = 0.01", fontsize=6, color=THR, ha="left", va="top")
    ax.axhline(tau, color=THR, lw=0.7, ls=":")
    ax.text(190, tau * 1.12, f"τ_prfB = {tau:.3g}", fontsize=6, color=THR, ha="right")
    for cls, (lab, col, mk) in CLASSES.items():
        xs, ys = [], []
        for run, r in res.items():
            n = r["candidates_by_class"].get(cls, 0)
            if n == 0 or (run, cls) not in minq:
                continue
            xs.append(n)
            ys.append(minq[(run, cls)])
            rows.append({"panel": "a", "run": run, "class": cls, "n_candidates": n,
                         "min_class_fdr": minq[(run, cls)]})  # fmt: skip
        ax.scatter(xs, ys, s=12, c=col, marker=mk, lw=0, alpha=0.85, label=lab)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1.5, 200)
    ax.set_ylim(0.004, 1.5)
    plain_log(ax, "x")
    plain_log(ax, "y")
    ax.set_xlabel("candidates in the class (count)")
    ax.set_ylabel("smallest class FDR (q) in the run")
    ax.set_title("No class reaches α = 0.01", loc="left")
    ax.legend(loc="upper left", bbox_to_anchor=(-0.02, -0.17), frameon=False, fontsize=6,
              ncol=2, title="candidate class (one point per run and class)",
              title_fontsize=6)  # fmt: skip
    ax.get_legend()._legend_box.align = "left"
    panel_label(ax, "a")

    # (b) the tau_prfB tier per run
    ax = fig.add_subplot(gs[0, 1])
    det = cam.detection()  # the genome order of Fig. 4
    runs = [r["run"] for r in det]
    names = [_name(r["dataset"]) for r in det]
    y = np.arange(len(runs))[::-1]
    ctl = np.array([res[r]["n_controls_passing"] for r in runs], dtype=float)
    match = np.array([res[r]["noncontrol_refseq_match"][0] for r in runs], dtype=float)
    nonc = np.array([res[r]["n_noncontrol_passing"] for r in runs], dtype=float)
    cand = np.array([res[r]["n_candidates"] for r in runs], dtype=float)
    ax.barh(y, ctl, color=OKABE["black"], height=0.7, label="positive control")
    ax.barh(y, match, left=ctl, color="#666666", height=0.7,
            label="non-control, RefSeq match")  # fmt: skip
    ax.barh(y, nonc - match, left=ctl + match, color="#CCCCCC", height=0.7,
            label="non-control, no RefSeq match")  # fmt: skip
    for yi, r, c in zip(y, runs, cand, strict=True):
        tot = res[r]["n_pass"]
        ax.text(tot + 2, yi, f"{tot} / {int(c)}", fontsize=6, va="center", color="#555555")
        rows.append(
            {
                "panel": "b",
                "run": r,
                "dataset": res[r]["dataset"],
                "n_candidates": int(c),
                "n_pass": tot,
                "controls": int(ctl[runs.index(r)]),
                "noncontrol_refseq_match": int(match[runs.index(r)]),
                "noncontrol": int(nonc[runs.index(r)]),
                "expected_false_gippa_class_fdr": res[r]["gippa_class_fdr_expected_false"],
            }
        )
    ax.set_yticks(y)
    ax.set_yticklabels(names, fontsize=6)
    ax.set_xlim(0, 160)
    ax.set_xlabel(f"candidates with class FDR ≤ τ_prfB = {tau:.3g} (count)")
    tot = pc["results"]["tau_prfB_total"]
    ax.set_title(f"Post hoc tier: {tot['n_pass']} of {tot['n_candidates']:,} pass; "
                 "GIPPA calls 0", loc="left")  # fmt: skip
    ax.legend(loc="upper left", bbox_to_anchor=(-0.02, -0.17), frameon=False, fontsize=6,
              ncol=1, title="bar labels: passing / all candidates in the run",
              title_fontsize=6)  # fmt: skip
    ax.get_legend()._legend_box.align = "left"
    panel_label(ax, "b")

    save(fig, d, "fig5_novel_orfs", rows)


if __name__ == "__main__":
    main()
