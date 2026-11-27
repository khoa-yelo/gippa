"""Figure 5: novel ORFs in the campaign: why GIPPA calls none at α = 0.01, and the post hoc
positive-control q tier (an analysis tier, not a GIPPA call or a change of α).

Reads validation/campaign/novel_pc/novel_pc_threshold.json and passers.tsv (route-B
novel_orfs.tsv aggregates of 17 runs; 1,581 candidates). (a) the smallest class FDR of each
(run, class) against the class size, with the 1/n floor; (b) passers of the τ_prfB tier per run,
split into positive controls and non-controls with or without a RefSeq match.
"""

from __future__ import annotations

from collections import defaultdict

import numpy as np

import _campaign as cam
from _common import OKABE, new_figure, out_dir, panel_label, parse_args, save

CLASSES = {"nterm_extension": ("N-terminal extension", OKABE["blue"], "o"),
           "intergenic": ("intergenic", OKABE["bluish_green"], "s"),
           "cds_outframe": ("CDS out of frame", OKABE["black"], "D"),
           "antisense": ("antisense", OKABE["reddish_purple"], "^")}  # fmt: skip
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

    fig = new_figure("double", 3.7)
    gs = fig.add_gridspec(1, 2, left=0.075, right=0.985, bottom=0.3, top=0.86, wspace=0.45,
                          width_ratios=[1.0, 1.15])  # fmt: skip

    # (a) smallest class FDR against class size
    ax = fig.add_subplot(gs[0, 0])
    nn = np.logspace(0.1, 3, 50)
    ax.plot(nn, 1 / nn, color="#333333", lw=0.6, ls="-")
    ax.text(150, 0.0052, "1/n (guide)", fontsize=6, color="#333333", ha="right")
    ax.axhline(ALPHA, color=OKABE["vermillion"], lw=0.7, ls="--")
    ax.text(1100, ALPHA * 1.12, "α = 0.01", fontsize=6, color=OKABE["vermillion"], ha="right")
    ax.axhline(tau, color=OKABE["sky_blue"], lw=0.7, ls=":")
    ax.text(1100, tau * 1.12, f"τ_prfB = {tau:.3g}", fontsize=6, color="#0072B2", ha="right")
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
    ax.set_xlim(1.5, 1200)
    ax.set_ylim(0.004, 1.5)
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
    runs = list(res)
    tpa = dict(cam.TPA)
    names = [f"T. pallidum {tpa[r]}" if r in tpa else res[r]["dataset"]
             .replace("E. coli K-12 (Phase 1 E2)", "E. coli K-12").replace(" CC9311", "")
             for r in runs]  # fmt: skip
    y = np.arange(len(runs))[::-1]
    ctl = np.array([res[r]["n_controls_passing"] for r in runs], dtype=float)
    match = np.array([res[r]["noncontrol_refseq_match"][0] for r in runs], dtype=float)
    nonc = np.array([res[r]["n_noncontrol_passing"] for r in runs], dtype=float)
    cand = np.array([res[r]["n_candidates"] for r in runs], dtype=float)
    ax.barh(y, ctl, color=OKABE["vermillion"], height=0.7, label="positive control")
    ax.barh(y, match, left=ctl, color=OKABE["blue"], height=0.7,
            label="non-control, RefSeq match")  # fmt: skip
    ax.barh(y, nonc - match, left=ctl + match, color="#BBBBBB", height=0.7,
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
    ax.set_xlim(0, 175)
    ax.set_xlabel(f"candidates with class FDR ≤ τ_prfB = {tau:.3g} (count)")
    tot = pc["results"]["tau_prfB_total"]
    ax.set_title(f"Post hoc tier: {tot['n_pass']} of {tot['n_candidates']:,} pass; "
                 "GIPPA calls 0", loc="left")  # fmt: skip
    ax.legend(loc="upper left", bbox_to_anchor=(-0.02, -0.17), frameon=False, fontsize=6,
              ncol=1, title="bar labels: passing / all candidates in the run",
              title_fontsize=6)  # fmt: skip
    ax.get_legend()._legend_box.align = "left"
    panel_label(ax, "b")

    fig.suptitle("Novel ORFs: none called at α = 0.01 in 17 genomes; a post hoc tier shows "
                 "where controls sit", x=0.01, ha="left", y=0.985)  # fmt: skip
    save(fig, d, "fig5_novel_orfs", rows)


if __name__ == "__main__":
    main()
