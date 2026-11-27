"""Supplementary Figure S1: region tests in the Nissle 1917 pilot, strain-database arm (A2).

Redraws the region-level rows of GIPPA's tables/regions_stats.tsv from the 1.1.1 fix rerun of
the pilot ($GIPPA_ROOT/preprint/pilot_render/a2/out, a copy of the rerun's context_stats
output). Aggregates per region only: no peptide and no per-gene value. (a) detection,
observed minus expected detected genes (circular-shift null) with its interval; (b) abundance,
Hodges-Lehmann shift of log2 iBAQ against the rest of the genome. Plasmid rows are replicon
contrasts with no null and no q-value; regions with no eligible gene are left out.
"""

from __future__ import annotations

import csv
import os
from pathlib import Path

from _common import OKABE, new_figure, out_dir, panel_label, parse_args, save

THR = "#444444"  # reference lines
TYPES = {"prophage": "prophage", "bgc": "BGC", "amr": "AMR genes", "plasmid": "plasmid"}
ORDER = ["prophage", "bgc", "amr", "plasmid"]


def _table() -> Path:
    root = os.environ.get("GIPPA_ROOT")
    if not root:
        raise SystemExit("GIPPA_ROOT is not set (source scripts/farmshare_env.sh)")
    return Path(root) / "preprint/pilot_render/a2/out/tables/regions_stats.tsv"


def _f(v: str) -> float | None:
    return None if v in ("", "NA") else float(v)


def main() -> None:
    args = parse_args(__doc__)
    d = out_dir(args)
    with open(_table()) as fh:
        recs = [r for r in csv.DictReader(fh, delimiter="\t") if r["subset"] == "all"]
    regions: dict[str, dict[str, dict]] = {}
    for r in recs:
        if int(r["n_eligible"]) == 0:
            continue
        regions.setdefault(r["region_id"], {})[r["stat"]] = r
    ids = sorted(regions, key=lambda k: (ORDER.index(regions[k]["o_minus_e"]["region_type"]),
                                         k))  # fmt: skip
    rows: list[dict] = []
    labels = []
    for k in ids:
        r = regions[k]["o_minus_e"]
        num = k.split("_")[-1]
        n, ne = int(r["n_genes"]), int(r["n_eligible"])
        size = f"{n} genes" if n == ne else f"{n} genes, {ne} eligible"
        name = TYPES[r["region_type"]] + ("" if r["region_type"] == "amr" else f" {num}")
        labels.append(f"{name} ({size})")

    fig = new_figure("double", 3.6)
    gs = fig.add_gridspec(1, 2, left=0.2, right=0.87, bottom=0.2, top=0.9, wspace=0.6)
    y = list(range(len(ids)))[::-1]
    for j, (stat, xlab, title) in enumerate(
        [
            ("o_minus_e", "detected genes: observed − expected", "Detection (O − E)"),
            ("hl_shift", "log2 iBAQ shift (Hodges–Lehmann)", "Abundance (HL shift)"),
        ]
    ):
        ax = fig.add_subplot(gs[0, j])
        ax.axvline(0, color=THR, lw=0.6, ls="--")
        for yi, k in zip(y, ids, strict=True):
            r = regions[k][stat]
            est, lo, hi, q = (_f(r[c]) for c in ("effect", "ci_low", "ci_high", "q_value"))
            contrast = r["method"] == "replicon_contrast"
            col = "#888888" if contrast else OKABE["blue"]
            if est is not None and lo is not None and hi is not None:
                ax.errorbar(est, yi, xerr=[[est - lo], [hi - est]], fmt="o", ms=3, color=col,
                            mfc="white" if contrast else col, lw=0.8, capsize=1.5)  # fmt: skip
            qtxt = "not tested" if contrast else f"q = {q:.2g}"
            ax.text(1.02, yi, qtxt, transform=ax.get_yaxis_transform(), fontsize=6,
                    va="center", color="#555555")  # fmt: skip
            rows.append({"panel": "ab"[j], "region_id": k, "region_type": r["region_type"],
                         "n_genes": r["n_genes"], "n_eligible": r["n_eligible"],
                         "observed": r["observed"], "expected": r["expected"], "stat": stat,
                         "effect": est, "ci_low": lo, "ci_high": hi, "p_value": _f(r["p_value"]),
                         "q_value": q, "method": r["method"], "n_perm": r["n_perm"]})  # fmt: skip
        ax.set_yticks(y)
        ax.set_yticklabels(labels if j == 0 else [], fontsize=6)
        ax.set_ylim(-0.6, len(ids) - 0.4)
        ax.set_xlabel(xlab)
        ax.set_title(title, loc="left")
        panel_label(ax, "ab"[j])
    ax.plot([], [], "o", color=OKABE["blue"], ms=3,
            label="circular-shift test (estimate, interval)")  # fmt: skip
    ax.plot([], [], "o", color="#888888", mfc="white", ms=3, label="replicon contrast (no null)")
    ax.legend(loc="upper left", bbox_to_anchor=(-1.5, -0.13), frameon=False, fontsize=6, ncol=2)
    save(fig, d, "figS1_pilot_regions", rows)


if __name__ == "__main__":
    main()
