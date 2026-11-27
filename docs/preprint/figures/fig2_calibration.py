"""Figure 2: simulation calibration of GIPPA 1.1.1 (full tier, d2bd17e4, Slurm array 2022338).

Reads only the harness records of the release-lineage calibration_full run
(``$GIPPA_ROOT/logs/calibration/full-d2bd17e42dcf-20261125-073329``: ``lanes/*/*.jsonl`` and the
``*.merged.jsonl`` files of the sharded suites). (a) QQ plots of null p-values; (b) every gated
level / FDR / FWER record against its bound, by status; (c) observed error against nominal level
for representative streams.
"""

from __future__ import annotations

import glob
import json
import math
import os
from pathlib import Path

import numpy as np
from scipy import stats

from _common import OKABE, new_figure, out_dir, panel_label, parse_args, plain_log, save

RUN = "full-d2bd17e42dcf-20261125-073329"
CRITERIA = {"level", "fdr", "fwer", "fdr_counts"}
LEVELS = (0.01, 0.05, 0.10)
#: one colour per method in panels a and c; line style separates streams of one method
METHOD = {"region": OKABE["blue"], "discordance": OKABE["bluish_green"],
          "autocorrelation": OKABE["reddish_purple"], "known": OKABE["vermillion"],
          "placement": OKABE["black"], "gene_fdr": OKABE["sky_blue"], "islands": "#777777",
          "novel": OKABE["black"]}  # fmt: skip
QQ = [  # (file under lanes/, stream, label, method, line style)
    ("regions-group0/regions.jsonl", "sim_expression.region_effect",
     "region test, O−E", "region", "-"),
    ("regions-group3/regions.jsonl", "sim_expression.ar1_lag1.rho0.9",
     "region test, AR(1) ρ = 0.9", "region", "--"),
    ("discordance-group1/discordance.jsonl", "discordance.gate.operon_null_1800",
     "operon discordance, gate null", "discordance", "-"),
    ("autocorrelation-group2/autocorrelation.jsonl", "sim_scenarios.tu_null.between_tu",
     "autocorrelation, TU null", "autocorrelation", ":"),
    ("discordance-group0/discordance.jsonl", "discordance.fallback_a.operon_null",
     "discordance, distance-model TUs\n(known failure, conservative)", "known", "-."),
]  # fmt: skip
REPS = [  # (suite, stream, criterion, label, method, line style, marker)
    ("placement", "placement.readthrough_stressed.shuffle", "fdr", "placement FDR (read-through)",
     "placement", "-", "s"),
    ("gene_fdr", "gene_fdr.picked.pi0.05", "fdr", "gene FDR, picked", "gene_fdr", "-", "o"),
    ("regions", "regions.region_tests.null_full_rho0.0.chr1_pooled.hl_shift", "level",
     "region test, HL shift (pooled)", "region", ":", "v"),
    ("novel_orfs", "novel_orfs.class_fdr.picked", "fdr", "novel-ORF class FDR (database)",
     "novel", "--", "^"),
]  # fmt: skip
STATUS_STYLE = {
    "pass": ("pass", OKABE["blue"], "o"),
    "unconfirmed": ("unconfirmed (over in 1 of 2 blocks)", OKABE["black"], "D"),
    "fail_flagged": ("fail, recorded only or flagged", OKABE["reddish_purple"], "s"),
    "fail_xfail": ("fail, known failure", OKABE["vermillion"], "^"),
}


def cal_dir() -> Path:
    root = os.environ.get("GIPPA_ROOT")
    if not root:
        raise SystemExit("GIPPA_ROOT is not set (source scripts/farmshare_env.sh)")
    return Path(root) / "logs" / "calibration" / RUN


def load(c: Path) -> list[dict]:
    """Last block of every (file, stream, criterion, level) record; lanes plus merged files."""
    files = sorted(glob.glob(str(c / "lanes" / "*" / "*.jsonl")))
    files += sorted(glob.glob(str(c / "*.merged.jsonl")))
    last: dict[tuple, dict] = {}
    for f in files:
        rel = os.path.relpath(f, c)
        with open(f) as fh:
            lines = fh.readlines()
        for line in lines:
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue
            if r.get("status") is None or "stream" not in r:
                continue
            r["_file"] = rel
            k = (rel, r["stream"], r.get("criterion"), r.get("level"))
            if k not in last or (r.get("block") or 0) >= (last[k].get("block") or 0):
                last[k] = r
    return list(last.values())


def status_key(r: dict) -> str | None:
    s = r["status"]
    if s == "pass":
        return "pass"
    if s == "unconfirmed":
        return "unconfirmed"
    if s == "fail":
        return "fail_xfail" if r.get("xfail") else "fail_flagged"
    return None  # vacuous, recorded


def main() -> None:
    args = parse_args(__doc__)
    d = out_dir(args)
    c = cal_dir()
    recs = load(c)
    rows: list[dict] = []
    fig = new_figure("double", 3.8)
    gs = fig.add_gridspec(1, 3, left=0.065, right=0.985, bottom=0.45, top=0.93, wspace=0.42)

    # (a) QQ of null p-values, -log10 scale, with a pointwise 95 % band for n = 1,000
    ax = fig.add_subplot(gs[0, 0])
    lim = 3.2
    n_band = 1000
    i = np.arange(1, n_band + 1)
    lo = -np.log10(stats.beta.ppf(0.975, i, n_band - i + 1))
    hi = -np.log10(stats.beta.ppf(0.025, i, n_band - i + 1))
    ex = -np.log10(i / (n_band + 1))
    ax.fill_between(ex, lo, hi, color="#DDDDDD", lw=0, label="pointwise 95 % band")
    ax.plot([0, lim], [0, lim], color="#333333", lw=0.6)
    for rel, stream, label, meth, ls in QQ:
        rec = next(r for r in recs if r["_file"] == f"lanes/{rel}" and r["stream"] == stream
                   and r.get("criterion") == "uniform")  # fmt: skip
        p = np.sort(np.asarray(rec["values"], dtype=float))
        n = p.size
        e = -np.log10(np.arange(1, n + 1) / (n + 1))
        o = -np.log10(np.clip(p, 1e-300, 1))
        ax.plot(e, np.minimum(o, lim), color=METHOD[meth], ls=ls, lw=1.0, label=label)
        rows.append({"panel": "a", "file": rec["_file"], "stream": stream, "n": n,
                     "ks_p": rec["observed"], "status": rec["status"]})  # fmt: skip
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_xlabel("expected −log10 p (uniform)")
    ax.set_ylabel("observed −log10 p")
    ax.set_title("Null p-values (QQ)", loc="left")
    ax.legend(loc="upper left", bbox_to_anchor=(-0.2, -0.2), frameon=False, fontsize=6,
              handlelength=2.0, title="null p-values, n = 1,000 each",
              title_fontsize=6)  # fmt: skip
    ax.get_legend()._legend_box.align = "left"
    panel_label(ax, "a", x=-0.22)

    # (b) every gated record: observed / nominal against bound / nominal
    ax = fig.add_subplot(gs[0, 1])
    floor = 0.02
    counts: dict[str, int] = {}
    for key in STATUS_STYLE:
        xs, ys = [], []
        for r in recs:
            if r.get("criterion") not in CRITERIA or status_key(r) != key:
                continue
            lev = r.get("level")
            if not lev or not any(math.isclose(lev, a, rel_tol=1e-9) for a in LEVELS):
                continue
            if r.get("observed") is None or r.get("bound") is None:
                continue
            xs.append(r["bound"] / lev)
            ys.append(max(r["observed"] / lev, floor))
            rows.append({"panel": "b", "file": r["_file"], "stream": r["stream"],
                         "criterion": r["criterion"], "level": lev, "observed": r["observed"],
                         "bound": r["bound"], "R": r.get("R"), "m": r.get("m"),
                         "status": key})  # fmt: skip
        counts[key] = len(xs)
        lab, col, mk = STATUS_STYLE[key]
        z = 3 if key == "pass" else 4
        ax.scatter(xs, ys, s=6 if key == "pass" else 12, c=col, marker=mk, lw=0, alpha=0.6
                   if key == "pass" else 0.9, label=f"{lab} ({len(xs):,})", zorder=z)  # fmt: skip
    ax.set_yscale("log")
    plain_log(ax)
    xs = np.linspace(0.8, 3.4, 200)
    ax.plot(xs, xs, color="#333333", lw=0.6)
    ax.text(3.15, 3.6, "observed = bound", fontsize=6, ha="right", va="bottom", color="#333333")
    ax.axhline(1.0, color="#999999", lw=0.5, ls=":")
    ax.text(3.15, 0.85, "observed = nominal", fontsize=6, ha="right", va="top", color="#777777")
    ax.set_xlim(0.9, 3.2)
    ax.set_ylim(floor * 0.8, 400)
    ax.set_xlabel("bound / nominal level")
    ax.set_ylabel("observed / nominal level")
    ax.set_title("Above the line: observed > bound", loc="left")
    total = sum(counts.values())
    ax.legend(loc="upper left", bbox_to_anchor=(-0.2, -0.2), frameon=False, fontsize=6,
              handletextpad=0.3, title=f"{total:,} records at α or q = 0.01, 0.05,\n"
              f"0.10; zero observed is drawn at {floor}.\nNot shown: GSEA arm null_b2, "
              "checked\nagainst 1.5α; 0.083 and 0.137 at α = 0.05\nand 0.10 exceed α + 2 SE",
              title_fontsize=6)  # fmt: skip
    ax.get_legend()._legend_box.align = "left"
    panel_label(ax, "b", x=-0.22)

    # (c) observed error against nominal level for representative streams
    ax = fig.add_subplot(gs[0, 2])
    ax.plot([0, 0.11], [0, 0.11], color="#333333", lw=0.6)
    for j, (suite, stream, crit, label, meth, ls, mk) in enumerate(REPS):
        pts = sorted(
            (r["level"], r["observed"], r["bound"])
            for r in recs
            if r.get("suite") == suite
            and r["stream"] == stream
            and r.get("criterion") == crit
            and r.get("level") is not None
            and any(math.isclose(r["level"], a, rel_tol=1e-9) for a in LEVELS)
        )
        if not pts:
            raise SystemExit(f"missing stream {stream}")
        lv, ob, bd = (np.array(v) for v in zip(*pts, strict=True))
        off = (j - 1.5) * 0.0012
        ax.plot(lv + off, ob, color=METHOD[meth], marker=mk, ms=2.6, lw=0.8, ls=ls, label=label)
        ax.scatter(lv + off, bd, color=METHOD[meth], marker="_", s=18, lw=0.8)
        for a, o, b in pts:
            rows.append({"panel": "c", "suite": suite, "stream": stream, "criterion": crit,
                         "level": a, "observed": o, "bound": b})  # fmt: skip
    ax.set_xlim(0, 0.11)
    ax.set_ylim(0, 0.13)
    ax.set_xticks(LEVELS)
    ax.set_xlabel("nominal level (α or q)")
    ax.set_ylabel("observed rate or mean FDP")
    ax.set_title("Observed vs nominal error", loc="left")
    ax.legend(loc="upper left", bbox_to_anchor=(-0.2, -0.2), frameon=False, fontsize=6,
              handlelength=2.0, title="short bar: the stream's bound (nominal + 2 SE)",
              title_fontsize=6)  # fmt: skip
    ax.get_legend()._legend_box.align = "left"
    panel_label(ax, "c", x=-0.22)

    save(fig, d, "fig2_calibration", rows)
    print(json.dumps(counts))


if __name__ == "__main__":
    main()
