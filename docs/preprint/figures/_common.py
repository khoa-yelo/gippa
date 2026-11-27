"""Shared helpers for the preprint figure scripts (docs/preprint/figures/).

Every figure script reads committed evidence files or run directories only, never writes into the
repository, and writes PNG (300 dpi), SVG, PDF and the plotted data as TSV to
``$GIPPA_ROOT/preprint/figures/<round>/``. Style: ``gippa.viz.style`` (GIPPA_SPEC.md §10.1):
DejaVu Sans, 7-pt text, Okabe-Ito categories, 3.5-in or 7.2-in widths. A deuteranopia twin
(``cvd_<name>.png``) is written next to each PNG for the colour-blind check.
"""

from __future__ import annotations

import argparse
import csv
import io
import os
from collections.abc import Sequence
from pathlib import Path

import numpy as np
from matplotlib.figure import Figure

from gippa.io.plot_files import read_png, stable_svg_ids, write_png
from gippa.viz import style
from gippa.viz.checklist import simulate_deuteranopia

REPO = Path(__file__).resolve().parents[3]
OKABE = style.OKABE_ITO
#: worktrees that hold the campaign evidence until it reaches main
WORKTREES = REPO.parent if REPO.parent.name == "worktrees" else REPO / "local" / "worktrees"


def parse_args(doc: str | None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=(doc or "").splitlines()[0] if doc else None)
    p.add_argument("--round", default="round1", help="output sub-directory (default round1)")
    p.add_argument("--outdir", type=Path, default=None, help="override the output directory")
    return p.parse_args()


def out_dir(args: argparse.Namespace) -> Path:
    if args.outdir is not None:
        d = Path(args.outdir)
    else:
        root = os.environ.get("GIPPA_ROOT")
        if not root:
            raise SystemExit("GIPPA_ROOT is not set (source scripts/farmshare_env.sh)")
        d = Path(root) / "preprint" / "figures" / args.round
    d.mkdir(parents=True, exist_ok=True)
    return d


def new_figure(width: str, height_in: float) -> Figure:
    style.apply()
    return style.new_figure(width, height_in)  # type: ignore[arg-type]


def panel_label(ax, letter: str, x: float | None = None, y: float = 1.04) -> None:
    """Bold panel letter just left of the axes' title (``x`` is ignored; kept for callers)."""
    pos = ax.get_position()
    ax.figure.text(
        pos.x0 - 0.012,
        pos.y1 + (pos.height * (y - 1.0)),
        letter,
        fontsize=9,
        fontweight="bold",
        va="bottom",
        ha="right",
    )


def save(fig: Figure, d: Path, name: str, rows: Sequence[dict] | None = None) -> list[Path]:
    """Write ``name``.png/.svg/.pdf (deterministic metadata), the data TSV and the CVD twin."""
    small = font_check(fig)
    if small:
        raise SystemExit(f"{name}: text below {MIN_PT} pt: {small[:8]}")
    paths = []
    png = d / f"{name}.png"
    fig.savefig(png, dpi=300, metadata={"Software": None})
    paths.append(png)
    buf = io.BytesIO()
    fig.savefig(buf, format="svg", metadata={"Date": None, "Creator": None})
    (d / f"{name}.svg").write_bytes(stable_svg_ids(buf.getvalue()))
    paths.append(d / f"{name}.svg")
    fig.savefig(d / f"{name}.pdf", metadata={"CreationDate": None, "Creator": None,
                                              "Producer": None})  # fmt: skip
    paths.append(d / f"{name}.pdf")
    rgba = read_png(png)
    sim = simulate_deuteranopia(rgba[..., :3])
    write_png(d / f"cvd_{name}.png", np.clip(np.concatenate([sim, rgba[..., 3:]], axis=-1), 0, 1))
    if rows:
        keys: list[str] = []
        for r in rows:
            keys += [k for k in r if k not in keys]
        with open(d / f"{name}.tsv", "w", newline="") as fh:
            w = csv.DictWriter(fh, keys, delimiter="\t", lineterminator="\n", restval="NA")
            w.writeheader()
            w.writerows(rows)
        paths.append(d / f"{name}.tsv")
    return paths


MIN_PT = 6.0


def font_check(fig: Figure) -> list[tuple[str, float]]:
    """Visible texts (titles, labels, ticks, legends, tables) smaller than ``MIN_PT``."""
    from matplotlib.text import Text

    fig.canvas.draw()
    out = []
    for t in fig.findobj(Text):
        s = t.get_text()
        if not t.get_visible() or not s.strip():
            continue
        if t.get_fontsize() < MIN_PT - 1e-9:
            out.append((s[:30], float(t.get_fontsize())))
        elif "$" in s and "^" in s:  # mathtext exponents render at about 70 % of the size
            out.append((s[:30], round(0.7 * float(t.get_fontsize()), 1)))
    return out


def plain_log(ax, axis: str = "y") -> None:
    """Log-axis tick labels as plain numbers (no mathtext exponents below 6 pt)."""
    from matplotlib.ticker import FuncFormatter, NullFormatter

    fmt = FuncFormatter(lambda v, _pos: f"{v:g}")
    target = ax.yaxis if axis == "y" else ax.xaxis
    target.set_major_formatter(fmt)
    target.set_minor_formatter(NullFormatter())


def wilson(k: int, n: int, z: float = 1.959963984540054) -> tuple[float, float, float]:
    """Share k/n with its Wilson 95 % interval."""
    if n == 0:
        return (float("nan"),) * 3
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return p, max(0.0, c - h), min(1.0, c + h)
