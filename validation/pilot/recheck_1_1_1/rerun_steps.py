# ruff: noqa: B905, E741, RUF005, SIM115  (copied as run for the recheck)
"""Re-execute gene_table -> detectability -> context_stats -> function of one freeze arm with the
gippa on PYTHONPATH, into DEST, then build DEST/out (original outputs symlinked, regenerated tables
overriding). Usage: rerun_steps.py ORIG_OUT DEST"""

import os
import re
import shlex
import subprocess
import sys
from pathlib import Path

orig, dest = Path(sys.argv[1]), Path(sys.argv[2])
dest.mkdir(parents=True, exist_ok=True)
trace = [l.rstrip("\n").split("\t") for l in open(orig / "logs/trace.tsv")]
hdr = trace[0]
rows = [dict(zip(hdr, r)) for r in trace[1:]]
wd = {r["tag"]: r["workdir"] for r in rows}
remap = {}
for tag in ("gene_table", "detectability", "context_stats", "function"):
    o = Path(wd[tag])
    n = dest / tag
    n.mkdir(exist_ok=True)
    for line in open(o / ".command.run"):
        m = re.search(r"mkdir -p (\S+) && ln -sfn (\S+) (\S+)$", line.strip())
        if not m:
            continue
        d, src, link = m.groups()
        for k, v in remap.items():
            if src.startswith(k + "/"):
                src = v + src[len(k) :]
        (n / d).mkdir(parents=True, exist_ok=True)
        lp = n / link
        if lp.is_symlink() or lp.exists():
            lp.unlink()
        lp.symlink_to(src)
        assert Path(src).exists(), src
    cmd = next(l for l in open(o / ".command.sh") if l.startswith("gippa _step"))
    argv = [sys.executable, "-c", "from gippa.cli import cli; cli()"] + shlex.split(cmd)[1:]
    env = dict(
        os.environ,
        GIPPA_FGSEA_PREFIX="/home/classes/bios/270/khoa/giap/local/micromamba/envs/gippa-fgsea",
        TMPDIR=str(n / "tmp"),
        PYTHONNOUSERSITE="1",
        PYTHONDONTWRITEBYTECODE="1",
    )
    (n / "tmp").mkdir(exist_ok=True)
    print(tag, "running", flush=True)
    r = subprocess.run(argv, cwd=n, env=env, capture_output=True, text=True)
    open(n / "step.log", "w").write(r.stdout + r.stderr)
    print(tag, "exit", r.returncode, flush=True)
    if r.returncode:
        sys.exit(r.returncode)
    remap[str(o)] = str(n)
# fake output dir
out = dest / "out"
out.mkdir(exist_ok=True)
for e in orig.iterdir():
    if e.name in ("tables", "work"):
        continue
    t = out / e.name
    if not (t.is_symlink() or t.exists()):
        t.symlink_to(e)
tb = out / "tables"
tb.mkdir(exist_ok=True)
new = {
    p.name: p for tag in ("context_stats", "function") for p in (dest / tag / "tables").glob("*")
}
for e in (orig / "tables").iterdir():
    t = tb / e.name
    if t.is_symlink() or t.exists():
        t.unlink()
    t.symlink_to(new.get(e.name, e))
print("overridden tables:", sorted(new), flush=True)
