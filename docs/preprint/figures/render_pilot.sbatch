#!/bin/bash
#SBATCH -J gippa-preprint-pilot-render
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH -o /home/classes/bios/270/khoa/giap/local/preprint/logs/%x-%j.out
# Re-render the GIPPA report plots of the Nissle pilot arms (K-12 reference, strain database) with
# the code of this worktree, from the 1.1.1 detection rerun (Slurm 2002006). The rerun outdirs
# symlink into the frozen pilot run, so they are copied (links resolved) to
# $GIPPA_ROOT/preprint/pilot_render/<arm>/out first; nothing in the freeze is written.
# Usage: sbatch --exclude=... docs/preprint/figures/render_pilot.sbatch WORKTREE
set -euo pipefail
WT=${1:?worktree path}
source "$WT/scripts/farmshare_env.sh"
SRC=$GIPPA_ROOT/scratch/a1-recheck/runs/fix
DST=$GIPPA_ROOT/preprint/pilot_render
PY=$UV_PROJECT_ENVIRONMENT/bin/python
export PYTHONPATH=$WT/src PYTHONDONTWRITEBYTECODE=1
for a in a1 a2; do
  rm -rf "$DST/$a"
  mkdir -p "$DST/$a"
  rsync -aL --exclude=work --exclude=plots --exclude=report "$SRC/$a/out/" "$DST/$a/out/"
  $PY -c 'from gippa.cli import cli; cli()' report "$DST/$a/out" --no-interactive \
    > "$DST/$a/report.log" 2>&1
  echo "$a report exit $?"
  ls "$DST/$a/out/plots" | head -50
done
$PY -c 'import gippa; print("gippa", gippa.__version__)'
