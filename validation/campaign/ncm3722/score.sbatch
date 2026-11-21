#!/bin/bash
#SBATCH --job-name=gippa-ncm-score
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=8000M
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%j.out
# Phase 1b step 6: gippa compare N1 N2 (exit code recorded), then score.py.
set -uo pipefail
export GIPPA_ROOT=/home/classes/bios/270/khoa/giap/local
PIN="$GIPPA_ROOT/worktrees/campaign-pin"
export GIPPA_REPO="$PIN"
source "$PIN/scripts/farmshare_env.sh"
export PYTHONPATH="$PIN/src"; export PATH="$UV_PROJECT_ENVIRONMENT/bin:$PATH"
D="$GIPPA_ROOT/campaign/ncm3722"; HERE="$GIPPA_ROOT/worktrees/campaign-ncm3722/validation/campaign/ncm3722"
cd "$PIN"
rm -rf "$D/cmp_n1n2"
gippa compare "$D/out_n1" "$D/out_n2" --labels ncm3722,mg1655 --outdir "$D/cmp_n1n2" > "$D/score/compare.log" 2>&1
echo "compare_exit=$?" | tee -a "$D/score/compare.log"
tail -5 "$D/score/compare.log"
"$UV_PROJECT_ENVIRONMENT/bin/python" -I "$HERE/score.py" "$D" 20261014
