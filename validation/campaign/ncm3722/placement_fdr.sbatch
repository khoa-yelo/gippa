#!/bin/bash
#SBATCH --job-name=gippa-ncm-te8
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=16000M
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%j.out
# Phase 1b check 1b-5 (placement_fdr.py) at the pinned commit.
#   sbatch --exclude=... --time=04:00:00 validation/campaign/ncm3722/placement_fdr.sbatch
set -euo pipefail
export GIPPA_ROOT=/home/classes/bios/270/khoa/giap/local
PIN="$GIPPA_ROOT/worktrees/campaign-pin"
export GIPPA_REPO="$PIN"
source "$PIN/scripts/farmshare_env.sh"
export PYTHONPATH="$PIN/src:$PIN/validation/proxies"
HERE="$GIPPA_ROOT/worktrees/campaign-ncm3722/validation/campaign/ncm3722"
D="$GIPPA_ROOT/campaign/ncm3722"; mkdir -p "$D/score"
cd "$PIN"
t0=$(date +%s)
"$UV_PROJECT_ENVIRONMENT/bin/python" "$HERE/placement_fdr.py" "$D/in_n2/peptides.tsv.gz" \
  "$PIN/tests/data/genomes/k12_mg1655.fna.gz" \
  "$GIPPA_ROOT/recordings/bakta_full_circ/k12_mg1655/bakta/k12_mg1655.json" "$D/score/placement_fdr.json"
echo "wall_s=$(( $(date +%s) - t0 )) job=${SLURM_JOB_ID:-none}"
