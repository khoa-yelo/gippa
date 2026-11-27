#!/bin/bash
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=64000M
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%j.out
# Campaign Phase 2, group A: run arms one after another inside one 16-CPU allocation, so the
# group never holds more than (lanes x 16) CPUs and nothing sits pending (agent cap 32 CPUs).
#   sbatch -J gippa-camp-p2a-L1 --exclude=... --time=47:00:00 lane.sbatch "ACC TAG R [REF]" ...
# Each argument is one arm.sh call; a failed arm is logged and the lane goes on.
set -uo pipefail
HERE=/home/classes/bios/270/khoa/giap/local/worktrees/campaign-p2a/validation/campaign/p2a
rc_all=0
for spec in "$@"; do
  s=$(date +%s)
  # shellcheck disable=SC2086
  bash "$HERE/arm.sh" $spec
  rc=$?
  echo "LANE arm='$spec' exit=$rc wall_s=$(( $(date +%s) - s )) job=${SLURM_JOB_ID:-none} node=$(hostname -s)"
  [ "$rc" -eq 0 ] || rc_all=1
done
exit $rc_all
