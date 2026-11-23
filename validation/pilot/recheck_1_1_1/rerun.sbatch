#!/bin/bash
#SBATCH -J gippa-a1recheck
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:30:00
#SBATCH -o /home/classes/bios/270/khoa/giap/local/scratch/a1-recheck/logs/%x-%j.out
set -uo pipefail
source /home/classes/bios/270/khoa/giap/scripts/farmshare_env.sh
S=/home/classes/bios/270/khoa/giap/local/scratch/a1-recheck
F=/home/classes/bios/270/khoa/giap/local/runs/pilot/freeze/fa49b8f
CTL=/home/classes/bios/270/khoa/giap/local/worktrees/pilot-fa49b8f145ea/src
PY=$UV_PROJECT_ENVIRONMENT/bin/python
export PYTHONDONTWRITEBYTECODE=1
for v in ctl fix; do
  src=$CTL; [ $v = fix ] && src=$S/src_fix
  for a in a1 a2; do
    PYTHONPATH=$src $PY $S/rerun_steps.py $F/out_$a $S/runs/$v/$a > $S/logs/$v-$a.log 2>&1 &
  done
done
wait
for v in ctl fix; do
  src=$CTL; [ $v = fix ] && src=$S/src_fix
  rm -rf $S/runs/$v/cmp
  PYTHONPATH=$src $PY -c 'from gippa.cli import cli; cli()' compare $S/runs/$v/a1/out $S/runs/$v/a2/out --outdir $S/runs/$v/cmp --labels A1,A2 > $S/logs/$v-cmp.log 2>&1
  echo "$v compare exit $?"
done
grep -h exit $S/logs/*-a?.log
