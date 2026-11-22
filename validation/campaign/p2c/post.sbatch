#!/bin/bash
#SBATCH --job-name=gippa-camp-p2c-post
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=16000M
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%j.out
# Standard recipe step 8 for one dataset: gippa compare out_a out_b --labels genes,sixframe, then
# gippa report --redact on both run directories, at the pinned commit (PIN, default rel111-pin).
#   sbatch -J gippa-camp-p2c-post-<ds> --time=02:00:00 --exclude=... post.sbatch DS
set -uo pipefail
DS="${1:?DS}"
export GIPPA_ROOT=/home/classes/bios/270/khoa/giap/local
PIN="${PIN:-$GIPPA_ROOT/worktrees/rel111-pin}"
source "$PIN/scripts/farmshare_env.sh"
export PYTHONPATH="$PIN/src"
export TMPDIR="$(mktemp -d /tmp/gippa.XXXXXX)"; trap 'rm -rf "$TMPDIR"' EXIT
GIPPA="$UV_PROJECT_ENVIRONMENT/bin/gippa"
D="$GIPPA_ROOT/campaign/$DS"
rc=0
rm -rf "$D/cmp_ab"
s=$(date +%s); "$GIPPA" compare "$D/out_a" "$D/out_b" --labels genes,sixframe --outdir "$D/cmp_ab" > "$D/cmp_ab.log" 2>&1; r=$?
echo "compare exit=$r wall_s=$(( $(date +%s) - s ))" | tee "$D/post_runtime.txt"; [ "$r" -eq 0 ] || rc=1
for o in out_a out_b; do
  ( cd "$D/$o/tables" && find . -type f -print0 | sort -z | xargs -0 sha256sum ) > "$D/$o.tables.pre_report.sha256"
  s=$(date +%s); "$GIPPA" report "$D/$o" --redact > "$D/$o/report_redact.log" 2>&1; r=$?
  ( cd "$D/$o/tables" && sha256sum -c --quiet "$D/$o.tables.pre_report.sha256" ) && same=yes || same=no
  echo "report $o --redact exit=$r wall_s=$(( $(date +%s) - s )) tables_unchanged=$same" | tee -a "$D/post_runtime.txt"
  [ "$r" -eq 0 ] || rc=1
done
echo "job=${SLURM_JOB_ID:-none} node=$(hostname -s)" >> "$D/post_runtime.txt"
exit $rc
