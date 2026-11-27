#!/bin/bash
#SBATCH --job-name=gippa-campaign-report
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=16000M
#SBATCH --time=02:00:00
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%j.out
# Standard recipe step 8: `gippa report OUT --redact` at the pinned commit, per run directory.
# Usage: sbatch report_redact.sbatch PIN OUT [OUT ...]
set -uo pipefail
source "${GIPPA_REPO:-/home/classes/bios/270/khoa/giap}/scripts/farmshare_env.sh"
PIN="$1"; shift
export TMPDIR="$(mktemp -d /tmp/gippa.XXXXXX)"; trap 'rm -rf "$TMPDIR"' EXIT
g() { PYTHONPATH="$PIN/src" "$UV_PROJECT_ENVIRONMENT/bin/python" -c 'from gippa.cli import cli; cli()' "$@"; }
rc=0
for o in "$@"; do
  s=$(date +%s); g report "$o" --redact >"$o/report_redact.log" 2>&1; r=$?
  echo "gippa report $o --redact: exit $r wall_s $(( $(date +%s) - s ))"; tail -3 "$o/report_redact.log"
  [ "$r" -eq 0 ] || rc=1
done
exit $rc
