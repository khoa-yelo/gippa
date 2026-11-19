#!/bin/bash
#SBATCH --job-name=gippa-campaign-demo
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=16000M
#SBATCH --time=02:00:00
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%j.out
# `gippa demo` and `gippa example fetch` once at the pinned commit (handover plan, Phase 1), so
# the Phase 1 report shows every GIPPA command. Usage: sbatch demo_example.sbatch PIN_WORKTREE OUT_DIR
set -uo pipefail
source "${GIPPA_REPO:-/home/classes/bios/270/khoa/giap}/scripts/farmshare_env.sh"
PIN="$1" OUT="$2"
export TMPDIR="$(mktemp -d /tmp/gippa.XXXXXX)"; trap 'rm -rf "$TMPDIR"' EXIT
mkdir -p "$OUT"; rm -rf "$OUT/demo" "$OUT/example"
echo "commit $(git -C "$PIN" rev-parse HEAD) dirty=$(git -C "$PIN" status --porcelain | wc -l)"
g() { PYTHONPATH="$PIN/src" "$UV_PROJECT_ENVIRONMENT/bin/python" -c 'from gippa.cli import cli; cli()' "$@"; }
cd "$OUT"
for c in "demo --outdir $OUT/demo" "example fetch --outdir $OUT/example"; do
  s=$(date +%s); g $c >"$OUT/${c%% *}.log" 2>&1; rc=$?
  echo "gippa $c: exit $rc wall_s $(( $(date +%s) - s ))"; tail -25 "$OUT/${c%% *}.log"
done
