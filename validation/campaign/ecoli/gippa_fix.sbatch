#!/bin/bash
#SBATCH --job-name=gippa-camp-efix
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=64000M
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%j.out
# Phase 1 recheck with 1.1.1 code (ADR-0078): rerun `gippa run` for arm E1 or E2 from the SAME
# converted inputs ($D/in_a, $D/in_b; DIA-NN and convert are not redone) into a new outdir
# $D/runs/out_<arm>_fix, refusing to overwrite. Same arguments as the 1.1.0 runs (all modules,
# default seed, --threads 14 --memory 54.GB), but the farmshare_local profile: every Nextflow task
# runs inside this 16-CPU allocation (as ncm3722/gippa_arm.sbatch).
#   PIN=<pinned worktree> sbatch --exclude=... --time=05:00:00 -J gippa-camp-e2fix gippa_fix.sbatch E2
set -euo pipefail
ARM="${1:?E1|E2}"
export GIPPA_ROOT=/home/classes/bios/270/khoa/giap/local
PIN="${PIN:?pinned worktree}"
export GIPPA_REPO="$PIN"
source "$PIN/scripts/farmshare_env.sh"
export PYTHONPATH="$PIN/src"
export PATH="$UV_PROJECT_ENVIRONMENT/bin:$PATH"   # tasks call `gippa _step`
export XDG_CACHE_HOME="$GIPPA_ROOT/.cache"
export TMPDIR="$(mktemp -d /tmp/gippa.XXXXXX)"; trap 'rm -rf "$TMPDIR"' EXIT
GIPPA="$UV_PROJECT_ENVIRONMENT/bin/gippa"
D="$GIPPA_ROOT/campaign/ecoli_k12"
case "$ARM" in
  E1) IN="$D/in_a" ;;
  E2) IN="$D/in_b" ;;
  *) echo "arm E1|E2" >&2; exit 2 ;;
esac
a="${ARM,,}"; OUTD="$D/runs/out_${a}_fix"
[ ! -e "$OUTD" ] || { echo "$OUTD exists; refusing to overwrite" >&2; exit 3; }
G="$D/genome.fna.gz"; J="$D/bakta/k12_mg1655.json"
cd "$PIN"
t0=$(date +%s)
set +e
"$GIPPA" run --profile farmshare_local,apptainer --genome "$G" --annotation "$J" \
  --peptides "$IN/peptides.tsv.gz" --db-dir "$GIPPA_DB_DIR" --threads 14 --memory 54.GB \
  --outdir "$OUTD"
rc=$?
t1=$(date +%s)
echo "arm=$ARM run_s=$((t1-t0)) run_exit=$rc node=$(hostname -s) job=${SLURM_JOB_ID:-none}" | tee "$D/runs/runtime_${a}_fix.txt"
sha256sum "$G" "$J" "$IN/peptides.tsv.gz" > "$D/runs/inputs_${a}_fix.sha256"
exit $rc
