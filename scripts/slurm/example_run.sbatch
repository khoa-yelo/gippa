#!/bin/bash
#SBATCH --job-name=gippa-example
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%j.out
# Quick-start example timing (feat/example-dataset): convert one route of a bundle built by
# scripts/example/build_example.py, then one `gippa run` with real tools on the local executor
# inside this allocation (profile farmshare_local,apptainer), so the wall time is what a user
# sees on one machine with $SLURM_CPUS_PER_TASK CPUs.
#
# Usage: sbatch --exclude=... [-c N] scripts/slurm/example_run.sbatch NAME ROUTE MODULES [extra gippa run args]
#   NAME     bundle name under $GIPPA_ROOT/example (e.g. nissle_r100k)
#   ROUTE    A (diann) | B (proteins+abundance) | C (cascadia) | D (generic)
#   MODULES  core | default (all four modules)
# Writes $GIPPA_ROOT/example/runs/<NAME>_<ROUTE>_<MODULES>_c<CPUS>_<jobid>/ with out/, timing.txt.
set -euo pipefail
export GIPPA_ROOT="${GIPPA_ROOT:-/home/classes/bios/270/khoa/giap/local}"
WT="${EXAMPLE_REPO:-/home/classes/bios/270/khoa/giap/local/worktrees/example-data}"
export GIPPA_REPO="$WT"
source "$WT/scripts/farmshare_env.sh"
export TMPDIR="$(mktemp -d /tmp/gippa.XXXXXX)"; trap 'rm -rf "$TMPDIR"' EXIT
export PYTHONPATH="$WT/src"
export PATH="$UV_PROJECT_ENVIRONMENT/bin:$PATH"   # tasks call `gippa _step` from PATH
GIPPA="$UV_PROJECT_ENVIRONMENT/bin/gippa"
NAME=$1; ROUTE=$2; MODULES=$3; shift 3
CPUS=${SLURM_CPUS_PER_TASK:-8}
B="$GIPPA_ROOT/example/$NAME"
R="$GIPPA_ROOT/example/runs/${NAME}_${ROUTE}_${MODULES}_c${CPUS}_${SLURM_JOB_ID:-local}"
mkdir -p "$R/inputs"
G="$B/genome/$NAME.fna.gz"
cd "$WT"
t0=$(date +%s)
case "$ROUTE" in
  A) "$GIPPA" convert --from diann "$B/A_diann/report.parquet" --out "$R/inputs/peptides.tsv" \
       --abundance-out "$R/inputs/abundance.tsv" --search-db same
     EV=(--peptides "$R/inputs/peptides.tsv") ;;
  B) EV=(--proteins "$B/B_proteins/proteins.faa" --abundance "$B/B_proteins/abundance.tsv" --abundance-type maxlfq) ;;
  C) "$GIPPA" convert --from cascadia "$B"/C_cascadia/WT_{1,2,3}.ssl --out "$R/inputs/peptides.tsv"
     EV=(--peptides "$R/inputs/peptides.tsv") ;;
  D) "$GIPPA" convert --from generic "$B/D_generic/peptides.tsv" --out "$R/inputs/peptides.tsv" \
       --source database --search-db same
     EV=(--peptides "$R/inputs/peptides.tsv") ;;
  *) echo "unknown route $ROUTE" >&2; exit 1 ;;
esac
t1=$(date +%s)
MOD=(); [[ "$MODULES" == core ]] && MOD=(--modules core)
mem_gb=$(( ${SLURM_MEM_PER_NODE:-32768} / 1024 - 2 ))
set +e
"$GIPPA" run --profile farmshare_local,apptainer --genome "$G" "${EV[@]}" "${MOD[@]}" \
  --db-dir "$GIPPA_DB_DIR" --threads "$CPUS" --memory "${mem_gb}.GB" --outdir "$R/out" "$@"
rc=$?
set -e
t2=$(date +%s)
{ echo "job=${SLURM_JOB_ID:-} node=$(hostname) cpus=$CPUS route=$ROUTE modules=$MODULES name=$NAME"
  echo "convert_s=$((t1-t0)) run_s=$((t2-t1)) total_s=$((t2-t0)) exit=$rc"
  echo "commit=$(git -C "$WT" rev-parse HEAD)"; echo "extra_args=${*:-none}" | tr " " "_"; } | tee "$R/timing.txt"
exit $rc
