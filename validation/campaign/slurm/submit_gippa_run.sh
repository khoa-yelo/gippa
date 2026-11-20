#!/bin/bash
# Submit one campaign `gippa run` through the pinned checkout's scripts/slurm/gippa_run.sbatch
# (standard recipe step 7: all modules, default seed, --threads 14 --memory 54.GB).
#   submit_gippa_run.sh NAME PIN OUTDIR GENOME BAKTA_JSON PEPTIDES [extra gippa run args...]
# GIPPA_REPO=PIN pins the code (gippa_run.sbatch also sets PYTHONPATH=PIN/src); GIPPA_ROOT and
# GIPPA_DB_DIR are exported so the job does not fall back to PIN/local. Prints the job ID.
set -euo pipefail
NAME="$1" PIN="$2" OUT="$3" GENOME="$4" JSON="$5" PEP="$6"; shift 6
source /home/classes/bios/270/khoa/giap/scripts/farmshare_env.sh
[ ! -e "$OUT" ] || { echo "$OUT exists" >&2; exit 2; }
export GIPPA_ROOT GIPPA_DB_DIR
GIPPA_REPO="$PIN" sbatch --parsable -J "gippa-run-$NAME" --time=2-00:00:00 \
  --exclude=barley-01,barley-02,barley-03,barley-04,rye-01,rye-02 \
  "$PIN/scripts/slurm/gippa_run.sbatch" --genome "$GENOME" --annotation "$JSON" \
  --peptides "$PEP" --db-dir "$GIPPA_DB_DIR" --threads 14 --memory 54.GB --outdir "$OUT" "$@"
