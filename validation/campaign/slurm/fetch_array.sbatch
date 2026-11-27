#!/bin/bash
#SBATCH --job-name=gippa-camp-fetch
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=2000M
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%A_%a.out
# Campaign Phase 2 prep: array task N fetches row N of validation/campaign/fetch_tasks.tsv
# (accession + selected raw files from run_selection.tsv) into $GIPPA_ROOT/campaign/<acc>/raw/
# with fetch_runs.sbatch (PRIDE https, resumable, size + deposit sha1 verified, sha256 recorded).
#   sbatch --exclude=... --time=36:00:00 --array=1-13%4 validation/campaign/slurm/fetch_array.sbatch
set -uo pipefail
REPO="${REPO:-/home/classes/bios/270/khoa/giap/local/worktrees/campaign-phase0}"
C="$REPO/validation/campaign"
source /home/classes/bios/270/khoa/giap/scripts/farmshare_env.sh
N="${SLURM_ARRAY_TASK_ID:?array task id}"
row=$(awk -F'\t' -v n="$N" 'NR > 1 && $1 == n' "$C/fetch_tasks.tsv")
[ -n "$row" ] || { echo "no task $N" >&2; exit 2; }
ACC=$(cut -f2 <<<"$row")
case "$ACC" in PXD[0-9][0-9][0-9][0-9][0-9][0-9]) ;; *) echo "bad accession" >&2; exit 2 ;; esac
read -r -a FILES <<<"$(cut -f5 <<<"$row")"
echo "task $N: $ACC ${#FILES[@]} files"
exec bash "$C/slurm/fetch_runs.sbatch" "$C/$ACC/files.tsv" "$GIPPA_ROOT/campaign/$ACC/raw" "${FILES[@]}"
