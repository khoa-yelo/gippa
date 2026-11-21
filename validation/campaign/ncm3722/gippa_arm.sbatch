#!/bin/bash
#SBATCH --job-name=gippa-ncm-run
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=64000M
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%j.out
# Phase 1b steps 3b/5: gippa convert + preflight --strict + gippa run for arm N1 or N2.
#   N1: the DIA-NN report, --search-db same, on NCM3722 (its Bakta JSON)
#   N2: the same report, --search-db other, on MG1655 U00096.3 (recorded Bakta JSON), --genome-relation related
# gippa run is scripts/slurm/gippa_run.sbatch's command with the venv entry point (no `uv run`) and
# the farmshare_local profile: every Nextflow task runs inside this 16-CPU allocation, so the arm
# never holds more than 16 CPUs (the campaign caps this agent at 48; the farmshare profile can
# hold 8 x 14). --threads 14 --memory 54.GB, all modules, default seed.
#   sbatch --exclude=... --time=36:00:00 -J gippa-ncm-run-n1 validation/campaign/ncm3722/gippa_arm.sbatch N1
# Rerun with other code (ADR-0078 fix): PIN=<pinned worktree> SUF=_fix writes in_<arm>_fix, preflight_<arm>_fix and
# out_<arm>_fix, never touching the earlier outputs.
set -euo pipefail
ARM="${1:?N1|N2}"
export GIPPA_ROOT=/home/classes/bios/270/khoa/giap/local
PIN="${PIN:-$GIPPA_ROOT/worktrees/campaign-pin}"
export GIPPA_REPO="$PIN"
source "$PIN/scripts/farmshare_env.sh"
export PYTHONPATH="$PIN/src"
export PATH="$UV_PROJECT_ENVIRONMENT/bin:$PATH"   # tasks call `gippa _step` (gippa_run.sbatch gets this from uv run)
export XDG_CACHE_HOME="$GIPPA_ROOT/.cache"
export TMPDIR="$(mktemp -d /tmp/gippa.XXXXXX)"; trap 'rm -rf "$TMPDIR"' EXIT
GIPPA="$UV_PROJECT_ENVIRONMENT/bin/gippa"
D="${D:-$GIPPA_ROOT/campaign/ncm3722}"
REPORT="$D/diann_a/report.parquet"
case "$ARM" in
  N1) SDB=same; G="$D/genome/ncm3722.fna.gz"; J="$D/bakta/out/ncm3722.json"; REL=() ;;
  N2) SDB=other; G="$PIN/tests/data/genomes/k12_mg1655.fna.gz"
      J="$GIPPA_ROOT/recordings/bakta_full_circ/k12_mg1655/bakta/k12_mg1655.json"; REL=(--genome-relation related) ;;
  *) echo "arm N1|N2" >&2; exit 2 ;;
esac
SUF="${SUF:-}"
a="${ARM,,}"; IN="$D/in_$a$SUF"; OUTD="$D/out_$a$SUF"; PF="$D/preflight_$a$SUF"
[ -z "$SUF" ] || [ ! -e "$OUTD" ] || { echo "$OUTD exists; refusing to overwrite" >&2; exit 3; }
mkdir -p "$IN"
cd "$PIN"
t0=$(date +%s)
"$GIPPA" convert --from diann "$REPORT" --out "$IN/peptides.tsv.gz" --abundance-out "$IN/abundance.tsv.gz" \
  --search-db "$SDB" 2>&1 | tee "$IN/convert.log"
t1=$(date +%s)
rm -rf "$PF"
"$GIPPA" preflight --genome "$G" --annotation "$J" --peptides "$IN/peptides.tsv.gz" \
  --outdir "$PF" --profile farmshare_local,apptainer --strict 2>&1 | tee "$IN/preflight.log"
pf=${PIPESTATUS[0]}
echo "preflight_exit=$pf"
[ "$pf" -eq 0 ] || exit "$pf"
t2=$(date +%s)
set +e
"$GIPPA" run --profile farmshare_local,apptainer --genome "$G" --annotation "$J" \
  --peptides "$IN/peptides.tsv.gz" --db-dir "$GIPPA_DB_DIR" --threads 14 --memory 54.GB "${REL[@]}" \
  --outdir "$OUTD"
rc=$?
t3=$(date +%s)
echo "arm=$ARM convert_s=$((t1-t0)) preflight_s=$((t2-t1)) run_s=$((t3-t2)) run_exit=$rc node=$(hostname -s) job=${SLURM_JOB_ID:-none}" | tee "$IN/runtime.txt"
sha256sum "$G" "$J" "$REPORT" "$IN/peptides.tsv.gz" "$IN/abundance.tsv.gz" > "$IN/inputs.sha256"
exit $rc
