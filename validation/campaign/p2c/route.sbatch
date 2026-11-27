#!/bin/bash
#SBATCH --job-name=gippa-camp-p2c-route
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=64000M
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%j.out
# Campaign Phase 2 group C: the standard recipe (handover plan) for one route of one dataset, in
# one 16-CPU allocation, steps in order:
#   0. raw check: sha256sum -c raw/SHA256SUMS; Bruker .d.tar archives unpacked once (flock) to raw_d/
#   3. gippa makedb (route a: called genes; route b: + --six-frame-orfs --min-orf-aa 20)
#   4. DIA-NN 2.3.1, two calls (p2c/diann_search.sbatch run as a script; flags = pilot diann_A2.cfg)
#   5-6. gippa convert --search-db same (twice, byte-identical) + preflight --strict (p2c/convert_preflight.sbatch)
#   7. gippa run, all modules, default seed, --threads 14 --memory 54.GB, profile
#      farmshare_local,apptainer: every Nextflow task runs inside this allocation (the group's cap
#      is 32 CPUs, so the farmshare profile's per-task jobs are not used), as Phase 1b gippa_arm.sbatch.
#
#   sbatch -J gippa-camp-p2c-<ds>-<route> --time=24:00:00 --exclude=... route.sbatch DS LABEL a|b
#
# DS = accession (data under $GIPPA_ROOT/campaign/DS), LABEL = genome/Bakta label
# ($D/genome/LABEL/LABEL.fna.gz, $D/bakta/LABEL/out/LABEL.json). SUF= suffixes in_/out_ dirs
# and refuses to overwrite an existing out dir. Code: PIN (default rel111-pin, GIPPA 1.1.1).
set -euo pipefail
DS="${1:?DS}" LABEL="${2:?LABEL}" R="${3:?a|b}"
case "$R" in a|b) ;; *) echo "route a|b" >&2; exit 2 ;; esac
export GIPPA_ROOT=/home/classes/bios/270/khoa/giap/local
PIN="${PIN:-$GIPPA_ROOT/worktrees/rel111-pin}"
REPO="${REPO:-$GIPPA_ROOT/worktrees/campaign-p2c}"
HERE="$REPO/validation/campaign/p2c"
export GIPPA_REPO="$PIN"
source "$PIN/scripts/farmshare_env.sh"
module load apptainer/1.4.4 >/dev/null 2>&1 || true
export PYTHONPATH="$PIN/src"
export PATH="$UV_PROJECT_ENVIRONMENT/bin:$PATH"   # Nextflow tasks call `gippa _step`
export XDG_CACHE_HOME="$GIPPA_ROOT/.cache"
export TMPDIR="$(mktemp -d /tmp/gippa.XXXXXX)"; trap 'rm -rf "$TMPDIR"' EXIT
GIPPA="$UV_PROJECT_ENVIRONMENT/bin/gippa"
D="$GIPPA_ROOT/campaign/$DS"
G="$D/genome/$LABEL/$LABEL.fna.gz"
J="$D/bakta/$LABEL/out/$LABEL.json"
SUF="${SUF:-}"
DB="$D/db_$R"; DN="$D/diann_$R"; IN="$D/in_$R$SUF"; OUT="$D/out_$R$SUF"
for f in "$G" "$J"; do [ -s "$f" ] || { echo "missing $f" >&2; exit 2; }; done
[ ! -e "$OUT" ] || { echo "$OUT exists; refusing to overwrite" >&2; exit 3; }
echo "ds=$DS label=$LABEL route=$R pin=$(git -C "$PIN" rev-parse HEAD) repo=$(git -C "$REPO" rev-parse HEAD) job=${SLURM_JOB_ID:-none} node=$(hostname -s)"
T0=$(date +%s)
RT="$D/runtime_$R$SUF.txt"; : > "$RT"
tick() { echo "$1_s=$(( $(date +%s) - $2 ))" | tee -a "$RT"; }

# 0. raw files
t=$(date +%s)
(
  flock 9
  if [ ! -f "$D/raw/.sha256_ok" ]; then
    (cd "$D/raw" && sha256sum -c SHA256SUMS) && touch "$D/raw/.sha256_ok"
  fi
  if ls "$D"/raw/*.d.tar >/dev/null 2>&1 && [ ! -f "$D/raw_d/.done" ]; then
    rm -rf "$D/raw_d"; mkdir -p "$D/raw_d"
    for a in "$D"/raw/*.d.tar; do tar -xf "$a" -C "$D/raw_d"; done
    touch "$D/raw_d/.done"
  fi
  if [ ! -s "$D/runs.txt" ]; then
    if [ -f "$D/raw_d/.done" ]; then find "$D/raw_d" -maxdepth 2 -type d -name '*.d' | sort > "$D/runs.txt"
    else ls -1 "$D"/raw/*.raw | sort > "$D/runs.txt"; fi
  fi
) 9>"$D/.raw.lock"
cat "$D/runs.txt"
tick raw "$t"

# 3. makedb
t=$(date +%s)
if [ ! -s "$DB/strain_db.faa" ]; then
  rm -rf "$DB"
  if [ "$R" = a ]; then "$GIPPA" makedb --annotation "$J" --genome "$G" --outdir "$DB"
  else "$GIPPA" makedb --annotation "$J" --genome "$G" --six-frame-orfs --min-orf-aa 20 --outdir "$DB"; fi
fi
FASTA="$DB/strain_db.faa"
echo "makedb records=$(grep -c '^>' "$FASTA") cont=$(grep -c '^>CONT_' "$FASTA" || true) bytes=$(stat -c %s "$FASTA")" | tee -a "$RT"
tick makedb "$t"

# 4. DIA-NN
t=$(date +%s)
if [ ! -s "$DN/report.parquet" ]; then
  GIPPA_REPO="$PIN" CAMPAIGN_REPO="$REPO" CFG="$HERE/diann_campaign.cfg" \
    bash "$HERE/diann_search.sbatch" "$DN" "$FASTA" "$D/runs.txt"
fi
tick diann "$t"

# 5-6. convert + preflight --strict
t=$(date +%s)
bash "$HERE/convert_preflight.sbatch" "$PIN" "$DN/report.parquet" "$IN" "$G" "$J" same
tick convert_preflight "$t"

# 7. gippa run
t=$(date +%s)
cd "$PIN"
set +e
"$GIPPA" run --profile farmshare_local,apptainer --genome "$G" --annotation "$J" \
  --peptides "$IN/peptides.tsv.gz" --db-dir "$GIPPA_DB_DIR" --threads 14 --memory 54.GB --outdir "$OUT"
rc=$?
set -e
tick gippa_run "$t"
echo "run_exit=$rc total_s=$(( $(date +%s) - T0 )) node=$(hostname -s) job=${SLURM_JOB_ID:-none}" | tee -a "$RT"
exit $rc
