#!/bin/bash
#SBATCH --job-name=gippa-ncm-diann
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=64000M
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%j.out
# Phase 1b step 3a: gippa makedb route A on the NCM3722 Bakta JSON (contaminants appended as
# CONT_ records), then DIA-NN 2.3.1 in the two calls of validation/pilot/diann_A2.cfg. The flags
# are read from that cfg (its --fasta lines replaced by the makedb FASTA), so they cannot drift;
# --report-decoys is kept; --met-excision/--reanalyse/--rt-profiling are absent there.
#   sbatch --exclude=... --time=12:00:00 validation/campaign/ncm3722/diann_n1.sbatch
set -euo pipefail
export GIPPA_ROOT=/home/classes/bios/270/khoa/giap/local
PIN="${PIN:-$GIPPA_ROOT/worktrees/campaign-pin}"
source /home/classes/bios/270/khoa/giap/scripts/farmshare_env.sh
export PYTHONPATH="$PIN/src"
GIPPA="$UV_PROJECT_ENVIRONMENT/bin/gippa"
D="${D:-$GIPPA_ROOT/campaign/ncm3722}"
CFG="$PIN/validation/pilot/diann_A2.cfg"
SIF="$GIPPA_ROOT/images/diann_2.3.1.sif"
PINSHA="$(cut -d' ' -f1 "$PIN/validation/pilot/diann_image.sha256")"
[ "$(sha256sum "$SIF" | cut -d' ' -f1)" = "$PINSHA" ] || { echo "image off pin" >&2; exit 2; }
T="${SLURM_CPUS_PER_TASK:-16}"
OUT="$D/diann_a"; mkdir -p "$OUT"
JSON="$D/bakta/out/ncm3722.json"
t0=$(date +%s)
if [ ! -s "$D/db_a/strain_db.faa" ]; then
  "$GIPPA" makedb --annotation "$JSON" --outdir "$D/db_a"
fi
echo "makedb_s=$(( $(date +%s) - t0 ))" > "$OUT/runtime.txt"
FASTA="$D/db_a/strain_db.faa"
grep -c '^>' "$FASTA"; grep -c '^>CONT_' "$FASTA" || true
flags() {  # [section] flags of the pilot cfg without --fasta and its value
  awk -v want="$1" '{sub(/#.*/, "")} /\[.*\]/ {gsub(/[^a-z0-9]/, ""); sec=$0; next}
    sec == want {for (i = 1; i <= NF; i++) { if ($i == "--fasta") {i++; continue}; print $i }}' "$CFG"
}
mapfile -t S1 < <(flags step1); mapfile -t S2 < <(flags step2)
RAWS=(); for r in "$D"/raw/*.raw; do RAWS+=(--f "$r"); done
[ "${#RAWS[@]}" -eq 6 ] || { echo "expected 3 raw files" >&2; exit 2; }
dn() { apptainer exec -B /home/classes -B /scratch -B /farmshare "$SIF" /diann-2.3.1/diann-linux "$@"; }
dn --version > "$OUT/diann_version.txt" 2>&1 || true
LIB="$OUT/n1.predicted.speclib"
t1=$(date +%s)
if [ ! -s "$LIB" ]; then
  printf '%s\n' "${S1[@]}" --fasta "$FASTA" --threads "$T" --out-lib "$OUT/n1.parquet" > "$OUT/command_step1.txt"
  dn "${S1[@]}" --fasta "$FASTA" --threads "$T" --out-lib "$OUT/n1.parquet"
fi
echo "step1_s=$(( $(date +%s) - t1 ))" >> "$OUT/runtime.txt"
t2=$(date +%s)
rm -rf "$OUT/quant"; mkdir -p "$OUT/quant"
printf '%s\n' "${S2[@]}" --fasta "$FASTA" "${RAWS[@]}" --threads "$T" --lib "$LIB" --temp "$OUT/quant" \
  --out "$OUT/report.parquet" --out-lib "$OUT/n1_empirical.parquet" > "$OUT/command_step2.txt"
dn "${S2[@]}" --fasta "$FASTA" "${RAWS[@]}" --threads "$T" --lib "$LIB" --temp "$OUT/quant" \
  --out "$OUT/report.parquet" --out-lib "$OUT/n1_empirical.parquet"
echo "step2_s=$(( $(date +%s) - t2 )) node=$(hostname -s) job=${SLURM_JOB_ID:-none} threads=$T" >> "$OUT/runtime.txt"
[ -s "$OUT/report.parquet" ] || { echo "no report" >&2; exit 1; }
sha256sum "$OUT/report.parquet" "$FASTA" "$CFG" "$SIF" > "$OUT/sha256.txt"
cat "$OUT/runtime.txt"
