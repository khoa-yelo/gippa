#!/bin/bash
#SBATCH --job-name=gippa-campaign-rawcheck
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=32000M
#SBATCH --time=03:00:00
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%j.out
# Phase 0 check "Thermo .raw on Linux" (campaign handover plan; ADR-0077 decision 2).
#
#   sbatch --exclude=barley-01,barley-02,barley-03,barley-04,rye-01,rye-02 \
#     validation/campaign/slurm/raw_check.sbatch FILE_NAME
#
# 1. Downloads FILE_NAME (one row of validation/campaign/PXD058422/files.tsv) over https to
#    $GIPPA_ROOT/campaign/ecoli_k12/raw/, verifies the deposit checksum (one retry), records sha256.
# 2. DIA-NN 2.3.1 (the verified read-only image of diann_image.sbatch) converts the .raw to .dia
#    (--convert): the direct test of whether DIA-NN reads Thermo .raw on Linux.
# 3. A short library-free search of that run against the K-12 Bakta proteins (same flags as
#    validation/pilot/diann_A2.cfg) shows the spectra are usable. Only the log is inspected; no
#    report row goes into git.
# Exit 0 both steps ok; 1 download/verify failed; 2 precondition; 10 + DIA-NN exit of a step.
set -uo pipefail
source "${GIPPA_REPO:-/home/classes/bios/270/khoa/giap}/scripts/farmshare_env.sh"
module load apptainer/1.4.4 >/dev/null 2>&1 || true

FILE="${1:?usage: raw_check.sbatch FILE_NAME}"
REPO="${CAMPAIGN_REPO:-$GIPPA_REPO/local/worktrees/campaign-phase0}"
MANIFEST="$REPO/validation/campaign/PXD058422/files.tsv"
SIF="${DIANN_SIF:-$GIPPA_ROOT/images/diann_2.3.1.sif}"
PIN="$(cut -d' ' -f1 "$GIPPA_REPO/validation/pilot/diann_image.sha256")"
DIANN=/diann-2.3.1/diann-linux
RAW_DIR="$GIPPA_ROOT/campaign/ecoli_k12/raw"
OUT="$GIPPA_ROOT/campaign/ecoli_k12/raw_check"
FAA="${CHECK_FASTA:-$GIPPA_ROOT/recordings/bakta_full_circ/k12_mg1655/bakta/k12_mg1655.faa}"
CONTAM="$GIPPA_DB_DIR/contaminants/contaminants.fasta"
T="${SLURM_CPUS_PER_TASK:-8}"
mkdir -p "$RAW_DIR" "$OUT"

case "$FILE" in */* | *..* | *[!A-Za-z0-9._-]*) echo "unsafe file name: $FILE" >&2; exit 2 ;; esac
row="$(awk -F'\t' -v f="$FILE" '$1 == f' "$MANIFEST")"
[ -n "$row" ] || { echo "$FILE is not in $MANIFEST" >&2; exit 2; }
URL="$(cut -f7 <<<"$row")" BYTES="$(cut -f4 <<<"$row")" CK="$(cut -f5 <<<"$row")"
case "$URL" in https://ftp.pride.ebi.ac.uk/*) ;; *) echo "not an ADR-0077 https URL: $URL" >&2; exit 2 ;; esac
[ "$(sha256sum "$SIF" | cut -d' ' -f1)" = "$PIN" ] || { echo "$SIF does not match the pin" >&2; exit 2; }
for f in "$FAA" "$CONTAM"; do [ -f "$f" ] || { echo "missing $f" >&2; exit 2; }; done

verify() {
  [ -f "$RAW_DIR/$FILE" ] && [ "$(stat -c %s "$RAW_DIR/$FILE")" = "$BYTES" ] || return 1
  case "$CK" in
    sha1:*) [ "$(sha1sum "$RAW_DIR/$FILE" | cut -d' ' -f1)" = "${CK#sha1:}" ] ;;
    md5:*) [ "$(md5sum "$RAW_DIR/$FILE" | cut -d' ' -f1)" = "${CK#md5:}" ] ;;
    *) true ;;
  esac
}
echo "== 1. download $URL ($BYTES bytes, $CK)"
if ! verify; then
  ok=0
  for attempt in 1 2; do
    curl --fail --location --silent --show-error --proto '=https' --proto-redir '=https' \
      --continue-at - --max-time 14400 --output "$RAW_DIR/$FILE" "$URL" && verify && { ok=1; break; }
    echo "attempt $attempt failed (curl or checksum)" >&2
    # keep a short partial file for the resume (--continue-at -); drop it only when it is
    # complete or too long and still fails the checksum
    if [ -f "$RAW_DIR/$FILE" ] && [ "$(stat -c %s "$RAW_DIR/$FILE")" -ge "$BYTES" ]; then
      verify || rm -f "$RAW_DIR/$FILE"
    fi
  done
  [ "$ok" -eq 1 ] || { echo "download failed after one retry: record the block (R6)" >&2; exit 1; }
fi
echo "verified: $CK"
echo "sha256: $(sha256sum "$RAW_DIR/$FILE" | cut -d' ' -f1)"

run() {  # label, DIA-NN args...
  local label="$1"; shift
  echo "== $label: apptainer exec $SIF $DIANN $*"
  local s; s=$(date +%s)
  apptainer exec --bind "$GIPPA_ROOT" "$SIF" "$DIANN" "$@" >"$OUT/$label.log" 2>&1
  local rc=$?
  echo "exit=$rc wall_s=$(( $(date +%s) - s ))"
  head -40 "$OUT/$label.log"; echo "..."; tail -25 "$OUT/$label.log"
  return $rc
}

echo "== 2. DIA-NN --convert"
rm -f "$RAW_DIR/$FILE.dia"
run convert --f "$RAW_DIR/$FILE" --convert --threads "$T" --verbose 1
RC1=$?
ls -la "$RAW_DIR" | grep -F "${FILE%.raw}" || true

echo "== 3. library-free search (two calls, diann_A2.cfg flags)"
COMMON=(--cut 'K*,R*' --missed-cleavages 1 --min-pep-len 7 --max-pep-len 30 --unimod4
  --var-mods 1 --var-mod UniMod:35,15.994915,M --min-pr-charge 1 --max-pr-charge 4
  --min-pr-mz 300 --max-pr-mz 1800 --min-fr-mz 200 --max-fr-mz 1800 --peptidoforms)
run predict --fasta-search --gen-spec-lib --predictor --fasta "$FAA" --fasta "$CONTAM" \
  "${COMMON[@]}" --out-lib "$OUT/k12.speclib" --threads "$T" --verbose 1
RC2=$?
LIB="$(ls "$OUT"/k12*.predicted.speclib 2>/dev/null | head -1)"
RC3=99
if [ "$RC2" -eq 0 ] && [ -n "$LIB" ]; then
  run search --f "$RAW_DIR/$FILE" --lib "$LIB" --fasta "$FAA" --fasta "$CONTAM" \
    --qvalue 0.01 --report-decoys --matrices --peptidoforms --out "$OUT/report.parquet" \
    --threads "$T" --verbose 1
  RC3=$?
fi
echo "== summary convert=$RC1 predict=$RC2 search=$RC3"
[ "$RC1" -eq 0 ] || exit $((10 + RC1 % 100))
[ "$RC3" -eq 0 ] || exit $((10 + RC3 % 100))
exit 0
