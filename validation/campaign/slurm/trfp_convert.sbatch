#!/bin/bash
#SBATCH --job-name=gippa-campaign-trfp
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=16000M
#SBATCH --time=04:00:00
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%j.out
# Thermo .raw -> indexed mzML with ThermoRawFileParser 1.4.5 (ADR-0077 decision 2), for Cascadia
# (arm E4). Vendor peak picking (the TRFP default: centroided MS1 and MS2), no gzip, so Cascadia
# reads <run>.mzML directly. Usage: sbatch trfp_convert.sbatch OUT_DIR RAW [RAW ...]
# Writes OUT_DIR/<run>.mzML, <run>.trfp.log and appends sha256 to OUT_DIR/SHA256SUMS.
set -uo pipefail
source "${GIPPA_REPO:-/home/classes/bios/270/khoa/giap}/scripts/farmshare_env.sh"
OUT="$1"; shift
case "$OUT" in "$GIPPA_ROOT"/campaign/*) ;; *) echo "OUT must be under \$GIPPA_ROOT/campaign" >&2; exit 2 ;; esac
ENV="$MAMBA_ROOT_PREFIX/envs/gippa-thermorawfileparser"
[ -d "$ENV" ] || { echo "missing env $ENV (build_trfp_env.sbatch)" >&2; exit 2; }
mkdir -p "$OUT"
micromamba run -p "$ENV" ThermoRawFileParser.sh --version
rc_all=0
for raw in "$@"; do
  run="$(basename "$raw" .raw)"
  s=$(date +%s)
  micromamba run -p "$ENV" ThermoRawFileParser.sh -i "$raw" -o "$OUT" -f 2 -L 1,2 >"$OUT/$run.trfp.log" 2>&1
  rc=$?
  echo "$run: exit $rc wall_s $(( $(date +%s) - s ))"; tail -3 "$OUT/$run.trfp.log"
  if [ "$rc" -ne 0 ] || [ ! -s "$OUT/$run.mzML" ]; then rc_all=1; continue; fi
  sha="$(sha256sum "$OUT/$run.mzML" | cut -d' ' -f1)"
  echo "$run.mzML $(stat -c %s "$OUT/$run.mzML") bytes sha256 $sha"
  { grep -v "  $run.mzML\$" "$OUT/SHA256SUMS" 2>/dev/null; printf '%s  %s\n' "$sha" "$run.mzML"; } | sort -k2 >"$OUT/SHA256SUMS.$$" && mv "$OUT/SHA256SUMS.$$" "$OUT/SHA256SUMS"
done
exit $rc_all
