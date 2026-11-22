#!/bin/bash







# DIA-NN 2.3.1 two-call search for one campaign arm (generalises scripts/slurm/diann_pilot.sbatch,
# which is hard-wired to PXD078459 and Bruker .d runs).
#
#   sbatch -J gippa-campaign-diann-e1a --time=08:00:00 --exclude=... \
#     validation/campaign/slurm/diann_search.sbatch OUT_DIR FASTAS RUN_LIST [REPORT_STEM]
#
# OUT_DIR   arm directory under $GIPPA_ROOT/campaign (library, report, provenance.json)
# FASTAS    colon-separated absolute FASTA paths (first: the gippa makedb strain_db.faa)
# RUN_LIST  file with one absolute raw path per line (Thermo .raw or Bruker .d); "-" = none
#           (step 1 only: predict the library and stop)
# REPORT_STEM  report name stem (default: report); a probe search uses another stem so that
#           it reuses the arm's library without touching the full report.
#
# Settings: validation/campaign/diann_campaign.cfg (= the flags of validation/pilot/diann_A2.cfg).
# Step 1 (library from the FASTA alone) is skipped when step1.done records the sha256 of the
# current library, cfg, FASTA files and image. Step 2 always starts fresh, with its own --temp so
# that concurrent arms never share .quant files (pilot job 1879914). DIA-NN runs only from the
# verified read-only image ($GIPPA_ROOT/images/diann_2.3.1.sif, pin validation/pilot/diann_image.sha256).
# provenance.json is written on every exit. Exit 0 ok; 1 DIA-NN failed; 2 precondition.
set -uo pipefail
source "${GIPPA_REPO:-/home/classes/bios/270/khoa/giap}/scripts/farmshare_env.sh"
module load apptainer/1.4.4 >/dev/null 2>&1 || true

[ $# -ge 3 ] || { echo "usage: diann_search.sbatch OUT_DIR FASTAS RUN_LIST [REPORT_STEM]" >&2; exit 2; }
OUT="$1" FASTAS="$2" RUN_LIST="$3" STEM="${4:-report}"
HERE="${CAMPAIGN_REPO:-$GIPPA_REPO/local/worktrees/campaign-p2a}/validation/campaign/p2a"
CFG="${CFG:-$HERE/diann_campaign.cfg}"
SIF="${DIANN_SIF:-$GIPPA_ROOT/images/diann_2.3.1.sif}"
PIN="$(cut -d' ' -f1 "$GIPPA_REPO/validation/pilot/diann_image.sha256")"
DIANN=/diann-2.3.1/diann-linux
T="${SLURM_CPUS_PER_TASK:-16}"
case "$OUT" in "$GIPPA_ROOT"/campaign/*) ;; *) echo "OUT_DIR must be under \$GIPPA_ROOT/campaign" >&2; exit 2 ;; esac
mkdir -p "$OUT"
LIB_STEM="$OUT/lib"
LIB="$LIB_STEM.predicted.speclib"
MARKER="$OUT/step1.done"
START=$(date +%s) S1=not_reached S1_RC="" S1_S="" S2=not_reached S2_RC="" S2_S="" NRUNS=0

sha() { sha256sum "$1" | cut -d' ' -f1; }
provenance() {
  local rc=$?
  trap - EXIT
  for f in ${FASTAS//:/ }; do printf '%s\t%s\n' "$f" "$(sha "$f" 2>/dev/null)"; done >"$OUT/.fastas.$STEM.tsv"
  P_RC=$rc P_START=$START P_S1=$S1 P_S1_RC=$S1_RC P_S1_S=$S1_S P_S2=$S2 P_S2_RC=$S2_RC P_S2_S=$S2_S \
  P_NRUNS=$NRUNS P_SIF=$SIF P_PIN=$PIN P_CFG=$CFG P_CFG_SHA="$(sha "$CFG" 2>/dev/null)" P_LIB=$LIB \
  P_COMMIT="$(git -C "$HERE" rev-parse HEAD 2>/dev/null)" P_RUNS=$RUN_LIST P_STEM=$STEM P_OUT=$OUT \
  python3 -I -c '
import json, os, time
e = os.environ.get
n = lambda k: (int(e(k)) if e(k) else None)
fas = dict(l.rstrip("\n").split("\t") for l in open(e("P_OUT") + "/.fastas." + e("P_STEM") + ".tsv"))
lib = e("P_LIB")
d = dict(job_id=e("SLURM_JOB_ID"), node=os.uname().nodename, exit_code=n("P_RC"),
         image=e("P_SIF"), image_sha256=e("P_PIN"), cfg=e("P_CFG"), cfg_sha256=e("P_CFG_SHA"),
         campaign_commit=e("P_COMMIT"), fasta_sha256=fas, run_list=e("P_RUNS"), n_runs=n("P_NRUNS"),
         report_stem=e("P_STEM"),
         step1=dict(status=e("P_S1"), exit=n("P_S1_RC"), wall_s=n("P_S1_S")),
         step2=dict(status=e("P_S2"), exit=n("P_S2_RC"), wall_s=n("P_S2_S")),
         wall_s=int(time.time()) - int(e("P_START")),
         library_bytes=(os.path.getsize(lib) if os.path.exists(lib) else None))
json.dump(d, open(e("P_OUT") + "/provenance_" + e("P_STEM") + ".json", "w"), indent=1)
' || echo "warning: provenance not written" >&2
  exit $rc
}
trap provenance EXIT
trap 'exit 143' TERM

[ "$(sha "$SIF")" = "$PIN" ] || { echo "$SIF does not hash to the pin" >&2; exit 2; }
[ $((0$(stat -c %a "$SIF") & 0222)) -eq 0 ] || { echo "$SIF is writable" >&2; exit 2; }
[ -f "$CFG" ] || { echo "missing $CFG" >&2; exit 2; }
FA=()
for f in ${FASTAS//:/ }; do
  case "$f" in /*) ;; *) echo "FASTA not absolute: $f" >&2; exit 2 ;; esac
  [ -s "$f" ] || { echo "missing FASTA $f" >&2; exit 2; }
  FA+=(--fasta "$f")
done
RAW=()
if [ "$RUN_LIST" != "-" ]; then
  while IFS= read -r r; do
    [ -n "$r" ] || continue
    [ -e "$r" ] || { echo "missing run $r" >&2; exit 2; }
    RAW+=(--f "$r"); NRUNS=$((NRUNS + 1))
  done <"$RUN_LIST"
  [ "$NRUNS" -gt 0 ] || { echo "no runs in $RUN_LIST" >&2; exit 2; }
fi

section() {  # flags of [step1]/[step2]
  set -f
  awk -v want="[$1]" '/^\[/ {s=$0; next} s==want {sub(/#.*/, ""); if (NF) print}' "$CFG" | tr -s ' \n' '\n'
  set +f
}
mapfile -t F1 < <(section step1)
mapfile -t F2 < <(section step2)
diann() { apptainer exec -B /home/classes -B /scratch -B /farmshare "$SIF" "$DIANN" "$@"; }
diann --version 2>&1 | grep -m1 . >"$OUT/diann_version.txt" || true
grep -q "DIA-NN 2.3.1 Academia" "$OUT/diann_version.txt" || { echo "image banner is not DIA-NN 2.3.1 Academia" >&2; exit 2; }

INPUTS="$( { echo "cfg $(sha "$CFG")"; echo "image $PIN"; for f in ${FASTAS//:/ }; do echo "fasta $(sha "$f") $f"; done; } )"
if [ -f "$MARKER" ] && [ -f "$LIB" ] && [ "$(head -1 "$MARKER")" = "$(sha "$LIB")" ] \
   && [ "$(sed 1,2d "$MARKER")" = "$INPUTS" ]; then
  echo "step 1: library matches step1.done; skipped ($(sed -n 2p "$MARKER"))"
  S1=skipped
else
  echo "step 1: predicting the library"
  rm -f "$MARKER" "$LIB" "$LIB_STEM.parquet"
  S1=ran; s=$(date +%s)
  printf '%s\n' apptainer exec "$SIF" "$DIANN" "${F1[@]}" "${FA[@]}" --threads "$T" --out-lib "$LIB_STEM.parquet" >"$OUT/command_step1.txt"
  diann "${F1[@]}" "${FA[@]}" --threads "$T" --out-lib "$LIB_STEM.parquet" >"$OUT/step1.log" 2>&1
  S1_RC=$?; S1_S=$(( $(date +%s) - s ))
  tail -5 "$OUT/step1.log"
  [ "$S1_RC" -eq 0 ] && [ -s "$LIB" ] || { echo "step 1 failed ($S1_RC)" >&2; exit 1; }
  { sha "$LIB"; echo "job_id=${SLURM_JOB_ID:-} wall_s=$S1_S"; echo "$INPUTS"; } >"$MARKER"
fi
echo "library: $(stat -c %s "$LIB") bytes"
grep -E -i "precursors|proteins|genes" "$OUT/step1.log" 2>/dev/null | tail -5 || true
[ "$RUN_LIST" != "-" ] || { echo "step 1 only"; exit 0; }

echo "step 2: searching $NRUNS run(s)"
Q="$OUT/quant_$STEM"; rm -rf "$Q"; mkdir -p "$Q"
rm -f "$OUT/$STEM".* "$OUT/${STEM}_empirical".*
S2=ran; s=$(date +%s)
printf '%s\n' apptainer exec "$SIF" "$DIANN" "${F2[@]}" "${FA[@]}" "${RAW[@]}" --threads "$T" --lib "$LIB" --temp "$Q" --out "$OUT/$STEM.parquet" --out-lib "$OUT/${STEM}_empirical.parquet" >"$OUT/command_step2_$STEM.txt"
diann "${F2[@]}" "${FA[@]}" "${RAW[@]}" --threads "$T" --lib "$LIB" --temp "$Q" \
  --out "$OUT/$STEM.parquet" --out-lib "$OUT/${STEM}_empirical.parquet" >"$OUT/step2_$STEM.log" 2>&1
S2_RC=$?; S2_S=$(( $(date +%s) - s ))
tail -8 "$OUT/step2_$STEM.log"
[ "$S2_RC" -eq 0 ] && [ -s "$OUT/$STEM.parquet" ] || { echo "step 2 failed ($S2_RC)" >&2; exit 1; }
sha256sum "$OUT/$STEM.parquet" >"$OUT/$STEM.parquet.sha256"
[ -z "${SLURM_JOB_ID:-}" ] || seff "$SLURM_JOB_ID" 2>/dev/null | tail -4 || true
exit 0
