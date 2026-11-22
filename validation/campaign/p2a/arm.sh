#!/bin/bash
# Campaign Phase 2, group A: one arm of the standard recipe, inside an existing 16-CPU allocation
# (called by lane.sbatch; never on a login node).
#
#   arm.sh ACC TAG ROUTE            ROUTE = a | b  (standard recipe steps 4-7)
#   arm.sh ACC TAG x REF_TAG        cross-strain arm (PXD065022 T6): TAG's route-A DIA-NN report,
#                                   --search-db other, on REF_TAG's genome, --genome-relation related
#   arm.sh ACC TAG te7              TE7 rerun of route A (out_a_te7)
#
# Layout: $GIPPA_ROOT/campaign/ACC/TAG/{runs.txt, diann_R, in_R, out_R}; databases in
# ACC/db/TAG/db_R (prep.sbatch). DIA-NN: diann_search.sh (pilot flags); convert + preflight
# --strict: convert_preflight.sh (twice, byte-identical, exit 0 required); gippa run: all modules,
# default seed, --threads 14 --memory 54.GB, profile farmshare_local,apptainer (tasks stay inside
# this allocation). A finished step (report.parquet + provenance exit 0; convert_provenance.json
# all ok; run_manifest.json exit_code 0) is skipped on a rerun. Exit 0 ok, else the failing step's.
set -uo pipefail
export GIPPA_ROOT=/home/classes/bios/270/khoa/giap/local
PIN="${PIN:-$GIPPA_ROOT/worktrees/rel111-pin}"
source /home/classes/bios/270/khoa/giap/scripts/farmshare_env.sh
export CAMPAIGN_REPO="${CAMPAIGN_REPO:-$GIPPA_ROOT/worktrees/campaign-p2a}"
HERE="$CAMPAIGN_REPO/validation/campaign/p2a"
ACC="$1" TAG="$2" R="$3" REF="${4:-}"
C="$GIPPA_ROOT/campaign/$ACC"
D="$C/$TAG"; mkdir -p "$D"
genome() { echo "$C/genome/$1/$1.fna.gz"; }
json() { echo "$C/bakta/$1/out/$1.json"; }
PY="$UV_PROJECT_ENVIRONMENT/bin/python"
log() { echo "[$(date +%F' '%T)] $ACC/$TAG/$R: $*"; }

runs_list() {  # runs.txt: the selected .d folders or .raw files of this TAG
  local pat
  case "$TAG" in
    tpa_dal1) pat=Dalas ;; tpa_grady) pat=Grady ;; tpa_haitib) pat=HAITI ;; tpa_madras) pat=Madras ;;
    tpa_mexicoa) pat=Mexiko ;; tpa_phil1) pat=PHIL ;; tpa_ss14) pat=SS14 ;; *) pat="" ;;
  esac
  if [ ! -d "$C/raw_d" ]; then ls "$C"/raw/*.raw
  elif [ -n "$pat" ]; then ls -d "$C"/raw_d/*.d | grep -- "_DIA_${pat}_"
  else ls -d "$C"/raw_d/*.d; fi
}

diann() {  # route
  local r="$1" out="$D/diann_$1" fas
  if [ -s "$out/report.parquet" ] && grep -q '"exit_code": 0' "$out/provenance_report.json" 2>/dev/null; then
    log "DIA-NN $r done, skipped"; return 0; fi
  runs_list > "$D/runs.txt"
  log "DIA-NN $r on $(wc -l < "$D/runs.txt") runs"
  fas="$C/db/$TAG/db_$r/strain_db.faa"
  [ "$ACC" = PXD065022 ] && fas="$fas:$C/host/rabbit_host.fasta"
  bash "$HERE/diann_search.sh" "$out" "$fas" "$D/runs.txt"
}

convert() {  # route report search_db genome json in_dir
  local in="$6"
  if [ -s "$in/convert_provenance.json" ] && "$PY" -I -c 'import json,sys; d=json.load(open(sys.argv[1])); sys.exit(0 if d["convert_exit"]==0 and d["deterministic"] and d["preflight_strict_exit"]==0 else 1)' "$in/convert_provenance.json"; then
    log "convert/preflight done, skipped"; return 0; fi
  log "convert + preflight --strict ($3)"
  bash "$HERE/convert_preflight.sh" "$PIN" "$2" "$in" "$4" "$5" "$3"
}

run() {  # in_dir out_dir genome json [extra args]
  local in="$1" out="$2" g="$3" j="$4"; shift 4
  if [ -s "$out/run_manifest.json" ] && "$PY" -I -c 'import json,sys; sys.exit(0 if json.load(open(sys.argv[1])).get("exit_code")==0 else 1)' "$out/run_manifest.json"; then
    log "gippa run done, skipped"; return 0; fi
  [ ! -e "$out" ] || mv "$out" "$out.failed.$(date +%s)"
  log "gippa run -> $out"
  (
    export GIPPA_REPO="$PIN" PYTHONPATH="$PIN/src"
    export PATH="$UV_PROJECT_ENVIRONMENT/bin:$PATH"   # Nextflow tasks call `gippa _step`
    export XDG_CACHE_HOME="$GIPPA_ROOT/.cache"
    TMPDIR="$(mktemp -d /tmp/gippa.XXXXXX)"; export TMPDIR
    cd "$PIN" || exit 2
    s=$(date +%s)
    "$UV_PROJECT_ENVIRONMENT/bin/gippa" run --profile farmshare_local,apptainer --genome "$g" --annotation "$j" \
      --peptides "$in/peptides.tsv.gz" --db-dir "$GIPPA_DB_DIR" --threads 14 --memory 54.GB "$@" --outdir "$out"
    rc=$?
    echo "run_s=$(( $(date +%s) - s )) run_exit=$rc node=$(hostname -s) job=${SLURM_JOB_ID:-none}" > "$in/runtime_$(basename "$out").txt"
    rm -rf "$TMPDIR"
    exit $rc
  )
}

case "$R" in
  a|b)
    diann "$R" || exit $?
    convert "$R" "$D/diann_$R/report.parquet" same "$(genome "$TAG")" "$(json "$TAG")" "$D/in_$R" || exit $?
    run "$D/in_$R" "$D/out_$R" "$(genome "$TAG")" "$(json "$TAG")" || exit $?
    ;;
  te7)
    run "$D/in_a" "$D/out_a_te7" "$(genome "$TAG")" "$(json "$TAG")" || exit $?
    ;;
  x)
    [ -n "$REF" ] || { echo "x needs REF_TAG" >&2; exit 2; }
    convert x "$D/diann_a/report.parquet" other "$(genome "$REF")" "$(json "$REF")" "$D/in_x_$REF" || exit $?
    run "$D/in_x_$REF" "$D/out_x_$REF" "$(genome "$REF")" "$(json "$REF")" --genome-relation related || exit $?
    ;;
  *) echo "route a|b|te7|x" >&2; exit 2 ;;
esac
log "ok"
