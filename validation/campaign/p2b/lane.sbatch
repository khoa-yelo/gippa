#!/bin/bash
#SBATCH --job-name=gippa-camp-p2b-lane
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=64000M
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%j.out
# Campaign Phase 2 group B: the standard recipe for one dataset, every step inside this one
# 16-CPU allocation (the group's share is 32 CPUs, so two lanes at most), resumable: a step
# whose output is complete is skipped, so a resubmission continues where a failed job stopped.
#
#   sbatch -J gippa-camp-p2b-<acc> --time=... --exclude=... lane.sbatch ACC [STEP ...]
#
# Steps (default: all, in order): makedb lib diann_a diann_b conv_a conv_b run_a run_b post te7
#   makedb   gippa makedb route A and route B (--six-frame-orfs --min-orf-aa 20)
#   lib      DIA-NN step 1 only (in-silico library) for A and B: library size before any search
#   diann_X  DIA-NN step 2 on the selected runs (validation/campaign/p2b/diann_search.sbatch)
#   conv_X   gippa convert (twice, byte-identical) + preflight --strict (convert_preflight.sbatch)
#   run_X    gippa run, all modules, default seed, --threads 14 --memory 54.GB,
#            profile farmshare_local,apptainer (every Nextflow task inside this allocation)
#   post     gippa compare out_a out_b; gippa report --redact on both
#   te7      route A rerun into out_a_te7 and diff -r of tables/ (TE7)
# Outputs: $GIPPA_ROOT/campaign/<ACC>/{db_a,db_b,diann_a,diann_b,in_a,in_b,out_a,out_b,cmp_ab}
# and lane_steps.tsv (step, exit, wall seconds, job). A complete out_X is never overwritten; an
# incomplete one is moved aside to out_X.failed.<epoch>.
set -uo pipefail
ACC="${1:?ACC}"; shift
STEPS=("$@"); [ ${#STEPS[@]} -gt 0 ] || STEPS=(makedb lib diann_a diann_b conv_a conv_b run_a run_b post)
export GIPPA_ROOT=/home/classes/bios/270/khoa/giap/local
PIN="${PIN:-$GIPPA_ROOT/worktrees/rel111-pin}"
HERE=/home/classes/bios/270/khoa/giap/local/worktrees/campaign-p2b/validation/campaign/p2b
export GIPPA_REPO="$PIN" CAMPAIGN_REPO=/home/classes/bios/270/khoa/giap/local/worktrees/campaign-p2b
source "$PIN/scripts/farmshare_env.sh"
export GIPPA_ROOT=/home/classes/bios/270/khoa/giap/local
export PYTHONPATH="$PIN/src"
export PATH="$UV_PROJECT_ENVIRONMENT/bin:$PATH"   # Nextflow tasks call `gippa _step`
export XDG_CACHE_HOME="$GIPPA_ROOT/.cache"
export TMPDIR="$(mktemp -d /tmp/gippa.XXXXXX)"; trap 'rm -rf "$TMPDIR"' EXIT
GIPPA="$UV_PROJECT_ENVIRONMENT/bin/gippa"
case "$ACC" in
  PXD049476) L=vnat_atcc14048 ;;
  PXD065346) L=cdif_630derm ;;
  PXD072227) L=syn_cc9311 ;;
  PXD060688) L=mxan_dk1622 ;;
  *) echo "unknown dataset $ACC" >&2; exit 2 ;;
esac
D="$GIPPA_ROOT/campaign/$ACC"
G="$D/genome/$L/$L.fna.gz"
J="$D/bakta/$L/out/$L.json"
EXTRA=""; [ "$ACC" = PXD065346 ] && EXTRA=":$D/spike/spike.faa"
echo "pin $(git -C "$PIN" rev-parse HEAD) dirty=$(git -C "$PIN" status --porcelain --untracked-files=no | wc -l) node=$(hostname -s) job=${SLURM_JOB_ID:-}"
ls "$D"/raw/*.raw > "$D/runs.txt"
log() { printf '%s\t%s\t%s\t%s\n' "$1" "$2" "$3" "${SLURM_JOB_ID:-}" >> "$D/lane_steps.tsv"; }
step() {  # name command...
  local n="$1"; shift; local s; s=$(date +%s)
  echo "=== $n $(date -Is)"; "$@"; local rc=$?
  log "$n" "$rc" $(( $(date +%s) - s )); echo "=== $n exit $rc"
  [ "$rc" -eq 0 ] || exit "$rc"
}
done_run() { [ -f "$1/run_manifest.json" ] && python3 -I -c 'import json,sys; m=json.load(open(sys.argv[1])); sys.exit(0 if m.get("exit_code")==0 else 1)' "$1/run_manifest.json"; }
makedb() {
  [ -s "$D/db_a/strain_db.faa" ] || "$GIPPA" makedb --annotation "$J" --genome "$G" --outdir "$D/db_a" || return 1
  [ -s "$D/db_b/strain_db.faa" ] || "$GIPPA" makedb --annotation "$J" --genome "$G" --six-frame-orfs --min-orf-aa 20 --outdir "$D/db_b" || return 1
  for r in a b; do echo "db_$r records=$(grep -c '^>' "$D/db_$r/strain_db.faa") cont=$(grep -c '^>CONT_' "$D/db_$r/strain_db.faa") bytes=$(stat -c %s "$D/db_$r/strain_db.faa")"; done
}
diann() {  # route runlist
  local r="$1" rl="$2"
  if [ "$rl" != "-" ] && [ -s "$D/diann_$r/report.parquet" ] && [ -s "$D/diann_$r/report.parquet.sha256" ]; then echo "diann_$r done"; return 0; fi
  bash "$HERE/diann_search.sbatch" "$D/diann_$r" "$D/db_$r/strain_db.faa$EXTRA" "$rl"
}
conv() {
  local r="$1"
  if [ -f "$D/in_$r/convert_provenance.json" ] && python3 -I -c 'import json,sys; d=json.load(open(sys.argv[1])); sys.exit(0 if d["preflight_strict_exit"]==0 and d["deterministic"] else 1)' "$D/in_$r/convert_provenance.json"; then echo "conv_$r done"; return 0; fi
  bash "$HERE/convert_preflight.sbatch" "$PIN" "$D/diann_$r/report.parquet" "$D/in_$r" "$G" "$J" same
}
run() {  # route outdir
  local r="$1" o="$2"
  if done_run "$o"; then echo "$o complete"; return 0; fi
  [ ! -e "$o" ] || mv "$o" "$o.failed.$(date +%s)"
  ( cd "$PIN" && "$GIPPA" run --profile farmshare_local,apptainer --genome "$G" --annotation "$J" \
      --peptides "$D/in_$r/peptides.tsv.gz" --db-dir "$GIPPA_DB_DIR" --threads 14 --memory 54.GB --outdir "$o" )
}
post() {
  rm -rf "$D/cmp_ab"
  "$GIPPA" compare "$D/out_a" "$D/out_b" --labels genes,sixframe --outdir "$D/cmp_ab" > "$D/compare.log" 2>&1
  echo "compare_exit=$?"; tail -5 "$D/compare.log"
  sha256sum "$D"/out_a/tables/* > "$D/out_a_tables.pre_report.sha256"
  for r in a b; do "$GIPPA" report "$D/out_$r" --redact > "$D/out_$r/report_redact.log" 2>&1; echo "report_$r exit $?"; done
  sha256sum -c --quiet "$D/out_a_tables.pre_report.sha256" && echo "tables unchanged by report"
}
te7() {
  run a "$D/out_a_te7" || return 1
  diff -r "$D/out_a/tables" "$D/out_a_te7/tables" > "$D/te7.diff" 2>&1; local rc=$?
  echo "te7_diff_exit=$rc files=$(ls "$D/out_a_te7/tables" | wc -l)"; return 0
}
cd "$D"
for s in "${STEPS[@]}"; do
  case "$s" in
    makedb) step makedb makedb ;;
    lib) step lib_a diann a -; step lib_b diann b - ;;
    diann_a) step diann_a diann a "$D/runs.txt" ;;
    diann_b) step diann_b diann b "$D/runs.txt" ;;
    conv_a) step conv_a conv a ;;
    conv_b) step conv_b conv b ;;
    run_a) step run_a run a "$D/out_a" ;;
    run_b) step run_b run b "$D/out_b" ;;
    post) step post post ;;
    te7) step te7 te7 ;;
    *) echo "unknown step $s" >&2; exit 2 ;;
  esac
done
echo "lane $ACC finished"
