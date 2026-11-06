#!/bin/bash
#SBATCH --job-name=gippa-example-routes
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --exclude=barley-01,barley-02,barley-03,barley-04,rye-01,rye-02
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%j.out
# M12-01 slow route test (docs/specs/M12-01-example-dataset.md section 5, AC 17): fetch the example
# bundle offline (`gippa example fetch --from-file`), run the four routes' printed commands with
# real tools on this allocation's 4 CPUs (profile farmshare_local,apptainer) and compare each run
# with the bundle's expected/<route>.json. Never in CI or the session allocation (R14: 4 CPUs).
#
# --exclude: barley-01..04 do not mount the home directory, rye-01/02 stall start-up or fail
# "user env retrieval" (scripts/slurm/nightly.sbatch).
#
# Usage: sbatch scripts/slurm/example_routes.sbatch <repo-or-worktree-dir> [pytest args...]
#   env GIPPA_EXAMPLE_BUNDLE   the bundle file (default $GIPPA_ROOT/example/gippa-example-<X.Y>.tar.gz)
#   env GIPPA_EXAMPLE_VERSION  X.Y when it is not the installed GIPPA's MAJOR.MINOR (adds --version)
#   env GIPPA_EXAMPLE_DIR      measure mode (release step 6): an unpacked build tree; no fetch, no
#                              expected/ comparison, only runs/measured/<route>.json
# Outputs: $GIPPA_ROOT/example/routes/<jobid>/gippa-example/runs/ with each route's run,
#   wall_times.tsv (recorded, not asserted) and measured/<route>.json (input for expected/).
# The last log line is "PYTEST_EXIT=<code>". Rebuilding the bundle with expected/ is not done here.
set -uo pipefail
WT="${1:?repo or worktree dir}"; shift || true
export GIPPA_ROOT="${GIPPA_ROOT:-/home/classes/bios/270/khoa/giap/local}"   # before farmshare_env.sh
export GIPPA_REPO="$WT"
source "$WT/scripts/farmshare_env.sh"
unset PYTHONHOME
export TMPDIR="$(mktemp -d /tmp/gippa.XXXXXX)"; trap 'rm -rf "$TMPDIR"' EXIT
export NUMBA_CACHE_DIR="$TMPDIR/numba"
cd "$WT"
export PYTHONPATH="$PWD/src:$PWD/tests/sim:$PWD/tests"
export PATH="$UV_PROJECT_ENVIRONMENT/bin:$PATH"   # tasks call `gippa _step` from PATH
mem_gb=$(( ${SLURM_MEM_PER_NODE:-16384} / 1024 - 2 ))
export GIPPA_EXAMPLE_RUN_ARGS="${GIPPA_EXAMPLE_RUN_ARGS:---profile farmshare_local,apptainer --memory ${mem_gb}.GB}"
export GIPPA_EXAMPLE_OUT="${GIPPA_EXAMPLE_OUT:-$GIPPA_ROOT/example/routes/${SLURM_JOB_ID:-local}}"
export GIPPA_EXAMPLE_REQUIRE=1   # a missing bundle or antiSMASH DB fails instead of skipping
mkdir -p "$GIPPA_EXAMPLE_OUT"
echo "worktree=$WT branch=$(git rev-parse --abbrev-ref HEAD) head=$(git rev-parse --short HEAD) node=$(hostname)"
echo "out=$GIPPA_EXAMPLE_OUT run_args=$GIPPA_EXAMPLE_RUN_ARGS"
s=$SECONDS; nextflow -version >/dev/null 2>&1; echo "warm-up: nextflow rc=$? $((SECONDS - s))s"
[ "$#" -eq 0 ] && set -- tests/e2e/test_m12_01_example_routes.py -v -rA
set -- "$@" --require-capabilities needs_apptainer,needs_bakta_db,needs_promotech_models
"$UV_PROJECT_ENVIRONMENT/bin/python" -m pytest -p no:cacheprovider "$@"
rc=$?
[ -f "$GIPPA_EXAMPLE_OUT/gippa-example/runs/wall_times.tsv" ] && cat "$GIPPA_EXAMPLE_OUT/gippa-example/runs/wall_times.tsv"
echo "PYTEST_EXIT=$rc"
exit $rc
