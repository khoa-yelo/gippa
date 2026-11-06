#!/bin/bash
#SBATCH --job-name=gippa-example-build
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --exclude=barley-01,barley-02,barley-03,barley-04,rye-01,rye-02,oat-01,oat-02,oat-03,oat-04,oat-05,oat-06
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%j.out
# The release build of the quick-start example bundle (M12-01 spec AC 18, release step 6):
#   1. build_example.py into $GIPPA_ROOT/example (the release archive and work tree);
#   2. a second build at the same commit into $TMPDIR; both archives must have one sha256;
#   3. --check against validation/example/bundles.tsv when it has this version's row, otherwise
#      against a TSV made from build 1's printed row (proves --check reproduces the archive).
# The printout of build 1 (sha256, size, db_sources.yaml row, bundles.tsv row with this job's ID
# in build_job) is in this log; the job ID stays out of the bundle (rebuilds are byte-identical).
#
# Usage: sbatch scripts/slurm/build_example.sbatch [build_example.py args, e.g. --expected DIR
#        --expected-commit SHA, --force, --allow-dirty]
# GIPPA_REPO (default the main checkout) is the checkout whose scripts/ and src/ are used.
# --exclude: barley-01..04 do not mount the home directory, rye-01/02 stall start-up or fail
# "user env retrieval" (scripts/slurm/nightly.sbatch), and the oat nodes are kept free.
set -euo pipefail
export GIPPA_REPO="${GIPPA_REPO:-/home/classes/bios/270/khoa/giap}"
# before the source: farmshare_env.sh defaults GIPPA_ROOT to $GIPPA_REPO/local (a worktree's)
export GIPPA_ROOT="${GIPPA_ROOT:-/home/classes/bios/270/khoa/giap/local}"
source "$GIPPA_REPO/scripts/farmshare_env.sh"
export TMPDIR="$(mktemp -d /tmp/gippa.XXXXXX)"; trap 'rm -rf "$TMPDIR"' EXIT
export PYTHONPATH="$GIPPA_REPO/src"
PY="$UV_PROJECT_ENVIRONMENT/bin/python"
cd "$GIPPA_REPO"
echo "build_example: job ${SLURM_JOB_ID:-local} on $(hostname), commit $(git rev-parse HEAD)"

# drop --force from the second build and --check (they never replace an archive)
again=()
for x in "$@"; do [[ "$x" == --force ]] || again+=("$x"); done

echo "== build 1 (--out, default $GIPPA_ROOT/example)"
"$PY" scripts/example/build_example.py "$@" | tee "$TMPDIR/build1.log"
sha1=$(awk '/^sha256: /{print $2}' "$TMPDIR/build1.log")
archive=$(awk '/^archive: /{print $2}' "$TMPDIR/build1.log")

echo "== build 2 ($TMPDIR/second)"
"$PY" scripts/example/build_example.py "${again[@]}" --out "$TMPDIR/second/work" \
  --bundle-out "$TMPDIR/second/rel" | tee "$TMPDIR/build2.log"
sha2=$(awk '/^sha256: /{print $2}' "$TMPDIR/build2.log")
if [[ -z "$sha1" || "$sha1" != "$sha2" ]]; then
  echo "build_example: FAIL: the two builds differ (${sha1:-none} vs ${sha2:-none})" >&2
  exit 1
fi
echo "build_example: two builds agree: sha256 $sha1 ($archive)"

echo "== --check"
version=$(basename "$archive" .tar.gz); version=${version#gippa-example-}
tsv=validation/example/bundles.tsv
if ! awk -F'\t' -v v="$version" 'NR > 1 && $1 == v {found = 1} END {exit !found}' "$tsv"; then
  tsv="$TMPDIR/bundles.tsv"
  { head -n 1 validation/example/bundles.tsv
    awk -F'\t' -v v="$version" 'NF == 10 && $1 == v' "$TMPDIR/build1.log"; } > "$tsv"
  echo "build_example: $version has no row in validation/example/bundles.tsv; checking against build 1's row"
fi
"$PY" scripts/example/build_example.py "${again[@]}" --check --bundles-tsv "$tsv"
echo "build_example: done; add the printed rows (build 1) to db_sources.yaml and bundles.tsv"
