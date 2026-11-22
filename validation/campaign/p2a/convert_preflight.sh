#!/bin/bash







# Standard recipe steps 5-6 for one arm (generalises scripts/slurm/diann_convert.sbatch):
# gippa convert --from diann (twice; the two outputs must be byte-identical), then
# gippa preflight --strict, which must exit 0 before any gippa run is submitted.
#
#   sbatch convert_preflight.sbatch PIN REPORT.parquet IN_DIR GENOME BAKTA_JSON [same|other]
#
# Writes IN_DIR/peptides.tsv.gz, abundance.tsv.gz, convert.log, preflight/, preflight.log and
# IN_DIR/convert_provenance.json (commit, sha256s, exit codes, wall times). Exit 0 all ok.
set -uo pipefail
source "${GIPPA_REPO:-/home/classes/bios/270/khoa/giap}/scripts/farmshare_env.sh"
PIN="$1" REPORT="$2" IN="$3" GENOME="$4" JSON="$5" SDB="${6:-same}"
export TMPDIR="$(mktemp -d /tmp/gippa.XXXXXX)"; trap 'rm -rf "$TMPDIR"' EXIT
case "$IN" in "$GIPPA_ROOT"/campaign/*) ;; *) echo "IN_DIR must be under \$GIPPA_ROOT/campaign" >&2; exit 2 ;; esac
for f in "$REPORT" "$GENOME" "$JSON"; do [ -s "$f" ] || { echo "missing $f" >&2; exit 2; }; done
mkdir -p "$IN"; rm -rf "$IN/preflight"
g() { PYTHONPATH="$PIN/src" "$UV_PROJECT_ENVIRONMENT/bin/python" -c 'from gippa.cli import cli; cli()' "$@"; }
s=$(date +%s)
g convert --from diann "$REPORT" --out "$IN/peptides.tsv.gz" --abundance-out "$IN/abundance.tsv.gz" \
  --search-db "$SDB" >"$IN/convert.log" 2>&1
C1=$?; C1_S=$(( $(date +%s) - s ))
g convert --from diann "$REPORT" --out "$TMPDIR/p.tsv.gz" --abundance-out "$TMPDIR/a.tsv.gz" \
  --search-db "$SDB" >"$TMPDIR/c2.log" 2>&1
C2=$?
DET=false
cmp -s "$IN/peptides.tsv.gz" "$TMPDIR/p.tsv.gz" && cmp -s "$IN/abundance.tsv.gz" "$TMPDIR/a.tsv.gz" && DET=true
tail -15 "$IN/convert.log"
s=$(date +%s)
g preflight --genome "$GENOME" --annotation "$JSON" --peptides "$IN/peptides.tsv.gz" \
  --outdir "$IN/preflight" --strict >"$IN/preflight.log" 2>&1
PF=$?; PF_S=$(( $(date +%s) - s ))
tail -25 "$IN/preflight.log"
echo "convert=$C1 convert2=$C2 deterministic=$DET preflight_strict=$PF"
python3 -I - "$IN/convert_provenance.json" <<PY
import json, sys, hashlib, os
def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()
d = {"job_id": os.environ.get("SLURM_JOB_ID"), "pin": "$PIN",
     "commit": os.popen("git -C $PIN rev-parse HEAD").read().strip(),
     "search_db": "$SDB", "convert_exit": $C1, "convert2_exit": $C2, "convert_wall_s": $C1_S,
     "deterministic": "$DET" == "true", "preflight_strict_exit": $PF, "preflight_wall_s": $PF_S,
     "sha256": {k: sha(p) for k, p in (("report", "$REPORT"), ("genome", "$GENOME"),
                ("annotation", "$JSON"), ("peptides", "$IN/peptides.tsv.gz"),
                ("abundance", "$IN/abundance.tsv.gz")) if os.path.exists(p)}}
json.dump(d, open(sys.argv[1], "w"), indent=1)
PY
[ "$C1" -eq 0 ] && [ "$C2" -eq 0 ] && [ "$DET" = true ] && [ "$PF" -eq 0 ]
