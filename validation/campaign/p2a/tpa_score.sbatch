#!/bin/bash
#SBATCH --job-name=gippa-camp-p2a-score
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=32000M
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%j.out
# PXD065022 scoring of the _fix2 runs (pin rel111c-pin, main 38e1cdba). Per strain TAG given:
#   gippa compare out_a_fix2 out_b_fix2 (when route B finished) -> TAG/cmp_ab_fix2
#   p2_metrics.py --host (route A, B, compare, answer loci)    -> score_fix2/TAG/p2.json
#   provenance.py a / b / x_tpa_dal1 (/ te7 for DAL-1)          -> score_fix2/TAG/prov_*.json
#   non-DAL-1: score.py (T6 truth vs cross-strain run)          -> score_fix2/TAG/t6/score.json
#              placement_fdr.py (pooled TE8 estimator)          -> score_fix2/TAG/placement_fdr.json
#   DAL-1: diff -r of tables/ against the TE7 rerun             -> score_fix2/tpa_dal1/te7_diff.txt
# Reads GIPPA outputs only; writes under $GIPPA_ROOT. Never edits GIPPA code or settings.
#   sbatch --exclude=... --time=02:00:00 tpa_score.sbatch tpa_dal1 tpa_grady ...
set -uo pipefail
export GIPPA_ROOT=/home/classes/bios/270/khoa/giap/local
source /home/classes/bios/270/khoa/giap/scripts/farmshare_env.sh
PIN=$GIPPA_ROOT/worktrees/rel111c-pin
HERE=$GIPPA_ROOT/worktrees/campaign-p2a/validation/campaign/p2a
KEY=$HERE/../PXD065022
C=$GIPPA_ROOT/campaign/PXD065022
PY="$UV_PROJECT_ENVIRONMENT/bin/python"
SUF=_fix2
export XDG_CACHE_HOME="$GIPPA_ROOT/.cache" MPLBACKEND=Agg
rc_all=0
step() { local name="$1"; shift; local s; s=$(date +%s); "$@"; local rc=$?
  echo "STEP $name exit=$rc wall_s=$(( $(date +%s) - s ))"; [ $rc -eq 0 ] || rc_all=1; }
for t in "$@"; do
  D=$C/$t S=$C/score$SUF/$t; mkdir -p "$S"
  A=$D/out_a$SUF B=$D/out_b$SUF
  gff=$(ls "$C"/genome/$t/GCF_*_genomic.gff.gz)
  extra=()
  if [ -s "$B/run_manifest.json" ]; then
    if [ ! -s "$D/cmp_ab$SUF/tables/compare_abundance.tsv" ]; then
      step "$t compare" env PYTHONPATH="$PIN/src" "$UV_PROJECT_ENVIRONMENT/bin/gippa" compare "$A" "$B" \
        --labels genes,sixframe --outdir "$D/cmp_ab$SUF"
    fi
    extra=(--b "$B" --cmp "$D/cmp_ab$SUF")
    step "$t prov b" "$PY" -I "$HERE/provenance.py" "$D" b "$S/prov_b.json" "$SUF"
  fi
  step "$t p2_metrics" env PYTHONPATH="$PIN:$PIN/src" "$PY" "$HERE/p2_metrics.py" --a "$A" \
    --pep-a "$D/in_a$SUF/peptides.tsv.gz" "${extra[@]}" --bakta-json "$C/bakta/$t/out/$t.json" \
    --genome "$C/genome/$t/$t.fna.gz" --refseq-gff "$gff" \
    --contaminants "$GIPPA_DB_DIR/contaminants/contaminants.fasta" --loci "$KEY/answer_loci.tsv" \
    --tag "$t" --host --out "$S/p2.json"
  step "$t prov a" "$PY" -I "$HERE/provenance.py" "$D" a "$S/prov_a.json" "$SUF"
  if [ "$t" = tpa_dal1 ]; then
    step "$t prov te7" "$PY" -I "$HERE/provenance.py" "$D" te7 "$S/prov_te7.json" "$SUF"
    diff -r "$A/tables" "$D/out_a${SUF}_te7/tables" > "$S/te7_diff.txt" 2>&1
    echo "TE7 diff exit=$? lines=$(wc -l < "$S/te7_diff.txt")"
    (cd "$A/tables" && find . -type f | sort | xargs sha256sum) > "$S/te7_tables_a.sha256"
    (cd "$D/out_a${SUF}_te7/tables" && find . -type f | sort | xargs sha256sum) > "$S/te7_tables_te7.sha256"
  else
    X=$D/out_x_tpa_dal1$SUF
    step "$t prov x" "$PY" -I "$HERE/provenance.py" "$D" x_tpa_dal1 "$S/prov_x.json" "$SUF"
    step "$t score" "$PY" -I "$HERE/score.py" "$C/truth/$t" "$A" "$X" "$S/t6" > "$S/t6.log" 2>&1
    step "$t placement_fdr" env PYTHONPATH="$PIN/src:$PIN/validation/proxies" "$PY" "$HERE/placement_fdr.py" \
      "$D/in_x_tpa_dal1$SUF/peptides.tsv.gz" "$C/genome/tpa_dal1/tpa_dal1.fna.gz" \
      "$C/bakta/tpa_dal1/out/tpa_dal1.json" "$S/placement_fdr.json"
  fi
done
exit $rc_all
