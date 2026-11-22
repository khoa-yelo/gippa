#!/bin/bash
#SBATCH --job-name=gippa-camp-p2a-truth
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=16000M
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%j.out
# PXD065022 T6 truth sets (ANSWER_KEY.md): each strain's Bakta proteins vs DAL-1's, with
# truth_set.py (Phase 1b code, unchanged: "ncm" = the strain, "mg" = the DAL-1 reference) and
# blastp from the pinned GIPPA BLAST image. Reads annotations only.
#   sbatch --exclude=... --time=02:00:00 truth.sbatch
set -euo pipefail
export GIPPA_ROOT=/home/classes/bios/270/khoa/giap/local
source /home/classes/bios/270/khoa/giap/scripts/farmshare_env.sh
HERE=$GIPPA_ROOT/worktrees/campaign-p2a/validation/campaign/p2a
PIN=$GIPPA_ROOT/worktrees/rel111-pin
PY="$UV_PROJECT_ENVIRONMENT/bin/python"
C=$GIPPA_ROOT/campaign/PXD065022
MJ=$C/bakta/tpa_dal1/out/tpa_dal1.json
DIG=$(grep -oP "blast@sha256:\K[0-9a-f]{64}" "$PIN/workflow/conf/containers.config" | head -1)
IMG="$GIPPA_ROOT/containers/quay.io-biocontainers-blast@sha256-$DIG.img"
[ -f "$IMG" ] || { echo "missing $IMG" >&2; exit 2; }
F="6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore qseq sseq"
for s in tpa_grady tpa_haitib tpa_madras tpa_mexicoa tpa_phil1 tpa_ss14; do
  T=$C/truth/$s; mkdir -p "$T"
  NJ=$C/bakta/$s/out/$s.json
  b() { apptainer exec -B "$T" "$IMG" "$@"; }
  "$PY" -I "$HERE/truth_set.py" proteins "$NJ" "$MJ" "$T"
  b blastp -version | head -1 > "$T/blast_version.txt"; echo "image sha256:$DIG" >> "$T/blast_version.txt"
  b makeblastdb -in "$T/ncm.faa" -dbtype prot -out "$T/ncm" >/dev/null
  b makeblastdb -in "$T/mg.faa" -dbtype prot -out "$T/mg" >/dev/null
  b blastp -query "$T/ncm.faa" -db "$T/mg" -seg no -evalue 1e-10 -max_hsps 1 -max_target_seqs 5 -num_threads 4 -outfmt "$F" -out "$T/n2m.tsv"
  b blastp -query "$T/mg.faa" -db "$T/ncm" -seg no -evalue 1e-10 -max_hsps 1 -max_target_seqs 5 -num_threads 4 -outfmt "$F" -out "$T/m2n.tsv"
  echo "== $s"; "$PY" -I "$HERE/truth_set.py" sites "$T" "$NJ" "$MJ"
  sha256sum "$NJ" "$MJ" "$T/truth_sites.tsv" > "$T/sha256.txt"
done
