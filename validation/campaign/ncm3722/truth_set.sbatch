#!/bin/bash
#SBATCH --job-name=gippa-ncm-truth
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=16000M
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%j.out
# Phase 1b step 4: truth set (truth_set.py) with BLAST+ blastp from the pinned GIPPA BLAST image.
#   sbatch --exclude=... --time=02:00:00 validation/campaign/ncm3722/truth_set.sbatch
set -euo pipefail
export GIPPA_ROOT=/home/classes/bios/270/khoa/giap/local
source /home/classes/bios/270/khoa/giap/scripts/farmshare_env.sh
HERE="${CAMPAIGN_REPO:-$GIPPA_ROOT/worktrees/campaign-ncm3722}/validation/campaign/ncm3722"
PIN="$GIPPA_ROOT/worktrees/campaign-pin"
PY="$UV_PROJECT_ENVIRONMENT/bin/python"
D="${D:-$GIPPA_ROOT/campaign/ncm3722}"; T="$D/truth"; mkdir -p "$T"
NJ="$D/bakta/out/ncm3722.json"; MJ="$GIPPA_ROOT/recordings/bakta_full_circ/k12_mg1655/bakta/k12_mg1655.json"
DIG=$(grep -oP "blast@sha256:\K[0-9a-f]{64}" "$PIN/workflow/conf/containers.config" | head -1)
IMG="$GIPPA_ROOT/containers/quay.io-biocontainers-blast@sha256-$DIG.img"
[ -f "$IMG" ] || { echo "missing $IMG" >&2; exit 2; }
b() { apptainer exec -B "$T" "$IMG" "$@"; }
"$PY" -I "$HERE/truth_set.py" proteins "$NJ" "$MJ" "$T"
b blastp -version | head -1 > "$T/blast_version.txt"; echo "image sha256:$DIG" >> "$T/blast_version.txt"
b makeblastdb -in "$T/ncm.faa" -dbtype prot -out "$T/ncm" >/dev/null
b makeblastdb -in "$T/mg.faa" -dbtype prot -out "$T/mg" >/dev/null
F="6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore qseq sseq"
b blastp -query "$T/ncm.faa" -db "$T/mg" -seg no -evalue 1e-10 -max_hsps 1 -max_target_seqs 5 -num_threads 4 -outfmt "$F" -out "$T/n2m.tsv"
b blastp -query "$T/mg.faa" -db "$T/ncm" -seg no -evalue 1e-10 -max_hsps 1 -max_target_seqs 5 -num_threads 4 -outfmt "$F" -out "$T/m2n.tsv"
"$PY" -I "$HERE/truth_set.py" sites "$T" "$NJ" "$MJ"
sha256sum "$NJ" "$MJ" "$T/truth_sites.tsv" > "$T/sha256.txt"
cat "$T/blast_version.txt"
