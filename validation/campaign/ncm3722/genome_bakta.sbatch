#!/bin/bash
#SBATCH --job-name=gippa-ncm-bakta
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=64000M
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%j.out
# Phase 1b step 1: NCM3722 genome GCF_001043215.1 from NCBI (MD5 against md5checksums.txt; both
# replicons tagged [topology=circular], NCBI esummary says circular/complete), then Bakta 1.12.2
# full DB 6.0 with the command of gippa.annotate.bakta_cmd.bakta_command (--keep-contig-headers)
# in the pinned image, as scripts/slurm/record_bakta_full.sbatch does.
#   sbatch --exclude=... --time=06:00:00 validation/campaign/ncm3722/genome_bakta.sbatch
# Env: PIN (pinned worktree), D (default $GIPPA_ROOT/campaign/ncm3722).
set -euo pipefail
export GIPPA_ROOT=/home/classes/bios/270/khoa/giap/local
PIN="${PIN:-$GIPPA_ROOT/worktrees/campaign-pin}"
source /home/classes/bios/270/khoa/giap/scripts/farmshare_env.sh
export PYTHONPATH="$PIN/src"
PY="$UV_PROJECT_ENVIRONMENT/bin/python"
GIPPA="$UV_PROJECT_ENVIRONMENT/bin/gippa"
D="${D:-$GIPPA_ROOT/campaign/ncm3722}"
A=GCF_001043215.1_ASM104321v1
URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/043/215/$A"
mkdir -p "$D/genome" "$D/bakta"
cd "$D/genome"
get() { curl --fail --silent --show-error --proto '=https' --retry 1 -o "$1" "$URL/$1"; }
get md5checksums.txt
for f in ${A}_genomic.fna.gz ${A}_genomic.gff.gz; do
  want=$(awk -v f="./$f" '$2 == f {print $1}' md5checksums.txt)
  [ -f "$f" ] && [ "$(md5sum "$f" | cut -d' ' -f1)" = "$want" ] || get "$f"
  got=$(md5sum "$f" | cut -d' ' -f1)
  [ "$got" = "$want" ] || { echo "MD5 mismatch $f: $got != $want" >&2; exit 1; }
  echo "md5 ok $f $got"
done
# both replicons circular (NCBI esummary topology=circular, esummary_topology.json)
zcat ${A}_genomic.fna.gz | awk '/^>/ {print $0 " [topology=circular]"; next} {print}' > ncm3722.fna
grep '^>' ncm3722.fna
gzip -9nf ncm3722.fna
sha256sum ${A}_genomic.fna.gz ncm3722.fna.gz > genome.sha256
GENOME="$D/genome/ncm3722.fna.gz"

cd "$PIN"
DB="$GIPPA_DB_DIR/bakta/db"
DIGEST=$(grep -oP "withName: 'BAKTA' \\{ container = '(docker://)?oschwengers/bakta@sha256:\\K[0-9a-f]{64}" \
  workflow/conf/containers.config)
IMG="$GIPPA_ROOT/containers/oschwengers-bakta@sha256-$DIGEST.img"
[ -f "$IMG" ] || { echo "missing $IMG" >&2; exit 1; }
grep -q '"type": "full"' "$DB/version.json" || { echo "full Bakta DB missing at $DB" >&2; exit 1; }
d="$D/bakta"
[ -f "$d/.done" ] && { echo "bakta done"; exit 0; }
rm -rf "$d/preflight" "$d/out"
printf '>placeholder_evidence\nMSTNPKPQRKTKRNTNRRPQDVKFPGG\n' > "$d/placeholder.faa"
"$GIPPA" preflight --genome "$GENOME" --proteins "$d/placeholder.faa" --outdir "$d/preflight"
w="$d/preflight/work/preflight"
CMD=(); while IFS= read -r a || [ -n "$a" ]; do CMD+=("$a"); done < <("$PY" - "$w/genome.norm.fa" "$d/out" ncm3722 "$DB" "$w/replicons.tsv" <<'PYEOF'
import sys
from pathlib import Path
from gippa.annotate.bakta_cmd import bakta_command
g, o, p, db, rep = sys.argv[1:]
for a in bakta_command(genome=Path(g), outdir=Path(o), prefix=p, db=Path(db), threads=16,
                       all_circular=False, complete=False, replicons_tsv=Path(rep), regions_gff3=None,
                       extra_args=[], translation_table=11):
    print(a)
PYEOF
)
printf '%q ' "${CMD[@]}" > "$d/command.txt"; echo >> "$d/command.txt"
cat "$d/command.txt"; cat "$w/replicons.tsv"
t0=$(date +%s)
apptainer exec --bind "$GIPPA_DB_DIR,$d" "$IMG" "${CMD[@]}"
echo "wall_s=$(( $(date +%s) - t0 )) cpus=16 node=$(hostname -s) job=${SLURM_JOB_ID:-none}" > "$d/runtime.txt"
{ apptainer exec "$IMG" bakta --version; echo "image sha256:$DIGEST"; cat "$DB/version.json"; } > "$d/versions.txt" 2>&1
(cd "$d" && find out -type f -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS)
touch "$d/.done"; cat "$d/runtime.txt"
