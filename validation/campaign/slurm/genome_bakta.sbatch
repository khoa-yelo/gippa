#!/bin/bash
#SBATCH --job-name=gippa-camp-bakta
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=64000M
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%A_%a.out
# Campaign Phase 2 prep: one genome per array task (row N of validation/campaign/genomes.tsv).
# NCBI RefSeq assembly FASTA + GFF (https, MD5 against md5checksums.txt), replicons tagged
# [topology=circular] where NCBI nuccore says circular (genomes.tsv column `circular`), then Bakta
# 1.12.2 full DB 6.0 with gippa.annotate.bakta_cmd.bakta_command (--keep-contig-headers) in the
# pinned image: the same steps as Phase 1b job 1998570 (validation/campaign/ncm3722/genome_bakta.sbatch).
#   sbatch --exclude=... --time=06:00:00 --array=1-16%2 validation/campaign/slurm/genome_bakta.sbatch
# Output: $GIPPA_ROOT/campaign/<dataset>/genome/<label>/ and .../bakta/<label>/
# Env: PIN (pinned worktree, default campaign-pin), REPO (worktree holding genomes.tsv).
set -euo pipefail
export GIPPA_ROOT=/home/classes/bios/270/khoa/giap/local
PIN="${PIN:-$GIPPA_ROOT/worktrees/campaign-pin}"
REPO="${REPO:-$GIPPA_ROOT/worktrees/campaign-phase0}"
source /home/classes/bios/270/khoa/giap/scripts/farmshare_env.sh
export PYTHONPATH="$PIN/src"
PY="$UV_PROJECT_ENVIRONMENT/bin/python"
GIPPA="$UV_PROJECT_ENVIRONMENT/bin/gippa"
N="${SLURM_ARRAY_TASK_ID:?array task id}"
row=$(awk -F'\t' -v n=$((N + 1)) 'NR == n' "$REPO/validation/campaign/genomes.tsv")
[ -n "$row" ] || { echo "no row $N" >&2; exit 2; }
DS=$(cut -f1 <<<"$row"); LABEL=$(cut -f2 <<<"$row"); ASM=$(cut -f3 <<<"$row")
URL=$(cut -f7 <<<"$row"); CIRC=$(cut -f9 <<<"$row")
case "$URL" in https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/*) ;; *) echo "bad url $URL" >&2; exit 2 ;; esac
case "$LABEL" in *[!a-z0-9_]*) echo "unsafe label" >&2; exit 2 ;; esac
A=$(basename "$URL")
echo "task $N: $DS $LABEL $A circular=$CIRC"
D="$GIPPA_ROOT/campaign/$DS"
G="$D/genome/$LABEL"; d="$D/bakta/$LABEL"
mkdir -p "$G" "$d"
cd "$G"
get() { curl --fail --silent --show-error --proto '=https' --retry 1 -o "$1" "$URL/$1"; }
get md5checksums.txt
for f in ${A}_genomic.fna.gz ${A}_genomic.gff.gz; do
  want=$(awk -v f="./$f" '$2 == f {print $1}' md5checksums.txt)
  [ -n "$want" ] || { echo "$f not in md5checksums.txt" >&2; exit 1; }
  [ -f "$f" ] && [ "$(md5sum "$f" | cut -d' ' -f1)" = "$want" ] || get "$f"
  got=$(md5sum "$f" | cut -d' ' -f1)
  [ "$got" = "$want" ] || { echo "MD5 mismatch $f: $got != $want" >&2; exit 1; }
  echo "md5 ok $f $got"
done
zcat ${A}_genomic.fna.gz | awk -v circ=",$CIRC," '
  /^>/ { id = substr($1, 2); if (index(circ, "," id ",")) print $0 " [topology=circular]"; else print; next }
  { print }' > "$LABEL.fna"
grep '^>' "$LABEL.fna"
gzip -9nf "$LABEL.fna"
sha256sum ${A}_genomic.fna.gz ${A}_genomic.gff.gz "$LABEL.fna.gz" > genome.sha256
GENOME="$G/$LABEL.fna.gz"

cd "$PIN"
DB="$GIPPA_DB_DIR/bakta/db"
DIGEST=$(grep -oP "withName: 'BAKTA' \\{ container = '(docker://)?oschwengers/bakta@sha256:\\K[0-9a-f]{64}" \
  workflow/conf/containers.config)
IMG="$GIPPA_ROOT/containers/oschwengers-bakta@sha256-$DIGEST.img"
[ -f "$IMG" ] || { echo "missing $IMG" >&2; exit 1; }
grep -q '"type": "full"' "$DB/version.json" || { echo "full Bakta DB missing at $DB" >&2; exit 1; }
[ -f "$d/.done" ] && { echo "bakta done"; exit 0; }
rm -rf "$d/preflight" "$d/out"
printf '>placeholder_evidence\nMSTNPKPQRKTKRNTNRRPQDVKFPGG\n' > "$d/placeholder.faa"
"$GIPPA" preflight --genome "$GENOME" --proteins "$d/placeholder.faa" --outdir "$d/preflight"
w="$d/preflight/work/preflight"
CMD=(); while IFS= read -r a || [ -n "$a" ]; do CMD+=("$a"); done < <("$PY" - "$w/genome.norm.fa" "$d/out" "$LABEL" "$DB" "$w/replicons.tsv" <<'PYEOF'
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
echo "wall_s=$(( $(date +%s) - t0 )) cpus=16 node=$(hostname -s) job=${SLURM_ARRAY_JOB_ID:-none}_${N}" > "$d/runtime.txt"
{ apptainer exec "$IMG" bakta --version; echo "image sha256:$DIGEST"; cat "$DB/version.json"; } > "$d/versions.txt" 2>&1
(cd "$d" && find out -type f -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS)
touch "$d/.done"; cat "$d/runtime.txt"
