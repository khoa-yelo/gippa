#!/bin/bash
#SBATCH --job-name=gippa-ncm-fetch
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=4000M
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%j.out
# Phase 1b step 2: the 3 PXD060434 wild-type exponential runs (run_selection.tsv, arm N1,N2) over
# https from ftp.pride.ebi.ac.uk (ADR-0077). PRIDE drops long transfers every 0.2-0.9 GB, so each
# file is resumed (--continue-at -) while every attempt makes progress (max 40 resumes; a resume
# that adds no byte counts as a failure, two in a row stop the file). The deposit has no checksum:
# the byte size must equal the PRIDE file list; sha256 goes to $D/raw/SHA256SUMS.
#   sbatch --exclude=... --time=08:00:00 validation/campaign/ncm3722/fetch_raw.sbatch
set -uo pipefail
export GIPPA_ROOT=/home/classes/bios/270/khoa/giap/local
source /home/classes/bios/270/khoa/giap/scripts/farmshare_env.sh
REPO="${CAMPAIGN_REPO:-$GIPPA_ROOT/worktrees/campaign-ncm3722}"
SEL="$REPO/validation/campaign/run_selection.tsv"
MAN="$REPO/validation/campaign/PXD060434/files.tsv"
RAW="${D:-$GIPPA_ROOT/campaign/ncm3722}/raw"
mkdir -p "$RAW"; cd "$RAW"
rc=0
for f in $(awk -F'\t' '$1 == "PXD060434" && $2 == "N1,N2" {print $4}' "$SEL"); do
  case "$f" in */* | *..* | *[!A-Za-z0-9._-]*) echo "unsafe $f" >&2; rc=2; continue ;; esac
  url=$(awk -F'\t' -v f="$f" '$1 == f {print $7}' "$MAN"); bytes=$(awk -F'\t' -v f="$f" '$1 == f {print $4}' "$MAN")
  case "$url" in https://ftp.pride.ebi.ac.uk/*) ;; *) echo "bad url $url" >&2; rc=2; continue ;; esac
  stalls=0
  for i in $(seq 1 40); do
    have=$(stat -c %s "$f" 2>/dev/null || echo 0)
    [ "$have" = "$bytes" ] && break
    [ "$have" -gt "$bytes" ] && { echo "$f larger than listed; removing" >&2; rm -f "$f"; have=0; }
    curl --fail --location --silent --show-error --proto '=https' --proto-redir '=https' \
      --continue-at - --max-time 7200 --output "$f" "$url"; c=$?
    now=$(stat -c %s "$f" 2>/dev/null || echo 0)
    echo "$(date -u +%FT%TZ) $f attempt $i curl=$c bytes $have -> $now of $bytes"
    if [ "$now" -le "$have" ]; then stalls=$((stalls + 1)); [ $stalls -ge 2 ] && break; else stalls=0; fi
  done
  if [ "$(stat -c %s "$f" 2>/dev/null || echo 0)" = "$bytes" ]; then
    grep -q "  $f\$" SHA256SUMS 2>/dev/null || sha256sum "$f" >> SHA256SUMS
    echo "ok $f $bytes"
  else
    echo "FAILED $f" >&2; rc=1
  fi
done
cat SHA256SUMS 2>/dev/null
exit $rc
