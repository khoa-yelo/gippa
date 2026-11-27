#!/bin/bash
#SBATCH --job-name=gippa-campaign-fetch
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=2000M
#SBATCH --time=12:00:00
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%j.out
# Resumable download of selected raw files of one campaign manifest (ADR-0077 + Amendment 1).
#
#   sbatch --time=12:00:00 --exclude=barley-01,barley-02,barley-03,barley-04,rye-01,rye-02 \
#     validation/campaign/slurm/fetch_runs.sbatch MANIFEST DEST_DIR FILE_NAME [FILE_NAME ...]
#
# MANIFEST is a validation/campaign/<accession>/files.tsv (file_name, bytes, checksum, url, ...).
# Each file lands as DEST_DIR/<basename>. PRIDE's HTTPS server drops long transfers every
# 0.2-0.9 GB (PROGRESS.md > Blocked, 2026-11-19), so one "download" is a chain of range
# requests (curl --continue-at -) on the same partial file: the transfer is resumed while it
# makes progress, and declared blocked after 3 consecutive resumes that add no byte (R6: one
# retry of a failed transfer, not of every dropped connection). Verified: byte size, and the
# deposit's checksum when the manifest has one (sha1:/md5:); sha256 always recorded in
# DEST_DIR/SHA256SUMS. Hosts: ftp.pride.ebi.ac.uk, panoramaweb.org and download.iprox.cn (iProX,
# ADR-0077 decision 1) over https only.
# Exit 0 all verified; 1 at least one file blocked or failed verification; 2 precondition.
set -uo pipefail
source "${GIPPA_REPO:-/home/classes/bios/270/khoa/giap}/scripts/farmshare_env.sh"

[ $# -ge 3 ] || { echo "usage: fetch_runs.sbatch MANIFEST DEST_DIR FILE..." >&2; exit 2; }
MANIFEST="$1" DEST="$2"; shift 2
[ -f "$MANIFEST" ] || { echo "no manifest $MANIFEST" >&2; exit 2; }
case "$DEST" in "$GIPPA_ROOT"/campaign/*) ;; *) echo "DEST must be under \$GIPPA_ROOT/campaign (R11)" >&2; exit 2 ;; esac
mkdir -p "$DEST"
MAX_STALLS=3 MAX_RESUMES=200

verify() {  # path bytes checksum
  [ -f "$1" ] && [ "$(stat -c %s "$1")" = "$2" ] || return 1
  case "$3" in
    sha1:*) [ "$(sha1sum "$1" | cut -d' ' -f1)" = "${3#sha1:}" ] ;;
    md5:*) [ "$(md5sum "$1" | cut -d' ' -f1)" = "${3#md5:}" ] ;;
    *) true ;;
  esac
}

failed=0
for FILE in "$@"; do
  row="$(awk -F'\t' -v f="$FILE" '$1 == f' "$MANIFEST")"
  [ -n "$row" ] || { echo "$FILE is not in $MANIFEST" >&2; failed=1; continue; }
  BYTES="$(cut -f4 <<<"$row")" CK="$(cut -f5 <<<"$row")" URL="$(cut -f7 <<<"$row")"
  case "$URL" in
    https://ftp.pride.ebi.ac.uk/* | https://panoramaweb.org/_webdav/* | https://download.iprox.cn/*) ;;
    *) echo "not an ADR-0077 https URL: $URL" >&2; failed=1; continue ;;
  esac
  NAME="$(basename "$FILE")"
  case "$NAME" in *[!A-Za-z0-9._+-]*) echo "unsafe file name: $NAME" >&2; failed=1; continue ;; esac
  OUTF="$DEST/$NAME"
  echo "== $NAME ($BYTES bytes, ${CK}) <- $URL"
  if verify "$OUTF" "$BYTES" "$CK"; then
    echo "already verified"
  else
    stalls=0 resumes=0 t0=$(date +%s)
    while :; do
      before=0; [ -f "$OUTF" ] && before=$(stat -c %s "$OUTF")
      [ "$before" -lt "$BYTES" ] || break
      curl --fail --location --silent --show-error --proto '=https' --proto-redir '=https' \
        --continue-at - --max-time 7200 --speed-limit 10000 --speed-time 300 \
        --user-agent gippa-campaign-fetch/1 --output "$OUTF" "$URL"
      rc=$?
      after=0; [ -f "$OUTF" ] && after=$(stat -c %s "$OUTF")
      resumes=$((resumes + 1))
      echo "  resume $resumes: curl exit $rc, $before -> $after bytes"
      if [ "$after" -le "$before" ]; then stalls=$((stalls + 1)); else stalls=0; fi
      [ "$stalls" -lt "$MAX_STALLS" ] && [ "$resumes" -lt "$MAX_RESUMES" ] || break
      [ "$rc" -eq 0 ] || sleep 15
    done
    echo "  transfer: $resumes request(s), $(( $(date +%s) - t0 )) s"
    if ! verify "$OUTF" "$BYTES" "$CK"; then
      echo "BLOCKED or corrupt: $NAME (size $( [ -f "$OUTF" ] && stat -c %s "$OUTF" || echo 0) of $BYTES; checksum $CK)" >&2
      # a full-length file that fails the checksum is discarded; a short partial is kept for a resume
      if [ -f "$OUTF" ] && [ "$(stat -c %s "$OUTF")" -ge "$BYTES" ]; then rm -f "$OUTF"; fi
      failed=1
      continue
    fi
  fi
  sha="$(sha256sum "$OUTF" | cut -d' ' -f1)"
  echo "verified ${CK}; sha256 $sha"
  { grep -v "  $NAME\$" "$DEST/SHA256SUMS" 2>/dev/null; printf '%s  %s\n' "$sha" "$NAME"; } \
    | LC_ALL=C sort -k2 >"$DEST/SHA256SUMS.tmp.$$" && mv "$DEST/SHA256SUMS.tmp.$$" "$DEST/SHA256SUMS"
done
exit "$failed"
