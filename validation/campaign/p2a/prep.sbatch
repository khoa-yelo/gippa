#!/bin/bash
#SBATCH --job-name=gippa-camp-p2a-prep
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=32000M
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%j.out
# Campaign Phase 2, group A: inputs that need no search (no spectrum is read here).
#   1. unpack every Bruker .d archive of PXD065022 and PXD055395 to $D/raw_d/<run>.d
#      (zip CRC test of every member; archives that hold the .d contents at their root, as the
#      PXD055395 ones do, are unpacked into <run>.d/; no member may be absolute or contain "..");
#   2. fetch the rabbit reference proteome (UniProt UP000001811, taxon 9986) over HTTPS from
#      ftp.uniprot.org, check its MD5 against the release's RELEASE.metalink, and write it with
#      headers CONT_RABIT_<accession> so that host peptides land in GIPPA's contaminant bucket
#      (input.contaminant_prefixes default includes CONT_; pre-registered in ANSWER_KEY.md);
#   3. gippa makedb route A and route B (--six-frame-orfs --min-orf-aa 20) for every genome.
#   sbatch --exclude=... --time=06:00:00 validation/campaign/p2a/prep.sbatch
set -uo pipefail
export GIPPA_ROOT=/home/classes/bios/270/khoa/giap/local
PIN="${PIN:-$GIPPA_ROOT/worktrees/rel111-pin}"
source /home/classes/bios/270/khoa/giap/scripts/farmshare_env.sh
export PYTHONPATH="$PIN/src"
GIPPA="$UV_PROJECT_ENVIRONMENT/bin/gippa"
C="$GIPPA_ROOT/campaign"
rc_all=0

echo "== 1. unpack .d archives"
unpack() {  # archive dest_dir
  python3 -I - "$1" "$2" <<'PY'
import sys, zipfile, pathlib, shutil
arc, dest = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2])
run = arc.name
for suf in (".zip", ".Zip", ".ZIP"):
    run = run.removesuffix(suf)
run = run.removesuffix(".d")
out = dest / f"{run}.d"
if (out / ".unpacked").exists():
    print(f"{run}: already unpacked"); sys.exit(0)
z = zipfile.ZipFile(arc)
names = z.namelist()
for n in names:
    if n.startswith("/") or ".." in pathlib.PurePosixPath(n).parts:
        sys.exit(f"unsafe member {n}")
bad = z.testzip()
if bad:
    sys.exit(f"corrupt member {bad}")
tops = {n.split("/", 1)[0] for n in names}
strip = None
if len(tops) == 1 and next(iter(tops)).endswith(".d"):
    strip = next(iter(tops)) + "/"
if out.exists():
    shutil.rmtree(out)
out.mkdir(parents=True)
n_files = 0
for info in z.infolist():
    name = info.filename
    if strip:
        name = name[len(strip):]
    if not name or name.endswith("/"):
        continue
    # a nested "<run>.d/<run>.d" zero-byte file exists in some PXD065022 archives; keep it
    target = out / name
    target.parent.mkdir(parents=True, exist_ok=True)
    with z.open(info) as src, open(target, "wb") as dst:
        shutil.copyfileobj(src, dst, 1 << 20)
    n_files += 1
if not (out / "analysis.tdf").exists() or not (out / "analysis.tdf_bin").exists():
    sys.exit(f"{run}: no analysis.tdf/analysis.tdf_bin after unpacking")
(out / ".unpacked").write_text(f"{n_files}\n")
print(f"{run}: {n_files} files")
PY
}
export -f unpack
for acc in PXD065022 PXD055395; do
  mkdir -p "$C/$acc/raw_d"
  ls "$C/$acc"/raw/*.d.[zZ]ip | xargs -P 4 -I{} bash -c 'unpack "$1" "$2"' _ {} "$C/$acc/raw_d" || rc_all=1
  echo "$acc: $(ls -d "$C/$acc"/raw_d/*.d | wc -l) .d folders"
done

echo "== 2. rabbit reference proteome"
H="$C/PXD065022/host"; mkdir -p "$H"
BASE=https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes
F=UP000001811_9986.fasta.gz
get() { curl --fail --location --silent --show-error --proto '=https' --proto-redir '=https' --max-time 1800 -o "$2" "$1"; }
ok=0
for attempt in 1 2; do
  get "$BASE/Eukaryota/UP000001811/RELEASE.metalink" "$H/RELEASE.metalink" \
    && get "$BASE/Eukaryota/UP000001811/$F" "$H/$F" \
    && get "$BASE/README" "$H/README" && { ok=1; break; }
  echo "rabbit fetch attempt $attempt failed" >&2
done
if [ "$ok" -eq 1 ]; then
  want=$(python3 -I - "$H/RELEASE.metalink" "$F" <<'PY'
import sys, re
t = open(sys.argv[1]).read()
m = re.search(r'<file name="' + re.escape(sys.argv[2]) + r'">.*?<hash type="md5">([0-9a-f]{32})</hash>', t, re.S)
print(m.group(1) if m else "")
PY
)
  got=$(md5sum "$H/$F" | cut -d' ' -f1)
  echo "md5 metalink=$want file=$got"
  if [ -n "$want" ] && [ "$want" = "$got" ]; then
    zcat "$H/$F" | awk '/^>/ {split($1, a, "|"); print ">CONT_RABIT_" a[2]; next} {print}' > "$H/rabbit_host.fasta"
    echo "records: $(grep -c '^>' "$H/rabbit_host.fasta")"
    grep -E "^Release|UP000001811" "$H/README" > "$H/release.txt"
    sha256sum "$H/$F" "$H/rabbit_host.fasta" > "$H/SHA256SUMS"
    echo "md5 $got (= RELEASE.metalink)" >> "$H/release.txt"
  else
    echo "rabbit proteome md5 mismatch or missing in metalink" >&2; rc_all=1
  fi
else
  rc_all=1
fi

echo "== 3. gippa makedb A and B"
mk() {  # D tag
  local D="$1" tag="$2" J="$1/bakta/$2/out/$2.json"
  for r in a b; do
    local out="$D/db/$tag/db_$r" s
    [ -s "$out/strain_db.faa" ] && { echo "$tag $r: exists"; continue; }
    s=$(date +%s)
    if [ "$r" = a ]; then "$GIPPA" makedb --annotation "$J" --outdir "$out" >"$D/db/$tag/makedb_$r.log" 2>&1
    else "$GIPPA" makedb --annotation "$J" --six-frame-orfs --min-orf-aa 20 --outdir "$out" >"$D/db/$tag/makedb_$r.log" 2>&1; fi
    echo "$tag $r: exit $? wall_s $(( $(date +%s) - s )) records $(grep -c '^>' "$out/strain_db.faa") bytes $(stat -c %s "$out/strain_db.faa")"
  done
}
for tag in tpa_dal1 tpa_grady tpa_haitib tpa_madras tpa_mexicoa tpa_phil1 tpa_ss14; do
  mkdir -p "$C/PXD065022/db/$tag"; mk "$C/PXD065022" "$tag" || rc_all=1; done
mkdir -p "$C/PXD055395/db/bthe_vpi5482"; mk "$C/PXD055395" bthe_vpi5482 || rc_all=1
mkdir -p "$C/PXD069004/db/ccre_na1000"; mk "$C/PXD069004" ccre_na1000 || rc_all=1
exit $rc_all
