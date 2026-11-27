"""Metadata-only file inventory of a Panorama Public folder (campaign, ADR-0077 Amendment 1).

PXD066241 is hosted at Panorama Public (panoramaweb.org), which ProteomeXchange lists without a
file list. This walks the folder's WebDAV tree (``<folder>/@files/?method=JSON``, LabKey's JSON
listing) and writes ``<out-dir>/<ACCESSION>/files.tsv`` with the columns of ``pxd_inventory.py``.
Stdlib only; run as ``python3 -I``. No raw file is downloaded. Panorama publishes no checksums,
so ``checksum`` is NA and integrity rests on TLS plus the byte size; the fetch job records sha256.

    panorama_inventory.py PXD066241 --folder '/Panorama Public/2025/<lab> - <title>' \
        --record-dir $GIPPA_ROOT/campaign/PXD066241/record --out-dir validation/campaign

R6 / ADR-0077 Amendment 1: https only, TLS verified, panoramaweb.org only; one retry per request.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

BASE = "https://panoramaweb.org"
USER_AGENT = "gippa-campaign-inventory/1"
NA = "NA"
FILE_COLUMNS = (
    "file_name",
    "file_type",
    "raw_format",
    "bytes",
    "checksum",
    "checksum_source",
    "url",
    "acquisition",
    "acquisition_source",
    "sample",
    "condition",
)


def get_json(href: str) -> tuple[dict, bytes]:
    url = BASE + href + "?method=JSON"
    assert url.startswith(BASE + "/_webdav/"), url
    last: Exception | None = None
    for attempt in (1, 2):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=120) as r:
                body = r.read()
            return json.loads(body), body
        except Exception as e:
            last = e
            if attempt == 1:
                time.sleep(20)
    raise RuntimeError(f"{url}: {last}")


def walk(href: str, records: list[bytes]):
    data, body = get_json(href)
    records.append(body)
    for f in data["files"]:
        if f.get("collection"):
            yield from walk(f["href"], records)
        else:
            yield f


def classify(rel: str) -> tuple[str, str, str]:
    """(file_type, raw_format, acquisition) from the path; the name is the only source."""
    name = rel.rsplit("/", 1)[-1]
    low = name.lower()
    if low.endswith(".raw"):
        if re.search(r"_dia_\d+", low):
            acq = "DIA"
        elif "_prm" in low:
            acq = "PRM"
        elif rel.startswith("RawFiles/Library runs/"):
            acq = "DDA"  # spectral-library runs searched with pepXML next to them
        else:
            acq = NA
        return "raw", "thermo_raw", acq
    if low.endswith(
        (".sky", ".sky.zip", ".skyd", ".skyl", ".sky.view", ".blib", ".protdb", ".msp")
    ):
        return "skyline", NA, NA
    if low.endswith((".pepxml", ".pep.xml", ".mzml", ".parquet")):
        return "result", NA, NA
    return "other", NA, NA


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("accession")
    ap.add_argument("--folder", required=True)
    ap.add_argument("--record-dir", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    a = ap.parse_args(argv)
    root = "/_webdav" + urllib.parse.quote(a.folder) + "/%40files/"
    records: list[bytes] = []
    files = list(walk(root, records))
    a.record_dir.mkdir(parents=True, exist_ok=True)
    blob = b"\n".join(records)
    (a.record_dir / "panorama_listing.jsonl").write_bytes(blob)
    utc = dt.datetime.now(dt.UTC).strftime("%Y-%m-%dT%H:%M:%SZ")
    (a.record_dir / "panorama_listing.sha256").write_text(
        f"{hashlib.sha256(blob).hexdigest()}  panorama_listing.jsonl  retrieved={utc}\n"
    )
    rows = []
    prefix = urllib.parse.unquote(root)
    for f in files:
        path = urllib.parse.unquote(f["href"])
        rel = path[len(prefix) :]
        ftype, fmt, acq = classify(rel)
        rows.append(
            {
                "file_name": rel,
                "file_type": ftype,
                "raw_format": fmt,
                "bytes": str(f.get("size", NA)),
                "checksum": NA,
                "checksum_source": NA,
                "url": BASE + f["href"],
                "acquisition": acq,
                "acquisition_source": "file_name" if acq != NA else NA,
                "sample": rel.rsplit("/", 1)[-1].rsplit(".", 1)[0] if ftype == "raw" else NA,
                "condition": NA,
            }
        )
    rows.sort(key=lambda r: r["file_name"])
    out = a.out_dir / a.accession / "files.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, FILE_COLUMNS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"{a.accession}: {len(rows)} files -> {out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
