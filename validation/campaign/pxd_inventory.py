"""Metadata-only file inventory of one ProteomeXchange dataset (campaign Phase 0, ADR-0077).

Generalises ``validation/pilot/pxd_fetch.py`` / ``px_record.py`` (hard-wired to PXD078459) to
any campaign accession. Stdlib only; run as ``python3 -I``. No raw file is downloaded: only the
ProteomeXchange XML, the PRIDE v3 project and file-list JSON, the deposit's own ``checksum.txt``
and SDRF (small text files), or, for iProX, an HTTPS ``HEAD`` per file for its size.

    pxd_inventory.py ACCESSION --record-dir DIR --out-dir validation/campaign \
        [--accessions validation/campaign/accessions.tsv]

Writes ``<out-dir>/<ACCESSION>/files.tsv`` and replaces the accession's row of
``<out-dir>/datasets.tsv``. The fetched records, with their sha256 and UTC retrieval time, stay
under ``--record-dir`` (``$GIPPA_ROOT/campaign/<ACCESSION>/record``; R11).

R6 / ADR-0077: https only, TLS verified (the default ``ssl`` context), official hosts only
(proteomecentral.proteomexchange.org, www.ebi.ac.uk/pride, ftp.pride.ebi.ac.uk over https,
download.iprox.cn). Every fetch is tried twice (one retry); a second failure is recorded as a
block in the dataset row (status ``blocked``) and the script moves on.

Exit codes: 0 inventory written (status ok or needs-decision); 1 blocked (row still written).
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import io
import json
import re
import ssl
import sys
import time
import urllib.error
import urllib.request
import xml.etree.ElementTree as ET
from collections.abc import Callable, Sequence
from dataclasses import dataclass, field
from pathlib import Path, PurePosixPath, PureWindowsPath

PX_URL = "https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID={acc}&outputMode=XML"
PRIDE_API = "https://www.ebi.ac.uk/pride/ws/archive/v3/projects/{acc}"
ALLOWED_HOSTS = (
    "proteomecentral.proteomexchange.org",
    "www.ebi.ac.uk",
    "ftp.pride.ebi.ac.uk",
    "download.iprox.cn",
    "www.iprox.cn",
)
USER_AGENT = "gippa-campaign-inventory/1"
RETRY_DELAY_S = 20
PAGE_SIZE = 100
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
DATASET_COLUMNS = (
    "accession",
    "phase",
    "organism_strain",
    "repository",
    "public",
    "n_raw_files",
    "total_raw_gb",
    "raw_format",
    "instrument",
    "acquisition",
    "sdrf_present",
    "status",
    "note",
)

Getter = Callable[[str], bytes]
Header = Callable[[str], tuple[int, int | None]]


class Blocked(Exception):
    """A source stayed unreachable after one retry."""


# ------------------------------------------------------------------ network (https only)


def _check_url(url: str) -> None:
    m = re.match(r"^https://([^/]+)/", url)
    if not m or m.group(1) not in ALLOWED_HOSTS:
        raise ValueError(f"refused (R6 / ADR-0077: https on an allowed host only): {url}")


def _open(url: str, method: str = "GET"):
    _check_url(url)
    req = urllib.request.Request(url, method=method, headers={"User-Agent": USER_AGENT})
    return urllib.request.urlopen(req, timeout=120, context=ssl.create_default_context())


def _twice(what: str, fn: Callable[[], object]) -> object:
    err: Exception | None = None
    for attempt in (1, 2):
        try:
            return fn()
        except urllib.error.HTTPError as e:
            if e.code in (401, 403, 404):  # not public or not there: a retry does not help
                raise Blocked(f"{what}: HTTP {e.code}") from e
            err = e
        except (urllib.error.URLError, TimeoutError, OSError) as e:
            err = e
        if attempt == 1:
            time.sleep(RETRY_DELAY_S)
    raise Blocked(f"{what}: {err} (after one retry)")


def http_get(url: str) -> bytes:
    def go() -> bytes:
        with _open(url) as r:
            return bytes(r.read())

    return _twice(url, go)  # type: ignore[return-value]


def http_head(url: str) -> tuple[int, int | None]:
    def go() -> tuple[int, int | None]:
        with _open(url, "HEAD") as r:
            n = r.headers.get("Content-Length")
            return r.status, int(n) if n and n.isdigit() else None

    return _twice(url, go)  # type: ignore[return-value]


def https_of(uri: str) -> str:
    """The https form of a record URI on the same official host (iProX: the .cn mirror the
    pilot used)."""
    u = re.sub(r"^(ftp|http)://", "https://", uri)
    return u.replace("https://download.iprox.org/", "https://download.iprox.cn/")


# ------------------------------------------------------------------ records


def _local(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


@dataclass
class PxRecord:
    accession: str
    repository: str
    announced: str
    title: str
    species: list[str]
    instruments: list[str]
    files: list[tuple[str, str, str]]  # (name, uri, cv name of the uri)
    links: list[str]
    text: str  # description, keywords: for the dataset-level DIA hint


def parse_px(xml: bytes) -> PxRecord:
    root = ET.fromstring(xml)
    summ = next((e for e in root.iter() if _local(e.tag) == "DatasetSummary"), None)
    repo = summ.get("hostingRepository", NA) if summ is not None else NA
    species, instruments, files, links, text = [], [], [], [], []
    for e in root.iter():
        tag = _local(e.tag)
        if tag == "Species":
            for p in e:
                if p.get("name") == "taxonomy: scientific name":
                    species.append(p.get("value", ""))
        elif tag == "Instrument":
            instruments += [p.get("name", "") for p in e if _local(p.tag) == "cvParam"]
        elif tag == "DatasetFile":
            for p in e:
                if (p.get("name") or "").lower().endswith("uri"):
                    files.append((e.get("name", ""), p.get("value", ""), p.get("name", "")))
                    break
        elif tag == "FullDatasetLink":
            links += [f"{p.get('name')}={p.get('value')}" for p in e]
        elif tag in ("Description", "Keyword", "keywordList"):
            text.append(e.text or "")
        elif tag == "cvParam" and (p_v := e.get("value")) and "keyword" in (e.get("name") or ""):
            text.append(p_v)
    # a file may be listed twice (the iProX XML repeats the list): keep the first
    seen, uniq = set(), []
    for f in files:
        if f[0] not in seen:
            seen.add(f[0])
            uniq.append(f)
    return PxRecord(
        root.get("id", NA),
        repo,
        summ.get("announceDate", NA) if summ is not None else NA,
        summ.get("title", "") if summ is not None else "",
        species,
        instruments,
        uniq,
        links,
        " ".join(text),
    )


# ------------------------------------------------------------------ classification

RAW_SUFFIX = (
    (".d.zip", "bruker_d"),
    (".d.tar.gz", "bruker_d"),
    (".d.tar", "bruker_d"),
    (".d", "bruker_d"),
    (".raw", "thermo_raw"),
    (".wiff2", "sciex_wiff"),
    (".wiff", "sciex_wiff"),
    (".wiff.scan", "sciex_wiff"),
    (".mzml", "mzml"),
    (".mzml.gz", "mzml"),
    (".mgf", "mgf"),
)
_ACQ = (
    (re.compile(r"(?i)(?<![a-z])(dia[-_]?pasef|dia|swath)(?![a-z])"), "DIA"),
    (re.compile(r"(?<=[a-z])DIA(?![a-z])"), "DIA"),  # e.g. "TimsDIA"
    (re.compile(r"(?i)(?<![a-z])(dda[-_]?pasef|dda|pasef_dda)(?![a-z])"), "DDA"),
    (re.compile(r"(?i)(?<![a-z])prm(?![a-z])"), "PRM"),
)


def raw_format(name: str, category: str, bruker: bool) -> str:
    low = name.lower()
    for suf, fmt in RAW_SUFFIX:
        if low.endswith(suf):
            return fmt
    if category == "RAW" and low.endswith((".zip", ".tar", ".tar.gz", ".7z", ".rar")):
        return "bruker_d_archive?" if bruker else "archive?"
    return NA


def file_type(name: str, category: str, fmt: str) -> str:
    low = name.lower()
    if "sdrf" in low:
        return "sdrf"
    if low in ("checksum.txt", "md5.txt", "md5sum.txt") or low.endswith(".md5"):
        return "checksum"
    if category == "RAW" or fmt in ("bruker_d", "thermo_raw", "sciex_wiff"):
        return "raw"
    if fmt in ("mzml", "mgf"):
        return "peak"
    return {"SEARCH": "search", "RESULT": "result", "PEAK": "peak", "FASTA": "fasta"}.get(
        category, "other"
    )


def acquisition_from_name(name: str) -> str:
    stem = PurePosixPath(name).name
    for pat, mode in _ACQ:
        if pat.search(stem):
            return mode
    return NA


# ------------------------------------------------------------------ SDRF and checksum.txt


def xlsx_rows(data: bytes) -> list[list[str]]:
    """The first sheet of an .xlsx file (stdlib zip + XML; no formulas, no styles)."""
    import zipfile

    ns = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
    z = zipfile.ZipFile(io.BytesIO(data))
    shared = []
    if "xl/sharedStrings.xml" in z.namelist():
        for si in ET.fromstring(z.read("xl/sharedStrings.xml")).iter(f"{ns}si"):
            shared.append("".join(t.text or "" for t in si.iter(f"{ns}t")))
    sheet = sorted(n for n in z.namelist() if n.startswith("xl/worksheets/sheet"))[0]
    rows = []
    for row in ET.fromstring(z.read(sheet)).iter(f"{ns}row"):
        vals: dict[int, str] = {}
        for c in row.iter(f"{ns}c"):
            ref = re.match(r"[A-Z]+", c.get("r", "A"))
            col = 0
            for ch in ref.group(0) if ref else "A":
                col = col * 26 + ord(ch) - 64
            v, t = c.find(f"{ns}v"), c.get("t")
            if t == "inlineStr":
                text = "".join(x.text or "" for x in c.iter(f"{ns}t"))
            elif v is None:
                text = ""
            else:
                text = shared[int(v.text or 0)] if t == "s" else v.text or ""
            vals[col - 1] = text
        rows.append([vals.get(i, "") for i in range(max(vals) + 1)] if vals else [])
    return rows


def parse_sdrf(data: bytes) -> list[dict[str, str]]:
    if data[:2] == b"PK":  # an SDRF deposited as .xlsx
        table = [r for r in xlsx_rows(data) if any(r)]
        if not table:
            return []
        head = [h.strip().lower() for h in table[0]]
        return [{h: (r[i].strip() if i < len(r) else "") for i, h in enumerate(head) if h}
                for r in table[1:]]  # fmt: skip
    text = data.decode("utf-8-sig", errors="replace").replace("\r\n", "\n").replace("\r", "\n")
    first = text.splitlines()[0] if text else ""
    delim = max(("\t", ";", ","), key=first.count)
    rows = list(csv.DictReader(io.StringIO(text), delimiter=delim))
    return [{k.strip().lower(): (v or "").strip() for k, v in r.items() if isinstance(k, str)}
            for r in rows]  # fmt: skip


SDRF_FILE_KEYS = (
    "comment[data file]",
    "data file",
    "datafile",
    "file name",
    "raw file",
    "raw file name",
)


def sdrf_index(rows: Sequence[dict[str, str]]) -> dict[str, dict[str, str]]:
    """file basename -> its SDRF row."""
    out: dict[str, dict[str, str]] = {}
    for r in rows:
        key = next((k for k in r if k in SDRF_FILE_KEYS), None)
        if key and r[key]:
            out[run_stem(PurePosixPath(r[key].replace("\\", "/")).name)] = r
    return out


def run_stem(name: str) -> str:
    """A file name without its raw/archive suffixes: SDRFs name runs with or without them."""
    low = name.lower()
    for suf in (
        ".d.tar.gz",
        ".d.zip",
        ".d.tar",
        ".tar.gz",
        ".tar",
        ".zip",
        ".7z",
        ".d",
        ".raw",
        ".wiff.scan",
        ".wiff2",
        ".wiff",
        ".mzml.gz",
        ".mzml",
    ):
        if low.endswith(suf):
            return name[: -len(suf)]
    return name


def _varying(rows: Sequence[dict[str, str]], keys: Sequence[str]) -> list[str]:
    return [k for k in keys if len({r.get(k, "") for r in rows}) > 1]


def sdrf_fields(row: dict[str, str], all_rows: Sequence[dict[str, str]]) -> tuple[str, str, str]:
    """(acquisition, sample, condition) of one SDRF row."""
    acq = NA
    for k, v in row.items():
        if "acquisition" in k:
            vl = v.lower()
            acq = (
                "DIA" if "dia" in vl or "independent" in vl or "swath" in vl
                else "DDA" if "dda" in vl or "dependent" in vl
                else "PRM" if "prm" in vl or "parallel reaction" in vl
                else v or NA
            )  # fmt: skip
            break
    if acq == NA:  # a non-standard SDRF: a DIA/DDA token in any value (e.g. "DIA, Spike-in")
        for k, v in row.items():
            if k not in ("comment[data file]", "data file", "file name") and v:
                for pat, mode in _ACQ:
                    if pat.search(v):
                        acq = mode
                        break
            if acq != NA:
                break
    sample = row.get("source name") or NA
    keys = list(row)
    factors = _varying(all_rows, [k for k in keys if k.startswith("factor value[")])
    chars = _varying(all_rows, [k for k in keys if k.startswith("characteristics[")])
    plain = _varying(
        all_rows,
        [k for k in keys if not k.startswith(("comment[", "factor", "characteristics"))
         and k not in ("source name", "assay name", "data file", "acquisition", "file name",
                       "associated search(es)")],
    )  # fmt: skip
    use = factors or chars or plain
    cond = "; ".join(f"{k}={row[k]}" for k in use if row.get(k)) or NA
    return acq, sample, cond


def parse_checksums(data: bytes) -> dict[str, str]:
    """basename -> ``<algorithm>:<hex>`` from a deposit ``checksum.txt``."""
    out = {}
    for line in data.decode("utf-8", errors="replace").splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        parts = re.split(r"\t|\s{2,}", line.strip())
        if len(parts) < 2:
            parts = line.strip().rsplit(" ", 1)
        a, b = parts[0].strip(), parts[-1].strip()
        path, hexd = (a, b) if re.fullmatch(r"[0-9a-fA-F]{32,64}", b) else (b, a)
        if not re.fullmatch(r"[0-9a-fA-F]{32,64}", hexd):
            continue
        alg = {32: "md5", 40: "sha1", 64: "sha256"}.get(len(hexd))
        if alg:
            out[PureWindowsPath(path.lstrip("*")).name] = f"{alg}:{hexd.lower()}"
    return out


# ------------------------------------------------------------------ inventory


@dataclass
class Inventory:
    accession: str
    repository: str = NA
    public: str = NA
    instruments: list[str] = field(default_factory=list)
    organism: list[str] = field(default_factory=list)
    files: list[dict[str, str]] = field(default_factory=list)
    sdrf_names: list[str] = field(default_factory=list)
    blocks: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)
    text: str = ""


class Recorder:
    """Saves every fetched record with its sha256 and UTC retrieval time."""

    def __init__(self, record_dir: Path, get: Getter) -> None:
        self.dir, self.get = record_dir, get
        record_dir.mkdir(parents=True, exist_ok=True)
        self.info: list[str] = []

    def fetch(self, url: str, name: str) -> bytes:
        data = self.get(url)
        (self.dir / name).write_bytes(data)
        sha = hashlib.sha256(data).hexdigest()
        now = dt.datetime.now(dt.UTC).strftime("%Y-%m-%dT%H:%M:%SZ")
        self.info.append(f"{name}\t{url}\t{now}\t{sha}")
        (self.dir / "records.info").write_text("name\turl\tretrieved_utc\tsha256\n"
                                               + "\n".join(self.info) + "\n")  # fmt: skip
        return data


def _pride(acc: str, rec: Recorder, inv: Inventory) -> None:
    proj = json.loads(rec.fetch(PRIDE_API.format(acc=acc), "pride_project.json"))
    inv.instruments = [i.get("name", "") for i in proj.get("instruments", [])] or inv.instruments
    inv.organism = [o.get("name", "") for o in proj.get("organisms", [])] or inv.organism
    inv.text += " " + " ".join(
        str(proj.get(k, "")) for k in ("sampleProcessingProtocol", "dataProcessingProtocol",
                                       "keywords", "experimentTypes", "projectTags")
    )  # fmt: skip
    listing: list[dict] = []
    page = 0
    while True:
        url = f"{PRIDE_API.format(acc=acc)}/files?pageSize={PAGE_SIZE}&page={page}"
        chunk = json.loads(rec.fetch(url, f"pride_files.{page}.json"))
        if not chunk:
            break
        listing += chunk
        if len(chunk) < PAGE_SIZE:
            break
        page += 1
    bruker = any("tims" in i.lower() or "bruker" in i.lower() for i in inv.instruments)
    for f in listing:
        name = f.get("fileName", "")
        cat = (f.get("fileCategory") or {}).get("value", "")
        ftp = next((loc["value"] for loc in f.get("publicFileLocations", [])
                    if loc.get("value", "").startswith("ftp://")), "")  # fmt: skip
        fmt = raw_format(name, cat, bruker)
        ck = (f.get("checksum") or "").strip().lower()
        alg = {32: "md5", 40: "sha1", 64: "sha256"}.get(len(ck))
        inv.files.append({
            "file_name": name,
            "file_type": file_type(name, cat, fmt),
            "raw_format": fmt,
            "bytes": str(f.get("fileSizeBytes", NA)),
            "checksum": f"{alg}:{ck}" if alg else NA,
            "checksum_source": "pride_api" if alg else NA,
            "url": https_of(ftp) if ftp else NA,
        })  # fmt: skip


def _iprox(px: PxRecord, inv: Inventory, head: Header) -> None:
    for name, uri, cv in px.files:
        cat = "RAW" if "raw file" in cv.lower() else "SEARCH" if "search" in cv.lower() else ""
        fmt = raw_format(name, cat, False)
        url = https_of(uri)
        size: int | None = None
        try:
            _, size = head(url)
        except (Blocked, ValueError) as e:
            inv.blocks.append(f"HEAD {name}: {e}")
        inv.files.append({
            "file_name": name, "file_type": file_type(name, cat, fmt), "raw_format": fmt,
            "bytes": str(size) if size is not None else NA, "checksum": NA,
            "checksum_source": NA, "url": url,
        })  # fmt: skip
    inv.notes.append("iProX: no record checksums (sha256 recorded at download, as in the pilot)")


def inventory(acc: str, record_dir: Path, get: Getter = http_get,
              head: Header = http_head) -> Inventory:  # fmt: skip
    inv = Inventory(acc)
    rec = Recorder(record_dir, get)
    try:
        px = parse_px(rec.fetch(PX_URL.format(acc=acc), f"{acc}.xml"))
    except Blocked as e:
        inv.blocks.append(f"ProteomeXchange record: {e}")
        return inv
    inv.repository, inv.instruments, inv.organism = px.repository, px.instruments, px.species
    inv.text = px.title + " " + px.text
    inv.public = "yes" if px.announced != NA else "no"
    try:
        if px.repository == "PRIDE":
            _pride(acc, rec, inv)
        elif px.repository == "iProX":
            _iprox(px, inv, head)
        else:
            inv.notes.append(
                f"hosted at {px.repository} ({'; '.join(px.links)}): not an ADR-0077 source, "
                "no file list fetched"
            )
            return inv
    except Blocked as e:
        inv.blocks.append(f"{px.repository} file list: {e}")
        return inv
    _enrich(inv, rec)
    return inv


def _deposit_checksums(inv: Inventory, rec: Recorder) -> None:
    """Fill missing raw checksums from the deposit's ``checksum.txt`` (PRIDE writes one per
    submission; the v3 file list does not always include it)."""
    raws = [g for g in inv.files if g["file_type"] == "raw"]
    if not raws or all(g["checksum"] != NA for g in raws):
        return
    urls = [f["url"] for f in inv.files if f["file_type"] == "checksum" and f["url"] != NA]
    dirs = sorted({g["url"].rsplit("/", 1)[0] for g in raws if g["url"].startswith("https://")})
    urls += [f"{d}/checksum.txt" for d in dirs if f"{d}/checksum.txt" not in urls]
    by_name = {f["file_name"]: f for f in inv.files}
    for i, url in enumerate(urls):
        try:
            sums = parse_checksums(rec.fetch(url, f"checksum.{i}.txt"))
        except (Blocked, ValueError) as e:
            inv.notes.append(f"no deposit checksum.txt ({str(e).rsplit(': ', 1)[-1]})")
            continue
        for n, c in sums.items():
            g = by_name.get(n)
            if g and g["checksum"] == NA:
                g["checksum"], g["checksum_source"] = c, "deposit:checksum.txt"


def _norm_acq(v: str) -> str:
    return "top-down" if v.lower().startswith("top-down") else v


def _enrich(inv: Inventory, rec: Recorder) -> None:
    """Deposit checksum.txt for missing checksums; SDRF for acquisition, sample, condition.

    Precedence for the acquisition mode: the submitter's SDRF, then a DIA/DDA/PRM token in the
    file name, then a PRIDE community-annotated SDRF (a dataset-wide template, often one value
    for every file)."""
    _deposit_checksums(inv, rec)
    tables: dict[bool, tuple[dict[str, dict[str, str]], list[dict[str, str]]]] = {
        False: ({}, []), True: ({}, [])}  # fmt: skip
    for f in inv.files:
        if f["file_type"] == "sdrf" and f["url"] != NA and not f["file_name"].endswith(".zip"):
            try:
                rows = parse_sdrf(rec.fetch(f["url"], f"sdrf__{f['file_name']}"))
            except (Blocked, ValueError) as e:
                inv.blocks.append(f"{f['file_name']}: {e}")
                continue
            inv.sdrf_names.append(f["file_name"])
            idx, all_rows = tables["community" in f["file_name"].lower()]
            all_rows += rows
            idx.update(sdrf_index(rows))
    conflicts = 0
    for f in inv.files:
        acq, src, sample, cond = NA, NA, NA, NA
        if f["file_type"] == "raw":
            name = f["file_name"]
            token = acquisition_from_name(name)
            for community in (False, True):
                idx, all_rows = tables[community]
                row = idx.get(run_stem(name))
                if not row:
                    continue
                a, s, c = sdrf_fields(row, all_rows)
                if sample == NA:
                    sample = s
                if cond == NA:
                    cond = c
                if acq == NA and a != NA and (not community or token == NA):
                    # a community SDRF is a dataset-wide template: its mode is a hint only
                    acq = _norm_acq(a) + ("?" if community else "")
                    src = "sdrf_community" if community else "sdrf"
            if acq == NA and token != NA:
                acq, src = token, "file_name"
            if token != NA and acq not in (NA, token) and src == "sdrf":
                conflicts += 1
        f.update(acquisition=acq, acquisition_source=src, sample=sample, condition=cond)
    if conflicts:
        inv.notes.append(
            f"{conflicts} raw files: SDRF acquisition differs from the file-name token"
        )


# ------------------------------------------------------------------ outputs


def _clean(v: str) -> str:
    return re.sub(r"[\t\r\n]+", " ", v).strip() or NA


def write_tsv(path: Path, cols: Sequence[str], rows: Sequence[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, cols, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({c: _clean(str(r.get(c, NA))) for c in cols})


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def dataset_row(inv: Inventory, meta: dict[str, str]) -> dict[str, str]:
    raws = [f for f in inv.files if f["file_type"] == "raw"]
    sized = [int(f["bytes"]) for f in raws if f["bytes"].isdigit()]
    fmts = sorted({f["raw_format"] for f in raws})
    modes = sorted({f["acquisition"] for f in raws if f["acquisition"] != NA})
    unknown = sum(f["acquisition"] == NA for f in raws)
    if modes:
        acq = modes[0] if len(modes) == 1 else "mixed (" + ", ".join(modes) + ")"
        if unknown:
            acq += f"; {unknown} files undetermined"
    else:
        hint = re.search(r"(?i)dia[- ]?pasef|data[- ]independent|\bDIA\b|SWATH", inv.text)
        acq = "DIA (record text; per-file unknown)" if hint else NA
    missing_ck = sum(f["checksum"] == NA for f in raws)
    notes = list(inv.notes)
    if raws and missing_ck:
        notes.append(f"{missing_ck}/{len(raws)} raw files without a record checksum")
    if raws and len(sized) < len(raws):
        notes.append(f"{len(raws) - len(sized)} raw files without a size")
    if inv.blocks:
        status = "blocked"
        notes = [f"BLOCKED: {b}" for b in inv.blocks] + notes
    elif not raws:
        status = "needs-decision"
        notes.append("no raw file list")
    else:
        status = "ok"
    if meta.get("note"):
        notes.append(meta["note"])
    if meta.get("status") and not inv.blocks:  # a fetch block is never overridden
        status = meta["status"]
    return {
        "accession": inv.accession,
        "phase": meta.get("phase", NA),
        "organism_strain": meta.get("organism_strain") or "; ".join(inv.organism) or NA,
        "repository": inv.repository,
        "public": inv.public,
        "n_raw_files": str(len(raws)),
        "total_raw_gb": f"{sum(sized) / 1e9:.1f}" if sized else NA,
        "raw_format": ", ".join(fmts) or NA,
        "instrument": "; ".join(dict.fromkeys(inv.instruments)) or NA,
        "acquisition": acq,
        "sdrf_present": "yes: " + ", ".join(inv.sdrf_names) if inv.sdrf_names else "no",
        "status": status,
        "note": "; ".join(notes) or NA,
    }


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("accession")
    ap.add_argument("--record-dir", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument(
        "--accessions", type=Path, help="accession, phase, organism_strain, status, note"
    )
    a = ap.parse_args(argv)
    acc = a.accession.strip().upper()
    if not re.fullmatch(r"PXD\d{6}", acc):
        print(f"not a ProteomeXchange accession: {acc}", file=sys.stderr)
        return 2
    meta = (
        next((r for r in read_tsv(a.accessions) if r.get("accession") == acc), {})
        if a.accessions
        else {}
    )
    inv = inventory(acc, a.record_dir)

    def nat(f: dict[str, str]) -> list[int | str]:
        return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", f["file_name"])]

    write_tsv(a.out_dir / acc / "files.tsv", FILE_COLUMNS, sorted(inv.files, key=nat))
    row = dataset_row(inv, meta)
    table = a.out_dir / "datasets.tsv"
    rows = [r for r in read_tsv(table) if r.get("accession") != acc] + [row]
    order = [r["accession"] for r in read_tsv(a.accessions)] if a.accessions else []
    rows.sort(key=lambda r: (order.index(r["accession"]) if r["accession"] in order else 999,
                             r["accession"]))  # fmt: skip
    write_tsv(table, DATASET_COLUMNS, rows)
    print("\t".join(f"{k}={row[k]}" for k in DATASET_COLUMNS))
    return 1 if row["status"] == "blocked" else 0


if __name__ == "__main__":
    sys.exit(main())
