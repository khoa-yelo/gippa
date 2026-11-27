"""Archaeal annotation check of the handover plan (datasets 9 and 10; M. hyorhinis rule).

Compares the Bakta 1.12.2 CDS calls (CDS + sORF features of the Bakta JSON) with the RefSeq
annotation (answer key, never a GIPPA input):

- CDS count: Bakta CDS + sORF over RefSeq CDS (non-pseudo), pass when within 10 % (0.90-1.10);
- stop agreement: share of RefSeq CDS (non-pseudo, one per ID; multi-part CDS by their 3' end)
  whose stop codon (contig, strand, 3' end coordinate) is the stop of a Bakta CDS/sORF; pass when
  >= 90 %. Wilson 95 % CI reported.

Aggregates only (counts and shares) are written as JSON on stdout.

    annotation_check.py BAKTA_JSON REFSEQ_GFF_GZ
"""

from __future__ import annotations

import gzip
import json
import math
import sys
from collections import defaultdict


def wilson(k: int, n: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if n == 0:
        return (math.nan, math.nan)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def refseq_cds(path: str) -> list[tuple[str, str, int]]:
    parts: dict[str, list[tuple[str, str, int, int]]] = defaultdict(list)
    with gzip.open(path, "rt") as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "CDS":
                continue
            attr = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            if attr.get("pseudo") == "true":
                continue
            parts[attr.get("ID", f"{f[0]}:{f[3]}")].append((f[0], f[6], int(f[3]), int(f[4])))
    out = []
    for segs in parts.values():
        contig, strand = segs[0][0], segs[0][1]
        stop = max(s[3] for s in segs) if strand == "+" else min(s[2] for s in segs)
        out.append((contig, strand, stop))
    return out


def bakta_cds(path: str) -> list[tuple[str, str, int]]:
    d = json.load(open(path))  # noqa: SIM115
    out = []
    for x in d["features"]:
        if x["type"] not in ("cds", "sorf"):
            continue
        stop = x["stop"] if x["strand"] == "+" else x["start"]
        out.append((x["sequence"], x["strand"], int(stop)))
    return out


def main() -> int:
    bj, gff = sys.argv[1:3]
    ref, bak = refseq_cds(gff), bakta_cds(bj)
    bset = set(bak)
    k = sum(r in bset for r in ref)
    n = len(ref)
    ratio = len(bak) / n
    lo, hi = wilson(k, n)
    res = {
        "refseq_cds_nonpseudo": n,
        "bakta_cds_plus_sorf": len(bak),
        "count_ratio": round(ratio, 4),
        "count_within_10pct": 0.9 <= ratio <= 1.1,
        "refseq_cds_with_bakta_stop": k,
        "stop_share": round(k / n, 4),
        "stop_share_ci95": [round(lo, 4), round(hi, 4)],
        "stop_share_ge_0.90": k / n >= 0.90,
        "bakta_cds_with_refseq_stop": sum(b in set(ref) for b in bak),
        "pass": (0.9 <= ratio <= 1.1) and k / n >= 0.90,
    }
    print(json.dumps(res, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
