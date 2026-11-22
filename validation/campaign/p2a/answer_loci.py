"""Answer-key loci: RefSeq (or GenBank) coordinates of named genes and the Bakta CDS/sORF that
shares their stop codon (same contig, strand and 3' end). Reads annotations only, never spectra.

    python -I answer_loci.py BAKTA_JSON SPEC.tsv OUT.tsv

SPEC.tsv columns: key, label, contig, start, end, strand (1-based closed; from the RefSeq GFF or
the GenBank feature table, written by hand into the spec). OUT.tsv adds bakta_locus,
bakta_start, bakta_end, bakta_gene, match (same_stop | overlap_same_strand | none).
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path


def main(bakta: Path, spec: Path, out: Path) -> None:
    d = json.loads(bakta.read_text())
    feats = [f for f in d["features"] if f["type"] in ("cds", "sorf")]
    rows = list(csv.DictReader(spec.open(), delimiter="\t"))
    with out.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow([*rows[0].keys(), "bakta_locus", "bakta_start", "bakta_end", "bakta_gene", "match"])
        for r in rows:
            a, b, s = int(r["start"]), int(r["end"]), r["strand"]
            stop = b if s == "+" else a
            same = [f for f in feats if f["sequence"] == r["contig"] and f["strand"] == s
                    and (f["stop"] if s == "+" else f["start"]) == stop]
            ov = [f for f in feats if f["sequence"] == r["contig"] and f["strand"] == s
                  and f["start"] <= b and f["stop"] >= a]
            f, m = (same[0], "same_stop") if same else ((ov[0], "overlap_same_strand") if ov else (None, "none"))
            w.writerow([*r.values(), *(["", "", "", ""] if f is None else
                       [f["locus"], f["start"], f["stop"], f.get("gene") or ""]), m])


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]))
