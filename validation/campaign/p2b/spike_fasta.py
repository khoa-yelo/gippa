"""Build the PXD065346 spike-in FASTA: the six H. volcanii small proteins of the paper.

The paper (PMC12850539) spikes HVO_0758, HVO_2212, HVO_2753, HVO_2922, HVO_2983 and HVO_A0101
into the C. difficile lysate. Their translations are taken from the RefSeq DS2 annotation
(GCF_000025685.1, the PXD066008 prep download: genome + GFF, MD5-checked against NCBI) by the
GFF ``old_locus_tag``, translated with table 11, and written as ``>SPIKE_<old_locus_tag>``.
Output lives under $GIPPA_ROOT (not committed); this script prints aggregates only.

    spike_fasta.py GENOME.fna.gz GENOMIC.gff.gz OUT.faa
"""

from __future__ import annotations

import gzip
import hashlib
import sys

from Bio import SeqIO
from Bio.Seq import Seq

WANT = ("HVO_0758", "HVO_2212", "HVO_2753", "HVO_2922", "HVO_2983", "HVO_A0101")


def main(genome: str, gff: str, out: str) -> int:
    seqs = {r.id: r.seq for r in SeqIO.parse(gzip.open(genome, "rt"), "fasta")}
    gene_old: dict[str, str] = {}
    cds: dict[str, tuple[str, int, int, str]] = {}
    with gzip.open(gff, "rt") as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            a = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            if f[2] == "gene" and a.get("old_locus_tag") in WANT:
                gene_old[a["ID"]] = a["old_locus_tag"]
            elif f[2] == "CDS" and a.get("Parent") in gene_old:
                cds[gene_old[a["Parent"]]] = (f[0], int(f[3]), int(f[4]), f[6])
    with open(out, "w") as w:
        for tag in WANT:
            if tag not in cds:
                print(f"{tag}\tmissing", file=sys.stderr)
                continue
            c, s, e, st = cds[tag]
            nt = seqs[c][s - 1 : e]
            nt = nt.reverse_complement() if st == "-" else nt
            aa = str(Seq(str(nt)).translate(table=11, to_stop=True))
            aa = "M" + aa[1:]
            w.write(f">SPIKE_{tag}\n{aa}\n")
            print(f"{tag}\t{c}:{s}-{e}{st}\t{len(aa)} aa\tsha256:{hashlib.sha256(aa.encode()).hexdigest()[:16]}")
    return 0


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:4]))
