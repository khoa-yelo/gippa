"""Phase 1b truth set: amino-acid differences between NCM3722 and MG1655 Bakta proteins.

Method (pre-registered in validation/campaign/ecoli/ANSWER_KEY_1b.md, "Truth set"):
blastp both ways (-seg no -evalue 1e-10 -max_hsps 1), reciprocal best hits by bitscore with
identity >= 90 % and HSP coverage >= 80 % of both proteins; every aligned non-gap column with
different residues is a site; I<->L columns, gap columns, residues outside the HSP, unpaired
proteins are counted only. Sites get MG1655 and NCM3722 codon coordinates (1-based lowest base).

    python -I truth_set.py proteins NCM_JSON MG_JSON OUTDIR        # writes ncm.faa, mg.faa, genes.tsv
    python -I truth_set.py sites OUTDIR NCM_JSON MG_JSON           # reads n2m.tsv, m2n.tsv (blast)

Writes under OUTDIR (never in git): truth_sites.tsv, orthologs.tsv, counts.json.
"""

from __future__ import annotations

import csv
import json
import sys
from collections import Counter
from pathlib import Path

COLS = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "qlen",
    "slen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "qseq",
    "sseq",
]


def proteins(ncm_json: Path, mg_json: Path, out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    with (out / "genes.tsv").open("w", newline="") as g:
        w = csv.writer(g, delimiter="\t", lineterminator="\n")
        w.writerow(["genome", "locus", "contig", "start", "stop", "strand", "length_bp", "aa_len"])
        for tag, path in (("ncm", ncm_json), ("mg", mg_json)):
            d = json.loads(path.read_text())
            lens = {s["id"]: int(s["length"]) for s in d["sequences"]}
            with (out / f"{tag}.faa").open("w") as fa:
                for f in d["features"]:
                    if f["type"] not in ("cds", "sorf") or not f.get("aa"):
                        continue
                    fa.write(f">{f['locus']}\n{f['aa']}\n")
                    w.writerow(
                        [
                            tag,
                            f["locus"],
                            f["sequence"],
                            f["start"],
                            f["stop"],
                            f["strand"],
                            lens[f["sequence"]],
                            len(f["aa"]),
                        ]
                    )


def read_blast(path: Path) -> dict[str, dict[str, str]]:
    best: dict[str, dict[str, str]] = {}
    with path.open() as fh:
        for line in fh:
            r = dict(zip(COLS, line.rstrip("\n").split("\t"), strict=False))
            q = r["qseqid"]
            if q not in best or float(r["bitscore"]) > float(best[q]["bitscore"]):
                best[q] = r
    return best


def codon(g: dict[str, str], pos: int) -> int:
    """1-based lowest base of residue ``pos`` (1-based) of gene ``g`` (circular wrap)."""
    start, stop, n = int(g["start"]), int(g["stop"]), int(g["length_bp"])
    c = start + 3 * (pos - 1) if g["strand"] == "+" else stop - 3 * pos + 1
    return (c - 1) % n + 1


CODE = dict(
    zip(
        (a + b + c for a in "TCAG" for b in "TCAG" for c in "TCAG"),
        "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        strict=False,
    )
)
COMP = str.maketrans("ACGT", "TGCA")


def residue(seqs: dict[str, str], contig: str, start: int, strand: str) -> str:
    s = seqs[contig]
    c = (s + s[:2])[start - 1 : start + 2].upper()
    return CODE.get(c if strand == "+" else c.translate(COMP)[::-1], "X")


def sites(out: Path, ncm_json: Path, mg_json: Path) -> None:
    seqs = {
        s["id"]: s["nt"]
        for p in (ncm_json, mg_json)
        for s in json.loads(p.read_text())["sequences"]
    }
    genes = {
        (r["genome"], r["locus"]): r
        for r in csv.DictReader((out / "genes.tsv").open(), delimiter="\t")
    }
    n2m, m2n = read_blast(out / "n2m.tsv"), read_blast(out / "m2n.tsv")
    cnt: Counter[str] = Counter()
    cnt["ncm_proteins"] = sum(1 for k in genes if k[0] == "ncm")
    cnt["mg_proteins"] = sum(1 for k in genes if k[0] == "mg")
    rows, orth = [], []
    for q, r in n2m.items():
        s = r["sseqid"]
        if m2n.get(s, {}).get("sseqid") != q:
            cnt["ncm_no_rbh"] += 1
            continue
        qlen, slen, _ln = int(r["qlen"]), int(r["slen"]), int(r["length"])
        qcov = (int(r["qend"]) - int(r["qstart"]) + 1) / qlen
        scov = (int(r["send"]) - int(r["sstart"]) + 1) / slen
        if float(r["pident"]) < 90 or qcov < 0.8 or scov < 0.8:
            cnt["rbh_below_limits"] += 1
            continue
        gn, gm = genes[("ncm", q)], genes[("mg", s)]
        cnt["pairs"] += 1
        cnt[f"pairs_ncm_{gn['contig']}"] += 1
        cnt["outside_hsp_residues"] += (qlen - (int(r["qend"]) - int(r["qstart"]) + 1)) + (
            slen - (int(r["send"]) - int(r["sstart"]) + 1)
        )
        qi, si = int(r["qstart"]) - 1, int(r["sstart"]) - 1
        nd = 0
        for a, b in zip(r["qseq"], r["sseq"], strict=False):
            if a != "-":
                qi += 1
            if b != "-":
                si += 1
            if a == "-" or b == "-":
                cnt["gap_columns"] += 1
                continue
            if a == b:
                continue
            if {a, b} == {"I", "L"}:
                cnt["il_columns"] += 1
                continue
            nd += 1
            rows.append(
                {
                    "mg_locus": s,
                    "mg_contig": gm["contig"],
                    "mg_strand": gm["strand"],
                    "mg_aa_pos": si,
                    "mg_codon_start": codon(gm, si),
                    "ref_aa": b,
                    "alt_aa": a,
                    "ncm_locus": q,
                    "ncm_contig": gn["contig"],
                    "ncm_strand": gn["strand"],
                    "ncm_aa_pos": qi,
                    "ncm_codon_start": codon(gn, qi),
                }
            )
        cnt["pairs_identical"] += int(nd == 0 and float(r["pident"]) == 100.0 and qlen == slen)
        orth.append({"ncm_locus": q, "mg_locus": s, "pident": r["pident"], "n_sub": nd})
    kept = []
    for x in rows:  # each site's codons translate to its REF (MG1655) and ALT (NCM3722) residues
        ok_m = residue(seqs, x["mg_contig"], x["mg_codon_start"], x["mg_strand"]) == x["ref_aa"]
        ok_n = residue(seqs, x["ncm_contig"], x["ncm_codon_start"], x["ncm_strand"]) == x["alt_aa"]
        if ok_m and not ok_n and x["ncm_aa_pos"] == 1 and x["alt_aa"] == "M":
            # Bakta's initiator M over a non-ATG start: a start-call difference, same codon
            cnt["initiator_start_differences"] += 1
            continue
        cnt["codon_check_failed"] += int(not (ok_m and ok_n))
        kept.append(x)
    rows = kept
    cnt["sites"] = len(rows)
    cnt["proteins_with_sites"] = len({x["mg_locus"] for x in rows})
    for x in rows:
        cnt[f"sites_ncm_{x['ncm_contig']}"] += 1
    cnt["ncm_no_hit"] = cnt["ncm_proteins"] - len(n2m)
    for name, data in (("truth_sites.tsv", rows), ("orthologs.tsv", orth)):
        with (out / name).open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(data[0]), delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(data)
    (out / "counts.json").write_text(json.dumps(dict(sorted(cnt.items())), indent=2) + "\n")
    print(json.dumps(dict(sorted(cnt.items()))))


if __name__ == "__main__":
    if sys.argv[1] == "proteins":
        proteins(Path(sys.argv[2]), Path(sys.argv[3]), Path(sys.argv[4]))
    else:
        sites(Path(sys.argv[2]), Path(sys.argv[3]), Path(sys.argv[4]))
