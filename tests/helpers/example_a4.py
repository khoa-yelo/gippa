"""A synthetic stand-in for the pilot's six-frame arm A4 (route E of M12-01; spec
docs/specs/M12-01-example-dataset.md, amendment 2026-11-06: subset rules 3E, 7 and 9, AC 20-21c).

Nothing here is real data (R11): every sequence, ID and number is synthetic. The files follow the
paths ``build_example.sources()`` reads for route E (section 3):

- ``$GIPPA_ROOT/runs/pilot/sixframe-fa49b8f/A4/A4_report.parquet`` (DIA-NN 2.3.1 rows),
- ``$GIPPA_ROOT/runs/pilot/sixframe-fa49b8f/db/db_id_map.tsv`` (``gippa makedb`` columns),
- ``$GIPPA_ROOT/runs/pilot/sixframe-fa49b8f/db/strain_db.faa`` (plain FASTA, as makedb writes it).

``write_a4(root, chrom, start, end)`` writes them for a chromosome sequence ``chrom`` and the
0-based half-open region ``[start, end)``. The ORF records of ``strain_db.faa`` are the code-11
translation of ``chrom`` at each ORF's interval on its strand (some with the first residue read
as ``M``, which rule 9 tolerates); ``nterm_extension`` records are not (rule 9 does not check
them). ``plant_orfs`` first makes every ORF interval G/C-only, so every frame of either strand is
stop-free there (G/C codons never code a stop) and the records hold no ``*``.

The entries and rows are the spec's rule 9 examples, generalised to any region (``S`` = start):
with the default region ``[1965000, 2165000)`` the IDs are exactly the spec's. Each row has its
own ``Stripped.Sequence`` (``Row.peptide``), so a test finds it after the build.

``make_pilot_e(root, start, end, length)`` writes a minimal routes A-D pilot (pass 3 and the
Cascadia arm, the files ``build_example.sources()`` reads) plus A4 on one planted chromosome, for
builds with ``--start``/``--end``.
"""

from __future__ import annotations

import gzip
import random
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

CHROM = "NZ_CP058217.1"
RUNS = ("WT_1", "WT_2", "WT_3")
A4_DIR = Path("runs") / "pilot" / "sixframe-fa49b8f"
REPORT_REL = A4_DIR / "A4" / "A4_report.parquet"
MAP_REL = A4_DIR / "db" / "db_id_map.tsv"
FAA_REL = A4_DIR / "db" / "strain_db.faa"
ORF_COLUMNS = ("Protein.Ids", "Protein.Group", "Protein.Names")  # ORF IDs in these (rule 9)
PERSONAL = "/home/users/alice/raw"  # File.Name holds an absolute path (rule 8, R12)
ID_MAP_COLUMNS = ("db_id", "kind", "locus_tag", "class", "contig", "start", "end", "strand")

_AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
_CODE = {a + b + c: _AA[16 * i + 4 * j + k] for i, a in enumerate("TCAG")
         for j, b in enumerate("TCAG") for k, c in enumerate("TCAG")}  # fmt: skip


def translate(nt: str) -> str:
    return "".join(_CODE[nt[i : i + 3]] for i in range(0, len(nt) - 2, 3))


def revcomp(nt: str) -> str:
    return nt.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def region_contig(start: int, end: int) -> str:
    """The bundle genome's contig ID (``build_example.py``; spec rule 9's ``C``)."""
    return f"nissle1917_{CHROM.replace('.', '_')}_{start}_{end}"


def orf_id(cls: str, contig: str, s1: int, e1: int, strand: str) -> str:
    """``GIPPAORF_{class}_{contig}_{s1}_{e1}_{p|m}`` (1-based closed, makedb's ``OrfRecord``)."""
    return f"GIPPAORF_{cls}_{contig}_{s1}_{e1}_{'p' if strand == '+' else 'm'}"


# ---------------------------------------------------------------- entries


@dataclass(frozen=True)
class Entry:
    """One ``db_id_map.tsv`` record (1-based closed ``s1..e1`` on ``NZ_CP058217.1``)."""

    key: str
    kind: str  # cds | sorf | orf
    s1: int
    e1: int
    strand: str
    cls: str | None = None  # the ORF class (kind orf)
    locus_tag: str | None = None
    name: str | None = None  # CDS/sORF ID; ORFs take orf_id() unless set (a malformed ID)
    first_m: bool = False  # the strain_db record reads the first codon as M

    @property
    def db_id(self) -> str:
        if self.kind == "orf" and self.name is None:
            assert self.cls is not None
            return orf_id(self.cls, CHROM, self.s1, self.e1, self.strand)
        assert self.name is not None
        return self.name

    def in_region(self, start: int, end: int) -> bool:
        """Rule 3E: ``start - 1 >= S`` and ``end <= E``, either strand."""
        return self.s1 - 1 >= start and self.e1 <= end

    def rewritten(self, start: int, end: int) -> str:
        """Rule 9: ORF IDs in region coordinates; CDS and sORF IDs unchanged."""
        if self.kind != "orf":
            return self.db_id
        assert self.cls is not None
        return orf_id(self.cls, region_contig(start, end), self.s1 - start, self.e1 - start,
                      self.strand)  # fmt: skip

    def parsed(self, start: int, end: int) -> tuple[str, str, int, int, str]:
        """What ``gippa.novel.fdr.parse_orf_id`` must give for the rewritten ID (rule 9)."""
        assert self.cls is not None
        return (self.cls, region_contig(start, end), self.s1 - start - 1, self.e1 - start,
                self.strand)  # fmt: skip


def entries(start: int, end: int) -> dict[str, Entry]:
    """The rule-9 examples of the spec (``S = start``, ``E = end``), keyed by a short name.

    With the default region the ORF IDs are the spec's literally:
    ``O_AS`` = ``GIPPAORF_antisense_NZ_CP058217.1_1965101_1965400_m``, ``O_OUT`` =
    ``..._cds_outframe_NZ_CP058217.1_2000002_2000301_p``, ``O_XEND`` = ``..._intergenic_..._
    2164900_2165120_p``, ``O_FAR`` = ``..._antisense_..._10_309_m``, ``O_FIRST`` =
    ``..._intergenic_..._1965001_1965300_p``, ``O_LAST`` = ``..._intergenic_..._2164701_2165000_p``.
    """
    n = end - start
    mid = 35_001 if n > 40_000 else n // 2  # S + 35002 in the default region
    cds_x = Entry("CDS_X", "cds", start + 2_001, start + 2_600, "+", name="CDS_X",
                  locus_tag="CDS_X")  # fmt: skip
    out = [
        cds_x,
        Entry(
            "CDS_S", "sorf", start + 3_001, start + 3_090, "-", name="SORF_Y", locus_tag="SORF_Y"
        ),
        Entry("CDS_FAR", "cds", 1_001, 1_600, "+", name="CDS_FAR", locus_tag="CDS_FAR"),
        Entry("O_AS", "orf", start + 101, start + 400, "-", cls="antisense", first_m=True),
        Entry("O_OUT", "orf", start + mid + 1, start + mid + 300, "+", cls="cds_outframe"),
        Entry("O_XEND", "orf", end - 100, end + 120, "+", cls="intergenic"),
        Entry("O_XSTART", "orf", start, start + 299, "-", cls="antisense"),
        Entry("O_FAR", "orf", 10, 309, "-", cls="antisense"),
        Entry("O_FIRST", "orf", start + 1, start + 300, "+", cls="intergenic"),
        Entry("O_LAST", "orf", end - 299, end, "+", cls="intergenic", first_m=True),
        Entry(
            "O_NT",
            "orf",
            start + mid + 1_001,
            start + mid + 1_300,
            "-",
            cls="nterm_extension",
            locus_tag="CDS_X",
        ),
    ]
    return {e.key: e for e in out}


# ---------------------------------------------------------------- rows


@dataclass(frozen=True)
class Row:
    """One DIA-NN precursor (written once per run): ``ids`` are entry keys or literal IDs."""

    peptide: str
    ids: tuple[str, ...]
    decoy: int
    kept: bool  # rule 3E
    group: tuple[str, ...] | None = None  # Protein.Group (default: the same as ids)


ROWS: tuple[Row, ...] = (
    Row("CDSTARGETK", ("CDS_X",), 0, True),  # region CDS, target: IDs unchanged
    Row("ASDECOYK", ("O_AS",), 1, True),  # region ORF only, decoy, minus strand
    Row("ASTARGETK", ("O_AS",), 0, True),  # the same ORF, target
    Row("MIXTARGETK", ("CDS_X", "O_OUT"), 0, True),  # CDS + ORF: second ID rewritten, order kept
    Row("MIXDECOYK", ("O_OUT", "CDS_X"), 1, True, ("O_OUT",)),  # ORF first; a shorter group
    Row("XENDK", ("O_XEND",), 0, False),  # crosses the region end
    Row("XENDDECOYK", ("O_XEND",), 1, False),
    Row("XSTARTK", ("CDS_X", "O_XSTART"), 0, False),  # crosses the region start
    Row("FARORFK", ("CDS_X", "O_FAR"), 0, False),  # region CDS + out-of-region ORF
    Row("FARCDSK", ("CDS_FAR",), 0, False),  # out-of-region CDS
    Row("CONTORFK", ("Cont_P1", "O_AS"), 0, False),  # contaminant + region ORF
    Row("CONTDECOYK", ("Cont_P1", "O_AS"), 1, False),
    Row("FIRSTBASEK", ("O_FIRST",), 0, True),  # ORF starts at the first region base
    Row("LASTBASEK", ("O_LAST",), 0, True),  # ORF ends at the last region base
    Row("NTERMK", ("O_NT",), 0, True),  # nterm_extension: rewritten, no translation check
    Row("SORFDECOYK", ("CDS_S",), 1, True),  # region sORF, decoy (no ORF ID)
)
KEPT = frozenset(r.peptide for r in ROWS if r.kept)
REWRITTEN_KEYS = ("O_AS", "O_OUT", "O_FIRST", "O_LAST", "O_NT")  # distinct ORF entries kept


def _ids(keys: Iterable[str], ent: Mapping[str, Entry]) -> list[str]:
    return [ent[k].db_id if k in ent else k for k in keys]


def report_rows(rows: Sequence[Row], ent: Mapping[str, Entry]) -> list[dict[str, Any]]:
    """The DIA-NN columns of the route A fixture plus ``Protein.Names`` (a third ID column)."""
    out: list[dict[str, Any]] = []
    for r, run in enumerate(RUNS):
        for j, row in enumerate(rows):
            ids = ";".join(_ids(row.ids, ent))
            group = ";".join(_ids(row.group or row.ids, ent))
            out.append({
                "File.Name": f"{PERSONAL}/{run}.mzML", "Run": run, "Protein.Group": group,
                "Protein.Ids": ids, "Protein.Names": ids, "Genes": "",
                "Stripped.Sequence": row.peptide, "Modified.Sequence": row.peptide,
                "Precursor.Id": f"{row.peptide}2", "Precursor.Charge": 2,
                "Precursor.Normalised": 1000.0 * (j + 1) + 10 * r,
                "Q.Value": 0.001 * (j + 1), "Global.Q.Value": 0.0005 * (j + 1),
                "PEP": 0.002 * (j + 1), "PG.MaxLFQ": 5000.0 + 100 * j + r, "Decoy": row.decoy,
            })  # fmt: skip
    return out


# ---------------------------------------------------------------- writing


def plant_orfs(chrom: str, ents: Iterable[Entry], seed: int = 6) -> str:
    """``chrom`` with every ORF interval (clipped to the chromosome) replaced by G/C bases."""
    rng = random.Random(seed)
    nt = list(chrom)
    for e in ents:
        if e.kind == "orf":
            a, b = max(0, e.s1 - 1), min(len(nt), e.e1)
            nt[a:b] = rng.choices("GC", k=b - a)
    return "".join(nt)


def record(e: Entry, chrom: str, rng: random.Random) -> str:
    """The ``strain_db.faa`` sequence of an entry."""
    if e.kind != "orf" or e.cls == "nterm_extension":
        return "M" + "".join(rng.choices("ACDEFGHIKLMNPQRSTVWY", k=60))
    nt = chrom[e.s1 - 1 : e.e1]
    aa = translate(revcomp(nt) if e.strand == "-" else nt)
    return "M" + aa[1:] if e.first_m else aa


def _mutate(aa: str, i: int = 4) -> str:
    return aa[:i] + ("W" if aa[i] != "W" else "Y") + aa[i + 1 :]


@dataclass(frozen=True)
class A4:
    report: Path
    id_map: Path
    strain_db: Path
    entries: Mapping[str, Entry]
    rows: tuple[Row, ...]


def write_a4(
    root: Path,
    chrom: str,
    start: int,
    end: int,
    *,
    rows: Sequence[Row] = ROWS,
    ents: Mapping[str, Entry] | None = None,
    map_shift: str | None = None,
    bad_translation: str | None = None,
    map_strand: str | None = None,
    long_record: str | None = None,
    column_types: Mapping[str, str] | None = None,
    extras: Mapping[str, str] | None = None,
) -> A4:
    """The three A4 files under ``root`` for the chromosome sequence ``chrom``.

    ``map_shift=KEY`` writes that entry's map start 3 bases after its ID's (rule 9 check 1);
    ``bad_translation=KEY`` changes residue 5 of its record (rule 9 check 3); ``map_strand=KEY``
    writes the other strand in its map row; ``long_record=KEY`` appends one residue to its
    record; ``column_types={column: "dictionary" | "string_view"}`` stores report columns in that
    Arrow type; ``extras={path under the A4 dir: text}`` writes PROVENANCE's command and config
    files (``A4/command_step1.txt``, ``jobs/diann_A4.cfg``, ...). The defaults write none of
    these."""
    ent = dict(ents if ents is not None else entries(start, end))
    rng = random.Random(4)
    report, id_map, faa = root / REPORT_REL, root / MAP_REL, root / FAA_REL
    for p in (report, id_map, faa):
        p.parent.mkdir(parents=True, exist_ok=True)
    with id_map.open("w", encoding="utf-8") as fh:
        fh.write("\t".join(ID_MAP_COLUMNS) + "\n")
        for e in ent.values():
            s1 = e.s1 + 3 if e.key == map_shift else e.s1
            strand = ("-" if e.strand == "+" else "+") if e.key == map_strand else e.strand
            cells = (e.db_id, e.kind, e.locus_tag or "", e.cls or "", CHROM, s1, e.e1, strand)
            fh.write("\t".join(str(c) for c in cells) + "\n")
    text = []
    for e in ent.values():
        seq = record(e, chrom, rng)
        if e.key == bad_translation:
            seq = _mutate(seq)
        if e.key == long_record:
            seq += "W"
        text.append(f">{e.db_id}\n" + "\n".join(seq[i : i + 60] for i in range(0, len(seq), 60)))
    faa.write_text("\n".join(text) + "\n", encoding="utf-8")
    table = pa.Table.from_pandas(pd.DataFrame(report_rows(rows, ent)), preserve_index=False)
    for col, kind in (column_types or {}).items():
        i = table.schema.get_field_index(col)
        arr = table.column(i).combine_chunks()
        new = arr.dictionary_encode() if kind == "dictionary" else arr.cast(pa.string_view())
        table = table.set_column(i, pa.field(col, new.type), new)
    pq.write_table(table, report)
    for rel, text in (extras or {}).items():
        path = root / A4_DIR / rel
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8")
    return A4(report, id_map, faa, ent, tuple(rows))


def with_rows(
    rows: Sequence[Row],
    *,
    decoys: bool = True,
    targets: bool = True,
    orf_rows: bool = True,
    ents: Mapping[str, Entry] | None = None,
) -> tuple[Row, ...]:
    """``rows`` without decoys, targets or rows naming an ORF entry (rule 7 variants)."""
    orf_keys = {k for k, e in (ents or entries(10_000, 40_000)).items() if e.kind == "orf"}
    keep = [r for r in rows if (decoys or not r.decoy) and (targets or r.decoy)
            and (orf_rows or not set(r.ids) & orf_keys)]  # fmt: skip
    return tuple(keep)


def bad_class(e: Entry) -> Entry:
    """An ORF entry whose class ``parse_orf_id`` does not know (rule 9 check 2 fails)."""
    return replace(e, cls="weird_class")


# ---------------------------------------------------------------- a minimal routes A-D pilot


def _fasta_gz(path: Path, records: list[tuple[str, str]]) -> None:
    text = "".join(f">{n}\n" + "\n".join(s[i : i + 60] for i in range(0, len(s), 60)) + "\n"
                   for n, s in records)  # fmt: skip
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(gzip.compress(text.encode(), mtime=0))


def make_pilot_e(
    root: Path, start: int, end: int, length: int, *, seed: int = 1917, **a4: Any
) -> tuple[str, A4]:
    """Routes A-D pass-3 and Cascadia files (two region proteins, a target and a decoy row per
    run, one placeable Cascadia peptide) and A4 (``write_a4(**a4)``) on one chromosome of
    ``length`` bp with the ORFs planted; (chromosome, A4)."""
    rng = random.Random(seed)
    ents = a4.get("ents") or entries(start, end)
    chrom = plant_orfs("".join(rng.choices("ACGT", k=length)), ents.values())
    p3, a3 = root / "pilot" / "pass3", root / "pilot" / "A3_cascadia_cpu"
    _fasta_gz(p3 / "genome" / "nissle1917_wt.fna.gz", [(CHROM, chrom)])
    prot = [("PA_1", start + 5_001, start + 5_600, "+", "cds"),
            ("PA_2", start + 6_001, start + 6_300, "-", "sorf")]  # fmt: skip
    (p3 / "db_a2").mkdir(parents=True, exist_ok=True)
    with (p3 / "db_a2" / "db_id_map.tsv").open("w", encoding="utf-8") as fh:
        fh.write("db_id\tcontig\tstart\tend\tstrand\tkind\n")
        for pid, s, e, strand, kind in prot:
            fh.write(f"{pid}\t{CHROM}\t{s}\t{e}\t{strand}\t{kind}\n")
    _fasta_gz(p3 / "db_a2" / "strain_db.faa.gz",
              [(p[0], "M" + "".join(rng.choices("ACDEFGHIKLMNPQRSTVWY", k=80)))
               for p in prot])  # fmt: skip
    a_rows = [Row("APROTONEK", ("PA_1",), 0, True), Row("APROTTWOK", ("PA_2",), 0, True),
              Row("KENOTORPA", ("PA_1",), 1, True)]  # fmt: skip
    (p3 / "A2").mkdir(parents=True, exist_ok=True)
    a_table = pd.DataFrame(report_rows(a_rows, {})).drop(columns=["Protein.Names"])
    pq.write_table(pa.Table.from_pandas(a_table, preserve_index=False),
                   p3 / "A2" / "A2_report.parquet")  # fmt: skip
    o = ents["O_OUT"]  # planted, plus strand: its frame is stop-free
    pep = translate(chrom[o.s1 - 1 : o.s1 - 1 + 36])
    for run in RUNS:
        d = a3 / run
        d.mkdir(parents=True, exist_ok=True)
        (d / f"{run}.ssl").write_text(
            "file\tscan\tcharge\tsequence\tscore-type\tscore\n"
            f"{PERSONAL}/{run}.mzML\t100\t2\t{pep}\tUNKNOWN\t0.9\n", encoding="utf-8")  # fmt: skip
    return chrom, write_a4(root, chrom, start, end, **a4)
