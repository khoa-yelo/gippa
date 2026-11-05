"""The text files of the example bundle and the release printout (M12-01, spec section 2).

Pure functions: each takes the facts it needs and returns a string. ``build_example.py`` writes
them; nothing here reads or writes files. The route commands come from
``gippa/data/example_routes.yaml`` (the single source for the printed commands, this README and the
docs page); ``{dir}`` is written as ``gippa-example`` and ``{out}`` as ``gippa-example/runs``.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

PXD = "PXD078459"
IPX = "IPX0017169000"
ASSEMBLY = "GCF_021559835.1"
GENBANK = "GCA_021559835.1"
PXD_RECORD = (
    "https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD078459&outputMode=XML"
)
PXD_RECORD_SHA256 = "5d3f16c3137fb11ecae7f79ad02063d8e3008b21380784fcfc8aa77b7fb533f9"
RUN_FILES = {
    "WT_1": "MLA1914689_WT_P1_1_Y-D7_1_13119.d.zip",
    "WT_2": "MLA1914690_WT_P2_2_Y-D8_1_13120.d.zip",
    "WT_3": "MLA1914691_WT_P3_3_Y-D9_1_13121.d.zip",
}
TOOLS = {"diann": "2.3.1 Academia", "cascadia": "0.0.7"}
REPO_URL = "https://github.com/khoa-yelo/gippa"
ASSET_HOST = "release-assets.githubusercontent.com"
TSV_COLUMNS = ("version", "file", "size_bytes", "sha256", "url", "release_tag", "build_commit",
               "build_job", "expected_commit", "asset_digest")  # fmt: skip
DOC_DIR = "gippa-example"


def archive_name(version: str) -> str:
    return f"gippa-example-{version}.tar.gz"


def asset_url(version: str) -> str:
    return f"{REPO_URL}/releases/download/v{version}.0/{archive_name(version)}"


def fill(cmd: str, directory: str = DOC_DIR) -> str:
    """A route command with ``{out}`` -> ``{dir}/runs`` and ``{dir}`` -> ``directory``."""
    return " ".join(cmd.split()).replace("{out}", "{dir}/runs").replace("{dir}", directory)


def readme(version: str, region: str, routes: Sequence[Mapping[str, Any]], expected: bool) -> str:
    blocks = []
    for route in routes:
        cmds = "\n".join(fill(str(c)) for c in route["commands"])
        blocks.append(f"### {route['title']}\n\n```bash\n{cmds}\n```\n")
    exp = (
        "`expected/` holds the key numbers of each route (`A.json` to `D.json`, with a tolerance "
        "per key) and `expected/README.md` shows them as a table, with the file and column each "
        "comes from."
        if expected
        else "This build has no `expected/` directory (the numbers are added at release prep)."
    )
    return f"""# GIPPA example bundle {version}

A real, cut-down dataset for trying GIPPA: one region of the *E. coli* Nissle 1917 chromosome,
`{region}` (1-based closed; 200 kb), and the peptides of the ProteomeXchange dataset {PXD}
(runs WT_1 to WT_3) that fall in it. The same region is given through four input routes: a DIA-NN
report (A), a protein FASTA with an abundance table (B), Cascadia de novo results (C) and a generic
peptide table (D). Nothing was re-searched; `SELECTION.json` holds the selection rule and the
sha256 of every source and output, and `PROVENANCE.md` says where each file comes from. The data
licence is in `LICENSE-DATA.md` (the data are not MIT-licensed).

Fetch it with `gippa example fetch`; the commands below assume the default `--outdir
gippa-example`. The routes need `gippa setup --no-genomad` (Bakta light and the antiSMASH
databases) and a container runtime (Apptainer or Docker).

## Routes

{chr(10).join(blocks)}
The runs use four faster real-tool settings: `--modules core,bgc`, the Bakta light database,
TransTermHP for terminators and Promoter Calculator off. A run with the defaults takes much
longer (see the GIPPA docs, "Try the example").

## Expected numbers

{exp}

## Caveats

- The numbers are for a region, not the genome: a peptide shared by a region gene and an
  out-of-region paralog looks unique here, and genome-wide statistics see few genes. Do not
  compare them with the full pilot.
- Route A's detection tier is `mayu`, not `decoy_picked`: on one region too few decoys pair with a
  gene.
- Route C keeps only the Cascadia rows whose peptide occurs in the region's six-frame translation,
  so every row is placeable: route C is a format demo, not a result.
- Route B needs GIPPA issue #536 fixed; before it, `gippa run --proteins` writes no gene table.
- The q-values, PEPs and MaxLFQ values of routes A, B and D are DIA-NN's, computed on the
  whole-proteome search (4,685 strain proteins plus contaminants, three runs). The subset does not
  recompute them.
- A Bakta full database, or an installed UniOP, can change the gene and operon numbers.
"""


def provenance(
    version: str,
    region: Mapping[str, Any],
    sources: Mapping[str, Mapping[str, str]],
    counts: Mapping[str, Any],
    build: Mapping[str, str],
    rewritten: Mapping[str, Sequence[str]],
) -> str:
    src = "\n".join(f"| `{k}` | `{v['path']}` | `{v['sha256']}` |" for k, v in sources.items())
    runs = "\n".join(f"- {r}: `{f}`" for r, f in RUN_FILES.items())
    a, c = counts["A"], counts["C"]
    c_rows = ", ".join(f"{r} {v['rows_kept']} of {v['rows_total']}" for r, v in c.items())
    scrub = "; ".join(f"`{k}`: {', '.join(f'`{c}`' for c in v)}" for k, v in rewritten.items())
    return f"""# Provenance of the GIPPA example bundle {version}

## Raw data

- ProteomeXchange {PXD} (iProX {IPX}); record {PXD_RECORD}, retrieved 2026-10-24, sha256
  `{PXD_RECORD_SHA256}`. Instrument: timsTOF HT (diaPASEF).
- Runs:
{runs}

## Genome

- *E. coli* Nissle 1917, {ASSEMBLY} (GenBank twin {GENBANK}), chromosome `{region["contig"]}`.
- Region: `{region["one_based_closed"]}` (1-based closed) = `[{region["start_0based"]},
  {region["end_exclusive"]})` (0-based half-open), {region["length_bp"]} bp, written as one linear
  contig `{region["output_contig_id"]}` (no topology tag).

## Routes

- A: DIA-NN {TOOLS["diann"]} report rows of the pilot's A2 search whose proteins all lie in the
  region; {a["target_rows_kept"]} target and {a["decoy_rows_kept"]} decoy rows of
  {a["rows_total"]}.
- B: the protein groups of route A's target rows with a MaxLFQ value;
  {counts["B"]["abundance_rows"]} proteins, sequences from the pilot's A2 strain database.
- C: Cascadia {TOOLS["cascadia"]} rows whose peptide occurs in the region's six-frame translation:
  {c_rows}.
- D: route A's rows grouped by (stripped peptide, run); {counts["D"]["rows"]} rows.
- Paths shortened to base names (R12): {scrub or "none"}.

## Sources (paths relative to `$GIPPA_ROOT`)

| Source | Path | sha256 |
|---|---|---|
{src}

## Build

- Script: `scripts/example/build_example.py` at commit `{build["commit"] or "unknown"}`.
- Bundle {version}, built with GIPPA {build["gippa_version"]}.
- `expected/`: {build["expected"]}.

`SELECTION.json` (next to this file) holds the selection rule verbatim, the counts and the sha256
of every source and output.
"""


LICENSE_DATA = f"""# Data licence

The GIPPA software is MIT-licensed. This licence does **not** cover the data files of this bundle.

The files under `A_diann/`, `B_proteins/`, `C_cascadia/` and `D_generic/` are derived from the
ProteomeXchange dataset **{PXD}** (iProX **{IPX}**), announced 2026-05-16. The dataset record
states no licence, and its publication is pending. The files are redistributed for the sole
purpose of demonstrating and testing GIPPA, with attribution to the depositors of {PXD}. The
content is unchanged: rows were selected and file paths shortened. The GIPPA authors grant no
licence to these files. Cite {PXD} if you use them. If you are a rights holder and object, open an
issue and the bundle will be withdrawn.

`genome/` holds a region of the *E. coli* Nissle 1917 assembly {ASSEMBLY} (GenBank {GENBANK}),
which is INSDC data with no use restrictions by NCBI. `expected/`, `README.md`, `PROVENANCE.md`,
`SELECTION.json`, `bundle.json` and `CHECKSUMS.sha256` are written by GIPPA tooling and are
MIT-licensed, like GIPPA.
"""


def db_sources_row(version: str, sha256: str, size: int, job: str) -> str:
    """The ``example_<X>_<Y>`` row of ``db_sources.yaml`` (spec section 1 > Pin)."""
    key = "example_" + version.replace(".", "_")
    job_note = f"build job {job}" if job else "local build"
    note = f"M12-01: the example bundle of GIPPA {version} (ADR-0074); never installed by setup"
    return f"""{key}:  # {note}
  version: "{version}"
  url: {asset_url(version)}
  checksum: sha256:{sha256}  # validation/example/bundles.tsv; {job_note}
  size_bytes: {size}
  path: example/{archive_name(version)}  # unused by setup; must be safe and relative
  tool: null
  check: sha256
  size_gib: 1
  redirect_hosts: [{ASSET_HOST}]"""


def bundles_tsv_row(
    version: str, sha256: str, size: int, commit: str, job: str, expected_commit: str
) -> str:
    """One ``validation/example/bundles.tsv`` row (``asset_digest`` empty until the upload)."""
    vals = (version, archive_name(version), str(size), f"sha256:{sha256}", asset_url(version),
            f"v{version}.0", commit, job, expected_commit, "")  # fmt: skip
    return "\t".join(vals)
