# Try the example

The example is real data, cut small enough to run in minutes. It shows the four ways into GIPPA:
a DIA-NN report, a protein FASTA with an abundance table, de novo results and a generic peptide
table. Every annotation tool runs live on it; nothing is precomputed.

## What the example is

The bundle holds a 200 kb region of the *E. coli* Nissle 1917 chromosome
(`NZ_CP058217.1:1965001-2165000`, 1-based and closed; assembly GCF_021559835.1) and the results
of GIPPA's pilot for the ProteomeXchange dataset **PXD078459** (three wild-type runs, WT_1 to WT_3),
cut to that region: the DIA-NN 2.3.1 report rows and the Cascadia 0.0.7 de novo results. Nothing
was searched again. The region holds the colibactin *pks* island and the yersiniabactin
high-pathogenicity island, so the `bgc` module has a gene cluster to find. It has no plasmid and
no prophage.

The archive is about 0.5 MB. One bundle is published with each minor release of GIPPA (1.1, 1.2,
…) as a GitHub release file; patch releases use their minor's bundle. It is not part of the
installed package. The data comes with licence terms of its own: see
[Licence and citation](#licence-and-citation) before you use it.

## Before you start

- GIPPA installed (`pip install gippa`, see the [Quick start](quickstart.md)).
- Databases and images from `gippa setup --no-genomad`. The routes need the Bakta light database
  (the default) and the antiSMASH databases; `gippa setup --minimal` leaves the antiSMASH
  databases out, so it is not enough.
- Apptainer or Docker, and Java 17 or newer, on `PATH`.

## Fetch the bundle

```bash
gippa example fetch
gippa example fetch --outdir <dir>
```

`gippa example fetch` downloads the bundle of your GIPPA version, checks its size and sha256
against the value pinned in GIPPA, unpacks it into `gippa-example` (or `--outdir <dir>`) and
checks every file against the bundle's `CHECKSUMS.sha256`. It then prints the commands of the four
routes with `<dir>` filled in as an absolute path. The output directory must be absent or empty.
Running the command again on a finished fetch prints `already fetched` and the commands, with no
network access; your outputs under `<dir>/runs` are left alone.

`--version X.Y` picks another bundle (this is the bundle version, not `gippa --version`). Without
it, GIPPA uses the bundle of its own `MAJOR.MINOR`; when there is none yet, it uses the newest
older bundle of the same major version and warns that the expected numbers were made with that
version.

### Without internet access

If the machine that runs GIPPA cannot reach GitHub, download the file on a machine that can (the
URL is in the error message) and copy it over:

```bash
gippa example fetch --from-file gippa-example-1.1.tar.gz --outdir <dir>
```

`--from-file` checks and unpacks the file in place and never opens a network connection. The file
itself is not changed.

### Exit codes

| Code | When |
|---|---|
| 0 | fetched and verified, or already fetched |
| 1 | your input: an unknown `--version`; the output directory exists and is not empty; it cannot be created; the `--from-file` file is missing, not a file, or not the pinned bundle. A file that is another bundle names it: `this is example bundle 1.2; use --version 1.2` |
| 2 | the environment: no bundle for this GIPPA major version; a network or HTTP error (a 404 adds that the bundle may not be published yet); a downloaded file whose size or sha256 differs from the pin; an unsafe or malformed archive; a file that does not match `CHECKSUMS.sha256`; a full disk or quota |
| 5 | an unexpected error |

A failed fetch leaves no partial directory behind. See [all exit codes](../user-guide/run.md#exit-codes).

## The four routes

The commands below assume the default bundle directory, `gippa-example`, in the current
directory. Each route writes its outputs to `gippa-example/runs/<route>`. They are the commands
`gippa example fetch` prints; the printed ones use the absolute path of your bundle directory, so
copy those if you used `--outdir`. Run each route's commands in order; the routes are
independent.

### Route A: DIA-NN report

```bash
gippa convert --from diann gippa-example/A_diann/report.parquet --out gippa-example/runs/A/peptides.tsv --search-db same
gippa preflight --genome gippa-example/genome/nissle_r200k.fna.gz --peptides gippa-example/runs/A/peptides.tsv --set 'modules=[core,bgc]' --set annotation.bakta.db=light --set annotation.terminators.tool=transtermhp --set annotation.promoter_calculator.enabled=false
gippa run --genome gippa-example/genome/nissle_r200k.fna.gz --peptides gippa-example/runs/A/peptides.tsv --genome-relation same --outdir gippa-example/runs/A --modules core,bgc --set annotation.bakta.db=light --terminators transtermhp --set annotation.promoter_calculator.enabled=false --threads 4
```

### Route B: protein FASTA and abundance table

```bash
gippa preflight --genome gippa-example/genome/nissle_r200k.fna.gz --proteins gippa-example/B_proteins/proteins.faa --abundance gippa-example/B_proteins/abundance.tsv --set 'modules=[core,bgc]' --set annotation.bakta.db=light --set annotation.terminators.tool=transtermhp --set annotation.promoter_calculator.enabled=false
gippa run --genome gippa-example/genome/nissle_r200k.fna.gz --proteins gippa-example/B_proteins/proteins.faa --abundance gippa-example/B_proteins/abundance.tsv --abundance-type maxlfq --genome-relation same --outdir gippa-example/runs/B --modules core,bgc --set annotation.bakta.db=light --terminators transtermhp --set annotation.promoter_calculator.enabled=false --threads 4
```

### Route C: Cascadia de novo results (a format demo)

```bash
gippa convert --from cascadia gippa-example/C_cascadia/WT_1.ssl gippa-example/C_cascadia/WT_2.ssl gippa-example/C_cascadia/WT_3.ssl --out gippa-example/runs/C/peptides.tsv
gippa preflight --genome gippa-example/genome/nissle_r200k.fna.gz --peptides gippa-example/runs/C/peptides.tsv --set 'modules=[core,bgc]' --set annotation.bakta.db=light --set annotation.terminators.tool=transtermhp --set annotation.promoter_calculator.enabled=false
gippa run --genome gippa-example/genome/nissle_r200k.fna.gz --peptides gippa-example/runs/C/peptides.tsv --genome-relation same --outdir gippa-example/runs/C --modules core,bgc --set annotation.bakta.db=light --terminators transtermhp --set annotation.promoter_calculator.enabled=false --threads 4
```

### Route D: a generic peptide table

```bash
gippa convert --from generic gippa-example/D_generic/peptides.tsv --out gippa-example/runs/D/peptides.tsv --source database --search-db same
gippa preflight --genome gippa-example/genome/nissle_r200k.fna.gz --peptides gippa-example/runs/D/peptides.tsv --set 'modules=[core,bgc]' --set annotation.bakta.db=light --set annotation.terminators.tool=transtermhp --set annotation.promoter_calculator.enabled=false
gippa run --genome gippa-example/genome/nissle_r200k.fna.gz --peptides gippa-example/runs/D/peptides.tsv --genome-relation same --outdir gippa-example/runs/D --modules core,bgc --set annotation.bakta.db=light --terminators transtermhp --set annotation.promoter_calculator.enabled=false --threads 4
```

Route D's table was made from route A's rows, as an engine without a GIPPA converter would write
it. To make such inputs from your own spectra, see [From raw spectra](from-spectra.md).

## How long it takes

About 5 min per route for `gippa run` on 4 CPUs (route B about 8.5 min), after `gippa setup` and with the container images
already pulled (setup and image pulls are not counted). `gippa convert` and `gippa preflight` take
seconds. Measured with the commands above:

| Route | `gippa run` (pipeline time) |
|---|---|
| A | 5:00 to 5:05 |
| B | 8:31 |
| C | 5:07 |
| D | 4:59 |

More CPUs do not make it faster: the steps run one after another, and each one spends most of
its few seconds starting up.

### Why the example is faster than a default run

The commands use four settings that keep every tool real but drop the steps whose time does not
shrink with a smaller genome:

| Default | Example setting | Why |
|---|---|---|
| Bakta full database | `--set annotation.bakta.db=light` | Bakta with the full database took 26 min on 55 kb alone; with the light database about 1 min |
| TerminatorNet | `--terminators transtermhp` | TerminatorNet took 4.5 to 4.8 min per strand of 100 kb on 4 CPUs; TransTermHP takes about 1 s |
| all modules | `--modules core,bgc` | geNOMAD (`mobile`) takes about 5 min whatever the size; `function` finds no COG, GO or KEGG IDs in the light database |
| Promoter Calculator on | `--set annotation.promoter_calculator.enabled=false` | 428 s on 200 kb; Promotech still runs |

With the defaults the example is much slower. Measured on route A:

| Region | Settings | `gippa run` (pipeline time) |
|---|---|---|
| 55 kb | `core`, Bakta full database | not finished; Bakta alone took 26 min |
| 100 kb | all modules, Bakta light database | 13:45 |
| 100 kb | `core`, Bakta light database | 8:51 |
| 200 kb | `core`, Bakta light database, TransTermHP, Promoter Calculator on | 10:29 |

A run on the 200 kb region with every default (full database, all modules) was not measured.
For a whole genome with the defaults, see the [quick start](quickstart.md).

## What results to expect

| Route | Detection tier | Genes | Detected |
|---|---|---|---|
| A | `mayu` | 158 | 99 (gene-level FDR 0.012 at q 0.01) |
| B | `protein_list` | 158 | 106 |
| C | `two_peptide` | | 4 (from 76 converted rows) |
| D | `mayu` | see `expected/D.json` | see `expected/D.json` |

- **Routes A and D are `mayu`, not `decoy_picked`.** A 200 kb region keeps too few decoys: 37
  decoy rows remain, but only 4 of 10 decoy peptides pair with a gene, below the 80 % that
  `decoy_picked` needs, so the run warns `decoy_pairing_low` and uses the `mayu` tier. Route D
  carries the same decoys as route A.
- **Route B** has no peptides, so its tier is `protein_list`.
- **Route C is a format demo, not a result.** The bundle keeps only Cascadia rows that match the
  region exactly, which inflates the placement rate; 4 genes are detected.
- **Region effects.** The q-values, PEPs and MaxLFQ values come from DIA-NN's search of the whole
  proteome; the subset does not recompute them. A peptide shared with a gene outside the region
  looks unique here. Do not compare these numbers with whole-genome runs.
- **Empty on this region:** genomic islands, variants and novel ORFs (0 each), enrichment and GSEA
  (the light database has no categories), and autocorrelation (`not_estimable`). Route A finds 46
  operons (4 flagged) and one gene cluster with 30 of its 71 genes detected, against 42.4
  expected.
- Preflight info rows such as a short genome or small replicons are expected on a region; they are
  not faults.

`gippa-example/expected/` holds each route's numbers (`A.json` to `D.json`) and their tolerances, and
`expected/README.md` shows them as a table with the file and column each one comes from. Other
machines, thread counts or container runtimes may differ within those tolerances. A Bakta full
database or an installed UniOP can change the gene and operon numbers.

To read the four runs and compare them, follow the tutorial
[One Nissle region, four routes](../tutorials/example-routes.md).

## Licence and citation

- **No licence is granted for the PXD078459-derived files** (`A_diann/`, `B_proteins/`,
  `C_cascadia/`, `D_generic/`). The dataset record states no licence. The files are provided only
  to demonstrate and test GIPPA, with attribution to the depositors; rows were selected and file
  paths shortened, nothing else was changed.
- **Cite PXD078459** (iProX **IPX0017169000**) if you use them.
- **Withdrawn on request.** If you hold rights to the data and object, open an issue and the
  bundle will be withdrawn.
- `genome/` is a region of the INSDC assembly GCF_021559835.1 (GenBank GCA_021559835.1), which
  NCBI places under no use restrictions. The bundle's own files (`expected/`, the README and the
  checksums) are MIT-licensed, like GIPPA.

The tools the routes run keep their own terms ([Licences](../user-guide/licences.md)). The example
needs neither of the two restricted resources: geNOMAD (non-commercial use only) belongs to the
`mobile` module, which the example switches off, so `gippa setup --no-genomad` is enough; and
KEGG (free for academic use, others need a KEGG licence) is used only by `function`, which is off
too, so `gippa setup --no-genomad --no-kegg` also works. The GPL tools the routes run (Bakta,
Promotech, TransTermHP) and antiSMASH (AGPL) run as separate programs in their own containers;
GIPPA does not link to them.
