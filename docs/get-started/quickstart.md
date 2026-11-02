# Quick start

This page takes a new user from nothing to a finished GIPPA run: install, set up, check, run, and
which settings to change. Each step links to the page with the details.

Paths written in angle brackets, such as `<db-dir>`, are yours to choose. The notes next to them
say what each one needs.

## The short version

```bash
pip install gippa
gippa demo                                   # about 2 minutes; prints the report path
gippa setup --db-dir <db-dir> --profile apptainer
gippa run --genome genome.fna --peptides peptides.tsv --outdir out --profile apptainer
```

## 1. Install

```bash
python3 -m venv <venv-dir>                   # any directory; Python 3.11 or newer
source <venv-dir>/bin/activate
pip install gippa                            # on glibc < 2.28 (CentOS 7): pip install --prefer-binary gippa
gippa --version
```

`gippa run` also needs these on `PATH` ([Install › Prerequisites](install.md#prerequisites)):

- **Java 17** or newer, for Nextflow;
- **a container runtime**: Apptainer (no root needed on most clusters) or Docker;
- **Nextflow** 25.04.8 or newer. If none is on `PATH`, `gippa setup` installs the pinned 26.10.0
  launcher for you.

## 2. Check the install: `gippa demo`

```bash
gippa demo
```

The demo runs GIPPA on a small synthetic genome bundled in the package. It needs no Nextflow, no
containers, no databases and no network. It takes about 2 minutes and ends with the path of the
report, for example `gippa demo: done in 48 s; report: <cwd>/gippa-demo/report/report.html`.
The numbers describe a simulated genome, so do not read them as biology.

## 3. Set up the databases and images (once)

```bash
gippa setup --db-dir <db-dir> --profile apptainer     # or --profile docker
```

- `<db-dir>`: a directory with about 50 GB free. On a cluster, it must be on a filesystem that
  the compute nodes can read, and home directories are often too small.
- Setup prints a plan with every item, its size and the free space, then asks `Proceed? [Y/n]`.
  The default databases take about 20 GB and the images 15–25 GB. Most of the time is spent
  downloading.
- Set `NXF_APPTAINER_CACHEDIR=<image-dir>` first to keep the Apptainer images somewhere other
  than the default cache.
- On a cluster, run setup as a batch job, not on a login node.
- Setup records `<db-dir>`, so `gippa run`, `gippa preflight` and `gippa makedb` find it without
  `--db-dir`.

Licences and options to leave things out:

| Flag | Leaves out | Then run with |
|---|---|---|
| `--no-genomad` | geNOMAD, which is for non-commercial use only | `--modules core,bgc,function` |
| `--no-kegg` | KEGG, which is free for academic use only | (function runs without KEGG) |
| `--minimal` | geNOMAD and antiSMASH, for small disks | `--modules core,function` |
| `--bakta full` | adds the full Bakta database, about 84 GB more | |

UniOP, the operon predictor, is built on your machine only after you agree, at the prompt or with
`--accept-uniop-terms`. Without it, runs use a distance-based operon model.
[Set up the databases](setup-databases.md) lists every flag; [Licences](../user-guide/licences.md)
gives the terms.

## 4. Prepare your inputs

GIPPA needs **one assembled bacterial genome** (FASTA) and **one table of identified peptides**,
or a protein list. It does not search spectra itself.

| You have | Do this |
|---|---|
| a peptide table from a search engine | pass it with `--peptides` (format: [Inputs](../user-guide/inputs.md#peptide-table-peptidestsv)) |
| DIA-NN output | `gippa convert --from diann report.parquet --out peptides.tsv`, then `--peptides peptides.tsv` |
| Cascadia, Casanovo or InstaNovo (de novo) | `gippa convert --from cascadia …`, then `--peptides`, ideally with `--reference-peptides` |
| a protein list (FASTA) | pass it with `--proteins` |
| raw spectra and no table yet | [From raw spectra to a GIPPA run](from-spectra.md): `gippa makedb`, DIA-NN, `gippa convert` |

- Mark each circular record with ` [topology=circular]` in its FASTA header, or pass
  `--set annotation.bakta.complete=true` for a complete assembly.
- Keep the search engine's decoy rows: with them GIPPA estimates gene-level FDR.
- If you already have a Bakta annotation, pass `--annotation genome.bakta.json` and Bakta is not run.

## 5. Check, then run

```bash
gippa preflight --genome genome.fna --peptides peptides.tsv --profile apptainer   # seconds
gippa run --genome genome.fna --peptides peptides.tsv --outdir <out-dir> --profile apptainer
```

Preflight checks the inputs and the installed tools, and exits with code 0 when the run can
start. A run takes hours, not minutes: a standard E. coli run took about 3–4 h on 16 threads and
64 GB ([Runtimes](../user-guide/runtimes.md)).

**On a Slurm cluster**, submit the run as one job. With `--profile apptainer`, every task runs
inside that job's allocation:

```bash
#!/bin/bash
#SBATCH --job-name=gippa-run
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=8:00:00
source <venv-dir>/bin/activate
gippa run --genome genome.fna --peptides peptides.tsv --outdir <out-dir> \
    --threads 16 --memory 64.GB --profile apptainer
```

`--profile slurm,apptainer` instead submits each task as its own Slurm job, on the cluster's
default partition. It is not tested yet, and every node then needs the same `gippa` on `PATH`
([Install › Execution profiles](install.md#execution-profiles)). If a run stops, rerun the same
command with `--resume`.

The report is `<out-dir>/report/report.html`, one self-contained file. Copy the whole output
directory to your computer to open it, so that its links to the tables keep working.
[Outputs](../user-guide/outputs.md) describes every file.

## 6. The settings you can change

Every setting has a default, so you only give the ones you change. Give them as flags, in a YAML
file passed with `--config`, or one at a time with `--set key=value`. The order of precedence is:
defaults, then the config file, then flags, then `--set`. A misspelled key stops the run with exit
code 1. Every run writes the values it used to `<out-dir>/resolved_config.yaml`.

```yaml
# my.yaml: only the keys you change
threads: 16
modules: [core, mobile, bgc, function]
genome_relation: same
abundance:
  type: ibaq
annotation:
  terminators:
    tool: transtermhp
```

```bash
gippa run --genome genome.fna --peptides peptides.tsv --outdir <out-dir> --config my.yaml \
    --set seed=7
```

The common settings:

| Flag | Config key | Options (default first) | What it means |
|---|---|---|---|
| `--modules` | `modules` | `core,mobile,bgc,function` | which analyses run. `core` (annotation, mapping, the gene table, detection, operons, variants, novel ORFs) is always on; `mobile` adds prophages and plasmids (geNOMAD), `bgc` biosynthetic gene clusters (antiSMASH), `function` COG, GO and KEGG enrichment |
| `--genome-relation` | `genome_relation` | `related`, `same` | `same` when the proteomics comes from this exact strain; `related` for a related strain |
| `--abundance-type` | `abundance.type` | `unknown`, `ibaq`, `top3`, `maxlfq`, `intensity` | what the abundance values in your input are |
| `--abundance-scale` | `abundance.scale` | `auto`, `log2`, `linear` | whether those values are already log-transformed |
| `--abundance-missing` | `abundance.missing` | `median`, `exclude`, `min`, `zero`, `mean` | how a gene with no abundance value is filled |
| `--promoters` | `annotation.promoters.mode` | `targeted`, `genome` | predict promoters only upstream of genes, or across the whole genome |
| `--terminators` | `annotation.terminators.tool` | `terminatornet`, `transtermhp` | terminator predictor; TerminatorNet took most of the pilot runs' time |
| `--add-bactermfinder` | `annotation.terminators.bactermfinder` | off | also run BacTermFinder |
| `--sensitivity` | `output.sensitivity` | `basic`, `full`, `off` | how many re-runs check that the conclusions are robust |
| `--publication` | | off | 10,000 permutations and bootstrap draws instead of the default, for final numbers; slower |
| `--parquet` | `output.parquet` | off | also write Parquet copies of the large tables |
| `--seed` | `seed` | `20261014` | the master seed; the same seed and inputs give the same results |
| `--threads` | `threads` | `16` | the most threads one task uses |
| `--memory` | | `64.GB` | the memory the workflow may use |
| `--profile` | | `auto` | where tools run: `apptainer`, `docker`, `conda`, `slurm`, `standard`, combined with commas; `auto` detects one and prints its choice |
| `--strict` | | off | treat warnings as errors (exit code 3) |

Error-rate settings: `input.q_value_max` (default `0.01`) is the q-value above which a peptide is
not evidence, and `stats.alpha.identification` (default `0.01`) is the gene-level FDR for calling a
gene detected. The other levels are under `stats.alpha` (`placement`, `novel_orfs`, `variants`,
`context`, `enrichment`, `operon_flags`). GIPPA's error rates were validated at the defaults.

[Running GIPPA](../user-guide/run.md) explains every option, and the
[configuration reference](../reference/configuration.md) lists every key.

## If this fails

!!! failure "No container runtime found"
    `gippa run` exits with code 2 and prints `no container runtime or conda found: install
    Apptainer (no root needed on most Linux systems) or Docker`. Install one, then run
    `gippa setup --images --profile apptainer` with the same `<db-dir>`, so that it pulls the
    images.

!!! failure "Java older than 17"
    Nextflow needs Java 17 or newer. Check with `java -version`. Without root, for example:
    `micromamba create -n gippa-java -c conda-forge openjdk=21`.

!!! failure "Not enough disk"
    Setup stops before it downloads anything, with exit code 2 and a message such as
    `needs 31.4 GiB, 12.0 GiB free on /home/u: pass --db-dir (or set GIPPA_DB_DIR) on a larger
    disk`. Point `--db-dir` at a larger filesystem, or install less with `gippa setup --minimal`.

## Next

- [From raw spectra to a GIPPA run](from-spectra.md): DIA-NN, six-frame ORFs and de novo.
- [Inputs](../user-guide/inputs.md): peptide tables, protein lists, abundance and annotations.
- [Interpreting results](../user-guide/interpreting-results.md): what the q-values mean.
