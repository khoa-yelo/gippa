# Quick start

From nothing to a finished GIPPA run. Paths in angle brackets, such as `<db-dir>`, are yours to
choose. Expect about 5 min to install, about 1.5 h for setup on 4 CPUs and 3–4 h per run
([measured times](install.md#how-long-each-step-takes)).

## 1. Install and check

Use Python 3.12 (3.11 also works). If `python3 --version` is older, get 3.12 first
([Install › Python](install.md#python)).

```bash
python3.12 -m venv <venv-dir> && source <venv-dir>/bin/activate
pip install --prefer-binary gippa   # prefer wheels: no numpy/numba builds from source
gippa demo                          # about 2 min on bundled synthetic data; prints the report path
```

`gippa run` also needs Java 17 or newer and Apptainer or Docker on `PATH`. `gippa setup` installs
Nextflow if none is found ([Install](install.md#prerequisites)).

## 2. Set up databases and images (once)

```bash
gippa setup --db-dir <db-dir> --profile apptainer --yes --accept-uniop-terms   # or --profile docker
```

`--yes` skips the "Proceed?" question after the download plan. `--accept-uniop-terms` means you
decide that your use of UniOP, the operon predictor, is permitted: it has no licence file.
`<db-dir>` needs about 50 GB free where your compute nodes can read it; on a cluster, run setup as
a batch job. Flags you may change ([all flags](setup-databases.md#flags)):

- `--no-genomad` for commercial use (geNOMAD is non-commercial only), then run with
  `--modules core,bgc,function`;
- `--no-kegg` outside academia without a KEGG licence, then run with `--modules core,mobile,bgc`;
- `--bakta full`: functions for more genes, about 84 GB more; `--minimal` for small disks.

## 3. Run

GIPPA takes one assembled bacterial genome and a table of identified peptides. It does not search
spectra: convert DIA-NN output with `gippa convert --from diann report.parquet --out
peptides.tsv`, or see [From raw spectra](from-spectra.md). Keep the decoy rows.

```bash
gippa preflight --genome genome.fna --peptides peptides.tsv --profile apptainer   # seconds
gippa run --genome genome.fna --peptides peptides.tsv --outdir <out-dir> --profile apptainer
```

A run takes about 3–4 h on 16 threads and 64 GB. On Slurm, submit it as one job:

```bash
#!/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=8:00:00
source <venv-dir>/bin/activate
gippa run --genome genome.fna --peptides peptides.tsv --outdir <out-dir> --profile apptainer
```

The report is `<out-dir>/report/report.html`. Rerun with `--resume` if a run stops.

## 4. Settings you can change

Use a flag, a YAML file (`--config my.yaml`) or `--set key=value`; later ones win. A misspelled key
is an error, and every run writes what it used to `<out-dir>/resolved_config.yaml`.

```yaml
# my.yaml: only what you change
genome_relation: same
abundance:
  type: ibaq
```

| Flag | Options (default first) | Meaning |
|---|---|---|
| `--modules` | `core,mobile,bgc,function` | analyses to run: `core` always; `mobile` prophages and plasmids, `bgc` gene clusters, `function` COG/GO/KEGG |
| `--genome-relation` | `related`, `same` | `same` if the proteomics comes from this exact strain |
| `--abundance-type` | `unknown`, `ibaq`, `top3`, `maxlfq`, `intensity` | what your abundance values are |
| `--abundance-missing` | `median`, `exclude`, `min`, `zero`, `mean` | how genes without a value are filled |
| `--terminators` | `terminatornet`, `transtermhp` | terminator predictor; TerminatorNet is one of the two slowest steps |
| `--promoters` | `targeted`, `genome` | promoters upstream of genes only, or genome-wide |
| `--sensitivity` | `basic`, `full`, `off` | robustness re-runs |
| `--publication` | off | 10,000 permutations for final numbers; slower |
| `--threads`, `--memory` | `16`, `64.GB` | resources |
| `--seed` | `20261014` | same seed and inputs give the same results |
| `--profile` | `auto` | `apptainer`, `docker`, `conda`, `slurm`; `auto` detects and prints its choice |

Error rates: `input.q_value_max` (`0.01`, peptide evidence) and `stats.alpha.identification`
(`0.01`, gene-level FDR). They were validated at these defaults.
All options: [Running GIPPA](../user-guide/run.md); all keys: [configuration reference](../reference/configuration.md).
