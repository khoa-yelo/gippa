# Install

GIPPA has two parts:

- the Python package, which provides the `gippa` command and bundles the Nextflow workflow, its
  container pins (`workflow/conf/containers.config`, `workflow/conf/images.json`) and the database
  pins (`src/gippa/data/db_sources.yaml`);
- the tools it runs as separate programs (Bakta, miniprot, geNOMAD, antiSMASH and others). Each
  tool runs from a container image pinned by digest; UniOP, which GIPPA may not redistribute, is
  built on your machine by `gippa setup` ([UniOP](#uniop)).

`gippa preflight`, `gippa convert`, `gippa compare` and `gippa report` need only the Python package.
`gippa setup` needs no Nextflow (it installs the pinned launcher when none is on `PATH`). `gippa
run` needs Nextflow, Java 17 or newer, a container runtime, the tool images and the databases.

## Tested platforms

The full pipeline is tested on Linux x86-64 with Slurm and Apptainer. Every other row is as stated:

--8<-- "docs/user-guide/platforms.md:platforms"

The cell words (tested, smoke-tested, untested) and the CI jobs are explained on
[Platforms](../user-guide/platforms.md).

## Prerequisites

| Requirement | Version | Where it is pinned |
|---|---|---|
| Python | Python ≥ 3.11 | `requires-python = ">=3.11"` in `pyproject.toml`; 3.12 is the version to use ([Python](#python)) |
| Nextflow | 26.10.0 is the pinned, tested version; at least 25.04.8 (advice: GIPPA does not check the version) | the launcher pin `nextflow` in `src/gippa/data/db_sources.yaml` |
| Java | 17 or newer, which Nextflow needs; tested with OpenJDK 21 | checked by `gippa run` (`src/gippa/runner/prereqs.py`) |
| Docker or Apptainer | the `docker` CLI with access to a daemon, Docker 20.10 or later (podman 3.0 or later through podman-docker: best effort, untested), or `apptainer` | see [Platforms](../user-guide/platforms.md) |

### Python

Use Python 3.12: the release package is built and smoke-tested with it. CI also tests 3.11 and
3.13; any Python 3.11 or newer is accepted. When the system `python3` is older, get 3.12 one of
these ways, then make the venv with it:

```bash
module avail python; module load python/3.12        # clusters: the module name varies by site
uv python install 3.12 && uv venv --seed -p 3.12 <venv-dir>    # uv, no root; --seed adds pip
micromamba create -n gippa -c conda-forge python=3.12          # or conda; then activate it
```

**Nextflow.** When `nextflow` is not on `PATH`, `gippa setup` downloads the official
self-contained launcher `nextflow-26.10.0-dist` from the `nextflow-io/nextflow` GitHub release
`v26.10.0`, checks its pinned sha256 and installs it as `DIR/bin/nextflow` (recorded in
`DIR/db_manifest.json` as item `nextflow`). The `-dist` launcher downloads nothing more at its
first run, so it also works on nodes without internet access. `--no-nextflow` skips it. When a
`nextflow` is on `PATH`, setup names its version and path and downloads nothing. `gippa run` uses
`nextflow` on `PATH`, else `DIR/bin/nextflow`, and exits with code 2 when it finds neither.

**Java** is never downloaded. `gippa run` runs `$JAVA_HOME/bin/java -version` (else `java` on
`PATH`) and exits with code 2 before Nextflow starts when Java is missing or older than 17;
`gippa setup` only warns, since the database downloads need no Java. Install Java from
conda-forge (never the `defaults` channel), for example with micromamba:

```bash
micromamba create -n gippa-java -c conda-forge openjdk=21
micromamba activate gippa-java
java -version
```

To install Nextflow yourself instead, use the same environment:
`micromamba create -n gippa-nf --override-channels -c conda-forge -c bioconda nextflow=26.10.0 openjdk=21`.

GIPPA supports genetic codes 11 (the default), 4 (Mycoplasma and Spiroplasma) and 25 (candidate
division SR1 and Gracilibacteria). Preflight refuses any other `annotation.translation_table`, and
warns (`genome.genetic_code`) when the genome's open reading frames suggest TGA does not act the
way the configured code says (it tells 11 from 4/25, not 4 from 25)
([Inputs › Genetic code](../user-guide/inputs.md#genetic-code)).

## The five steps

```bash
# 1. the Python package (Python 3.12; 3.11 also works)
pip install --prefer-binary gippa
gippa --version

# 2. Nextflow, Java and Docker or Apptainer on PATH (Prerequisites above)

# 3. pinned images, databases and the local UniOP build
export NXF_APPTAINER_CACHEDIR=/big/disk/gippa-images     # optional, apptainer only (below)
gippa setup --db-dir /path/to/db --profile apptainer \
  --accept-download-plan --install-uniop-from-git      # or: --profile docker; no prompts

# 4. check the inputs and the installed tools
GIPPA_DB_DIR=/path/to/db gippa preflight --genome genome.fna --peptides peptides.tsv \
    --profile apptainer                                  # or: --profile docker

# 5. run
gippa run --genome genome.fna --peptides peptides.tsv --db-dir /path/to/db \
    --outdir out --profile apptainer                     # or: --profile docker
```

`pip install --prefer-binary` makes pip choose the newest release that has a wheel for your
system over a newer one that it would have to build from source. Use it everywhere. On older Linux
(glibc < 2.28, such as CentOS or RHEL 7) a plain `pip install gippa` fails, because the newest
numpy, numba and pyarrow have no wheels there. With the flag, on a CentOS 7 cluster (glibc 2.17,
Python 3.12), pip chose numpy 2.2.6, numba 0.67.0, pyarrow 20.0.0 and lightgbm 3.3.5, and
`gippa demo` passed. These older versions are not part of GIPPA's tested matrix
([#502](https://github.com/khoa-yelo/gippa/issues/502)).

**When `pip install` fails:**

- `ERROR: Could not find a version that satisfies the requirement gippa (from versions: none)`:
  the Python is older than 3.11. Make the venv with Python 3.12 ([Python](#python)).
- pip builds numba or numpy from source and stops with `NumPy requires GCC >= 10.3` (or another
  compiler error): rerun with `pip install --prefer-binary gippa`
  ([#502](https://github.com/khoa-yelo/gippa/issues/502)).

`uv tool install gippa` works in place of `pip install --prefer-binary gippa`. The extra
`gippa[docs]` installs MkDocs Material for building this site. **Bioconda:** the recipe is
pending (M11-10); until it is merged, use PyPI.

**Step 3.** `gippa setup` pulls the images, then downloads the databases. Every item is on by
default; each has an opt-out such as `--no-antismash` (the flags are on
[Set up the databases](setup-databases.md)). Before downloading, setup prints the licence terms
of geNOMAD (non-commercial use only) and KEGG (academic use only).

`gippa setup` asks at most two questions, and only when stdin is a terminal. No other `gippa`
command asks anything. To run setup without questions, pass the flags in the last column:

| Question | When | Flag that answers it |
|---|---|---|
| `gippa setup: Proceed? [Y/n]` (default yes) | after the download plan, before the first download or pull | `--accept-download-plan`; in a batch job (no terminal) it is not asked and setup goes on |
| `Build UniOP for your own use? [y/N]` (default no) | after `Proceed?`, when setup builds UniOP: with `--profile` when the modules need it, or with `--uniop` | `--install-uniop-from-git` (yes: you decide your use is permitted) or `--no-uniop` (skip UniOP); in a batch job without either, UniOP is skipped ([UniOP](#uniop)) |

- With `--profile`, setup pulls the pinned images that the configuration's modules need;
  `--images all` pulls every pinned image and `--no-images` none. Without `--profile` the images
  are skipped with a warning. `--config` and `--set` select the configuration, as for
  `gippa run`.
- `--profile apptainer` runs `apptainer pull <cache>/<name>.img docker://<ref>`, where `<name>` is
  the file name Nextflow looks for ([ADR-0059](../adr/ADR-0059-bare-image-refs.md)).
  `--profile docker` runs `docker pull <ref>` and checks that the pulled digest equals the pin.
  A comma-joined profile is accepted for the images, but the UniOP step needs one runtime:
  `--profile docker,conda` exits 2 with `ambiguous UniOP runtime for docker,conda: use one of
  conda | docker | apptainer`.
- Each image is recorded as an `image:<tool>` item in `DIR/db_manifest.json`. A failed pull or a
  digest mismatch (`pulled <d1>, pinned <d2>`) is exit code 2; the other images are still pulled.
- **A failed Apptainer pull.** On some unprivileged Apptainer builds (seen with
  apptainer 1.5.5 on one cluster), `mksquashfs` can segfault on a first pull: the message ends in
  `mksquashfs command failed: exit status 139`. Setup retries a failed pull once automatically.
  Setting `APPTAINER_CACHEDIR` keeps the layers between attempts, so a rerun can reuse them;
  without it each attempt pulls from scratch. On that cluster the failure followed the node type,
  not memory: it repeated on Intel Broadwell nodes and setup passed on AMD ones. If it keeps
  failing, run setup on a different node type (with Slurm, `sbatch --constraint=<feature>`). A
  missing bakta image also makes `bakta_light` skip, so rerun setup after such a failure. Large
  pulls print no progress: the antiSMASH image (4.6 GiB) can take several minutes
  ([#504](https://github.com/khoa-yelo/gippa/issues/504)).
- **The DB directory's nearest existing parent must be writable.** Otherwise setup stops with a
  clean error (exit code 1)
  ([#503](https://github.com/khoa-yelo/gippa/issues/503)).
- **The Apptainer cache** is one directory: `NXF_APPTAINER_CACHEDIR`, else
  `DIR/containers` (`setup --cache-dir` overrides it for that
  `setup` command). Every reader uses that same directory:
    - setup's image pull and UniOP build;
    - setup's Bakta and antiSMASH downloaders, which run inside the pulled images;
    - preflight under `--profile` (`gippa preflight --profile` and `gippa run`);
    - Nextflow: `gippa run` exports `NXF_APPTAINER_CACHEDIR` to the cache it checked when the
      variable is not set (`apptainer.cacheDir`).

  `DIR/containers` therefore works with no variable set. Export `NXF_APPTAINER_CACHEDIR` to put
  the images on another disk, or to share one cache with Nextflow runs outside `gippa`; keep it
  set for setup, preflight and run. A `--cache-dir` other than that directory is seen only by
  that `setup` command: preflight then reports the images `not in the apptainer cache`.
- On a machine with Apptainer and no Docker, the Bakta (and, with `--antismash`, antiSMASH)
  downloaders run inside the tool image: setup pulls the images first, so one `setup --images`
  run installs those databases too.
- Compute nodes without network: run step 3 where the network is. Preflight checks the cache,
  not the registry. To download the databases and images yourself, see
  [Download the databases by hand](manual-databases.md).

**Step 4.** `gippa preflight` validates the inputs in seconds and reads the installed versions from
`$GIPPA_DB_DIR/db_manifest.json`. With `--profile`, it also checks that every tool is installed
under that profile, exactly as `gippa run --profile` does before Nextflow starts (next section):
the same profile names and the DB directory from `$GIPPA_DB_DIR`. A missing tool is exit code 2; for missing images the message says `run: gippa setup
--images --profile apptainer`. Without `--profile`, no tool is checked.

**Step 5.** The default `--profile` is `auto`: `gippa run` chooses a profile from what it finds on
the host and prints its choice ([Execution profiles](#execution-profiles)); pass `--profile` to
choose yourself.
`--memory` defaults to `64.GB`; set it to what the machine has. `--threads` defaults to 16, capped to
the CPUs available unless you set it or a Slurm profile runs the tasks. See
[Run](../user-guide/run.md).

### UniOP

UniOP, the operon predictor, has no licence file, so GIPPA never redistributes it: there is no
UniOP image. `setup --images` (or `setup --uniop` alone) fetches the authors' archive at the pinned
commit from GitHub, checks its sha256, and builds it on your machine. Nothing is pushed anywhere.

- **Terms.** Setup prints this notice whenever it builds or checks a local build:
  "UniOP (github.com/hongsua/UniOP @3ed14dfc) is fetched from its authors and built on this
  machine; it has no licence file, so gippa does not distribute it and whether your use is
  permitted is yours to decide." On a terminal it asks `Build UniOP for your own use? [y/N]`; the
  default is no. For scripts and batch jobs, pass `--install-uniop-from-git`. Without either, nothing
  is fetched: `setup --images` warns `image:uniop: terms not accepted; rerun with
  --install-uniop-from-git` and exits 0; `setup --uniop` exits 2.
- **Runtime.** micromamba when a profile enables conda, else Apptainer, else Docker. Apptainer
  runs `apptainer build --fakeroot` (plain `apptainer build` as root), and retries once as
  `apptainer build --fakeroot --ignore-fakeroot-command`. The build is named in
  `DIR/local_tools.config`, which `gippa run` passes to Nextflow.
- **When the build fails** (no network, no Docker daemon access, no fakeroot), setup records the
  reason and exits 2 (once the terms are accepted, under `--images` and `--uniop` alike). The ways
  out:
    - run `setup` again on a machine that can build, with the same DB directory and cache on a
      filesystem the run host sees (preflight checks the run host: the Docker image ID or the
      `.sif` sha256 must equal the recorded build);
    - switch profile, as the Apptainer message says: `use --profile docker or conda, or set
      annotation.uniop.enabled=auto`;
    - or run without UniOP. `annotation.uniop.enabled` is `auto` by default: with UniOP missing,
      preflight warns `uniop not installed (unlicensed, not distributed); operons use the distance
      model; run: gippa setup --uniop --install-uniop-from-git --profile <profile> --db-dir <db_dir>`
      and the run goes on. `--set annotation.uniop.enabled=true` makes a missing UniOP an error
      (exit 2) instead.
- **What the distance model loses.** Without UniOP, the operon pairs and the TUs come from the
  intergenic-distance model of [ADR-0022](../adr/ADR-0022-uniop-terminatornet-cutoffs.md)
  (linked when d ≤ 102 bp). On the RegulonDB K-12 pairs it has precision 0.828 [0.805, 0.847]
  and recall 0.916 [0.901, 0.931], against UniOP's precision 0.90 at its cut-off. The run carries
  the flag `operon_source_distance_model` and the manifest `operon_source: distance_model`. Per
  [ADR-0061](../adr/ADR-0061-distance-model-fallback-scope.md), Moran's I rows under the TU-block
  null and the discordance `both_detected_t` rows carry `tu_source_distance_model`: Moran's I is
  anti-conservative at short lags (rejection 0.34 at α = 0.05 in simulation), and discordance is
  conservative (fewer pairs are found). Under `--publication`, `auto` with UniOP missing is an
  error; set `annotation.uniop.enabled: false` to publish with the distance model.
- `gippa compare` warns `compare.operon_source` when two runs used different operon sources.

## How long each step takes

Single measurements, not averages. Unless the row says otherwise, they were taken on 2026-11-03
in one Slurm job with 4 CPUs (Intel Xeon Gold 6148) and 16 GB, Apptainer 1.4.4 and Python 3.12,
writing to a network filesystem. Downloads depend on your network and the source servers; image
pulls also depend on CPUs (Apptainer converts every image with `mksquashfs`).

| Step | Time | Notes |
|---|---|---|
| `pip install --prefer-binary gippa` | 3 min 35 s | fresh venv, no pip cache |
| `gippa demo` | 1 min 55 s | 92 s in a second 4-CPU job on the same day |
| `gippa preflight` | 2 s; 10 s with `--profile` | demo inputs (100 kb genome, 1,685 peptides) |
| `gippa setup`: the 11 pinned images | 58 min | 8.7 GB of images; the antiSMASH image alone took 11 min 30 s. With 18 CPUs on another node (2026-10-30): 7 min 15 s |
| `gippa setup`: UniOP | 4 min 56 s | 13 s to fetch the source, 4 min 43 s to build |
| `gippa setup`: Bakta light database | not measured | 1.3 GB download |
| `gippa setup`: geNOMAD database | about 3 min | the 803 MiB archive at 4.5 MiB/s (2026-10-14, 4 CPUs); unpacking not timed |
| `gippa setup`: antiSMASH databases | 9 min 46 s | antiSMASH's downloader (2026-10-14, 4 CPUs) |
| `gippa setup`: Promotech model | 1 min 46 s | from the GitHub LFS copy |
| `gippa setup`: GO, KEGG, COG, contaminants | 11 s | all four together |
| `gippa setup`: Nextflow launcher | 2 s | |
| `gippa setup --bakta full` (not the default) | 3 h 04 min | about 84 GB: download and unpack (2026-10-14, 4 CPUs) |
| `gippa run`, E. coli K-12, all modules | 3 h 41 min | mean of 2 runs, 16 threads and 64 GB ([Runtimes](../user-guide/runtimes.md)) |

In total: about 5 min to install and check; for a default `gippa setup` on 4 CPUs, 1 h 20 min
for the measured rows plus the Bakta light download (not measured), so plan for about 1.5 h; and
3–4 h per run.

## Execution profiles

`gippa run --profile` takes one or more of these Nextflow profiles, joined by commas: `slurm`,
`apptainer`, `conda`, `docker`, `standard`, at most one executor and one container runtime. The
default is `auto`, which prints its choice:

- with `sbatch` on `PATH`: `slurm,<runtime>`, where the runtime is the one `gippa setup --images`
  recorded in `db_manifest.json` when it is available here, else Apptainer (or Singularity), else
  conda; with only Docker available it chooses `docker`, and the tasks run on the launching host;
- without `sbatch`: the recorded runtime when available, else Docker, Apptainer (or Singularity)
  or conda, in that order;
- with no container runtime and no conda it stops and asks you to install Apptainer or Docker.

The `farmshare` and `farmshare_local` profiles are site profiles used by the developers (see the
[developer guide](../developer.md)); they are not meant for users.

| Profile | What it does | Status |
|---|---|---|
| `apptainer` | runs every pinned tool from its image in the Apptainer cache; binds the DB directory read-only | tested (Linux x86-64, with Slurm) |
| `docker` | the same images through Docker: DB directory bound read-only, tasks run as your user (`-u $(id -u):$(id -g)`); never pulls (`docker run --pull=never`): images come from `gippa setup --images --profile docker` | smoke-tested (`docker-smoke`); see [Platforms](../user-guide/platforms.md) |
| `conda` | lets Nextflow create each tool's environment from `workflow/envs/*.yml`; preflight still requires the prebuilt environments of the conda-only tools | not tested on its own |
| `farmshare` | developers' site profile ([developer guide](../developer.md)) | not for users |
| `farmshare_local` | developers' site profile ([developer guide](../developer.md)) | not for users |
| `slurm` | generic Slurm executor: the cluster's default partition, no site-specific options; the workflow reads the queue, account, extra `sbatch` options, job cap and time limit from the `params.slurm_*` options, but no `gippa run` flag sets them yet in 1.0.0 (`slurm_queue`, `slurm_account`, `slurm_options`, `slurm_max_jobs`, `slurm_max_time`); combine it with one container runtime, e.g. `slurm,apptainer` | not tested |
| `standard` | Nextflow's local executor, with no container or Conda settings | not tested |

### Run on a cluster

The simplest route is one batch job: request the CPUs and memory for the whole run (for example
16 CPUs and 64 GB) and run `gippa run --profile apptainer` inside it. Every task then runs inside
that job's allocation, and nothing else has to be set up on the other nodes. Do not start
`gippa run` on a login node: the Nextflow head process runs on the host where you start it, for
the whole run.

`--profile slurm,apptainer` (which `auto` chooses where `sbatch` is on `PATH`) instead submits
each task as its own Slurm job, on the cluster's default partition. It is not tested yet (table
above), and it needs more set-up:

- the analysis steps (`GIPPA_STEP`) run `gippa` on the host, not in a container, so every node
  that runs tasks needs the same Python environment with the same `gippa` version on `PATH`, for
  example one venv on a shared filesystem;
- the DB directory and the Apptainer cache (`NXF_APPTAINER_CACHEDIR`) must be on a filesystem
  that every node sees;
- 1.0.0 has no command-line flag to choose a partition or account. `gippa run` has no `-c`
  option for your own Nextflow config, but Nextflow also reads `$NXF_HOME/config` (default
  `~/.nextflow/config`), where `process.queue` and `process.clusterOptions` can be set.

## What `gippa run` needs before preflight passes

In run mode, and in `gippa preflight --profile`, preflight checks that every tool the run will use
is installed under `--profile` (`missing_tools` in `src/gippa/preflight/environment.py`). A
missing tool is an environment error: `gippa run` exits with code 2 before Nextflow starts, and
the message names the tool and the reason. Only `gippa setup` pulls images: `gippa run` and
`gippa preflight` use the images already on the machine, or stop and name the `setup` command
that fetches them.

| Tool kind | `docker` | `apptainer` (no conda profile) | with `conda` (or a site profile) |
|---|---|---|---|
| pinned in `containers.config` | present iff `docker image inspect <ref>` succeeds (runs never pull); Docker < 20.10 is refused | present iff its `.img` is in the cache | with `apptainer` in the list: as in that column; with `docker` in the list (no `apptainer`): as in the `docker` column; else not checked |
| no image pin (conda-only) | missing: `conda-only tool under a container profile without conda` | same | needs the prebuilt env `$GIPPA_ROOT/micromamba/envs/gippa-<tool>/` |
| UniOP | the local build in `local_tools.config` (image ID checked) | the local `.sif` (sha256 checked) | site profiles: the prebuilt env; plain `conda`: the micromamba build |
| fgsea (`function` with `--abundance`) | prebuilt env, else its `images.json` image | prebuilt env, else that image in the cache | prebuilt env |
| in-process tools | never missing | never missing | never missing |

A prebuilt env counts only with its marker file `.gippa-built`, holding the SHA-256 of
`workflow/envs/<tool>.yml`. Preflight also reads `DIR/db_manifest.json` for the databases
([Set up the databases](setup-databases.md)). The antiSMASH image alone is 4.9 GB, and all 11
images 8.7 GB ([Set up the databases › Layout and sizes](setup-databases.md#layout-and-sizes)).

## What `gippa setup` needs

`gippa setup` fetches most items itself over HTTPS (the Nextflow launcher too, when none is on
`PATH`) and needs no Nextflow, no Java and no checkout. Two
downloaders run inside a tool environment or image: `bakta_db download` (which also runs
`amrfinder_update`) and antiSMASH's database downloader (`--antismash`). For each,
`ToolRunner` (`src/gippa/setupdb/run.py`) takes the first of:

1. a micromamba environment `gippa-bakta` or `gippa-antismash` under `$MAMBA_ROOT_PREFIX/envs/`;
2. the pinned image in the Apptainer cache (as in step 3), run with `apptainer exec`;
3. the pinned image through `docker run -u <uid>:<gid>`, when `docker` is on `PATH`.

When none is found, the item stays `absent` with the message `no gippa-<tool> environment, no
pinned <tool> image … run gippa setup --images --profile apptainer|docker`, and setup exits 2.

## Licences

GIPPA is MIT-licensed; the tools keep their own licences ([Licences](../user-guide/licences.md)).
geNOMAD (module `mobile`) is for non-commercial use only: for commercial use, run with
`--modules core,bgc,function` and install with `--no-genomad`. KEGG (`--kegg`) is free for
academic use only. The GPL and AGPL tools (Bakta, Promotech, Promoter Calculator, OSTIR,
TransTermHP, BacTermFinder; antiSMASH) run as separate programs in their own images or
environments; GIPPA never imports or vendors their code.

## Check the installation

```bash
gippa --version
gippa setup --db-dir /path/to/db --check          # verify what is installed, download nothing
```

The [CLI reference](../reference/cli.md) lists every command and option. To move to a newer
release, see [Upgrade](upgrade.md).
