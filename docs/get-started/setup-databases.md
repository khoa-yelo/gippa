# Set up the databases

`gippa setup` downloads the databases and models into one directory and records each item in
`DIR/db_manifest.json`. Run it once per installation. Then pass the directory to `gippa run` with
`--db-dir DIR`, or set `GIPPA_DB_DIR=DIR`. Without internet access, or to reuse a shared copy,
see [Download the databases by hand](manual-databases.md).

## The command for the default modules

`gippa run` runs all four modules by default: `core`, `mobile` (geNOMAD), `bgc` (antiSMASH) and
`function` (COG, GO, KEGG). `gippa setup` installs everything they need by default:

```bash
gippa setup --db-dir /path/to/db --profile apptainer \
  --accept-download-plan --install-uniop-from-git      # or --profile docker
```

The last two flags answer setup's two questions in advance, so it also runs as a batch job:
`--accept-download-plan` the "Proceed?" after the download plan, and `--install-uniop-from-git` the UniOP build
(UniOP has no licence file; the flag means you decide your use is permitted). Leave them out to be
asked ([Install › Step 3](install.md#the-five-steps)).

That is the Bakta light database, the geNOMAD database, the antiSMASH databases, the Promotech
models, GO, KEGG, COG, the contaminant library and, with `--profile`, the container images. Each
item has an opt-out: `--no-genomad`, `--no-antismash`, `--no-go`, `--no-kegg`, `--no-cog`,
`--no-contaminants`, `--no-images`. `gippa setup` is the only step that pulls images: `gippa run`
never pulls one (under Docker its tasks run with `docker run --pull=never`) and stops when one is
missing.

`--minimal` is for small disks: it leaves out geNOMAD and antiSMASH and pulls only the images of
the `core` and `function` modules. It is recorded as `setup.set: minimal` in `db_manifest.json`;
a later `gippa run` on that directory without `--modules` runs `core,function` and says so.

Before any download, setup prints the licence terms of the planned items:

- `geNOMAD (module mobile) is for non-commercial use only; for commercial use: gippa setup
  --no-genomad, gippa run --modules core,bgc,function`;
- `KEGG is free for academic use only; for other use: gippa setup --no-kegg`.

Each run records the terms of what it used in `run_manifest.json` (`licence_notes`), and the
report's Methods section shows them. See [Licences](../user-guide/licences.md).

To use the Bakta full database, add `--bakta full`. Runs annotated with the light database carry the
run flag `bakta_db_light`. With `gippa run --annotation bakta.json`, Bakta is not run. Preflight
still checks that a Bakta database is installed, though, so do not use `--bakta skip` for a run.

Before the first `gippa setup`, read [Install › What `gippa setup` needs](install.md#what-gippa-setup-needs):
the Bakta and antiSMASH downloaders run in a tool environment or image that `gippa setup` does not
install.

On a shared cluster, run `gippa setup` as a batch job: `--bakta full` alone downloads about 84 GB.

## Flags

Every flag of `gippa setup`, in plain words. The [CLI reference](../reference/cli.md#gippa-setup)
has the same flags in short form. Sizes are the unpacked files (a compressing filesystem stores
less). Times are single measurements on one cluster, mostly with 4 CPUs (conditions in
[Install › How long each step takes](install.md#how-long-each-step-takes)); they depend on your
network and CPUs.

"Preflight stops" means: `gippa run` (and `gippa preflight --profile`) exits with code 2 before
Nextflow starts and names what is missing. Every module's databases are checked, so when you leave
an item out, also leave out the module that needs it with `gippa run --modules`.

### Where and how

| Flag | Default | What it does | Turn it off or change it when |
|---|---|---|---|
| `--db-dir DIR` | `$GIPPA_DB_DIR`, else the directory the last setup recorded, else `$XDG_DATA_HOME/gippa/db` or `~/.local/share/gippa/db` | where everything is installed; recorded, so `gippa run`, `preflight` and `makedb` find it without the flag | always set it to a filesystem with about 50 GB free that your compute nodes can read |
| `--profile` | `auto` | the container runtime for the images and the UniOP build: `apptainer` or `docker` (a comma-joined Nextflow profile is accepted); `auto` detects one and prints it | set it when `auto` picks the wrong runtime; it must match the `--profile` of your runs |
| `--images [all]` / `--no-images` | on when `--profile` names a runtime | pulls the pinned container images the configured modules need (`all`: every pinned image). `gippa run` never pulls: a missing image stops preflight. Size and time: 11 images, 8.7 GB; 58 min on 4 CPUs, 7 min 15 s on 18 CPUs | `--no-images` when the images are already in the cache, or you run with the `conda` profile |
| `--cache-dir DIR` | `$NXF_APPTAINER_CACHEDIR`, else `DIR/containers` | the Apptainer image cache and UniOP build directory, for this setup command only | rarely: set `NXF_APPTAINER_CACHEDIR` instead, so preflight and run see the same cache ([Install › Step 3](install.md#the-five-steps)) |
| `--nextflow` / `--no-nextflow` | on | installs the pinned Nextflow 26.10.0 launcher as `DIR/bin/nextflow`, only when no `nextflow` is on `PATH`. 47 MiB; 2 s | `--no-nextflow` when you install Nextflow yourself |
| `--accept-download-plan` | off | answers the "Proceed?" question after the download plan with yes; does not accept any licence terms; in a batch job (no terminal) it is never asked | use it in scripts and batch jobs |
| `--check` | off | downloads nothing and asks nothing: verifies what is installed and prints each item's status | after setup, or to find what is missing |
| `--config FILE` / `--set key=value` | the defaults | the run configuration, as for `gippa run`; setup reads `modules` from it to choose the images and whether UniOP is needed | when you will run fewer modules, for example `--set modules=core,function` |
| `--log-level` | `INFO` | how much setup prints | for debugging |

### What it downloads

| Flag | Default | Downloads (size; time) | Terms | Needed by; if you leave it out | Turn it off when |
|---|---|---|---|---|---|
| `--bakta` `light`, `full` or `skip` | `light` | Bakta database v6.0 (Zenodo record 14916843) for genome annotation; `bakta_db download` also runs `amrfinder_update`. Light: 4.2 GB (1.3 GB archive), time not measured. Full: about 84 GB, 3 h 04 min | see the Zenodo record | module `core`, always: preflight requires a Bakta database even with `--annotation bakta.json`, so `skip` stops preflight. The full database assigns functions to more genes; runs with the light one carry the flag `bakta_db_light`, and the report says the light database assigns fewer functions | use `full` when gene functions matter (the `function` module's GO, KEGG and COG sets); keep `light` for small disks and quick tests |
| `--genomad` / `--no-genomad` | on | geNOMAD database v1.9 (Zenodo record 14886553): 1.5 GB; the 803 MiB download took about 3 min (unpacking not timed) | **non-commercial use only** | module `mobile` (prophages and plasmids). Without it, run `--modules core,bgc,function` | commercial use, or small disks (`--minimal`) |
| `--antismash` / `--no-antismash` | on | antiSMASH 8.0.4 databases, through antiSMASH's own downloader in its image (skipped when the image already holds them): 10.1 GB, 9 min 46 s | antiSMASH is AGPL-3.0 and runs as a separate program | module `bgc` (biosynthetic gene clusters). Without it, run `--modules core,mobile,function` | small disks (`--minimal`), or no interest in gene clusters |
| `--promotech-models` / `--no-promotech-models` | on | the Promotech RF-HOT promoter model (`www.cs.mun.ca`, else the authors' GitHub LFS copy): 4.8 GB (the 1.0 GB archive and the 3.8 GB model in it); 1 min 46 s | Promotech is GPL-3.0 and runs as a separate program | module `core` (promoter predictions), always: without it preflight stops | only when you put the model in place by hand ([Download the databases by hand](manual-databases.md)) |
| `--go` / `--no-go` | on | `go-basic.obo` and `goslim_generic.obo` (release 2026-05-19): 32 MB; 11 s together with KEGG, COG and the contaminants | CC BY 4.0 | module `function` (GO enrichment). Without it, run `--modules core,mobile,bgc` | when you do not run `function` |
| `--kegg` / `--no-kegg` | on | KEGG KO-to-pathway links and pathway names from `rest.kegg.jp` (unversioned: the retrieval date and SHA-256 are recorded): 2.3 MB; see GO for the time | **free for academic use only**; others need a KEGG licence | module `function`, optional: without KEGG, preflight warns (`environment.kegg`) and `function` runs GO and COG and skips only the KEGG pathway analyses; the report says "KEGG not installed: pathway enrichment skipped" | non-academic use without a KEGG licence |
| `--cog` / `--no-cog` | on | the COG 2020 category tables `cog-20.def.tab` and `fun-20.tab` (NCBI FTP): 0.4 MB; see GO for the time | public (NCBI FTP) | module `function` and the `cog_category` column of `genes.tsv`. Without it, run `--modules core,mobile,bgc` | when you do not run `function` |
| `--contaminants` / `--no-contaminants` | on | the universal protein contaminant FASTA (commit `bddf6da9`): 0.2 MB; see GO for the time | no licence file: downloaded from its authors, never redistributed | contaminant labels on peptides in mapping and novel-ORF calls, and `gippa makedb`. Preflight does not require it: without it, peptides from common contaminants are not labelled | rarely |
| `--minimal` | off | everything above except geNOMAD and antiSMASH, and only the images of the `core` and `function` modules | as above | a later `gippa run` on this directory without `--modules` runs `core,function` and says so | small disks |

### UniOP

UniOP, the operon predictor, has no licence file, so GIPPA never redistributes it: setup fetches
its source from the authors (pinned commit) and builds it on your machine ([Install › UniOP](install.md#uniop)).

| Flag | Default | What it does | If you leave UniOP out |
|---|---|---|---|
| `--install-uniop-from-git` | off: setup asks `Build UniOP for your own use? [y/N]` on a terminal, and skips UniOP in a batch job | answers that question with yes: **you decide that your use of UniOP is permitted**; UniOP is fetched from its authors' GitHub repository (hongsua/UniOP) at the pinned commit, sha256-checked and built during the normal setup. Fetch and build: 4 min 56 s; the Apptainer image is 0.3 GB | operons and transcription units come from the intergenic-distance model; the run carries the flag `operon_source_distance_model` |
| `--no-uniop` | off | does not build UniOP under `--images`, and asks nothing | as above |
| `--uniop` | off | installs only UniOP (no databases, no other images); add `--install-uniop-from-git` to skip the question | — |

## Layout and sizes

The paths are relative to `--db-dir`, from `src/gippa/data/db_sources.yaml`. Sizes are the
unpacked files of an installed copy (measured 2026-10-15 and 2026-11-03); a compressing filesystem
stores less. Setup checks free space against the larger `size_gib` reserves of `db_sources.yaml`.

| Item | Path | Size |
|---|---|---|
| Bakta light DB v6.0 | `bakta/db-light` | 4.2 GB (the 1.3 GB archive can be deleted) |
| Bakta full DB v6.0 | `bakta/db` | about 84 GB |
| geNOMAD DB v1.9 | `genomad_db` | 1.5 GB (803 MiB archive) |
| antiSMASH databases | `antismash` | 10.1 GB |
| Promotech RF-HOT | `promotech` | 4.8 GB (`RF-HOT.zip` 1.0 GB and the 3.8 GB `RF-HOT.model`) |
| GO | `go/go-basic.obo`, `go/goslim_generic.obo` | 32 MB |
| KEGG | `kegg/ko_pathway.tsv`, `kegg/pathway_names.tsv` | 2.3 MB |
| COG 2020 | `cog/cog-20.def.tab`, `cog/fun-20.tab` | 0.4 MB |
| Contaminants | `contaminants/contaminants.fasta` | 0.2 MB |
| Nextflow launcher 26.10.0 | `bin/nextflow` | 47 MiB (49,084,616 bytes) |

The Apptainer images take a separate 8.7 GB in the image cache for the 11 pinned images, plus
0.3 GB for the UniOP build (measured 2026-11-03).

## `db_manifest.json`

`gippa setup` writes `DIR/db_manifest.json`, next to the items. Each item gets a record with its
name, version, URL, checksum, size and retrieval date. The schema is
`src/gippa/schemas/db_manifest.schema.json`. `gippa run` copies the database versions and checksums
into the run's `run_manifest.json`. Preflight reads the manifest to check that the database each
module needs is installed.

KEGG has no release version. For KEGG, the manifest records the retrieval date and the SHA-256 of
the downloaded bytes ([ADR-0017](../adr/ADR-0017-unversioned-sources.md)). The contaminant
library is pinned: commit `bddf6da9`, with a fixed SHA-256 in `db_sources.yaml`.

## Failures

`gippa setup` exits with code 2 in two cases: a source cannot be reached, or a checksum does not
match. The message names the host, and the item stays absent. GIPPA never falls back to an
unofficial source. One item has a second official source: when `www.cs.mun.ca` cannot be reached,
the Promotech RF-HOT model is fetched from the Promotech authors' GitHub LFS copy on
`media.githubusercontent.com`, checked against the same SHA-256, and the manifest records the URL
that served it ([ADR-0069](../adr/ADR-0069-promotech-model-github-lfs-fallback.md)). A checksum
mismatch is final: setup never tries another source after one. When both Promotech sources fail,
the message names both hosts and the directory to fill by hand: put `RF-HOT.zip` and its unzipped
`RF-HOT.model` in `DIR/promotech/`, then rerun `gippa setup`.

Run the command again later; `--check` shows which items are present. To download an item
yourself and have setup adopt it, see [Download the databases by hand](manual-databases.md).

## What `gippa setup` does not download

- **RegulonDB.** Users never download RegulonDB. GIPPA's developers use RegulonDB v12 only to
  calibrate cut-offs on E. coli K-12, and GIPPA does not redistribute it.
- **Search-engine software and spectra.** GIPPA never searches spectra. You run DIA-NN or a de novo
  tool yourself, under your own licence (see [Licences](../user-guide/licences.md)).
- **Genomes.** You provide the genome FASTA.

## Licences of the downloaded data

geNOMAD is for non-commercial use only. KEGG is free for academic use; other users need a licence
from KEGG. Read [Licences](../user-guide/licences.md) before you install either.
