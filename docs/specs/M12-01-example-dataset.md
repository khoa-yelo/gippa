# M12-01: a real, cut-down example dataset and `gippa example fetch`

Milestone M12 (1.1.0) · depends on: M10 pilot (PXD078459 arms A2 and A3, `validation/pilot/`),
M11-15 (`gippa demo`, `HttpsFetcher`, `release-assets` redirect rows), M11-04/M11-05
(`release.yml`, `docs/release.md`), M11-07/M11-08 (site, nav and docs tests) · estimate ≈ 300 code
lines (CLI + loader + reference generator) + ≈ 450 test lines + the build script (owned by
validation, branch `feat/example-dataset`) + docs, split into nine PRs (table at the end) · roles:
spec-writer (this file), lead (ADR-0074, routing), validation (`scripts/example/`,
`scripts/slurm/build_example.sbatch`, the measurement), test-author (`tests/unit/test_m12_01_*.py`,
`tests/e2e/test_m12_01_example_routes.py`, the nav and docs tests), implementer
(`src/gippa/example.py`, `src/gippa/cli.py`, `src/gippa/data/`, `scripts/docs/gen_reference.py`),
docs (`docs/get-started/`, `docs/tutorials/`, `docs/release.md`) · source: Khoa, 2026-11-03 (in
session).

## Source: Khoa's decisions, 2026-11-03 (in session)

> "we should provide a real but subset dataset that user can run on and get result quick, say
> within 5 mins; this example should show different routes of using the software, through diann
> output, through peptides faa and abundance tsv, through denovo output etc..."

Further decisions in the same session (Khoa, 2026-11-03):

- The data is a subset of the M10 pilot: E. coli Nissle 1917, PXD078459. The pilot's "no
  row-level DIA-NN content" rule (Khoa 2026-10-23, `validation/pilot/README.md` › Conversion) is
  **lifted for this subset only**. Everything else of the pilot stays aggregate-only.
- Real tools only: no precomputed annotation layers. The data is a small genome region.
- Hosting is a GitHub release asset: not inside the wheel, and not on Zenodo.

## Spec sections touched

- §5.1 command line and exit codes (a new command group `gippa example`).
- §3.2 download sources: "GitHub releases" is already a pinned source; the bundle is GIPPA's own
  release asset. R6 applies: pinned sha256, official host only, no mirror, TLS on.
- R1 scope: the example is not in §2 and not excluded by §16. It is a new user-facing command, so
  **ADR-0074** records it (R10, R13): the command, the hosting, the lifted row-level rule and the
  data licence statement. Khoa accepts the ADR; the spec does not change.
- R11: the bundle is never committed. Nothing row-level from it goes into the repository, the
  docs site or `tests/data/`. R12: the bundle holds no personal paths (see "Scrubbing").
- §9.2 / ADR-0060 (semantic versioning): a patch release gives byte-identical outputs, a minor may
  change numbers. This decides the bundle versioning below.
- §6.5: coordinates are 0-based half-open internally. Human-facing files of the bundle state both
  conventions explicitly.

## What already exists (reuse, do not duplicate)

| Piece | Where | Reused for |
|---|---|---|
| HTTPS download, redirect allowlist, body cap, `SourceError` | `gippa.setupdb.run.HttpsFetcher` | the bundle download, unchanged |
| per-row redirect hosts and size cap | `db_sources.yaml` `redirect_hosts`, `size_gib`; `_redirect_hosts`, `_body_cap` look the URL up in `SOURCES` | the bundle row (see "Pin") |
| GitHub release asset precedent | `db_sources.yaml › nextflow` (`redirect_hosts: [release-assets.githubusercontent.com]`) | same host, same redirect |
| sha256 of a file | `gippa.setupdb.items.file_digest` | archive and per-file checks |
| safe unpack (no absolute names, no `..`, links only inside, `filter="data"`) | `gippa.setupdb.run._unpack` | the tarball. The implementer may rename it public (`unpack_archive`) and keep `_unpack` as an alias; no copy |
| "exists → exit 1", `_first_missing`, `_out_path` | `gippa.cli.demo` | outdir handling |
| exit codes | `gippa.exitcodes.ExitCode` | unchanged |

`gippa setup` builds its plan from named rows (`items.plan`), so an extra row in
`db_sources.yaml` is never installed by setup. `test_setupdb.py`'s invariant over every row (https,
host, sha256 or md5, relative `path`, `tool` null or known) already holds for a bundle row.

## 1. `gippa example fetch`

```
gippa example fetch [--outdir DIR] [--version X.Y] [--from-file FILE]
```

`gippa example` is a Click group with one subcommand, `fetch`. `gippa example` alone prints the
group help (exit 0). The shared options (`--config`, `--set`, `--log-level`) are accepted and
ignored, as for `demo`.

| Option | Default | Meaning |
|---|---|---|
| `--outdir DIR` | `gippa-example` | where the bundle is unpacked; must be absent or empty (see below) |
| `--version X.Y` | the installed GIPPA's `MAJOR.MINOR` (selection rule below) | which bundle. `X.Y.Z` is accepted and read as `X.Y` |
| `--from-file FILE` | none | offline: verify and unpack a bundle file already on disk, with no network access |

`--version` is the bundle version, not `gippa --version`; the subcommand option shadows nothing
because the root `--version` is only on the `gippa` group. The help text says "bundle version".

### Steps

1. **Select the pin.** Read the `example_<X>_<Y>` rows of `gippa/data/db_sources.yaml` (one per
   bundle; "Pin" below). Selection without `--version`: the row whose version equals the
   installed `MAJOR.MINOR`; else the newest row with the same major and a lower minor, with the
   warning `gippa example: no bundle for <installed>; using <X.Y> (expected numbers were made
   with GIPPA <X.Y> and may differ, see Upgrade)`; else exit 2 `no example bundle for GIPPA
   <MAJOR>.x`. A dev or rc version (`1.1.0.dev3`, `1.1.0rc1`) counts as its `MAJOR.MINOR`.
   With `--version`: that row, else exit 1 listing the pinned versions. A pin older than the
   installed minor gives the same warning.
2. **Check the outdir.** Absent, or an empty directory: continue. Holding a previous fetch of the
   **same** bundle (its `bundle.json` `sha256` equals the pin and every file verifies against
   `CHECKSUMS.sha256`): print `already fetched`, print the commands, exit 0, no network. Anything
   else: exit 1 `<DIR> exists and is not this example; use --outdir or remove it`.
3. **Get the archive** into a private staging directory `DIR.partial-<random>` created next to
   `DIR` (same filesystem, so the final rename is atomic).
   - Default: `HttpsFetcher().fetch(row.url, staging/<file name>)`.
   - `--from-file FILE`: copy nothing; hash `FILE` in place. The fetcher is never constructed.
4. **Verify the archive.** Size equal to the pin's `size_bytes` (a cheap, clear first message), then
   `file_digest(archive, "sha256") == row.checksum`.
5. **Unpack** with `_unpack` into the staging directory. The archive must hold exactly one
   top-level directory, `gippa-example-<X.Y>/`; anything else is a mismatch.
6. **Verify the content.** Every line of `CHECKSUMS.sha256` (`sha256sum` format, relative paths)
   matches its file, and no regular file exists that the list does not name (except
   `CHECKSUMS.sha256` itself).
7. **Move into place.** Rename `gippa-example-<X.Y>/` to `DIR` (or move its entries into an
   existing empty `DIR`), remove the staging directory and the downloaded archive (the
   `--from-file` file is never touched).
8. **Print** to stdout: one line `gippa example <X.Y>: <DIR>` and then, per route, a heading and
   its commands, with `DIR` resolved to an absolute path. The command text comes from
   `gippa/data/example_routes.yaml` (below), never from the downloaded bundle. A last line says
   that the routes need `gippa setup` (at least `--minimal`) and a container runtime.
   Progress and warnings go to stderr, prefixed `gippa example:`.

On every failure after step 3 the staging directory is removed. A failed fetch leaves no `DIR`
that it created and no partial file (as `HttpsFetcher.fetch` already does for its own file).

### Exit codes

| Code | When |
|---|---|
| 0 | fetched and verified; or already fetched (same bundle) |
| 1 | invalid input: unknown `--version`; `DIR` not empty and not this bundle; `DIR` cannot be created; `--from-file` missing, not a file, or not the pinned bundle (wrong size or sha256) |
| 2 | environment: no bundle for this GIPPA major; network error or HTTP error (`SourceError`), including 404; the downloaded archive's size or sha256 differs from the pin; unsafe or malformed archive; `CHECKSUMS.sha256` mismatch; disk full or quota (`DISK_FULL` errnos) |
| 5 | an unexpected exception (the standard internal-error message) |

A wrong `--from-file` is the user's input (1); a wrong download is the source's fault (2), the
same as `gippa setup`'s `mismatch`. When a `--from-file` archive matches **another** pin, the
message says so: `this is example bundle 1.2; use --version 1.2`.

### No network

`HttpsFetcher` raises `SourceError(host, reason)` for DNS, connect, TLS, timeout, HTTP status and
truncated bodies. The command prints, on stderr:

```
gippa example: cannot download <url> (<host>: <reason>)
gippa example: on a machine with internet access, download that URL, copy the file here and run
  gippa example fetch --from-file gippa-example-<X.Y>.tar.gz --outdir <DIR>
```

and exits 2. It does not retry (R6's "one later retry" is the user's rerun) and never tries
another host. A 404 adds: `the bundle may not be published yet (GIPPA <X.Y> release)`.

### Pin (version pinning)

- **One bundle per GIPPA minor.** Its version is `X.Y`. It is built at the release-prep commit of
  `X.Y.0` and attached to the GitHub release `vX.Y.0` as `gippa-example-X.Y.tar.gz`. Patch releases
  `X.Y.Z` reuse it: by ADR-0060 a patch keeps the goldens byte-identical, so `expected/` still
  holds; a patch that changes example numbers is a minor by that same rule. A minor whose inputs
  are unchanged still gets a new bundle, because `expected/` changes with the minor.
- The pin is a row in `src/gippa/data/db_sources.yaml`, so `HttpsFetcher` finds its redirect hosts
  and size cap without any change:

  ```yaml
  example_1_1:                    # M12-01: the example bundle of GIPPA 1.1 (ADR-0074); never installed by setup
    version: "1.1"
    url: https://github.com/khoa-yelo/gippa/releases/download/v1.1.0/gippa-example-1.1.tar.gz
    checksum: sha256:<64 hex>     # validation/example/bundles.tsv; build job <id>
    size_bytes: <n>
    path: example/gippa-example-1.1.tar.gz   # unused by setup; must be safe and relative
    tool: null
    check: sha256
    size_gib: 1
    redirect_hosts: [release-assets.githubusercontent.com]
  ```

  `size_bytes` is a new optional key of the loader (`DbItem.size_bytes: int = 0`, ignored by
  setup). The package keeps every older bundle row it can still run (same major), so
  `--version 1.1` works from GIPPA 1.2.
- `src/gippa/data/example_routes.yaml` holds the route commands with the placeholders `{dir}`
  (the absolute outdir) and `{out}` (`{dir}/runs`). It is the single source for the printed
  commands, the bundle README and the docs page (tests below keep them equal).

## 2. The bundle

### Layout

```
gippa-example-<X.Y>/
├── README.md              what this is, the routes and their commands, the expected numbers, caveats
├── PROVENANCE.md          sources, sha256s, region, subset rules, build commit and Slurm job
├── LICENSE-DATA.md        the licence statement (below)
├── CHECKSUMS.sha256       sha256 of every other file, sorted by path
├── bundle.json            bundle version, GIPPA version and commit used for expected/, region, counts, tool and DB versions
├── genome/
│   └── region.fna.gz      the Nissle region (and, if kept, one whole plasmid)
├── A_diann/
│   └── report.parquet     DIA-NN 2.3.1 report rows of the region (targets and decoys)
├── B_proteins/
│   ├── proteins.faa       detected region proteins (sequences from validation/pilot/db_a2)
│   └── abundance.tsv      protein-group MaxLFQ of those proteins
├── C_cascadia/
│   ├── WT_1.ssl
│   ├── WT_2.ssl
│   └── WT_3.ssl           Cascadia 0.0.7 rows of the region
├── D_generic/
│   └── peptides.tsv       a generic peptide table as a user of another engine would write it
└── expected/
    ├── A.json  B.json  C.json  D.json
    └── README.md          the same numbers as a table, with the file and column each comes from
```

The archive is a deterministic tar.gz: members sorted, mode 0644/0755, uid/gid 0, empty
user/group names, mtime 0, gzip header without name and with mtime 0. Two builds from the same
inputs give the same sha256.

### Routes (one input set each, same region)

| Route | Input | Commands (about three) | Expected tier (to be measured) |
|---|---|---|---|
| A | `A_diann/report.parquet` | `gippa convert --from diann … --out {out}/A/peptides.tsv --search-db same`; `gippa preflight`; `gippa run --peptides …` | `decoy_picked` if ≥ 80 % of the region decoys pair, else `mayu` |
| B | `B_proteins/proteins.faa` + `abundance.tsv` | `gippa preflight --proteins … --abundance …`; `gippa run --proteins … --abundance … --abundance-type maxlfq` | `protein_list` |
| C | `C_cascadia/*.ssl` | `gippa convert --from cascadia … --out {out}/C/peptides.tsv`; `gippa preflight`; `gippa run --peptides …` | `two_peptide`; `c_assumed` (no reference peptides) |
| D | `D_generic/peptides.tsv` | `gippa convert --from generic … --out {out}/D/peptides.tsv --source database --search-db same`; `gippa preflight`; `gippa run --peptides …` | `mayu` (q-values, no decoys) |

Every `gippa run` uses `--genome {dir}/genome/region.fna.gz --genome-relation same --modules core
--outdir {out}/<route>` and the thread and memory values chosen by the measurement (default
`--threads 4 --memory 16.GB`). No route uses `--annotation`, `--strict` or `--resume`: Bakta and
every other core tool run live (Khoa: real tools only).

**Route D's table.** Built from the route-A subset, as an engine without a GIPPA converter would
write it: columns `peptide`, `modified_peptide` (mass-bracket notation such as `M[+15.995]`, so
the generic converter's stripping is exercised), `protein_ids`, `q_value`, `intensity`, `run`; no
decoy rows and no `is_decoy` column. It therefore shows the column-to-tier table of
`from-spectra.md` route D: q-values give `mayu`, not `decoy_picked`.

### Subset rules (requirements on the build script)

The build script (below) implements these; its internals are not specified here.

1. **Region.** One contiguous slice of chromosome CP171242.1 (GCA_043228065.1), not crossing its
   origin. Both ends fall in intergenic gaps of the pilot's Bakta annotation (no gene cut), so no
   gene is partial by construction. Written as a linear contig (no topology tag), ID
   `CP171242.1_<start0>_<end0>` (0-based half-open, matches `[A-Za-z0-9._-]+`), description
   `region <start1>-<end1> of CP171242.1 (1-based closed)`. Optionally (measurement decides) one
   whole Nissle plasmid with `[topology=circular]`, to exercise a circular replicon. Region choice
   (size, place) is TBD until measured; when two candidates score equal, the lower start wins.
2. **Region proteins.** A pilot protein (`validation/pilot/db_a2/db_id_map.tsv`) belongs to the
   region when its CDS or sORF lies wholly inside it, on either strand.
3. **DIA-NN rows (A).** Keep a row of the A2 pass-1 report (`A2_report.parquet`, sha256
   `733a3884…`) when any ID in `Protein.Ids`, after stripping the decoy prefix, is a region
   protein. Targets and decoys use the **same** rule, so the target/decoy balance of the region is
   not biased. All columns are kept, apart from the scrubbing below.
4. **Protein table (B).** The region proteins in protein groups that pass q ≤ 0.01 in the A2
   report (target, non-contaminant), with sequences from `db_a2/strain_db.faa.gz`;
   `abundance.tsv` = those groups' rows of the A2 `--abundance-out` table (`maxlfq`, the median of
   `PG.MaxLFQ` over runs).
5. **Cascadia rows (C).** Rows of the three A3 `.ssl` files whose peptide the pilot's A3 run placed
   inside the region (`mapping/peptide_hits.tsv`, any accepted locus in the region), plus a seeded
   random sample of the rows that run left unplaced, at the region's share of the genome length.
   Without the second part every kept row would be placeable, and the example's placement rate
   would be inflated. The seed is fixed in the script and recorded in `PROVENANCE.md`.
6. **No empty route.** The build fails if any route keeps zero rows or zero proteins, or if route
   A keeps no decoy.
7. **Scrubbing (R12).** Columns that hold absolute paths (`File.Name` in the DIA-NN report, the
   file column of the `.ssl`) are rewritten to the base file name (for example
   `MLA1914689_WT_P1_1_Y-D7_1_13119.d`). No other value is changed. `PROVENANCE.md` lists each
   rewritten column.

### `PROVENANCE.md` (content)

- The PXD078459 record (URL, retrieval date, record sha256 `5d3f16c3…`), the three runs WT_1–WT_3
  with their file names, the instrument.
- The genome: GCA_043228065.1, NCBI md5 and the sha256 of `tests/data/genomes/nissle1917.fna.gz`;
  the region in both coordinate conventions; the plasmid, if kept.
- Route A: DIA-NN 2.3.1 Academia, the image sha256, `validation/pilot/diann_A2.cfg`, the input
  report's sha256, row counts kept (targets, decoys).
- Route B: the A2 convert table sha256s (`a2_abundance.tsv.gz`), protein counts.
- Route C: Cascadia 0.0.7, the weights sha256, the three `.ssl` sha256s, rows kept per run (placed,
  sampled), the seed.
- Route D: how it was derived from route A.
- The build: script path and commit, Slurm job ID, GIPPA version and commit that produced
  `expected/`, the Bakta DB (light 6.0) and image digests used, and that UniOP was not installed.

### Licence statement (`LICENSE-DATA.md`, default until Khoa answers Q-K1)

> The GIPPA software is MIT-licensed. This licence does **not** cover the data files of this
> bundle.
>
> The files under `A_diann/`, `B_proteins/`, `C_cascadia/` and `D_generic/` are derived from the
> ProteomeXchange dataset PXD078459 (iProX IPX0017169000). The dataset record states no licence.
> They are redistributed unchanged in content (rows selected; file paths shortened) for the sole
> purpose of demonstrating and testing GIPPA, with attribution to the depositors of PXD078459.
> No licence to them is granted by the GIPPA authors. Cite PXD078459 if you use them. If you are a
> rights holder and object, open an issue and the bundle will be withdrawn.
>
> `genome/` holds a region of GenBank GCA_043228065.1 (INSDC data, no use restrictions by NCBI).
> `expected/`, `README.md`, `PROVENANCE.md`, `bundle.json` and `CHECKSUMS.sha256` are written by
> GIPPA tooling and are MIT, as GIPPA.

### `expected/` (key numbers)

One JSON file per route, with a `tolerance` per key, made by running the four routes with the
release's GIPPA at the release-prep commit (Bakta light DB 6.0, no UniOP, default seed).

| Key | Source | Tolerance |
|---|---|---|
| `status` | `run_manifest.json › status` | exact |
| `detection_tier` | `run_manifest.json › detection_tier` | exact |
| `layers_missing` | `run_manifest.json › layers_missing` (as a set) | exact |
| `convert_rows` | data rows of the converted `peptides.tsv` (A, C, D) | exact (conversion is deterministic) |
| `n_genes` | rows of `tables/genes.tsv` | abs ≤ max(1, 1 %) |
| `n_detected` | `genes.tsv` rows with `detected = true` | abs ≤ max(1, 2 %) |
| `n_detected_2pep` | `detected_2pep = true` | abs ≤ max(1, 2 %) |
| `n_peptides_placed` | distinct `peptide` in `mapping/peptide_hits.tsv` | abs ≤ max(1, 2 %) |
| `n_abundance` | `genes.tsv` rows with `abundance_raw` present and `abundance_imputed = false` (A, B, D) | abs ≤ max(1, 2 %) |
| `warnings` | `run_manifest.json › warnings` codes | informational, not compared |
| `wall_time_s` | the measured wall time of the route | informational, not compared |

R7 makes a run byte-identical for the same inputs, seed, versions and images, so on FarmShare the
numbers should match exactly. The tolerances are for users on other machines, thread counts or
runtimes (`docker` versus `apptainer`), where tool tie orders may differ. `expected/README.md`
shows the same numbers as a table and says where each one is in the outputs. It says that a Bakta
full DB, or an installed UniOP, can change the gene and operon numbers.

### Caveats the README states

- The numbers are for a region, not the genome: shared peptides of a region gene and an
  out-of-region paralog look unique here, and genome-wide statistics see few genes. Do not compare
  them with the pilot.
- Route C keeps placed rows plus a share of unplaced rows (rule 5); its precision is not the
  pilot's.
- Preflight info rows expected on the region (for example `genome.small_replicons`, a short total
  length or N50) are listed, so users do not read them as faults.

## 3. How the bundle is built

`scripts/example/build_example.py` (validation agent, branch `feat/example-dataset`; internals not
designed here) reads the pilot files under `$GIPPA_ROOT/pilot/` and the committed
`validation/pilot/db_a2/`, applies the subset rules above, writes the layout above
deterministically, and writes or checks `validation/example/bundles.tsv`. It runs through
`scripts/slurm/build_example.sbatch` (job name `gippa-example-build`, R14) and leaves the bundle,
read-only, at `$GIPPA_ROOT/example/gippa-example-<X.Y>.tar.gz`. Requirements on it:

- **Reproducible.** Two builds at the same commit give the same archive sha256 (the sbatch builds
  twice and compares).
- `--check` mode: rebuild and compare with the recorded sha256 in `bundles.tsv`.
- `expected/` is filled from the slow route test's outputs (section 5), not computed by hand.

`validation/example/bundles.tsv` (committed; the "release asset manifest"): one row per bundle,
columns `version`, `file`, `size_bytes`, `sha256`, `url`, `release_tag`, `build_commit`,
`build_job`, `expected_commit`, `asset_digest` (the GitHub API `digest` of the uploaded asset;
empty until Khoa uploads it).

### Region, size and modules: TBD until measured

The 5-minute target (Khoa) is the wall time of **one route's commands** (convert, preflight,
run) from start to `report/report.html`, after `gippa setup --minimal` and with the images already
pulled. Measured on a FarmShare `wheat` node, 4 threads and 16 GB, Apptainer, no `--resume`, as a
Slurm job. Setup and image pulls are not counted, and the docs say so.

| Route | Region (kb) | Genes | convert | preflight | run | Total | ≤ 5 min? | Job | Commit |
|---|---|---|---|---|---|---|---|---|---|
| A | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD |
| B | TBD | TBD | n/a | TBD | TBD | TBD | TBD | TBD | TBD |
| C | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD |
| D | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD |
| `gippa demo --full` (SIM100K, the fixed-cost yardstick) | 108 | 85 | n/a | TBD | TBD | TBD | n/a | TBD | TBD |

Measure `gippa demo --full` first: it is the same pipeline on 108 kb, so it bounds the fixed cost
(Nextflow start, container start, Bakta DB load) before any region is chosen. Also record the
Bakta task's own time from `logs/trace.tsv`.

**Fallback if real tools cannot meet 5 min**, in this order, each re-measured:

1. a smaller region, down to a floor of about 40 genes in the region, so that each route still
   detects genes and route A still has decoys;
2. drop the plasmid, if one was kept;
3. faster **real-tool** settings for the example commands: `--terminators transtermhp`,
   `--sensitivity off` (`--promoters targeted` is already the default);
4. if the fixed cost alone (step 0: `demo --full`, or the Bakta task) is above 5 min, publish the
   measured time instead ("about N min") and ask Khoa (Q-K3). Precomputed layers or a Bakta JSON
   are **not** a fallback (Khoa: real tools only).

The docs quote only measured times.

## 4. Docs

- **New page `docs/get-started/example.md`, "Try the example"**, in the Get started nav right after
  the quick start. Content: what the bundle is (one paragraph, PXD078459, Nissle region, licence
  link); prerequisites (`pip install gippa`, `gippa setup --minimal`, a container runtime); `gippa
  example fetch`, and `--from-file` for offline compute nodes; routes A–D, about three commands
  each, taken from `example_routes.yaml`; how to compare with `expected/README.md`; the measured
  time table; a pointer to [From raw spectra](from-spectra.md) for how such inputs are made. No
  peptide sequence or other row-level content on the page (the lifted rule covers the bundle, not
  the site); the M11-08 guard (`PEPTIDE_LIKE`) covers this page too.
- **Quick start.** In section 1, after `gippa demo`, one line:
  `gippa example fetch   # real Nissle data, four input routes, about N min each after setup ([Try the example](example.md))`.
- **`docs/tutorials/ecoli-k12.md`** is a stub in the nav. Default (Q-K4): replace it with
  `docs/tutorials/example-routes.md`, "Tutorial: one Nissle region, four routes": run the four
  routes, read `genes.tsv` and the report per route, see how the tier changes with the input
  (`decoy_picked`/`mayu`, `protein_list`, `two_peptide`), and compare the four runs with
  `gippa compare {out}/A {out}/B {out}/C {out}/D --labels A,B,C,D`. Aggregates only. The K-12
  tutorial moves to `docs/backlog.md`. `tutorials/nissle-pilot.md` stays a stub (out of scope).
- **`from-spectra.md`** and **`user-guide/convert.md`**: one sentence each linking the example as a
  way to try routes A, C and D on real files.
- **CLI reference.** `scripts/docs/gen_reference.py render_cli` recurses into groups and renders
  `## gippa example fetch` with its option table; `docs/reference/cli.md` is regenerated.
- **`docs/release.md`.** New steps (section 5, item 4).
- Changelog fragment `changelog.d/<PR>.feature.md` in the implementation PR.

## 5. Tests and acceptance criteria

### Unit tests (network mocked; `tests/unit/test_m12_01_example_fetch.py`)

A fake fetcher (the `Fetcher` protocol) writes a small tar.gz built in `tmp_path`, and the test
monkeypatches the selected row (url, checksum, size) to that archive's values. No socket is opened
(the tests also monkeypatch `http.client.HTTPSConnection` to raise).

### CI-safe pin test (`tests/unit/test_m12_01_example_pin.py`)

Every `example_*` row of `db_sources.yaml` has exactly one row in
`validation/example/bundles.tsv` with the same version, sha256, size and URL, and vice versa; the
URL is `https://github.com/khoa-yelo/gippa/releases/download/v<X.Y>.0/gippa-example-<X.Y>.tar.gz`;
when `asset_digest` is filled it equals the sha256. No network. With no bundle rows yet the test
passes on the empty set; the pin PR adds the assertion that the current minor has a row.

A `needs_net(github.com)` test (not in CI, run in release step 14) downloads the live asset and
checks the pin.

### Slow route test (`tests/e2e/test_m12_01_example_routes.py`)

Markers `e2e`, `slow`, `needs_apptainer`, `needs_bakta_db`, `needs_promotech_models`. Submitted with
`scripts/slurm/example_routes.sbatch` (job `gippa-example-routes`, 4 CPUs, 16 GB, R14), never in
CI or the session allocation. It runs `gippa example fetch --from-file
$GIPPA_ROOT/example/gippa-example-<X.Y>.tar.gz` (so the offline path is exercised), then the
printed commands of each route, and compares each run with `expected/<route>.json` by the
tolerances there. It writes the per-route wall times to a TSV next to the outputs; times are
recorded, not asserted (a shared node is noisy; the table in section 3 is filled from them).

### Acceptance criteria

| # | Input | Expected | Test |
|---|---|---|---|
| 1 | fetch, fake fetcher serves the pinned archive, `--outdir` absent | exit 0; `DIR` holds the layout; every `CHECKSUMS.sha256` line verifies; stdout has the four route headings and their commands with the absolute `DIR`; no staging dir left | unit |
| 2 | fake fetcher serves bytes with another sha256 | exit 2; stderr `checksum mismatch … pinned sha256:<pin>`; no `DIR`, no staging dir | unit |
| 3 | fake fetcher raises `SourceError("github.com", "timed out")` | exit 2; stderr names the URL and the `--from-file` command; no `DIR` | unit |
| 3b | `SourceError` with `HTTP 404` | exit 2; message adds "may not be published yet" | unit |
| 4 | `--from-file` with the pinned archive; the fetcher constructor raises if called | exit 0; same tree as #1; `FILE` unchanged (bytes and mtime) | unit |
| 5 | `--from-file` with an archive that matches a second (monkeypatched) pin `1.2` | exit 1; stderr `this is example bundle 1.2; use --version 1.2` | unit |
| 5b | `--from-file` empty file, or a directory, or missing | exit 1 (empty: size mismatch; missing: Click's error) | unit |
| 6 | `DIR` exists and is empty | exit 0, fetched into it | unit |
| 6b | `DIR` holds a previous fetch of the same bundle; fetcher raises if called | exit 0, `already fetched`, commands printed | unit |
| 6c | `DIR` holds any other file | exit 1, nothing written | unit |
| 7 | `--version 9.9` | exit 1; stderr lists the pinned versions | unit |
| 7b | installed `1.2.0`, only `1.1` pinned | uses `1.1`, warning on stderr, exit 0 | unit |
| 7c | installed `2.0.0`, only `1.x` pinned | exit 2 `no example bundle for GIPPA 2.x` | unit |
| 7d | installed `1.1.0rc1` / `1.1.0.dev3`, `1.1` pinned | uses `1.1`, no warning | unit |
| 7e | `--version 1.1.4` | read as `1.1` | unit |
| 8 | archive with a `../x` member, an absolute member, or a symlink leaving the tree (pin set to its sha256) | exit 2 `unsafe member …`; nothing extracted outside staging; staging removed | unit |
| 8b | archive with two top-level directories, or a file at top level | exit 2 | unit |
| 9 | a member altered after `CHECKSUMS.sha256` was written; or an extra unlisted file | exit 2 naming the file | unit |
| 10 | fetcher raises `OSError(ENOSPC)` | exit 2 with the disk-full advice; nothing left | unit |
| 11 | the implementation calls `HttpsFetcher`, `file_digest` and `_unpack` (or its public alias) from `gippa.setupdb`; no second HTTPS or tar-extraction code | AST/import check: `src/gippa/example.py` imports no `http.client`, `urllib.request`, `tarfile` or `zipfile` | unit |
| 12 | `items.plan()` for every `SetupOptions` combination | contains no `example_*` item | unit |
| 13 | `db_sources.yaml` example rows and `validation/example/bundles.tsv` | equal sets (CI-safe pin test above) | unit |
| 14 | `example_routes.yaml` | every command starts with `gippa`, names only options in `docs/reference/cli.md`, uses only `{dir}`/`{out}`; route A has `--from diann`, B `--proteins` and `--abundance`, C `--from cascadia`, D `--from generic`; no `--annotation`, `--strict`, `--resume` or `--executor` | unit |
| 15 | `docs/get-started/example.md` | its fenced route commands equal `example_routes.yaml`'s, with `{dir}` written as `gippa-example`; page in nav after the quick start; no `PEPTIDE_LIKE` match | unit (docs) |
| 16 | `docs/reference/cli.md` | has `## gippa example fetch` with `--outdir`, `--version`, `--from-file`; `gen_reference.py --check` passes | unit (docs) |
| 17 | the slow route test on the published bundle | all four routes exit 0, `status: completed`, every compared key within tolerance | e2e (Slurm) |
| 18 | `build_example.py` twice at one commit | same sha256; `--check` passes against `bundles.tsv` | the build script's own tests + the sbatch |
| 19 | `docs/release.md` | has the build, attach and verify steps below | unit (docs) |

### Release procedure (`docs/release.md`)

- **Step 6 (release-prep, agent), minor or major:** `sbatch scripts/slurm/build_example.sbatch`;
  `sbatch scripts/slurm/example_routes.sbatch` on the new bundle; fill `expected/` and rebuild
  (the archive sha changes with `expected/`; the sbatch does both in one job); add the
  `example_<X>_<Y>` row to `db_sources.yaml` and the row to `bundles.tsv` (empty `asset_digest`);
  fill the measured table on the docs page. **Patch:** run `example_routes.sbatch` on the minor's
  bundle; a difference outside tolerance means the release is not a patch (ADR-0060).
- **Step 13b (Khoa), right after `release.yml` creates the release:**
  `gh release upload v<X.Y>.0 $GIPPA_ROOT/example/gippa-example-<X.Y>.tar.gz --repo khoa-yelo/gippa`,
  then `gh release view v<X.Y>.0 --json assets` and check the asset's `digest` equals the pin.
  Only Khoa uploads (publishing is Khoa's). Patch releases attach nothing.
- **Step 14 (smoke):** add `gippa example fetch --outdir "$smoke/example"` (exit 0) and the
  `needs_net` pin test.
- **Step 15 (post-release PR):** record `asset_digest` in `bundles.tsv`.
- Zenodo archives the source zipball, not release assets, so the bundle does not reach Zenodo
  (Khoa: not on Zenodo). Step 13b checks that the Zenodo record has no bundle file.

## 6. Edge cases

| Case | Expected |
|---|---|
| Linear region of a circular chromosome | written without a topology tag, so it is linear; nothing wraps; boundaries are intergenic, so no `partial_excluded` gene from the cut. Preflight's `genome.topology` info row says linear |
| Plasmid kept (circular) | tagged `[topology=circular]`; a gene or peptide may cross its origin; `expected/` covers it |
| Region must not cross the CP171242.1 origin | the build refuses such a region (it would need an unwrapped interval and a rotated contig) |
| Minus-strand genes at a region end | containment (rule 2) does not depend on strand; a minus-strand gene ending at the boundary gap is wholly inside or wholly outside |
| Shared peptides across the boundary | a row is kept when any of its proteins is in the region; its out-of-region IDs stay in `Protein.Ids` and map to no gene (they vote nowhere). README caveat |
| Decoy pairing on the region | decoys pair through the peptide vote of `protein_gene_map` (pilot locus tags are not the region's Bakta locus tags). If < 80 % pair, route A's tier is `mayu` with `decoy_pairing_low`; `expected/A.json` records whichever is measured. The lead decides whether a `mayu` route A is acceptable (Q-L2) |
| Empty route after subsetting | build fails (rule 6); never shipped |
| Empty `--from-file` | exit 1, size mismatch |
| Empty existing `--outdir` | allowed |
| Two pins with equal version | impossible: the YAML key is unique; the loader rejects two rows with the same `version` (ValueError at import, caught by the unit tests) |
| Region ties in the build | equal scores: the lower start wins |
| Installed GIPPA newer than every pin | newest pin of the same major, with the warning; another major: exit 2 |
| The asset is not uploaded yet (between `release.yml` and step 13b) | 404, exit 2 with "may not be published yet"; step 13b follows the release immediately |
| User has Bakta full DB, or UniOP | runs fine; gene or operon numbers can differ; README says so |
| User passes `--strict` | the region's preflight warnings give exit 3; the printed commands never use it |
| Bundle on a compute node without internet | `--from-file`, documented on the page and in the no-network message |
| Disk full in the middle of unpacking | exit 2; staging removed; `DIR` not created |

## 7. Open questions

### For Khoa (each with the reversible default used until he answers)

| # | Question | Default |
|---|---|---|
| Q-K1 | The data licence. PXD078459 states no licence. Use the statement above, or contact the depositors for permission or a licence first? Also: is redistributing DIA-NN 2.3.1 Academia **output** rows acceptable to you under its licence? | the statement above (no licence granted, attribution, demonstration and testing use, withdraw on request); DIA-NN output is data, not software, so it is shipped. Reversible: delete the asset and the pin row |
| Q-K2 | Bundle cadence: one per minor (above), or one per release? | one per minor; patches reuse it |
| Q-K3 | What "5 min" means and what to do if real tools cannot reach it | one route's convert + preflight + run, after setup, 4 threads, 16 GB; if unmet after the fallbacks, publish the measured time, still no precomputed layers |
| Q-K4 | Replace the empty `tutorials/ecoli-k12.md` with the example tutorial? | yes; the K-12 tutorial goes to the backlog |
| Q-K5 | The asset upload is a manual step by you (13b), since the bundle is built on FarmShare and `release.yml` cannot see it. Fine, or should the release wait on an agent-built draft? | manual step 13b by Khoa |
| Q-K6 | Keep one whole Nissle plasmid (circular) in the region bundle? | yes, if the measured time allows |
| Q-K7 | Paths in `File.Name` and the `.ssl` file column are shortened to base names. Acceptable? | yes |

### For the lead

| # | Question | Default |
|---|---|---|
| Q-L1 | ADR-0074 (R1, R10): the command, hosting, the lifted row-level rule, the licence statement | the lead writes it as PR 2; Khoa accepts before PR 5 merges |
| Q-L2 | Is a route A at tier `mayu` (decoy pairing below 80 % on the region) acceptable, or must the region be chosen so that `decoy_picked` holds? | prefer regions where it holds; accept `mayu` if no region under the time target gives it, and say so in the README |
| Q-L3 | `_unpack` is private in `gippa.setupdb.run`; make a public alias? | yes, `unpack_archive` with `_unpack` kept |
| Q-L4 | Group rendering in `gen_reference.py` and the `test_m11_08_docs.py` option-drift parser (`## gippa <name>` assumes one word) both need changes for `gippa example fetch` | yes; the test-author updates the parser in PR 4 |
| Q-L5 | The nav and page lists in `test_m11_07_mkdocs_config.py` and `test_m11_08_docs.py` are exact; adding `get-started/example.md` and renaming the tutorial needs test edits | test-author PR 6a, before the docs PR |
| Q-L6 | An M12 index (`docs/specs/M12-index.md`)? | not now; this is M12's first issue |

## 8. Work split (≤ ~400 lines each, in order)

| # | PR | Role | Contents | Code / test lines | Depends on |
|---|---|---|---|---|---|
| 1 | `docs(spec): M12-01 example dataset and gippa example` | spec-writer | this file | 0 | — |
| 2 | `docs(adr): ADR-0074 example bundle` | lead | the ADR; Khoa accepts | 0 | 1 |
| 3 | `feat(example): build script` (`feat/example-dataset`) | validation + test-author | `scripts/example/build_example.py`, `scripts/slurm/build_example.sbatch`, `validation/example/bundles.tsv` (header only), tests on a synthetic pilot fixture (subset rules 1–7, determinism, scrubbing) | ≈ 300 / 300 | 1 |
| 4 | `test(example): gippa example fetch` | test-author | `test_m12_01_example_fetch.py` (AC 1–12), `test_m12_01_example_pin.py` (AC 13), `test_m12_01_example_routes_yaml.py` (AC 14), parser update for groups (Q-L4); red | 0 / ≈ 400 | 1 |
| 5 | `feat(cli): gippa example fetch` | implementer | `src/gippa/example.py`, the `example` group in `cli.py`, `DbItem.size_bytes`, `example_routes.yaml`, the `unpack_archive` alias, `gen_reference.py` group rendering, regenerated `cli.md`, changelog fragment | ≈ 250 / 0 | 2, 4 |
| 6 | `test(example): slow route test` | test-author | `tests/e2e/test_m12_01_example_routes.py`, `scripts/slurm/example_routes.sbatch` | ≈ 40 / ≈ 150 | 5 |
| 7 | `chore(example): measure and pin bundle 1.1` | validation | the runs of section 3, the filled time table (this spec and the page), `expected/`, the `example_1_1` row, the `bundles.tsv` row, the "current minor has a row" assertion (test-author commit) | ≈ 20 / ≈ 10 | 3, 5, 6 |
| 8a | `test(docs): example page and tutorial` | test-author | nav lists, `TUTORIALS`, AC 15, 16, 19; red | 0 / ≈ 80 | 4 |
| 8b | `docs: try the example, tutorial, quick start` | docs | `get-started/example.md`, `tutorials/example-routes.md` (replacing `ecoli-k12.md`), `mkdocs.yml`, quick start, from-spectra and convert links, backlog entry | 0 (prose) | 7, 8a |
| 9 | `docs(release): build and attach the example bundle` | docs | `docs/release.md` steps 6, 13b, 14, 15 | 0 (prose) | 7 |

PR 7 is part of the 1.1.0 release-prep (step 6), since `expected/` must come from the release
commit. Before PR 7 merges, `gippa example fetch` on a 1.1 dev install exits 2 (`no example bundle
for GIPPA 1.x` while no row exists), which the unit tests cover with monkeypatched rows.
