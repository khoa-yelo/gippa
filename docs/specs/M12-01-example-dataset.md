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
session); measurements: `scripts/example/README.md` on `feat/example-dataset` @0698587f
(2026-11-03) · blocked in part by #536 (route B).

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

Khoa's answers to this spec's questions (2026-11-03, in session):

- The four faster real-tool settings: "accept" (Bakta light DB, TransTermHP, Promoter Calculator
  off, `mobile` off; section 3).
- #536 (`gippa run` with `--proteins` and no `--peptides` writes no gene table): "fix", so route
  B works once that lands.
- The data licence (Q-K1): "use the data", with the default statement (section 2).

### Amendment 2026-11-06: route E, a six-frame database (Khoa, in session)

> "we should add the 6 frame db as part of the suite as well, as that is one of the routes mostly
> use by the user"

Further decisions in the same session (Khoa, 2026-11-06):

- Route E's data is the pilot's six-frame arm **A4**: the same PXD078459 runs WT_1–WT_3, searched
  by DIA-NN 2.3.1 against the `gippa makedb --six-frame-orfs` database, cut to the same 200 kb
  Nissle region. The 2026-11-03 lift of the "no row-level DIA-NN content" rule extends to **A4's
  rows inside the region**, and only to them.
- Route E ships in **1.1.0rc1**. The release-prep PR #572 waits for it; the 1.1 bundle is rebuilt
  and re-pinned with route E in it.

Every change this amendment makes is marked "(route E, 2026-11-06)" where it lands: the layout,
the routes table, subset rules 3E, 7 and 9, PROVENANCE, the licence statement, `expected/`, the
caveats, section 4, section 5 (slow test and AC 14, 15, 17, 20–24), the edge cases, the open
questions Q-L9 to Q-L13 and Q-K8, and work-split rows 10–12. Estimate of the addition: ≈ 150 code
lines (build script, `example_routes.yaml`, help text) + ≈ 250 test lines + prose.

**What A4 is** (read-only look at `$GIPPA_ROOT/runs/pilot/sixframe-fa49b8f`, 2026-11-06; counts
only):

| Item | Value |
|---|---|
| database | `db/strain_db.faa` from `gippa makedb --six-frame-orfs` (GIPPA 1.0.0, `--min-orf-aa 20`, no contaminants, code 11): 4,613 cds, 72 sorf, 143,068 orf records; the 4,685 cds/sorf records are byte-identical to A2's (`db_counts.json`) |
| ORF IDs | `GIPPAORF_{class}_{contig}_{start1}_{end1}_{p\|m}`, classes `antisense`, `cds_outframe`, `intergenic`, `nterm_extension`; in `db_id_map.tsv` they have `kind = orf`, `class` set; `locus_tag` NA except `nterm_extension` (the full-DB Bakta gene it extends) |
| search | DIA-NN 2.3.1 Academia, `jobs/diann_A4.cfg` = `diann_A2.cfg` with only `--fasta` changed (plus the 381 contaminant records), job 1918156; `A4/command_step1.txt`, `command_step2.txt`, `diann_version.txt`, `provenance.json` (image sha256 `1946452d…`) |
| report | `A4/A4_report.parquet`: 91,076 rows (901 decoys), 71 columns |
| genome-wide run | `gippa convert --from diann … --search-db same` (`jobs/convert_A4.sbatch`), then `gippa run --annotation <full-DB Bakta JSON> --peptides a4.tsv.gz` (job 1921228): tier `decoy_picked`, 109 novel-ORF candidates (all `picked_decoy`, paired share 1.0), **0 called**: each class had 8–64 target wins against the 100 needed for a q floor ≤ 0.01 (`docs/user-guide/makedb.md`, PROGRESS › A4) |

**In the region** (same rule as the build will apply; measured 2026-11-06, aggregates only):

| Item | Count |
|---|---|
| region entries of `db_id_map.tsv` (contained, either strand) | 153 cds, 1 sorf, 5,678 orf (3,023 antisense, 2,101 cds_outframe, 447 intergenic, 107 nterm_extension) |
| entries crossing a region end (dropped) | 9 orf, 2 cds |
| contained ORFs whose stop codon falls outside the region | 0 |
| A4 rows kept by rule 3E | **2,218** (2,181 targets, **37 decoys**); per run 749 / 730 / 739 |
| kept rows with at least one ORF ID | 63 (24 with ORF IDs only, of which 16 are decoys) |
| distinct region ORF entries in the kept rows (targets and decoys) | 35 |
| rows that rule A alone would keep on A4 (no ORF IDs) | 2,155 |
| columns holding ORF IDs | `Protein.Ids`, `Protein.Group` |
| `E_sixframe/report.parquet` (zstd, all 71 columns) | ≈ 387 kB, so the bundle grows from 0.48 MB to **≈ 0.9 MB** (estimate; parquet barely compresses further in the tar.gz) |

**Why route E needs an ID rewrite (rule 9).** `gippa run` recognises six-frame entries from the
IDs alone: `gippa.novel.class_call.makedb_id_map` parses every `GIPPAORF_…` ID of the peptide
table with `gippa.novel.fdr.parse_orf_id`, and `picked_class_fdr` pairs a target ORF with the entry
of the same `(contig, strand, stop)`. No FASTA, `db_id_map.tsv` or makedb manifest is read, and no
extra flag is needed (`--search-db same` as for route A). The bundle's genome, however, is the
region written as the contig `nissle1917_NZ_CP058217_1_1965000_2165000` with coordinates starting
at 0 (`build_example.py`). A4's IDs name `NZ_CP058217.1` and chromosome coordinates, so unchanged
IDs never match a candidate's key: every candidate would get `no_makedb_pair` and route E would
show nothing of the six-frame workflow. Rule 9 therefore rewrites each region ORF ID to the
region's contig and coordinates (Q-L9).

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
   **same** bundle: print `already fetched`, print the commands, exit 0, no network. Anything
   else: exit 1 `<DIR> exists and is not this example; use --outdir or remove it`.
   *Amendment (2026-11-04, from the #542 review):* a re-fetch is detected with the marker
   `DIR/.gippa-example.json`, written at unpack time and holding the archive sha256, the bundle
   version and the GIPPA version, because an archive cannot hold its own sha256 (so `bundle.json`
   cannot be compared with the pin). The marker must name the selected pin (archive sha256 and
   bundle version), and every file listed in `CHECKSUMS.sha256` must verify. A re-fetch checks
   only the listed files: unlisted entries, such as the routes' outputs under `DIR/runs/`, are the
   user's and are ignored. The marker is not listed in `CHECKSUMS.sha256`.
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
   that the routes need `gippa setup --no-genomad` (Bakta light and the antiSMASH databases) and
   a container runtime.
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
├── SELECTION.json         the selection rule, counts, sha256 of every source and output (written by the build script)
├── genome/
│   └── <name>.fna.gz      the Nissle region, one linear contig
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
├── E_sixframe/
│   └── report.parquet     DIA-NN 2.3.1 report rows of the region from the six-frame search A4 (targets and decoys; ORF IDs in region coordinates, rule 9)
└── expected/
    ├── A.json  B.json  C.json  D.json  E.json
    └── README.md          the same numbers as a table, with the file and column each comes from
```

*Route E, 2026-11-06:* `E_sixframe/` holds only the report. The bundle ships **no** six-frame
FASTA and no `db_id_map.tsv` subset: `gippa convert` and `gippa run` do not read them (the IDs
carry class and coordinates; `nterm_extension` locus tags come from the run's own gene model,
`makedb.build.nterm_entries`). PROVENANCE names the database by its sha256 and command (Q-L11).

The bundle is **0.48 MB** (measured, routes A–D); with route E it is **≈ 0.9 MB** (estimated from
the region subset; the rebuild in work-split row 12 measures it). The archive is a deterministic tar.gz: members sorted, mode 0644/0755, uid/gid 0, empty
user/group names, mtime 0, gzip header without name and with mtime 0. Two builds from the same
inputs give the same sha256.

### Routes (one input set each, same region)

| Route | Input | Commands (about three) | Tier (measured at 200 kb) |
|---|---|---|---|
| A | `A_diann/report.parquet` | `gippa convert --from diann … --out {out}/A/peptides.tsv --search-db same`; `gippa preflight`; `gippa run --peptides …` | **`mayu`**: 158 genes, 99 detected, gene FDR 0.012 at q 0.01. 37 decoy rows survive, but only 4 of 10 decoy peptides pair with a gene, below the 80 % `decoy_picked` threshold (warning `decoy_pairing_low`) |
| B | `B_proteins/proteins.faa` + `abundance.tsv` | `gippa preflight --proteins … --abundance …`; `gippa run --proteins … --abundance … --abundance-type maxlfq` | `protein_list`, **after #536**. In GIPPA 1.0.0 `workflow/main.nf` runs MAPPING only with `inputs.peptides`, so this route runs ANNOTATE and REPORT and writes no gene table |
| C | `C_cascadia/*.ssl` | `gippa convert --from cascadia … --out {out}/C/peptides.tsv`; `gippa preflight`; `gippa run --peptides …` | `two_peptide`, `c_assumed`; 76 rows, **4 genes detected**. A **format demo**, not a result |
| D | `D_generic/peptides.tsv` | `gippa convert --from generic … --out {out}/D/peptides.tsv --source database --search-db same`; `gippa preflight`; `gippa run --peptides …` | the same decoys as A, so expected `mayu`; the measured value goes into `expected/D.json` (PR 7) |
| E (route E, 2026-11-06) | `E_sixframe/report.parquet` | `gippa convert --from diann … --out {out}/E/peptides.tsv --search-db same`; `gippa preflight`; `gippa run --peptides …` (no extra flag; see below) | expected **`mayu`**, as A: decoys of ORF entries that map to no gene leave the gene-level decoy counts (M7-12 decision 28), so the 21 non-ORF decoys pair through peptide votes as in A (Q-L10). **Novel ORFs:** a few `picked_decoy` candidates, **0 called** (below) |

**Route E's commands** (`example_routes.yaml`, title `Route E: DIA-NN report, six-frame
database`), with the `gippa run` settings of every route:

```
gippa convert --from diann {dir}/E_sixframe/report.parquet --out {out}/E/peptides.tsv --search-db same
gippa preflight --genome {dir}/genome/nissle_r200k.fna.gz --peptides {out}/E/peptides.tsv --set 'modules=[core,bgc]' --set annotation.bakta.db=light --set annotation.terminators.tool=transtermhp --set annotation.promoter_calculator.enabled=false
gippa run --genome {dir}/genome/nissle_r200k.fna.gz --peptides {out}/E/peptides.tsv --genome-relation same --outdir {out}/E --modules core,bgc --set annotation.bakta.db=light --terminators transtermhp --set annotation.promoter_calculator.enabled=false --threads 4
```

They are route A's with `A_diann` → `E_sixframe` and `{out}/A` → `{out}/E`. The genome-wide A4
run used `--annotation`; route E does not (no route does, AC 14): live Bakta light runs as in A,
and the six-frame IDs do not depend on Bakta's locus tags. The database itself is not rebuilt by
the route; the example page shows, in prose, the `gippa makedb --six-frame-orfs` command that made
it (PROVENANCE has the exact one).

**Route E's novel ORFs.** Novel ORFs run in `core` for every route; route E is the one whose
decoys are class-tagged makedb entries, so its candidates take the `picked_decoy` regime and pair
with their entries (rule 9). **No call is possible on the region:** a class's q floor is
`1 / n_target_wins`, and the region's kept rows hold 35 distinct ORF entries (targets and decoys
together), so every class has far fewer than the 100 target wins that a floor ≤ α = 0.01 needs.
Genome-wide, A4 called 0 of 109 for the same reason. `n_novel_called = 0` is therefore exact
(`expected/` below), and route E demonstrates the workflow (IDs, regime, class summaries, the
`floor_above_alpha` flag), not a discovery.

Every `gippa run` uses `--genome {dir}/genome/<name>.fna.gz --genome-relation same --outdir
{out}/<route> --threads 4` and the example settings Khoa accepted:

```
--modules core,bgc --set annotation.bakta.db=light --terminators transtermhp \
  --set annotation.promoter_calculator.enabled=false
```

No route uses `--annotation`, `--strict` or `--resume`. Bakta, antiSMASH, Promotech, TransTermHP,
miniprot and the mappers run live (Khoa: real tools only). `function` is off: the Bakta light DB
carries no COG, GO or KEGG IDs, so the module would give 0 rows. `bgc` needs the antiSMASH
databases, which `gippa setup --minimal` leaves out, so the prerequisite is `gippa setup
--no-genomad` (Bakta light is the default) and not `--minimal`.

**Route A's tier is `mayu`, and the README says so (Q-L2, resolved: accept `mayu`).** Live Bakta
writes new locus tags, so decoys pair only through peptide votes. Adding the search FASTA
(`--proteins`) pairs every decoy (tested at 100 kb), but the `decoy_picked` q floor is then 1/24 >
0.01 and 0 genes are detected. With this few decoys, `mayu` is the honest tier.

**FDR caveat (README).** The q-values, PEPs and MaxLFQ values of routes A, B and D are DIA-NN's,
computed on the whole-proteome search (4,685 strain proteins plus contaminants, three runs). The
subset does not recompute them. *Route E, 2026-11-06:* route E's are DIA-NN's on the whole
six-frame search A4 (147,753 strain records plus 381 contaminants, the same three runs); a
29-fold larger database changes the scores, so route E's gene numbers differ from route A's.

**Route D's table.** Built from route A's rows, grouped by (stripped peptide, run), as an engine
without a GIPPA converter would write it. The columns are `peptide`, `modified_peptide` (ProForma,
`(UniMod:n)` written as `[UNIMOD:n]`, so the generic converter's stripping is exercised),
`protein_ids`, `q_value`, `pep`, `is_decoy`, `intensity` and `run`.

### Subset rules (as `build_example.py` implements them)

The script's docstring holds the rule and copies it verbatim into the bundle's `SELECTION.json`,
with the sha256 of every source and output file. In summary:

1. **Region.** `NZ_CP058217.1:1965001-2165000` (1-based closed; `[1965000, 2165000)` 0-based
   half-open; 200 kb), the chromosome of GCF_021559835.1 (ASM2155983v1; GenBank twin
   GCA_021559835.1), the pilot's pass-3 assembly (`validation/pilot/README.md` › pass 3). It holds
   the colibactin *pks* island and the yersiniabactin HPI. antiSMASH calls both as one
   NRPS/T1PKS/NRP-metallophore region over the first 117 kb. It is written as one linear contig
   (no topology tag) and does not cross the origin. **No plasmid** (Q-K6, resolved): the region is
   a piece of the chromosome. It holds no prophage, because geNOMAD is too slow to include.
2. **Region proteins.** Records of the pass-3 A2 strain database (`db_id_map.tsv`, 1-based closed)
   on `NZ_CP058217.1` with `start - 1 >= region start` and `end <= region end` (CDS and sORF), on
   either strand.
3. **DIA-NN rows (A).** A row is kept iff **every** ID in `Protein.Ids` is a region protein.
   The same rule applies to target and decoy rows (DIA-NN decoys carry the target IDs with
   `Decoy = 1`), so the target/decoy balance of the region is not biased by the rule. Rows with a
   contaminant (`Cont_`) or out-of-region ID are dropped. All 71 columns are kept.

   **3E. DIA-NN rows (E; route E, 2026-11-06).** Source: A4's `A4/A4_report.parquet` and
   `db/db_id_map.tsv` under `$GIPPA_ROOT/runs/pilot/sixframe-fa49b8f/` (read in place, Q-L13).
   **Region entries** are the `db_id_map.tsv` records on `NZ_CP058217.1` with
   `start - 1 >= region start` and `end <= region end`, of kind `cds`, `sorf` **or `orf`**, on
   either strand (rule 2, with `orf` added; an ORF is placed by its map coordinates, never by
   parsing its ID). A row is kept iff **every** ID in `Protein.Ids` (split on `;`) is a region
   entry. Target and decoy rows get the same rule (A4's decoys carry the target IDs with
   `Decoy = 1`, as A2's). Rows with a contaminant (`Cont_`) ID, an out-of-region protein or an
   out-of-region ORF (including one that crosses a region end) are dropped. All 71 columns are
   kept; rule 9 then rewrites the ORF IDs. Route A is not changed: it still reads A2.
4. **Protein table (B).** The protein groups of route A's target rows with a non-missing, non-zero
   `PG.MaxLFQ`; value = median over runs (the `gippa convert --abundance-out` rule). A group `A;B`
   is one row per member ID. The FASTA holds each member's sequence from the strain database.
5. **Cascadia rows (C).** A row is kept iff its stripped sequence (I and L equated) is an exact
   substring of a stop-free stretch of the region's six-frame translation (code 11). Every kept
   row is therefore placeable, which inflates the placement rate. The README says so, and calls
   route C a format demo.
6. **Route D.** As described above.
7. **No empty route.** The build fails if any route keeps zero rows or zero proteins, or if route
   A keeps no decoy. *Route E, 2026-11-06:* also if route E keeps no target row, no decoy row, or
   no row with a region ORF ID (measured: 2,181, 37 and 63).
8. **Scrubbing (R12).** The `.ssl` `file` column is reduced to the mzML base name, whose stem is
   the run name. Any DIA-NN column holding an absolute path is shortened the same way (route E
   included), and `SELECTION.json` lists each rewritten column. No other value is changed,
   except by rule 9.
9. **ORF IDs in region coordinates (route E, 2026-11-06; Q-L9).** In `E_sixframe/report.parquet`
   every ORF ID (`GIPPAORF_{class}_NZ_CP058217.1_{s1}_{e1}_{p|m}`, in `Protein.Ids` and
   `Protein.Group`, and in any other string column that holds one) becomes
   `GIPPAORF_{class}_{C}_{s1 - S}_{e1 - S}_{p|m}`, where `C` is the bundle genome's contig ID
   (`nissle1917_NZ_CP058217_1_1965000_2165000`) and `S` the 0-based region start (1965000): the
   region-frame 1-based closed coordinates of the same interval. Class, strand and the ID order
   within a cell are unchanged; CDS and sORF IDs are unchanged (they carry no coordinates).
   Checks, each a build failure:
   - the ID's parsed coordinates equal its `db_id_map.tsv` row (the ID is the map's);
   - the rewritten ID parses with `gippa.novel.fdr.parse_orf_id` to `(class, C, s1 - S - 1,
     e1 - S, strand)`;
   - for every class but `nterm_extension`, the region sequence at the rewritten interval,
     translated with code 11 on its strand, equals the ID's record in A4's `strain_db.faa`
     except the first residue (a start codon read as M);
   - no ORF ID naming `NZ_CP058217.1` is left in the file.
   `SELECTION.json` records the rule, the columns rewritten and the number of distinct IDs
   rewritten (35 expected). The converted table (`{out}/E/peptides.tsv`) carries the rewritten
   IDs, so a user's ORF IDs in route E name the bundle's contig, as the outputs do.

   Examples (synthetic IDs and rows, `S = 1965000`, `C` as above):

   | Input (A4 row's `Protein.Ids`, map coordinates 1-based closed) | Expected |
   |---|---|
   | `CDS_X` (region CDS), target | kept; IDs unchanged |
   | `GIPPAORF_antisense_NZ_CP058217.1_1965101_1965400_m` (map 1965101–1965400), decoy | kept; becomes `GIPPAORF_antisense_C_101_400_m`; `parse_orf_id` gives `(antisense, C, 100, 400, -)` |
   | `CDS_X;GIPPAORF_cds_outframe_NZ_CP058217.1_2000002_2000301_p`, target | kept; second ID becomes `GIPPAORF_cds_outframe_C_35002_35301_p`, order kept |
   | `GIPPAORF_intergenic_NZ_CP058217.1_2164900_2165120_p` (crosses the region end) | dropped |
   | `CDS_X;GIPPAORF_antisense_NZ_CP058217.1_10_309_m` (ORF outside the region) | dropped |
   | `Cont_P1;GIPPAORF_antisense_NZ_CP058217.1_1965101_1965400_m` | dropped |
   | `GIPPAORF_intergenic_NZ_CP058217.1_1965001_1965300_p` (ORF starts at the first region base) | kept; becomes `GIPPAORF_intergenic_C_1_300_p` |
   | `GIPPAORF_intergenic_NZ_CP058217.1_2164701_2165000_p` (ends at the last region base) | kept; becomes `GIPPAORF_intergenic_C_199701_200000_p` |

### `PROVENANCE.md` (content)

- The PXD078459 record (URL, retrieval date, record sha256 `5d3f16c3…`), the three runs WT_1–WT_3
  with their file names, the instrument.
- The genome: GCF_021559835.1 (`NZ_CP058217.1`; GenBank twin GCA_021559835.1), the sha256 of the
  pass-3 genome file; the region in both coordinate conventions.
- Route A: DIA-NN 2.3.1 Academia, the image sha256, `validation/pilot/diann_A2.cfg`, the input
  report's sha256, row counts kept (targets, decoys).
- Route B: the A2 convert table sha256s (`a2_abundance.tsv.gz`), protein counts.
- Route C: Cascadia 0.0.7, the weights sha256, the three `.ssl` sha256s, rows kept per run.
- Route D: how it was derived from route A.
- Route E (route E, 2026-11-06): the database: `gippa makedb --six-frame-orfs` at GIPPA 1.0.0
  (`makedb_manifest.json`: `--min-orf-aa 20`, no contaminants, code 11; the sha256 of
  `strain_db.faa`, `db_id_map.tsv` and the Bakta JSON input, the counts 4,613 / 72 / 143,068);
  the search: DIA-NN 2.3.1 Academia (`A4/diann_version.txt`), the same image sha256 as route A,
  `jobs/diann_A4.cfg` (sha256; `diann_A2.cfg` with only `--fasta` changed) and the two commands
  of `A4/command_step1.txt` and `command_step2.txt` with their paths shortened to base names
  (R12), the Slurm job 1918156, the input report's sha256; row counts kept (targets, decoys,
  rows with an ORF ID) and the rule-9 rewrite count. It says the genome-wide run of A4 called 0
  novel ORFs.
- The build: script path and commit, Slurm job ID, GIPPA version and commit that produced
  `expected/`, the Bakta DB (light 6.0) and image digests used, the run settings, and that UniOP was
  not installed.
- `SELECTION.json` (the rule and all sha256s) sits next to it.

### Licence statement (`LICENSE-DATA.md`; Khoa, 2026-11-03: "use the data")

> The GIPPA software is MIT-licensed. This licence does **not** cover the data files of this
> bundle.
>
> The files under `A_diann/`, `B_proteins/`, `C_cascadia/`, `D_generic/` and `E_sixframe/` are
> derived from the ProteomeXchange dataset **PXD078459** (iProX **IPX0017169000**), announced
> 2026-05-16. The dataset record states no licence, and its publication is pending. The files are
> redistributed for the sole purpose of demonstrating and testing GIPPA, with attribution to the
> depositors of PXD078459. The content is unchanged: rows were selected, file paths shortened,
> and in `E_sixframe/` the identifiers of six-frame ORFs (made by GIPPA, not by the depositors)
> renumbered to the region's coordinates. The GIPPA
> authors grant no licence to these files. Cite PXD078459 if you use them. If you are a rights
> holder and object, open an issue and the bundle will be withdrawn.
>
> `genome/` holds a region of the *E. coli* Nissle 1917 assembly GCF_021559835.1 (GenBank
> GCA_021559835.1), which is INSDC data with no use restrictions by NCBI. `expected/`, `README.md`,
> `PROVENANCE.md`, `SELECTION.json`, `bundle.json` and `CHECKSUMS.sha256` are written by GIPPA
> tooling and are MIT-licensed, like GIPPA.

*Route E, 2026-11-06:* the statement above adds `E_sixframe/` and the ORF-ID clause to Khoa's
accepted 2026-11-03 text; Khoa confirms the wording (Q-K8), and ADR-0074 gets the same amendment
(Q-L12).

### `expected/` (key numbers)

One JSON file per route, with a `tolerance` per key, made by running the routes (five, with
route E) with the
release's GIPPA at the release-prep commit, with the example settings (Bakta light DB 6.0, no
UniOP, default seed).

| Key | Source | Tolerance |
|---|---|---|
| `status` | `run_manifest.json › status` | exact |
| `detection_tier` | `run_manifest.json › detection_tier` | exact |
| `layers_missing` | `run_manifest.json › layers_missing` (as a set) | exact |
| `convert_rows` | data rows of the converted `peptides.tsv` (A, C, D, E) | exact (conversion is deterministic) |
| `n_genes` | rows of `tables/genes.tsv` | abs ≤ max(1, 1 %) |
| `n_detected` | `genes.tsv` rows with `detected = true` | abs ≤ max(1, 2 %) |
| `n_detected_2pep` | `detected_2pep = true` | abs ≤ max(1, 2 %) |
| `n_peptides_placed` | distinct `peptide` in `mapping/peptide_hits.tsv` | abs ≤ max(1, 2 %) |
| `n_abundance` | `genes.tsv` rows with `abundance_raw` present and `abundance_imputed = false` (A, B, D, E) | abs ≤ max(1, 2 %) |
| `n_novel_candidates` | **E only**: data rows of `tables/novel_orfs.tsv` (`#` lines skipped) | abs ≤ max(1, 2 %) |
| `n_novel_called` | **E only**: `tables/novel_orfs.tsv` rows with `called = true` | exact (0: no class can reach 100 target wins on the region, section 2 › Routes) |
| `n_novel_unpaired` | **E only**: rows of `novel/novel_class_fdr.parquet` with `reason = no_makedb_pair` | abs ≤ max(1, 2 %) |
| `novel_classes` | **E only**: `run_manifest.json › novel_orfs › classes`, per entry `regime`, `class`, `estimated`, `n_targets`, `n_target_wins`, `n_decoys`, `q_floor`, `floor_above_alpha` | informational, not compared |
| `warnings` | `run_manifest.json › warnings` codes | informational, not compared |
| `wall_time_s` | the measured wall time of the route | informational, not compared |

*Route E, 2026-11-06:* `expected/E.json` holds every key of the table (all apply to E); A–D keep their key sets (the four novel keys are not required for them). `n_novel_unpaired`
is compared, not asserted to be 0: live Bakta light can call a gene differently from the full-DB
annotation A4 was built on, and a peptide of such a gene then makes a candidate with no makedb
entry. AC 22 checks that the rewrite works.

R7 makes a run byte-identical for the same inputs, seed, versions and images, so on FarmShare the
numbers should match exactly. The tolerances are for users on other machines, thread counts or
runtimes (`docker` versus `apptainer`), where tool tie orders may differ. `expected/README.md`
shows the same numbers as a table and says where each one is in the outputs. It says that a Bakta
full DB, or an installed UniOP, can change the gene and operon numbers.

### Caveats the README states

- The numbers are for a region, not the genome: shared peptides of a region gene and an
  out-of-region paralog look unique here, and genome-wide statistics see few genes. Do not compare
  them with the pilot.
- Route A's tier is `mayu`, not `decoy_picked`: too few decoys pair on a region (above).
- Route C keeps only rows that match the region exactly (rule 5) and detects 4 genes: a format
  demo. Its precision is not the pilot's.
- Route B needs #536 (until then it writes no gene table).
- The q-values come from the whole-proteome search (FDR caveat above).
- Empty or not estimable on this region: islands (0), enrichment and GSEA (the light DB has no
  categories), autocorrelation (`not_estimable`), variants (0) and novel ORFs (0). Non-empty:
  46 operons with 4 significant operon flags, and the BGC contrast (bgc_1: 71 genes, 30 detected
  against 42.4 expected).
- **The defaults take longer.** The example uses four faster real-tool settings. A run with the
  defaults is slower; the measured default times are in section 3 and are repeated in the README.
- Preflight info rows expected on the region (for example `genome.small_replicons`, a short total
  length or N50) are listed, so users do not read them as faults.
- *Route E, 2026-11-06:* **Route E shows the six-frame workflow, not discoveries.** Expect a few
  novel-ORF candidates and **no call**: a class needs at least 100 target wins over its decoys
  before its q floor can reach 0.01, and the kept rows name 35 ORF entries (of 5,678 in the region). The genome-wide
  A4 run called 0 of 109 for the same reason (8–64 wins per class); `docs/user-guide/makedb.md`
  explains the floor. The novel-ORF line of the caveat above ("novel ORFs (0)") is for routes
  A–D.
- *Route E, 2026-11-06:* route E's ORF IDs name the bundle's region contig and region
  coordinates (rule 9), not `NZ_CP058217.1`; add 1,965,000 to get chromosome coordinates.
- *Route E, 2026-11-06:* route E's tier is expected to be `mayu`, like route A's: its ORF decoys
  pair no gene, and its other decoys pair through peptide votes.

## 3. How the bundle is built

`scripts/example/build_example.py` (validation agent, branch `feat/example-dataset`; internals not
designed here) reads the pilot's pass-3 files under `$GIPPA_ROOT/pilot/`, applies the subset rules
above, and writes `$GIPPA_ROOT/example/<name>/` with `SELECTION.json`. It runs in seconds.
`scripts/example/example_run.sbatch` converts one route and runs `gippa run` with
`-profile farmshare_local,apptainer`. `summarize_runs.py` prints the wall and per-process times.
Still to add (PR 3 or 7): the release layout above, the archive, and the writing or checking of
`validation/example/bundles.tsv`. The bundle is left read-only at
`$GIPPA_ROOT/example/gippa-example-<X.Y>.tar.gz`. Requirements on the script:

- **Reproducible.** Two builds at the same commit give the same archive sha256 (the sbatch builds
  twice and compares).
- `--check` mode: rebuild and compare with the recorded sha256 in `bundles.tsv`.
- `expected/` is filled from the slow route test's outputs (section 5), not computed by hand.
- *Route E, 2026-11-06:* `sources()` gains `a4_report` (`runs/pilot/sixframe-fa49b8f/A4/A4_report.parquet`),
  `a4_db_id_map` (`…/db/db_id_map.tsv`) and `a4_strain_db` (`…/db/strain_db.faa`, read only for
  the rule-9 translation check), all under `$GIPPA_ROOT` so the `_rel` guard holds. The build
  writes `E_sixframe/report.parquet` (rules 3E, 8, 9), adds route E to the rule-7 check, to
  `SELECTION.json`, `bundle.json` counts, README, PROVENANCE and LICENSE-DATA (`bundle_text.py`),
  and to the per-route sizing tars. It still runs in seconds (the A4 map is 15 MB, the report
  11 MB).

`validation/example/bundles.tsv` (committed; the "release asset manifest"): one row per bundle,
columns `version`, `file`, `size_bytes`, `sha256`, `url`, `release_tag`, `build_commit`,
`build_job`, `expected_commit`, `asset_digest` (the GitHub API `digest` of the uploaded asset;
empty until Khoa uploads it).

### Measured times (2026-11-03, wheat nodes, warm image cache)

The 5-minute target (Khoa) is the wall time of one route's `gippa run`, after setup and with the
images already pulled. Setup and image pulls are not counted, and the docs say so. Settings: the
200 kb region and the example settings above. Source: `scripts/example/README.md` @0698587f.

| Route | Region | CPUs (granted) | Nextflow span | Slurm elapsed | ≤ 5 min? | Top steps (s) |
|---|---|---|---|---|---|---|
| A | 200 kb | 9 | 5:05 | 5:46 | about (Nextflow 5:00–5:07) | MAP_TOLERANT 72, Promotech 65, Bakta 62, plots 31, antiSMASH 30 |
| A | 200 kb | 5 | 5:00 | 5:47 | about | Bakta 64, MAP_TOLERANT 61, Promotech 59 |
| C | 200 kb | 9 | 5:07 | 5:45 | about | |
| D | 200 kb | 9 | 4:59 | 5:41 | yes | |
| B | 200 kb | 9 | 3:08 | 3:41 | yes, but only ANNOTATE and REPORT run until #536 | |

`gippa convert` and `gippa preflight` take seconds and are not in the table. About 35 tasks run
in sequence along the critical path, and each `gippa _step` task takes about 10 s, mostly start-up.
So 4 and 8 CPUs give the same time, and so do 100 and 200 kb. The docs and printed commands use
`--threads 4`. The docs give the time as "about 5 min (Nextflow), about 6 min as a Slurm job".
Route B is re-measured after #536 (PR 7). *Route E, 2026-11-06:* not measured yet; expected
close to route A (same steps, a few more novel-ORF windows for Promotech and OSTIR). Work-split
row 12 measures it and fills the row.

**The four faster settings (accepted by Khoa, 2026-11-03).** Each one replaces a step whose time
does not fall when the region gets smaller:

| Default | Measured cost | Example setting |
|---|---|---|
| Bakta full DB (`psc.dmnd`, 46 GB) | 26 min on 55 kb: DIAMOND reads the whole database | `--set annotation.bakta.db=light` (1.3 GB): about 1 min |
| TerminatorNet | 4.5–4.8 min per strand of 100 kb on 4 CPUs; it sets the critical path | `--terminators transtermhp`: 1 s |
| geNOMAD (`mobile`) | 4.6–5.4 min, roughly fixed (100 kb) | module off (`--modules core,bgc`) |
| Promoter Calculator | 100 s on 100 kb, 428 s on 200 kb (one task) | `--set annotation.promoter_calculator.enabled=false`; Promotech still runs |

`function` is off too, for content rather than time: 13 s, but 0 rows with the light DB.

**Measured times with the defaults** (route A; the README and docs page repeat them, so users
know the defaults take longer):

| Region | Settings | Nextflow span | Slurm elapsed |
|---|---|---|---|
| 55 kb | `core`, Bakta full DB, other defaults | not finished: Bakta alone took 1,563 s (26 min) before a PATH bug in the first sbatch stopped the job at 26:51 | n/a |
| 100 kb | default modules, Bakta light | 13:45 | 14:14 |
| 100 kb | `core`, Bakta light | 8:51 | 9:28 |
| 200 kb | `core`, Bakta light, TransTermHP (Promoter Calculator on) | 10:29 | 11:07 |

A full-default run on the 200 kb region (full DB, all modules) was not measured. From the rows
above, expect well over 30 min. The docs say "not measured", not a number.

## 4. Docs

- **New page `docs/get-started/example.md`, "Try the example"**, in the Get started nav right after
  the quick start. Content: what the bundle is (one paragraph, PXD078459, Nissle region, licence
  link); prerequisites (`pip install gippa`, `gippa setup --no-genomad`, a container runtime); `gippa
  example fetch`, and `--from-file` for offline compute nodes; routes A–D, about three commands
  each, taken from `example_routes.yaml`; what the four faster settings change, and the
  measured default times; route A's `mayu` tier, route B's #536 note (until it lands) and route C as
  a format demo; how to compare with `expected/README.md`; the measured time table; a pointer to [From raw spectra](from-spectra.md) for how such inputs are made. No
  peptide sequence or other row-level content on the page (the lifted rule covers the bundle, not
  the site); the M11-08 guard (`PEPTIDE_LIKE`) covers this page too.
- **Quick start.** In section 1, after `gippa demo`, one line:
  `gippa example fetch   # real Nissle data, four input routes, about 5 min each after setup ([Try the example](example.md))`.
- **`docs/tutorials/ecoli-k12.md`** is a stub in the nav. Default (Q-K4): replace it with
  `docs/tutorials/example-routes.md`, "Tutorial: one Nissle region, four routes": run the four
  routes, read `genes.tsv` and the report per route, see how the tier changes with the input
  (`mayu`, `protein_list`, `two_peptide`), and compare the four runs with
  `gippa compare {out}/A {out}/B {out}/C {out}/D --labels A,B,C,D`. Aggregates only. The K-12
  tutorial moves to `docs/backlog.md`. `tutorials/nissle-pilot.md` stays a stub (out of scope).
- **`from-spectra.md`** and **`user-guide/convert.md`**: one sentence each linking the example as a
  way to try routes A, C and D on real files.
- **CLI reference.** `scripts/docs/gen_reference.py render_cli` recurses into groups and renders
  `## gippa example fetch` with its option table; `docs/reference/cli.md` is regenerated.
- **`docs/release.md`.** New steps (section 5, item 4).
- Changelog fragment `changelog.d/<PR>.feature.md` in the implementation PR.

*Route E, 2026-11-06:*

- **Example page:** a section `### Route E: DIA-NN report, six-frame database` holding exactly
  route E's three commands from `example_routes.yaml` (AC 15, per route); one paragraph on when a
  user picks a six-frame database (looking for unannotated ORFs), the `gippa makedb
  --six-frame-orfs` command that built it shown in prose and pointing to
  `../user-guide/makedb.md`; the route E caveats of section 2 (few candidates, no call, the
  floor; IDs in region coordinates; tier `mayu`). "four routes" becomes "five routes" on the page,
  in the quick-start line (`five input routes`) and in the `gippa example` group help (`cli.py`,
  then `docs/reference/cli.md` regenerated). `from-spectra.md`'s "four routes" are another set
  (spectra to `gippa run`) and stay.
- **Tutorial** (`docs/tutorials/example-routes.md`): title `# Tutorial: one Nissle region, five
  routes`; it runs route E too, reads `{out}/E/tables/novel_orfs.tsv` and the report's novel-ORF
  section (aggregates: candidates per class, called 0, the floor), and compares the five runs with
  `gippa compare {out}/A {out}/B {out}/C {out}/D {out}/E --labels A,B,C,D,E --outdir
  {out}/compare` (`--outdir` is required). Still aggregates only (`PEPTIDE_LIKE`); no ORF ID is
  printed on the page (they are row-level content).
- **`docs/release.md`** step 14: "prints the five routes".

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

Markers `e2e`, `slow`, `needs_apptainer`, `needs_bakta_db`, `needs_promotech_models`; the antiSMASH
databases must be installed (`bgc`). Submitted with
`scripts/slurm/example_routes.sbatch` (job `gippa-example-routes`, 4 CPUs, 16 GB, R14), never in
CI or the session allocation. It runs `gippa example fetch --from-file
$GIPPA_ROOT/example/gippa-example-<X.Y>.tar.gz` (so the offline path is exercised), then the
printed commands of each route, and compares each run with `expected/<route>.json` by the
tolerances there. Route B is part of the test only once #536 has landed; until then it is
listed as a known gap in the test's docstring, not skipped silently. It writes the per-route wall times to a TSV next to the outputs; times are
recorded, not asserted (a shared node is noisy; the table in section 3 is filled from them).

*Route E, 2026-11-06:*

- `tests/helpers/example_routes.py`: `ROUTES = ("A", "B", "C", "D", "E")`; `CONVERTED` and
  `ABUNDANCE` gain E (route E converts, and its peptide intensities give `n_abundance` as in A);
  a `NOVEL = frozenset("E")` set adds `n_novel_candidates`, `n_novel_called`, `n_novel_unpaired`
  and `novel_classes` to `required_keys("E")` only, with the tolerances of section 2;
  `measure()` reads them from `tables/novel_orfs.tsv`, `novel/novel_class_fdr.parquet` and
  `run_manifest.json`. `tests/helpers/example_m12_01.py`'s fixture bundle gains
  `E_sixframe/report.parquet` and `expected/E.json`.
- `tests/e2e/test_m12_01_example_routes.py`: `TIERS["E"] = "mayu"` (Q-L10), and for route E the
  AC 22 assertion. The route runs in the same job; `example_routes.sbatch` keeps 4 CPUs, 16 GB and
  its 2 h limit (five routes of about 6 min each fit).

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
| 17 | the slow route test on the published bundle | all four routes exit 0, `status: completed`, every compared key within tolerance; route A `detection_tier = mayu`, route C `two_peptide`, route B `protein_list` (after #536) | e2e (Slurm) |
| 17b | `example_routes.yaml` route commands | every `gippa run` carries the four example settings, `--modules core,bgc`, `--threads 4` | unit |
| 18 | `build_example.py` twice at one commit | same sha256; `--check` passes against `bundles.tsv` | the build script's own tests + the sbatch |
| 19 | `docs/release.md` | has the build, attach and verify steps below | unit (docs) |

*Route E, 2026-11-06:* AC 1 reads "five route headings"; AC 14 adds "route E has `--from diann`,
`--search-db same`, input `{dir}/E_sixframe/report.parquet`, `--out {out}/E/peptides.tsv`, and
the same three-command sequence and settings as A" (AC 17b covers E through the route list); AC 15
adds the route E section; AC 17 reads "all five routes … route E `detection_tier = mayu`". New
criteria:

| # | Input | Expected | Test |
|---|---|---|---|
| 20 | build on a synthetic A4 fixture: rows with (a) region CDS + region ORF, (b) a region ORF only, target and decoy, (c) an ORF crossing the region end, (d) an ORF outside, (e) a `Cont_` ID with a region ORF, (f) a region CDS + an out-of-region ORF | `E_sixframe/report.parquet` keeps exactly (a) and (b), target and decoy alike; all 71 (fixture: all) columns kept | build tests (test-author) |
| 21 | same fixture, the ORF IDs of the kept rows | every ORF ID in `Protein.Ids` and `Protein.Group` is rewritten (rule 9): `parse_orf_id` gives the region contig and `start1 - S - 1`, `end1 - S`, class and strand unchanged; CDS IDs and ID order unchanged; no `NZ_CP058217.1` ORF ID left; `SELECTION.json` lists the columns and the count | build tests |
| 21b | a fixture whose map coordinates disagree with an ID, or whose region translation differs from the `strain_db.faa` record (beyond the first residue) | the build fails naming the ID | build tests |
| 21c | a fixture where route E keeps no decoy, or no ORF row | the build fails (rule 7) naming route E | build tests |
| 22 | the slow route test, route E | exit 0, `status: completed`, tier `mayu`; `n_novel_called = 0`; if `n_novel_candidates > 0`, then `n_novel_unpaired < n_novel_candidates` and `run_manifest.json › novel_orfs › classes` has at least one `picked_decoy` entry (the rewrite makes candidates pair); every compared key of `expected/E.json` within tolerance | e2e (Slurm) |
| 23 | `expected/E.json` (`required_keys("E")`) | holds the base keys including `convert_rows` and `n_abundance`, plus the four novel keys; `required_keys` of A–D unchanged | unit (helpers, `test_m12_01_example_routes_compare.py`) |
| 24 | the tutorial | title "five routes"; mentions `{out}/E` and `novel_orfs.tsv`; exactly one `gippa compare` with `{out}/A … {out}/E`, `--labels A,B,C,D,E` and `--outdir`; no `PEPTIDE_LIKE` match and no `GIPPAORF_` string | unit (docs) |

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
| Linear region of a circular chromosome | written without a topology tag, so it is linear; nothing wraps. A gene that Bakta calls across a region end is `partial` and `partial_excluded`; `expected/` counts them as measured. Preflight's `genome.topology` info row says linear |
| No circular replicon in the bundle | no plasmid (Q-K6); origin-crossing is not exercised by the example (the demo and the unit tests cover it) |
| Region must not cross the `NZ_CP058217.1` origin | the build refuses such a region (it would need an unwrapped interval and a rotated contig) |
| Minus-strand genes at a region end | containment (rule 2) does not depend on strand: a protein on either strand is in the region only if both of its ends are |
| Shared peptides across the boundary | a DIA-NN row is kept only when **every** protein ID is in the region (rule 3), so a peptide shared with an out-of-region paralog is dropped. Within the region, peptides look as unique as they are in the region genome. README caveat |
| Decoy pairing on the region | decoys pair through the peptide vote of `protein_gene_map` (live Bakta writes new locus tags): 4 of 10 pair, so route A is `mayu` with `decoy_pairing_low` (Q-L2 resolved: accepted and stated) |
| Route B before #536 | preflight accepts it; the run writes no gene table. The docs say so until #536 lands; the slow test adds route B then |
| Empty route after subsetting | build fails (rule 6); never shipped |
| Empty `--from-file` | exit 1, size mismatch |
| Empty existing `--outdir` | allowed |
| Two pins with equal version | impossible: the YAML key is unique; the loader rejects two rows with the same `version` (ValueError at import, caught by the unit tests) |
| Region ties in the build | equal scores: the lower start wins |
| Installed GIPPA newer than every pin | newest pin of the same major, with the warning; another major: exit 2 |
| The asset is not uploaded yet (between `release.yml` and step 13b) | 404, exit 2 with "may not be published yet"; step 13b follows the release immediately |
| User has Bakta full DB, or UniOP | the example sets `annotation.bakta.db=light`, so a full-only install fails preflight with the missing light DB (exit 2); with UniOP installed the operon numbers can differ. README says so |
| No antiSMASH databases (`setup --minimal`) | `bgc` fails preflight (exit 2) naming them; the page's prerequisite is `setup --no-genomad` |
| User passes `--strict` | the region's preflight warnings give exit 3; the printed commands never use it |
| Bundle on a compute node without internet | `--from-file`, documented on the page and in the no-network message |
| Disk full in the middle of unpacking | exit 2; staging removed; `DIR` not created |
| *Route E:* ORF overlapping a CDS (`cds_outframe`, `antisense`; 5,124 of the 5,678 region ORFs) | kept by its own coordinates (rule 3E). A row shared by the CDS and the ORF is kept when both are region entries; the peptide places on the CDS gene as in A, and its ORF ID leaves the decoy counts only if it maps to no gene (decision 28). 39 kept rows mix CDS and ORF IDs |
| *Route E:* ORF crossing a region end (9 entries) | not a region entry, so every row naming it is dropped (target or decoy), even when its other IDs are in the region; the region contig cannot hold its interval |
| *Route E:* ORF whose closing stop codon lies outside the region although the ORF is contained | measured 0 on this region. Had there been one: on the linear region it is a 3'-open candidate; `_orf_key` and `_stop` both wrap by the region length, so it still pairs. Not a build failure |
| *Route E:* minus-strand ORF | containment (rule 3E) uses `start - 1 >= S` and `end <= E` on either strand; rule 9 shifts both ends by `S` and keeps the `m` tag; the translation check reverse-complements |
| *Route E:* `nterm_extension` entry (107 in the region) | its ID is rewritten like the others; the run pairs it by locus tag through its own `nterm_entries`, which only works when live Bakta light calls the same gene start as the full-DB annotation A4 was built on. Otherwise the candidate is `no_makedb_pair` and counts in `n_novel_unpaired`. Not checked by translation (upstream part + junction) |
| *Route E:* a decoy whose only IDs are ORF entries mapping to no gene (16 kept) | no gene-level decoy row (decision 28); it still counts as a makedb decoy for its class in the novel-ORF FDR |
| *Route E:* `Cont_` with a region ORF in one cell | dropped (rule 3E), as in A |
| *Route E:* ties | the build has no ranking; rule 3E is a filter. In the run, picked ties go to the decoy (`_class_rank`), unchanged |
| *Route E:* empty route | the build fails (rule 7: no target, no decoy or no ORF row) |
| *Route E:* bundle size | `E_sixframe/report.parquet` ≈ 387 kB (zstd); the bundle ≈ 0.9 MB, under the pin's `size_gib: 1` cap with room |
| *Route E:* a user's own six-frame genome-wide run vs route E | numbers differ (region, live Bakta light, rewritten IDs); the README says not to compare them with the pilot |

## 7. Open questions

### For Khoa

Resolved on 2026-11-03: Q-K1 (licence: "use the data", statement in section 2), Q-K3 (the four
faster settings: "accept"), Q-K6 (no plasmid: the region is a chromosome piece), and #536 ("fix").
Still open, each with its reversible default:

| # | Question | Default |
|---|---|---|
| Q-K2 | Bundle cadence: one per minor (above), or one per release? | one per minor; patches reuse it |
| Q-K4 | Replace the empty `tutorials/ecoli-k12.md` with the example tutorial? | yes; the K-12 tutorial goes to the backlog |
| Q-K5 | The asset upload is a manual step by you (13b), because the bundle is built on FarmShare and `release.yml` cannot see it. Fine? | manual step 13b by Khoa |
| Q-K7 | File paths in the data are shortened to base names. Acceptable? | yes |
| Q-K8 (route E, 2026-11-06) | The licence statement adds `E_sixframe/` and "the identifiers of six-frame ORFs (made by GIPPA, not by the depositors) renumbered to the region's coordinates". Fine? | yes, the amended text of section 2 |

### For the lead

| # | Question | Default |
|---|---|---|
| Q-L1 | ADR-0074 (R1, R10): the command, hosting, the lifted row-level rule, the licence statement | the lead writes it as PR 2; Khoa accepts before PR 5 merges |
| Q-L2 | Route A at `mayu` | **resolved**: accept `mayu` and say so (README, docs page) |
| Q-L7 | **Genome accession.** Khoa's answer names GCA_043228065.1, but the measured region is on `NZ_CP058217.1` of GCF_021559835.1 (GenBank twin GCA_021559835.1), the pilot's pass-3 assembly, and `build_example.py` reads that one. This spec, the licence statement and PROVENANCE use the assembly the data actually comes from | GCF_021559835.1 / GCA_021559835.1; confirm with Khoa that the answer meant "the Nissle genome, INSDC" |
| Q-L8 | `bgc` needs the antiSMASH databases (about 10 GB), so the prerequisite is `gippa setup --no-genomad`, not `--minimal`. Keep `bgc` in the example, or drop it so that `--minimal` is enough? | keep `bgc` (it shows the pks/HPI contrast), and document `--no-genomad` |
| Q-L3 | `_unpack` is private in `gippa.setupdb.run`; make a public alias? | yes, `unpack_archive` with `_unpack` kept |
| Q-L4 | Group rendering in `gen_reference.py` and the `test_m11_08_docs.py` option-drift parser (`## gippa <name>` assumes one word) both need changes for `gippa example fetch` | yes; the test-author updates the parser in PR 4 |
| Q-L5 | The nav and page lists in `test_m11_07_mkdocs_config.py` and `test_m11_08_docs.py` are exact; adding `get-started/example.md` and renaming the tutorial needs test edits | test-author PR 6a, before the docs PR |
| Q-L6 | An M12 index (`docs/specs/M12-index.md`)? | not now; this is M12's first issue |
| Q-L9 (route E) | **ORF IDs vs the region contig.** `gippa run` pairs a candidate with its makedb entry by `(contig, strand, stop)` parsed from the ID, and the bundle's contig is `nissle1917_NZ_CP058217_1_1965000_2165000` with region coordinates. Unchanged A4 IDs pair nothing. Rewrite the IDs in the bundle (rule 9), or ship them unchanged and accept that every candidate is `no_makedb_pair`? (A third way, a run option for a coordinate offset, is a code change outside M12-01.) | rewrite (rule 9); reversible: drop rule 9 and the ORF-ID clause of the licence |
| Q-L10 (route E) | Route E's tier: `mayu` (as A) is a prediction from decision 28 and route A, not a measurement | `TIERS["E"] = "mayu"`; if work-split row 12 measures another tier, validation records it in `expected/E.json`, and the lead amends this spec and `TIERS` in the same PR |
| Q-L11 (route E) | Ship a region subset of A4's `strain_db.faa` or `db_id_map.tsv`? `gippa run` and the report need neither | no; PROVENANCE names them by sha256, and the page shows the makedb command |
| Q-L12 (route E) | ADR-0074 records the lifted row-level rule and the licence statement; route E extends both | the lead adds a dated amendment to ADR-0074 before row 11 merges |
| Q-L13 (route E) | A4 lives at `$GIPPA_ROOT/runs/pilot/sixframe-fa49b8f/`, not under `$GIPPA_ROOT/pilot/` with the other sources | read it in place (it is under `$GIPPA_ROOT`, so the `_rel` guard holds; SELECTION records the sha256s); copy it to `$GIPPA_ROOT/pilot/A4/` later if the runs tree is pruned |

## 8. Work split (≤ ~400 lines each, in order)

| # | PR | Role | Contents | Code / test lines | Depends on |
|---|---|---|---|---|---|
| 1 | `docs(spec): M12-01 example dataset and gippa example` | spec-writer | this file | 0 | — |
| 2 | `docs(adr): ADR-0074 example bundle` | lead | the ADR; Khoa accepts | 0 | 1 |
| 3 | `feat(example): build script` (`feat/example-dataset`) | validation + test-author | `scripts/example/build_example.py`, `scripts/slurm/build_example.sbatch`, `validation/example/bundles.tsv` (header only), tests on a synthetic pilot fixture (subset rules 1–7, determinism, scrubbing) | ≈ 300 / 300 | 1 |
| 4 | `test(example): gippa example fetch` | test-author | `test_m12_01_example_fetch.py` (AC 1–12), `test_m12_01_example_pin.py` (AC 13), `test_m12_01_example_routes_yaml.py` (AC 14), parser update for groups (Q-L4); red | 0 / ≈ 400 | 1 |
| 5 | `feat(cli): gippa example fetch` | implementer | `src/gippa/example.py`, the `example` group in `cli.py`, `DbItem.size_bytes`, `example_routes.yaml`, the `unpack_archive` alias, `gen_reference.py` group rendering, regenerated `cli.md`, changelog fragment | ≈ 250 / 0 | 2, 4 |
| 6 | `test(example): slow route test` | test-author | `tests/e2e/test_m12_01_example_routes.py`, `scripts/slurm/example_routes.sbatch` | ≈ 40 / ≈ 150 | 5 |
| 7 | `chore(example): pin bundle 1.1` | validation | re-measure at the release-prep commit (route B after #536), `expected/`, the `example_1_1` row, the `bundles.tsv` row, the "current minor has a row" assertion (test-author commit) | ≈ 20 / ≈ 10 | 3, 5, 6, #536 |
| 8a | `test(docs): example page and tutorial` | test-author | nav lists, `TUTORIALS`, AC 15, 16, 19; red | 0 / ≈ 80 | 4 |
| 8b | `docs: try the example, tutorial, quick start` | docs | `get-started/example.md`, `tutorials/example-routes.md` (replacing `ecoli-k12.md`), `mkdocs.yml`, quick start, from-spectra and convert links, backlog entry | 0 (prose) | 7, 8a |
| 9 | `docs(release): build and attach the example bundle` | docs | `docs/release.md` steps 6, 13b, 14, 15 | 0 (prose) | 7 |
| 10 | `docs(spec): M12-01 route E, six-frame database` (route E, 2026-11-06) | spec-writer | this amendment | 0 | — |
| 11a | `test(example): route E` | test-author | build tests on a synthetic A4 fixture (AC 20–21c); `tests/helpers/example_routes.py` (E in `ROUTES`, `CONVERTED`, `ABUNDANCE`, the novel keys, AC 23) and its tests; `tests/helpers/example_m12_01.py` fixture (E files); `test_m12_01_example_routes_yaml.py` (AC 14 for E, `ROUTE_IDS`); `test_m12_01_example_docs.py` (AC 15 for E, AC 24, five routes); `tests/e2e/test_m12_01_example_routes.py` (`TIERS["E"]`, AC 22); red | 0 / ≈ 250 | 10 |
| 11b | `feat(example): route E, six-frame database` | implementer | `scripts/example/build_example.py` (sources, rules 3E, 7, 9, SELECTION, counts) and `bundle_text.py` (README, PROVENANCE, LICENSE-DATA, caveats); `src/gippa/data/example_routes.yaml` (route E); the `gippa example` help text in `cli.py` and regenerated `docs/reference/cli.md`; `docs/get-started/example.md`, `docs/tutorials/example-routes.md`, the quick-start line, `docs/release.md` step 14 (the YAML and the page must change together: AC 15 compares them); changelog fragment | ≈ 150 / 0 (+ prose) | 11a, Q-L12 |
| 12 | `chore(example): rebuild and re-pin bundle 1.1 with route E` | validation | `sbatch scripts/slurm/build_example.sbatch`, then `example_routes.sbatch` on the new bundle (five routes), fill `expected/` (E.json included) and rebuild; the `example_1_1` row and the `bundles.tsv` row with the new sha256 and size; the measured route E time in section 3's table and on the docs page; folded into or merged before the release-prep PR #572 | ≈ 20 / 0 | 11b |

*Route E, 2026-11-06:* rows 10–12 come before the 1.1.0rc1 release-prep (#572), which waits for
row 12: the 1.1 bundle and its pin are made once, with route E in it.

PR 7 is part of the 1.1.0 release-prep (step 6), since `expected/` must come from the release
commit. Before PR 7 merges, `gippa example fetch` on a 1.1 dev install exits 2 (`no example bundle
for GIPPA 1.x` while no row exists), which the unit tests cover with monkeypatched rows.
