# ADR-0074: A real example dataset, fetched by `gippa example`, from a subset of the pilot
- Status: accepted-by-khoa (2026-11-03 and 2026-11-04, in session: "Nissle pilot subset",
  "GitHub release asset", "Real tools only, small region", "accept" (the four faster settings),
  "use the data", and "yes" to the M12-01 spec)
- Date: 2026-11-04
- Spec sections: §13 (M10 pilot data), §14.9; `docs/specs/M12-01-example-dataset.md`
- Needs Khoa: yes (a new command, a new distribution channel for data, and an exception to the
  pilot data rule)

## Context
GIPPA 1.0.0 ships `gippa demo` (synthetic SIM100K, precomputed layers) but no real input a user
can run through `gippa run`. Khoa asked for "a real but subset dataset that user can run on and
get result quick, say within 5 mins", showing the routes through DIA-NN output, protein FASTA plus
abundance, de novo output and a generic table.

Three standing rules are touched:
- the pilot data rule (PROGRESS.md › Answered by Khoa, 2026-10-24): no row-level DIA-NN content is
  committed or distributed;
- R6 (official sources only) and the wheel's contents (code and pinned metadata, no data);
- the source licence: PXD078459 (iProX IPX0017169000, public since 2026-05-16) states no licence,
  and its publication is pending.

## Decision
1. **Data.** The example is a 200 kb region of the E. coli Nissle 1917 chromosome,
   `NZ_CP058217.1:1965001-2165000` (assembly GCF_021559835.1, the pilot's pass-3 genome; INSDC),
   with the pilot's own DIA-NN 2.3.1 and Cascadia results for PXD078459 WT_1–WT_3, cut to that
   region. Nothing is re-searched. The pilot data rule is lifted **for this subset only**.
2. **Routes.** A: DIA-NN `report.parquet`; B: protein FASTA plus abundance TSV (needs #536,
   merged as #537); C: Cascadia `.ssl` (a format demo: 4 genes detected); D: a generic peptide TSV.
3. **Real tools only.** No precomputed annotation layers. The example command uses four
   non-default settings that keep every step real but fit the time target:
   `--modules core,bgc --set annotation.bakta.db=light --terminators transtermhp
   --set annotation.promoter_calculator.enabled=false`, 4 threads. Measured: 5:00–5:07 of
   Nextflow time for routes A, C and D on a warm install. The README states that the defaults are
   slower and gives the measured default times.
4. **Licence statement.** The code is MIT. For the PXD078459-derived files no licence is granted;
   they are provided for demonstration and testing, with citation of PXD078459 and iProX
   IPX0017169000, and are withdrawn on request. The genome region is INSDC data.
5. **Distribution.** One bundle per GIPPA minor version, attached by Khoa to the GitHub release
   `vX.Y.0` (patch releases reuse it). It is pinned by sha256 in `src/gippa/data/db_sources.yaml`
   and fetched by `gippa example fetch [--outdir DIR] [--version X.Y] [--from-file FILE]`, which
   reuses the `gippa setup` downloader and verifier. The bundle is not in the wheel and not on
   Zenodo.
6. **Honest results.** The bundle's `expected/` records each route's tier and key counts. Route
   A and D report `mayu` (too few decoys pair in a 200 kb region for `decoy_picked`); the README
   says so, and that the q-values come from the whole-proteome search.

## Evidence (before → after: metrics, tests, simulations)
- Feasibility runs (Slurm 1930684–1931447, `scripts/example/README.md` on branch
  `feat/example-dataset`): 200 kb, recommended settings, route A 5:05 (9 CPUs) and 5:00 (5 CPUs);
  C 5:07; D 4:59. With defaults: 13:45 at 100 kb (all modules), 10:29 at 200 kb (Promoter
  Calculator on), and the full Bakta DB alone 26 min at 55 kb.
- Route A at 200 kb: 158 genes, 99 detected, gene-level FDR 0.012; the BGC contrast (30 of 71
  cluster genes detected against 42.4 expected) and 46 operons (4 flagged).
- Route B on the #536 fix (Slurm 1931828): exit 0, 8:31, 158 genes, 106 detected, `protein_list`.
- Bundle: 0.48 MB for all routes.
- PXD078459: ProteomeXchange record announced 2026-05-16; raw files download without a login.

## Consequences
- A new command (`gippa example`) and a new release step: Khoa uploads the bundle with
  `gh release upload`, because the bundle is built on FarmShare from pilot results that CI cannot
  see.
- The example does not exercise geNOMAD, TerminatorNet, Promoter Calculator, the full Bakta DB or
  the `function` module; the default pipeline stays covered by the weekly end-to-end suite.
- If the PXD078459 authors ask, the bundle is withdrawn and replaced with a public K-12 dataset
  (backlog).

## Amendment, 2026-11-06: route E (six-frame database)

Khoa (in session, 2026-11-06): "we should add the 6 frame db as part of the suite as well, as that
is one of the routes mostly use by the user". The bundle gains route E, a region subset of the
pilot's six-frame DIA-NN search (arm A4, the same PXD078459 runs WT_1–WT_3, DIA-NN 2.3.1 against
`gippa makedb --six-frame-orfs`). Khoa extended the row-level lift to A4's rows inside the region;
everything else of the pilot stays aggregate-only. Six-frame ORF identifiers, which GIPPA made and
the depositors did not, are renumbered to the region's coordinates so that `gippa run` pairs them
with the bundle genome (spec M12-01 rule 9). Route E demonstrates the six-frame workflow; with 35
distinct ORF entries in the region no novel ORF can reach q ≤ 0.01, so `expected/E.json` pins 0
calls. The 1.1 bundle is rebuilt and re-pinned before 1.1.0rc1.
