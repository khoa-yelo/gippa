# Quick-start example dataset (draft, targets 1.1.0)

A cut-down subset of the M10 pilot (`validation/pilot/README.md`, pass 3): *E. coli* Nissle 1917
(GCF_021559835.1, chromosome `NZ_CP058217.1`) and PXD078459 runs WT_1 to WT_3. Nothing is
re-searched: `build_example.py` subsets the pilot's DIA-NN 2.3.1 A2 report and the Cascadia `.ssl`
files to one contiguous region. The selection rule is in that script's docstring and is copied
into every bundle's `SELECTION.json`, with the sha256 of each source and output file.

**FDR caveat.** The q-values, PEPs and MaxLFQ values in routes A, B and D come from DIA-NN. DIA-NN
computed them on the whole-proteome search (4,685 strain proteins plus contaminants, three runs).
The subset does not recompute them.

| Script | Purpose |
|---|---|
| `build_example.py` | writes `$GIPPA_ROOT/example/<name>/`: genome, routes A to D, `SELECTION.json`, `bundles/*.tar.gz` |
| `../slurm/example_run.sbatch` | converts one route, then runs one `gippa run` on the local executor in its allocation (`farmshare_local,apptainer`) |
| `summarize_runs.py` | prints wall time, per-process trace times and result counts for each run |

## Recommended configuration (measured 2026-11-03)

- Region: `NZ_CP058217.1:1965001-2165000` (200 kb; `--start 1965000 --end 2165000`). It holds
  the colibactin *pks* island and the yersiniabactin HPI. antiSMASH calls both as one
  NRPS/T1PKS/NRP-metallophore region over the first 117 kb.
- Run options:
  `--modules core,bgc --set annotation.bakta.db=light --terminators transtermhp --set annotation.promoter_calculator.enabled=false`
- CPUs: 8 (Slurm granted 9). The run takes about 5 min wall.

Each non-default option replaces one step whose time does not fall when the region gets smaller:

| Default | Measured cost (wheat nodes) | Replaced by |
|---|---|---|
| Bakta full DB (`psc.dmnd`, 46 GB) | 26 min on 55 kb: DIAMOND reads the whole database | Bakta light DB (1.3 GB): about 1 min |
| TerminatorNet | 4.5 to 4.8 min per strand of 100 kb on 4 CPUs, which sets the critical path | TransTermHP: 1 s |
| geNOMAD (`mobile`) | 4.6 to 5.4 min, apparently fixed (100 kb) | module off |
| Promoter Calculator | 100 s on 100 kb, 428 s on 200 kb (one task) | switched off; Promotech still runs |
| `function` | 13 s, but the Bakta light DB carries no COG, GO or KEGG IDs, so it gives 0 rows | module off |

Bakta and antiSMASH stay live, and Bakta's Pyrodigal warns when a region is under 100 kb. The
region contains no prophage: geNOMAD is too slow to include.

## Measured runs (route A unless noted; local executor; wheat nodes)

| Region | Options | CPUs (granted) | Slurm elapsed | Nextflow span | Top steps (s) |
|---|---|---|---|---|---|
| 55 kb | core, full DB, defaults | 9 | (failed at 26:51, PATH bug in the first sbatch) | n/a | Bakta 1563 |
| 55 kb | core, light | 9 | 7:33 | 6:39 | TerminatorNet 147/strand, Promotech 79, MAP_TOLERANT 78, Bakta 65 |
| 100 kb | core, light | 9 | 9:28 | 8:51 | TerminatorNet 285/strand, Promoter Calculator 99, MAP_TOLERANT 96 |
| 100 kb | default modules, light | 9 | 14:14 | 13:45 | geNOMAD 275, TerminatorNet 268/strand, Promoter Calculator 96 |
| 100 kb | core, light, TransTermHP | 9 | 5:17 | 4:43 | Promoter Calculator 103, MAP_TOLERANT 63, Bakta 48 |
| 100 kb | core, light, TransTermHP | 5 | 7:34 | 6:29 | Promoter Calculator 134, Promotech 92, Bakta 89 |
| 100 kb | default modules, light, TransTermHP | 9 | 10:58 | 10:18 | geNOMAD 326, Promoter Calculator 100 |
| 100 kb | core,bgc,function, light, TransTermHP | 9 | 5:56 | 5:16 | Promoter Calculator 109, MAP_TOLERANT 72, Promotech 69 |
| 200 kb | core, light, TransTermHP | 9 | 11:07 | 10:29 | Promoter Calculator 428 |
| 200 kb | core,bgc,function, light, TransTermHP, no PC | 9 | 5:35 | 5:02 | MAP_TOLERANT 68, Promotech 63, Bakta 57 |
| **200 kb** | **core,bgc, light, TransTermHP, no PC** | 9 | **5:46** | **5:05** | MAP_TOLERANT 72, Promotech 65, Bakta 62, plots 31, antiSMASH 30 |
| **200 kb** | same | 5 | **5:47** | **5:00** | Bakta 64, MAP_TOLERANT 61, Promotech 59 |
| 200 kb | same, route C (Cascadia) | 9 | 5:45 | 5:07 | |
| 200 kb | same, route D (generic) | 9 | 5:41 | 4:59 | |
| 200 kb | same, route B (proteins + abundance) | 9 | 3:41 | 3:08 | annotation and report only (below) |

With the recommended options, about 35 tasks run in sequence along the critical path. Each
`gippa _step` task takes about 10 s, mostly start-up. As a result, 4 and 8 CPUs give the same
time, and so do 100 and 200 kb.

## Result sanity at 200 kb (route A)

- Genes and detection: 158 genes and 99 detected. The tier is `mayu`, with gene-level FDR 0.012 at
  q 0.01.
- Decoys: 37 decoy rows survive, but only 4 of 10 decoy peptides pair with a gene. The cause is
  that live Bakta writes new locus tags, so pairing has to go through peptide votes; the run
  therefore falls back from `decoy_picked`. Adding the search FASTA (`--proteins`) pairs every
  decoy (tested at 100 kb). The `decoy_picked` q floor is then 1/24 > 0.01, and 0 genes are
  detected. With this few decoys, `mayu` is the honest tier.
- BGC contrast: bgc_1 covers 71 genes, with 30 detected against 42.4 expected.
- Operons: 46 operons and 4 significant operon flags.
- Empty or not estimable: islands (0), enrichment and GSEA (light DB, no categories),
  autocorrelation (`not_estimable`), variants (0) and novel ORFs (0).
- Route C: 76 Cascadia rows match the region exactly. Only 4 genes are detected (`two_peptide`
  tier, `c_assumed`), so this route is a format demo.
- Route B: GIPPA 1.0.0 runs only ANNOTATE and REPORT when no peptide table is given
  (`workflow/main.nf` gates MAPPING on `inputs.peptides`). Proteins with an abundance table and no
  peptides therefore produce no gene table, although preflight accepts them.
