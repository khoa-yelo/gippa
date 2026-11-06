# Tutorial: one Nissle region, five routes

This tutorial runs the five routes of the example on one 200 kb region of the *E. coli* Nissle
1917 genome. Each route feeds a different kind of input into GIPPA. You then read `genes.tsv` and
the report of each run, see how the detection tier changes with the input, and compare the five
runs with `gippa compare`.

The bundle, its licence terms and the prerequisites are described in
[Try the example](../get-started/example.md). Read that page first. This tutorial quotes only
counts and tiers, never table rows.

## 1. Fetch the bundle

```bash
gippa example fetch
```

This unpacks the bundle into `gippa-example` in the current directory. The commands below use that
path. Every run writes to its own directory under `gippa-example/runs`.

## 2. Run the five routes

Run each route's commands in order. The routes do not depend on each other. Each `gippa run` takes
about 5 min on 4 CPUs (route B about 8.5 min) after `gippa setup`.

### Route A: a DIA-NN report

The DIA-NN report is converted to GIPPA's peptide table, checked and run. The output goes to
`gippa-example/runs/A`.

```bash
gippa convert --from diann gippa-example/A_diann/report.parquet --out gippa-example/runs/A/peptides.tsv --search-db same
gippa preflight --genome gippa-example/genome/nissle_r200k.fna.gz --peptides gippa-example/runs/A/peptides.tsv --set 'modules=[core,bgc]' --set annotation.bakta.db=light --set annotation.terminators.tool=transtermhp --set annotation.promoter_calculator.enabled=false
gippa run --genome gippa-example/genome/nissle_r200k.fna.gz --peptides gippa-example/runs/A/peptides.tsv --genome-relation same --outdir gippa-example/runs/A --modules core,bgc --set annotation.bakta.db=light --terminators transtermhp --set annotation.promoter_calculator.enabled=false --threads 4
```

### Route B: a protein list with abundances

Route B has no peptides: a protein FASTA and a MaxLFQ abundance table. There is nothing to
convert. The output goes to `gippa-example/runs/B`.

```bash
gippa preflight --genome gippa-example/genome/nissle_r200k.fna.gz --proteins gippa-example/B_proteins/proteins.faa --abundance gippa-example/B_proteins/abundance.tsv --set 'modules=[core,bgc]' --set annotation.bakta.db=light --set annotation.terminators.tool=transtermhp --set annotation.promoter_calculator.enabled=false
gippa run --genome gippa-example/genome/nissle_r200k.fna.gz --proteins gippa-example/B_proteins/proteins.faa --abundance gippa-example/B_proteins/abundance.tsv --abundance-type maxlfq --genome-relation same --outdir gippa-example/runs/B --modules core,bgc --set annotation.bakta.db=light --terminators transtermhp --set annotation.promoter_calculator.enabled=false --threads 4
```

### Route C: Cascadia de novo results

The three Cascadia files are converted together. The output goes to `gippa-example/runs/C`.

```bash
gippa convert --from cascadia gippa-example/C_cascadia/WT_1.ssl gippa-example/C_cascadia/WT_2.ssl gippa-example/C_cascadia/WT_3.ssl --out gippa-example/runs/C/peptides.tsv
gippa preflight --genome gippa-example/genome/nissle_r200k.fna.gz --peptides gippa-example/runs/C/peptides.tsv --set 'modules=[core,bgc]' --set annotation.bakta.db=light --set annotation.terminators.tool=transtermhp --set annotation.promoter_calculator.enabled=false
gippa run --genome gippa-example/genome/nissle_r200k.fna.gz --peptides gippa-example/runs/C/peptides.tsv --genome-relation same --outdir gippa-example/runs/C --modules core,bgc --set annotation.bakta.db=light --terminators transtermhp --set annotation.promoter_calculator.enabled=false --threads 4
```

### Route D: a generic peptide table

Route D's table was made from route A's rows, in the form any other search engine can write. The
output goes to `gippa-example/runs/D`.

```bash
gippa convert --from generic gippa-example/D_generic/peptides.tsv --out gippa-example/runs/D/peptides.tsv --source database --search-db same
gippa preflight --genome gippa-example/genome/nissle_r200k.fna.gz --peptides gippa-example/runs/D/peptides.tsv --set 'modules=[core,bgc]' --set annotation.bakta.db=light --set annotation.terminators.tool=transtermhp --set annotation.promoter_calculator.enabled=false
gippa run --genome gippa-example/genome/nissle_r200k.fna.gz --peptides gippa-example/runs/D/peptides.tsv --genome-relation same --outdir gippa-example/runs/D --modules core,bgc --set annotation.bakta.db=light --terminators transtermhp --set annotation.promoter_calculator.enabled=false --threads 4
```

### Route E: a DIA-NN report from a six-frame database

The same runs as route A, searched against a database that also holds every six-frame ORF of the
genome (`gippa makedb --six-frame-orfs`). It is converted and run exactly like route A. The output
goes to `gippa-example/runs/E`.

```bash
gippa convert --from diann gippa-example/E_sixframe/report.parquet --out gippa-example/runs/E/peptides.tsv --search-db same
gippa preflight --genome gippa-example/genome/nissle_r200k.fna.gz --peptides gippa-example/runs/E/peptides.tsv --set 'modules=[core,bgc]' --set annotation.bakta.db=light --set annotation.terminators.tool=transtermhp --set annotation.promoter_calculator.enabled=false
gippa run --genome gippa-example/genome/nissle_r200k.fna.gz --peptides gippa-example/runs/E/peptides.tsv --genome-relation same --outdir gippa-example/runs/E --modules core,bgc --set annotation.bakta.db=light --terminators transtermhp --set annotation.promoter_calculator.enabled=false --threads 4
```

## 3. Read `genes.tsv` per route

Each run writes one row per gene to `tables/genes.tsv`. Three columns matter here:

- `detection_tier`: the tier the run used (the same on every row);
- `detected`: whether the gene is detected in that tier;
- `detected_2pep`: whether the gene has at least two distinct accepted peptides, in every tier.

This prints the tier and the counts of each run:

```bash
python3 - <<'EOF'
import csv

for route in "ABCDE":
    path = f"gippa-example/runs/{route}/tables/genes.tsv"
    with open(path, encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    detected = sum(r["detected"] == "true" for r in rows)
    two_pep = sum(r["detected_2pep"] == "true" for r in rows)
    print(route, rows[0]["detection_tier"], len(rows), detected, two_pep)
EOF
```

The tier is also in `run_manifest.json` (`detection_tier`). Compare your numbers with
`gippa-example/expected/README.md`, which lists each route's numbers and their tolerances.

## 4. See how the tier changes with the input

GIPPA picks the strongest tier the input allows
([Detection tiers](../user-guide/interpreting-results.md#detection-tiers)). On this region:

| Route | Input | Tier | Why |
|---|---|---|---|
| A | DIA-NN report | `mayu` | the peptides carry q-values; too few decoys pair with a gene for `decoy_picked` (warning `decoy_pairing_low`) |
| B | protein FASTA and abundance table | `protein_list` | proteins only, no peptides |
| C | Cascadia results | `two_peptide` | peptide sequences only: a gene needs two distinct peptides |
| D | generic peptide table | `mayu` | the same rows and decoys as route A |
| E | DIA-NN report, six-frame database | `mayu` | as route A: the six-frame decoys pair with no gene |

- **A and D.** On the whole genome, DIA-NN's decoys would give `decoy_picked`. A 200 kb region
  keeps too few of them, so GIPPA warns and falls back to `mayu`. The gene q-value is still
  estimated. Route A detects 99 of 158 genes at a gene-level FDR of 0.012.
- **B.** A gene is detected when an input protein maps to it. Route B detects 106 of 158 genes.
  The gene error is `unestimated`: it comes from the search engine, not from GIPPA.
- **C.** A gene with one peptide is `single_peptide_unestimated` and counts as not detected.
  Route C detects 4 genes. It is a format demo, not a result: the bundle keeps only Cascadia rows
  that match the region exactly.

These numbers are for a region. Peptides shared with genes outside it look unique here. Do not
compare them with whole-genome runs.

### Route E's novel ORFs

Route E is the one route whose peptides can name six-frame ORFs. Its candidates are in
`gippa-example/runs/E/tables/novel_orfs.tsv`, one row per candidate, with its class and
`called`. This counts them per class:

```bash
python3 - <<'EOF'
import csv
from collections import Counter

path = "gippa-example/runs/E/tables/novel_orfs.tsv"
with open(path, encoding="utf-8") as fh:
    rows = list(csv.DictReader((ln for ln in fh if not ln.startswith("#")), delimiter="\t"))
print(len(rows), "candidates;", sum(r["called"] == "true" for r in rows), "called")
print(Counter(r["class"] for r in rows))
EOF
```

Expect a few candidates and 0 called. A class's q-value cannot go below a floor of one over its
number of target wins, so a floor of 0.01 needs at least 100 wins; the region holds a few dozen
six-frame entries in all. `run_manifest.json` (`novel_orfs` › `classes`) lists each class with its
regime (`picked_decoy`), its counts, its `q_floor` and `floor_above_alpha`. The report's novel-ORF
section shows the same class summaries.

## 5. Read the report per route

Open `report/report.html` of each run in a browser, for example
`gippa-example/runs/A/report/report.html`. It is one self-contained file.

- **Summary** has a genes card: the number of genes detected, out of all genes, and how the
  tier estimates the error. The `mayu` runs (A and D) give the Mayu gene-level FDR and the count
  with at least two distinct peptides. Route B reads "protein list: no FDR estimated" and route C
  "two-peptide rule: no FDR estimated". Routes A, C and D also have a card for the peptides
  placed; route B has one for the proteins.
- **Gene detection** shows only the columns the tier estimates. Route B (`protein_list`) shows
  the detected count only. Route C (`two_peptide`) also shows the count detected with two
  peptides.
- **Operons** and the gene cluster (module `bgc`) do not depend on the route.
- Several sections are empty on this region: islands, variants, novel ORFs (route E has a few
  candidates and no call), enrichment and autocorrelation. The [example page](../get-started/example.md#what-results-to-expect) says why.

See [Plots and report](../user-guide/plots-and-report.md) for every section.

## 6. Compare the five runs

The five runs share one genome and one gene model (the same Bakta database and settings), so
`gippa compare` can put them side by side:

```bash
gippa compare gippa-example/runs/A gippa-example/runs/B gippa-example/runs/C gippa-example/runs/D gippa-example/runs/E --labels A,B,C,D,E --outdir gippa-example/runs/compare
```

It writes to `gippa-example/runs/compare`:

- `tables/compare_detection.tsv`: which genes each run detects, and how many genes are detected in
  exactly each set of runs;
- `tables/compare_abundance.tsv`: abundance agreement for each pair of runs, on the genes with a
  measured abundance in both. Route C has no abundance;
- `tables/compare_stability.tsv`: the overlap of flagged sets, such as operon flags;
- `report/compare.html`: the comparison page with plot V16.

Route D was made from route A's rows, so the two runs are a check on the generic converter.
Route E searched the same spectra against a much larger database, so A against E shows what the
six-frame search changes at the gene level. See
[Comparing runs](../user-guide/compare.md) for the columns.

## Next steps

- Make such inputs from your own spectra: [From raw spectra](../get-started/from-spectra.md).
- Run a whole genome with the defaults: [Quick start](../get-started/quickstart.md).
