# Campaign disk and CPU budget (Phase 0)

Derived from `run_selection.tsv` (sizes from the PRIDE / iProX file lists in
`<accession>/files.tsv`) on 2026-11-19. Priority 1 is the first pass of every dataset (one
condition, 3–16 runs); priority 2 is the rest of T. pallidum (the other strains); priority 3 are
reserve proposals, not budgeted. Rows marked needs-decision in `run_selection.tsv` are included
at their default selection.

## Pilot numbers used (handover plan)

| Step | Pilot wall time | CPUs assumed | CPU-h per unit |
|---|---|---|---|
| Bakta 1.12.2, full DB | 1 h 43 min per genome | 16 (`record_bakta_full.sbatch`) | 27.5 per genome |
| DIA-NN 2.3.1, route A (called genes) | 14–42 min per 3 runs (both calls) | 16 (`diann_pilot.sbatch`) | 1.2–3.7 per run |
| DIA-NN 2.3.1, route B (six-frame ORFs) | **unmeasured** | 16 | **unmeasured** (library several times larger; the demo E2 search measures it first) |
| `gippa run` | 1.6–4.5 h per run | 16 (`--threads 14` + head) | 26–72 per run |

Bakta time is taken per genome as measured on the pilot strain; it scales roughly with genome
size (1.1 Mb T. pallidum to 9.1 Mb M. xanthus), so the per-dataset numbers are order-of-magnitude.
The K-12 MG1655 annotation is reused (job 1849096), so the demo needs no Bakta run.

## Per dataset

| Dataset | Accession | Priority | Runs | GB | Bakta genomes | Bakta CPU-h | DIA-NN A CPU-h | gippa runs | gippa CPU-h |
|---|---|---|---|---|---|---|---|---|---|
| E. coli MG1655 demo (E1, E2, E3 placeholder) | PXD058422 | 1 | 4 | 7.2 | 0 | 0 | 5–15 | 3 | 77–216 |
| E. coli NCM3722 (N1, N2) | PXD060434 | 1 | 3 | 2.9 | 1 | 27 | 4–11 | 2 | 51–144 |
| T. pallidum DAL-1 | PXD065022 | 1 | 3 | 11.6 | 1 | 27 | 4–11 | 2 | 51–144 |
| B. thetaiotaomicron VPI-5482 | PXD055395 | 1 | 16 | 10.4 | 1 | 27 | 20–60 | 2 | 51–144 |
| C. crescentus NA1000 | PXD069004 | 1 | 3 | 5.3 | 1 | 27 | 4–11 | 2 | 51–144 |
| V. natriegens ATCC 14048 | PXD049476 | 1 | 4 | 9.6 | 1 | 27 | 5–15 | 2 | 51–144 |
| C. difficile 630Δerm | PXD065346 | 1 | 3 | 15.2 | 1 | 27 | 4–11 | 2 | 51–144 |
| Synechococcus CC9311 | PXD072227 | 1 | 5 | 6.0 | 1 | 27 | 6–19 | 2 | 51–144 |
| M. tuberculosis H37Rv | PXD057335 | 1 | 3 | 9.2 | 1 | 27 | 4–11 | 2 | 51–144 |
| M. xanthus DK1622 | PXD060688 | 1 | 4 | 13.1 | 1 | 27 | 5–15 | 2 | 51–144 |
| H. volcanii H26 | PXD066008 | 1 | 5 | 13.6 | 1 | 27 | 6–19 | 2 | 51–144 |
| S. islandicus E233S | PXD079213 | 1 | 3 | 4.7 | 1 | 27 | 4–11 | 2 | 51–144 |
| T. pallidum, 6 further strains + isolates 171, 176 | PXD065022 | 2 | 24 | 81.6 | 6 | 165 | 30–90 | 12 | 307–864 |
| **Total, priority 1** | | | 56 | **108.8** | 11 | 302 | 70–209 | 25 | 640–1800 |
| **Total, priorities 1 + 2** | | | 80 | **190.4** | 17 | 467 | 100–299 | 37 | 947–2664 |

gippa runs per dataset: route A and route B (2); the demo adds E3 (route A on the knockouts);
NCM3722 is N1 (route A) and N2 (the same report on MG1655). Not counted: route B DIA-NN searches
(unmeasured), the optional second round with `--patch`, the TE7 determinism rerun (one per phase,
one more `gippa run`), Cascadia (arm E4, GPU hours, unmeasured), ThermoRawFileParser conversions
if needed, and the T. pallidum cross-strain runs on DAL-1.

Not yet sized: the E3 knockout runs (PXD066241 at Panorama Public, needs-decision) and the reserve
swaps (priority 3: PXD065462 0.8 GB, PXD073896 4.5 GB, PXD080037 12.9 GB).

## Against the caps

- **Disk (~500 GB).** Raw downloads are 108.8 GB for priority 1 and 190.4 GB with priority 2.
  The pilot fetch keeps 3x the download free (archive + unpacked `.d` + DIA-NN working files),
  i.e. about 330 GB for priority 1, within the cap; priorities 1 + 2 at 3x (about 570 GB) exceed
  it unless the T. pallidum archives are deleted once unpacked and searched, or the strains run
  in batches. Six-frame route-B libraries add an unmeasured amount per dataset.
- **CPU.** Priority 1 is roughly 1,000–2,300 CPU-h before route B, about 8–18 h of the R14 cap
  (128 CPUs) if fully parallel; in practice the 16-CPU jobs run 4–8 at a time.
- **Route B is the unknown.** Six-frame databases are several times larger than the called-gene
  database; DIA-NN library prediction and search time grow with it. Measure on the demo (E2) before
  budgeting route B for Phase 2.
