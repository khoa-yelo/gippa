# Phase 0 check: does DIA-NN 2.3.1 on Linux read Thermo `.raw`?

**Verdict: yes.** The pinned DIA-NN 2.3.1 Academia image reads the Orbitrap Eclipse `.raw` of the
demo dataset directly on Linux, both in conversion mode and in a library-free search. The image
ships Thermo's RawFileReader (`ThermoFisher.CommonCore.RawFileReader.dll`, `RawWrapper.dll`) and
the .NET 8.0.22 runtime (`/usr/share/dotnet`). ThermoRawFileParser (ADR-0077 decision 2) is
**not** needed for the DIA-NN searches; it is still needed for Cascadia's mzML input (arm E4).
Nothing was built.

## Input

| Item | Value |
|---|---|
| File | `231128_2465_EcUM2_FTFT_Pepsep25_90min_4ul_DIA_2.raw` (PXD058422; E1 run: wild type, LB, S-Trap, DIA; the smallest of the four selected E1 runs) |
| Source | `https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/09/PXD058422/…` (TLS verified) |
| Size | 1,594,236,621 bytes (matches the PRIDE file list) |
| Checksum | SHA-1 `c514bab65061425aab82d1dff56ae5d27aa4675d` from the deposit's `checksum.txt`: verified |
| sha256 | `74d304a3bcde67ae9623b4d80b49c164a6cccac53fa419254cdbeb7e2b78fe94` |
| Location | `$GIPPA_ROOT/campaign/ecoli_k12/raw/` (not in git) |
| Image | `$GIPPA_ROOT/images/diann_2.3.1.sif`, sha256 `1946452d…a6381` (= `validation/pilot/diann_image.sha256`), read-only copy |
| Proteins (search only) | K-12 Bakta 1.12.2 recording `k12_mg1655.faa` (job 1849096) + `$GIPPA_DB_DIR/contaminants/contaminants.fasta` |

## Jobs

| Job | What | Result |
|---|---|---|
| 1998296 | download, two attempts | failed: `curl: (56) OpenSSL SSL_read: unexpected eof` twice; the script discarded the partial file between attempts (bug, fixed) |
| 1998349 | download, two resumed attempts | failed at 1.48 of 1.59 GB (same TLS EOF); partial kept |
| 1998427 | resume of the last 0.12 GB, then the DIA-NN steps (`--qos=dev`, 8 CPUs, 30 GB) | **COMPLETED, exit 0**, 12 min 25 s on wheat-08 |

The PRIDE HTTPS server drops long transfers every 0.2–0.9 GB at 0.5–1 MB/s here; the campaign
downloads must resume partial files across attempts (PROGRESS.md › Blocked).

## Commands (from `validation/campaign/slurm/raw_check.sbatch`)

```bash
sbatch --qos=dev --mem=30000M --exclude=barley-01,barley-02,barley-03,barley-04,rye-01,rye-02 \
  validation/campaign/slurm/raw_check.sbatch 231128_2465_EcUM2_FTFT_Pepsep25_90min_4ul_DIA_2.raw

# step 2: conversion mode (reads the .raw, writes <file>.raw.dia)
apptainer exec --bind $GIPPA_ROOT $GIPPA_ROOT/images/diann_2.3.1.sif /diann-2.3.1/diann-linux \
  --f $GIPPA_ROOT/campaign/ecoli_k12/raw/231128_2465_EcUM2_FTFT_Pepsep25_90min_4ul_DIA_2.raw \
  --convert --threads 8 --verbose 1
# step 3a: library prediction (diann_A2.cfg [step1] flags)
apptainer exec … diann-linux --fasta-search --gen-spec-lib --predictor --fasta k12_mg1655.faa \
  --fasta contaminants.fasta --cut 'K*,R*' --missed-cleavages 1 --min-pep-len 7 --max-pep-len 30 \
  --unimod4 --var-mods 1 --var-mod UniMod:35,15.994915,M --min-pr-charge 1 --max-pr-charge 4 \
  --min-pr-mz 300 --max-pr-mz 1800 --min-fr-mz 200 --max-fr-mz 1800 --peptidoforms \
  --out-lib $GIPPA_ROOT/campaign/ecoli_k12/raw_check/k12.speclib --threads 8 --verbose 1
# step 3b: search of the .raw (diann_A2.cfg [step2] flags)
apptainer exec … diann-linux --f …DIA_2.raw --lib k12.predicted.speclib --fasta k12_mg1655.faa \
  --fasta contaminants.fasta --qvalue 0.01 --report-decoys --matrices --peptidoforms \
  --out $GIPPA_ROOT/campaign/ecoli_k12/raw_check/report.parquet --threads 8 --verbose 1
```

## Exit codes and log excerpts (aggregates only)

| Step | Exit | Wall | Log |
|---|---|---|---|
| 2, `--convert` | 0 | 6 s | `1 files will be processed` · `[0:00] Loading run …DIA_2.raw` · `Finished`; wrote `…DIA_2.raw.dia` (190,399,160 bytes) |
| 3a, predict | 0 | 3 min 38 s | `743434 precursors generated` · `Saving the library to …k12.predicted.speclib` |
| 3b, search | 0 | 5 min 6 s | see below |

Search log (DIA-NN 2.3.1 Academia, 8 threads):

```
[0:04] Loading run …/231128_2465_EcUM2_FTFT_Pepsep25_90min_4ul_DIA_2.raw
[0:09] 1438 MS1 and 57520 MS2 scans in 1438 (inferred) and 1438 (encoded) cycles, 679152 precursors in range
[0:29] Recommended MS1 mass accuracy setting: 2.6 ppm
[0:52] Optimised mass accuracy: 8 ppm
[5:02] Number of IDs at 0.01 FDR: 37578
[5:02] Number of genes identified at 1% FDR: 1229 (precursor-level), 1157 (protein-level)
[5:06] Protein groups with global q-value <= 0.01: 2460
Finished
```

Only these counts are recorded; the report, matrices and library stay under
`$GIPPA_ROOT/campaign/ecoli_k12/raw_check/` and are not committed.

## Consequences for the campaign

- The eight Thermo datasets can be searched from `.raw` with the pinned image, as the Bruker
  datasets are from `.d`; no conversion step and no new tool for the DIA-NN routes.
- The pilot DIA-NN wrapper (`diann_pilot.sbatch`) only accepts `.d` runs from an
  `acquisition.json`; the campaign wrapper must accept `.raw` (and the `.d.zip`/`.d.tar`/`.d.Zip`
  archive variants seen in the inventory).
- Throughput on one 90-min Eclipse run, 8 threads: 5 min search after a 3.6 min library
  prediction for the K-12 Bakta proteins plus contaminants (4,740 library isoforms; route A size). Route B (six-frame) is unmeasured.
- Arm E4 (Cascadia) still needs mzML: ThermoRawFileParser per ADR-0077 when E4 starts.
