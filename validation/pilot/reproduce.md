<!-- reproduce: status=final freeze_commit=fa49b8f145ea4598a5f896f9d0b1956912fc5df3 a3_branch=cascadia -->
# Pilot reproduction runbook (GIPPA_SPEC.md §15.2 steps 4-9; M10-14)

This is the recorded command (§3.9) for every step of the pilot, in execution order: the inputs
it reads, the job it runs as and the job that ran it. The format is fixed by
docs/specs/M10-14-pilot-reproducibility.md ("The runbook format"), and
`validation/pilot/reproduce_check.py static` checks it on every commit
(`tests/unit/test_m10_14_reproduce.py`). **Status: final.** The freeze commit F is in the header. S1-S7 carry their earlier evidence; S8-S14 ran at F, except the entrapment rerun, which ran from main (Evidence table).
The freeze uses the pass-3 wild-type reference (Khoa, 2026-10-25; README "Pass-3 runs"), not the
committed chassis genome the M10-14 commands name. A pinned detached worktree held F, so the main
checkout stayed free to move (README "Pass-3 runs" did the same). R11 audit: `tests/unit/test_m10_14_r11_audit.py`
passes on this tree.

Every command runs on FarmShare with the environment sourced, from the worktree pinned at F:

```bash
export GIPPA_ROOT=/home/classes/bios/270/khoa/giap/local
export GIPPA_REPO=$GIPPA_ROOT/worktrees/pilot-fa49b8f145ea   # git worktree add --detach, at F
source "$GIPPA_REPO/scripts/farmshare_env.sh"
X="--exclude=barley-01,rye-02 --export=ALL,XDG_CACHE_HOME=$GIPPA_ROOT/.cache,GIPPA_REPO=$GIPPA_REPO,GIPPA_ROOT=$GIPPA_ROOT"
PY="$UV_PROJECT_ENVIRONMENT/bin/python"
P=$GIPPA_ROOT/pilot/inputs
F7=$(git -C "$GIPPA_REPO" rev-parse --short=7 HEAD)
R=$GIPPA_ROOT/runs/pilot/freeze/$F7
```

The job classes are the scripts' `#SBATCH` headers unless the command overrides them. No step
runs on a login node, and none runs `uv sync` or `pip`.

## Inputs

A `.d` run has no single hash, so each run is listed twice: through its archive (the hash in the
committed manifest; the archives were deleted after unpacking, so `verify` reports them `absent`)
and through its `analysis.tdf` (M10-02 `provenance.json`). `analysis.tdf_bin` is not hashed (L11).

| id | path | sha256 | recorded_in |
|---|---|---|---|
| genome | tests/data/genomes/nissle1917.fna.gz | fe74801661e688a525eef7a7a6be6210659ebc9a2e61691efb1f10b65ea623ee | tests/data/genomes/PROVENANCE.md;tests/data/recorded/bakta/1.12.2-full/nissle/PROVENANCE.md |
| bakta_json | tests/data/recorded/bakta/1.12.2-full/nissle/nissle1917.json.gz | a012798400eaa71587dc7954473a201ef046b4907bdfa31936a2927a1cb6c97c | validation/pilot/db_a2/makedb_manifest.json;tests/data/recorded/bakta/1.12.2-full/nissle/PROVENANCE.md |
| strain_db | validation/pilot/db_a2/strain_db.faa.gz | 80dd01f3eeacf209131c7e14ebc65c12174aad89d9a9fd5a285704843faa887e | $GIPPA_ROOT/pilot/A2/provenance.json |
| cfg_a1 | validation/pilot/diann_A1.cfg | 0546a8905bcf092675991a793a45ba46f8e93a875834546587c07cc7d97ae608 | $GIPPA_ROOT/pilot/A1/provenance.json |
| cfg_a2 | validation/pilot/diann_A2.cfg | c0d145e87e582e450295517ac393d2c86093587605dd3da66c6f6a923c1f98e4 | $GIPPA_ROOT/pilot/A2/provenance.json |
| acquisition | validation/pilot/inputs/acquisition.json | a7af58ae2a31392f199dbca2259b899ab75f4a613ad020005ccf322e286bacb3 | $GIPPA_ROOT/pilot/A1/provenance.json;$GIPPA_ROOT/pilot/A2/provenance.json |
| proteome_a1 | $GIPPA_DATA_DIR/proteomes/UP000000625.fasta | a8052a607a98d2444c9d5dd685bd1370150df8c5570c3116811544bbf3d093a2 | $GIPPA_ROOT/pilot/A1/provenance.json |
| contaminants | $GIPPA_DB_DIR/contaminants/contaminants.fasta | c0782146be6dad6ecf53e23e5a21a7ed4fea9b4802cb01e8d00091404ec66926 | $GIPPA_ROOT/pilot/A1/provenance.json;$GIPPA_ROOT/pilot/A2/provenance.json |
| diann_image | $GIPPA_ROOT/images/diann_2.3.1.sif | 1946452d56d2fb0f66e05ced473df88e53484e28040f1ecbd83c4a43405a6381 | validation/pilot/diann_image.sha256;validation/pilot/README.md |
| px_record | $GIPPA_DATA_DIR/pxd078459/record/PXD078459.xml | 5d3f16c3137fb11ecae7f79ad02063d8e3008b21380784fcfc8aa77b7fb533f9 | validation/pilot/README.md;$GIPPA_DATA_DIR/pxd078459/record/PXD078459.xml.info |
| wt1_zip | $GIPPA_DATA_DIR/pxd078459/MLA1914689_WT_P1_1_Y-D7_1_13119.d.zip | 58d09a41d6642f35f0cb25e88e4a715fbfe0c30f7b4f839c9fd30601bd99f013 | validation/pilot/inputs/pxd078459_files.tsv;$GIPPA_DATA_DIR/pxd078459/SHA256SUMS |
| wt2_zip | $GIPPA_DATA_DIR/pxd078459/MLA1914690_WT_P2_2_Y-D8_1_13120.d.zip | 05d6d6499ac7ec3cf1af10d7eeabc28d741e826d81f4768667386828554e3102 | validation/pilot/inputs/pxd078459_files.tsv;$GIPPA_DATA_DIR/pxd078459/SHA256SUMS |
| wt3_zip | $GIPPA_DATA_DIR/pxd078459/MLA1914691_WT_P3_3_Y-D9_1_13121.d.zip | c45618952c709de617934ca63eb29b9329286f41397c00341671056cdda1cd11 | validation/pilot/inputs/pxd078459_files.tsv;$GIPPA_DATA_DIR/pxd078459/SHA256SUMS |
| wt1_tdf | $GIPPA_DATA_DIR/pxd078459/WT_1.d/analysis.tdf | f4d5fd823d21d4d6c979bed74c74024d196e2393696bade787ae42b1bde66c5f | $GIPPA_ROOT/pilot/A1/provenance.json;$GIPPA_ROOT/pilot/A2/provenance.json |
| wt2_tdf | $GIPPA_DATA_DIR/pxd078459/WT_2.d/analysis.tdf | f77a283e669422688c6e7fe88d8d0d40c3921fa67ee0d02e4645d175abd71892 | $GIPPA_ROOT/pilot/A1/provenance.json;$GIPPA_ROOT/pilot/A2/provenance.json |
| wt3_tdf | $GIPPA_DATA_DIR/pxd078459/WT_3.d/analysis.tdf | 512187712c5a523b60155c0f746ac098826dfbd6c282a78567291eeb23babf9a | $GIPPA_ROOT/pilot/A1/provenance.json;$GIPPA_ROOT/pilot/A2/provenance.json |
| genome_wt | $GIPPA_ROOT/pilot/pass3/genome/nissle1917_wt.fna.gz | 1a654e6b83292353f15be4e7a7973b896d728761a4c9dfde533953f76e7961ad | validation/pilot/README.md |
| bakta_json_wt | $GIPPA_ROOT/pilot/pass3/recorded/bakta/1.12.2-full/nissle_wt/nissle1917_wt.json.gz | b680d455d8760a23476edf4d624dd5e816bba0df387098812354ea13adc7f223 | validation/pilot/README.md;$GIPPA_ROOT/pilot/pass3/recorded/bakta/1.12.2-full/nissle_wt/PROVENANCE.md |
| bakta_gbff_wt | $GIPPA_ROOT/pilot/pass3/recorded/bakta/1.12.2-full/nissle_wt/nissle1917_wt.gbff.gz | 1dfcdd287aca7bf6410ee0cede40ba77d7e8dc6e1b38e947f23994f4e15ac583 | $GIPPA_ROOT/pilot/pass3/recorded/bakta/1.12.2-full/nissle_wt/PROVENANCE.md |
| strain_db_wt | $GIPPA_ROOT/pilot/pass3/db_a2/strain_db.faa.gz | 20e36adfa049110906bc688d20ea83f8ca251f9afcc7add18d78ca31029417b8 | validation/pilot/README.md |

## Steps

### S1 resolve (M10-01)

Job: script=scripts/slurm/fetch_pxd078459.sbatch name=gippa-fetch-pxd078459 partition=normal qos=dev cpus=1 gpus=0 mem=2G time=00:30:00 constraint=- array=-

```bash
sbatch --qos=dev --cpus-per-task=1 --mem=2G --time=00:30:00 \
  "$GIPPA_REPO/scripts/slurm/fetch_pxd078459.sbatch" --resolve
```

Reads: -
Writes: $GIPPA_DATA_DIR/pxd078459/record/

The draft manifest was reviewed by hand and committed as
`validation/pilot/inputs/pxd078459_files.tsv` (README › Raw data). Job 1879851 was submitted with
the normal QOS (4 CPUs, 8G, 1 h). It ended with exit 2 after writing the record and the draft.

### S2 download (M10-01)

Job: script=scripts/slurm/fetch_pxd078459.sbatch name=gippa-fetch-pxd078459 partition=normal qos=normal cpus=2 gpus=0 mem=8G time=24:00:00 constraint=- array=-

```bash
sbatch "$GIPPA_REPO/scripts/slurm/fetch_pxd078459.sbatch"
```

Reads: px_record, wt1_zip, wt2_zip, wt3_zip
Writes: $GIPPA_DATA_DIR/pxd078459/

### S3 acquisition mode (M10-01)

Job: script=- name=gippa-acqmode partition=normal qos=dev cpus=2 gpus=0 mem=2G time=00:15:00 constraint=- array=-

```bash
srun -J gippa-acqmode -p normal --qos=dev -c 2 --mem=2G -t 00:15:00 \
  "$PY" -I "$GIPPA_REPO/validation/pilot/acquisition_mode.py" \
  "$GIPPA_DATA_DIR"/pxd078459/WT_1.d "$GIPPA_DATA_DIR"/pxd078459/WT_2.d "$GIPPA_DATA_DIR"/pxd078459/WT_3.d \
  --out "$TMPDIR/acquisition.json"
cmp "$TMPDIR/acquisition.json" "$GIPPA_REPO/validation/pilot/inputs/acquisition.json"
```

Reads: wt1_tdf, wt2_tdf, wt3_tdf
Writes: validation/pilot/inputs/acquisition.json (committed)

Job 1879882 ran it as a dev-QOS batch job. The `srun` above is the equivalent command. It writes
to `$TMPDIR` and compares, so a repeat never rewrites the committed file.

### S4 image install (M10-02)

Job: script=scripts/slurm/diann_image.sbatch name=gippa-diann-image partition=normal qos=normal cpus=1 gpus=0 mem=2G time=00:30:00 constraint=- array=-

```bash
sbatch "$GIPPA_REPO/scripts/slurm/diann_image.sbatch"
```

Reads: diann_image
Writes: $GIPPA_ROOT/images/diann_2.3.1.sif

Job 1879883 ended FAILED (exit 2) on the banner check fixed in e7edc457. The copy was verified by
hand afterwards (mode 444, the pin). S-check re-verifies it at F. The image is never re-pinned.

### S5 searches (M10-02)

Job: script=scripts/slurm/diann_pilot.sbatch name=gippa-diann-A1 partition=normal qos=normal cpus=16 gpus=0 mem=64000M time=24:00:00 constraint=- array=-

```bash
sbatch -J gippa-diann-A1 "$GIPPA_REPO/scripts/slurm/diann_pilot.sbatch" A1
# after the A1 job has ended (two jobs of one arm must never run at once):
sbatch -J gippa-diann-A2 "$GIPPA_REPO/scripts/slurm/diann_pilot.sbatch" A2
```

Reads: proteome_a1, contaminants, strain_db, cfg_a1, cfg_a2, acquisition, diann_image, wt1_tdf, wt2_tdf, wt3_tdf
Writes: $GIPPA_ROOT/pilot/A1/, $GIPPA_ROOT/pilot/A2/

**Not re-run at the freeze** (index decision 4). DIA-NN may not repeat byte for byte (M10-06 L3).
A later reader who repeats S5 should take the converted-table hashes below as the target. The
pass-1 jobs asked for `--mem=64G` and got 17 CPUs; the script now asks for 64000M (16 CPUs).

### S6 report check (M10-02)

Job: script=scripts/slurm/diann_report.sbatch name=gippa-diann-report-A1 partition=normal qos=normal cpus=4 gpus=0 mem=16000M time=01:00:00 constraint=- array=-

```bash
sbatch -J gippa-diann-report-A1 "$GIPPA_REPO/scripts/slurm/diann_report.sbatch" A1
sbatch -J gippa-diann-report-A2 "$GIPPA_REPO/scripts/slurm/diann_report.sbatch" A2
```

Reads: a1_report, a2_report
Writes: $GIPPA_ROOT/pilot/A1/check.json, $GIPPA_ROOT/pilot/A2/check.json, $GIPPA_ROOT/pilot/recorded/diann/2.3.1/

### S7 conversion (M10-03)

Job: script=scripts/slurm/diann_convert.sbatch name=gippa-diann-convert-A1 partition=normal qos=normal cpus=2 gpus=0 mem=8G time=00:30:00 constraint=- array=-

```bash
sbatch -J gippa-diann-convert-A1 "$GIPPA_REPO/scripts/slurm/diann_convert.sbatch" A1
sbatch -J gippa-diann-convert-A2 "$GIPPA_REPO/scripts/slurm/diann_convert.sbatch" A2
```

Reads: a1_report, a2_report
Writes: $GIPPA_ROOT/pilot/inputs/

S-check re-converts both reports at F and compares the hashes (example 7). A changed table is
re-converted with these commands at F, and the old Outputs row becomes `superseded` (example 8).

### S-check reproducibility check (M10-14)

Job: script=scripts/slurm/pilot_reproduce_check.sbatch name=gippa-pilot-repro-check partition=normal qos=normal cpus=2 gpus=0 mem=8G time=01:00:00 constraint=- array=-

```bash
sbatch "$GIPPA_REPO/scripts/slurm/pilot_reproduce_check.sbatch"
```

Reads: genome, bakta_json, strain_db, cfg_a1, cfg_a2, acquisition, proteome_a1, contaminants, diann_image, px_record, wt1_tdf, wt2_tdf, wt3_tdf, a1_report, a2_report
Writes: $GIPPA_ROOT/pilot/reproduce/

It runs at F, before S8, and must exit 0 (freeze gate item 4). It runs `static`, then `verify`
into `verify-<F7>.json`, then the re-conversion into `$TMPDIR`.

### S8 freeze A2 (M10-14)

Job: script=scripts/slurm/gippa_run.sbatch name=gippa-run-freeze-a2 partition=normal qos=normal cpus=2 gpus=0 mem=8G time=2-00:00:00 constraint=- array=- children=nextflow:normalx14cpu

```bash
G=$GIPPA_ROOT/pilot/pass3/genome/nissle1917_wt.fna.gz
J=$GIPPA_ROOT/pilot/pass3/recorded/bakta/1.12.2-full/nissle_wt/nissle1917_wt.json.gz   # next to it: nissle1917_wt.gbff.gz, so antiSMASH runs
P3=$GIPPA_ROOT/pilot/pass3/inputs
sbatch -J gippa-run-freeze-a2 $X "$GIPPA_REPO/scripts/slurm/gippa_run.sbatch" --genome $G --annotation $J \
  --peptides $P3/a2.tsv.gz --db-dir $GIPPA_DB_DIR --threads 14 --memory 54.GB --outdir $R/out_a2
```

Reads: genome_wt, bakta_json_wt, bakta_gbff_wt, a2_table_wt
Writes: $R/out_a2

Each arm ran after the previous head job ended. Each attempt is in `$R/ATTEMPTS.tsv`.

### S9 freeze A1 (M10-14)

Job: script=scripts/slurm/gippa_run.sbatch name=gippa-run-freeze-a1 partition=normal qos=normal cpus=2 gpus=0 mem=8G time=2-00:00:00 constraint=- array=- children=nextflow:normalx14cpu

```bash
# after S8's head job has ended (any state):
sbatch -J gippa-run-freeze-a1 $X "$GIPPA_REPO/scripts/slurm/gippa_run.sbatch" --genome $G --annotation $J \
  --peptides $P/a1.tsv.gz --db-dir $GIPPA_DB_DIR --threads 14 --memory 54.GB --outdir $R/out_a1
```

Reads: genome_wt, bakta_json_wt, bakta_gbff_wt, a1_table
Writes: $R/out_a1

### S10 freeze TE7 replicate (M10-14)

Job: script=scripts/slurm/gippa_run.sbatch name=gippa-run-freeze-a2-rep partition=normal qos=normal cpus=2 gpus=0 mem=8G time=2-00:00:00 constraint=- array=- children=nextflow:normalx14cpu

```bash
sbatch -J gippa-run-freeze-a2-rep $X "$GIPPA_REPO/scripts/slurm/gippa_run.sbatch" --genome $G --annotation $J \
  --peptides $P3/a2.tsv.gz --db-dir $GIPPA_DB_DIR --threads 14 --memory 54.GB --outdir $R/out_a2_rep
```

Reads: genome_wt, bakta_json_wt, bakta_gbff_wt, a2_table_wt
Writes: $R/out_a2_rep

### S11 freeze A3 Cascadia (M10-14)

Job: script=scripts/slurm/gippa_run.sbatch name=gippa-run-freeze-a3 partition=normal qos=normal cpus=2 gpus=0 mem=8G time=2-00:00:00 constraint=- array=- children=nextflow:normalx14cpu

```bash
sbatch -J gippa-run-freeze-a3 $X "$GIPPA_REPO/scripts/slurm/gippa_run.sbatch" --genome $G --annotation $J \
  --peptides $P/a3.tsv.gz --reference-peptides $P3/a2.tsv.gz --db-dir $GIPPA_DB_DIR --threads 14 --memory 54.GB --outdir $R/out_a3
```

Reads: genome_wt, bakta_json_wt, bakta_gbff_wt, a3_table, a2_table_wt
Writes: $R/out_a3

### S12 compare (M10-04, M10-10)

Job: script=scripts/slurm/pilot_compare.sbatch name=gippa-pilot-compare partition=normal qos=normal cpus=2 gpus=0 mem=8000M time=01:00:00 constraint=- array=-

```bash
sbatch -J gippa-pilot-compare --cpus-per-task=2 --mem=8000M $X "$GIPPA_REPO/scripts/slurm/pilot_compare.sbatch" $R/out_a1 $R/out_a2 $R/cmp_a1_a2
sbatch -J gippa-pilot-compare --cpus-per-task=2 --mem=8000M $X "$GIPPA_REPO/scripts/slurm/pilot_compare.sbatch" $R/out_a1 $R/out_a2 $R/out_a3 $R/cmp_a1_a2_a3
```

Reads: genome_wt
Writes: $R/cmp_a1_a2, $R/cmp_a1_a2_a3, $R/evidence/

The second compare ran under the job name `gippa-pilot-compare3`.

### S13 evidence (M10-04)

Job: script=- name=gippa-pilot-evidence partition=normal qos=dev cpus=1 gpus=0 mem=4000M time=00:30:00 constraint=- array=-

```bash
cd "$GIPPA_REPO"; export PYTHONPATH="$GIPPA_REPO/src:$GIPPA_REPO"
"$PY" validation/pilot/run_evidence.py arm $R/out_a2 --arm A2 --job <S8 job> --expect-commit <F> --expect-input genome=<sha> \
  --expect-input annotation=<sha> --expect-input peptides=<sha> --out $R/evidence/a2.json
# the same for out_a1, out_a2_rep and out_a3 (--expect-input reference_peptides=<sha> --a2 $R/out_a2)
```

Reads: genome_wt, bakta_json_wt, a1_table, a2_table_wt, a3_table
Writes: $R/evidence/

### S14 TE and metrics (M10-05, M10-06, M10-07)

Job: script=- name=gippa-pilot-te partition=normal qos=dev cpus=4 gpus=0 mem=16000M time=01:00:00 constraint=- array=-

```bash
M=$GIPPA_ROOT/pilot/metrics/freeze-$F7
"$PY" validation/pilot/te_real.py compute --a2 $R/out_a2 --a1 $R/out_a1 --peptides-a1 $P/a1.tsv.gz --peptides-a2 $P3/a2.tsv.gz \
  --contaminants $GIPPA_DB_DIR/contaminants/contaminants.fasta --out $M
"$PY" validation/pilot/te_real.py te7 --ref $R/out_a2 --rep $R/out_a2_rep --repo $GIPPA_REPO --out $M
# Entrapment: the first attempt at F failed its guard (tool bug); the rerun that counts ran with GIPPA_REPO at the main
# checkout, commit db3a95f1 (fix #522; src/ and workflow/ identical to F), not at F. A deviation; see REPORT.md section 10.
sbatch $X,PEP_A1=$P/a1.tsv.gz,GENOME=$G,STRAIN_DB=$GIPPA_ROOT/pilot/pass3/db_a2/strain_db.faa.gz,ANNOTATION=$J,OUT_A1=$R/out_a1,METRICS_DIR=$M "$GIPPA_REPO/scripts/slurm/entrapment.sbatch"
"$PY" validation/pilot/metrics.py compute --a1 $R/out_a1 --a2 $R/out_a2 --cmp $R/cmp_a1_a2 --cmp3 $R/cmp_a1_a2_a3 --a3 $R/out_a3 \
  --peptides-a1 $P/a1.tsv.gz --peptides-a2 $P3/a2.tsv.gz --peptides-a3 $P/a3.tsv.gz --convert-a1 $P/a1.convert.json \
  --convert-a2 $P3/a2.convert.json --convert-a3 $P/a3.convert.json --bakta-json $P3/nissle1917_wt.bakta.json --genome $G \
  --k12-fasta $GIPPA_DATA_DIR/proteomes/UP000000625.fasta --strain-db $GIPPA_ROOT/pilot/pass3/db_a2/strain_db.faa.gz \
  --te-real $M/te_real.json --te7 $M/te7.json --entrapment $M/entrapment.json --a3-feasibility $GIPPA_ROOT/pilot/pass3/a3_feasibility.json \
  --baseline benchmarks/baseline.json --benchmark-adr docs/adr/ADR-0072-benchmark-baseline-1-0-0.md --sacct $R/sacct.tsv --out $M
```

Reads: genome_wt, a1_table, a2_table_wt, a3_table, strain_db_wt, contaminants, proteome_a1
Writes: $GIPPA_ROOT/pilot/metrics/freeze-fa49b8f/

### S-bench standard-run benchmark (M10-12)

Job: script=scripts/slurm/benchmark_standard.sbatch name=gippa-bench partition=normal qos=normal cpus=16 gpus=0 mem=64000M time=22:00:00 constraint=CPU_SKU:6148 array=0-2%1

```bash
sbatch $X,BENCH_COMMIT=$F,BENCH_ID=2026-10-31_$F7 "$GIPPA_REPO/scripts/slurm/benchmark_standard.sbatch"
```

Reads: -
Writes: $GIPPA_ROOT/bench/2026-10-31_fa49b8f/

All four cells of M10-12's matrix ran, in one array at F. Array tasks 0 and 1 (repeats 1 and 2, all 8 runs ok) are the evidence; task 2
(repeat 3) was cancelled before it started, on Khoa's decision (ADR-0072). The baseline was written by main's
`scripts/dev/benchmark_summary.py --rounds 2 --adr ADR-0072` (#523): that option is not in the script at F, and the numbers it reads are from F.

### S15 clean-up (M10-14)

Job: script=scripts/slurm/pilot_cleanup.sbatch name=gippa-pilot-cleanup partition=normal qos=normal cpus=2 gpus=0 mem=8G time=02:00:00 constraint=- array=-

```bash
sbatch "$GIPPA_REPO/scripts/slurm/pilot_cleanup.sbatch"            # dry run: plan.tsv only
CLEAN=1 sbatch "$GIPPA_REPO/scripts/slurm/pilot_cleanup.sbatch"    # after reading the plan
```

Reads: -
Writes: $GIPPA_ROOT/pilot/cleanup/

Not run for this report: the clean-up is the lead's step, and K2 is still open. It runs once, last. It deletes the Nextflow work directories and the mzML files. Until Khoa
answers K2, the DIA-NN reports and libraries are listed as "kept pending K2" and stay. After the
clean-up, M10-03's `slow` full-report test still finds the reports, and the 1,000-row sample test
keeps running.

## A3

### A3-cascadia (status: active)

Active when `validation/pilot/metrics/a3_feasibility.json` has verdict `pass`. The Cascadia
searches and their conversion are not re-run. S11 runs `gippa run` on `$P/a3.tsv.gz` with
`--reference-peptides $P3/a2.tsv.gz` into `$R/out_a3`, and S12 compares A1, A2 and A3. The
commands come from M10-09 and M10-10.

### A3-fallback (status: inactive)

Active when the verdict is `fail`, or when an E1-E3 outcome sent A3 to M10-11. Casanovo is not
re-run. S11 runs `gippa run` on the UTI89 genome without `--reference-peptides` into
`$R/out_a3_uti89`, and S12 compares only A1 and A2. The inputs join the Inputs table from
M10-11a's manifests.

## Verify a reproduction

Run the S8 command at F into a new outdir. Then run M10-06's
`te_real.py te7 --ref $R/out_a2 --rep <new> --repo $GIPPA_REPO --out <dir>`. "Reproduced" means
TE7 `pass`: every gate file is byte-identical, so the gate digest matches the Outputs row. A
`known_lapack` diagnosis on another CPU class is the documented non-portability (M4-20 decision
9). It does not count as reproduced.

## Evidence

Pass-1 jobs (`sacct`). The `commit` column is the commit the job ran at, where it was recorded.

| step | label | job | state | exit | elapsed | maxrss | node | commit |
|---|---|---|---|---|---|---|---|---|
| S1 | resolve (first try) | 1879840 | superseded | 2:0 | 00:00:07 | 14M | wheat-02 | - |
| S1 | resolve | 1879851 | FAILED | 2:0 | 00:00:07 | 15M | wheat-08 | - |
| S2 | download | 1879859 | COMPLETED | 0:0 | 00:23:21 | 5.7G | wheat-07 | - |
| S3 | acquisition mode | 1879882 | COMPLETED | 0:0 | 00:00:07 | 12M | wheat-05 | - |
| S4 | image install | 1879883 | FAILED | 2:0 | 00:00:16 | 654M | wheat-05 | - |
| S5 | A1 step 1 | 1879913 | COMPLETED | 0:0 | 00:08:52 | 8.0G | barley-03 | 26adc18c40e6eceb7f30b42f84e27a491d606b82 |
| S5 | A2 step 1 | 1879914 | superseded | 1:0 | 00:08:06 | 8.1G | barley-04 | 26adc18c40e6eceb7f30b42f84e27a491d606b82 |
| S5 | A1 | 1879955 | COMPLETED | 0:0 | 00:14:21 | 6.6G | rye-02 | 81272fc0c68dd45fc4f82309d9eb3e6c770ee0df |
| S5 | A2 | 1879956 | COMPLETED | 0:0 | 00:14:21 | 6.4G | rye-02 | 81272fc0c68dd45fc4f82309d9eb3e6c770ee0df |
| S6 | A1 | 1880002 | COMPLETED | 0:0 | 00:00:11 | 203M | wheat-01 | - |
| S6 | A2 | 1880003 | COMPLETED | 0:0 | 00:00:07 | 203M | wheat-08 | - |
| S7 | A1 (before acquisition.json) | 1880208 | superseded | 0:0 | 00:00:18 | - | wheat-01 | b051c2b75f1daf341dcf13dc491b4cbbee402903 |
| S7 | A2 (before acquisition.json) | 1880209 | superseded | 0:0 | 00:00:20 | - | wheat-05 | b051c2b75f1daf341dcf13dc491b4cbbee402903 |
| S7 | A1 | 1880219 | COMPLETED | 0:0 | 00:00:17 | 266M | wheat-01 | fdf586be0b1310a484ac51c8483f2e3b9b6d663c |
| S7 | A2 | 1880220 | COMPLETED | 0:0 | 00:00:17 | 290M | wheat-05 | fdf586be0b1310a484ac51c8483f2e3b9b6d663c |
| S-check | at F (draft runbook) | 1914946 | superseded | 0:0 | 00:00:23 | 867208K | wheat-01 | fa49b8f145ea4598a5f896f9d0b1956912fc5df3 |
| S-check | at F | 1922538 | COMPLETED | 0:0 | 00:00:11 | 219332K | barley-03 | fa49b8f145ea4598a5f896f9d0b1956912fc5df3 |
| S8 | A2 freeze | 1914948 | COMPLETED | 0:0 | 02:50:55 | 1497620K | wheat-07 | fa49b8f145ea4598a5f896f9d0b1956912fc5df3 |
| S9 | A1 freeze | 1915264 | COMPLETED | 0:0 | 02:36:36 | 1446348K | wheat-03 | fa49b8f145ea4598a5f896f9d0b1956912fc5df3 |
| S10 | A2 replicate | 1915683 | COMPLETED | 0:0 | 02:49:15 | 1142356K | wheat-03 | fa49b8f145ea4598a5f896f9d0b1956912fc5df3 |
| S11 | A3 freeze | 1916086 | COMPLETED | 0:0 | 01:35:10 | 1177696K | wheat-03 | fa49b8f145ea4598a5f896f9d0b1956912fc5df3 |
| S12 | compare A1, A2 | 1916505 | COMPLETED | 0:0 | 00:00:21 | 278716K | wheat-05 | fa49b8f145ea4598a5f896f9d0b1956912fc5df3 |
| S12 | compare A1, A2, A3 | 1916515 | COMPLETED | 0:0 | 00:00:20 | 316712K | wheat-06 | fa49b8f145ea4598a5f896f9d0b1956912fc5df3 |
| S13 | rule set A | 1916525 | COMPLETED | 0:0 | 00:00:10 | 145148K | wheat-07 | fa49b8f145ea4598a5f896f9d0b1956912fc5df3 |
| S14 | TE2-TE5, TE7 | 1916553 | COMPLETED | 0:0 | 00:00:37 | 393068K | wheat-01 | fa49b8f145ea4598a5f896f9d0b1956912fc5df3 |
| S14 | entrapment, first attempt (tool bug, fixed in #522) | 1916565 | superseded | 1:0 | 00:01:04 | 742340K | wheat-05 | fa49b8f145ea4598a5f896f9d0b1956912fc5df3 |
| S14 | entrapment rerun at main, not F (#522) | 1917692 | COMPLETED | 0:0 | 00:00:37 | 744392K | barley-04 | db3a95f16342efbdc3a7cdd6ed06b199818df595 |
| S14 | metrics | 1921926 | COMPLETED | 0:0 | 00:00:27 | 391484K | wheat-01 | fa49b8f145ea4598a5f896f9d0b1956912fc5df3 |
| S-bench | r1 | 1916711_0 | COMPLETED | 0:0 | 13:56:48 | 65405340K | wheat-01 | fa49b8f145ea4598a5f896f9d0b1956912fc5df3 |
| S-bench | r2 | 1916711_1 | COMPLETED | 0:0 | 13:57:34 | 65133744K | wheat-01 | fa49b8f145ea4598a5f896f9d0b1956912fc5df3 |

## Outputs

The DIA-NN reports (kept pending K2) and the converted tables. The freeze `report.html` files,
`compare.html`, `metrics.json` and the per-arm gate digests join in 14b.

| id | path | sha256 | recorded_in | kept |
|---|---|---|---|---|
| a1_report | $GIPPA_ROOT/pilot/A1/A1_report.parquet | f217a40fc0777e99783f3d410743bc80d301290dcdb219a6c69574ee6a4926d7 | validation/pilot/README.md;$GIPPA_ROOT/pilot/A1/outputs.sha256 | yes |
| a2_report | $GIPPA_ROOT/pilot/A2/A2_report.parquet | 733a3884a15664f13385b241e48186e538308993ff2d5a8cc3ec9ff212866d36 | validation/pilot/README.md;$GIPPA_ROOT/pilot/A2/outputs.sha256 | yes |
| a1_table | $GIPPA_ROOT/pilot/inputs/a1.tsv.gz | 81818184f0cd222f091ae6b835ff72cc5743b6ce7e6403e7a2cb1740bbf24f3a | validation/pilot/README.md;$GIPPA_ROOT/pilot/inputs/a1.convert.json | yes |
| a1_abundance | $GIPPA_ROOT/pilot/inputs/a1_abundance.tsv.gz | 20284d4234055b26cf779dc405548087e4efc9bf1007dc3bebc031c5b6c74c54 | validation/pilot/README.md;$GIPPA_ROOT/pilot/inputs/a1.convert.json | yes |
| a2_table | $GIPPA_ROOT/pilot/inputs/a2.tsv.gz | 1eb5e1e9a8855fc8ac62bed4d60f2715f9270f2b96f40206dd35a3a754adf14f | validation/pilot/README.md;$GIPPA_ROOT/pilot/inputs/a2.convert.json | yes |
| a2_abundance | $GIPPA_ROOT/pilot/inputs/a2_abundance.tsv.gz | c61122f3afc9aa31986de8bbfef1c542719627d3d632ec9113e7719582b41c19 | validation/pilot/README.md;$GIPPA_ROOT/pilot/inputs/a2.convert.json | yes |
| a2_table_wt | $GIPPA_ROOT/pilot/pass3/inputs/a2.tsv.gz | d7f1ddf5d2badb0038b50f55a6d99b43be61ccc9e2c8fadff9f5a7ab0c405fbf | validation/pilot/README.md;$GIPPA_ROOT/pilot/pass3/inputs/a2.convert.json | yes |
| a2_abundance_wt | $GIPPA_ROOT/pilot/pass3/inputs/a2_abundance.tsv.gz | 3970344b9c0edc0fcf756b267063a1d1c30c05f6b6bdb6654d973eae1e6f888e | validation/pilot/README.md;$GIPPA_ROOT/pilot/pass3/inputs/a2.convert.json | yes |
| a2_report_wt | $GIPPA_ROOT/pilot/pass3/A2/A2_report.parquet | c274bbb028b386c67d72410fd5ae48a6eaf0edeb7ac2a8030ab8bbf73ea6bc06 | validation/pilot/README.md;$GIPPA_ROOT/pilot/pass3/A2/outputs.sha256 | yes |
| a3_table | $GIPPA_ROOT/pilot/inputs/a3.tsv.gz | 743396d090db4f7b0e1075d63d68c4a81b903b8b6fd3a8b61d80fedcac1bf0b8 | validation/pilot/README.md;$GIPPA_ROOT/pilot/inputs/a3.convert.json | yes |
| freeze_report_a2 | $GIPPA_ROOT/runs/pilot/freeze/fa49b8f/out_a2/report/report.html | 5def12406f78772f11425af3c77e4696058daa6133a1e799001eef1252ae900f | - | yes |
| gate_digest_a2 | $GIPPA_ROOT/runs/pilot/freeze/fa49b8f/out_a2/run_manifest.json | edf95a8b5843d41b097cf4254bad9609a86a69a8754fddd95352e1335bdb953b | - | yes |
| freeze_report_a1 | $GIPPA_ROOT/runs/pilot/freeze/fa49b8f/out_a1/report/report.html | 0f5dd331489bcf2d9dfaf38bc78a56977aebd3028046fb0977f1d4640914ff0b | - | yes |
| gate_digest_a1 | $GIPPA_ROOT/runs/pilot/freeze/fa49b8f/out_a1/run_manifest.json | d4862323a820393f629361df1108b2c66aaa5862d610e688e074f86d6ebd8f0f | - | yes |
| freeze_report_a2_rep | $GIPPA_ROOT/runs/pilot/freeze/fa49b8f/out_a2_rep/report/report.html | b8c7671d6cf294dd0cd18730f2c1e98ba2a985df54b8d9471cde1cc2c7164b62 | - | yes |
| gate_digest_a2_rep | $GIPPA_ROOT/runs/pilot/freeze/fa49b8f/out_a2_rep/run_manifest.json | edf95a8b5843d41b097cf4254bad9609a86a69a8754fddd95352e1335bdb953b | - | yes |
| freeze_report_a3 | $GIPPA_ROOT/runs/pilot/freeze/fa49b8f/out_a3/report/report.html | 21d838d35e232cd5f68f1aeac2677b49669df40d8ec31ce84470e16847b1a5de | - | yes |
| gate_digest_a3 | $GIPPA_ROOT/runs/pilot/freeze/fa49b8f/out_a3/run_manifest.json | f2c2ec768527d47e8e2f51efdfd9e114e7572abdbda138f1815a3470ec40989d | - | yes |
| freeze_compare2_html | $GIPPA_ROOT/runs/pilot/freeze/fa49b8f/cmp_a1_a2/report/compare.html | 5f3c63ba080a552f92591135d8ca0f4440c9512cf0ccdbcd5aac76988d8aad99 | - | yes |
| freeze_compare3_html | $GIPPA_ROOT/runs/pilot/freeze/fa49b8f/cmp_a1_a2_a3/report/compare.html | d9168a0efae3ef93a8c45492ed360a4c2b46afc98800ad73710e6c4ac8fa896c | - | yes |
| freeze_metrics | $GIPPA_ROOT/pilot/metrics/freeze-fa49b8f/metrics.json | a92f7eadda3da0c6876d4334c3d51f04ffc291c90539c73307e6fa6633918d7d | - | yes |
