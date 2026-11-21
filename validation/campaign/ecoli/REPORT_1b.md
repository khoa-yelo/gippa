# Phase 1b report: NCM3722 strain-variant arm (N1, N2)

Checks pre-registered in `ANSWER_KEY_1b.md` (commit 18a47b54, branch `validation/campaign-ncm3722`)
before any download or search. Aggregates only: no peptide sequences or DIA-NN rows. Every share
has a Wilson 95 % interval (optimistic, because peptides and sites are treated as independent).
Runs on 2026-11-19/20 at `a2a95c76` (v1.1.0); N1 and N2 rerun on 2026-11-21 at `7d4c4bfd` (1.1.1 code, ADR-0078 fix). Scripts: `validation/campaign/ncm3722/`;
`$D = $GIPPA_ROOT/campaign/ncm3722`.

## Verdicts

| # | Check | Result | Verdict |
|---|---|---|---|
| 1b-1 | Recall | 0/11 covered truth sites = 0.00 [0.00, 0.26] (any confidence and `high`); eligible stratum 0/10 [0.00, 0.28]; prior ≥ 0.8 not met | reported: no recall |
| 1b-2 | Precision of `high` ALT | 0 `high` ALT rows, 0 `mixed`; `variants.tsv` empty | unestimated (n = 0) |
| 1b-3 | `high` ALT calls not in the truth set | 0 (mass_equivalent 0, start/stop 0, unexplained 0) | none |
| 1b-4 | N1 peptides unplaced on MG1655, by N1 replicon | chromosome 63/24,594 = 0.26 % [0.20, 0.33]; plasmid F 242/242 = 100 % [98.4, 100]; both replicons 0/6 [0, 39]; all 305/24,842 = 1.23 % [1.10, 1.37]; priors met | reported |
| 1b-5 | Placement FDR ≤ 1 % (TE8 analogue) | pooled TE8 estimator 0.35 % (24,538 accepted, all exact); run `mapping_fdr.tsv` stratum 0: 5.6e-5 (24,539 targets / 86 decoys); strata 1 and 2 accept nothing | pass (error unestimated); exact placements only |
| 1b-6 | N1 vs N2 | `gippa compare` exit 1 ("different genomes"), as pre-registered. Ortholog join (4,235 pairs): Spearman 0.9988 [0.9963, 1.000], n = 2,280 (1,000 gene bootstraps, seed 20261014); N1-detected orthologs kept 2,284/2,284 [0.998, 1]; ~~detection Jaccard 0.961~~ (1.1.0, superseded: 1.000 with the fix, see Rerun); islands 0 vs 0 (Jaccard undefined); ~~operon flags 10 vs 9, 8 shared, Jaccard 0.73~~ (superseded: 10 vs 8, 8 shared, 0.80); prior ρ ≥ 0.9 met | reported |

## Why recall is zero (exploratory)

- 11 of 31 truth sites are covered by 12 distinct N1 evidence peptides.
- 11 of the 12 reach placement on MG1655 as 1- or 2-difference tolerant hits; the 12th is 7 aa
  (`too_short`). None is exact.
- All 11 are rejected (`rejected_placement_fdr`). The 1-difference stratum has 13 targets (11 of
  them truth peptides) against 7 reverse decoys; pooled over length classes its FDR is 0.62. The
  2-difference stratum has 1 target against 4 decoys.
- Even with zero decoys, (0 + 1)/13 = 0.077 > α = 0.01: a stratum with fewer than about 100
  tolerant targets cannot accept anything. This is the same floor as
  `novel_orfs.class_q_floor_above_alpha`, but the run gives no warning for it (its only warning is
  `novel_orfs.decoys_not_makedb`).
- So 1b-5 passes on exact hits only, and GIPPA 1.1.0 at defaults cannot report substitutions
  between strains this close.

## Truth set

- Method: blastp 2.16.0+ (pinned image `sha256:9dfc69f9…959c`) on the Bakta 1.12.2 full-DB CDS and
  sORF proteins of both genomes; reciprocal best hits with ≥ 90 % identity and ≥ 80 % HSP coverage of
  both proteins (`blastp -seg no -evalue 1e-10 -max_hsps 1`). Job 1998777, run before any N2 output
  existed.
- NCM3722: 4,498 proteins (4,423 chromosome, 75 plasmid F); MG1655: 4,373.
- 4,235 pairs kept, all chromosomal; 4,193 identical.
- Unpaired NCM3722 proteins: chromosome 93 hit but not kept, 95 no hit; plasmid 6 hit but not kept,
  69 no hit.
- 31 substitution sites in 26 proteins. Not scored: 1 I↔L column, 1 gap column, 528 residues
  outside HSPs.
- **Deviation from the pre-registration:** the codon check removed 2 raw sites (33 → 31) before any
  N2 output existed. In both, the NCM3722 residue was Bakta's initiator M over a non-ATG start whose
  bases equal the MG1655 codon. Counts: `validation/campaign/ncm3722/truth_counts.json`.

## N1 baseline (not gated)

- Preflight `--strict` exited 0 on both arms (91,603 rows, decoy share 0.99 %).
- Both runs completed with exit 0, `git_dirty` false, all modules; only warning
  `novel_orfs.decoys_not_makedb`.
- `decoy_picked` tier: 2,341 genes at q ≤ 0.01, estimated FDR 0.0099 (2,315/4,423 chromosomal,
  26/75 plasmid F; 2,138 with ≥ 2 peptides).
- 25,222 exact placements, 99.6 % `cds_inframe`.
- ~~N2 detects 2,399 genes~~ (1.1.0, superseded), including about 94 orthologs N1 does not detect. **These are a GIPPA defect, not
  biology** (see Correction below). With the ADR-0078 fix N2 detects 2,305 genes (see Rerun).

## Command log

All `gippa` commands ran from `local/worktrees/campaign-pin` (`a2a95c76`) with
`PYTHONPATH=<pin>/src` and the venv entry points. Every sbatch used
`--exclude=barley-01..04,rye-01,rye-02` and an explicit `--time`.

- **1998570 genome + Bakta** (21 min 18 s, 16 CPUs): NCBI GCF_001043215.1 `genomic.fna.gz`, MD5
  `b3a7430b1beda7a8541a6953667cc28e` = `md5checksums.txt`; both replicons tagged
  `[topology=circular]` (NCBI esummary: circular, complete); `bakta_command` argv with
  `--keep-contig-headers`, image `sha256:1e24ae2f…7f18`, DB 6.0 full. 4,419 CDS, 79 sORFs,
  5 pseudogenes.
- **1998571 fetch** (1 h 01 min, 1 CPU): PRIDE HTTPS, resumable; sizes match the PRIDE list (no
  deposit checksum). sha256: FRH_049 `13a859dc…6f25`, FRH_056 `0b85bef9…7d94`, FRH_061
  `62573261…2668`.
- **1998776 makedb + DIA-NN** (13 min 49 s: predict 121 s, search 689 s; 16 CPUs):
  `gippa makedb --annotation $D/bakta/out/ncm3722.json --outdir $D/db_a` (4,874 records, 381
  `CONT_`); DIA-NN 2.3.1 with the `[step1]`/`[step2]` flags of `diann_A2.cfg`. About 28k precursors (28,025–28,761)
  and 1,040–1,067 genes per run at 1 %; 2,301 protein groups.
- **1998777 truth set** (1 min 08 s).
- **1998866 N1** (2 h 13 min 31 s, 16 CPUs): `gippa convert --from diann $D/diann_a/report.parquet
  --search-db same`; `gippa preflight --strict` (exit 0); `gippa run --profile
  farmshare_local,apptainer --genome $D/genome/ncm3722.fna.gz --annotation
  $D/bakta/out/ncm3722.json --peptides $D/in_n1/peptides.tsv.gz --db-dir $GIPPA_DB_DIR --threads 14
  --memory 54.GB --outdir $D/out_n1`.
- **1998867 N2** (2 h 18 min 46 s): same report converted with `--search-db other`; genome
  `tests/data/genomes/k12_mg1655.fna.gz` (U00096.3); recorded Bakta JSON (job 1849096);
  `--genome-relation related`.
- **1998868 placement FDR** (1 min 18 s): TE8 chain on 24,843 N2 evidence peptides, then
  `pooled_decoy_fdr`.
- **1999152 compare + score** (9 s): `gippa compare … --labels ncm3722,mg1655` (exit 1), then
  `score.py`.
- **Failed first attempts 1998861 and 1998862:** `gippa: command not found` inside Nextflow tasks
  (venv `bin/` not on PATH); fixed and resubmitted as 1998866/1998867.
- **Deviation:** each run used the `farmshare_local` profile inside one 16-CPU job instead of the
  `farmshare` Slurm profile, to stay within this arm's 48-CPU share of the R14 cap.
- **Provenance:** `$D/n1/provenance.json`, `$D/n2/provenance.json`. DIA-NN image
  `sha256:1946452d…a6381`.

## Decisions recorded (Khoa, 2026-11-20)

- No cross-genome `gippa compare`: the ortholog join in `score.py` stays campaign tooling.
- The ~94 extra N2 orthologs were investigated (below).
- Open: how to handle the tolerant-placement FDR floor (warning, issue/ADR, or documentation).

## Correction: N2's gene-level FDR is anti-conservative (GIPPA defect, 2026-11-20)

Investigated at Khoa's request and confirmed by an independent statistician review (aggregates in
the lead session's scratchpad `n2extra/`, `n2review/`):
- The 94 extra genes have identical peptide evidence in N1 and N2 (same unique/shared peptides, same
  best peptide q); only their gene q differs (N1 0.010–0.026, N2 0.0026–0.0051).
- Cause: in `genetable/detection.py`, with no protein hits and protein IDs that are not locus tags of
  the run's annotation, a protein is mapped to a gene only through votes of **accepted** target
  peptides, so `pair_decoys` leaves unpaired exactly the decoys whose protein has no accepted target
  — the would-be decoy wins. Unpaired decoys enter only the paired share f, so the 1/f correction of
  ADR-0009 cannot repair the loss.
- N1 vs N2: decoys paired/unpaired 301/3 vs 244/60; decoy wins 62 vs 9 (53 vs 0 on genes without an
  accepted target). Pairing N1 by vote only reproduces N2 (2,439 genes, max q 0.0050).
- Estimated realised FDR of the N2 list: about 0.026–0.028 (target–decoy estimates) at a nominal
  0.01 (reported max q 0.005).
- Consequence for this report: the N2 gene counts and the "kept 2,284/2,284" and detection-Jaccard
  figures of 1b-6 are inflated by this defect; the abundance correlation on shared genes is not
  affected. N1 and the K-12 arms E1/E2 (locus-tag protein IDs) are unaffected.
- A GitHub issue is drafted (fix options: pair decoys independently of target acceptance, or count
  unpaired decoys as wins; plus a new attack null). Phase 2 cross-strain runs wait for the fix.

## Rerun with the fix (1.1.1 code, 7d4c4bfd)

These figures supersede the 1.1.0 N2 gene counts and the 1b-6 detection figures above. The 1.1.0
numbers above are kept for the record.

- Code: the pinned detached worktree `local/worktrees/fix-pin` at `7d4c4bfd` (PR #607, ADR-0078
  option A), `PYTHONPATH=<fix-pin>/src`, venv `bin/` on PATH inside Nextflow tasks. The manifests
  record `git_commit` 7d4c4bfd and `git_dirty` false. The version string still reads 1.1.0 because
  the fix is not released yet.
- Command: the Command log's `gippa_arm.sbatch` with `PIN=<fix-pin> SUF=_fix`. The new `SUF` option
  writes `in_<arm>_fix`, `preflight_<arm>_fix` and `out_<arm>_fix`, and refuses an existing outdir,
  so the 1.1.0 outputs are untouched. Convert, preflight and run arguments are unchanged. The
  input sha256s (genome, Bakta JSON, DIA-NN report, converted peptides and abundance) are identical
  to the 1.1.0 runs.
- **2001968 N2** (`gippa-camp-n2fix`): 3 h 26 min 37 s, 16 CPUs, wheat-01, MaxRSS 39.3 GB.
  **2001969 N1** (`gippa-camp-n1fix`): 2 h 30 min 10 s, 16 CPUs, wheat-03, MaxRSS 24.0 GB.
  - Both ran with `--time=04:00:00` and `--exclude=barley-01..04,rye-01,rye-02`.
  - Deviation, as before: the `farmshare_local` profile ran every task inside one 16-CPU job.
    Total GIPPA CPUs during the reruns: 33.
  - Preflight `--strict` exited 0 and the run exited 0 on both arms, with all modules. The only
    warning was `novel_orfs.decoys_not_makedb`.
- Score: `score.py` is unchanged (seed 20261014). It ran on `$D/score_fix/<pair>/`, a directory
  that symlinks `truth`, `out_n1` and `out_n2`, so `$D/score/` is untouched.
  - It ran against the 1.1.0 N1 (as pre-registered) and against the rerun N1. Both give identical
    figures.
  - `gippa compare` still exits 1 ("different genomes").
- Decoy counts were recomputed from the run outputs with the fix-pin `genetable.detection`
  functions. The recomputed detections and max q equal `genes.tsv`.

| | N1 1.1.0 | N1 fix | N2 1.1.0 (superseded) | N2 fix |
|---|---|---|---|---|
| Genes at q ≤ 0.01 (`decoy_picked`) | 2,341 | 2,341 | ~~2,399~~ | **2,305** |
| Max gene q of the list | 0.0099 | 0.0098 | ~~0.0051~~ | **0.0100** |
| `n_detected_2pep` summary row (own criterion, not a subset of the q list) | 2,138 | 2,138 | 2,113 | 2,113 |
| Decoys q ≤ 0.01: gene-paired / own-unit rows | 301 / 3 | 301 / 3 (2 own units) | 244 / 60 (unpaired, dropped from D) | 244 / 60 (57 own units) |
| `decoy_paired_share` | 0.990 | 0.990 | 0.803 | 0.803 |
| `decoy_win_weight` | 1.01 | 1 (retired) | 1.246 | 1 (retired) |
| Decoy wins in genes (on genes without accepted target) | 62 (53) | 62 (53) | 9 (0) | 9 (0) |
| Decoy wins incl. own units | n/a | 64 | n/a | 66 |
| Target wins (genes with finite q) | 2,439 | 2,439 | 2,399 | 2,399 |

- N2 now has 2,305 genes at max q 0.0100, down from 2,399 at max q 0.0051. That is 94 fewer
  genes, the size of the extra set in the Correction. This matches the statistician's
  function-level check (2,305 genes, max q 0.0100).
- N1 keeps the same 2,341 genes. Its max q moves from 0.00992 to 0.00982 because 1/f (f = 0.990)
  is retired.
- The 57 own units of N2 count as decoy wins in full. N2's total decoy wins (66) are now in line
  with N1's (64). In 1.1.0, N2 lost 53 of its wins.
- 1b-6, recomputed (ortholog join, 4,235 pairs):
  - N1-detected orthologs kept in N2: 2,284/2,284 = 1.00 [0.998, 1] (unchanged).
  - Detection Jaccard: **1.000** (was 0.961). The N2 list now detects exactly the orthologs N1
    detects.
  - Spearman: 0.9988 [0.9963, 1.000], n = 2,280 (unchanged, as expected).
  - Operon flags: 10 vs 8, 8 shared, Jaccard 0.80 (was 10 vs 9, 0.73).
  - Islands: 0 vs 0.
  - The prior ρ ≥ 0.9 is still met. The verdict stays "reported".
- 1b-1 to 1b-5 are unchanged: the score.json keys for recall, variants, unplaced peptides and
  mapping FDR are identical.
- Provenance: `$D/n1_fix/provenance.json`, `$D/n2_fix/provenance.json`; aggregate copy
  `validation/campaign/ecoli/provenance/N1N2_fix.json` (commit, image digests, input sha256s, job
  IDs, wall times, decoy counts, score.json sha256s).
