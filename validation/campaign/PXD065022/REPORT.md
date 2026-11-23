# PXD065022 report: Treponema pallidum, 7 strains (Phase 2, group A, dataset 1), interim: blocked

Checks pre-registered in `ANSWER_KEY.md` (commit `854c143c`, before any search). Aggregates only.
Machine-readable numbers so far: `metrics_partial.json` (`p2a/tpa_partial.py`).
`$D` = `$GIPPA_ROOT/campaign/PXD065022`. **No `gippa run` has completed for any strain**, so every
check that reads GIPPA output is open.

## Status

| Step | State |
|---|---|
| Raw files | 21 runs (3 per strain) unpacked to `.d` (prep job 2005458, zip CRC test). `PHIL_A` failed the CRC test (corrupt member `analysis.tdf_bin`); one refetch (job 2006167) passed. Old copy sha256 `d7a07bde…` kept as `.corrupt_crc`; new sha256 `cd2e9aa9c51c3f283095579e72d6fc0aa592b39c014aa1c5012453b87820004f`. `PREP.md` / `files.tsv` updated to the new sha256 |
| Host | UniProt UP000001811 (release 2026_04), 21,205 sequences, MD5 `97516cae5c749aeb11cf5887864ab31a` = `RELEASE.metalink`; headers `CONT_RABIT_<acc>` |
| makedb | all 7 strains, A: 1,353-1,357 records (≈ 573 kB); B: 31,751-31,759 records (3.7 MB) (job 2005458) |
| DIA-NN route A | done: DAL-1, Grady, Haiti B, Madras, Mexico A (jobs 2006168, 2009201); Philadelphia 1 and SS14 running (L6, job 2013337) |
| DIA-NN route B | done: DAL-1 (2006168), Grady (2011438); Haiti B, Madras, Mexico A, Philadelphia 1, SS14 running (L4, job 2011438) |
| convert + preflight | `convert` exit 0, byte-identical twice; `preflight --strict` exits 3 on every strain because of the warning `genome.gaps` alone (assembly N runs). Preflight without `--strict` exits 0 (deviation below) |
| `gippa run` at `bfc74ff6` | DAL-1 B: exit 4, `ANNOTATE:PROMOTECH` crash on 40 windows with N → **#611** (fixed in #612, main `cff0fd68`) |
| `gippa run` at `cff0fd68` (`rel111b-pin`, new outdirs `out_*_fix`) | DAL-1 A, DAL-1 B, Grady A: exit 4 after 25-31 min, `ANNOTATE:PROMOTER_CALCULATOR` `KeyError` on a non-ACGT base in TU upstream sequences (DAL-1: IUPAC `R`; Grady: `N`) → **#614** (open). Every T. pallidum assembly has N runs (64-324 N; DAL-1 also 1 `R`), so all 14 runs are blocked. The run lane was cancelled (job 2012859) |

## Numbers available now

| # | Check | State | Number |
|---|---|---|---|
| T1 | Host share ≥ 50 % in every strain (distinct accepted route-A target peptides with only `CONT_` proteins) | 5 of 7 strains measured, all above 50 %; verdict open until Philadelphia 1 and SS14 are in | DAL-1 16,160/24,510 = 65.9 % [65.3, 66.5]; Grady 9,371/12,117 = 77.3 % [76.6, 78.1]; Haiti B 23,018/32,302 = 71.3 % [70.8, 71.7]; Madras 11,612/16,532 = 70.2 % [69.5, 70.9]; Mexico A 11,586/17,153 = 67.5 % [66.8, 68.2]. Rabbit-only share 55.0-64.3 %; contaminant-library-only 6.9-9.4 % |
| T6 | Truth sets vs DAL-1 (blastp RBH ≥ 90 % identity, ≥ 80 % coverage; job 2005467, built before any cross-strain run) | built | pairs / sites / proteins with sites: Grady 962/218/96, Haiti B 962/82/51, Madras 966/26/21, Mexico A 961/247/93, Philadelphia 1 964/221/96, SS14 964/224/97 |
| — | DIA-NN route A (5.30 M precursors: rabbit + strain; 16 CPUs) | | protein groups at global q ≤ 0.01: DAL-1 4,322, Grady 3,341, Haiti B 4,702, Madras 3,915, Mexico A 3,726 (rabbit + treponemal). Precursors at 1 % per run 3,080-33,704. Step 2: 2 h 50 min-3 h 40 min per strain |
| — | DIA-NN route B (6.20 M precursors) | | DAL-1: 4,257 protein groups, step 2 5 h 00 min; Grady: 3,286, 3 h 26 min |

Open: S1-S11, T2-T5, T6 scoring and T7 (TE7 on DAL-1 route A).

## Resume once #614 is fixed

New pin, then for every strain and route:
`sbatch --export=ALL,ARM_SCRIPT=arm_v3.sh,PIN=<pin> -J gippa-camp-p2a-L* --time=47:00:00 validation/campaign/p2a/lane_v2.sbatch "PXD065022 <tag> a" ...`.
DIA-NN is skipped when its report exists. Convert, preflight and run write to `in_<r>_fix` /
`out_<r>_fix`; set `SUF` to a new suffix if `_fix` is taken. Then run the cross-strain arms
(`"PXD065022 <tag> x tpa_dal1"`), TE7 (`"PXD065022 tpa_dal1 te7"`), `p2a/score.py` per strain and
`p2a/p2_metrics.py --host`.

## Deviations and problems

- **Preflight strict waiver.** `--strict` turns any warning into exit 3. The only warning is
  `genome.gaps` (N runs in the RefSeq assemblies themselves, which no permitted input change can
  remove). `p2a/convert_preflight.sh` then reruns preflight without `--strict`, requires exit 0 and
  records `preflight_strict_exit` 3, `preflight_exit` 0 and `strict_waived` `genome.gaps` in
  `convert_provenance.json`. Any other warning or error still stops the arm.
- **Host as contaminant.** Rabbit proteins enter DIA-NN as `CONT_RABIT_*` records, so GIPPA's default
  contaminant prefixes bucket them; no GIPPA setting was changed (pre-registered).
- **Bugs:** #611 (Promotech, fixed), #614 (Promoter Calculator, open). Both were opened by this agent.
- DIA-NN search cost is dominated by the rabbit library (5.3 M precursors): about 3-5 h per strain
  and route on 16 CPUs, against 12-63 min for the other two datasets.

The full report replaces this one when the runs complete.
