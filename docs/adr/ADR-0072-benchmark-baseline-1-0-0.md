# ADR-0072: The 1.0.0 benchmark: two repeats, measured runtimes, and the §13.2 range waived
- Status: accepted-by-khoa (2026-10-31 and 2026-11-01, in session: "Accept measured runtime (Recommended)"; "Cut repeat 3")
- Date: 2026-11-01
- Spec sections: §13.2, §15.2, §15.6; M10-12, M10-15 bullet 4; `scripts/release/readiness.py` rows R7, R8
- Needs Khoa: yes (it supersedes a spec number and a repeat count)

Waives readiness R7

This is the decision record M10-12 calls ADR-M10-D.

## Context
§13.2 estimates a K-12 standard run on 16 threads and 64 GB at about 35 min (range 16-55 min). §15.2
asks for the median of 3 repeats. The benchmark at commit `fa49b8f1` ran on wheat-01, array 1916711.
Each repeat took about 14 h for its four cells, and later release steps had to wait for the array.

## Decision
1. **The standard run** is `scripts/slurm/benchmark_standard.sbatch`: 16 threads, 64 GB, one node
   class (`wheat`), `farmshare_local,apptainer`, Bakta's full database for the standard variant,
   on the K-12 and `sim5m` inputs the script names. Both use simulated peptide tables (`gippasim`).
2. **Two repeats, not three, for 1.0.0.** Repeat 3 was cancelled after the first repeats showed the
   range was missed by far more than the repeats differ. Each cell's `median_s` is the mean of its 2
   runs, recorded as `rounds: 2` in `benchmarks/baseline.json`. Later baselines return to 3.
3. **The measured runtimes replace the §13.2 estimates in the docs.** R7 is waived for 1.0.0, and the
   16-55 min range is not met. No bound changes.

## Evidence
| Cell | Repeat 1 | Repeat 2 | Mean | Spread | Peak memory, both repeats (GiB) |
|---|---|---|---|---|---|
| K-12, standard | 3:28:22 | 3:52:43 | 3:40:32 | 11.0 % | 62.4, 49.2 |
| K-12, core | 3:01:36 | 3:00:51 | 3:01:13 | 0.4 % | 26.4, 37.9 |
| `sim5m`, standard | 4:02:58 | 3:27:01 | 3:44:59 | 16.0 % | 62.4, 62.1 |
| `sim5m`, core | 3:22:26 | 3:35:45 | 3:29:05 | 6.4 % | 39.6, 41.2 |

All 8 runs exited 0 with no invalidating flag. `sim5m` standard is above the 15 % noisy threshold of
`ci/gates.yaml`, and the summary marks it noisy. With 2 repeats a spread is a weak estimate, and the
true variation may be larger: the cells differ from each other by more than the first K-12 core
pair suggested. The node was busy: its median allocation was 60-80 of 80 CPUs.

Two steps take about 92 % of summed task time in K-12 standard's repeat 1: `promoter_calculator` and
`terminatornet_chunk`. That is one run's breakdown, summed across parallel tasks.

Not established: why peak memory differs so much between repeats of one cell, and whether the 62 GiB
peaks reflect a limit. Treat 64 GB as the tested allocation, not as a measured requirement.

## Consequences
- The runtime promise in §13.2 is not kept in 1.0.0; the docs state the measured times. Users on
  similar hardware should expect about 3-4 h on 16 threads.
- The waiver applies to 1.0.0 only. A later release either speeds up the two slow steps or changes
  the claim, with 3 repeats.
