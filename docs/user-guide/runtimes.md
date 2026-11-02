# Runtimes

Measured wall-clock times of a `gippa run` from `scripts/slurm/benchmark_standard.sbatch` (M10-12).
These numbers replace the estimates of the specification. The source is
`benchmarks/baseline.json › standard_run`, and the decision record is
[ADR-0072](../adr/ADR-0072-benchmark-baseline-1-0-0.md).

**This 1.0.0 baseline has 2 repeats per run, not 3.** The third repeat was cancelled (ADR-0072), so
each time below is the mean of 2 runs, and the spread is (slowest − fastest) / mean. With 2 repeats a
spread is a weak estimate, and the true variation may be larger.

| Run | Mean of 2 runs (h:mm) | Range of 2 runs | Spread | Peak memory, both runs (GiB) | Threads / memory | Node class | Date | Commit |
|---|---|---|---|---|---|---|---|---|
| E. coli K-12, standard (all modules, live Bakta full DB) | 3:41 | 3:28:22-3:52:43 | 11.0 % | 62.4, 49.2 | 16 / 64 GB | wheat | 2026-10-31 | fa49b8f |
| E. coli K-12, `--modules core` (recorded Bakta JSON) | 3:01 | 3:00:51-3:01:36 | 0.4 % | 26.4, 37.9 | 16 / 64 GB | wheat | 2026-10-31 | fa49b8f |

How to read the table:

- Plan for about 3-4 h on 16 threads on similar hardware. The specification's first estimates were
  not met, and ADR-0072 waives that range for 1.0.0.
- 16 threads are 16 hardware threads (8 physical cores) of an Intel Xeon Gold 6148 node, shared
  with other users' jobs; the node's median allocation was 60-80 of 80 CPUs.
- Peak memory is the run's Slurm MaxRSS (`sacct`, sampled every 60 s). Why it differs so much
  between repeats is not established. Each run's Slurm limit was 64000 MiB (62.5 GiB), and the 62.4 GiB
  peaks sit just under it, so they may reflect that limit rather than a need; this is not established
  either. Treat 64 GB as the tested allocation, not as a measured requirement.
- Both runs used simulated peptide tables. The `core` run differs from `standard` in two ways, fewer
  modules and no Bakta run, so it is not "standard minus modules".
- The baseline's simulated 5-Mb genome is the speed-loop reference, not a user number, so it is not
  in this table. Its times are in ADR-0072 and `benchmarks/baseline.json`; its standard run is marked
  noisy (spread 16.0 %, above the 15 % threshold).

What affects the runtime, qualitatively:

- `--publication` raises permutations and bootstrap draws from 1,000 to 10,000;
- `--promoters genome` predicts promoters genome-wide instead of in upstream windows;
- Bakta runs unless `--annotation bakta.json` is given, and the full database is slower than the
  light one;
- the optional modules: antiSMASH (`bgc`) and geNOMAD (`mobile`).

Nextflow writes per-task times and peak memory to `OUT/logs/trace.tsv` and
`OUT/logs/nf_report.html`. The manifest records per-step time and peak memory.
