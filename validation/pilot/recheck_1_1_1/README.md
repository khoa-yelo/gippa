# Pilot addendum: 1.1.1 detection recheck (ADR-0078)

This is the full addendum to `validation/pilot/REPORT.md` for the ADR-0078 fix (PR #607, merged as
92fbf3af). `REPORT.md` is rendered from `metrics.json` and `report_notes.md`. Its render guard rejects
notes that hold a number of 4 or more digits not found in `metrics.json`, so the section there,
"Addendum: 1.1.1 detection recheck", gives only the changes and verdicts. The values with 4 or more
digits are in this file and in `../metrics/metrics_1_1_1.json`. `REPORT.md` and `metrics.json` keep
the 1.0.0 freeze values, and the values below supersede them for 1.1.1. Only aggregates are recorded
here: no sequence, no gene list, no per-gene value.


GIPPA 1.1.0 left out of the picked gene FDR any decoy whose own unit had no q-passing target, so on the `decoy_picked` route such decoys never competed. ADR-0078 (PR #607) lets every q-passing decoy compete. Every pilot `decoy_picked` arm was recomputed with the 1.1.0 code (a2a95c76) and with the fix (7d4c4bfd). Then gene_table, detectability, context_stats, function and compare were rerun on the freeze outputs with the freeze code and only `detection.py` swapped: Slurm 2001992 is the control (it reproduces the recorded tables byte for byte), and Slurm 2002006 is the fix. Inputs, scripts and aggregate outputs: `validation/pilot/recheck_1_1_1/`. The superseding values, under the key paths of `metrics.json`: `validation/pilot/metrics/metrics_1_1_1.json`.

The 1.0.0 values above stay as recorded, and `metrics.json` is unchanged. For a 1.1.1 reading, the values below supersede them. No threshold changed. No verdict changed: H1 (not_supported, "mostly" fails), H3 (supported), TE3, TE4, and the plasmid and prophage coverage.

Gene detections at q ≤ 0.01, `decoy_picked` route (decoys paired by gene; own units: unpaired decoys that compete as their own unit):

| arm | genes | max q | decoys paired | own units | decoy wins |
|---|---|---|---|---|---|
| pass 1 A1 (pre-390 and #390) | 2,670 → 2,642 | 0.0022 → 0.0098 | 189 | 0 → 33 | 4 → 37 |
| pass 1 A2 (pre-390 and #390) | 3,039 → 3,034 | 0.0096 → 0.0086 | 261 | 0 → 2 | 41 → 43 |
| pass 3 A1, freeze A1 | 2,673 → 2,645 | 0.0022 → 0.0098 | 189 | 0 → 33 | 4 → 37 |
| pass 3 A2, freeze A2, freeze A2 replicate | 3,057 → 3,051 | 0.0095 → 0.0098 | 323 | 0 → 2 | 49 → 51 |
| A4 (six-frame, exploratory) | 2,965 → 2,965 | 0.0055 → 0.0067 | 168 | 0 → 4 | 15 → 19 |

A1 searches K-12 protein accessions, so 38 of its decoys have no gene pair; 33 of them now compete as their own units. That moves A1's maximum q from 0.0022 to 0.0098, still below 0.01.

Superseded report values (1.0.0 → 1.1.1). Rows of metric 3 not listed here are unchanged.

| report value | 1.0.0 | 1.1.1 |
|---|---|---|
| Metric 3, A1 all | 0.5705 [0.5563, 0.5847] | 0.5646 [0.5503, 0.5787] |
| Metric 3, A1 chromosome | 0.5719 [0.5577, 0.586] | 0.5659 [0.5517, 0.5801] |
| Metric 3, A1 bgc | 0.3108 [0.2567, 0.3705] | 0.3028 [0.2493, 0.3623] |
| Metric 3, A1 k12_shared | 0.7111 [0.6964, 0.7254] | 0.7036 [0.6888, 0.718] |
| Metric 3, A2 all | 0.6525 [0.6388, 0.666] | 0.6512 [0.6375, 0.6647] |
| Metric 3, A2 chromosome | 0.6522 [0.6384, 0.6657] | 0.6509 [0.6371, 0.6644] |
| Metric 3, A2 strain_specific | 0.3461 [0.316, 0.3775] | 0.345 [0.3149, 0.3764] |
| Metric 3, A2 k12_shared | 0.7289 [0.7145, 0.7429] | 0.7276 [0.7131, 0.7416] |
| Metric 4, Spearman rho | 0.978 [0.9731, 0.9826] | 0.978 [0.9728, 0.9827] |
| Metric 4, shared genes | 2635 | 2613 |
| Metric 4, Bland–Altman bias (log2) | -0.195 | -0.194 |
| Metric 5, enriched categories Jaccard | 0.563 | 0.588 |
| H1, A2-only genes b / A1-only genes c | 417 / 33 | 434 / 28 |
| H1, A2-only share of discordant genes | 0.9267 [0.8988, 0.9473] | 0.9394 [0.9138, 0.9577] |
| H1, coverage difference | 0.082 | 0.0867 |
| H1, TU-level units (mixed, left out) | 275 against 31 (1) | 290 against 24 (3) |
| H1, p_tu | 2.5e-50 | 1.8e-59 |
| H1, unit effect | 0.8987 [0.8598, 0.9277] | 0.9236 [0.8888, 0.9481] |
| H1, p_gene | 4.6e-86 | 5.2e-95 |
| H1, share of A2-only genes in plasmid, prophage or BGC | 0.2542 [0.2148, 0.2981] (106 of 417) | 0.2442 [0.2062, 0.2868] (106 of 434) |
| H1, strain-specific share of A2-only genes | 0.7602 [0.717, 0.7987] | 0.7281 [0.6844, 0.7678] |
| H3, rho criterion | 0.978 [0.9731, 0.9826], holds | 0.978 [0.9728, 0.9827], holds |
| H3, enriched categories n_a / n_b / n_shared / union / Jaccard | 437 / 285 / 260 / 462 / 0.563, holds | 405 / 292 / 258 / 439 / 0.588, holds |
| Read these with the numbers, H1 | 106 of the 417 A2-only genes | 106 of the 434 A2-only genes |

H1's base rate (0.08538) and its failed clause are unchanged, so H1 stays not_supported. H3's rho and enriched-category criteria still hold, so H3 stays supported.

## Files

| file | what |
|---|---|
| `recheck.py` | recomputes `decoy_picked` gene detections (aggregates) with whatever gippa is on `PYTHONPATH` |
| `out_v110.jsonl`, `out_fix.jsonl` | its output per run, 1.1.0 code (a2a95c76) and the fix (7d4c4bfd) |
| `rerun_steps.py` | reruns gene_table → detectability → context_stats → function of one freeze arm |
| `rerun.sbatch`, `rerun_fix.sbatch` | Slurm 2001992 (control, freeze code) and 2002006 (freeze code with the fix's `detection.py`), each followed by `gippa compare` |
| `derived.py` | metric 3 (every subset), H1, metric 4, metric 5 and H3 from the rerun tables, with the freeze `validation/pilot/metrics.py` functions |
| `derived_ctl.json`, `derived_fix.json` | its output. The control equals `metrics.json` on every key it covers. |
| `build_metrics.py` | writes `../metrics/metrics_1_1_1.json` (superseding values, per-arm recheck and the list of changed keys). It asserts that the control equals `metrics.json`. |

The scripts are copied as they ran, with scratch paths under `local/scratch/a1-recheck` of the main checkout.
There are two changes: `ruff format` and a file-level `noqa` line, and `derived.py` now writes every
metric 3 subset, not 6. The values of the 6 are unchanged. Per-gene detection lists stayed in scratch and are not committed.
