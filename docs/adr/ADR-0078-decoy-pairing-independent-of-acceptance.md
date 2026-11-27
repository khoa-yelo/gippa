# ADR-0078: Decoys count in the picked gene FDR whatever the targets' acceptance
- Status: accepted-by-khoa (2026-11-20, in session: "Fix now"; fix in 1.1.1, option A)
- Date: 2026-11-20
- Spec sections: §7.3.2 (decision 4, picked target–decoy tier), §8.1, §8.5; amends ADR-0009 (1/f
  retired from the picked formula); keeps ADR-0047, M4-12 decision 12, M7-12 decisions 23/28/29
- Needs Khoa: yes (changes a statistical rule and the meaning of the `decoy_win_weight` row)

## Context
GIPPA 1.1.0 computes the picked gene FDR as (D(t)/f + 1)/T(t) (ADR-0009). Here f is the share of
decoys passing q <= `input.q_value_max` whose stripped protein IDs all map to one gene. Unpaired
decoys enter f only and are dropped from D. ADR-0009 assumes they are lost **at random with
respect to the outcome**.

`protein_gene_map` maps a protein ID in three ways: (1) it is a C2 locus tag; (2) C12 rank-1 hits;
(3) a vote of **accepted** unique target peptides. When protein_hits is empty and the IDs are not
locus tags, only the vote applies. Then a decoy pairs only when its protein has an accepted target.
The decoys left unpaired are those whose protein has none. Those are exactly the decoys that would
have won their gene (a gene with no accepted target loses to any passing decoy). The loss depends
on the outcome, so 1/f does not repair it. The KS check in docs/methods/detection-tiers.md does
not see it: the q of unpaired decoys is not shifted, their partners are.

Campaign NCM3722 (PXD060434, DIA-NN, 304 decoys at q <= 0.01; issue draft
`local/review/campaign/ISSUE_decoy_pairing.md`). The same table was run twice:
- N1: matched genome, `--search-db same`, locus-tag IDs;
- N2: MG1655, `--search-db other --genome-relation related`, NCM3722 IDs, vote route.

N2 reports 2,399 genes at max q 0.0051. Its realised FDR, estimated with honest decoys, is about
0.026–0.028 at nominal 0.01. Pairing N1 by vote alone reproduces the N2 result (2,439 genes,
max q 0.0050).

Affected: every run on the vote route. That covers `--search-db other`, and `--search-db same`
with non-locus-tag IDs (e.g. UniProt accessions) and no protein hits. Locus-tag runs are not
affected (N1, campaign E1/E2).

## Decision
Option A. The rule changes how decoys are counted. It does not change how targets are accepted
or mapped.

1. **Every q-passing decoy with protein IDs competes.** `pair_decoys` keeps its row set
   (decisions 12, 28/29 and ADR-0047 below are unchanged). Each row gets a competition **unit**:
   - the gene, when every counted stripped ID maps to one and the same gene (as today, "paired");
   - otherwise the decoy's **own unit**. Its key is the sorted, de-duplicated list of the decoy's
     counted stripped IDs, joined by `;`. Counted IDs are those left after removing the longest
     decoy prefix, the contaminant IDs (decision 12) and the unmapped `GIPPAORF_` IDs
     (decision 28/29).
   Decoys with the same key share one unit, and its best decoy counts. Own units form a separate
   namespace: an own unit never takes a target, even if its key equals a C2 locus tag. Own units
   contribute to D only, never to T, and never appear in the gene table.
2. **Picked FDR without 1/f.** Units = eligible genes ∪ own units. Per unit the best target
   (genes only) competes with the best decoy, and a tie goes to the decoy. q is
   FDR(t) = (D(t) + 1)/T(t), made monotone, where D counts decoy wins in genes and in own units.
   An own unit has no target, so it is a decoy win whenever it exists. `picked_fdr_rescaled` is
   called with paired_share = 1. Decoys paired to a partial (non-eligible) gene leave the
   competition with that gene, as today. `q_floor` = 1/T is unchanged.
3. **ADR-0009's 1/f is retired from the formula.** Random decoy loss (unmappable IDs) now enters
   D in full as own units. That is conservative: each such decoy is counted as a win, where 1/f
   only counted it in expectation. Decoys with **no** protein IDs stay outside the competition
   and outside f, as today; f was never applied to them.
4. **The tier guard stays.** `decoy_paired_share` = gene-paired rows / all rows, as today. It
   still drives the 80 % guard (`decoy_pairing_low` → MAYU), so no input changes tier. It now
   measures how well the decoys fit the gene model; validity no longer depends on it.
5. **Summary rows.** `n_decoys_paired` and `decoy_paired_share` are unchanged.
   `n_decoys_unpaired` keeps its count; these decoys now compete in own units. `decoy_win_weight`
   is 1.0 whenever there is at least one pairing row, and NA otherwise. The row is kept for schema
   stability and documented as retired.
6. **Ranking key.** S = `shared_key_columns(accepted, <all pairing rows>)`: every competing
   decoy, own units included, in `detect` and in `choose_tier`. Before, only paired rows counted.
7. **Vote map unchanged.** `protein_gene_map` keeps voting with accepted unique peptides.
   `Detection.protein_gene_map` and the abundance step are unchanged. The lead asked whether
   pairing should vote with all C11 unique placements. No:
   - Under rule 1 the map no longer decides whether a decoy counts. It only decides whether the
     decoy competes in a gene or alone, and alone is always a win, which is conservative.
   - Votes from placements that failed acceptance (failed `placement_q`, tolerant route in an
     exact run, non-assigning locus class) are the noisier ones. A wrong pairing lets a decoy lose
     to an unrelated target, which is anti-conservative.
   - On N1 and N2 both maps give identical results (Evidence).
8. **Contaminants (decision 12), unchanged.** All-contaminant decoys get no row; a mixed decoy is
   unitised on its non-contaminant IDs. The known contaminant-overlap attack
   (`test_detection_tiers_attack_inputs.py`) is outside this ADR.
9. **ORF entries (M7-12 decisions 23/28/29), unchanged.** An all-`GIPPAORF_` decoy whose IDs map
   to no gene gets no row; it belongs to the M7-11 class FDR. One whose entry is voted to a gene
   pairs there.
   - Known residual: an `nterm_extension` entry whose in-gene targets have no accepted peptide maps
     to no gene, so its decoy is dropped (the same selection, in small form). The gene's own
     `rev_<locus_tag>` decoy still competes.
   - Backlog: map `GIPPAORF_nterm_extension` entries to their parent gene from makedb coordinates
     (database side), then pair them there.

Assumptions, stated in the methods note:
- (a) A decoy whose protein maps to no gene, or to several, is a null for some unit in which no
  accepted target would beat it. Counting it as a win can only overstate D.
- (b) A gene pairing made by votes of accepted peptides is correct, or at worst groups the decoy
  with a target of the same protein.
- (c) Decoys without protein IDs are rare or absent; an input with many of them is outside the
  rule's guarantee, as before.

## Evidence (before → after: metrics, tests, simulations)
Campaign NCM3722, same inputs, `decoy_picked` at α 0.01. Scripts:
`local/scratch/adr-0078/r1.py` (load), `r2.py` (diagnosis), `r3.py` (this rule, both vote maps).

| run | rule | decoys gene-paired / own unit | weight | decoy wins | genes q <= 0.01 | max q |
|---|---|---|---|---|---|---|
| N1 (`same`, locus tags) | 1.1.0 | 301 / 3 dropped | 1/0.990 | 62 | 2,341 | 0.0099 |
| N1 | ADR-0078 | 301 / 3 own units | 1 | 64 | 2,341 (same set) | 0.0098 |
| N2 (`other`, vote route) | 1.1.0 | 244 / 60 dropped | 1/0.803 | 9 | 2,399 | 0.0051 |
| N2 | ADR-0078 | 244 / 60 own units | 1 | 66 | 2,305 (94 fewer, none new) | 0.0100 |

- N1: the detected set is unchanged. Max |Δq| over all genes is 0.0007. Two of the three own
  units win, but below the last call; the dropped 1/0.99 weight offsets them. The rule's effect on
  locus-tag runs is nil here, and conservative in general.
- N2: the 94 lost genes are the ones N1 rejects (N1 q 0.010–0.026). With all C11 unique
  placements as the vote source, both runs give the same numbers (vote map: 2,389 proteins either
  way).

Simulation: `detection_tiers_calsim` replicate, pi_wrong 0.05, every protein ID re-labelled to an
accession not in C2, protein_hits empty. Mean FDP ± SE at q = 0.01; "fixed" is a reference
implementation of rules 1–2 (`local/scratch/adr-0078/a2.py`, `a3.py`; seeds
`[20261120, 78, …, rep]`).

| input | R | gene-paired share | 1.1.0 | ADR-0078 |
|---|---|---|---|---|
| 4,500 genes, all IDs re-labelled | 50 | 0.895 (min 0.864) | 0.01275 ± 0.00051 | 0.01006 ± 0.00038 |
| 2,000 genes, all IDs re-labelled | 50 | 0.894 (min 0.836) | 0.01220 ± 0.00072 | 0.01002 ± 0.00062 |
| 2,000 genes, locus tags (baseline) | 40 | 1.000 | 0.01019 ± 0.00070 | identical |
| 2,000 genes, 20 % random loss (ADR-0009) | 40 | 0.802 | 0.00965 ± 0.00071 | 0.00705 ± 0.00057 |

At q = 0.05 (4,500 genes): 0.0530 ± 0.0005 under 1.1.0 against 0.0494 ± 0.0007 under ADR-0078.
q = 0.10 is saturated by the peptide q cap: both are equal and far below the bound. The tier is
`decoy_picked` in every replicate.

## Spec for tests
Notation: C2 genes G1–G4, non-partial, all eligible. protein_hits is empty. Prefixes are the
defaults (`rev_`, `CON__`). Key S = {q}: targets and decoys carry q only, pep and score NA.
Accepted unique targets, one peptide per gene, sit on proteins P1–P4 voted to G1–G4, with q 0.001,
0.002, 0.003 and 0.004. Every decoy below passes q <= 0.01. Expected gene q comes from
`picked_fdr_rescaled(..., 1.0)` over the units.

Unit tests (test-author, `tests/unit/genetable/`):
1. **Target-less decoy, vote route (the defect).**
   - Decoys: `[rev_P2]` q 0.0025 and `[rev_P9]` q 0.0015. P9 has no target.
   - `pair_decoys` gives 2 rows: `gene` G2 / NA and `unit` G2 / `P9`.
   - Gene q is 0.5 for G1–G4; 1.1.0 gives 0.25 (f 0.5, D 0).
   - Summary: `n_decoys_paired` 1, `n_decoys_unpaired` 1, `decoy_paired_share` 0.5,
     `decoy_win_weight` 1.0. The tier falls to MAYU at the default guard, so either run picked
     with `min_identity` = `min_coverage` = 0.5 or test `_picked_q` / the unit frame directly.
2. **Multi-gene and partly mapped decoys.**
   - Decoys: `[rev_P2]` q 0.0025, `[rev_P1;rev_P2]` q 0.0015 (unit `P1;P2`), and
     `[rev_P7;rev_P1]` q 0.0035 (P7 unmapped; unit `P1;P7`, sorted).
   - Gene q: G1, G2, G3 0.666667; G4 0.75. 1.1.0 gives 0.25 for all.
3. **Same key twice.** Decoys `[rev_P9]` q 0.0015 and `[rev_P9]` q 0.0030 form one unit P9, and
   only its best decoy counts. Gene q is 0.5 for all genes, as in example 1.
4. **Contaminants unchanged.** `[rev_CON__X]` gets no row. `[rev_P9;rev_CON__X]` gets a row with
   unit `P9`. `[rev_P2;rev_CON__X]` pairs to G2.
5. **ORF entries unchanged.** `[rev_GIPPAORF_intergenic_c1_100_200_p]` (unmapped) gets no row.
   An `nterm_extension` entry voted to G1 pairs to G1 with unit G1. The existing
   test_detection_orf_decoys*.py pins stay green.
6. **Locus-tag run unchanged.**
   - IDs are C2 locus tags; decoys `[rev_G2]` q 0.0025 and `[rev_G3]` q 0.0005.
   - Gene q: G1, G2, G4 0.666667; G3 NA (a decoy win). Identical to 1.1.0, byte for byte, for any
     input where every decoy pairs.
7. **Namespace.** An own unit never takes a target. Unit `P1;P2` stays separate even when a C2
   gene is named `P1;P2`.
8. **ADR-0047 unchanged.** A decoy with q 0.02 gets no row when some target carries q. When no
   target carries q, it gets a row.
9. **Random loss.** `[rev_LOST_0]` gets unit `LOST_0` and counts as a decoy win.
10. **Ranking key over all rows.** When only an own-unit decoy has NA score, score leaves S and
    `decoy_unscored` warns (before, only paired rows were read).
11. **Partial gene.** A decoy paired to a partial gene leaves the competition, as in 1.1.0.

Calibration attack (statistician-specified, test-author writes;
`tests/calibration/attacks/test_detection_tiers_attack_decoy_pairing.py`):
- **Data.** The `gene_fdr` replicate of `tests/helpers/detection_tiers_calsim.py` at pi_wrong 0.05.
  It runs at **4,500 genes in both modes** (fast and full), because 2,000 genes lack power in the
  fast mode.
  - Every target protein ID `GSIM_<k>` and its decoy `rev_GSIM_<k>` is re-labelled to
    `ACC<perm(k):06d>` / `rev_ACC<...>`, with a permutation drawn from the replicate rng. C11
    placements (`gene`) stay as they are: targets still place by sequence.
  - protein_hits is empty (`cs.protein_hits`), so the map is vote-only, as in N2.
  - Sim stream `gene_fdr.attack_decoy_pairing.sim.reid100` (+ `.full`). Claim stream
    `gene_fdr.attack_decoy_pairing.reid100`.
- **Replicates.** `C.replicates("gene_fdr", full)` (50 fast / 300 full), ADR-0015 two-block
  confirmation.
- **Criterion.** `fdr` (mean FDP <= q + `se_multiplier` · SE, k = 2 from tolerances.yaml) at each
  q in `cs.LEVELS`, with detection rerun at `stats.alpha.identification` = q. No new tolerance
  key; nothing looser than the existing gene_fdr claims.
- **Assertion.** The tier is `decoy_picked` in every replicate (share >= 0.86 in 50/50
  simulated). If a replicate falls to MAYU, the test must fail, not skip.
- **Teeth.** Required before the implementer starts:
  - On 1.1.0 the fast claim at q = 0.01 must be a confirmed fail. Expected about
    0.0128 ± 0.0005 against a bound of about 0.0110, roughly 3.4 SE over the bound in each block.
  - q = 0.05 is expected to fail too.
  - Mark the claims `xfail(strict=True)` until the fix lands.
- **Fixed-code expectation.** About 0.0101 ± 0.0004 at q = 0.01. The chance of a false
  two-block failure is about 5e-4.
- **Regression.** The existing `test_detection_tiers_attack_decoy_loss.py` must keep passing. It
  becomes more conservative: 0.0071 against 0.0097 at 20 % loss, and the tier assertions are
  unchanged.

## Consequences
- 1.1.1 patch, src/gippa/genetable/detection.py only:
  - `pair_decoys` gains a `unit` column;
  - `_picked_q` competes over genes ∪ own units with weight 1;
  - S is read over all pairing rows;
  - the summary rows change as in rule 5.
  - Changelog `fixed` entry.
- Vote-route runs lose anti-conservative calls (N2: −94 genes). Locus-tag runs keep their calls.
  Their q may rise slightly when unmappable decoys win high in the ranking.
- docs/methods/detection-tiers.md is corrected in this ADR's PR: the cross-strain paragraph, the
  decoy pairing and the picked formula. The implementer's PR rewords the report block
  `decoy_picked` ("divided by the share of decoys that paired" goes) and regenerates the snippets
  with `scripts/dev/sync_report_snippets.py`, because the snippets live under src/.
- The 1.1.0 pilot arm A1 (UniProt K-12 database) used the vote route. Re-check its detection
  counts after the fix.
- Backlog: database-side pairing for `GIPPAORF_nterm_extension` (rule 9). The KS diagnostic named
  in the methods note is withdrawn as a guard for this failure; it cannot see it.
