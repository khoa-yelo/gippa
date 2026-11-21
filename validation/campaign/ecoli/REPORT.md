# Phase 1 report: E. coli K-12 MG1655 demonstration (arms E1-E4)

Demo and discovery campaign, Phase 1 (handover plan of 2026-11-17). The pre-registered answer
key is `ANSWER_KEY.md` (committed before any search, not edited). This report holds aggregates
only: counts, shares, q-values, locus tags and 1-based coordinates. No peptide sequence and no
DIA-NN row is in git. Machine-readable numbers: `metrics_phase1.json` (written by
`phase1_metrics.py`); per-arm provenance: `provenance/E1-E4.json` (copies of
`$GIPPA_ROOT/campaign/ecoli_k12/<arm>/provenance.json`).

All shares carry a Wilson 95 % CI in brackets. "unestimated" means GIPPA or this report gives
no error estimate for that number.

## Verdicts

| # | Check | Verdict | Number (error estimate) | Cause / note |
|---|---|---|---|---|
| 1 | Mapping, TE1 (≥ 99.5 % of E1 peptides exact in frame) | **pass** | 31,071 / 31,071 distinct accepted E1 peptides exact in frame at their own gene: 100 % [99.988, 100] | 38 peptides with a non-ATG N-terminus, all mapped; 23 peptides in `rna_overlap` at their own gene were masked by the pilot's te_proxies rule. Placement FDR: not applicable to E1 (`--search-db same`, exact hits accepted as `same_db`; `mapping_fdr.tsv` holds no accepted placement rows) |
| 2 | Expression, TE3 (ribosomal ≥ 90 % detected; median abundance in the top quartile) | **pass** | 53 / 53 ribosomal proteins detected, 100 % [93.2, 100]; median iBAQ log2 29.25 vs Q3 22.97 (percentile rank 0.989) | §7.1.2 set: 59 / 59 [93.9, 100] |
| 3 | Stoichiometry, TE4 | **fail** (1 of 2 pairs) | log2(RplL/RplJ) = −0.33 (iBAQ; top3 −0.12), outside [1, 3]; log2(AtpA/AtpD) = −0.30 (top3 −0.21), inside [−1, 1]. Error unestimated | Same failure as the Nissle pilot (waived for 1.1.0, ADR-0076); the answer key named it a known risk. Cause not established. Candidates: L7 is N-acetylated L12, and N-terminal acetylation is not a search modification here, so part of the L7/L12 signal is not searched; the 12 kDa protein has 7 observable peptides (rplJ 8), so iBAQ rests on few peptides. Not a GIPPA bug as far as this run shows |
| 4 | Contaminants, TE5 | **pass** | 196 / 196 trypsin (4) and keratin (192) peptides in the contaminant bucket, 100 % [98.1, 100]; 0 placed on the genome | Every E1 unmapped peptide is a contaminant (217 `no_hit`, 59 `too_short`, 1 `rejected_placement_fdr`) |
| 5 | Genome enrichment: cryptic prophages detected less than detectability predicts | **fail** | 5 geNomad prophage regions, all observed < expected (O/E 0.14-0.51); q ≤ 0.05 in 2 of 5 (DLP12 q 0.036, Rac q 0.013); e14 q 0.056, region at 1.63 Mb q 0.090, region at 2.46 Mb q 0.082 | Rule fixed before reading (in `phase1_metrics.py`): every region O < E **and** q ≤ 0.05. Direction holds in 5 / 5; significance does not in 3 / 5 (small regions: 18-46 genes; per-region circular-shift test). CP4-6 and CP4-57 were not called as prophage regions by geNomad, so they were not tested. No island (`islands.tsv` empty) |
| 6 | Operons: ribosomal and ATP synthase partners concordant | **pass** | 37 adjacent pairs among 62 ribosomal-protein/atp genes (33 OP, 4 NO; 36 both detected); 0 with within-operon q ≤ 0.05 (min q 0.78); 0 discordant `operon_flags` rows | One `conflict_terminator_in_operon` flag (atpA-atpG) with a concordant abundance pattern (q 0.996). Median |Δlog2| 0.73 |
| 7 | Crl | **pass** (junction part not testable) | E1: Bakta crl (KOIAIB_00263, 258,734-259,006 +) detected, 7 peptides, gene q 4.0e-4. E2: 1 peptide in 258,676-258,733 (`nterm_extension` candidate NOVEL_U00096.3_258662_p); 2 peptides in fragment 1 (257,829-257,899), candidate NOVEL_U00096.3_257829_p, class `intergenic`, 257,829-257,912 +, RBS and promoter supported. Not called (class q floor, row 15) | The junction peptide cannot be identified by a database search: no database (routes A or B, both built from the IS1-carrying reference) contains it, so it can never reach `unmapped.tsv`. E4 (de novo) did not place a peptide across the junction either |
| 8 | PrfB | **pass** | E1: Bakta prfB (KOIAIB_03083) detected, 13 peptides. E2: 4 peptides in 3,036,066-3,036,206 (−) as `nterm_extension` (NOVEL_U00096.3_3035184_m, flag `no_start_codon`); 1 peptide in the 25-codon leader (3,036,208-3,036,282) as an `intergenic` candidate (NOVEL_U00096.3_3036205_m, 25 aa). Neither called | As expected for the upstream peptides (N-terminal extension). The other-frame leader landed in an `intergenic` candidate rather than the `cds_outframe` class the answer key named; it was not called, so the verdict holds. GIPPA reports where the peptides land, not a frameshift |
| 9 | Novel ORFs vs RefSeq (share of called ORFs matching a gene or pseudogene) | **not testable** | 0 ORFs called (124 candidates) | No call to compare. Exploratory, all 124 candidates: 23 / 124 overlap a RefSeq CDS or pseudogene on the same strand by ≥ 50 % (18.5 % [12.7, 26.3]); by class: nterm_extension 7/11, intergenic 8/21, cds_outframe 8/30, antisense 0/62. Candidates with ≥ 2 peptides: 3 / 6. These shares are low by construction: `phase1_metrics.py` counts every candidate without a start codon as a non-match, which includes the prfB and infC `nterm_extension` candidates that lie inside RefSeq prfB and infC (candidate table below) |
| 10 | Selenoproteins | **pass** | FdnG (KOIAIB_01574), FdoG (KOIAIB_04147, Bakta names it fdnG) and FdhF (KOIAIB_04360) all detected; peptides on both sides of the Sec TGA in each (fdnG 5 / 23, fdoG 49 / 11, fdhF 1 / 3, genomic left / right); 0 peptides across the TGA | makedb writes Sec as U in the protein, so the Sec tryptic peptide is not searched; this is consistent with 0 across |
| 11 | Negative controls (E3) | **dropped** (Khoa, 2026-11-20) | 0 knockout runs identified | Not pursued further (Suppl. sample key not chased). PXD066241 (Panorama Public) holds 18 DIA runs without a sample key. The 8 whole-lysate-depth runs all carry Lon (22-38 peptides) and ClpP (2-7 peptides) and no CAT, TetR(B) or TetA(B) cassette peptide (one APH(3')-Ia peptide in 2 runs), so none is Δlon::tet or ΔclpP::cm; the 10 other runs are 30 kDa MWCO flow-through. The paper measured ΔmntS by PRM only. Strain assignment used cassette proteins only; Lon/ClpP counts were read afterwards as a check |
| 12 | Variants, same strain, TE2 (E4) | **fail** | 3 `high` ALT calls of 295 ALT sites, 0 at start/stop codons → 3 violations. Error unestimated (`class_error_unestimated`: DIA-NN writes no spectrum IDs, so no spectrum pairing) | carA (KOIAIB_00034) codon 30,668 + aa 340, swap, `mass_equivalent`; pckA (KOIAIB_03612) codon 3,534,432 + aa 539, swap, `mass_equivalent`; rpmG (KOIAIB_03867) codon 3,811,352 − aa 22, T>S SNV, call `mixed` (1 REF, 2 ALT). Same genome as the spectra, so all three are de novo errors or modifications, not strain variants. Two are mass-equivalent residue swaps, which TE2 does not exempt. Candidates for review below; not a confirmed GIPPA bug |
| 13 | Modification-like calls (E4) | **pass** | 249 / 295 ALT sites flagged `mass_equivalent` (170 swap, 41 substitution, 38 block; 247 candidate, 2 high), listed separately in `variants.tsv` | All 249 carry the flag; GIPPA does not count them as strain variants. Row 12's two mass-equivalent `high` calls are among them |
| 14 | Route B cost | **reported** (no threshold) | E2 keeps 2,648 / 2,691 E1 genes, 98.4 % [97.9, 98.8]; 21 genes detected only in E2. DB: 4,749 records / 1.77 MB (A) vs 135,019 / 14.1 MB (B, + 130,270 ORFs). Library: 743,151 vs 4,088,653 precursors (182 MB vs 1.02 GB). Wall times below | Abundance agreement E1/E2: Spearman 0.996 [0.995, 0.997] over 2,634 shared genes (1,758 TUs); Bland-Altman bias 0.09 log2 |
| 15 | Class q floor | **pass** | E2 warns `novel_orfs.class_q_floor_above_alpha`; classes have 11-62 candidates (< 100); minimum class FDR 0.14 (intergenic, nterm_extension), 0.50 (cds_outframe), 0.67 (antisense) | Alpha 0.01 not loosened; hence 0 calls and no second round |
| 16 | Determinism, TE7 | **pass** | E1 rerun (job 1999366, same inputs, seed 20261014): `tables/` byte-identical, 22 / 22 files (`diff -r`) | `gippa report --redact` on E1 afterwards did not change `tables/` (sha256 checked) |

Counts: **10 pass** (row 7's junction sub-check not testable), **3 fail, 1 not testable, 1 dropped (row 11, E3), 1 reported** (row 14 has no pass threshold).

## Key numbers

| | E1 (route A) | E2 (route B) | E4 (route C, de novo) |
|---|---|---|---|
| DIA-NN precursors at 1 % (per run) | 35,342-40,070 | 34,455-38,634 (one-run probe: 33,762) | Cascadia: 1,502-1,652 `.ssl` rows per run; 6,389 converted, 6,362 ≥ 7 aa |
| DIA-NN protein groups, global q ≤ 0.01 | 2,649 | 2,606 | - |
| GIPPA genes detected (cds + sORF, 4,368) | 2,691 (tier `decoy_picked`) | 2,669 (`decoy_picked`) | 187 (`two_peptide`; 202 single-peptide, unestimated) |
| Gene-level FDR at q ≤ 0.01 (decoy estimate) | 0.0093 | 0.0097 | unestimated (no decoys) |
| Detected with ≥ 2 peptides | 2,455 | 2,403 | 187 |
| Novel-ORF candidates / called | - | 124 / 0 | 4 / 0 |
| Placement FDR | not applicable (`same_db`) | not applicable (`same_db`) | decoy counts: 0 differences 700 targets / 5 decoys accepted; 1 difference 314 / 1 accepted; 2 differences 68 / 0, none accepted |

## Candidate discoveries for Khoa's review (not confirmed)

E2 candidates with ≥ 2 peptides (none called; class FDRs as listed), with the RefSeq feature they
overlap (answer-key annotation):

| ORF | class | coordinates (strand) | peptides | best peptide q | class FDR | RefSeq | Reading |
|---|---|---|---|---|---|---|---|
| NOVEL_U00096.3_257829_p | intergenic | 257,829-257,912 (+) | 2 | 3.3e-7 | 0.14 | crl fragment 1 (pseudogene b0240) | expected (answer key): the crl+ strain expresses the part the IS1 splits off in the reference |
| NOVEL_U00096.3_3035184_m | nterm_extension | ..-3,035,184 (−), 334 aa | 4 | 3.0e-7 | 0.14 | prfB b2891 | expected: prfB upstream of the Bakta start |
| NOVEL_U00096.3_1800096_m | nterm_extension | ..-1,800,096 (−), 155 aa | 2 | 3.6e-7 | 0.14 | infC b1718 (RefSeq 1,800,096-1,800,638; Bakta 1,800,096-1,800,530) | known: Bakta starts infC 36 codons late (infC uses the AUU start); exploratory |
| NOVEL_U00096.3_3560288_m | intergenic | 3,560,288-3,560,605 (−) | 2 | 3.6e-7 | 0.50 | glpR pseudogene b3423 | known pseudogene, translated; exploratory |
| NOVEL_U00096.3_4081770_m | cds_outframe | 4,081,770-4,081,850 (−), 26 aa | 2 | 3.1e-3 | 1.0 | inside fdoI b3892, other frame | unexplained; class FDR 1.0, likely false |
| NOVEL_U00096.3_1547575_m | antisense | 1,547,575-1,548,744 (−), 389 aa | 2 | 2.5e-3 | 0.86 | antisense to fdnG | unexplained; class FDR 0.86, likely false |

E4 `high` ALT calls (row 12): carA aa 340 (30,668 +), pckA aa 539 (3,534,432 +), rpmG aa 22
(3,811,352 −), 2 ALT peptides each. On the strain's own genome these are best read as de novo
sequencing errors or modifications; listed for review, not as variants.

## Route A vs route B

| | Route A (E1) | Route B (E2) | B / A |
|---|---|---|---|
| makedb records | 4,290 CDS + 78 sORF + 381 CONT_ | + 130,270 six-frame ORFs (≥ 20 aa) | 28x records |
| FASTA size | 1,765,715 bytes | 14,084,025 bytes | 8.0x |
| makedb wall time (allocation) | 5 s | 11 s | |
| predicted library | 743,151 precursors, 182,445,389 bytes | 4,088,653 precursors, 1,016,383,009 bytes | 5.5x |
| DIA-NN step 1 (16 CPUs) | 2 min 01 s (job 1998573) | 9 min 41 s (1998574) | 4.8x |
| DIA-NN step 2, 4 runs (16 CPUs) | 13 min 14 s (1998769) | 32 min 30 s (1998857); one-run probe 7 min 37 s (1998796) | 2.5x |
| DIA-NN peak memory (step 1 / step 2) | 5.2 / 7.8 GiB | 7.4 / 14.5 GiB | |
| genes detected by GIPPA | 2,691 | 2,669 | E2 keeps 98.4 % [97.9, 98.8] of E1's |
| `gippa run` wall time (14 threads, 54 GB) | 2 h 32 min (1998814; TE7 rerun 2 h 23 min) | 2 h 16 min (1998879) | |

The six-frame search cost 1.6 % of route A's detected genes and about 2.8x the DIA-NN time in total (2.5x for the search step); GIPPA
time was unchanged.

## Command log (every GIPPA command, pinned commit `a2a95c76`)

All from the detached worktree `local/worktrees/campaign-pin` (tag v1.1.0, `git_dirty` false in
every manifest) with `PYTHONPATH=<pin>/src`; `$D` = `$GIPPA_ROOT/campaign/ecoli_k12`.

```bash
gippa demo --outdir $D/cmds/demo                                   # job 1998577: exit 0, 69 s
gippa example fetch --outdir $D/cmds/example                        # job 1998577: exit 0, 1 s
gippa makedb --annotation $D/bakta/k12_mg1655.json --genome $D/genome.fna.gz --outdir $D/db_a          # allocation, 5 s
gippa makedb --annotation $D/bakta/k12_mg1655.json --genome $D/genome.fna.gz \
  --six-frame-orfs --min-orf-aa 20 --outdir $D/db_b                                                    # allocation, 11 s
gippa convert --from diann $D/diann_a/report.parquet --out $D/in_a/peptides.tsv.gz \
  --abundance-out $D/in_a/abundance.tsv.gz --search-db same                    # job 1998810 (twice, byte-identical)
gippa preflight --genome $D/genome.fna.gz --annotation $D/bakta/k12_mg1655.json \
  --peptides $D/in_a/peptides.tsv.gz --outdir $D/in_a/preflight --strict       # job 1998810: exit 0
gippa convert --from diann $D/diann_b/report.parquet --out $D/in_b/peptides.tsv.gz \
  --abundance-out $D/in_b/abundance.tsv.gz --search-db same                    # job 1998878 (twice, byte-identical)
gippa preflight ... --peptides $D/in_b/peptides.tsv.gz --outdir $D/in_b/preflight --strict   # job 1998878: exit 0
gippa run --profile farmshare,apptainer --genome $D/genome.fna.gz --annotation $D/bakta/k12_mg1655.json \
  --peptides $D/in_a/peptides.tsv.gz --db-dir $GIPPA_DB_DIR --threads 14 --memory 54.GB --outdir $D/runs/out_e1      # 1998814
gippa run ... --peptides $D/in_b/peptides.tsv.gz ... --outdir $D/runs/out_e2                                      # 1998879
gippa run ... --peptides $D/in_a/peptides.tsv.gz ... --outdir $D/runs/out_e1_te7                                  # 1999366 (TE7)
gippa compare $D/runs/out_e1 $D/runs/out_e2 --labels genes,sixframe --outdir $D/runs/cmp_e1_e2   # allocation: exit 0, 9 s
gippa report $D/runs/out_e1 --redact                                           # job 1999542: exit 0, 66 s
gippa report $D/runs/out_e2 --redact                                           # job 1999542: exit 0, 88 s
gippa convert --from cascadia $D/cascadia/*/*.ssl --out $D/in_c/peptides.tsv.gz   # allocation: exit 0, 2 s (twice, byte-identical)
gippa preflight ... --peptides $D/in_c/peptides.tsv.gz --strict                # exit 3: only warning short_peptide (27 rows < 7 aa)
gippa preflight ... --peptides $D/in_c/peptides_min7.tsv.gz --strict           # exit 0 (the 27 rows GIPPA drops anyway removed first)
gippa run ... --peptides $D/in_c/peptides_min7.tsv.gz --reference-peptides $D/in_a/peptides.tsv.gz \
  ... --outdir $D/runs/out_e4                                                  # 2000505
```

`gippa run` was submitted through `validation/campaign/slurm/submit_gippa_run.sh`, which calls the
pinned `scripts/slurm/gippa_run.sbatch` with `GIPPA_REPO=<pin>`: profile `farmshare,apptainer`, all
modules (`modules_off` empty), default seed 20261014, no `--abundance`, no `--strict`. Every run:
116 tasks, 0 failed, exit 0. Not run: the E2 second round (`makedb --patch`; no ORF called) and any
E3 `gippa run` (no knockout run set).

## Deviations and problems

- **Contaminants.** `gippa makedb` appends the 381-protein contaminant library as `CONT_` records by
  default, so DIA-NN got one FASTA, not the pilot's extra `--fasta contaminants.fasta` (which would
  have duplicated them). The DIA-NN flags are otherwise the pilot's `diann_A2.cfg`
  (`validation/campaign/diann_campaign.cfg`, diffed identical).
- **TE5 family IDs.** makedb drops the contaminant descriptions, so the pilot's family rule was
  applied to `contaminants.fasta` and the IDs given makedb's `CONT_` prefix (`phase1_metrics.py`).
- **Route C preflight.** `preflight --strict` exited 3 on the Cascadia table because of the
  `short_peptide` warning alone (27 of 6,389 rows below 7 aa, which GIPPA drops itself); those rows
  were removed before preflight, which then exited 0. The pilot's A3 run had the same warning (62
  rows) and ran without `--strict`.
- **Genome file.** The run genome is `tests/data/genomes/k12_mg1655.fna.gz` (U00096.3, sequence
  identical to GCF_000005845.2), because the reused Bakta JSON names its contig U00096.3. The
  unstripped JSON and its GBFF were copied side by side, so antiSMASH ran (`layers_missing` empty in
  E1/E2).
- **E1 warning** `novel_orfs.decoys_not_makedb` (route A has no ORF decoys): expected.
- **PRIDE downloads** dropped every 0.2-0.9 GB; the resumable fetch (`fetch_runs.sbatch`) finished
  each E1 file in 2-3 range requests (1 h 45 min for 3 files). Panorama served 22 GB in 26 min.
- **E3** used one combined DIA-NN search (8 runs, DB A + 5 UniProt cassette proteins), 28 min.

## Wall times (Slurm)

| step | jobs | elapsed |
|---|---|---|
| E1 fetch (3 runs, 5.6 GB) | 1998566 | 1 h 45 min |
| E3 fetch (8 runs, 22 GB) | 1998567 | 26 min |
| DIA-NN E1 / E2 / E3 | 1998573+1998769 / 1998574+1998857 / 1998659 | 15 min / 42 min / 28 min |
| convert + preflight E1 / E2 | 1998810 / 1998878 | 33 s / 30 s |
| ThermoRawFileParser 1.4.5, 4 runs | 1998581, 1998770 | 3 min + 13 min (544-606 MB mzML each) |
| Cascadia 0.0.7, 1 L40S GPU per run | 1998593, 1998968-70 | 1 h 00-03 min per run |
| `gippa run` E1 / E2 / TE7 / E4 | 1998814 / 1998879 / 1999366 / 2000505 | 2 h 32 / 2 h 16 / 2 h 23 / 1 h 51 |

## Provenance

| item | value |
|---|---|
| GIPPA | `a2a95c76b41665464002d01a4d52c11ccb7e38db` (v1.1.0), clean |
| DIA-NN | 2.3.1 Academia, image sha256 `1946452d56d2fb0f66e05ced473df88e53484e28040f1ecbd83c4a43405a6381` |
| Cascadia | 0.0.7, weights sha256 per `validation/pilot/a3_assets.tsv` (`scripts/slurm/denovo_pilot.sbatch`) |
| ThermoRawFileParser | 1.4.5 (bioconda), lock `envs/locks/thermorawfileparser.explicit.txt` |
| Bakta | 1.12.2, DB 6.0 full, image `sha256:1e24ae2f…`, recording job 1849096 (reused) |
| genome | `k12_mg1655.fna.gz` sha256 `4bfd7163f6f0de6394a2131cbfdbe75100daa55cbee227ec3e403a44640ac03c` |
| Bakta JSON / GBFF | `ae0ccaf186b0…` / `c958d300ab43…` |
| E1 raw | 4 PXD058422 runs, SHA-1 verified against the deposit; sha256 in `provenance/E1.json` |
| E3 raw | 8 PXD066241 runs (no published checksum; TLS + size); sha256 in `$D/raw_e3/SHA256SUMS` |
| peptides E1 / E2 / E4 | `50dbd6b32c37…` / `dfb8e80e2580…` / `8c697dbf3a14…` (full sha256 in `provenance/`) |
| container digests | in `provenance/E*.json` › `gippa_run.tools` |
