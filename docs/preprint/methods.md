# Methods

<!-- Draft Methods section for the GIPPA 1.1.0 paper. Every paragraph ends with an HTML comment
naming its sources (repo-relative paths; `defaults.yaml` = src/gippa/config/defaults.yaml).
Where a methods note, ADR or the spec disagrees with the code, this text follows the code of
release 1.1.0 (release code e43683c7; the v1.1.0 tag is on a2a95c76, which adds only the readiness record). Results are reported in the Results section. -->

## 1. Overview and design

GIPPA (genome integrated prokaryotic proteomics analysis) takes one assembled bacterial genome and proteomics evidence from any search engine (identified peptides, protein-level abundances, or de novo peptides), and returns one gene table keyed by Bakta locus tags. Around that table it runs genome-context statistics: detection and abundance of annotated regions, autocorrelation along the replicon, abundance islands, the origin-to-terminus gradient, operon concordance, functional enrichment and sequence determinants of abundance. It also calls amino-acid variants and novel open reading frames (ORFs) from peptides that do not fit the annotation. Each inferential claim carries an error estimate (a q-value, FDR, FWER, confidence interval or an explicit "unestimated" flag). GIPPA does not search spectra. It inherits the identification error of the search engine and estimates two further error sources itself: the error of placing a peptide on the genome, and the error of calling a gene detected.

<!-- src: CITATION.cff:3-11; GIPPA_SPEC.md:2001-2005; docs/adr/ADR-0071-r2-mapping-fdr-waiver-1-0-0.md (Consequences); docs/preprint/outline.md §2 item 3 -->

The analysis has six stages: (i) preflight checks of inputs and environment; (ii) annotation layers (Bakta genes; mobile elements, biosynthetic gene clusters, CRISPR arrays and resistance genes; operons, promoters, terminators and ribosome binding sites); (iii) peptide- and protein-to-genome mapping with a placement FDR; (iv) the gene table with gene-level detection FDR, abundance and a cross-fitted detectability model; (v) the context, operon, variant, novel-ORF, function and determinant modules; and (vi) tables, plots and a self-contained HTML report. Version 1 analyses one sample per run. It has no tests across samples. `gippa compare` describes agreement between runs on the same genome.

<!-- src: src/gippa/cli.py:256,484,588,1316,1379,1442,1576; workflow/subworkflows/local/ (12 subworkflows); docs/preprint/outline.md §6 "Single sample per run" -->

Default significance levels are set in one configuration block (Table M1). Permutation and bootstrap counts (`stats.permutations`, `stats.bootstrap`) default to 1,000; `--publication` raises these two and the determinants null to 10,000, while the de novo precision bootstrap (1,000), fgsea's `nPermSimple` (1,000) and the sensitivity reruns (B = 199) stay fixed. Defaults below are those of `src/gippa/config/defaults.yaml` in release 1.1.0. The global seed is 20261014. Mathematical symbols are local to the section that defines them.

**Table M1.** Significance levels in use (`stats.alpha`; a further key, `determinants` = 0.05, is not read by the 1.1.0 code).

| Key | Default | Used for |
|---|---|---|
| `identification` | 0.01 | gene-level detection q |
| `placement` | 0.01 | peptide placement q |
| `novel_orfs` | 0.01 | class FDR of novel ORFs |
| `variants` | 0.01 | placement q of variant-supporting peptides |
| `context` | 0.05 | region tests, autocorrelation, islands |
| `enrichment` | 0.05 | ORA/GSEA significance and GO grouping |
| `operon_flags` | 0.10 | operon discordance and missed operons |

<!-- src: src/gippa/config/defaults.yaml:2,94-103; src/gippa/config/__init__.py:18,81-92; src/gippa/config/__init__.py:83-93; src/gippa/steps/mapping.py:441; src/gippa/context/sensitivity.py:47-48 -->

## 2. Inputs, preflight and genetic codes

**Genome.** The genome is a FASTA file (gzip accepted). The first header token is the contig ID. A header with `circular=true` or `[topology=circular]` marks a circular replicon. Characters outside `[A-Za-z0-9._-]` in contig IDs are replaced by `_` (with a warning); IDs that collide after mapping, and non-IUPAC characters such as `U` and `-`, are errors. Internally all intervals are 0-based and half-open; a feature crossing the origin of a circular replicon is stored unwrapped (end > length), and `crosses_origin` is derived from it. Writers convert once: TSV, GFF3 and VCF are 1-based closed, BED is 0-based half-open.

<!-- src: GIPPA_SPEC.md:767-772,938-942; docs/adr/ADR-0005-m1-input-decisions.md item 1; src/gippa/model/coords.py:1-3,37; src/gippa/preflight/genome.py:213-214 -->

**Peptides.** The generic peptide table has the columns `peptide, modified_peptide, protein_ids, q_value, pep, score, is_decoy, intensity, run, charge, spectrum_id, source, search_db`; only `peptide` is required. Rows with q or PEP outside [0, 1], negative or non-finite intensity, or charge < 1 are dropped and reported. Peptides that are not upper-case A–Z, or are shorter than 7 or longer than 40 residues, are dropped with a warning. Rows above the input q cut-off (`input.q_value_max` = 0.01) are kept in the table, because they inform FDR estimation, but are never evidence. Decoy status defaults from protein-ID prefixes (`rev_`, `REV_`, `DECOY_`, `decoy_`, `XXX_`, `##`); a row is a contaminant when all its protein IDs carry a contaminant prefix (`Cont_`, `CON__`, `CONT_`, `cRAP`). Rows are then aggregated to one row per stripped peptide: minimum q and PEP, maximum score, decoy only if all rows are decoy, the union of protein IDs, and an intensity equal to the per-run sum over charges and modified forms, then the median across runs. Rows without a `source` value take a table default: database search output when any row has a q or PEP, otherwise de novo output. `search_db` is resolved for the whole table: `same` when at least 95 % of target peptides with engine protein IDs name a Bakta locus tag of the run, otherwise `other`. A per-row `source` or `search_db` value overrides the table default. One `(run, spectrum_id)` used with different peptides is an error (rows dropped); repeated `(run, spectrum_id, peptide, charge)` rows keep the first; a peptide with both database and de novo rows counts as database.

<!-- src: src/gippa/io/peptides.py:32-46,181-255,305-312,354-405; src/gippa/steps/mapping.py:95-108,151,165-176; defaults.yaml:9-14; src/gippa/io/peptides.py:223-301; src/gippa/steps/mapping.py:94-106,176-199; src/gippa/io/peptides.py:309 -->

**Proteins and abundances.** Protein FASTA records are routed by length: records shorter than 7 residues are dropped, 7–40 residues follow the peptide path, and longer records follow the protein path (Section 5.3). UniProt `sp|ACC|NAME` IDs reduce to the accession. A protein abundance table has one ID column and one value column. Its scale is detected as log2 when any value is negative or all values are ≤ 50, else linear; the user can override it. More than 20 % missing values raises a warning. The declared abundance type (iBAQ, Top3, MaxLFQ, intensity or unknown; default unknown) decides whether protein-route values count as molar (Section 7).

<!-- src: src/gippa/io/protein_fasta.py:21-22,69-77; GIPPA_SPEC.md:774-778,812; src/gippa/io/abundance.py:42-43,141-148,195; defaults.yaml:25,81 -->

**Converters.** `gippa convert` turns DIA-NN 2.x reports (`report.parquet` or TSV), Cascadia 0.0.7 `.ssl` files and generic tables into the peptide table. For DIA-NN the stripped sequence, protein IDs, global q-value (else run q-value), PEP, decoy status, `Precursor.Normalised` (else `Precursor.Quantity`; zero read as missing), run and charge are carried over; an optional protein table holds, per protein group, the median over runs of `PG.MaxLFQ`. For Cascadia, the score is carried over without a q-value or PEP, because the score is not calibrated; Cascadia's own score cut-off is counted, not applied. Modifications are stripped by removing bracketed groups, flanking residues and mass deltas. Converters for Casanovo and InstaNovo are listed in the command-line interface but are not implemented in 1.1.0 and exit with code 5 and the message "not implemented".

<!-- src: src/gippa/io/converters/__init__.py:20-24; src/gippa/io/converters/diann.py:1-12,37-47,133-154; src/gippa/io/converters/cascadia.py:1-31; src/gippa/io/converters/strip.py:1-30; src/gippa/cli.py:1366-1399 -->

**Preflight.** `gippa preflight` (also the first step of `gippa run`) validates the genome, proteins, peptides, abundances, configuration and environment, and writes `preflight.json` and `preflight.md`. Errors stop the run; warnings stop it only with `--strict` (exit code 3). Genome warnings fire at more than 10 contigs, total length below 100 kb, N50 below 100 kb or a largest contig below 1 Mb, GC outside 0.25–0.75, and runs of at least 10 N, and runs of any other IUPAC ambiguity code. Exit codes distinguish invalid input (1), missing tools or databases (2), warnings under `--strict` (3), a failed step (4) and internal errors (5).

<!-- src: src/gippa/preflight/*.py (rule IDs); defaults.yaml:17-25; src/gippa/exitcodes.py:8-14; GIPPA_SPEC.md:976-1022; src/gippa/preflight/genome.py:161-174,257-275; src/gippa/preflight/genome.py:71,87,271 -->

**Genetic codes.** `annotation.translation_table` accepts 11 (default), 4 or 25; other values are refused. The three tables differ only at TGA (stop, Trp, Gly). One table applies to the whole run and is passed to every tool: Bakta receives `--translation-table` when the table is not 11, miniprot uses it at indexing, and tblastn's `-db_gencode` is derived from it. The six-frame index records its table, and an index or Bakta JSON built with another table is refused. A preflight check measures the genome length (bp, union over six frames) covered by stop-free open frames of at least 300 codons (no start codon required; records treated as linear), with TGA read as stop ($b_{11}$) and as sense ($b_4$). It warns when $b_{11}/b_4 < 0.55$ at table 11, or $b_{11}/b_4 \ge 0.55$ at table 4 or 25, and only when at least 20 such ORFs exist with TGA read as sense. The warning never changes results. Codes 4 and 25 cannot be told apart from the genome; the user chooses from taxonomy.

<!-- src: src/gippa/model/genetic_code.py:24; defaults.yaml:29; docs/adr/ADR-0067-genetic-codes-4-11-25.md decisions 1,6; src/gippa/annotate/bakta_cmd.py:58; src/gippa/mapping/miniprot.py:29-31; src/gippa/mapping/tblastn.py:100-116; src/gippa/mapping/sixframe.py:342-352; src/gippa/preflight/genetic_code.py:33-37; docs/methods/genetic-code.md:5-21; src/gippa/preflight/genetic_code.py:117-148; src/gippa/preflight/genetic_code.py:11,92-121 -->

## 3. Annotation layers

**Genes.** The gene model comes from Bakta 1.12.2 with database 6.0 [Schwengers et al. 2021]. A Bakta JSON supplied by the user is used as is. A GFF3 or GenBank file from another annotator is passed to Bakta with `--regions`, so that its CDS are kept and Bakta adds the rest. Bakta runs with `--complete` when every contig is circular. With `db: auto`, the full database is used when installed, else the light one; `gippa setup` installs the light database by default, and the pilot and the weekly real-genome runs used the full database.

<!-- src: src/gippa/annotate/tools.py:91-93; src/gippa/data/db_sources.yaml:6-21; src/gippa/annotate/bakta_cmd.py:18-33,54-61; defaults.yaml:28; src/gippa/cli.py:183; GIPPA_SPEC.md:814-817,1031 -->

**Regions.** Prophages come from geNOMAD 1.12.0 [Camargo et al. 2024] (database 1.9; `end-to-end --cleanup --splits 8`). Plasmids are replicons that Bakta marks as plasmids or that geNOMAD calls plasmid. Biosynthetic gene clusters (BGCs) come from antiSMASH 8.0.4 [Blin et al. 2025], run with `--genefinding-tool none --minimal` on the Bakta genes. CRISPR regions are Bakta's CRISPR arrays plus `cas` genes within 1,000 bp. Antimicrobial resistance (AMR) regions are gene sets, one per contig, formed from the AMRFinderPlus [Feldgarden et al. 2021] hits that Bakta reports. A gene belongs to a region when at least 50 % of its bases lie inside it. Prophage genes receive a descriptive function subset (integrase, repressor, structural, lysis, cargo) from keyword rules.

<!-- src: src/gippa/annotate/regions.py:1-8,45,124,147-228; src/gippa/annotate/genomad_cmd.py:13-16; src/gippa/annotate/antismash_cmd.py:16-32; src/gippa/annotate/tools.py:113-114; defaults.yaml:31 -->

**Operons and transcription units.** UniOP (git commit 3ed14dfc) [Su et al. 2024] gives a probability for each pair it scores. Other adjacent same-strand pairs, or all pairs when UniOP is unavailable, receive a distance model $p = 1/\{1 + \exp[(d - d_0)/s]\}$, where $d$ is the intergenic distance in bp. Transcription units (TUs) are maximal runs of same-strand genes linked by $p \ge 0.5$. The parameters $d_0 = 102.687$ bp and $s = 47.160$ bp were fitted on RegulonDB 14.5.0 [Salgado et al. 2024] transcription units of *E. coli* K-12 MG1655 (U00096.3), so that $p \ge 0.5$ when $d \le 102$ bp. The threshold 0.5 is a fixed value chosen by the project owner; the rule "smallest threshold with precision ≥ 0.90" gave 0.529. On the K-12 labelled pairs, UniOP precision at 0.5 is 0.896 and the distance model's 0.828, below the 0.90 target. When UniOP is unavailable, the run carries the flag `operon_source_distance_model` (with a preflight warning), and TU-block Moran's I rows and both-detected discordance rows carry `tu_source_distance_model`; individual pairs that a working UniOP does not score use the distance model without a flag. TUs are the resampling and grouping unit of most statistics below.

<!-- src: src/gippa/annotate/operons.py:1-7,65-71,145-153; src/gippa/data/calibration.yaml (uniop_pair_threshold, operon_distance); defaults.yaml:32; docs/adr/ADR-0022-uniop-terminatornet-cutoffs.md:38-40,256-273; docs/adr/ADR-0061-distance-model-fallback-scope.md decision 5; src/gippa/operons/flags.py:254-269; src/gippa/annotate/operons.py:24,83-95; src/gippa/operons/flags.py:59,259-275 -->

**Promoters, terminators and ribosome binding sites.** Promotech (RF-HOT model, git commit 56251ad9) [Chevez-Guardado and Peña-Castillo 2021] scores 40-nt windows whose 3′ ends lie from −200 to +20 nt around each start codon; a "promoter" call needs a score ≥ 0.5. Promoter Calculator 1.2.4 [LaFleur et al. 2022] supplies predicted transcription rates, used as features but never as a promoter call. TerminatorNet (git commit 72180309) [Tjaden 2026] is the primary terminator predictor; a call needs confidence ≥ 0.9999999999961204, the maximum-F1 cut-off on RegulonDB 14.5.0 K-12 terminators. Its precision there (0.226) is a lower bound, because RegulonDB is incomplete; the cut-off carries no precision guarantee. TransTermHP 2.09 [Kingsford et al. 2007] runs as a cross-check that never enters a call. When TransTermHP is chosen as the primary tool instead, its confidence (0–100) is divided by 100 and the same cut-off applies, so only confidence-100 hits count as calls. BacTermFinder (git commit 02487044) [Taheri Ghahfarokhi and Peña-Castillo 2025] is off by default. OSTIR 1.1.3 [Roots et al. 2021], with ViennaRNA 2.6.4 [Lorenz et al. 2011], computes translation initiation rates (TIR) from −50 to +50 nt around starts, with the anti-Shine–Dalgarno sequence taken from the 3′ end of the strain's 16S rRNA (default `ACCTCCTTA`, flag `anti_sd_default`, when no 16S is annotated or its most common 3′ 9-mer lacks CCTCC). Because these layers predict σ70 promoters and intrinsic terminators only, a missing call is never read as absence, and a missing layer raises a `layer_missing:<layer>` flag rather than a "no".

<!-- src: defaults.yaml:34-45; src/gippa/annotate/promoters.py:1-10,240-256; src/gippa/operons/pairs.py:175-195; src/gippa/operons/flags.py:9-12,89-104; src/gippa/data/calibration.yaml (terminatornet_yes); docs/adr/ADR-0022-uniop-terminatornet-cutoffs.md:41,50-52; src/gippa/steps/terminators.py:1-7; workflow/envs/ostir.yml:1,10; src/gippa/annotate/bakta_layers.py:20-21,40-52; workflow/envs/ostir.yml:1,8; src/gippa/annotate/terminators.py:145; src/gippa/operons/flags.py:96-99; src/gippa/novel/evidence.py:335 -->

**Codon usage.** CAI, tAI, ENC, Fop, GC3 and codon-pair bias are computed in-process with codon-bias 0.5.0 on codons 2 to n−1, in the run's genetic code, plus the AT content and the number of rare codons (genome-wide synonymous share < 0.10) in codons 2–10. CAI and Fop use ribosomal-protein genes as the reference set [Sharp and Li 1987]; tAI uses tRNA gene copy numbers [dos Reis et al. 2004]; ENC follows [Wright 1990].

<!-- src: src/gippa/annotate/codon_usage.py:1-10; src/gippa/annotate/tools.py:117-125; defaults.yaml:46; src/gippa/annotate/codon_usage.py:35,93-111 -->

## 4. Strain search database (`gippa makedb`)

From a Bakta JSON (`--annotation`; makedb does not run Bakta), `gippa makedb` writes an optional strain database for the user's own spectrum search. It contains Bakta CDS and sORF proteins, optional patched proteins from a previous run, optional six-frame ORFs and, unless switched off, the universal contaminant library. It contains no decoys; decoys come from the search engine. Six-frame ORFs are maximal stop-to-stop segments per frame (split at ambiguous codons) of at least 20 codons (`--min-orf-aa`), translated in the run's code. ORFs in frame with an annotated gene are dropped; each such gene with at least 7 stop-free upstream codons instead gets one N-terminal-extension entry (upstream segment plus the first 30 codons of the gene). The other ORFs keep their locus class (intergenic, out-of-frame or antisense) in a class-tagged ID (`GIPPAORF_<class>_…`). These tags let the novel-ORF module pair engine decoys with targets (Section 13). GIPPA builds no FDR on the database itself.

<!-- src: src/gippa/makedb/build.py:1-15,45-47; src/gippa/makedb/__init__.py:4; src/gippa/cli.py:1312,1329-1351; docs/methods/target-decoy.md:120-127,168-175; src/gippa/novel/fdr.py:53; src/gippa/cli.py:1306,1331-1333 -->

## 5. Peptide-to-genome mapping and placement FDR

### 5.1 Exact matching

All six frames of every replicon are translated in the run's genetic code. Ambiguous codons become X, which never matches, and I is folded to L in both the index and the queries. For circular replicons the first 126 nt ($3 \times (40 + 2)$) are appended, so that peptides of up to 40 residues, plus a two-difference window, can span the origin; hits are mapped back modulo the length. An Aho–Corasick automaton [Aho and Corasick 1975] of all distinct query peptides of 7–40 residues (non-decoy, q ≤ 0.01 or no q) scans every target frame and the decoy frames of the peptide's own source (Section 5.4), and every locus is kept. An initiator methionine at peptide position 0 is credited at the run's own annotated non-ATG CDS or sORF starts (GTG, TTG) on this route only.

<!-- src: src/gippa/mapping/sixframe.py:33-42,244-245; src/gippa/mapping/exact.py:1-12; src/gippa/steps/mapping.py:111-131; defaults.yaml:30,58; docs/adr/ADR-0068-initiator-methionine-non-atg-starts.md -->

### 5.2 Tolerant matching

Peptides of 10–40 residues with no exact target hit are searched with a difference budget $k$: one difference at 10–11 residues and two at 12–40. Each of the following costs one difference: a substitution; a swap of adjacent residues; a mass-equivalent block exchange ({N}↔{G,G}, {Q}↔{A,G}, {G,T}↔{A,S}, {A,V}↔{G,L}, {A,N}↔{G,Q}); a near-equal exchange ({W}↔{A,D}, {W}↔{E,G}, {R}↔{G,V}); and a ±1-nt frameshift split with at least 4 residues on each side. I/L is free. Insertions and deletions are not modelled. The search is a lossless q-gram filter, with the largest $q \in \{5, 4, 3\}$ that gives a positive threshold $t = (L - q + 1) - k(q + 1)$, followed by a banded dynamic-programming verifier (band ±2). A target stop may be consumed only at the annotated stop codon of a CDS or sORF, as a substitution; peptides shorter than 10 residues that span a stop are never placed. The decoy search runs for every searched peptide. Under the default `genome_relation = related`, accepted tolerant hits count as evidence; under `same`, only exact hits do.

<!-- src: src/gippa/mapping/tolerant.py:1-17; defaults.yaml:6,65-71; src/gippa/mapping/diffops.py:30-41; src/gippa/mapping/qgram.py:45-46; src/gippa/mapping/tolerant_dp.py:46; docs/methods/placement-fdr.md:835-847; docs/adr/ADR-0058-annotated-stop-substitution.md; src/gippa/genetable/evidence.py:142-163 -->

### 5.3 Proteins longer than 40 residues

Long proteins are aligned with miniprot 0.18 [Li 2023] (`-S -N 5 --outs=0.95`, genetic code at indexing). An alignment is assigned to a same-strand Bakta CDS or sORF when identity, query coverage and reciprocal overlap are all ≥ 0.80; the largest overlap wins. Proteins left unassigned are searched with tblastn from BLAST+ 2.16.0 [Camacho et al. 2009] (`-seg no -evalue 1e-5`), collinear HSPs are chained (genome gap ≤ 300 nt, query overlap ≤ 10 residues) and the same thresholds apply. Assignments are labelled identical (identity and coverage 1), near (identity ≥ 0.95) or divergent. Frameshifts and internal stops are flags only. This route has no decoys and carries no placement FDR.

<!-- src: src/gippa/mapping/miniprot.py:29-36; src/gippa/mapping/protein_assign.py:1-9,28,72-78; src/gippa/mapping/tblastn.py:1-28,100-133; defaults.yaml:59-64; src/gippa/steps/mapping.py:3-4 -->

### 5.4 Genome decoys and the placement FDR

The placement FDR is the expected share of accepted peptide placements that are chance hits of a wrong peptide. It does not estimate the identification error of the input, which is inherited from the search engine through the q ≤ 0.01 filter. Decoy genomes are built per input source. For database search output, every translated frame is reversed; for de novo output, residues are shuffled (Fisher–Yates, seeded) within each stop-free segment in windows of 50 residues. Shuffling is used for de novo input because de novo errors often reverse a stretch of sequence, which would inflate hits on reversed decoys. Both decoys keep frame composition, stop and X geometry, and the number of candidate positions. Decoy hits are counted, never placed.

<!-- src: docs/methods/placement-fdr.md:28-105; defaults.yaml:72; src/gippa/mapping/decoys.py:26-89; GIPPA_SPEC.md:2001-2005 -->

Each peptide is ranked by a chance-hit score $E_i = N \cdot \Pr(\le m \text{ mismatches})$, where $N$ is the number of six-frame codon positions, $m$ is the number of differences of its best hit, and the probability is a Poisson-binomial with per-position mismatch probability $1 - f_{a_j}$ from the genome's six-frame residue frequencies. $E_i$ is a ranking score, not an FDR estimate. Each peptide counts once at its best target cost $m_T$ and, separately, at its best decoy cost $m_D$. Lists are formed per decoy type and stratum $m \in \{0, 1, 2\}$, with length classes pooled. Within a list,

$$\widehat{\mathrm{FDR}}_m(e) = \frac{W_m(e) + c}{T_m(e)}, \qquad W_m(e) = \sum_{j \in \mathcal{D}_m:\, E_j \le e} \hat\pi_j ,$$

where $T_m(e)$ is the number of target peptides with $E \le e$, $\mathcal{D}_m$ the decoy hits, $\hat\pi_j$ the estimated probability that peptide $j$ is wrong, and $c$ a +1-type correction: $c = 1$ at $m \ge 1$; at $m = 0$, $c$ is the largest class estimate $\hat\pi_\ell$ (database lists) or the largest $\hat\pi$ of the searched peptides (de novo lists). Ties count against the target. The q-value is the running minimum over $e \ge E_i$, clipped to [0, 1], and a placement is accepted at q ≤ 0.01. With constant $\hat\pi$ the estimator reduces to $\hat\pi(D + 1)/T$ [Barber and Candès 2015; He et al. 2015]. Because q ≥ c/T, a list with c = 1 needs at least 100 targets before its first decoy to accept anything at 0.01.

<!-- src: src/gippa/mapping/chance.py:1-6,52-74; docs/methods/placement-fdr.md:221-233,535-570; src/gippa/mapping/placement_fdr.py:1-7,131-170; docs/adr/ADR-0010-placement-fdr-weighted-decoys.md; defaults.yaml:97; docs/methods/target-decoy.md:28; src/gippa/mapping/placement_fdr.py:494-515; src/gippa/mapping/denovo_precision.py:169,197-205 -->

The weight $\hat\pi$ depends on source and stratum. For database peptides at $m = 0$, it is estimated per length class $\ell$ ([7], [8], [9], [10–11], [12–14], [15–19], [20–40]) as

$$\hat\pi_\ell = \min\!\left\{1, \frac{U_\ell}{1 - r_\ell}\right\},$$

where $U_\ell$ is the one-sided 95 % Clopper–Pearson [Clopper and Pearson 1934] upper bound of the share of class-$\ell$ peptides with no exact target hit, and $r_\ell$ is the share with an exact decoy hit ($\hat\pi_\ell = 1$ if $r_\ell = 1$). Engine PEPs and q-values are not used here. They are read only for a warning (`input_prefiltered`) when fewer peptides miss the genome than PEP predicts. For de novo peptides at $m = 0$, $\hat\pi$ is one minus the lower end (2.5 % quantile) of the 95 % percentile-bootstrap interval of the de novo precision of the peptide's score bin (Section 13), and 1 when that bin is not estimable. At $m \ge 1$, $\hat\pi = 1$. Database peptides from a same-strain search (`search_db = same`) with an exact target hit receive placement q = 0 (flag `same_db`) and leave the stratum-0 lists; for such input GIPPA adds no estimate of identification error. Contaminants are excluded from the stratum-0 $\hat\pi$ population. Peptides with at most 3 distinct residues are flagged `low_complexity`, without changing their q.

<!-- src: src/gippa/mapping/placement_fdr.py:9-13,54-57,177-310,487-503; defaults.yaml:73; docs/adr/ADR-0019-stratum0-pi-hat-exact-share.md; docs/adr/ADR-0011-low-complexity-placement-flag.md; docs/methods/denovo-precision.md ("π̂ for M3-10"); docs/adr/ADR-0071-r2-mapping-fdr-waiver-1-0-0.md -->

*Calibration status.* List-level simulation suites test mean FDP ≤ q + 2 SE at q ∈ {0.01, 0.05, 0.10}. On a 100-kb simulated genome, the genome-level checks are smoke tests (no false peptides in most lists). Evidence for strata $m \ge 1$ rests on the list-level suites, and stratum-0 evidence rests on a stressed 1-Mb arm and checks of the estimate. The contaminant attack `flat_c200` exceeds its bound at q = 0.01 and is kept as a strict expected failure (contaminant weight undecided), and gate M3 was recorded as not met and shipped as a disclosed limit. Neither decoy is a null for wrong peptides that share sequence order with the genome's proteins, so chance hits of protein-like wrong peptides are under-counted.

<!-- src: docs/methods/placement-fdr.md:172-180,700-733,777-818; docs/preprint/outline.md C11; PROGRESS.md:155-165; tests/calibration/attacks (test_attack_contaminants_stays_at_level[flat_c200-full]) -->

### 5.5 Contaminants and unmapped sequences

Peptides that place nowhere are matched exactly (I = L) against the universal contaminant library [Frankenfield et al. 2022] (381 entries, fetched by `gippa setup` at a pinned commit with a SHA-256 check). Unassigned proteins are matched with blastp at ≥ 90 % identity over ≥ 50 % of the query. Remaining rows whose engine IDs all carry a contaminant prefix are labelled contaminants. Every other unplaced sequence is written to `unmapped.tsv` with a reason (`no_hit`, `rejected_placement_fdr`, `below_identity`, `too_short`, `contaminant_id`). Contaminant peptides that do place on a gene and pass the evidence rules count toward detection, with a `contaminant_evidence` flag.

<!-- src: src/gippa/data/db_sources.yaml:95-101; src/gippa/mapping/contaminants.py:1-7,15-16,33-35,174-182; src/gippa/genetable/evidence.py:8-10 -->

## 6. Gene table, detection tiers and gene-level FDR

**Locus classes and evidence.** Each placement takes the highest-priority locus class: `rna_overlap` > `pseudogene` > `sorf` > `cds_inframe` > `nterm_extension` > `cterm_extension` > `cds_outframe` > `antisense` > `intergenic`. Uniqueness is `unique` (one gene), `multi_identical` (several genes with identical protein sequence) or `shared` (different genes). A peptide is evidence for a gene when all of the following hold: input q ≤ 0.01 (or no q); placement q ≤ 0.01 or `same_db`; a gene-assigning class (`cds_inframe`, `sorf`, `nterm_extension`, `cterm_extension`, or `pseudogene`); uniqueness `unique` or `multi_identical`; and the exact route, unless `genome_relation = related`. Shared peptides count toward neither detection nor abundance. Peptides are counted as distinct stripped sequences.

<!-- src: src/gippa/mapping/loci.py:7-17; src/gippa/genetable/contract.py:95-104; src/gippa/genetable/evidence.py:28,142-150,228-232; docs/adr/ADR-0013 (cterm_extension) -->

**Tiers.** The detection tier is chosen once per run. `decoy_picked` applies when decoy rows carry protein IDs, at least 80 % of decoys pair to a gene (the larger of the protein identity and coverage cut-offs, 0.80 by default), and targets and decoys share a ranking column. Otherwise `mayu` applies when accepted targets carry q or PEP (a pairing share below that cut-off always falls back here, with the warning `decoy_pairing_low`). `two_peptide` applies to peptides without q or PEP, and `protein_list` to protein-only input. Partial genes are excluded (`partial_excluded`). A gene is detected at gene q ≤ 0.01 in the two estimated tiers. The number of non-partial genes with at least two distinct accepted peptides (`detected_2pep`) is reported in every tier.

<!-- src: src/gippa/genetable/detection.py:52,222-263,362-394,423-438; defaults.yaml:96; src/gippa/genetable/evidence.py:302; src/gippa/genetable/detection.py:240-242,366-373,417-419 -->

**Picked gene FDR.** In the `decoy_picked` tier, the picked protein FDR [Savitski et al. 2015] is adapted to genes. A decoy pairs with a gene when every one of its prefix-stripped protein IDs maps to that gene; decoys whose IDs are all contaminants, and decoys of `makedb` ORF entries that map to no gene, get no row. Decoys above the input q cut-off leave the competition when the targets carry q-values. For each gene, the best accepted target competes with the best paired decoy on a lexicographic key of the columns both sides carry (q ascending, PEP ascending, score descending); ties go to the decoy. With $T(t)$ target-won and $D(t)$ decoy-won genes at threshold $t$ and $f$ the paired share of decoys,

$$\widehat{\mathrm{FDR}}(t) = \frac{D(t)/f + 1}{T(t)},$$

made monotone into gene q-values; decoy winners get no q. The q floor $1/T$, with $T$ the number of target-won genes, is reported, with a warning when it exceeds α. The $1/f$ factor assumes that decoy loss in pairing is independent of score.

<!-- src: src/gippa/genetable/detection.py:48-50,149-298,367-372,392-394; src/gippa/stats/tda.py:117-143; docs/adr/ADR-0009-picked-fdr-paired-share.md; docs/adr/ADR-0047; docs/methods/detection-tiers.md:84-89 -->

**MAYU gene FDR.** In the `mayu` tier, a closed form of MAYU [Reiter et al. 2009] is used. False peptides land on gene $g$ with probability proportional to $a_g$, its number of observable peptides (floor 1). Genes fall into 10 equal-count size bins. For bin $b$ with $N_b$ genes, $n_b$ of them detected at threshold $t$, and mean weight $w_b = \bar a_b / \sum_g a_g$:

$$p_b = 1 - (1 - w_b)^{F(t)}, \qquad \mathrm{FP}_b = \frac{(N_b - n_b)\,p_b}{1 - p_b}, \qquad \widehat{\mathrm{FDR}}(t) = \frac{\sum_b \mathrm{FP}_b}{\sum_b n_b}.$$

$F(t)$, the expected number of false accepted peptides, is the sum of the PEPs of the distinct accepted peptides with q ≤ t (when all carry a PEP), else $t \cdot |P(t)|$, with $P(t)$ the set of those peptides. With PEP-only input, peptide q is the running mean of sorted PEPs. MAYU is an expectation with no +1 term and no finite-sample guarantee; it is anti-conservative when wrong hits spread more evenly than $a_g$ or when PEPs are optimistic.

<!-- src: src/gippa/stats/mayu.py:1-7,63-160; src/gippa/genetable/detection.py:313-323,400-412; defaults.yaml:77; docs/methods/mayu-gene-fdr.md:43-47,64-121 -->

**Unestimated tiers.** In the `two_peptide` tier, a gene is detected with at least two distinct accepted peptides. Genes with exactly one are `single_peptide_unestimated` and counted as detected only in a sensitivity rerun. In 1.1.0 this threshold is fixed at 2 in the gene table; the configuration key `detection.min_peptides_sequences_only` is read only by the detection-limit model. In the `protein_list` tier, a gene is detected when a protein is assigned to it. Neither tier has a gene q.

<!-- src: src/gippa/genetable/detection.py:430-437; src/gippa/genetable/evidence.py:258; src/gippa/steps/gene_table.py:94; src/gippa/operons/detection_limit.py:392 -->

*Calibration status.* Picked and MAYU tiers pass their simulation suites at q = 0.01, including random loss of 10–20 % of decoys in the picked tier and 30 % loss, which usually falls back to MAYU (the claim holds for whichever tier is chosen). In the simulator, the levels 0.05 and 0.10 are weak or vacuous: simulated peptides are capped at q ≤ 0.01, so the loosest reachable gene FDR is about 0.053; and at q = 0.05 with half of the target q-values missing, the FDP is bimodal with a mean slightly above q. The simulator places false hits on genes in proportion to their observable peptides, which is MAYU's own model, so the unattacked MAYU suite partly assumes what it tests. Loss of the best-scoring decoys is outside the claim (losing the best 10 % or 20 % gives an FDP of about 1.5× or 2× q).

<!-- src: docs/methods/detection-tiers.md:84-89,150-165; docs/preprint/outline.md C12; docs/methods/detection-tiers.md:84-87,150-160; docs/methods/detection-tiers.md:207-213,296-301; tests/calibration/attacks/test_detection_tiers_attack_decoy_loss.py:12-16 -->

## 7. Abundance

Abundance estimates the relative molar amount per gene in one sample, up to a run constant. GIPPA performs no normalisation across runs or samples. Peptide intensities are used when any abundance-eligible accepted peptide has one; otherwise protein-table values are used; the routes are never mixed. Abundance-eligible peptides are `unique` and not in an N- or C-terminal extension, with finite positive intensity. The default summary is iBAQ [Schwanhäusser et al. 2011], $\mathrm{iBAQ}_g = \sum_j I_j / \max(1, n^{\mathrm{obs}}_g)$, where $n^{\mathrm{obs}}_g$ is the number of distinct in-silico tryptic peptides of 7–30 residues with no missed cleavage, cleaving after K/R including before P. Top3 [Silva et al. 2006] (mean of the up to three largest intensities) is the alternative. Protein-route values join through the protein-to-gene map; they are molar only when declared iBAQ or Top3, otherwise relative; on a relative scale the operon discordance and missed-operon tests report `unestimated` with flag `relative_scale`, while the other context and enrichment statistics still use the values. With two proteins on one gene, the larger log2 value is kept (`two_proteins_one_gene`). Linear values are log2-transformed (non-positive values become missing); log2-scale protein values are used as given. Detected genes without a value are imputed (median of detected non-partial genes by default; alternatives `mean`, `min`, `zero` (min − 1) and `exclude`) and flagged. Undetected genes are never imputed; an observed value is reported for them but not used by the context statistics. Abundance itself has no test, p-value or interval.

<!-- src: docs/methods/abundance.md:6-8,25-47; src/gippa/genetable/abundance.py:11-13,58-62,86-89,101-105,145-238,288-325; src/gippa/genetable/covariates.py:43-47,103; defaults.yaml:79-87; src/gippa/operons/discordance.py:65,133; src/gippa/operons/missed.py:45,152; src/gippa/genetable/abundance.py:183-196; src/gippa/genetable/abundance.py:7-9,101-106; src/gippa/context/regions.py:117 -->

## 8. Detectability model and detection limit

**Detectability.** The probability that a gene is detected, given its physicochemical properties only, is modelled for eligible (non-partial) genes. The covariates are log(1 + observable peptides), log length, GRAVY, pI (Biopython ProtParam [Cock et al. 2009]) and the number of transmembrane segments (Kyte–Doolittle window 19, mean ≥ 1.6, runs merged, capped at 10) [Kyte and Doolittle 1982]. Each covariate enters as a cubic spline with 5 quantile knots and linear extrapolation (`tm_segments` enters linearly when it has fewer than 5 distinct values), followed by L2-penalised logistic regression. The penalty $C \in \{0.01, 0.1, 1, 10\}$ is chosen by inner 3-fold TU-grouped log loss. Predictions are cross-fitted: TUs are assigned to 5 folds with a seeded permutation, and each gene's `p_detect` comes from the model not trained on its TU. With fewer than 50 detected or 50 undetected genes, or fewer TUs than folds, `p_detect` is the overall rate (flag `detectability_constant`). Reliability diagnostics (binned calibration, Brier score, recalibration slope and intercept) are descriptive. The calibration gate requires a mean cross-fitted slope in [0.8, 1.2].

<!-- src: src/gippa/genetable/detectability.py:40-42,199-215,357; src/gippa/genetable/covariates.py:110; defaults.yaml:89-92; docs/methods/detectability.md:9-82; tests/calibration/tolerances.yaml:4 -->

**Detection limit and P(miss).** For operon pairs in which exactly one partner has at least $k_{\min}$ observed peptides (Section 11), a peptide-level model gives the probability of detection as a function of abundance: $\mathrm{logit}\,p_j = b_0 + b_a a + b_L L_j + b_G \mathrm{GRAVY}_j$ over the observable peptides $j$ of a gene with log2 abundance $a$. A gene is detected when at least $k_{\min}$ peptides are seen (Poisson-binomial), with $k_{\min} = 2$ in the `two_peptide` tier and 1 otherwise. Because training genes are selected on detection, the model is fitted by unpenalised conditional (truncated) likelihood [Huggins 1989] and cross-fitted by TU folds. A replicate-based variant exists in the code (it is not identified with two runs) but is not wired into the 1.1.0 pipeline, so P(miss) always uses the peptide-level model. The calibration gate requires a mean gene-level calibration slope in [0.8, 1.2]; calibration-in-the-large is recorded but not gated in 1.1.0 (the tolerance key exists, but no test reads it).

<!-- src: src/gippa/operons/detection_limit.py:1-15,55-59,362-392; docs/methods/detection-limit.md:3-4,29-123; docs/adr/ADR-0030 decision 1; tests/calibration/tolerances.yaml:6; tests/calibration/test_detection_limit_citl.py:1-10; tests/calibration/test_detection_limit_calibration.py:87-98; src/gippa/steps/operons.py:26-27,196; src/gippa/operons/detection_limit.py:537-547 -->

## 9. Region tests

**Null hypothesis.** For each prophage, BGC, CRISPR and AMR region, GIPPA tests whether the region's genes are detected more often, or are more or less abundant, than genes at random positions on the same replicon. The null is that the region's position is exchangeable with every other position. Genes are indexed $0, \dots, N - 1$ in genome order among eligible genes of one replicon. A region spread over several replicons gives one test per replicon (`split_contigs`).

<!-- src: docs/methods/region-tests.md:10-17,58-70; src/gippa/context/regions.py:40,91-115,541-552 -->

**Detection statistic.** With detection indicators $d_i$ and cross-fitted probabilities $p_i$ of the region's genes, the reported effect is $S = \sum_i (d_i - p_i)$ (observed minus expected). The test statistic is studentised and recomputed at every shift:

$$T_{\mathrm{OE}} = \frac{\sum_i (d_i - p_i)}{\sqrt{\sum_i p_i (1 - p_i)}}.$$

<!-- src: src/gippa/context/regions.py:458-475; docs/adr/ADR-0018-region-tests-confounding.md:18-19 -->

**Abundance statistic.** Log2 abundance is first residualised per replicon by ordinary least squares on log length, log(1 + observable peptides) and GRAVY, fitted on detected genes. The effect is the Hodges–Lehmann (HL) shift, the median of all pairwise differences between residuals of the region's detected genes and the replicon's other detected genes. The test statistic is $T_{\mathrm{HL}} = \mathrm{HL}/\sqrt{V}$, with a variance that uses no region labels:

$$V = \frac{S_{\mathrm{in,in}}}{m^2} + \frac{S_{\mathrm{out,out}}}{m'^2} - \frac{2 S_{\mathrm{in,out}}}{m m'}, \qquad S_{AB} = \sum_{i \in A,\, j \in B} \hat\gamma\big(\delta(i,j)\big),$$

where $m$ and $m'$ are the numbers of detected genes inside and outside, $\delta(i, j)$ is the gene-index distance (circular on circular replicons), and $\hat\gamma(h)$ is the residual autocovariance at lag $h = 0, \dots, 30$, set to 0 from the first lag $h \ge 1$ at which it is ≤ 0. If $V \le 0$ at the observed or any null shift, the HL p-value is NA (`hl_variance_nonpositive`); if $\hat\gamma(0)$ is below a rounding floor, $V$ and $T$ are set to 0 at every shift. The HL effect and its interval are reported on the residual log2 scale.

<!-- src: src/gippa/context/regions.py:45-51,118-145,271-398,476-495; docs/methods/region-tests.md:26-56; docs/adr/ADR-0018-region-tests-confounding.md decision 1, Amendments 1-2 -->

**Null distribution and p-values.** The null shifts gene indices circularly, $i \mapsto (i + k) \bmod N$ for $k = 1, \dots, N - 1$; on linear contigs only non-wrapping shifts are used. Shifts are exhaustive up to 20,000 genes per replicon, which covers bacterial replicons, so the p-value is an exact randomisation p with floor $1/(1 + K)$; above that, 1,000 shifts are sampled. The shift keeps region size, member spacing and the replicon's autocorrelation, and changes only position. Two-sided p-values follow the +1 rule [Phipson and Smyth 2010],

$$p = \frac{1 + \#\{b : |T_b - \bar T| \ge |T_{\mathrm{obs}} - \bar T|\}}{1 + K},$$

with $\bar T$ the mean of observed and null statistics pooled; float ties count as exceedances. Intervals are 95 % percentile TU-block bootstraps (B = 1,000): O − E resamples the region's TUs, and HL resamples TUs inside and outside as two strata. With few TUs these intervals under-cover; calibration checks them against an expected small-sample coverage, not 0.95. BH is applied per (region type, statistic) over all regions of the run, at q ≤ 0.05. Regions with fewer than 3 eligible genes are flagged `small_region`, and those with fewer than 19 shifts (`few_shifts`) cannot reach p ≤ 0.05.

<!-- src: src/gippa/context/regions.py:42-43,401-403,454-464,469-471,486-493,556-560; src/gippa/stats/permutation.py:27-121; docs/methods/permutation-and-shifts.md:7-122; defaults.yaml:104-105,110; docs/adr/ADR-0018 decision 3; docs/preprint/results.md §1.2 (region tests: TU-block CIs against c(n)) -->

**Plasmids.** Plasmids are not shift-tested. Each plasmid gets a descriptive replicon contrast: O − E for detection and the HL shift of raw log2 abundance against chromosomal replicons, with p and q set to NA and the label `dosage_plus_regulation`. With user-supplied replicon depths, a `dosage_adjusted` row subtracts log2 of the depth ratio, treating depth as known. Region tests, autocorrelation, islands and the origin–terminus fits use the imputed abundance of detected genes (Section 7); among the genome-context statistics, only the plasmid contrast excludes imputed values.

<!-- src: src/gippa/context/regions_table.py:9,51-55,214; docs/methods/region-tests.md:194-222; src/gippa/genetable/abundance.py:215-238; src/gippa/context/regions.py:101,112,476-478; docs/methods/region-tests.md:198-200 -->

*Calibration status.* The region tests passed their full-tier simulation suites and attacks (including all-linear replicons and autocorrelated abundance), in one prophage stream within tolerance rather than at level (Results). The HL shift is relative to random positions on the replicon and partly reflects position along the origin-to-terminus dosage gradient, which the simulator does not model.

<!-- src: docs/methods/region-tests.md:130-143,178-181; validation/release/1.1.0/readiness.json rows R3,R4; docs/preprint/outline.md C12; docs/preprint/results.md §1.2 (prophage FDR 0.063 at q = 0.05) -->

## 10. Autocorrelation, abundance islands and the origin–terminus gradient

**Autocorrelation.** For genes with a value, centred log2 abundance $y$ is detrended twice: $y_w = y - \mathrm{LOESS}_{0.30}(y)$, then $x = y_w - \mathrm{LOESS}_{0.10}(y_w)$ over the ordinals of genes with values, with circular padding (local-linear, tricube, no robustness iterations). Moran's I [Moran 1950] at lag $k$ is

$$I_k = \frac{W_k^{-1} \sum_{(i,j):\,|i-j| = k} z_i z_j}{n^{-1} \sum_i z_i^2},$$

over the $W_k$ lag-$k$ pairs of genes with values, for lags 1–20, with $z_i = x_i - \bar x$ and $n$ the number of genes with values. Lags are counted along all eligible genes (wrapping on circular replicons), and pairs are formed only between genes with values; with 30 or fewer genes with values, Moran's I rows are NA (`few_values`). For detection, the reported value is the both-detected join count at lag $k$, but the test statistic is $S_k = \sum r_i r_j$ over lag-$k$ pairs, with Pearson residuals $r_i = (d_i - p_i)/\sqrt{p_i(1 - p_i)}$ from cross-fitted `p_detect`. Two nulls are used, each with B = 1,000: a gene-label permutation (values permuted among genes with values; $d$ and $p$ permuted together), and a TU-block permutation (maximal same-TU runs placed in random order, internal order kept). Each null replicate re-applies the 0.10 smoother. P-values are two-sided permutation p-values as above; BH runs over the 20 lags within each (replicon, statistic, null). Rows have no confidence interval (`ci_unestimated`).

<!-- src: src/gippa/context/autocorrelation.py:7-23,76-81,143-180,249-352; src/gippa/stats/autocorr.py; docs/methods/autocorrelation.md:19-185; docs/adr/ADR-0020-autocorrelation-redetrend-residual-joincount.md decisions 1-4, Amendment 1; defaults.yaml:104,111; src/gippa/context/autocorrelation.py:127-150,279-318; src/gippa/stats/autocorr.py:163-175 -->

*Calibration status.* With exact TUs, both nulls passed their gates at lags 1 and 20 and for any BH rejection; one exact-TU attack with large, heavy-tailed TUs passed only by two-block confirmation (Results). Under realistic TU-boundary errors, the TU-block Moran's I is anti-conservative at short lags (it passes at lag 20). Every TU-block Moran's I row therefore carries the flag `tu_block_tu_error`, and short-lag TU-block Moran's I results are not presented as calibrated. The join-count statistic and the gene-label null are unaffected.

<!-- src: docs/adr/ADR-0061-distance-model-fallback-scope.md:21-82; docs/methods/autocorrelation.md:141-147; src/gippa/operons/flags.py:254-269; src/gippa/steps/context_stats.py:197; tests/calibration/attacks/test_attack_operon_distance_fallback.py:114-117; docs/preprint/results.md §1.2 (autocorrelation: exact-TU attack unconfirmed) -->

**Abundance islands.** Within each replicon, eligible genes are ranked by log2 abundance with undetected genes tied at the lowest average rank, and converted to Blom normal scores [Blom 1958] $\Phi^{-1}\{(r - 3/8)/(n + 1/4)\}$. TU scores $z_t$ are the mean scores of their genes, ordered along the replicon and detrended by LOESS (span 0.10, circular padding). For every window of $w = 2, \dots, 20$ consecutive TUs (wrapping on circular replicons), $S_w = \sum_{t \in w} z_t / \sqrt{n_w}$. Significance uses a two-sided single-step maxT procedure [Westfall and Young 1993]: $M_b = \max_w |S_w^{(b)}|$ in null replicate $b$, and

$$p_w = \frac{1 + \#\{b : M_b \ge |S_w|\}}{1 + B}.$$

The null replicates are circular block bootstraps [Politis and Romano 1992] of each replicon's detrended TU series, with B = 1,000 and block length $\ell' = \min(20, 2\hat\ell + 1, T)$. Here $\hat\ell$ is the first lag at which Moran's I of the detrended TU scores has a raw two-sided permutation p > 0.05, and $T$ is the number of TUs. $M_b$ is the maximum over all scanned replicons, so one FWER covers the whole run. Windows with $p_w \le 0.05$ are significant; overlapping same-direction windows merge into one island, reported by its best window. The island effect is the mean detrended score of that window, with a TU-bootstrap interval that is subject to selection (winner's curse). A replicon with 30 or fewer scored TUs (so that the 0.10 LOESS span covers fewer than 4 points) is not scanned (`few_values`). At normal detection rates the low direction cannot reach significance, so the absence of low islands is not evidence of silencing.

<!-- src: src/gippa/stats/scores.py:19; src/gippa/context/islands.py:11-32,98-107,336-348,379,406-481,610-620; src/gippa/stats/scan.py:160-194; docs/methods/islands.md:31-157,209-270; docs/adr/ADR-0021-island-oriter-block-length.md decision 1; defaults.yaml:105,112 -->

*Calibration status.* The island FWER is claimed within a stated dependence regime, not for every genome. Within it (simulated AR(1) background with ρ ≤ 0.5 at every size measured, and a 4,500-gene chromosome at ρ = 0.9), the FWER at α = 0.05 stayed within the gate bound. Outside it, where strong TU autocorrelation is resolved only to a short block length ($\ell' \le 7$), the FWER is breached, as it is when a strong origin–terminus gradient combined with left censoring makes the local scale heterogeneous. These breaches are recorded and not flagged at run time; the run report gives $\hat\ell$, $\ell'$ and the lag-1 Moran's I per replicon, as diagnostics, so users can compare them with the calibrated cells ($\ell' \le 7$ at strong autocorrelation lies outside); no rule turns them into an inside/outside verdict.

<!-- src: docs/methods/islands.md:159-207,301-326; docs/adr/ADR-0021-island-oriter-block-length.md Amendment 1; docs/preprint/outline.md C13; docs/methods/islands.md:175-201,327-331; docs/methods/context_tables.md:63-89; src/gippa/report/content.py:236-249 -->

**Origin–terminus gradient.** For circular chromosomes, oriC comes from Bakta (several calls are merged by circular mean if all lie within 0.05 of the chromosome length, else no row) and the terminus is placed opposite. Each TU's replichore position runs from 0 (origin) to 1 (terminus). TU mean log2 abundance is regressed on position by OLS (pooled, right and left replichores; at least 3 TUs). Intervals are 95 % percentile circular block bootstraps of each fit's own TU sequence; $\hat\ell$ is estimated as in the island scan but once per replicon from the residuals of the pooled fit, and $\ell' = \min(20, 2\hat\ell + 1, T_s)$ uses each fit's TU count $T_s$. These rows are descriptive: no p-value, no q-value, no family. The interval's simulated coverage (measured with B = 199 instead of the default 1,000) is a few points below nominal. In 1.1.0 this coverage is recorded, not gated: the calibration test warns when pooled coverage falls below 0.90 (a full-arm gate at non-coverage ≤ `level_bound(0.07, 500)`, proposed in the methods note, is not implemented). The replichore rows cover less than the pooled row. The interval is not a test: with no gradient it excludes 0 in 8.3 % of simulated replicates.

<!-- src: src/gippa/context/oriter.py:72-95; docs/methods/oriter.md:6-151,194-204; src/gippa/context/oriter.py:22-27; docs/methods/oriter.md:126-180; tests/calibration/test_oriter_calibration.py:11-13,46-47,119-133; tests/calibration/test_oriter_calibration.py:46 -->

## 11. Operons

**Pair layers and flags.** Adjacent same-strand pairs are operon pairs (OP, $p \ge 0.5$) or non-operon pairs (NO). Each pair is annotated with an internal promoter call (Promotech window score ≥ 0.5 between the genes), a primary terminator call (TerminatorNet, or TransTermHP when chosen as the tool), predicted transcription rate, and the abundance difference $\Delta = a_{\mathrm{down}} - a_{\mathrm{up}}$ (log2). Pairs with a significant abundance pattern (discordant OP or co-abundant NO pair, q ≤ 0.10) receive the first matching flag. A tested OP pair with both partners detected that is not significant but has a terminator call between its genes is flagged `conflict_terminator_in_operon`. All other pairs get no flag row. The flags are `internal_promoter_supported`, `internal_terminator_supported`, `unexplained_discordance`, `conflict_terminator_in_operon`, `novel_operon_candidate`, `coabundant_despite_terminator`, or `other`.

<!-- src: src/gippa/operons/flags.py:44-104; src/gippa/operons/pairs.py:175-195; docs/methods/operons.md:22-26; src/gippa/operons/flags.py:75-82,107-125,161-172; src/gippa/annotate/terminators.py:24-25; src/gippa/operons/flags.py:64-71,121,172 -->

**Discordance (both partners detected).** The model is fitted on usable OP pairs: both genes eligible and detected, with measured (not imputed) abundance. Their $\Delta$ values are modelled as Student-t with location $\mu$, degrees of freedom $\nu \in [2, 50]$ and scale $\sigma_i = \sigma_0 \exp(\gamma z_i)$, $\gamma \in [-2, 2]$, where $z_i$ is the standardised pair-mean abundance; parameters are fitted by maximum likelihood (Nelder–Mead). A constant-σ model replaces it when pair-mean abundance varies little (SD < 0.1 log2) or the four-parameter fit fails; if the constant fit also fails, every OP pair is left untested (`too_few_op_pairs`). Every eligible OP pair with both partners detected, imputed ones included, is then tested against that fit with $p_i = 2\,\bar F_t\{|\Delta_i - \hat\mu|/\hat\sigma_i;\, \hat\nu\}$ ($z$ clamped to the fitted range), so the fitted pairs are also tested. At least 20 usable OP pairs are required. A relative abundance scale disables this test, the missed-operon family and P(miss) (`relative_scale`). BH runs over all tested pairs of the run at q ≤ 0.10.

<!-- src: src/gippa/operons/discordance.py:53-66,309; src/gippa/stats/tdist.py:33,41-47,239,334; docs/methods/operons.md:43-69,164-211; defaults.yaml:102; src/gippa/stats/tdist.py:273-306; src/gippa/operons/discordance.py:288-291; src/gippa/operons/discordance.py:100-110,292-307 -->

*Calibration status.* Discordance holds level only under a narrowed claim: smooth, monotone dependence of noise on abundance. A step change in noise with abundance fails, and that case is kept as a strict expected failure outside the claim. Power at realistic effects is near zero, so an empty discordance list is expected. Under the distance-model TU fallback the test is conservative (its null p-values are not uniform), and those rows carry `tu_source_distance_model`.

<!-- src: docs/adr/ADR-0055-operon-discordance-abundance-noise.md Amendment 1, :113-126; docs/methods/operons.md:83-125,1309-1322; docs/preprint/outline.md C14b -->

**One partner detected: P(miss).** For an OP pair in which exactly one partner $v$ has at least $k_{\min}$ observed peptides (and no t test applies), P(miss) is the probability, under the detection-limit model (Section 8) fitted out-of-fold for $v$'s TU, that the unseen partner would be missed given $v$'s abundance, averaged over within-operon abundance ratios. P(miss) is a probability, not a p-value. It receives no q, never enters the BH family and never makes a discordance call. Without peptide hits or a usable detection-limit model, the rows stay untested (`no_detection_limit_model`). Under protein-level detectability heterogeneity it is anti-conservative.

<!-- src: docs/methods/operons.md:189-200; src/gippa/steps/operons.py:19-24; docs/adr/ADR-0056-pmiss-outside-the-bh-family.md decisions 1-4, :46-104; src/gippa/operons/discordance.py:66; src/gippa/operons/discordance.py:331-333; src/gippa/steps/operons.py:13-27; src/gippa/operons/discordance.py:10-11,293-297; src/gippa/steps/operons.py:19-22 -->

**Missed operons.** Target pairs are NO pairs with both genes detected and finite abundance, scored by $s = -|\Delta|$. Decoy pairs, which cannot share a transcript, are convergent pairs and "sandwich" pairs (→ ← →) with both genes detected. Within distance bins $(-\infty, 0), [0, 50), [50, 150), [150, 300), [300, \infty)$ bp,

$$\widehat{\mathrm{FDR}}(t) = \frac{\pi_0 \sum_b N_{\mathrm{NO},b}\,\big(D_b(t) + 1\big)/\big(N_{\mathrm{decoy},b} + 1\big)}{T(t)}, \qquad \pi_0 = 1,$$

taking the element-wise maximum over decoy sets, then made monotone and clipped to [0, 1]. Calls require q ≤ 0.10, at least 20 scored NO pairs (partial genes excluded) and at least one usable decoy set; pairs with imputed partners are kept and flagged `abundance_imputed`. A kernel local FDR (Silverman bandwidth) is used for ranking only.

<!-- src: src/gippa/stats/decoy_fdr.py:139,208,367; src/gippa/operons/missed.py:41-48,150,222; docs/methods/operons.md:290-332,434-447; defaults.yaml:114-117; src/gippa/stats/decoy_fdr.py:114-136; src/gippa/operons/missed.py:183-188,223-241 -->

*Calibration status.* With the shipped decoys (convergent plus sandwich), the missed-operon family made no call in the gate or in any attack, so its calibration pass is vacuous. FDR control was shown only in a convergent-only arm that is not shipped. The procedure is near-powerless by construction. Co-imputed pairs (both partners median-filled, giving $\Delta = 0$) are a known risk: a pre-implementation simulation of a large genome with full co-imputation produced false calls, while the recorded attack M-1 made no call.

<!-- src: docs/methods/operons.md:357-371,482-551,1179-1277; docs/preprint/outline.md C14a; docs/methods/operons.md:504-533,1240-1267; tests/calibration/attacks/test_attack_operons.py:40-42,189-212 -->

## 12. Variants

**Pileup and calls.** Variants are called per genome codon from distinct I/L-folded peptide sequences with placement q ≤ 0.01 that are not `shared`. Residue alignments are recomputed by banded dynamic programming, and the recomputed cost must equal the mapper's difference count. At each codon, each peptide counts once as REF, IL, or one ALT key (substituted residue; or the token of a swap or block exchange at its first codon). A codon is ALT if it has ALT support and no REF, `mixed` if both, REF if only REF, and `no_call` otherwise. Frameshift events form their own rows. An ALT row is `high` when at least 2 distinct peptides support that ALT, else `candidate`. Variant-bearing peptides can only come from the tolerant route, so peptides shorter than 10 residues never carry a variant. Ambiguous letters (B, J, X, Z), other U or O residues, and a credited initiator M at an annotated non-ATG start give no category; the exceptions are U over a TGA codon that is a stop in the run's table and O over a TAG stop codon, which count as ALT (read-through). Under `genome_relation = same`, a database peptide searched against another strain keeps only its match and I/L operations.

<!-- src: src/gippa/variants/pileup.py:3-27,57-73,115-145,203-206,253-286,328-357; src/gippa/variants/align.py:1-14; defaults.yaml:65-71,99,120; src/gippa/variants/pileup.py:228-260; src/gippa/model/genetic_code.py:65-70; src/gippa/variants/pileup.py:20-22,115-131,252-254; src/gippa/variants/flags.py:181-196; src/gippa/variants/pileup.py:247-253 -->

**Flags.** Flags describe a call but never change it, with one exception. `mass_equivalent` marks an ALT whose mass matches REF, or REF plus a listed modification, within 0.05 Da; every swap is mass-equivalent, and so is any substitution that a supporting peptide also aligns as a mass-equivalent 2→1 block (ADR-0049). For substitutions (other than U/O read-through), `snv` or `multi_nt` records whether a single-base change explains the ALT. `strain_difference` marks every ALT under `genome_relation = related`, so it is not evidence of a strain difference. Other flags are `database_derived`, `start_codon`, `stop_codon` and `vcf_unrepresentable`. The exception is `putative_start`: a substitution to M at a non-annotated GTG or TTG codon, supported only at peptide position 0 and (inside a gene) not preceded by K/R, is capped at `candidate`. No modification is ever named. Outputs are a per-codon pileup, a variant table and a VCF 4.3 file holding substitutions and equal-length two-codon events.

<!-- src: src/gippa/variants/flags.py:6-9,38-73,86-96,150-202; src/gippa/variants/starts.py:29-31,79-105; src/gippa/steps/variants.py:100-105,153-155,209-255; src/gippa/variants/vcf_records.py:1-6; defaults.yaml:121; docs/adr/ADR-0049-m7-variant-precision.md decision 1 -->

**Class error rates.** When reference identifications are supplied, de novo reads that support an ALT are paired with database identifications of the same spectrum (target, q ≤ 0.01). The class error rate is the share of pairs whose peptides differ after I/L folding, per class (alt type × mass-equivalence). It estimates the error share among ALT-supporting de novo reads that a database search also identifies, not the false-ALT rate; it is not an upper bound and never changes confidence. Without paired spectra it is NA (`class_error_unestimated`).

<!-- src: src/gippa/variants/error_rates.py:46-122; src/gippa/steps/variants.py:147-153; docs/methods/variants.md:46-95,131-133 -->

*Calibration status.* Gate M7 requires F1 ≥ 0.90 for `high` ALT calls (and for planted novel ORFs). The variant part of gate M7 is not met in 1.1.0, and the release ships it as a disclosed limitation. A narrower scoring scope (unflagged `high` calls against identifiable truth) was accepted in an ADR, but no full-tier result under that scope is recorded. Variant calls should be read as candidates.

<!-- src: tests/calibration/tolerances.yaml:9; docs/adr/ADR-0049-m7-variant-precision.md (Context, Evidence); PROGRESS.md:199-208; docs/preprint/outline.md C16 and §6 -->

## 13. Novel ORFs and de novo precision

**Candidates.** Peptides that place outside annotated reading frames form novel-ORF candidates in seven classes, in priority order: `pseudogene` > `sorf` > `nterm_extension` > `cterm_extension` > `cds_outframe` > `antisense` > `intergenic`. A peptide is excluded when it, or its I/L twin, also maps in frame to a CDS or overlaps an RNA gene, is a contaminant, or maps to more than one novel locus. An ORF groups all placements that share (contig, strand, stop codon), extends in frame to the previous stop, and lists candidate starts (ATG, GTG, TTG) upstream of the first peptide. Two starts are reported: the longest ORF and the start with the highest OSTIR TIR.

<!-- src: src/gippa/novel/candidates.py:27-46,80-103,212-279; src/gippa/steps/novel.py:290-322; defaults.yaml:30 -->

**Evidence rules.** RBS support requires the TIR at the chosen start to reach the 0.10 quantile of TIR over detected non-partial genes (at least 20 such genes; otherwise `rbs_threshold_unestimated`, and N-terminal extensions and intergenic ORFs cannot be called). An N-terminal extension needs a start upstream of the annotated start and RBS support. An intergenic ORF needs a start, at least 10 codons and RBS support. Out-of-frame and antisense ORFs need at least 2 distinct peptides. A C-terminal extension needs a peptide that spans the annotated stop. Promoter (Promotech ≥ 0.5 in targeted windows) and downstream terminator evidence (a primary terminator call within 300 nt after the stop; BacTermFinder calls excluded) enter ranking only, never the call.

<!-- src: src/gippa/novel/evidence.py:7-8,28-31,41-129,192-204,346-369; defaults.yaml:36-38,44,124,126; docs/adr/ADR-0012 decisions 1-5 -->

**Class FDR.** Error control is per class; classes and sources are never pooled. In the `picked_decoy` regime, the run contains engine decoys of class-tagged `makedb` entries (four tagged classes: N-terminal extension, intergenic, out-of-frame, antisense). Each target ORF competes with its paired decoy entry on the best database-peptide q (engine score as tie-break; ties to the decoy), and $\widehat{\mathrm{FDR}}(t) = (D(t)/f + 1)/T(t)$ as in Section 6, with $f$ the paired share of the class's decoys; a class with $f < 0.80$ or with unscored decoys is unestimated. An ORF that loses its pick gets class FDR 1. In the `denovo_subgroup` regime, peptides of a class are ordered by $E$ and $\widehat{\mathrm{FDR}}_{\mathrm{class}}(e) = (W_{\mathrm{class}}(e) + c)/T_{\mathrm{class}}(e)$ as in Section 5.4, with the ORF q taken as the minimum peptide q. That choice is anti-conservative at ORF level. Otherwise the class FDR is unestimated. An ORF is called when its class rule holds and q ≤ 0.01; when q is unestimated, the class rule plus at least 2 distinct peptides is required (`class_fdr_unestimated`). In the database regime, the pseudogene, sORF and C-terminal-extension classes have no tagged `makedb` decoys and are always unestimated (`untagged_class`); the de novo regime estimates them like any other class. The q floor (1/T or c/T) means a class needs about 100 candidates to be callable at 0.01; a floor above α raises a warning. The flag `denovo_only_unsupported` marks ORFs supported only by de novo peptides with no further RBS support. Called ORFs are ranked by peptide count, best peptide q, start codon, TIR, promoter score, terminator confidence and length.

<!-- src: src/gippa/novel/fdr.py:53,108-185,310,372-430; src/gippa/novel/class_call.py:105-106,148,165-177; src/gippa/novel/evidence.py:60-63,218-232,304-310,355-375; src/gippa/stats/tda.py:117-140; defaults.yaml:98,125; docs/methods/novel-orfs.md:31-65,102-211,265-288; docs/adr/ADR-0010; src/gippa/novel/fdr.py:410-432; docs/methods/novel-orfs.md:64; src/gippa/novel/fdr.py:1-12,216-243; src/gippa/novel/evidence.py:41-51,68-85,103-125 -->

**De novo precision.** The genome itself estimates the precision of de novo peptides, as the correct share among distinct de novo peptides per score decile $b$:

$$\widehat{\mathrm{prec}}_b = \sum_\ell w_{b\ell}\, \frac{m_{b\ell} - r_{b\ell}}{\rho_\ell - r_{b\ell}},$$

clipped to [0, 1]. Here $m$ is the share with an exact target match, $r$ the share with an exact shuffled-decoy match, $\rho_\ell$ the exact-match share of reference database peptides (q ≤ 0.01) in length class $\ell$ ($\rho_\ell = 0.99$, flagged `c_assumed`, without reference peptides), and $w_{b\ell}$ the share of bin-$b$ peptides in class $\ell$. Classes with $\rho_\ell - r_{b\ell} < 0.05$ are dropped and the remaining weights renormalised; a bin with no class kept is unreported (its $\hat\pi = 1$). Intervals are 2.5/97.5 percentiles of 1,000 multinomial bootstrap resamples with $\rho_\ell$ held fixed. The interval ignores the sampling error of $\rho_\ell$ and dependence between nested peptides, and it is anti-conservative when reference peptides match less well than the de novo peptides.

<!-- src: src/gippa/mapping/denovo_precision.py:25-27,63-100,169-184; docs/methods/denovo-precision.md:24-150; src/gippa/mapping/denovo_precision.py:9-13 -->

*Calibration status.* The novel-ORF F1 criterion (≥ 0.90) passes, and so do the class-FDR checks in the database regime. In the de novo regime the ORF-level q (the minimum peptide q) is anti-conservative: it fails the attack with multi-peptide false ORFs at all three levels, and the attacks with exact systematic de novo errors (out-of-frame and antisense ORFs at q = 0.01 and 0.05, intergenic ORFs at 0.01, and attack GC at q = 0.01 in the full tier). All of these are kept as strict expected failures.

<!-- src: tests/calibration/tolerances.yaml:9; PROGRESS.md:204-207; docs/adr/ADR-0050-denovo-orf-q-exact-errors.md; tests/calibration/attacks/test_attack_novel_class_fdr.py:32-34,286-289; tests/calibration/test_novel_orfs_gate.py:24-28,223-226; PROGRESS.md:211-212; docs/adr/ADR-0050-denovo-orf-q-exact-errors.md:24-25; tests/calibration/test_m7_gate_shared_peptide.py:308-317 -->

## 14. Functional enrichment

**Collections.** GO [Gene Ontology Consortium 2023] (go-basic and goslim_generic, release 2026-05-19), KEGG pathways [Kanehisa et al. 2025] (REST API, unversioned; retrieval date and SHA-256 recorded) and COG 2020 functional categories [Galperin et al. 2021] are read from Bakta full-database cross-references. GO terms are propagated over `is_a` and `part_of`, with alternative IDs mapped to primary and obsolete terms dropped, giving the go_bp, go_mf and go_cc collections (namespace roots removed). KEGG global and overview maps (map011xx, map012xx, map01310, map01320) are excluded. Only eligible genes are members.

<!-- src: src/gippa/function/collections.py:5-7,26,116-192; src/gippa/data/db_sources.yaml:53-94; defaults.yaml:129; docs/methods/function-collections.md:5-29; docs/adr/ADR-0023; src/gippa/data/db_sources.yaml:66-79; src/gippa/function/collections.py:5-6 -->

**Detection over-representation (ORA).** The unit is the TU. A TU is detected if any eligible gene is detected, and belongs to a category if any eligible gene does. The raw TU detection probability is $q_t^{\mathrm{raw}} = 1 - \prod_g (1 - p_g)$ over its genes, which is recalibrated by a cross-fitted logistic model $\mathrm{logit}\, q_t = \beta_0 + \beta_1 \mathrm{logit}\, q_t^{\mathrm{raw}} + \beta_2 \log n_t$ on the TU folds (skipped, with a flag, when it cannot be fitted). Under the null, TUs are independent with probabilities $q_t$. The test uses the exact conditional law of the number $K$ of detected TUs in category $C$, given the total number $N$ of detected TUs:

$$\Pr(K = j \mid N) = \frac{\mathrm{PB}_C(j)\,\mathrm{PB}_{\bar C}(N - j)}{\sum_i \mathrm{PB}_C(i)\,\mathrm{PB}_{\bar C}(N - i)},$$

where $\mathrm{PB}$ denotes Poisson-binomial probabilities computed by dynamic programming. Over- and under-representation tails are reported for categories of 5–500 TUs. For GO over-representation, the elim procedure [Alexa et al. 2006] removes the TUs of a term with elim p < 0.01 from its ancestors. BH runs per (collection, direction). A gene-level Wallenius tail with goseq-style weighting [Young et al. 2010] is reported as `p_classic`; it is outside BH. On committed golden data, its tails match an exact Wallenius reference to a relative error of 10⁻⁶. The flexible ridge recalibration of ADR-0037 and the far-tail Wallenius urn recursion of ADR-0045, both accepted, are not implemented in 1.1.0.

<!-- src: src/gippa/function/ora.py:4-58,113-136,144-298,301-352,365-433; src/gippa/function/elim.py:1-5; src/gippa/stats/wallenius.py:1-10,24,62-90; defaults.yaml:101,130; docs/adr/ADR-0029; docs/adr/ADR-0035; docs/methods/function.md:318-513,665-668; docs/backlog.md:735-737; tests/golden/test_golden_wallenius.py:1-12; tests/calibration/tolerances.yaml:10; docs/adr/ADR-0037-ora-flexible-tu-recalibration.md:1-3 -->

*Calibration status.* ORA p-values are not calibrated: the recalibrated $q_t$ is misfit in several simulated cells (large TUs, high detection), while the conditional law itself is at level with oracle $q_t$. Every ORA row carries `ora_uncalibrated`, and the report states that these p-values are not discovery claims.

<!-- src: src/gippa/function/ora.py:49-52,413; docs/methods/function.md:17-26,589-663; docs/preprint/outline.md C17 -->

**Abundance GSEA.** Gene set enrichment [Subramanian et al. 2005] runs in R 4.5.3 with fgsea 1.36.2 `fgseaMultilevel` [Korotkevich et al. 2021] as a separate process (`sampleSize = 101`, `eps = 0`, `nPermSimple = 1000`, `scoreType = "std"`). The ranked unit is the TU, scored by the mean log2 abundance of its detected genes (imputed values included), centred on the median. Sets of 10–500 TUs are tested; a seed from the step's stream is passed to R. BH runs per collection with both signs of enrichment in one family. A significant GO set with significant descendants is re-tested without their genes and flagged `redundant` if it no longer passes. The claim is restricted: p-values are valid under TU exchangeability; sets that are clustered on the genome (`genomic_cluster`) and sets whose detection probability opposes their enrichment (`detection_confounded`) are flagged and carry no level claim; results are conditional on detection; NES has no interval. Runs use the R engine. Calibration uses two engines: the single-set level claims (TU null, operon-packed sets, ties), the gated structure attacks and the BH-family checks are run on R fgsea, whereas the flag checks (`detection_confounded`, `genomic_cluster`) the restricted-claim expected failures and the TU-fallback (distance-model, UniOP-like) GSEA arms run only on `fgsea_simple` (a Python port of fgsea that remains the default engine of the GSEA function when no engine is passed), as does the original TU-level suite (re-run on R for its level claims). The `detection_confounded` flag is a screen, not a guarantee. In mixed-confounding simulations on the port, unflagged rows still rejected at about 1.3α, 1.15α and 1.07α at α = 0.01, 0.05 and 0.10, and at about 1.1α under an unconfounded null; on that evidence a statistician decision gates the unflagged-row check at 1.5α, not at α. In the 1.1.0 full tier, this port-engine check passed within 1.5α, with one arm above α (Results); the screen has not been re-measured on the R engine. `genomic_cluster` is calibrated only for clustering up to φ = 0.5 (port). Flagged rows stay in the BH family; their effect on unflagged rows (spillover) has not been measured on the R engine, so no FDR or FWER claim is made for unflagged rows in a collection that contains flagged rows.

<!-- src: workflow/envs/fgsea.yml; workflow/bin/gippa_fgsea.R:2,8,105-109; src/gippa/steps/function.py:17-19,125-151; src/gippa/function/gsea.py:11-13,65-69,192-229,276-277,399-407,484-531,560-561; defaults.yaml:131-132; docs/methods/function-gsea.md:7-60; docs/adr/ADR-0025, ADR-0026, ADR-0028, ADR-0038; docs/methods/function-gsea.md:52,56; src/gippa/function/gsea.py:403-407; docs/methods/function-gsea.md:56; docs/methods/function-gsea.md:55-59; tests/calibration/attacks/test_attack_gsea_flags.py:367,384; docs/preprint/results.md §1.2 (GSEA null_b2); src/gippa/function/gsea.py:5-7,444,545; src/gippa/steps/function.py:17-19; tests/calibration/test_gsea_r_engine_calibration.py:1-12; tests/calibration/attacks/test_attack_gsea_r_engine.py:1-28; tests/calibration/test_gsea_tu_calibration.py:53-60,93-96; tests/calibration/attacks/test_attack_gsea_flags.py:72-75,116-123,366-386; tests/calibration/attacks/test_attack_gsea_structure.py:52,61; tests/calibration/attacks/test_attack_operon_distance_fallback.py:448-457 -->

**Category summary.** For COG, KEGG and GO slim categories, observed minus expected detected genes ($\sum p_{\mathrm{detect}}$) is reported with a 95 % TU-block bootstrap interval (B = 1,000), without a test. Intervals from fewer than 10 TUs are NA (`ci_few_units`); at 10–30 TUs they under-cover.

<!-- src: src/gippa/function/summary.py:26; docs/methods/function.md:687-722 -->

## 15. Determinants of detection and abundance

Model rows are detected genes with measured (not imputed) log2 abundance, excluding ribosomal proteins; at least 100 rows in 25 TUs are required. Thirty-one features form five blocks: transcription (promoter scores, predicted transcription rates, promoters between TU start and gene, TU position, distance to TU start), translation (OSTIR TIR, Shine–Dalgarno ΔG, spacer, start codon, mRNA ΔG), codon usage (CAI, tAI, ENC, Fop, GC3, codon-pair bias, AT content and rare codons at codons 2–10), protein (length, mass, pI, GRAVY, N-terminal residue, transmembrane segments) and context (replicon type, replichore position, prophage and BGC membership, local GC deviation). Univariate associations are Spearman ρ with abundance on the model rows and AUROC for detection over all eligible genes (ribosomal genes excluded for CAI only), with TU-block bootstrap intervals that are flagged `ci_unestimated`.

<!-- src: src/gippa/determinants/models.py:37-38,66-82,123-170; src/gippa/determinants/features.py:31-84; docs/methods/determinants.md:9-35; src/gippa/determinants/models.py:150-170 -->

Two models are fitted with 5-fold × 5-repeat cross-validation, with TU-grouped, size-balanced folds: an elastic net (l1 ratio ∈ {0.1, 0.5, 0.9}, 100 penalties, inner 3-fold TU-grouped CV, median imputation and standardisation) [Zou and Hastie 2005] and LightGBM (500 trees, learning rate 0.03, 15 leaves, minimum 20 samples per leaf, row and column subsampling 0.8, deterministic mode) [Ke et al. 2017]. Performance is the pooled out-of-fold R² per repeat (mean and SD). Block importance is the drop in R² when a block is removed, refitted on the same folds; the elastic-net value is the primary attribution, and the LightGBM value and per-gene SHAP contributions (LightGBM `pred_contrib`) [Lundberg and Lee 2017; Lundberg et al. 2020] are secondary checks. Genes whose residual from the repeat-averaged out-of-fold LightGBM prediction exceeds 2 residual SDs are flagged `abundance_not_explained_by_sequence`. A permutation null permutes whole TU response vectors among TUs with the same number of rows (TUs with a unique row count stay in place) and refits one CV repeat on the first repeat's folds. The observed statistic is the mean R² over the five repeats, which makes the null conservative for the elastic net. It gives $p = (1 + b)/(1 + B)$ with B = 1,000 per model, without BH.

<!-- src: src/gippa/determinants/models.py:175-245,278,355-405; src/gippa/determinants/attribution.py:5,39-50,234-237; src/gippa/determinants/null.py:34-37,143-192; defaults.yaml:134-140; docs/adr/ADR-0031, ADR-0033; src/gippa/determinants/null.py:1-7; docs/methods/determinants.md:68-90; src/gippa/determinants/attribution.py:232-235 -->

*Calibration status.* The elastic-net ΔR² attribution meets its bias gate (|bias| ≤ 0.05 for the planted block; noise blocks within 0.05 for both models). LightGBM's planted-block ΔR² is recorded, not gated, and overestimates it (about 0.98 against a true 0.80 in the full tier). The permutation p-value is not calibrated: it fails attacks with dependence between TUs and with detection-driven selection. Every permutation result carries `null_uncalibrated`, and drop-block ΔR² intervals carry `ci_unestimated` because repeats share data.

<!-- src: tests/calibration/tolerances.yaml:7; src/gippa/determinants/null.py:37,207; docs/methods/determinants.md:83-95; docs/adr/ADR-0034-determinants-null-uncalibrated.md; docs/adr/ADR-0031; docs/preprint/outline.md C18; tests/calibration/test_determinants_models_calibration.py:13-18; CAL/lanes/determinants-group1/determinants_m808.jsonl (planted_codon_delta_r2: lightgbm 0.976) -->

## 16. Sensitivity analyses

By default (`output.sensitivity: basic`), the regions, islands, operon, enrichment and determinant families are rerun under each alternative imputation policy (min and exclude) and, in the `two_peptide` tier, with single-peptide genes counted as detected. The `full` level adds a Top3 rerun of every family (when peptide intensities exist and Top3 is not already primary) and, for the operon family only, UniOP thresholds 0.3 and 0.7. Each rerun, and the recomputed primary, uses B = 199 and a common random-number stream per family. The single-peptide variant refits `p_detect` on the primary's folds; the other variants reuse the primary `p_detect`. Agreement is reported as the Jaccard index of significant sets, without an interval, and for determinants also as the agreement of dominant feature blocks.

<!-- src: src/gippa/context/sensitivity.py:10-21,50-53,88-94; docs/methods/sensitivity.md:15-46,248-252; defaults.yaml:84,150; src/gippa/context/sensitivity.py:16-21,47-48,95-117; src/gippa/steps/sensitivity.py:5-8,94; src/gippa/steps/sensitivity.py:4-6,86-98,118-120,152; src/gippa/context/sensitivity.py:96-123 -->

## 17. Comparing runs

`gippa compare` takes two or more output directories on the same genome (identical genome SHA-256 and locus-tag sets). It reports detection intersections, and for each pair of runs the Spearman ρ of measured log2 abundance on genes detected in both, with a 95 % TU-block bootstrap interval, and Bland–Altman bias and limits of agreement [Bland and Altman 1986]. Stability is the Jaccard index of five result sets (island genes, operon flags, enriched categories, called novel ORFs, ALT sites), described without a test. Runs with different output-schema major versions are refused unless forced.

<!-- src: src/gippa/compare/abundance.py:1-31,87-106; src/gippa/compare/stability.py:1-62; src/gippa/compare/versions.py:1-17; docs/methods/compare.md:7-36; docs/adr/ADR-0046 -->

## 18. Multiple-testing families

Every p- and q-value is reported, whether or not it passes. BH controls the FDR under independence or positive regression dependence [Benjamini and Hochberg 1995; Benjamini and Yekutieli 2001]; no π₀ is estimated, and no error rate is controlled across families. Table M2 lists the families.

**Table M2.** Error-control families in GIPPA 1.1.0.

| Claim | Family (scope) | Error rate | Level | Status |
|---|---|---|---|---|
| Identification (input) | the input search | search-engine q, used as given | 0.01 | not re-estimated by GIPPA |
| Peptide placement | one list per decoy type × stratum (0, 1, 2) | weighted target–decoy FDR | 0.01 | list-level calibrated; contaminant attack strict xfail; gate M3 not met |
| Gene detection | all eligible genes of the run | picked FDR, or MAYU (expected) FDR | 0.01 | calibrated at 0.01; `two_peptide`/`protein_list` unestimated |
| Novel ORFs | per class and per regime (database or de novo) | picked or weighted decoy FDR | 0.01 | de novo ORF q anti-conservative (strict xfails); three classes unestimated in the database regime |
| Variants | none (calls by peptide count) | placement q of supporting peptides only | 0.01 | gate M7 not met |
| ALT class error rate | per (ALT type, mass-equivalent) class | paired-spectrum error rate, not a bound | – | NA without paired spectra |
| De novo precision | per score bin | precision with percentile bootstrap interval | – | gated on bias (≤ 0.05) |
| Region tests | per (region type, statistic) | BH FDR | 0.05 | passes within tolerance (one prophage stream above q); intervals under-cover with few TUs |
| Plasmid contrasts | none | descriptive | – | – |
| Autocorrelation | 20 lags per (replicon, statistic, null) | BH FDR | 0.05 | TU-block Moran's I anti-conservative at short lags (flagged) |
| Abundance islands | all windows of all scanned replicons, both directions | FWER (single-step maxT) | 0.05 | within a stated dependence regime |
| Origin–terminus | none | descriptive interval | – | coverage slightly below nominal; recorded, not gated |
| Operon discordance | all both-detected OP pairs of the run | BH FDR | 0.10 | narrowed claim |
| P(miss) | none | probability, no q | – | outside BH; anti-conservative under heterogeneity |
| Missed operons | NO pairs of the run | decoy FDR | 0.10 | vacuous with shipped decoys |
| ORA | per (collection, direction) | BH FDR | 0.05 | `ora_uncalibrated` |
| GSEA | per collection, both signs | BH FDR | 0.05 | restricted claim (flags); unflagged rows screened at 1.5α (port engine); no FDR claim for unflagged rows beside flagged ones |
| Determinants | one p per model | permutation p, no BH | – | `null_uncalibrated` |
| Category summary, compare | none | intervals or descriptive | – | – |

<!-- src: docs/methods/multiple-testing.md:6-53; src/gippa/stats/multitest.py; src/gippa/context/regions.py:556-558; src/gippa/context/autocorrelation.py:348-352; src/gippa/context/islands.py:448-481; src/gippa/operons/discordance.py:309; src/gippa/operons/missed.py:150; src/gippa/function/ora.py:419-433; docs/methods/function-gsea.md:23-27; src/gippa/determinants/null.py:192; defaults.yaml:94-103 -->

Flags do not remove a row from its BH family: rows flagged `small_region`, `few_shifts`, `tu_block_tu_error`, `genomic_cluster` or `detection_confounded` stay in their families. Among p-values of BH-controlled modules, only P(miss), `p_classic` and `p_noelim` (the GO term's p before elim) are outside BH; the determinants permutation p has no family, and island windows use FWER (Table M2). No interval is computed for ORA, GSEA NES, discordance, missed-operon or determinant permutation rows, or for autocorrelation statistics; these rows carry `ci_unestimated` or report no interval.

<!-- src: src/gippa/context/regions.py:402-403,556-558; src/gippa/operons/flags.py:254-269; docs/methods/function-gsea.md:55-59; src/gippa/function/ora.py:366; src/gippa/function/gsea.py:632; src/gippa/operons/discordance.py:63; src/gippa/operons/missed.py:48,154; src/gippa/context/autocorrelation.py:19-20; src/gippa/variants/error_rates.py:1-7; src/gippa/mapping/denovo_precision.py:1-13; tests/calibration/tolerances.yaml:5; src/gippa/determinants/null.py:1-7; src/gippa/function/ora.py:72,367-369,419-425 -->

## 19. Report and outputs

Each run writes `README.txt`, `run_manifest.json`, `resolved_config.yaml`, and the directories `preflight`, `annotation`, `mapping`, `tables`, `plots`, `report`, `schemas` and `logs`; intermediates go to `intermediate/`. Tables are UTF-8 TSV with `NA` for missing values, snake_case names with units, floats written as `%.6g`, and a Frictionless data-package descriptor. Rows of the hypothesis-test tables (regions, autocorrelation, islands, operon pairs, ORA, GSEA) carry `stat, effect, ci_low, ci_high, p_value, q_value, method, null_model, n_perm`; the determinants permutation p is a row of a long-format model table. Thirty-four tables are documented in a generated reference. The genes table holds detection status, tier, gene q, peptide counts by uniqueness, `p_detect`, abundance and covariates; outputs say "detected", never "expressed".

<!-- src: workflow/modules/local/gippa_step/main.nf:56-60; docs/adr/ADR-0075-outputs-revision-1-1.md decision 3; src/gippa/io/tsv.py:33-56; GIPPA_SPEC.md:1297,2011-2053; docs/reference/tables.md:1-6; src/gippa/schemas/tables/genes.schema.json -->

The report is one self-contained HTML file (Jinja2, Plotly.js inlined, no network access when viewed). Its sections are summary, inputs and QC, mapping, gene detection, genome context, operons, variants, novel ORFs, function, determinants, sensitivity, methods and provenance. The run report holds fifteen plots (V01–V03, V05–V08, V10–V15, V17, V18), written as SVG, PNG (300 dpi) and PDF with their data tables; eleven have interactive versions. `gippa compare` adds V16 with its own interactive view. Each figure caption includes a structured guide (question, one item per panel, how to read it, caveats). Plots use Okabe–Ito colours for categories and viridis or cividis for continuous scales.

<!-- src: src/gippa/report/build.py:70-88,168-226; src/gippa/report/content.py:610-870; src/gippa/viz/interactive.py:21-25; src/gippa/viz/style.py:22-25,87-88; defaults.yaml:142-146; GIPPA_SPEC.md:2200-2252; src/gippa/cli.py:1525-1539 -->

## 20. Implementation, containers and reproducibility

**Architecture.** GIPPA is an MIT-licensed Python package (Python ≥ 3.11) with a command-line interface (`setup`, `preflight`, `run`, `makedb`, `convert`, `compare`, `report`, `example`, `demo`). The wheel ships a Nextflow workflow [Di Tommaso et al. 2017] (Nextflow 26.10.0 on the development cluster) with 12 subworkflows. `gippa run` writes all parameters to a JSON file, so Nextflow never parses user input. Each Python step is one Nextflow process; external tools are separate processes. Profiles cover local, Slurm, conda (micromamba), Apptainer and Docker execution. Tasks that run out of memory or time are retried at most twice with more resources. In 1.1.0, Python steps run in the host Python environment; step container images are defined (ADR-0065) but not wired into the workflow.

<!-- src: pyproject.toml:2-12,65-77; src/gippa/cli.py:256-1629; workflow/nextflow.config:1-9,37-88; workflow/conf/base.config:13-20; docs/adr/ADR-0002-farmshare-env-pins.md:11-12; workflow/modules/local/gippa_step/main.nf:1,73-75; docs/adr/ADR-0065-gippa-step-images.md -->

**Tools and licences.** Third-party tools other than the MIT-licensed codon-bias library (imported in-process) run as separate programs from pinned conda environments or from container images pinned by SHA-256 digest (17 Nextflow processes pinned to 11 images, plus a digest-pinned fgsea image started by the function step; UniOP is built locally and never published). A CI check fails if the source tree imports any GPL, AGPL, non-commercial, unlicensed or process-only tool, or the R bridge rpy2, or calls Rscript outside the fgsea launcher. Table M3 lists versions and licences.

**Table M3.** External tools in GIPPA 1.1.0.

| Tool | Version | Licence | Role | How run |
|---|---|---|---|---|
| Bakta | 1.12.2 (DB 6.0) | GPL-3.0 | genes, CRISPR, AMR hits, oriC | container or conda |
| geNOMAD | 1.12.0 (DB 1.9) | non-commercial | prophages, plasmids | container or conda |
| antiSMASH | 8.0.4 | AGPL-3.0 | BGCs | container or conda |
| UniOP | git 3ed14dfc | no licence file | operon pair probabilities | conda or locally built image; never redistributed |
| Promotech (RF-HOT) | git 56251ad9 | GPL-3.0 | promoter windows | GIPPA-built image or conda |
| Promoter Calculator | 1.2.4 | GPL-3.0 | transcription rates | GIPPA-built image or conda |
| TerminatorNet | git 72180309 | MIT | terminators | GIPPA-built image or conda |
| TransTermHP | 2.09 | GPL-3.0 | terminator cross-check | container or conda |
| BacTermFinder (opt-in) | git 02487044 | GPL-3.0 | terminators | GIPPA-built image or conda |
| OSTIR | 1.1.3 (ViennaRNA 2.6.4) | GPL-3.0 | translation initiation rates | GIPPA-built image or conda |
| miniprot | 0.18 | MIT | protein alignment | container or conda |
| BLAST+ (tblastn, blastp) | 2.16.0 | public domain | protein alignment, contaminants | container or conda |
| R / fgsea | 4.5.3 / 1.36.2 | GPL-2/3 / MIT | GSEA | subprocess (conda env or image) |
| codon-bias | 0.5.0 | MIT | codon indices | in-process |

Principal Python libraries (locked for Python ≥ 3.12): numpy 2.5.4, scipy 1.18.1, pandas 3.0.6, numba 0.68.0, scikit-learn 1.9.1, lightgbm 4.7.0, statsmodels 0.15.0, pyahocorasick 2.3.1, biopython 1.88, goatools 1.6.5. geNOMAD is licensed for non-commercial use, and KEGG data for academic use; other users can switch off the `mobile` module or use GO and COG only.

<!-- src: src/gippa/annotate/tools.py:66-125; workflow/conf/containers.config:19-50; workflow/conf/images.json:1-65; workflow/envs/*.yml; NOTICE.md:3-23; scripts/ci/check_imports.py:14-27; scripts/ci/check_no_rscript.py; src/gippa/io/fgsea_r.py:1-13; uv.lock (numpy :2133, scipy :3403, pandas :2261, numba :2009, scikit-learn :3268, lightgbm :1511, statsmodels :3628, pyahocorasick :2517, biopython :137, goatools :996); docs/preprint/outline.md §5; workflow/conf/containers.config:36-54; workflow/conf/images.json; scripts/ci/check_no_rscript.py:1-20; src/gippa/annotate/codon_usage.py:20-26; src/gippa/annotate/tools.py:113-125; scripts/ci/check_imports.py:1-27 -->

**Databases.** `gippa setup` downloads and checks these databases: Bakta 6.0 (light or full; Zenodo 14916843), geNOMAD 1.9 (Zenodo 14886553), antiSMASH 8.0.4 databases, the Promotech model, GO 2026-05-19, KEGG (unversioned; retrieval date and SHA-256 recorded), COG 2020 and the contaminant library (381 entries); `bakta_db download` also fetches the AMRFinderPlus database current at install time (2026-08-07.1 in the pilot runs); GIPPA records its version but pins no checksum for it. Pinned checksums (MD5 as published for Bakta and geNOMAD, SHA-256 otherwise) are verified on download; antiSMASH's databases are checked by its own downloader, and KEGG, which is unversioned, is recorded by retrieval date and SHA-256.

<!-- src: src/gippa/data/db_sources.yaml:6-109; src/gippa/cli.py:183-251; src/gippa/data/db_sources.yaml:9,17,25,30-36,72-85; pilot freeze run manifests (databases.amrfinderplus, read by audit_R2_methods_C.md); src/gippa/setupdb/run.py:313-315,557-565 -->

**Determinism.** All randomness flows through `for_step(step, seed)`, a NumPy generator seeded by `SeedSequence([seed, crc32(step)])`, so adding a step never changes another step's draws. LightGBM fits are single-threaded, deterministic and row-wise, and parallel fits are collected in a fixed order; the tolerant mapper's output does not depend on thread count. The run manifest records input checksums, tool and database versions, the resolved configuration, the seed, per-step status and run flags, and a `provenance_digest` (SHA-256 of a canonical JSON core). Volatile fields (command line, paths, timestamps, executor and profile, resume state, per-step records) are excluded from comparisons. Two runs with the same seed and inputs are required to give byte-identical tables (true example TE7, Section 22). This covers tables, not the report or plots.

<!-- src: src/gippa/rng.py:14-16; src/gippa/determinants/models.py:234; src/gippa/mapping/tolerant.py:13-17; src/gippa/manifest.py:29-34,110-155,424-467; defaults.yaml:2; tests/integration/test_m10_06_sim_te7.py:1-44; validation/pilot/te7.py:1-10; src/gippa/determinants/models.py:259-266; src/gippa/manifest.py:30-34 -->

**Example data.** `gippa example fetch` downloads a 1.09 MB bundle (SHA-256 pinned, also attached to the GitHub release) holding two 200-kb regions of the wild-type Nissle 1917 chromosome (GCF_021559835.1) with region subsets of the pilot's DIA-NN (arm A2 and the post-freeze six-frame arm A4, Section 23) and Cascadia outputs, for five input routes (DIA-NN report, protein FASTA plus abundance, Cascadia `.ssl`, generic table, six-frame DIA-NN search); the user converts the search outputs. Example runs use real tools at faster settings (core and BGC modules, light Bakta database, TransTermHP terminators, Promoter Calculator off).

<!-- src: validation/example/bundles.tsv:2; src/gippa/data/db_sources.yaml:122-131; src/gippa/data/example_routes.yaml:4-37; docs/adr/ADR-0074 decisions 1,3; docs/adr/ADR-0074-example-dataset.md:24-25,68-84; docs/adr/ADR-0074-example-dataset.md (Amendment 2026-11-06); validation/pilot/REPORT.md:340 -->

## 21. Simulation and calibration

**Simulator.** gippasim, a simulator written for testing, generates a genome and a proteome with known truth. It plants TUs (geometric size, mean 1.8 genes), genes (log-normal length, median 280 residues), start codons (ATG 0.90, GTG 0.08, TTG 0.02), prophages with a GC shift, RBS motifs, and strain substitutions and indels. Abundance is generated per TU with replicon, region, island and AR(1) background effects, plus planted sequence determinants of known R². Tryptic peptides of 7–30 residues are detected by a logistic model. Search errors are drawn at rate $\pi_{\mathrm{wrong}}$ from a shuffled or other-strain proteome, and de novo errors from an error model; wrong peptides from non-coding six-frame translations are an opt-in third source, and an ori–ter gradient term is available. It uses its own random streams.

<!-- src: tests/sim/gippasim/simulate.py:93-129,604-613; tests/sim/gippasim/features.py:24; GIPPA_SPEC.md:2297-2341; tests/sim/gippasim/proteomics.py:95-131; tests/sim/gippasim/expression.py:1-18 -->

**Criteria.** Every statistical method has a calibration suite over $R$ independent replicates: null replicates for level, uniformity, FDR and FWER criteria, and planted-truth replicates for accuracy criteria. The constants come from a protected file (`tests/calibration/tolerances.yaml`): levels α ∈ {0.01, 0.05, 0.10} and multiplier $k = 2$. Four criteria are used. Level, FDR and FWER are one-sided, so conservative results pass and anti-conservative results fail. Uniformity is a two-sided Kolmogorov–Smirnov test, so it fails conservative as well as anti-conservative continuous p-value distributions:

- *Level:* the share of $p \le \alpha$ is at most $\alpha + k\sqrt{\alpha(1-\alpha)/n}$, with $n = R$, or $n = Rm$ for pooled checks with $m$ p-values per replicate.
- *Uniformity:* the two-sided Kolmogorov–Smirnov test against U(0, 1) gives p > 0.001 (refused when more than 5 % of values are tied; discrete p-values use the level criterion only).
- *FDR:* mean FDP ≤ $q + k\,\mathrm{sd}(\mathrm{FDP})/\sqrt{R}$ (sample SD; $R \ge 2$).
- *FWER:* the share of replicates with at least one false discovery is at most $\alpha + k\sqrt{\alpha(1-\alpha)/R}$.

Additional fixed criteria are: detectability slope in [0.8, 1.2]; de novo precision |bias| ≤ 0.05; detection-limit |calibration-in-the-large| ≤ 0.2 (a tolerance key not yet read by any test, so recorded only); determinants ΔR² |bias| ≤ 0.05; mapping recall ≥ 0.999; and F1 ≥ 0.90 for `high` variant calls and for planted novel ORFs. The harness also has pooled-level, pooled-F1 and FDR-from-counts variants, and a pooled check ends `vacuous` when its mean number of calls falls below a set floor. Non-finite or out-of-range inputs raise errors. Replicate $i$ of suite $s$ is seeded by `SeedSequence([20261014, crc32(s), i])`, so any shard reproduces any replicate.

<!-- src: tests/calibration/tolerances.yaml:1-11; tests/calibration/_criteria.py:1-60,177-273; docs/methods/calibration.md:17-24,165-180; GIPPA_SPEC.md:1956-1959; tests/calibration/_criteria.py:10-17,243-257; tests/calibration/_criteria.py:96,960-962; scripts/ci/check_protected.py:6; tests/calibration/_criteria.py:10-11,281-307 -->

**Two-block confirmation and replicate counts.** For checks run through the harness's `calibrate` gate, a failure on replicates $0, \dots, R-1$ triggers a second block $R, \dots, 2R-1$; the check fails only if both blocks fail, and is otherwise marked `unconfirmed`. This limits false alarms across many at-level checks. Fifteen suites run at a fast tier (nightly) and a full tier (one Slurm array task per suite, run on demand; Table M4). Release readiness requires both tiers to pass on the release commit or a reusable ancestor; known failures, in methods shipped with a flag or a narrowed claim, are kept as expected failures (xfail; next paragraph).

**Table M4.** Replicates per calibration suite (fast / full).

| Suite | R | Suite | R | Suite | R |
|---|---|---|---|---|---|
| placement | 50 / 300 | autocorrelation | 100 / 1,000 | gsea | 50 / 300 |
| gene_fdr | 50 / 300 | islands | 50 / 500 | determinants | 20 / 100 |
| detectability | 20 / 100 | discordance | 100 / 1,000 | novel_orfs | 50 / 300 |
| regions | 100 / 1,000 | missed_operons | 50 / 300 | denovo_precision | 50 / 300 |
| ora | 100 / 1,000 | variants | 20 / 100 | compare | 100 / 1,000 |

Some checks use larger R than their suite default (for example, compare streams with 10,000 or 20,000 replicates; Results).

<!-- src: tests/calibration/_criteria.py:27-46,669,1208-1227; docs/adr/ADR-0015-calibration-two-block-confirmation.md:12-16; tests/calibration/tolerances.yaml:12-19; scripts/slurm/calibration_full.sbatch:5,10-12; docs/adr/ADR-0042; docs/preprint/outline.md C3; scripts/release/readiness.py:339-369,698-714; validation/release/1.1.0/READINESS.md R3; tests/calibration/attacks/test_attack_operons_sigma.py:49; tests/calibration/attacks/test_attack_operon_distance_fallback.py:117; docs/preprint/results.md §1.2 (H1 and free-entrapment streams) -->

**Attack nulls.** Each method also faces attack nulls: null scenarios built to fool it, which it passes only by staying at level. The specification names six: autocorrelated abundance without islands; detectability confounded with region membership; gene sets packed into one operon; high-copy plasmids; de novo errors that reverse a stretch; and isobaric swaps everywhere. The attack suite (55 files) extends these per method. It also tests the harness: an oracle BH run at FDR 1.5q must be a confirmed failure. A known failing attack is kept as an expected failure (xfail) tied to an ADR or, for two attacks (the contaminant attack `flat_c200` and multi-peptide false ORFs), to an open design decision with no ADR. Most of these marks are strict: such a test fails the build if the method starts passing, so a fix must update the record.

<!-- src: GIPPA_SPEC.md:1980-1991; tests/calibration/attacks/ (55 test files); tests/calibration/attacks/test_attack_calibration_harness.py; docs/methods/calibration.md:158-161; tests/calibration/_criteria.py docstring (ADR-0057 xfail handling); PROGRESS.md Gate M3 and Gate M7 blocks; tests/calibration/attacks/test_attack_novel_class_fdr.py:25-34 -->

**Scope of the evidence.** All calibration claims rest on simulation with gippasim, a simulator built within the project (written in the implementer role and tested in the test-author role), so each claim holds only under its generative model; the attack nulls probe departures from that model but do not remove the dependence. The full validation set out in the specification was not done: the K-12 dataset PXD075568, an FDRBench entrapment benchmark and the K-12 co-abundance AUROC; the sealed EDL933 holdout (PXD068911) was not used, and no holdout claim is made. Operon and terminator cut-offs were calibrated on K-12 annotations only (RegulonDB). The only real-data evaluation is one Nissle 1917 sample (three runs) analysed in three gated arms (two DIA-NN searches and one de novo arm), plus an exploratory post-freeze six-frame arm (Section 23).

<!-- src: AGENTS.md (Conventions: the simulator is written by the implementer, tested by the test-author); docs/preprint/outline.md:139-145 (§6 Limitations: full validation not done; sealed holdout); docs/adr/ADR-0022-uniop-terminatornet-cutoffs.md (K-12 RegulonDB 14.5.0); scripts/release/readiness.py R10 (deferred) -->

## 22. True-example assertions and release readiness

**True examples.** Eight assertions, mostly on real *E. coli* genomes, check biology that any correct run must reproduce (Table M5).

**Table M5.** True-example assertions.

| ID | Assertion | Where evaluated |
|---|---|---|
| TE1 | ≥ 99.5 % of A2 peptides map exactly in frame | pilot (readiness R1); in-silico proxy |
| TE2 | no `high` ALT call except at `start_codon`/`stop_codon` sites | pilot; in-silico proxy (weekly, on GCA_043228065.1) |
| TE3 | ribosomal proteins ≥ 90 % detected, median abundance in the top quartile | pilot |
| TE4 | log2(RplL/RplJ) ∈ [1, 3]; log2(AtpA/AtpD) ∈ [−1, 1] | pilot |
| TE5 | trypsin and keratin peptides land in the contaminant bucket | pilot; simulation |
| TE6 | K-12: Bakta oriC next to mnmG | weekly K-12 run |
| TE7 | the same seed gives byte-identical tables | simulated 100-kb genome (integration test; ran and passed in the nightly integration job 1973700_0 on the release code); second pilot A2 run at the freeze commit |
| TE8 | UniProt K-12 peptides ≥ 99.5 % exact on K-12; Nissle placement FDR ≤ 1 % | in-silico proxy in the weekly K-12 job (Nissle side on GCA_043228065.1) |

<!-- src: GIPPA_SPEC.md:2369-2380; validation/proxies/te_proxies.py:67-68,280-317,402-470; validation/pilot/te_real.py:42-43,93,201,238,322; tests/contract/test_te6_oric.py; scripts/slurm/weekly.sbatch:19-21; tests/integration/test_m10_06_sim_te7.py; tests/data/genomes/PROVENANCE.md:8,20; scripts/slurm/weekly.sbatch:11,200-215; docs/preprint/results.md §3 (TE7); local/ci/runs/1973700_0/integration.txt (commit e43683c7; 967 passed, 84 skipped, none of them TE7; 0 failed); scripts/slurm/nightly.sbatch:207-214 -->

**Readiness.** A release script computes rows R1–R10 from the specification's acceptance criteria and housekeeping rows H1–H6. R1: A2 exact mapping ≥ 99.5 %. R2: decoy-estimated mapping FDR ≤ 1 %. R3: the full and nightly calibration tiers both pass (no failed test). R4: every attack null stays at level. R5: TE1–TE5, TE7 and TE8 pass (TE6 is checked in the weekly job, not in R5). R6: byte-identical tables (TE7) and passing integration runs. R7: median K-12 standard run time 16–55 min at 16 threads and 64 GB. R8: core-only run time (reported only). R9: pilot passes. R10: co-abundance (deferred). The H rows cover milestone gates, version agreement, changelog, digest-pinned containers, licences and CI status. A row may be waived only by an accepted ADR; a waived row is reported as waived, not as a pass, although the overall verdict counts it as non-failing. Calibration evidence from an earlier commit is reused only when every path changed since that commit lies on a fixed allow-list of non-statistical paths (documentation, release paperwork, CI and agent configuration, report and plotting code), apart from a one-hunk `__version__` change in `src/gippa/__init__.py`; 1.1.0 relies on that exception. For 1.1.0, three rows were waived. The release table cites the 1.0.0 waiver records (ADR-0070, ADR-0071, ADR-0072), which still apply because the script does not check a waiver's version; ADR-0076 extends their scope to 1.1.0: R2, because A2's peptides come from a same-strain database and are exempt from decoy-based placement estimation, so A2's mapping FDR was not measured; R5, for TE4's RplL/RplJ ratio (the waiver covers the whole row; the release lead checked by hand that only TE4 failed); and R7, runtime.

<!-- src: scripts/release/readiness.py:25,44-62,86-101; validation/release/1.1.0/READINESS.md:1-24; docs/adr/ADR-0076-waivers-1-1.md:20-22; docs/adr/ADR-0070, ADR-0071, ADR-0072; docs/preprint/outline.md R3 rule; scripts/release/readiness.py:39-55,392-451; scripts/release/readiness.py:45-52,698-715; docs/adr/ADR-0076-waivers-1-1.md (Context); scripts/release/readiness.py:65-66,186-203 -->

**Continuous integration.** Pull requests run lint, strict type checks (plus a macOS test job), licence-import checks, role and protected-file guards, and unit, property, golden, contract and snapshot tests (snapshot RMS ≤ 2), plus nf-test stub runs and a package build. Calibration and the coverage gate (≥ 0.90) run nightly on the cluster, together with integration tests, benchmarks (fail above 15 % slowdown) and a real run on a simulated genome. A weekly job runs K-12 MG1655 (U00096.3) and Nissle 1917 with every tool and the full Bakta database; its Nissle genome is the committed test assembly GCA_043228065.1 (an engineered chassis), not the wild-type assembly of the pilot.

<!-- src: .github/workflows/ci.yml:20-166; ci/gates.yaml:1-13; docs/adr/ADR-0052-pr-ci-under-20-minutes.md:18-24; scripts/slurm/nightly.sbatch:2-56; scripts/slurm/weekly.sbatch:2-37; tests/data/genomes/PROVENANCE.md:8,20; scripts/slurm/weekly.sbatch:11; .github/workflows/ci.yml:139-145 -->

## 23. Pilot study design

**Data.** The pilot used three wild-type *Escherichia coli* Nissle 1917 runs from PXD078459 (iProX IPX0017169000), acquired on a timsTOF HT; their diaPASEF acquisition was confirmed from the raw files. The reference genome was the wild-type assembly GCF_021559835.1 (chromosome NZ_CP058217.1 and plasmids pMUT1 NZ_CP058218.1 and pMUT2 NZ_CP058219.1), with all records marked circular. It was annotated with Bakta 1.12.2 and the full database 6.0. An earlier assembly named in the specification (GCA_043228065.1) proved to be an engineered chassis, so the pilot was re-run on the wild-type assembly. Nissle 1917 belongs to phylogroup B2 and K-12 to phylogroup A, so a K-12 database is a genuine mismatch.

<!-- src: validation/pilot/README.md:4,134,161-163,532-560; validation/pilot/inputs/pxd078459_files.tsv; GIPPA_SPEC.md:2603-2607; validation/pilot/REPORT.md "Deviations at the freeze" -->

**Arms.** Arm A1 (generic reference) searched the UniProt K-12 proteome UP000000625 plus contaminants. Arm A2 (matched strain) searched the Nissle Bakta proteins built by `gippa makedb` (4,613 CDS and 72 sORFs; no decoys), plus contaminants. Both searches used DIA-NN 2.3.1 Academia [Demichev et al. 2020] in place of the specified 2.7.0. They ran in two steps: library prediction from the FASTA files, then search of the runs. Settings were identical apart from the database: trypsin (K*, R*), one missed cleavage, peptides of 7–30 residues, fixed carbamidomethylation, up to one variable methionine oxidation, precursor charge 1–4, precursor m/z 300–1800, fragment m/z 200–1800, q ≤ 0.01, and decoy reporting. Peptidoforms and matrices output were enabled; initiator-methionine excision and match-between-runs were not used. DIA-NN ran from a checksum-pinned Apptainer image. Arm A3 (de novo, exploratory) ran Cascadia 0.0.7 [Sanders et al. 2025] on the same runs, with the checkpoint and settings fixed by a prior feasibility check. DIA-NN and Cascadia outputs were converted with `gippa convert`. GIPPA ran each arm with the default `genome_relation = related`; A3 received the A2 peptides as reference peptides for de novo precision. The A1 DIA-NN search and the Cascadia results do not depend on the Nissle genome and were not repeated after the genome change; the A2 database and search were rebuilt. A1, A2 and A3 were then compared with `gippa compare` (metrics on A1 versus A2). An exploratory six-frame arm (A4: DIA-NN 2.3.1 against a `gippa makedb --six-frame-orfs` database) ran after the freeze and is not part of the gate.

<!-- src: GIPPA_SPEC.md:2609-2613,2645-2647; validation/pilot/diann_A1.cfg; validation/pilot/diann_A2.cfg; docs/adr/ADR-0051-diann-version.md (Context); validation/pilot/README.md:400-430,560-604,877; validation/pilot/REPORT.md §2; validation/pilot/reproduce.md:7,183-230 (S8-S12) -->

**Freeze and provenance.** All arms, the comparison and the metrics ran at one frozen commit (fa49b8f1, 12 commits before the v1.0.0 tag; its run manifests report version 1.0.0), except that the A1 and A3 input tables are conversions made in earlier passes (pass 1 and pass 2) and the entrapment re-run (commit db3a95f1, identical `src/` and `workflow/`, fixed entrapment tool), recorded with job IDs in `validation/pilot/reproduce.md`. Release 1.1.0's readiness check reads these 1.0.0-era pilot metrics directly. ADR-0076 extends only the R2, R5 (TE4) and R7 waivers; no decision record covers reusing the passing rows. Between the freeze and 1.1.0, output columns were renamed and further code changed. Row-level inputs and outputs were not committed; only aggregates and checksums were.

<!-- src: validation/pilot/REPORT.md §1, §9, §10; docs/adr/ADR-0076-waivers-1-1.md; docs/adr/ADR-0075-outputs-revision-1-1.md; $GIPPA_ROOT/runs/pilot/freeze/fa49b8f/out_*/run_manifest.json (gippa_version, read by audit_R2_methods_C.md); git tag v1.0.0; validation/pilot/README.md:603; validation/pilot/REPORT.md:306-354; docs/adr/ADR-0076-waivers-1-1.md:1,27-31 -->

**Metrics and pass criteria.** The pilot report covered: identifications at 1 % FDR; mapping shares; detection coverage of strain CDS (including plasmid, prophage and BGC genes); A1–A2 Spearman ρ; Jaccard stability; novel ORFs and variants by class; runtime; and A1's free entrapment. For the entrapment, peptides absent from the Nissle translation cannot be true if the genome is complete. With $N_E$ such peptides, $N_T$ peptides present, and $r$ the in-silico ratio of K-12-only to shared tryptic sequences, computed at the same level (peptide or precursor), the combined estimator [Wen et al. 2025] is $\widehat{\mathrm{FDP}} = N_E(1 + 1/r)/(N_T + N_E)$. The pass criteria were: GIPPA runs end to end on A1, A2 and A3; at least 99.5 % of A2 peptides map exactly (TE1); A1 and A2 are compared; measured runtimes replace estimates; and TE2–TE5 are reported, with failures explained.

<!-- src: GIPPA_SPEC.md:2671-2694; validation/pilot/entrapment.py; docs/methods/free-entrapment.md; validation/pilot/metrics/metrics.json › m8.universe (r, r_prec) -->

**Pre-registered hypotheses.** Four hypotheses were fixed before the pilot ran: (H1) A2 detects more strain CDS than A1, mostly plasmid, prophage and BGC genes; (H2) A1's genome-inconsistent FDP is at least its entrapment FDP; (H3) context conclusions (islands, region tests, operon concordance) agree between A1 and A2 on shared genes; (H4) the genome-anchored de novo precision agrees with precision measured against A2 identifications. H2 required an entrapment benchmark (FDRBench) and was deferred to a full validation. H1 and H4 were evaluated. H3 could not be tested as stated (no islands or operon flags in either arm, and region tests were not compared); the substitute criteria used are reported in Results.

<!-- src: GIPPA_SPEC.md:2696-2703; docs/preprint/outline.md C25-C28; validation/pilot/REPORT.md §10; docs/preprint/results.md §4.4 (H3 not tested as pre-registered); validation/pilot/metrics/metrics.json › h3 -->

## 24. Data and code availability

The proteomics data are PXD078459 (ProteomeXchange [Deutsch et al. 2023], iProX [Chen et al. 2022]), cited by accession and not redistributed. Genomes: *E. coli* Nissle 1917 GCF_021559835.1 (pilot), the engineered Nissle chassis GCA_043228065.1 (weekly tests, in-silico proxies and committed test data only) and K-12 MG1655 GCA_000005845.2 (U00096.3); proteome: UniProt UP000000625 [UniProt Consortium 2025]. RegulonDB 14.5.0 was used for the K-12 operon and terminator calibration and is not redistributed, under its terms. GIPPA is MIT-licensed; the source code, specification, methods notes, ADRs, calibration suites and pilot aggregates are at https://github.com/khoa-yelo/gippa. GIPPA 1.1.0 is on PyPI (`pip install gippa==1.1.0`) and attached to the GitHub release v1.1.0 (wheel, sdist and example bundle; tagged and published 2026-11-12 UTC). An archival DOI: to be added. DIA-NN is proprietary (academic licence) and is not distributed with GIPPA.

<!-- src: CITATION.cff:12-18; docs/preprint/outline.md §5; validation/pilot/README.md:556; GIPPA_SPEC.md:2645-2653; docs/adr/ADR-0022 (RegulonDB 14.5.0); docs/adr/ADR-0051; PyPI project gippa 1.1.0 and GitHub release v1.1.0 (verified by audit_R1_methods_C.md F35-F36); tests/data/genomes/PROVENANCE.md:8 -->

## 25. Software development process

GIPPA was developed against a written specification (`GIPPA_SPEC.md`, read-only), with architecture decision records (ADRs) that override it when accepted; the project owner accepts any change to defaults, thresholds, interfaces, protected files or releases. Statistical methods have methods notes, reviewed before code. Features followed a test-first loop: issue specification, failing tests, the smallest passing implementation, then independent code review and statistical review. The statistical reviewer can block a change and owns the attack nulls. Most of the work was done by AI coding agents in separate roles (specification writer, test author, implementer, code reviewer, statistician, and others). Each commit records its role in an `Agent-Role` trailer. CI enforces that no role changes both source code and tests, that implementers do not edit existing tests, and that protected files (specification, defaults, calibration cut-offs and tolerances, CI gates) change only with an ADR the owner has accepted, named in a commit trailer (adding a new default key needs only an ADR reference).

<!-- src: AGENTS.md "Source of truth"; GIPPA_SPEC.md:2483-2578; .claude/agents/ (9 files); scripts/ci/check_roles.py; scripts/ci/check_protected.py; ci/protected.txt; .github/workflows/ci.yml:36-83; docs/preprint/outline.md §1 (AI-assistance statement: Khoa to write); ci/protected.txt; scripts/ci/check_protected.py:1-8 -->

## References (methods)

- Aho AV, Corasick MJ (1975). Efficient string matching: an aid to bibliographic search. *Commun ACM* 18:333–340.
- Alexa A, Rahnenführer J, Lengauer T (2006). Improved scoring of functional groups from gene expression data by decorrelating GO graph structure. *Bioinformatics* 22:1600–1607. doi:10.1093/bioinformatics/btl140
- Barber RF, Candès EJ (2015). Controlling the false discovery rate via knockoffs. *Ann Stat* 43:2055–2085.
- Benjamini Y, Hochberg Y (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. *J R Stat Soc B* 57:289–300. doi:10.1111/j.2517-6161.1995.tb02031.x
- Benjamini Y, Yekutieli D (2001). The control of the false discovery rate in multiple testing under dependency. *Ann Stat* 29:1165–1188.
- Bland JM, Altman DG (1986). Statistical methods for assessing agreement between two methods of clinical measurement. *Lancet* 327:307–310.
- Blin K, et al. (2025). antiSMASH 8.0: extended gene cluster detection capabilities and analyses of chemistry, enzymology, and regulation. *Nucleic Acids Res* 53(W1):W32–W38. doi:10.1093/nar/gkaf334
- Blom G (1958). *Statistical Estimates and Transformed Beta-Variables.* Wiley.
- Camacho C, et al. (2009). BLAST+: architecture and applications. *BMC Bioinformatics* 10:421.
- Camargo AP, et al. (2024). Identification of mobile genetic elements with geNomad. *Nat Biotechnol* 42(8):1303–1312. doi:10.1038/s41587-023-01953-y
- Chen T, et al. (2022). iProX in 2021: connecting proteomics data sharing with big data. *Nucleic Acids Res* 50(D1):D1522–D1527. doi:10.1093/nar/gkab1081
- Chevez-Guardado R, Peña-Castillo L (2021). Promotech: a general tool for bacterial promoter recognition. *Genome Biol* 22:318. doi:10.1186/s13059-021-02514-9
- Clopper CJ, Pearson ES (1934). The use of confidence or fiducial limits illustrated in the case of the binomial. *Biometrika* 26:404–413.
- Cock PJA, et al. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics* 25:1422–1423.
- Demichev V, et al. (2020). DIA-NN: neural networks and interference correction enable deep proteome coverage in high throughput. *Nat Methods* 17(1):41–44. doi:10.1038/s41592-019-0638-x
- Deutsch EW, et al. (2023). The ProteomeXchange consortium at 10 years: 2023 update. *Nucleic Acids Res* 51(D1):D1539–D1548. doi:10.1093/nar/gkac1040
- Di Tommaso P, et al. (2017). Nextflow enables reproducible computational workflows. *Nat Biotechnol* 35:316–319. doi:10.1038/nbt.3820
- dos Reis M, Savva R, Wernisch L (2004). Solving the riddle of codon usage preferences: a test for translational selection. *Nucleic Acids Res* 32:5036–5044. doi:10.1093/nar/gkh834
- Feldgarden M, et al. (2021). AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence. *Sci Rep* 11:12728. doi:10.1038/s41598-021-91456-0
- Frankenfield AM, et al. (2022). Protein contaminants matter: building universal protein contaminant libraries for DDA and DIA proteomics. *J Proteome Res* 21:2104–2113. (PMC10040255)
- Galperin MY, et al. (2021). COG database update: focus on microbial diversity, model organisms, and widespread pathogens. *Nucleic Acids Res* 49(D1):D274–D281. doi:10.1093/nar/gkaa1018
- Gene Ontology Consortium (2023). The Gene Ontology knowledgebase in 2023. *Genetics* 224(1):iyad031. doi:10.1093/genetics/iyad031
- He K, et al. (2015). A theoretical foundation of the target-decoy search strategy for false discovery rate control in proteomics. arXiv:1501.00537.
- Huggins RM (1989). On the statistical analysis of capture experiments. *Biometrika* 76:133–140.
- Kanehisa M, et al. (2025). KEGG: biological systems database as a model of the real world. *Nucleic Acids Res* 53(D1):D672–D677. doi:10.1093/nar/gkae909
- Ke G, et al. (2017). LightGBM: a highly efficient gradient boosting decision tree. *Adv Neural Inf Process Syst* 30.
- Kingsford CL, Ayanbule K, Salzberg SL (2007). Rapid, accurate, computational discovery of Rho-independent transcription terminators illuminates their relationship to DNA uptake. *Genome Biol* 8:R22.
- Korotkevich G, et al. (2021). Fast gene set enrichment analysis. *bioRxiv* 060012, version 3 (first posted 2016). doi:10.1101/060012
- Kyte J, Doolittle RF (1982). A simple method for displaying the hydropathic character of a protein. *J Mol Biol* 157:105–132. doi:10.1016/0022-2836(82)90515-0
- LaFleur TL, Hossain A, Salis HM (2022). Automated model-predictive design of synthetic promoters to control transcriptional profiles in bacteria. *Nat Commun* 13:5159. doi:10.1038/s41467-022-32829-5
- Li H (2023). Protein-to-genome alignment with miniprot. *Bioinformatics* 39:btad014.
- Lorenz R, et al. (2011). ViennaRNA Package 2.0. *Algorithms Mol Biol* 6:26.
- Lundberg SM, Erion G, Chen H, et al. (2020). From local explanations to global understanding with explainable AI for trees. *Nat Mach Intell* 2:56–67.
- Lundberg SM, Lee S-I (2017). A unified approach to interpreting model predictions. *Adv Neural Inf Process Syst* 30.
- Moran PAP (1950). Notes on continuous stochastic phenomena. *Biometrika* 37:17–23. doi:10.1093/biomet/37.1-2.17
- Phipson B, Smyth GK (2010). Permutation p-values should never be zero: calculating exact p-values when permutations are randomly drawn. *Stat Appl Genet Mol Biol* 9:39.
- Politis DN, Romano JP (1992). A circular block-resampling procedure for stationary data. In *Exploring the Limits of Bootstrap*, Wiley, 263–270.
- Reiter L, et al. (2009). Protein identification false discovery rates for very large proteomics data sets generated by tandem mass spectrometry. *Mol Cell Proteomics* 8:2405–2417. doi:10.1074/mcp.M900317-MCP200
- Roots CT, Lukasiewicz AJ, Barrick JE (2021). OSTIR: open source translation initiation rate prediction. *J Open Source Softw* 6(64):3362. doi:10.21105/joss.03362
- Salgado H, et al. (2024). RegulonDB v12.0: a comprehensive resource of transcriptional regulation in *E. coli* K-12. *Nucleic Acids Res* 52(D1):D255–D264. doi:10.1093/nar/gkad1072 (version 14.5.0 of the database was used)
- Sanders J, Wen B, Rudnick PA, et al. (2025). A transformer model for de novo sequencing of data-independent acquisition mass spectrometry data. *Nat Methods* 22:1447–1453. doi:10.1038/s41592-025-02718-y
- Savitski MM, Wilhelm M, Hahne H, Kuster B, Bantscheff M (2015). A scalable approach for protein false discovery rate estimation in large proteomic data sets. *Mol Cell Proteomics* 14:2394–2404. doi:10.1074/mcp.M114.046995
- Schwanhäusser B, et al. (2011). Global quantification of mammalian gene expression control. *Nature* 473:337–342. doi:10.1038/nature10098
- Schwengers O, et al. (2021). Bakta: rapid and standardized annotation of bacterial genomes via alignment-free sequence identification. *Microb Genom* 7:000685. (PMC8743544)
- Sharp PM, Li W-H (1987). The codon adaptation index — a measure of directional synonymous codon usage bias, and its potential applications. *Nucleic Acids Res* 15:1281–1295. doi:10.1093/nar/15.3.1281
- Silva JC, et al. (2006). Absolute quantification of proteins by LCMS^E: a virtue of parallel MS acquisition. *Mol Cell Proteomics* 5:144–156. doi:10.1074/mcp.M500230-MCP200
- Su H, Zhang R, Söding J (2024). UniOP: a universal operon prediction for high-throughput prokaryotic (meta-)genomic data using intergenic distance. *bioRxiv*. doi:10.1101/2024.11.11.623000
- Subramanian A, et al. (2005). Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. *Proc Natl Acad Sci USA* 102:15545–15550. doi:10.1073/pnas.0506580102
- Taheri Ghahfarokhi SMA, Peña-Castillo L (2025). BacTermFinder: a comprehensive and general bacterial terminator finder using a CNN ensemble. *NAR Genom Bioinform* 7(1):lqaf016. doi:10.1093/nargab/lqaf016
- Tjaden B (2026). TerminatorNet: comprehensive identification of intrinsic transcription terminators in bacteria. *Bioinformatics* 42(3):btag116. doi:10.1093/bioinformatics/btag116
- UniProt Consortium (2025). UniProt: the Universal Protein Knowledgebase in 2025. *Nucleic Acids Res* 53(D1):D609–D617. doi:10.1093/nar/gkae1010
- Wen B, Freestone J, Riffle M, MacCoss MJ, Noble WS, Keich U (2025). Assessment of false discovery rate control in tandem mass spectrometry analysis using entrapment. *Nat Methods* 22:1454–1463. doi:10.1038/s41592-025-02719-x
- Westfall PH, Young SS (1993). *Resampling-Based Multiple Testing.* Wiley.
- Wright F (1990). The 'effective number of codons' used in a gene. *Gene* 87:23–29. doi:10.1016/0378-1119(90)90491-9
- Young MD, Wakefield MJ, Smyth GK, Oshlack A (2010). Gene ontology analysis for RNA-seq: accounting for selection bias. *Genome Biol* 11:R14. doi:10.1186/gb-2010-11-2-r14
- Zou H, Hastie T (2005). Regularization and variable selection via the elastic net. *J R Stat Soc B* 67:301–320.

<!-- src: GIPPA_SPEC.md:2760-2819 (references); docs/methods/*.md (Phipson-Smyth, Politis-Romano, Huggins, Blom, Westfall-Young, Benjamini-Yekutieli, Barber-Candès, He et al.); digest_C §10 (implemented but not cited in repo: Zou-Hastie, Ke et al., Lundberg-Lee, Bland-Altman — added by the writer); audit_R2_methods_C.md reference verdicts (Crossref), Cascadia/BacTermFinder/antiSMASH 8 re-checked by web search -->
