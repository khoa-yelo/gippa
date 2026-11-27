# PXD066008: Phase 2 preparation record

Aggregates only. Data under `$GIPPA_ROOT/campaign/PXD066008/` (not committed). No GIPPA run yet.

## Raw files

Selection: H26 wild type, glucose synthetic medium, 5 runs (2 timsTOF Pro II, 3 timsTOF HT) (`run_selection.tsv`). Source: PRIDE over HTTPS (ADR-0077), resumable transfer, byte size checked against the PRIDE file list, deposit SHA-1 checked where the deposit has `checksum.txt`, sha256 recorded (also in `files.tsv`, column `sha256`). Fetch job 2001986 (one retry: 2003431, tasks 5, 8, 11, 12, 13), array task(s) 12.

| File | Bytes | Deposit checksum | Verified | sha256 |
|---|---|---|---|---|
| 230216_2305_AgSoppa_TimsDIA_PepSep15_17min_10x_2uL_13_Slot2-13_1_237.d.tar | 2960263680 | sha1:da219b5e404fe32029ad25de6c14d59ee061d671 | size + SHA-1 | cfd6c27d236245002696f222bc06e9acb58d3d96126e8c494c40f14f716443ca |
| 230216_2305_AgSoppa_TimsDIA_PepSep15_17min_10x_2uL_14_Slot2-14_1_239.d.tar | 3040278528 | sha1:1ba7ae913b4781ea74b25b819504a42bc3424b40 | size + SHA-1 | d6d128ab2f15c71a45908f35caab0114cacabfd65f3611fceb5faed81228d47f |
| 240118_2453_JoSoppa_DIA_Sample-13_Slot2-13_1_613.d.tar | 2493631488 | sha1:a99d0a4b93b304a944120943b2b8baa559b095f1 | size + SHA-1 | 11c217de55c9290c0d5c7a3fae249e28b328b95142f63ac969430fc5f261357e |
| 240118_2453_JoSoppa_DIA_Sample-14_Slot2-14_1_616.d.tar | 2598699520 | sha1:dd5e39d47774b2db76631eb688b621a2948a1144 | size + SHA-1 | 6e8668f8132b8f0aefec6e81e1b8ef58c610b52350dcb3d938e8dec5e3f2e574 |
| 240820_2549_JoSoppa_1_Slot2-19_1_2535.d.tar | 2508644352 | sha1:ce8bbf8b814b5742a6362eca9f7150f3eaef818d | size + SHA-1 | feaac018e117071e051c3e0066e5627f12488c28522c4e1b3773921b93b38cb2 |

Fetched 5/5 files, 13.6 GB.

## Genome and Bakta annotation

NCBI RefSeq assembly over HTTPS, MD5 checked against `md5checksums.txt`; replicons NCBI nuccore marks circular get `[topology=circular]`; Bakta 1.12.2, DB 6.0 full, pinned image and `gippa.annotate.bakta_cmd.bakta_command` argv (`--keep-contig-headers`) at the campaign pin, as Phase 1b job 1998570. Bakta job 2001987 (array, one task per genome).

### hvol_ds2: GCF_000025685.1 (ASM2568v1, Complete Genome, strain DS2)

- Replicons (accession:bp:NCBI topology): NC_013967.1:2847757:circular,NC_013968.1:85092:circular,NC_013965.1:6359:circular,NC_013964.1:437906:circular,NC_013966.1:635786:circular; tagged circular: NC_013967.1,NC_013968.1,NC_013965.1,NC_013964.1,NC_013966.1
- `GCF_000025685.1_ASM2568v1_genomic.fna.gz`: MD5 4197e319edfb2f801c2f1375b19392b8 (= md5checksums.txt, checked in job), sha256 7bb37b065dfeeea1802b24a08cdae400f08203d821a89896282dcb3adcb9ef91
- `GCF_000025685.1_ASM2568v1_genomic.gff.gz`: MD5 476dedf9933844600d8f090a8b413e02 (= md5checksums.txt, checked in job), sha256 c9f8b0fa5910361e0fe0bb045dd4a3e92358600800d14f5908546fb171ce5474
- Tagged genome `hvol_ds2.fna.gz`: sha256 94be40697a6f799d4f66c8ead03cfaeec62f78e065b60c7e8bddbcb71431fc56
- Bakta runtime: wall_s=1408 cpus=16 node=wheat-04 job=2001987_15
- Versions: bakta 1.12.2; image sha256:1e24ae2f435c46f876a322298ab96d994d15aadc2824e1f16a77fa9a1bc07f18
- Genome: 4012900 bp, 5 replicon(s), GC 65.5 %, coding density 86.3 %
- Bakta counts: CDSs 4023, sORFs 0, pseudogenes 0, hypotheticals 3488, tRNAs 53, rRNAs 6, ncRNAs 3, CRISPR arrays 4
- RefSeq annotation (answer key, not an input): 3897 CDS, 120 pseudogenes; Bakta CDS + sORF / RefSeq CDS = 1.032

## Notes

- Retry: first-pass resumed transfers produced a full-length file with a SHA-1 mismatch, discarded and fetched again once in job 2003431 (verified): 240118_2453_JoSoppa_DIA_Sample-13_Slot2-13_1_613.d.tar. Genome is DS2 (parent of H26) per the plan; H26 lacks pHV2 (NC_013965.1, 6,359 bp) and pyrE2. Bakta does not officially support archaea: 3,488 of 4,023 CDS are hypothetical; the plan's archaeal check (CDS count within 10 %: 1.032 here; ≥ 90 % of RefSeq CDS with the same stop) must pass before any run.
