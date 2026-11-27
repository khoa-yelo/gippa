# PXD060688: Phase 2 preparation record

Aggregates only. Data under `$GIPPA_ROOT/campaign/PXD060688/` (not committed). No GIPPA run yet.

## Raw files

Selection: DK1622 wild type, untreated, exponential (presumed), 4 replicates (`run_selection.tsv`). Source: PRIDE over HTTPS (ADR-0077), resumable transfer, byte size checked against the PRIDE file list, deposit SHA-1 checked where the deposit has `checksum.txt`, sha256 recorded (also in `files.tsv`, column `sha256`). Fetch job 2001986 (one retry: 2003431, tasks 5, 8, 11, 12, 13), array task(s) 11.

| File | Bytes | Deposit checksum | Verified | sha256 |
|---|---|---|---|---|
| DdiA-4-01.raw | 3240703211 | sha1:c204c1ddceeb2ea03f73758df10f9cd96983ddda | size + SHA-1 | 9c38f1f7a2880cb55427718e30a5c22988ecfad0b7d6ec857886c53aafeed6ed |
| DdiA-4-02.raw | 3323736144 | sha1:bffa12be741d3104deaac7a0bce6ffc07f622f4b | size + SHA-1 | b8790ccbc6abeb739a25ee1368b063748c3dfee82c201ae85d9c437f0c0ea327 |
| DdiA-4-03.raw | 3303039144 | sha1:4418a21200bb9b708fa9e3d35af059b9f41af00d | size + SHA-1 | bf3fff8ff77df2179013bf84ab66c7ad40a003abc4bc437b1e5b4f81ee22e367 |
| DdiA-4-04.raw | 3248336183 | sha1:89a6eeb72ce6baf0b69466d24d6d885a3e0d8771 | size + SHA-1 | a1bccf92cf16c8c6ff9b5c285d466e3c8d7ae9af308ea0428a58c1440e5abbf4 |

Fetched 4/4 files, 13.1 GB.

## Genome and Bakta annotation

NCBI RefSeq assembly over HTTPS, MD5 checked against `md5checksums.txt`; replicons NCBI nuccore marks circular get `[topology=circular]`; Bakta 1.12.2, DB 6.0 full, pinned image and `gippa.annotate.bakta_cmd.bakta_command` argv (`--keep-contig-headers`) at the campaign pin, as Phase 1b job 1998570. Bakta job 2001987 (array, one task per genome).

### mxan_dk1622: GCF_000012685.1 (ASM1268v1, Complete Genome, strain DK 1622)

- Replicons (accession:bp:NCBI topology): NC_008095.1:9139763:circular; tagged circular: NC_008095.1
- `GCF_000012685.1_ASM1268v1_genomic.fna.gz`: MD5 21612cf40348ac872b0d865b11f8df01 (= md5checksums.txt, checked in job), sha256 073ba0271fcb5142a367af4c9af97812086677c6cb7954cf01d0850b90fd9934
- `GCF_000012685.1_ASM1268v1_genomic.gff.gz`: MD5 40a1cf6f0586531525becf46c91b1034 (= md5checksums.txt, checked in job), sha256 d748f3face59025c0444cdc9133d904b0f69d5184734f1e72f5298cfb3431374
- Tagged genome `mxan_dk1622.fna.gz`: sha256 a36130a65c51933c5b18f1f44673773c1df710003dfb57268be853130e05b218
- Bakta runtime: wall_s=879 cpus=16 node=wheat-07 job=2001987_14
- Versions: bakta 1.12.2; image sha256:1e24ae2f435c46f876a322298ab96d994d15aadc2824e1f16a77fa9a1bc07f18
- Genome: 9139763 bp, 1 replicon(s), GC 68.9 %, coding density 90.8 %
- Bakta counts: CDSs 7326, sORFs 0, pseudogenes 8, hypotheticals 391, tRNAs 69, rRNAs 12, ncRNAs 6, CRISPR arrays 5
- RefSeq annotation (answer key, not an input): 7234 CDS, 72 pseudogenes; Bakta CDS + sORF / RefSeq CDS = 1.013

## Notes

- Retry: first-pass resumed transfers produced a full-length file with a SHA-1 mismatch, discarded and fetched again once in job 2003431 (verified): DdiA-4-03.raw.
