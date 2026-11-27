# PXD072227: Phase 2 preparation record

Aggregates only. Data under `$GIPPA_ROOT/campaign/PXD072227/` (not committed). No GIPPA run yet.

## Raw files

Selection: CC9311, one nutrient-stress treatment (N or P not labelled), 5 replicates (`run_selection.tsv`). Source: PRIDE over HTTPS (ADR-0077), resumable transfer, byte size checked against the PRIDE file list, deposit SHA-1 checked where the deposit has `checksum.txt`, sha256 recorded (also in `files.tsv`, column `sha256`). Fetch job 2001986 (one retry: 2003431, tasks 5, 8, 11, 12, 13), array task(s) 9.

| File | Bytes | Deposit checksum | Verified | sha256 |
|---|---|---|---|---|
| 231113_Garcia_DIA_2ug_23_9311_T1.raw | 1149210140 | none (no deposit checksum) | size | 84b33b766d3326a611a87b01a55d7709c8210a3d299016f378ac71f1d87aaa30 |
| 231113_Garcia_DIA_2ug_24_9311_T2.raw | 1230595300 | none (no deposit checksum) | size | 06ab67a00b70b574c9e38e5738771bb57cab6bd8db90b79bbf8da7ccb6ba60b7 |
| 231113_Garcia_DIA_2ug_25_9311_T3.raw | 1237520860 | none (no deposit checksum) | size | d46f7bd814d86245805c312aada2502e31ca9c0580be5f72c56df02c6b26ab2e |
| 231113_Garcia_DIA_2ug_26_9311_T4.raw | 1203011159 | none (no deposit checksum) | size | 1491cc26ba46e61c3d2e0d7a34c5ac3dfa44a00ba67a38c997c998fa13c552bc |
| 231113_Garcia_DIA_2ug_27_9311_T5.raw | 1207038226 | none (no deposit checksum) | size | e9d804a2be2da8c951289cd8106374df71c8502f722666f75b74a99447870162 |

Fetched 5/5 files, 6.0 GB.

## Genome and Bakta annotation

NCBI RefSeq assembly over HTTPS, MD5 checked against `md5checksums.txt`; replicons NCBI nuccore marks circular get `[topology=circular]`; Bakta 1.12.2, DB 6.0 full, pinned image and `gippa.annotate.bakta_cmd.bakta_command` argv (`--keep-contig-headers`) at the campaign pin, as Phase 1b job 1998570. Bakta job 2001987 (array, one task per genome).

### syn_cc9311: GCF_000014585.1 (ASM1458v1, Complete Genome, strain CC9311)

- Replicons (accession:bp:NCBI topology): NC_008319.1:2606748:circular; tagged circular: NC_008319.1
- `GCF_000014585.1_ASM1458v1_genomic.fna.gz`: MD5 860bb828b62d8fe98ed1fc58b3da0295 (= md5checksums.txt, checked in job), sha256 13e836a9752f81702a50b3ffc2952997e8179cd445dff294010f5b259e4a2970
- `GCF_000014585.1_ASM1458v1_genomic.gff.gz`: MD5 1091d39352b784900aba15b0ce397b18 (= md5checksums.txt, checked in job), sha256 9734e637576558b69bd33cd6871d2ef1b7d2a7547746a0791cf0493a0235e4f2
- Tagged genome `syn_cc9311.fna.gz`: sha256 675d870a5605c3a2658e71dc8118679b58101d034b139a79ef0654f9eaf0a726
- Bakta runtime: wall_s=606 cpus=16 node=wheat-07 job=2001987_12
- Versions: bakta 1.12.2; image sha256:1e24ae2f435c46f876a322298ab96d994d15aadc2824e1f16a77fa9a1bc07f18
- Genome: 2606748 bp, 1 replicon(s), GC 52.4 %, coding density 88.1 %
- Bakta counts: CDSs 2865, sORFs 0, pseudogenes 6, hypotheticals 477, tRNAs 44, rRNAs 6, ncRNAs 17, CRISPR arrays 0
- RefSeq annotation (answer key, not an input): 2699 CDS, 19 pseudogenes; Bakta CDS + sORF / RefSeq CDS = 1.062

## Notes

- CC9311 runs 23-27 (Khoa, question 21); the treatment (N or P stress) is not labelled in the deposit.
