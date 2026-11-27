# PXD065346: Phase 2 preparation record

Aggregates only. Data under `$GIPPA_ROOT/campaign/PXD065346/` (not committed). No GIPPA run yet.

## Raw files

Selection: 630Δerm total lysate + 0.5 pg/µg synthetic small-protein spike-in, 3 replicates (`run_selection.tsv`). Source: PRIDE over HTTPS (ADR-0077), resumable transfer, byte size checked against the PRIDE file list, deposit SHA-1 checked where the deposit has `checksum.txt`, sha256 recorded (also in `files.tsv`, column `sha256`). Fetch job 2001986 (one retry: 2003431, tasks 5, 8, 11, 12, 13), array task(s) 8.

| File | Bytes | Deposit checksum | Verified | sha256 |
|---|---|---|---|---|
| 221217_QE1_PM1_SK_JB_Cdiff_SpikeIn_05pg_ug_DIA_Ctrl_A.raw | 5379164712 | sha1:c1c0f7ecf7df5a646d7adbb0ab362ad64276d5b9 | size + SHA-1 | d727fda637da57de71db75250f5a374266f0a58b13b6bd9de08213708b5bb603 |
| 221219_QE1_PM1_SK_JB_Cdiff_SpikeIn_05pg_ug_DIA_Ctrl_B.raw | 5498035652 | sha1:ecf7972a5b756739dd94250889784acb7e52f81f | size + SHA-1 | c6036cd32cae72bc8459e745c2b85376dc790f9dba0b1515de9c2a4986c2058b |
| 230106_QE1_PM1_SK_JB_Cdiff_SpikeIn_05pg_ug_DIA_Ctrl_C.raw | 4309072956 | sha1:575f7e2b1217380dd7977a124182ac1f810669cf | size + SHA-1 | e523e3c8494d79feefdb96cdfe6871024417760064be2e99ada725267928506f |

Fetched 3/3 files, 15.2 GB.

## Genome and Bakta annotation

NCBI RefSeq assembly over HTTPS, MD5 checked against `md5checksums.txt`; replicons NCBI nuccore marks circular get `[topology=circular]`; Bakta 1.12.2, DB 6.0 full, pinned image and `gippa.annotate.bakta_cmd.bakta_command` argv (`--keep-contig-headers`) at the campaign pin, as Phase 1b job 1998570. Bakta job 2001987 (array, one task per genome).

### cdif_630derm: GCF_002080065.1 (ASM208006v1, Complete Genome, strain 630 delta erm)

- Replicons (accession:bp:NCBI topology): NZ_CP016318.1:4293048:circular,NZ_CP016319.1:7881:circular; tagged circular: NZ_CP016318.1,NZ_CP016319.1
- `GCF_002080065.1_ASM208006v1_genomic.fna.gz`: MD5 c19ed19053a12c3004fac62fe9b60fa6 (= md5checksums.txt, checked in job), sha256 abf952a92a35566453a7b8b011cbd532111d8172e87d691e49163d15da3d91bd
- `GCF_002080065.1_ASM208006v1_genomic.gff.gz`: MD5 efc7ee8023a8a58dd45be9fb4c33e93a (= md5checksums.txt, checked in job), sha256 be720b966227fe6adaee80ccffffe29feffceb4714abee8c6c8a080e504ead56
- Tagged genome `cdif_630derm.fna.gz`: sha256 49065d7c83c1d51f6f16acbafec642af68c91437fac5348a3858e0f4d323e95f
- Bakta runtime: wall_s=666 cpus=16 node=wheat-07 job=2001987_11
- Versions: bakta 1.12.2; image sha256:1e24ae2f435c46f876a322298ab96d994d15aadc2824e1f16a77fa9a1bc07f18
- Genome: 4300929 bp, 2 replicon(s), GC 29.1 %, coding density 84.7 %
- Bakta counts: CDSs 3820, sORFs 6, pseudogenes 2, hypotheticals 92, tRNAs 90, rRNAs 35, ncRNAs 9, CRISPR arrays 10
- RefSeq annotation (answer key, not an input): 3801 CDS, 49 pseudogenes; Bakta CDS + sORF / RefSeq CDS = 1.007

## Notes

- Retry: first-pass resumed transfers produced a full-length file with a SHA-1 mismatch, discarded and fetched again once in job 2003431 (verified): 230106_..._DIA_Ctrl_C.raw. Runs are the 0.5 pg/µg spike-in total lysate (Khoa, question 20): spiked synthetic small proteins are excluded from every claim.
