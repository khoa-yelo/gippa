# PXD049476: Phase 2 preparation record

Aggregates only. Data under `$GIPPA_ROOT/campaign/PXD049476/` (not committed). No GIPPA run yet.

## Raw files

Selection: parental (presumed), exponential (`run_selection.tsv`). Source: PRIDE over HTTPS (ADR-0077), resumable transfer, byte size checked against the PRIDE file list, deposit SHA-1 checked where the deposit has `checksum.txt`, sha256 recorded (also in `files.tsv`, column `sha256`). Fetch job 2001986 (one retry: 2003431, tasks 5, 8, 11, 12, 13), array task(s) 7.

| File | Bytes | Deposit checksum | Verified | sha256 |
|---|---|---|---|---|
| Vnat-1-1.raw | 2330482180 | sha1:7a78fbeb86d0fddb93092349f2a5e17a156c91d0 | size + SHA-1 | 6201de24dc10e2179590a71bf97ff11ddb2da4680cce186f120149ac565e2817 |
| Vnat-1-2.raw | 2391772933 | sha1:7366be064c9c951bbc8ffb8380618e0e7feabf53 | size + SHA-1 | d927f690a70ebb6c3bb1f7bbb6b0c26b85c9fa11fa0194d38f76cbc1102f6b18 |
| Vnat-1-3.raw | 2451875699 | sha1:68df2eec1b83fcf558fbd777d4e1637e9ce7ad65 | size + SHA-1 | 43ab7eadfa34a68daacb86ab244989a512199be75a04eba72bdab405de4d7923 |
| Vnat-1-4.raw | 2430403880 | sha1:6a40719e008e46082fc9903dcbc3108517e12c38 | size + SHA-1 | 235744f290e3fbe8310d4cfe0ce634f88d99794b353fbbf74459f73fdcd8f33e |

Fetched 4/4 files, 9.6 GB.

## Genome and Bakta annotation

NCBI RefSeq assembly over HTTPS, MD5 checked against `md5checksums.txt`; replicons NCBI nuccore marks circular get `[topology=circular]`; Bakta 1.12.2, DB 6.0 full, pinned image and `gippa.annotate.bakta_cmd.bakta_command` argv (`--keep-contig-headers`) at the campaign pin, as Phase 1b job 1998570. Bakta job 2001987 (array, one task per genome).

### vnat_atcc14048: GCF_001456255.1 (ASM145625v1, Complete Genome, strain ATCC 14048)

- Replicons (accession:bp:NCBI topology): NZ_CP009977.1:3248023:circular,NZ_CP009978.1:1927130:circular; tagged circular: NZ_CP009977.1,NZ_CP009978.1
- `GCF_001456255.1_ASM145625v1_genomic.fna.gz`: MD5 308331ba70b65ae5e3f3c15f65801a49 (= md5checksums.txt, checked in job), sha256 7384bd1c9412d7502d0bdb47950667f9262a62240c1816ea92b5d8d1d423068f
- `GCF_001456255.1_ASM145625v1_genomic.gff.gz`: MD5 8fa09cf175b8a3fe6ec6511258580cdd (= md5checksums.txt, checked in job), sha256 0a28822ad2ef6bb489601b7e6f61b82ea6ced3d26eda5fee9294b5a0c24328c2
- Tagged genome `vnat_atcc14048.fna.gz`: sha256 abffc6c36484a117c0646297a477452e3284bbae011f543d9c23a1f4377527cc
- Bakta runtime: wall_s=897 cpus=16 node=wheat-04 job=2001987_10
- Versions: bakta 1.12.2; image sha256:1e24ae2f435c46f876a322298ab96d994d15aadc2824e1f16a77fa9a1bc07f18
- Genome: 5175153 bp, 2 replicon(s), GC 45.1 %, coding density 87.3 %
- Bakta counts: CDSs 4553, sORFs 0, pseudogenes 2, hypotheticals 166, tRNAs 129, rRNAs 34, ncRNAs 23, CRISPR arrays 0
- RefSeq annotation (answer key, not an input): 4504 CDS, 19 pseudogenes; Bakta CDS + sORF / RefSeq CDS = 1.011

## Notes

- Vnat-1-1..4, assumed to be the parental strain (Khoa, question 18). Two chromosomes, both circular.
