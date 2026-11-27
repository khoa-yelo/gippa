# PXD079213: Phase 2 preparation record

Aggregates only. Data under `$GIPPA_ROOT/campaign/PXD079213/` (not committed). No GIPPA run yet.

## Raw files

Selection: E233S control strain, 3 replicates (`run_selection.tsv`). Source: PRIDE over HTTPS (ADR-0077), resumable transfer, byte size checked against the PRIDE file list, deposit SHA-1 checked where the deposit has `checksum.txt`, sha256 recorded (also in `files.tsv`, column `sha256`). Fetch job 2001986 (one retry: 2003431, tasks 5, 8, 11, 12, 13), array task(s) 13.

| File | Bytes | Deposit checksum | Verified | sha256 |
|---|---|---|---|---|
| cont_1.raw | 1580759346 | none (no deposit checksum) | size | d01cbd6c227e176ab3d6f95e75d58f484f9176af8177ce9938a6f0c5bddf15d8 |
| cont_2.raw | 1560198727 | none (no deposit checksum) | size | 4feb3514c01e70ca045c60ed71a51258b63deb273dee9d7246ef9938d8142c5f |
| cont_3.raw | 1564899673 | none (no deposit checksum) | size | 2a52fed961d990eac9db77dd8a732c06b68328a0744e9e4da95b03bb510f4bed |

Fetched 3/3 files, 4.7 GB.

## Genome and Bakta annotation

NCBI RefSeq assembly over HTTPS, MD5 checked against `md5checksums.txt`; replicons NCBI nuccore marks circular get `[topology=circular]`; Bakta 1.12.2, DB 6.0 full, pinned image and `gippa.annotate.bakta_cmd.bakta_command` argv (`--keep-contig-headers`) at the campaign pin, as Phase 1b job 1998570. Bakta job 2001987 (array, one task per genome).

### sisl_rey15a: GCF_000189555.1 (ASM18955v1, Complete Genome, strain REY15A)

- Replicons (accession:bp:NCBI topology): NC_017276.1:2522992:circular; tagged circular: NC_017276.1
- `GCF_000189555.1_ASM18955v1_genomic.fna.gz`: MD5 555f1d7e127fd5bf3db0506b7f30dd45 (= md5checksums.txt, checked in job), sha256 7e6fad1b125901f5e73c3228afd2b9b48779c700dce6de238ec7f8193f819c35
- `GCF_000189555.1_ASM18955v1_genomic.gff.gz`: MD5 50d0f20230b6f77724096c2ac6abe1b2 (= md5checksums.txt, checked in job), sha256 c60038dee2d6166ee0ab6a34386ce56996ef6516ee5ae6ac4456a0d6f28be2ad
- Tagged genome `sisl_rey15a.fna.gz`: sha256 5f5698b621c3e08a90c7c81c38e6e1dff6184ec93e8f1dc9c40defc5dd957369
- Bakta runtime: wall_s=1033 cpus=16 node=wheat-07 job=2001987_16
- Versions: bakta 1.12.2; image sha256:1e24ae2f435c46f876a322298ab96d994d15aadc2824e1f16a77fa9a1bc07f18
- Genome: 2522992 bp, 1 replicon(s), GC 35.3 %, coding density 89.3 %
- Bakta counts: CDSs 2808, sORFs 0, pseudogenes 1, hypotheticals 2664, tRNAs 46, rRNAs 3, ncRNAs 1, CRISPR arrays 2
- RefSeq annotation (answer key, not an input): 2578 CDS, 191 pseudogenes; Bakta CDS + sORF / RefSeq CDS = 1.089

## Notes

- iProX (download.iprox.cn, ADR-0077 decision 1); first pass refused by the fetch script's host allowlist (fixed), fetched in job 2003431. No deposit checksum: verified by byte size only. Genome is REY15A (parent of E233S, which lacks pyrEF and lacS). 2,664 of 2,808 Bakta CDS are hypothetical; the archaeal check (CDS ratio 1.089; stop match) must pass before any run.
