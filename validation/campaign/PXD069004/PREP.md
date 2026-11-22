# PXD069004: Phase 2 preparation record

Aggregates only. Data under `$GIPPA_ROOT/campaign/PXD069004/` (not committed). No GIPPA run yet.

## Raw files

Selection: NA1000 wild type, M2G (no GlcNAc) (`run_selection.tsv`). Source: PRIDE over HTTPS (ADR-0077), resumable transfer, byte size checked against the PRIDE file list, deposit SHA-1 checked where the deposit has `checksum.txt`, sha256 recorded (also in `files.tsv`, column `sha256`). Fetch job 2001986 (one retry: 2003431, tasks 5, 8, 11, 12, 13), array task(s) 6.

| File | Bytes | Deposit checksum | Verified | sha256 |
|---|---|---|---|---|
| AmiX-5-4.raw | 1851297444 | sha1:28600af86703e26a1b03d7bf23c4c0e90d95cadf | size + SHA-1 | 5dd1c5ef21226801322cec1b5625619424f902a6e44dce4451dfff8c5bbe6fe0 |
| AmiX-5-5.raw | 1732249783 | sha1:1ba341e4ed162b92fa390bb3cf66c12ecb1ce0bb | size + SHA-1 | 2c67e745e50763a85dfb2f12db302b595e1eb61832ebd83835629be503a72bb5 |
| AmiX-5-6.raw | 1735114935 | sha1:552548644e8ce0c6196ccba0016d7cc03d12f80f | size + SHA-1 | 8786a252d1c7343d7f602cf62179d69b7ed79756864c1170fb5b1098d58d7235 |

Fetched 3/3 files, 5.3 GB.

## Genome and Bakta annotation

NCBI RefSeq assembly over HTTPS, MD5 checked against `md5checksums.txt`; replicons NCBI nuccore marks circular get `[topology=circular]`; Bakta 1.12.2, DB 6.0 full, pinned image and `gippa.annotate.bakta_cmd.bakta_command` argv (`--keep-contig-headers`) at the campaign pin, as Phase 1b job 1998570. Bakta job 2001987 (array, one task per genome).

### ccre_na1000: GCF_000022005.1 (ASM2200v1, Complete Genome, strain NA1000)

- Replicons (accession:bp:NCBI topology): NC_011916.1:4042929:circular; tagged circular: NC_011916.1
- `GCF_000022005.1_ASM2200v1_genomic.fna.gz`: MD5 9362062e72a914fd10324a2f22326ea7 (= md5checksums.txt, checked in job), sha256 11d097d0cd1ca09dd126c92d563c4fb5e5e4692d979e76362325118b44ac5d9b
- `GCF_000022005.1_ASM2200v1_genomic.gff.gz`: MD5 2364b9578d6eba2fed2ffec1ef82244f (= md5checksums.txt, checked in job), sha256 d896d467681e528aa5f366d1c608c8a6ac00c2ab31f287c9b0106f3c33ba098b
- Tagged genome `ccre_na1000.fna.gz`: sha256 8c325e6358566d68947148db2294baa22be6f438db71a4fa90749479079d655f
- Bakta runtime: wall_s=696 cpus=16 node=wheat-07 job=2001987_9
- Versions: bakta 1.12.2; image sha256:1e24ae2f435c46f876a322298ab96d994d15aadc2824e1f16a77fa9a1bc07f18
- Genome: 4042929 bp, 1 replicon(s), GC 67.2 %, coding density 89.5 %
- Bakta counts: CDSs 3722, sORFs 2, pseudogenes 1, hypotheticals 159, tRNAs 51, rRNAs 6, ncRNAs 18, CRISPR arrays 0
- RefSeq annotation (answer key, not an input): 3886 CDS, 0 pseudogenes; Bakta CDS + sORF / RefSeq CDS = 0.958

## Notes

- NA1000 wild type in M2G. RefSeq has 3,886 CDS; Bakta has 3,722 + 2 sORFs (ratio 0.958).
