# PXD057335: Phase 2 preparation record

Aggregates only. Data under `$GIPPA_ROOT/campaign/PXD057335/` (not committed). No GIPPA run yet.

## Raw files

Selection: H37Rv, DMSO vehicle, 48 h 37 °C (`run_selection.tsv`). Source: PRIDE over HTTPS (ADR-0077), resumable transfer, byte size checked against the PRIDE file list, deposit SHA-1 checked where the deposit has `checksum.txt`, sha256 recorded (also in `files.tsv`, column `sha256`). Fetch job 2001986 (one retry: 2003431, tasks 5, 8, 11, 12, 13), array task(s) 10.

| File | Bytes | Deposit checksum | Verified | sha256 |
|---|---|---|---|---|
| 20230614_Exploris_IB_EXP57_DMSO_01.raw | 2986095367 | none (no deposit checksum) | size | 9324e626ab3c9e693cbd33b4e04795a2a12c72b9ba829762ebce9d0aa0d9de14 |
| 20230614_Exploris_IB_EXP57_DMSO_02.raw | 3081515553 | none (no deposit checksum) | size | 0b13121d2b755669d771c5bf7e9e64d2cc513569ddf317115966593b7479d11a |
| 20230614_Exploris_IB_EXP57_DMSO_03.raw | 3113059923 | none (no deposit checksum) | size | 4880cd9c1e3a824f639313037174391fd88e590980be54163f73ebdb0741ba5a |

Fetched 3/3 files, 9.2 GB.

## Genome and Bakta annotation

NCBI RefSeq assembly over HTTPS, MD5 checked against `md5checksums.txt`; replicons NCBI nuccore marks circular get `[topology=circular]`; Bakta 1.12.2, DB 6.0 full, pinned image and `gippa.annotate.bakta_cmd.bakta_command` argv (`--keep-contig-headers`) at the campaign pin, as Phase 1b job 1998570. Bakta job 2001987 (array, one task per genome).

### mtb_h37rv: GCF_000195955.2 (ASM19595v2, Complete Genome, strain H37Rv)

- Replicons (accession:bp:NCBI topology): NC_000962.3:4411532:linear; tagged circular: none
- `GCF_000195955.2_ASM19595v2_genomic.fna.gz`: MD5 f93726a0a6e57af2bf1173998fbeb546 (= md5checksums.txt, checked in job), sha256 35911742cea55cd556d54948f56e3efa757e40d44e18ebc8742cacc59b2a976e
- `GCF_000195955.2_ASM19595v2_genomic.gff.gz`: MD5 4d53755c7968f2509d772b12ba2c4136 (= md5checksums.txt, checked in job), sha256 fbd1a13d212e1d042059de7e44dfe5d7dbd96bcc71a53f2e5040ea0e50a4a2ed
- Tagged genome `mtb_h37rv.fna.gz`: sha256 b3bd16d538851d164443acae9ca3bbc1bf59869edaf3d7f2d9183860e37492c8
- Bakta runtime: wall_s=886 cpus=16 node=wheat-04 job=2001987_13
- Versions: bakta 1.12.2; image sha256:1e24ae2f435c46f876a322298ab96d994d15aadc2824e1f16a77fa9a1bc07f18
- Genome: 4411532 bp, 1 replicon(s), GC 65.6 %, coding density 91.0 %
- Bakta counts: CDSs 4078, sORFs 2, pseudogenes 5, hypotheticals 137, tRNAs 45, rRNAs 3, ncRNAs 21, CRISPR arrays 2
- RefSeq annotation (answer key, not an input): 3906 CDS, 30 pseudogenes; Bakta CDS + sORF / RefSeq CDS = 1.045

## Notes

- NCBI nuccore gives NC_000962.3 as linear (LOCUS line too), so the replicon is left untagged per the recipe, although the chromosome is circular biologically. Revisit before the run if origin-spanning features matter.
