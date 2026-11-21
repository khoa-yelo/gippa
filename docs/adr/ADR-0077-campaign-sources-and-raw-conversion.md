# ADR-0077: Download sources and Thermo raw conversion for the demo and discovery campaign
- Status: accepted-by-khoa (2026-11-19, in session: "Yes")
- Date: 2026-11-19
- Spec sections: §3.2 (allowed sources), R6 (official sources only), R11 (large files under
  `$GIPPA_ROOT`), §15 (pilot rules, reused)
- Needs Khoa: yes (adds download sources and a tool)

## Context
The demo and discovery campaign (handover plan of 2026-11-17) runs GIPPA on public DIA data:
E. coli K-12 MG1655 (PXD058422, PXD066241), NCM3722 (PXD060434) and ten further datasets.
Spec §3.2 allows ProteomeXchange and iProX only for the raw files of PXD078459. Most campaign
datasets are deposited at PRIDE. Eight of the twelve datasets are Thermo `.raw`; the pilot only ran
Bruker `.d`, and Cascadia (arm E4) needs mzML. NCBI (genomes, RefSeq GFF) and UniProt (host
proteomes, e.g. rabbit for T. pallidum) are already allowed by §3.2.

## Decision
1. **New sources**, for the campaign only (`validation/campaign/`, data under
   `$GIPPA_ROOT/campaign/`), on the same terms as the pilot (TLS verified, sha256 recorded in a
   per-dataset `files.tsv`, one retry, then record the block):
   - ProteomeXchange (`proteomecentral.proteomexchange.org`) for **any** campaign accession, for the
     dataset XML and file lists;
   - PRIDE (`www.ebi.ac.uk/pride/ws/archive/` for metadata and checksums, `ftp.pride.ebi.ac.uk`
     over HTTPS or FTP for files);
   - iProX (`www.iprox.cn` and the FTP host the dataset page names) for any campaign accession
     deposited there;
   - PubMed Central (`www.ncbi.nlm.nih.gov/pmc/`, `pmc.ncbi.nlm.nih.gov`) for the papers'
     supplementary tables that hold the answer keys (read, never redistributed).
2. **ThermoRawFileParser** (Apache-2.0, bioconda `thermorawfileparser`) is added as campaign
   tooling, run as a separate program in its own micromamba environment
   (`validation/campaign/envs/thermorawfileparser.yml`, version pinned there), only if DIA-NN 2.3.1
   on Linux cannot read `.raw` directly (Phase 0 check) and for Cascadia's mzML input. It is not a
   GIPPA dependency: nothing under `src/` or `workflow/` changes.
3. **Rules unchanged.** Aggregates only in git (no peptide sequences, no DIA-NN rows); downloads
   stay under `$GIPPA_ROOT/campaign/`; about 500 GB total before asking Khoa.

## Evidence (before → after: metrics, tests, simulations)
None yet; the Phase 0 inventory and the `.raw` test record the outcome.

## Consequences
- The Phase 0 inventory (`validation/campaign/datasets.tsv`, `<id>/files.tsv`) can start once this
  is accepted.
- No change to the GIPPA package, its tests or its licences.
- Each campaign report cites its accession and paper; raw files are never redistributed.

## Amendment 1 (2026-11-19, accepted-by-khoa in session: "Add both")
The Phase 0 inventory (`validation/campaign/`, PR #603) found two more sources:
- **Panorama Public** (`panoramaweb.org`, its PXD-linked project folders) for PXD066241, which
  holds the E3 knockout runs (ΔmntS, Δlon, ΔclpP) of the K-12 demonstration;
- **NCBI E-utilities** (`eutils.ncbi.nlm.nih.gov`) for the full text and supplementary-file
  links of the papers, because the PMC web pages answer automated clients with a reCAPTCHA.
Same terms as decision 1: TLS verified, sha256 in `files.tsv`, one retry, then record the block.

The Phase 0 check also found that the pinned DIA-NN 2.3.1 image reads Thermo `.raw` directly on
Linux (job 1998427, `validation/campaign/raw_check.md`), so ThermoRawFileParser is needed only for
Cascadia's mzML input (arm E4).
