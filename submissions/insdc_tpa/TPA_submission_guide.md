# INSDC / GenBank Third-Party Annotation (TPA) Submission Guide

This directory contains submission-ready files for the two DTR-verified, complete circularized bacteriophage genomes recovered from ENA freshwater metagenome assembly `ERZ18661299` (Lake Mekkojarvi, BioProject `PRJEB47226`, BioSample `SAMEA6872951`, raw paired-end reads `ERR4194644`).

## Files Included

- `CTG_000273_tpa.fsa`, `CTG_000273_tpa.tbl`, `CTG_000273_tpa.gff3`, `CTG_000273_tpa.gbk` (`37,043 bp` trimmed circularized genome; `34 bp` DTR removed; `47` complete CDSs, with linear `_48` and `_1` fused in-frame across the circular seam into a single 162-aa DNA repair protein)
- `CTG_000311_tpa.fsa`, `CTG_000311_tpa.tbl`, `CTG_000311_tpa.gff3`, `CTG_000311_tpa.gbk` (`34,122 bp` trimmed circularized genome; `55 bp` DTR removed; `51` CDSs, with `_51` wrapping across the circular seam)

## Primary Sequence & Read Evidence Citations

- **Primary ENA Assembly:** `ERZ18661299` (`ERZ18661299.273` and `ERZ18661299.311`)
- **Primary Raw Reads:** `ERR4194644` (`ERR4194644_1.fastq.gz`, `ERR4194644_2.fastq.gz`)
- **TPA Classification:** `TPA:experimental` / `TPA:inferential` (assemblies derived from `ERZ18661299`, circularity and DTR-seam continuity experimentally verified by mapping `ERR4194644` reads with `minimap2 -ax sr`, `MAPQ >= 20`, `NM <= 3`: `CTG_000273` at `29.72x` mean depth with `116` full-DTR reads; `CTG_000311` at `15.45x` mean depth with `3` full-DTR reads).

## Submitting via NCBI `table2asn` or BankIt

1. Create a submission template (`template.sbt`) at https://submit.ncbi.nlm.nih.gov/genbank/template/submission/ with your author/affiliation details.
2. Bundle each `.fsa` + `.tbl` pair using NCBI `table2asn` (or upload directly via BankIt / NCBI Submission Portal under **Third-Party Annotation**):
   ```bash
   table2asn -t template.sbt -i CTG_000273_tpa.fsa -f CTG_000273_tpa.tbl -o CTG_000273_tpa.sqn -Z -V b
   table2asn -t template.sbt -i CTG_000311_tpa.fsa -f CTG_000311_tpa.tbl -o CTG_000311_tpa.sqn -Z -V b
   ```
3. Reference `ERZ18661299` and `ERR4194644` in the TPA primary-accession field when prompted by BankIt or ENA Webin.
