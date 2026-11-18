# ICTV Taxonomic Proposal Summary: New Genus *Mekkojarvivirus* and Species *Mekkojarvivirus mekkojarviense* (`CTG_000311`)

**Proposer:** Niall Gillen ([Institution / Email — update prior to submission])  
**Subcommittee:** Bacterial Viruses Subcommittee (`Caudoviricetes`)  
**Companion Spreadsheet:** `ictv_taxoprop_CTG_000311.xlsx`  
**Note on Taxon Names:** *Mekkojarvivirus* and *Mekkojarvivirus mekkojarviense* are placeholder Latinized names formatted per ICTV binomial rules; customize as desired before formal submission, and insert the assigned GenBank/ENA TPA accession once `submissions/insdc_tpa/CTG_000311_tpa.*` is accepted.

## 1. Taxonomic Summary

We propose creating one new genus (**Mekkojarvivirus**) and one new species (**Mekkojarvivirus mekkojarviense**) directly within class *Caudoviricetes* to classify `CTG_000311` (`34,122 bp` circularized / `34,177 bp` DTR-flanked, `37.44%` GC, `51` ORFs), a complete, read-verified freshwater bacteriophage genome recovered from Lake Mekkojarvi metagenome `ERZ18661299` (BioProject `PRJEB47226`, BioSample `SAMEA6872951`, raw reads `ERR4194644`).

## 2. Phylogenomic and Proteomic Demarcation

* **Nucleotide Divergence:** `CTG_000311` shows `0.0%` whole-genome nucleotide coverage against NCBI Viral RefSeq and surveyed freshwater viromes, well below the standard `95%` ANI / `85%` coverage species demarcation threshold.
* **Whole-Proteome Tree (ViPTree v4.2):** Across `7,097` prokaryotic dsDNA reference viruses (session `u6lWCt_eU2Wt7KD86pfmNw`), `CTG_000311` falls below ViPTree's $S_G = 0.02$ related-genome threshold, with a maximum $S_G = 0.0174$ (`1.74%`) to *Parageobacillus* phage JGon-2020a (`CP063417`).
* **NCBI `ClusteredNR` Homology (`RID: FUJZE5Z7014`):** Five consecutive structural and packaging proteins (`CTG_000311_16` HK97-gp10, `_18` head-tail connector, `_21` HK97 MCP, `_22` prohead protease, `_23` portal) share conserved synteny with uncultured Candidate Phyla Radiation (*Patescibacteria*: *Candidatus Moraniibacteriota* `HBI50266.1`–`HBI50277.1` and *Parcubacteria* `MFA5777622.1`–`MFA5777624.1`) MAGs at `33.00%–55.65%` amino acid identity, while `TerL` (`_24`) hits *Candidatus Altimarinota* (`MHY2402454.1`, `49.35%` id, $E = 10^{-149}$). This level of amino acid divergence (`33%–56%`) across core structural proteins supports recognition as a distinct genus (and potentially a new family-level clade once additional CPR-associated genomes are deposited).

## 3. Genome Completeness, Read-Level Validation, and Architecture

* **Circularity & Read Support:** CheckV identified a `55 bp` direct terminal repeat (`100.0%` completeness). Mapping raw paired-end reads (`ERR4194644`) with `minimap2 -ax sr` (`MAPQ >= 20`, `NM <= 3`) across the circularized `34,122 bp` genome and `2 kb` junction-bridging reference confirmed `100.00%` base coverage (`3,541` mapped reads, `15.45x` mean depth, `15.0x` median, `10x` minimum seam depth), with `8` reads spanning the left seam, `9` spanning the right seam, and `3` single reads spanning the entire `55 bp` DTR with at least `15 bp` overhang on both flanks.
* **Functional Organization:** `CTG_000311` encodes a complete HK97-like virion morphogenesis, packaging, tail, lysis (`phrog_291` endolysin), and replication/recombination (`phrog_44` SSB, `phrog_72` NinB, `phrog_298`/`phrog_10477` exonucleases) suite, plus a three-gene queuosine/7-deazapurine modification cassette (`phrog_887` queuine tRNA-ribosyltransferase, `phrog_2197` GTP cyclohydrolase, `phrog_2076` *QueE*-like radical SAM protein).
