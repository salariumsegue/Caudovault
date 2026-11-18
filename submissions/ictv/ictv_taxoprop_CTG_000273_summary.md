# ICTV Taxonomic Proposal Summary: New Genus *Evovirus* and Species *Evovirus mekkojarviense* (`CTG_000273`)

**Proposer:** Niall Gillen (Indiana University Bloomington, ngillen@iu.edu)  
**Subcommittee:** Bacterial Viruses Subcommittee (`Caudoviricetes`)  
**Companion Spreadsheet:** `ictv_taxoprop_CTG_000273.xlsx`  
**Note on Taxon Names:** *Evovirus* and *Evovirus mekkojarviense* are placeholder Latinized names (referencing Lake Mekkojärvi in the Evo region, Finland) formatted per ICTV binomial rules; insert the assigned GenBank/ENA TPA accession once `SUB16642552` is accepted.

## 1. Taxonomic Summary

We propose creating one new genus (**Evovirus**) and one new species (**Evovirus mekkojarviense**) within class *Caudoviricetes* to classify `CTG_000273` (`37,043 bp` DTR-trimmed / `37,077 bp` linear, `34 bp` DTR, `44.62%` GC, `47` joined ORFs), a complete, read-verified freshwater bacteriophage assembly recovered from Lake Mekkojärvi metagenome `ERZ18661299` (BioProject `PRJEB47226`, BioSample `SAMEA6872951`, raw reads `ERR4194644`).

## 2. Phylogenomic and Database-Wide Demarcation

* **Nucleotide Divergence:** `CTG_000273` exhibits `0.0%` aligned bases against frozen NCBI Viral RefSeq and surveyed freshwater viromes, and at most `5%` query coverage against full NCBI `core_nt` (`RID: FVPFG1V3014`; *Polynucleobacter paneuropaeus* `CP049682.1` at `63.71%` id, $E = 2 \times 10^{-53}$; MAG phage `BK031899.1` at `4%` cover, `64.58%` id; *Bordetella* phage BPP-1 `NC_005357.1` at `3%` cover, `65.16%` id). These represent weak sequence affinity only rather than a definitive host assignment, and local CRISPR-spacer screening across `5,165` spacers in `ERZ18661299` yielded `0` hits.
* **Whole-Proteome Tree (ViPTree v4.2):** Across `7,097` prokaryotic dsDNA reference viruses (session `u6lWCt_eU2Wt7KD86pfmNw`), `CTG_000273` reaches a maximum $S_G = 0.0438$ (`4.38%`) to *Bordetella* phage BPP-1 (`NC_005357`), well below genus-level relatedness and supporting a distinct genus.

## 3. Read-Level Junction Validation and Genome Architecture

* **DTR Continuity & In-Frame Gene Fusion:** Mapping `ERR4194644` paired-end reads (`minimap2 v2.31 -ax sr`, `MAPQ >= 20`, `NM <= 3`) confirmed `100.00%` base coverage (`7,478` reads, `29.72x` mean / `29.0x` median depth, `23x` min / `87.44x` mean seam depth), with `157` left-seam, `150` right-seam, and `116` full-DTR spanning reads (`>= 15 bp` overhang). Joining across the `34 bp` DTR fuses linear ORFs `_48` and `_1` in-frame into a single complete `162-aa` DNA repair protein (`phrog_10276`). Physical virion topology (circular vs. linear DTR) remains unresolved.
* **Functional Architecture:** Encodes a complete head-packaging-tail-replication suite including portal (`phrog_1216`), major head (`phrog_1328`), head scaffolding (`phrog_3927`), tail (`phrog_11262`), `TerS` (`phrog_23444`), polymerase (`phrog_2825`), DnaD-like helicase loader (`phrog_2514`), RecT (`phrog_148`), DNA methyltransferase (`phrog_757`), and SbcC-like endonuclease (`phrog_77`).
