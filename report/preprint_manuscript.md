# Identification and Read-Level Validation of Two Complete, Previously Uncharacterized Caudoviricetes Assemblies (`CTG_000273` and `CTG_000311`) from Lake Mekkojarvi

**Authors:** Niall Gillen¹*  
**Affiliations:** ¹Indiana University Bloomington, Bloomington, IN, USA  
**Correspondence:** ngillen@iu.edu  

## Abstract

Freshwater ecosystems harbor vast uncharacterized viral diversity, yet standard nucleotide k-mer screening often fails to resolve evolutionary relationships for divergent bacteriophages. Using a reproducible hallmark-guided computational screen of 221,228 assembled contigs from public Lake Mekkojarvi metagenome `ERZ18661299` (BioProject `PRJEB47226`, BioSample `SAMEA6872951`, raw reads `ERR4194644`), we recovered 93 candidate tailed dsDNA phage (*Caudoviricetes*) contigs and focused full phylogenomic, database-wide, and raw-read validation on two complete, direct-terminal-repeat (DTR)-supported genomes: `CTG_000273` (`37,077 bp` linear / `37,043 bp` DTR-trimmed, `34 bp` DTR, `44.62%` GC, `48` linear / `47` joined ORFs) and `CTG_000311` (`34,177 bp` linear / `34,122 bp` DTR-trimmed, `55 bp` DTR, `37.44%` GC, `51` ORFs). Both assemblies showed `0.0%` aligned bases against frozen NCBI Viral RefSeq and surveyed freshwater virome snapshots. Against full NCBI `core_nt` (sensitive BLASTn, `RID: FVPFG1V3014`), `CTG_000273` yielded at most `5%` query coverage (`63.71%` identity to *Polynucleobacter paneuropaeus* `CP049682.1`, $E = 2 \times 10^{-53}$; `4%` cover at `64.58%` to MAG phage `BK031899.1`), and ViPTree v4.2 placed it at $S_G = 0.0438$ to *Bordetella* phage BPP-1 (`NC_005357`), consistent with a previously uncharacterized freshwater *Burkholderiales*/*Polynucleobacter*-associated phage. `CTG_000311` fell below ViPTree's $S_G = 0.02$ relatedness cutoff ($S_G = 0.0174$ to *Parageobacillus* phage JGon-2020a `CP063417`), but full NCBI `ClusteredNR` BLASTp (`RID: FUJZE5Z7014`) revealed conserved synteny (`33%–56%` amino acid identity) across five consecutive structural/packaging proteins with uncultured Candidate Phyla Radiation (*Patescibacteria*: *Candidatus Moraniibacteriota* and *Parcubacteria*) MAGs alongside a three-gene queuosine/7-deazapurine modification cassette (*QueE*, GTP cyclohydrolase, queuine tRNA-ribosyltransferase). Mapping raw paired-end reads (`ERR4194644`) across DTR-bridging references confirmed `100.00%` base coverage with zero dropouts for both `CTG_000273` (`29.72x` mean depth, `116` full-DTR reads, fusing terminal ORFs `_48` and `_1` in-frame into a `162-aa` DNA repair protein) and `CTG_000311` (`15.45x` mean depth, `8` left-seam, `9` right-seam, and `3` full-DTR reads). These results establish `CTG_000273` and `CTG_000311` as complete, terminally redundant, DTR-joined phage assemblies (with physical virion circularity vs. linear DTR packaging left unresolved).

## 1. Introduction

Bacteriophages of the class *Caudoviricetes* drive microbial mortality, horizontal gene transfer, and biogeochemical cycling across freshwater ecosystems. Despite rapid growth in viral metagenomics, a large fraction of assembled freshwater viral contigs—often termed "viral dark matter"—cannot be classified by standard nucleotide k-mer or alignment tools because dsDNA phage genomes accumulate rapid synonymous and structural divergence across distant lineages. Even protein-based viral trees built primarily from cultured isolates under-represent phages infecting uncultured bacterial phyla, notably the Candidate Phyla Radiation (CPR / *Patescibacteria*), whose ultra-small, largely uncultivated hosts are abundant in humic and oligotrophic freshwater lakes.

Here, we apply a reproducible, multi-stage computational screen to public freshwater metagenome `ERZ18661299` from humic Lake Mekkojarvi (BioProject `PRJEB47226`, BioSample `SAMEA6872951`, raw reads `ERR4194644`) and perform end-to-end functional annotation, whole-proteome placement, full NCBI `core_nt`/`ClusteredNR` novelty verification, and raw-read DTR junction validation on two complete, DTR-supported genomes: `CTG_000273` and `CTG_000311`.

## 2. Methods

### 2.1 Metagenome Preprocessing and Hallmark-Guided Screening
Assembled contigs from ENA metagenome `ERZ18661299` (`221,228` raw records) were filtered by length (`>= 5 kb` retained, `n = 3,388`; priority tier `>= 10 kb`, `n = 1,502`; `< 5 kb` archived, `n = 217,810`; `30` corrupt/duplicate records excluded). Open reading frames were called with Prodigal v2.6.3 (genetic code 11) and scanned against PHROGs v4 profile HMMs. Contigs encoding at least two independent viral hallmark categories (`major_capsid`, `portal`, `terminase_large`) or one hallmark plus at least three supporting phage genes without cellular housekeeping contamination were retained (`n = 93`), and CheckV was used to identify complete, DTR-supported assemblies, selecting `CTG_000273` (`ERZ18661299.273`) and `CTG_000311` (`ERZ18661299.311`) for complete validation.

### 2.2 Novelty, Whole-Proteome Phylogenetics, and Database-Wide Checks
Nucleotide similarity was evaluated against frozen NCBI Viral RefSeq and surveyed freshwater virome snapshots (`0.0%` aligned bases for both `CTG_000273` and `CTG_000311`, with nearest k-mer anchors `NC_049436` and `NC_070936`). Whole-proteome phylogenetic placement was performed with ViPTree v4.2 (session `u6lWCt_eU2Wt7KD86pfmNw`, Prodigal code 11) against `7,097` prokaryotic dsDNA reference viral genomes (`7,104` total entries). To test novelty beyond RefSeq, `CTG_000273` (`37,043 bp` DTR-trimmed) was queried against full NCBI `core_nt` using sensitive `blastn` (`RID: FVPFG1V3014`), and all six hallmark proteins of ViPTree orphan `CTG_000311` (`S_G < 0.02`) were queried against full NCBI `ClusteredNR` using `blastp` (`RID: FUJZE5Z7014`).

### 2.3 Raw-Read DTR Junction Verification
Redundant 3' DTRs were trimmed from `CTG_000273` (`34 bp` -> `37,043 bp`) and `CTG_000311` (`55 bp` -> `34,122 bp`), and `2 kb` junction-bridging references (`1,000 bp` 3' flank + `1,000 bp` 5' flank, placing each DTR seam at position `1001`) were constructed (`results/priority_7/dtr_junction_refs.fna`). Paired-end raw reads from ENA run `ERR4194644` (`2.44 GB` + `2.55 GB`) were mapped using `minimap2 -ax sr` (v2.31) and filtered with `samtools` (v1.25; primary alignments, `MAPQ >= 20`, `NM <= 3`). Reads spanning the left seam, right seam, and entire DTR with at least `15 bp` overhang on both flanks were quantified, and genome-wide coverage profiles were computed with `samtools depth -a`.

## 3. Results

### 3.1 Genome Architecture and Functional Annotation of `CTG_000273` and `CTG_000311`

| Metric | `CTG_000273` (`ERZ18661299.273`) | `CTG_000311` (`ERZ18661299.311`) |
| --- | --- | --- |
| Linear / DTR-trimmed length | `37,077 bp` / `37,043 bp` | `34,177 bp` / `34,122 bp` |
| DTR length & CheckV tier | `34 bp` (Complete, 100%) | `55 bp` (Complete, 100%) |
| GC (%) & Coding density | `44.62%` / `0.9261` | `37.44%` / `0.9411` |
| Predicted ORFs | `48` linear (`47` joined) | `51` |
| ViPTree v4.2 Max $S_G$ | `0.0438` (*Bordetella* phage BPP-1 `NC_005357`) | `0.0174` (*Parageobacillus* phage JGon-2020a `CP063417`) |
| Full NCBI `core_nt` / `ClusteredNR` | Max `5%` cover (`63.71%` id, `RID: FVPFG1V3014`) | `33%–56%` aa id to CPR MAGs (`RID: FUJZE5Z7014`) |
| Mapped reads (`MAPQ >= 20`) | `7,478` | `3,541` |
| Mean / Median depth | `29.72x` / `29.0x` (`100.00%` covered) | `15.45x` / `15.0x` (`100.00%` covered) |
| Seam-spanning reads (`>= 15 bp` overhang) | `157` L / `150` R / `116` full-DTR (min `23x`) | `8` L / `9` R / `3` full-DTR (min `10x`) |

* **`CTG_000273` (Figure 1):** Encodes a complete head-packaging-tail-replication architecture including portal (`phrog_1216`, bit score `360.49`), major capsid (`phrog_1328`, `275.02`), head scaffolding (`phrog_3927`), tail protein (`phrog_11262`), `TerS` (`phrog_23444`), polymerase (`phrog_2825`), DnaD-like helicase loader (`phrog_2514`), RecT (`phrog_148`), DNA methyltransferase (`phrog_757`), and SbcC-like endonuclease (`phrog_77`).
* **`CTG_000311` (Figure 2):** Encodes a complete HK97-like virion morphogenesis and lysis suite—portal (`phrog_2662`, `436.21`), prohead maturation protease (`phrog_7015`), major capsid (`phrog_15928`), Neck1/tail completion (`phrog_5`), tail tape measure (`phrog_3981`), endolysin (`phrog_291`), SSB (`phrog_44`), NinB (`phrog_72`), and exonucleases (`phrog_298`, `phrog_10477`)—alongside a three-gene queuosine/7-deazapurine modification cassette (`phrog_887` queuine tRNA-ribosyltransferase, `phrog_2197` GTP cyclohydrolase, `phrog_2076` *QueE*-like radical SAM protein).

### 3.2 Phylogenomic and Database-Wide Novelty Verification
* **`CTG_000273`:** Across full NCBI `core_nt` (`RID: FVPFG1V3014`), sensitive `blastn` recovered at most `5%` query coverage (`63.71%` identity to *Polynucleobacter paneuropaeus* `CP049682.1`, $E = 2 \times 10^{-53}$; `5%` at `62.83%` to *Bordetella* phage `PV330005.1`; `4%` at `64.58%` to MAG phage `BK031899.1`; `3%` at `65.16%` to *Bordetella* phage BPP-1 `NC_005357.1`), confirming it is far below the `85%` coverage / `95%` ANI species threshold and suggesting a freshwater *Burkholderiales*/*Polynucleobacter* host.
* **`CTG_000311`:** Sitting below ViPTree's $S_G = 0.02$ cutoff (`0.0174`), `CTG_000311`'s five consecutive structural/packaging proteins (`_16`, `_18`, `_21`, `_22`, `_23`) matched syntenic consecutive loci in uncultured *Candidatus Moraniibacteriota* (`HBI50266.1`–`HBI50277.1`, `33.00%–55.65%` aa identity) and *Parcubacteria* (`MFA5777622.1`–`MFA5777624.1`, `42.37%–53.21%` aa identity) MAGs in `ClusteredNR` (`RID: FUJZE5Z7014`), while `TerL` (`_24`) hit *Candidatus Altimarinota* (`MHY2402454.1`, `49.35%` id, $E = 10^{-149}$).

### 3.3 Raw-Read Mapping Confirms Terminally Redundant Assembly Continuity
Mapping `ERR4194644` paired-end reads against the trimmed assemblies and DTR-bridging junction references (Figure 3) confirmed `100.00%` base coverage with zero dropouts for both genomes. Notably, joining `CTG_000273` across its `34 bp` DTR seam fuses linear ORFs `CTG_000273_48` and `CTG_000273_1` in-frame into a single complete `162-aa` DNA repair protein (`phrog_10276`), directly supported by `116` full-DTR spanning reads.

## 4. Discussion

Together, our results show that `CTG_000273` and `CTG_000311` are complete, non-chimeric, previously uncharacterized *Caudoviricetes* assemblies recovered from public freshwater metagenomic data. `CTG_000311` in particular represents a novel genus/subfamily-level member of an uncultured CPR (*Patescibacteria*) HK97-like phage clade that appears as an isolated orphan in cultured-phage proteomic trees solely because CPR phages are absent from RefSeq. While short-read DTR-seam mapping proves terminally redundant assembly continuity for both `37,043 bp` (`CTG_000273`) and `34,122 bp` (`CTG_000311`), physical virion topology (circular vs. linear with direct terminal repeats) remains unresolved without long-read or physical topology validation.

## 5. Data and Code Availability

All analysis scripts, summary TSVs, genome maps (`results/figures/fig2_genome_map_CTG_000273.png`, `fig2_genome_map_CTG_000311.png`), DTR coverage plots (`results/figures/fig3_dtr_junction_coverage.png`), and INSDC Third-Party Annotation files (`submissions/insdc_tpa/`) are archived at https://github.com/salariumsegue/Caudovault. Primary metagenome data are available under ENA assembly `ERZ18661299`, BioProject `PRJEB47226`, BioSample `SAMEA6872951`, and raw read run `ERR4194644`.
