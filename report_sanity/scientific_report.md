> **SANITY RUN: synthetic control data only. No biological or scientific claims are made here.**

# Computational Screen for Candidate Tailed dsDNA Bacteriophages in an Assembled Metagenome

**Generated:** 2026-11-17 20:50 UTC  |  **Mode:** `sanity`

## 1. Research question

Can a reproducible computational workflow recover well-supported tailed dsDNA phage sequences (*Caudoviricetes*) from an assembled freshwater metagenome that lack close matches in specified public reference databases? Finding zero such candidates is a valid result.

## 2. Dataset and preprocessing

- Raw records evaluated: 3
- Retained (>=5 kb): 2 (priority tier >=10 kb: 2)
- Archived (<5 kb, not declared non-viral): 1
- Excluded (corrupt/duplicate), see `01_exclusion_log.tsv`: 0

## 3. Candidates

1 contig(s) met the criteria (>=2 independent hallmark categories, or 1 hallmark plus >=3 supporting genes, without cellular housekeeping contamination).

| Candidate | Original ID | Length (bp) | GC (%) | Coding density | Hallmarks | Closure | Closest ref | Approx WGS (%) | Interpretation |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CTG_000001 | CTRL_POS_TAILED_PHAGE_1 | 11106 | 49.604 | 0.8739 | major_capsid;portal;terminase_large | terminal-repeat-supported assembly (DTR) | REF_CTRL_POS_TAILED_PHAGE_1 | 99.46 | known_species_level_match_or_rediscovery |

## 4. Limitations and required follow-up

1. Without geNomad/CheckV, candidates are heuristic and completeness is `NA`, not an estimate.
2. Direct terminal repeats support a terminally redundant assembly, not proof of a circular virion genome.
3. Novelty is judged only against the frozen reference snapshots in the manifest. Sketch and exact-anchor methods detect only closely related genomes; divergent relatives need protein-level comparison (e.g. vConTACT/geNomad taxonomy) before any novelty claim.
4. Candidates are assembly-only; read mapping across any repeat junction and coverage uniformity are required before formal taxonomic proposals.
5. A candidate is 'provisionally novel relative to the snapshot', never 'confirmed new phage'.
