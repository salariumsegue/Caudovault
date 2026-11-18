# Agent Operating Instructions

## Research question
Can a reproducible computational workflow recover well-supported tailed dsDNA phage sequences
(*Caudoviricetes*) from an assembled freshwater metagenome that lack close matches in specified
public reference databases? A null result is a valid outcome.

## Scope
Purely computational analysis of existing public sequence data (assembled contigs + reference
databases). No wet-lab protocols, no work on live organisms.

## Rules
1. **Zero fabrication.** Never invent accessions, checksums, completeness values or taxonomy. Unmeasured = `NA`.
2. **Offline sandbox.** The sandbox cannot reach NCBI. Do not try to download databases mid-run; the user stages
   all inputs from their own machine. `src/00_preflight.py` reports what is present; log sources in
   `metadata/dataset_manifest.tsv`.
3. **Provenance.** Never overwrite `data/raw/` or `data/references/`. Every excluded contig is logged with a reason
   in `<results>/tables/01_exclusion_log.tsv`.
4. **Fixed HMM search space.** Always pass `Z` from config (`fixed_Z: 100000`) to `pyhmmer.hmmsearch`.
5. **Conservative evidence.**
   - One gene hitting two hallmark profiles counts as ONE independent hallmark.
   - Generic replication enzymes never qualify a contig as viral on their own.
   - Without geNomad, label candidates `heuristic_*_candidate`, never "validated virus".
   - Without CheckV, completeness is `NA`, never `0%`.
   - DTRs = "terminal-repeat-supported assembly", not "circular genome".
   - Short conserved-gene identity is not whole-genome identity; report aligned fractions.
6. **Novelty wording.** "Provisionally novel relative to the frozen snapshot" only. With no references staged the
   result is `not_assessed_no_references_staged`, never "novel".
7. **Sanity vs real.** `bash run_all.sh sanity` uses `config.sanity.yaml` and a synthetic panel (`results_sanity/`).
   It must show 5/5 PASS in `control_validation.tsv` before trusting a real run. Sanity output is never
   scientific evidence. `bash run_all.sh real` uses `config.yaml` and refuses to run without real PHROG HMMs.

## Execution order
`run_all.sh` runs: 00_preflight -> (00_make_sanity_panel | 01_build_phrog_hmms) -> 02_preprocess_contigs ->
03_identify_candidates -> 04_assess_quality -> 05_assess_novelty -> 06_annotate_and_plot -> 07_compile_results_and_report.

## Known limits to state in any report
- Sketch / exact-anchor novelty only sees close relatives (roughly >85-90% nucleotide identity). Divergent relatives
  need protein-level comparison; BLAST+ improves nucleotide sensitivity when installed.
- Candidates are assembly-only; no read-mapping validation of junctions.
