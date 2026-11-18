#!/usr/bin/env bash
# Usage: bash run_all.sh [auto|sanity|real]
#   sanity : synthetic control panel only (config.sanity.yaml -> results_sanity/)
#   real   : real assembly in data/raw/assembly.fasta.gz (config.yaml -> results/)
#   auto   : real if the assembly is staged, otherwise sanity
set -euo pipefail
cd "$(dirname "$0")"

MODE="${1:-auto}"
if [ "$MODE" = "auto" ]; then
  if [ -f data/raw/assembly.fasta.gz ] || [ -f data/raw/assembly.fasta ] || [ -f data/raw/assembly.fna.gz ]; then MODE=real; else MODE=sanity; fi
fi

if [ "$MODE" = "sanity" ]; then
  export PHAGE_CONFIG=config.sanity.yaml
else
  export PHAGE_CONFIG=config.yaml
fi
echo "=== Mode: $MODE (config: $PHAGE_CONFIG) ==="

echo "=== Stage 0: preflight + smoke test ==="
python3 src/00_preflight.py

if [ "$MODE" = "sanity" ]; then
  echo "=== Stage 0b: synthetic control panel + control HMMs ==="
  python3 src/00_make_sanity_panel.py
else
  echo "=== Stage 1: PHROG MSAs -> HMMs ==="
  python3 src/01_build_phrog_hmms.py
fi

echo "=== Stage 2: contig QC, dedup, length split ==="
python3 src/02_preprocess_contigs.py
echo "=== Stage 3: gene calling + hallmark HMM screen ==="
python3 src/03_identify_candidates.py
echo "=== Stage 4: terminal repeats + CheckV (if present) ==="
python3 src/04_assess_quality.py
echo "=== Stage 5: novelty vs frozen references ==="
python3 src/05_assess_novelty.py
echo "=== Stage 6: annotation + figures ==="
python3 src/06_annotate_and_plot.py
echo "=== Stage 7: results table + report ==="
python3 src/07_compile_results_and_report.py
echo "=== Done ($MODE). See results*/tables/main_candidate_results.tsv and report*/scientific_report.md ==="
