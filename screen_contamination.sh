#!/usr/bin/env bash
# Host/prophage screen of the shortlist with geNomad. Shortlist is ~2 MB, so this is light.
set -euo pipefail
cd "$(dirname "$0")"
source .venv/bin/activate 2>/dev/null || true
command -v mmseqs >/dev/null || brew install mmseqs2
python -c "import genomad" 2>/dev/null || pip install genomad
[ -d data/genomad_db/genomad_db ] || genomad download-database data/genomad_db
genomad end-to-end --cleanup --splits 4 results/sequences/shortlist_candidates.fna results/genomad data/genomad_db/genomad_db
python src/08_merge_genomad.py results/genomad results
