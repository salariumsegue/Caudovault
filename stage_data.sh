#!/usr/bin/env bash
# Downloads every input the real run needs, into this project. Run on your Mac, not in a sandbox.
#   bash stage_data.sh                 # default assembly (Lake Mekkojarvi, ~78 MB) + references
#   bash stage_data.sh ERZ18634931     # a specific ENA assembly accession
# Re-running is safe: curl -C - resumes, and existing complete files are skipped.
set -euo pipefail
cd "$(dirname "$0")"

ERZ="${1:-ERZ18661299}"            # default: lake Mekkojarvi MJ110629B, ~78 MB
INPHARED_DATE="9Aug2025"
PHAROKKA_URL="https://zenodo.org/api/records/21755221/files/pharokka_v1.11.0_databases.tar.gz/content"
PHAROKKA_MD5="143bb375ddb0b0653e5cb5671f4a7629"

mkdir -p data/raw data/references/{phrogs_hmm,refseq,inphared,pharokka}

dl () {  # dl <url> <dest>
  if [ -s "$2" ]; then echo "  already present: $2"; else
    echo "  downloading -> $2"; curl -L -C - --fail --retry 3 -o "$2" "$1"
  fi
}

echo "[1/4] ENA assembly $ERZ"
PREFIX="$(echo "$ERZ" | cut -c1-6)"
dl "https://ftp.sra.ebi.ac.uk/vol1/sequence/${PREFIX}/${ERZ}/contig.fa.gz" "data/raw/assembly.fasta.gz"

echo "[2/4] Pharokka database bundle (~735 MB; includes PHROGs HMMs + annotation)"
dl "$PHAROKKA_URL" "data/references/pharokka/pharokka_db.tar.gz"
echo "  verifying md5..."
GOT=$(md5 -q data/references/pharokka/pharokka_db.tar.gz 2>/dev/null || md5sum data/references/pharokka/pharokka_db.tar.gz | cut -d' ' -f1)
if [ "$GOT" != "$PHAROKKA_MD5" ]; then echo "  MD5 MISMATCH: got $GOT, expected $PHAROKKA_MD5"; exit 1; fi
echo "  md5 ok"
if [ ! -d data/references/pharokka/db ]; then
  echo "  extracting..."; mkdir -p data/references/pharokka/db
  tar -xzf data/references/pharokka/pharokka_db.tar.gz -C data/references/pharokka/db --strip-components=1
fi

echo "[3/4] INPHARED phage genomes (RefSeq subset) for the novelty comparison"
dl "https://millardlab-inphared.s3.climb.ac.uk/${INPHARED_DATE}_refseq_genomes.fa.gz" \
   "data/references/refseq/refseq_phages.fna.gz"

echo "[4/4] Recording sizes and checksums in the manifest log"
{
  echo -e "file\tbytes\tsha256\tdownloaded_utc"
  for f in data/raw/assembly.fasta.gz data/references/pharokka/pharokka_db.tar.gz data/references/refseq/refseq_phages.fna.gz; do
    [ -s "$f" ] || continue
    SZ=$(wc -c < "$f" | tr -d ' ')
    SH=$(shasum -a 256 "$f" | cut -d' ' -f1)
    echo -e "${f}\t${SZ}\t${SH}\t$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  done
} > metadata/staged_files.tsv
cat metadata/staged_files.tsv

cat <<EOF

Staged. Assembly accession: $ERZ  (record this in metadata/dataset_manifest.tsv)
Next:
  source .venv/bin/activate
  python3 src/01b_prepare_pharokka_phrogs.py --pharokka-db data/references/pharokka/db
  bash run_all.sh real
EOF
