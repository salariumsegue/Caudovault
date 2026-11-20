# Freshwater tailed-phage candidate screen

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24408019.svg)](https://doi.org/10.5281/zenodo.24408019)

Offline-ready Python pipeline that screens assembled metagenomic contigs (>=5 kb, priority >=10 kb) for candidate
tailed bacteriophages, checks terminal-repeat closure, and compares candidates to frozen reference snapshots.



## Setup
```bash
python3 -m venv .venv && source .venv/bin/activate
pip install biopython pandas numpy pyrodigal pyhmmer sourmash matplotlib seaborn networkx pyyaml
bash run_all.sh sanity     # synthetic controls; expect 5/5 PASS in results_sanity/tables/control_validation.tsv
```
Optional but worth it: `brew install blast` — real nucleotide alignment instead of the exact-31mer fallback.

## Real run
```bash
bash stage_data.sh                      # downloads assembly + PHROGs + reference genomes (~900 MB total)
source .venv/bin/activate
python3 src/01b_prepare_pharokka_phrogs.py --pharokka-db data/references/pharokka/db
bash run_all.sh real
```
`stage_data.sh ERZ18634931` stages a different ENA assembly instead of the default.

### Verified sources (checked 2026-11-17)
| Input | Source | Notes |
| --- | --- | --- |
| Assembly contigs | `https://ftp.sra.ebi.ac.uk/vol1/sequence/<ERZ prefix>/<ERZ>/contig.fa.gz` | ENA; browse with the ENA portal API |
| PHROGs HMMs + annotation | Pharokka db v1.11.0, Zenodo record 21755221 (735 MB, md5 `143bb375ddb0b0653e5cb5671f4a7629`) | ships all 38,880 PHROGs as a PyHMMER-ready `.h3m` |
| Phage reference genomes | `https://millardlab-inphared.s3.climb.ac.uk/9Aug2025_refseq_genomes.fa.gz` | INPHARED RefSeq subset |

`phrogs.lmge.uca.fr` currently serves an expired TLS certificate, so the Pharokka bundle is used instead of
building HMMs from PHROG MSAs. If that site is fixed, `src/01_build_phrog_hmms.py` still works from MSAs.

### Candidate freshwater assemblies (ENA, sizes verified 2026-11-17)
| Accession | Sample | Gzipped |
| --- | --- | --- |
| ERZ18661299 | lake Mekkojarvi MJ110629B (Finland) | 78 MB |
| ERZ18628969 | lake Alinen Mustajarvi AM4 (Finland) | 112 MB |
| ERZ18641798 | lake Keskinen Rajajarvi KR3 (Finland) | 130 MB |
| ERZ18624489 | lake Alinen Mustajarvi AM-2014-D7 | 176 MB |
| ERZ18652741 | lake Loclat Loc090116_1m (Switzerland) | 210 MB |
| ERZ18634931 | lake Erken Erken-D3 (Sweden) | 225 MB |

All from ENA study PRJEB47226. Find more:
```
https://www.ebi.ac.uk/ena/portal/api/search?result=analysis&query=assembly_type%3D%22primary%20metagenome%22%20AND%20tax_eq(1647806)&fields=analysis_accession,study_accession,sample_title,generated_bytes,generated_ftp&format=tsv&limit=100
```
(`tax_eq(1647806)` = lake water metagenome; `tax_eq(449393)` = freshwater metagenome.)

## Runtime expectation
The slow stage is 03: the hallmark/supporting PHROG profiles searched against every predicted protein. On 8 cores
expect roughly tens of minutes for a ~78 MB assembly and a few hours for a ~225 MB one. Start small, confirm the
whole thing runs, then scale up.

## Layout
`config.yaml` (real) | `config.sanity.yaml` (synthetic) | `src/` stages 00-07 | `data/` inputs |
`results*/` tables, sequences, figures | `report*/scientific_report.md` | `logs*/`

## Honesty rules the code enforces
Completeness is `NA` without CheckV, never 0%. Candidates are `heuristic_*_candidate` without geNomad, never
"validated virus". With no references staged, novelty is `not_assessed`, never "novel". DTRs mean a
terminal-repeat-supported assembly, not a circular genome. See `CLAUDE.md`.
