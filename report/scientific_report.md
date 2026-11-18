# Computational Screen for Candidate Tailed dsDNA Bacteriophages in an Assembled Metagenome

**Generated:** 2026-11-17 21:57 UTC  |  **Mode:** `sample`

## 1. Research question

Can a reproducible computational workflow recover well-supported tailed dsDNA phage sequences (*Caudoviricetes*) from an assembled freshwater metagenome that lack close matches in specified public reference databases? Finding zero such candidates is a valid result.

## 2. Dataset and preprocessing

- Raw records evaluated: 221228
- Retained (>=5 kb): 3388 (priority tier >=10 kb: 1502)
- Archived (<5 kb, not declared non-viral): 217810
- Excluded (corrupt/duplicate), see `01_exclusion_log.tsv`: 30

## 3. Candidates

93 contig(s) met the criteria (>=2 independent hallmark categories, or 1 hallmark plus >=3 supporting genes, without cellular housekeeping contamination).

| Candidate | Original ID | Length (bp) | GC (%) | Coding density | Hallmarks | Closure | Closest ref | Approx WGS (%) | Interpretation |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CTG_000043 | ERZ18661299.43 | 88860 | 34.508 | 0.9165 | major_capsid | no_terminal_repeat_detected | NC_029119 | 0.0 | no_close_match_in_frozen_references_provisional_candidate |
| CTG_000051 | ERZ18661299.51 | 82224 | 33.716 | 0.9035 | terminase_large | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_000058 | ERZ18661299.58 | 77643 | 63.333 | 0.916 | major_capsid;terminase_large | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_000065 | ERZ18661299.65 | 74764 | 34.641 | 0.8946 | major_capsid;portal;terminase_large | no_terminal_repeat_detected | NC_031057 | 0.02 | no_close_match_in_frozen_references_provisional_candidate |
| CTG_000111 | ERZ18661299.111 | 61836 | 54.326 | 0.9042 | major_capsid;portal | terminal-repeat-supported assembly (DTR) | NC_041892 | 0.0 | no_close_match_in_frozen_references_provisional_candidate |
| CTG_000122 | ERZ18661299.122 | 57818 | 39.789 | 0.9257 | major_capsid;portal;terminase_large | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_000129 | ERZ18661299.129 | 54726 | 42.0 | 0.9003 | major_capsid;portal;terminase_large | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_000137 | ERZ18661299.137 | 52985 | 50.598 | 0.9281 | major_capsid;portal | no_terminal_repeat_detected | NC_031120 | 0.0 | no_close_match_in_frozen_references_provisional_candidate |
| CTG_000145 | ERZ18661299.145 | 51844 | 37.032 | 0.9294 | major_capsid;terminase_large | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_000153 | ERZ18661299.153 | 49642 | 46.351 | 0.8901 | major_capsid;portal;terminase_large | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_000156 | ERZ18661299.156 | 49093 | 40.199 | 0.9637 | major_capsid;portal;terminase_large | no_terminal_repeat_detected | NC_023595 | 0.0 | no_close_match_in_frozen_references_provisional_candidate |
| CTG_000169 | ERZ18661299.169 | 46864 | 39.507 | 0.9128 | major_capsid | no_terminal_repeat_detected | NC_022762 | 0.0 | no_close_match_in_frozen_references_provisional_candidate |
| CTG_000191 | ERZ18661299.191 | 45202 | 35.899 | 0.9783 | portal;terminase_large | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_000194 | ERZ18661299.194 | 44832 | 31.347 | 0.8832 | terminase_large | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_000253 | ERZ18661299.253 | 38787 | 65.006 | 0.8763 | major_capsid;portal | terminal-repeat-supported assembly (DTR) | NC_023701 | 0.0 | no_close_match_in_frozen_references_provisional_candidate |
| CTG_000273 | ERZ18661299.273 | 37077 | 44.623 | 0.9261 | major_capsid;portal | terminal-repeat-supported assembly (DTR) | NC_049436 | 0.0 | no_close_match_in_frozen_references_provisional_candidate |
| CTG_000281 | ERZ18661299.281 | 36497 | 42.577 | 0.9081 | major_capsid;portal;terminase_large | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_000286 | ERZ18661299.286 | 36386 | 38.963 | 0.9396 | major_capsid;portal;terminase_large | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_000296 | ERZ18661299.296 | 35483 | 53.395 | 0.9442 | major_capsid;portal;terminase_large | no_terminal_repeat_detected | NC_061432 | 0.0 | no_close_match_in_frozen_references_provisional_candidate |
| CTG_000298 | ERZ18661299.298 | 35356 | 39.34 | 0.9647 | major_capsid;portal;terminase_large | terminal-repeat-supported assembly (DTR) | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_000311 | ERZ18661299.311 | 34177 | 37.437 | 0.9411 | major_capsid;portal | terminal-repeat-supported assembly (DTR) | NC_070936 | 0.0 | no_close_match_in_frozen_references_provisional_candidate |
| CTG_000316 | ERZ18661299.316 | 34007 | 41.012 | 0.9494 | major_capsid;portal | no_terminal_repeat_detected | NC_031909 | 0.0 | no_close_match_in_frozen_references_provisional_candidate |
| CTG_000337 | ERZ18661299.337 | 32547 | 48.997 | 0.9682 | major_capsid | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_000348 | ERZ18661299.348 | 31879 | 32.836 | 0.9719 | major_capsid | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_000349 | ERZ18661299.349 | 31871 | 44.551 | 0.9046 | major_capsid;portal | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_000363 | ERZ18661299.363 | 31043 | 36.44 | 0.9739 | major_capsid;portal | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_000364 | ERZ18661299.364 | 31035 | 47.35 | 0.905 | portal;terminase_large | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_000424 | ERZ18661299.424 | 28187 | 65.12 | 0.9273 | major_capsid;portal;terminase_large | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_000431 | ERZ18661299.431 | 27909 | 55.631 | 0.924 | major_capsid;portal | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_000498 | ERZ18661299.498 | 25177 | 44.791 | 0.9053 | portal;terminase_large | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_000564 | ERZ18661299.564 | 22932 | 41.797 | 0.9519 | terminase_large | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_000573 | ERZ18661299.573 | 22554 | 46.475 | 0.8845 | major_capsid;terminase_large | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_000581 | ERZ18661299.581 | 22443 | 41.532 | 0.886 | portal;terminase_large | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_000600 | ERZ18661299.600 | 22027 | 42.648 | 0.9634 | major_capsid | no_terminal_repeat_detected | NC_022968 | 0.0 | no_close_match_in_frozen_references_provisional_candidate |
| CTG_000623 | ERZ18661299.623 | 21504 | 49.205 | 0.948 | major_capsid | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_000635 | ERZ18661299.635 | 21181 | 56.905 | 0.9688 | major_capsid;portal | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_000647 | ERZ18661299.647 | 20880 | 52.313 | 0.9388 | major_capsid;portal;terminase_large | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_000648 | ERZ18661299.648 | 20838 | 55.01 | 0.8979 | major_capsid | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_000672 | ERZ18661299.672 | 20490 | 49.276 | 0.8793 | major_capsid;portal;terminase_large | no_terminal_repeat_detected | NC_019544 | 0.0 | no_close_match_in_frozen_references_provisional_candidate |
| CTG_000692 | ERZ18661299.692 | 20060 | 47.164 | 0.9653 | major_capsid | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_000699 | ERZ18661299.699 | 19982 | 61.541 | 0.9773 | major_capsid;portal;terminase_large | no_terminal_repeat_detected | NC_054664 | 0.0 | no_close_match_in_frozen_references_provisional_candidate |
| CTG_000703 | ERZ18661299.703 | 19955 | 53.731 | 0.9186 | major_capsid;portal | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_000708 | ERZ18661299.708 | 19912 | 58.898 | 0.9489 | major_capsid;terminase_large | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_000710 | ERZ18661299.710 | 19843 | 47.18 | 0.9018 | portal;terminase_large | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_000725 | ERZ18661299.725 | 19512 | 36.946 | 0.913 | major_capsid;portal | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_000736 | ERZ18661299.736 | 19197 | 50.867 | 0.9821 | major_capsid;portal;terminase_large | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_000767 | ERZ18661299.767 | 18419 | 62.247 | 0.9707 | major_capsid;terminase_large | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_000780 | ERZ18661299.780 | 18161 | 36.276 | 0.9613 | major_capsid;portal;terminase_large | no_terminal_repeat_detected | NC_020855 | 0.03 | no_close_match_in_frozen_references_provisional_candidate |
| CTG_000816 | ERZ18661299.816 | 17339 | 62.051 | 0.964 | major_capsid;portal | no_terminal_repeat_detected | NC_026606 | 0.0 | no_close_match_in_frozen_references_provisional_candidate |
| CTG_000826 | ERZ18661299.826 | 17129 | 59.686 | 0.9315 | major_capsid;terminase_large | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_000853 | ERZ18661299.853 | 16736 | 36.574 | 0.9549 | terminase_large | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_000961 | ERZ18661299.961 | 15067 | 44.678 | 0.9004 | major_capsid;terminase_large | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_001013 | ERZ18661299.1013 | 14401 | 63.422 | 0.9699 | major_capsid | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_001016 | ERZ18661299.1016 | 14343 | 36.534 | 0.926 | major_capsid | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_001038 | ERZ18661299.1038 | 14028 | 51.012 | 0.9281 | major_capsid;portal;terminase_large | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_001043 | ERZ18661299.1043 | 13989 | 42.869 | 0.9629 | major_capsid;portal;terminase_large | no_terminal_repeat_detected | NC_021559 | 0.0 | no_close_match_in_frozen_references_provisional_candidate |
| CTG_001095 | ERZ18661299.1095 | 13334 | 41.092 | 0.9234 | portal;terminase_large | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_001110 | ERZ18661299.1110 | 13139 | 55.204 | 0.9794 | major_capsid;portal;terminase_large | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_001127 | ERZ18661299.1127 | 12904 | 58.269 | 0.9392 | major_capsid;portal | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_001164 | ERZ18661299.1164 | 12561 | 47.608 | 0.8811 | portal;terminase_large | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_001251 | ERZ18661299.1251 | 11879 | 36.468 | 0.9718 | portal | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_001301 | ERZ18661299.1301 | 11492 | 54.714 | 0.9052 | major_capsid | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_001323 | ERZ18661299.1323 | 11271 | 64.457 | 0.9731 | major_capsid;portal | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_001335 | ERZ18661299.1335 | 11207 | 45.908 | 0.9313 | major_capsid;portal | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_001348 | ERZ18661299.1348 | 11129 | 36.697 | 0.96 | major_capsid | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_001416 | ERZ18661299.1416 | 10623 | 54.58 | 0.9561 | major_capsid | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_001426 | ERZ18661299.1426 | 10509 | 49.358 | 0.9801 | major_capsid;portal | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_001434 | ERZ18661299.1434 | 10479 | 46.407 | 0.9071 | major_capsid;portal;terminase_large | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_001539 | ERZ18661299.1539 | 9785 | 36.755 | 0.9323 | portal | no_terminal_repeat_detected | NC_020845 | 0.05 | no_close_match_in_frozen_references_provisional_candidate |
| CTG_001623 | ERZ18661299.1623 | 9355 | 38.707 | 0.9395 | portal;terminase_large | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_001795 | ERZ18661299.1795 | 8516 | 54.004 | 0.922 | major_capsid;portal;terminase_large | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_001836 | ERZ18661299.1836 | 8309 | 50.506 | 0.944 | portal | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_001873 | ERZ18661299.1873 | 8135 | 39.828 | 0.9888 | portal;terminase_large | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_001884 | ERZ18661299.1884 | 8122 | 40.483 | 0.9733 | major_capsid;portal;terminase_large | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_001942 | ERZ18661299.1942 | 7904 | 41.005 | 0.9552 | major_capsid;portal;terminase_large | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_002003 | ERZ18661299.2003 | 7713 | 56.022 | 0.9796 | portal;terminase_large | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_002010 | ERZ18661299.2010 | 7696 | 55.6 | 0.9585 | major_capsid;portal | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_002065 | ERZ18661299.2065 | 7524 | 44.019 | 0.9478 | major_capsid;portal | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_002075 | ERZ18661299.2075 | 7482 | 40.631 | 0.9566 | portal;terminase_large | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_002123 | ERZ18661299.2123 | 7334 | 59.531 | 0.9182 | major_capsid;terminase_large | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_002155 | ERZ18661299.2155 | 7272 | 48.804 | 0.971 | portal;terminase_large | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_002181 | ERZ18661299.2181 | 7206 | 61.31 | 0.9288 | major_capsid | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_002183 | ERZ18661299.2183 | 7201 | 54.631 | 0.8611 | portal;terminase_large | no_terminal_repeat_detected | NC_070773 | 0.0 | no_close_match_in_frozen_references_provisional_candidate |
| CTG_002419 | ERZ18661299.2419 | 6580 | 43.86 | 0.9919 | portal;terminase_large | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_002496 | ERZ18661299.2496 | 6386 | 39.132 | 0.9956 | portal;terminase_large | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_002578 | ERZ18661299.2578 | 6228 | 48.205 | 0.8565 | portal;terminase_large | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_002653 | ERZ18661299.2653 | 6073 | 52.84 | 0.7206 | portal;terminase_large | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_002657 | ERZ18661299.2657 | 6061 | 52.45 | 0.902 | major_capsid | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_002731 | ERZ18661299.2731 | 5929 | 43.819 | 0.9973 | portal;terminase_large | no_terminal_repeat_detected | NC_055877 | 0.0 | no_close_match_in_frozen_references_provisional_candidate |
| CTG_002833 | ERZ18661299.2833 | 5738 | 51.586 | 0.9582 | portal;terminase_large | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_003146 | ERZ18661299.3146 | 5288 | 54.784 | 0.9314 | major_capsid | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_003214 | ERZ18661299.3214 | 5196 | 67.071 | 0.9577 | portal;terminase_large | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |
| CTG_003216 | ERZ18661299.3216 | 5193 | 41.055 | 0.9623 | major_capsid;portal | no_terminal_repeat_detected | none | NA | no_sketch_or_anchor_match_in_references_provisional_candidate |

## 4. ViPTree Proteomic Tree & Clade Placement (Priority 7 Contigs)

To resolve evolutionary relationships beyond nucleotide k-mers, the 7 priority contigs (`CTG_000058`, `CTG_000137`, `CTG_000145`, `CTG_000169`, `CTG_000273`, `CTG_000296`, `CTG_000311`) were analyzed with **ViPTree v4.2** ([session `u6lWCt_eU2Wt7KD86pfmNw`](https://www.genome.jp/viptree/u/u6lWCt_eU2Wt7KD86pfmNw)) against `7,097` prokaryotic dsDNA reference viral genomes (`7,104` total rows in `set.ordered.lst.w.sim.pdf`, Prodigal genetic code 11). Full results are archived in `results/priority_7/viptree_7contigs_summary.tsv`.

- **Five contigs exceed ViPTree's $S_G > 0.02$ (`2.0%`) related-genome cutoff**, placing them into defined viral clades where nucleotide k-mers (`sourmash` / exact 31-mers) yielded `0.0%` WGS similarity:
  - **`CTG_000296`** (`35,483 bp`, `53.40%` GC, tree order `#6061`): Strongest hit at **$S_G = 0.2290$** (`22.90%`) to *Ralstonia phage YO374_1* (`ON170342`), followed by *Ralstonia phage YO108_2* (`ON170191`, `0.1661`) and *Klebsiella phage ST405-OXA48phi1.1* (`MK388859`, `0.1419`) — a **novel genus** within the Mu-like transposable *Pseudomonadota* phage clade.
  - **`CTG_000169`** (`46,864 bp`, `39.51%` GC, tree order `#4237`): Top hits to *Shewanella phage SppYZU01* (`KY622015`, **$S_G = 0.0909$**) and *Acinetobacter phage Scuro* (`OR348818`, `0.0901`), placing it in an **aquatic *Gammaproteobacteria* phage clade**.
  - **`CTG_000058`** (`77,643 bp`, `63.33%` GC, tree order `#6433`): Top hits to *Streptomyces phage Gilgamesh* (`NC_073461`, **$S_G = 0.0865$**) and freshwater *Microcystis phage Mwe-JY05* (`PP395638`, `0.0795`), consistent with a **high-GC *Actinomycetota* / cyanophage-allied lineage**.
  - **`CTG_000137`** (`52,985 bp`, `50.60%` GC, tree order `#4407`): Top hits to *Burkholderia phage phiE131* (`MH809532`, **$S_G = 0.0858$**) and *Aeromonas phage 62AhydR11PP* (`MH179474`, `0.0824`), bridging *Anamaviridae* and *Hirszfeldviridae*. Notably, `CTG_000137` and `CTG_000169` are the only pair among the 7 priority contigs with cross-contig proteomic similarity (**$S_G = 0.0427$**), corroborating five shared syntenic proteins in `results/priority_7/p7_vs_all93.tsv`.
  - **`CTG_000273`** (`37,077 bp`, DTR-complete, `44.62%` GC, tree order `#1007`): Top hit to *Bordetella phage BPP-1* (`NC_005357`, **$S_G = 0.0438$**), branching directly adjacent to *Pelagibacter phage HTVC106P* (`MN698244`, `0.0248`) and *Microcystis phage MJing1* (`OQ540923`, `0.0243`).
- **Two contigs are true "orphan" lineages below ViPTree's $S_G = 0.02$ cutoff**:
  - **`CTG_000311`** (`34,177 bp`, **DTR-complete**, `37.44%` GC, tree order `#5795`): Maximum $S_G$ across all `7,097` reference viruses is **$S_G = 0.0174$** (`1.74%`) to *Parageobacillus phage JGon-2020a* (`CP063417`), sitting at a deep basal node preceding *Gordonia*/*Rhodococcus* phages (`SG <= 0.0132`).
  - **`CTG_000145`** (`51,844 bp`, `37.03%` GC, tree order `#2217`): Maximum $S_G$ is **$S_G = 0.0049$** (`0.49%`) to *Bacillus phage vB_BceH_LBC2* (`PV449003`), sitting on an isolated branch with `0.0` similarity to immediate tree neighbors.
- **NCBI `ClusteredNR` BLASTp Resolution of the Two ViPTree Orphans (`RID: FUJZE5Z7014`)**: Querying all 11 hallmark proteins from `CTG_000145` (5 queries) and `CTG_000311` (6 queries) against NCBI `ClusteredNR` (`results/priority_7/ncbi_clusterednr_orphans_summary.tsv`) confirmed that both ViPTree orphans belong to uncultured environmental/MAG phage lineages absent from cultured RefSeq/ViPTree databases while remaining novel at the genus/subfamily level (`~33–69%` aa identity):
  - **`CTG_000311` (uncultured CPR / *Patescibacteria* HK97-like phage lineage):** Five structural/packaging loci (`CTG_000311_16`, `_18`, `_21`, `_22`, `_23`) match syntenic consecutive loci in uncultured *Candidatus Moraniibacteriota* (`HBI50266.1`–`HBI50277.1`, `33.00%–55.65%` id) and *Parcubacteria* (`MFA5777622.1`–`MFA5777624.1`, `42.37%–53.21%` id), while `TerL` (`CTG_000311_24`) hits *Candidatus Altimarinota* (`MHY2402454.1`, `49.35%` id, `E = 1e-149`).
  - **`CTG_000145` (uncultured Mu-like `gpT` environmental MAG prophage lineage):** `CTG_000145_42` hits Mu-like prophage major head subunit `gpT` in *Saccharimonadales* (`HEY1645038.1`, `66.88%` id, `100%` cover, `E = 7e-155`), with `TerL` (`_31`, `HEY1644814.1`, `58.07%` id, `E = 0.0`) and portal (`_37`, `HEY1645263.1`, `46.27%` id, `E = 0.0`) matching the same *Saccharimonadales*/*Sulfuriferula*/*Bacteroidia*/*Occallatibacter* environmental prophage cluster.
- **Correction of initial k-mer anchors**: Single 31-mer/k21 hits in Section 3 (e.g., `CTG_000137` -> `NC_031120` at $S_G = 0.0058$; `CTG_000273` -> `NC_049436` at $S_G = 0.0$; `CTG_000311` -> `NC_070936` at $S_G = 0.0$) were k-mer noise superseded by the whole-proteome ViPTree placement above.

## 5. Limitations and required follow-up

1. Without geNomad/CheckV, candidates are heuristic and completeness is `NA`, not an estimate (subsequent CheckV analysis confirmed DTR-complete status for `CTG_000273` and `CTG_000311`).
2. Direct terminal repeats support a terminally redundant assembly, and raw paired-end read mapping (`ERR4194644`, `minimap2 -ax sr`, `MAPQ >= 20`, `NM <= 3`, `>= 15 bp` overhang; `results/priority_7/dtr_junction_verification.tsv`, `results/figures/fig3_dtr_junction_coverage.png`) has now confirmed continuous circular/terminally redundant junction coverage across both DTR-complete genomes: `CTG_000273` (`29.72x` mean depth, `100.00%` bases covered, `157` left-seam, `150` right-seam, and `116` full-DTR spanning reads; min seam depth `23x`) and `CTG_000311` (`15.45x` mean depth, `100.00%` bases covered, `8` left-seam, `9` right-seam, and `3` full-DTR spanning reads; min seam depth `10x`).
3. ViPTree (`7,097` prokaryotic dsDNA reference genomes) establishes novelty relative to cultured/RefSeq dsDNA phages, and the follow-up NCBI `ClusteredNR` BLASTp check (`RID: FUJZE5Z7014`) resolved both ViPTree orphans (`CTG_000311` and `CTG_000145`) to uncultured MAG/CPR prophage clades at `33–69%` aa identity.
4. Both `CTG_000311` and `CTG_000145` remain provisionally novel genus/subfamily-level lineages relative to `ClusteredNR` and ViPTree v4.2, with `CTG_000273` and `CTG_000311` now validated at the raw-read level for INSDC/ICTV submission.
