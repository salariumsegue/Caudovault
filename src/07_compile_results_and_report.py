#!/usr/bin/env python3
"""
src/07_compile_results_and_report.py
Merges evidence, quality and novelty into main_candidate_results.tsv ('NA' for unmeasured values)
and writes report/scientific_report.md. In sanity mode it also validates the control panel and
exits non-zero if any control behaves wrongly.
"""

import logging
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from common import load_cfg, rdir, setup_logging


def read_or_empty(p):
    try:
        return pd.read_csv(p, sep="\t")
    except (FileNotFoundError, pd.errors.EmptyDataError):
        return pd.DataFrame()


def merge_unique(left, right):
    if right.empty or "contig_id" not in right.columns:
        return left
    return left.merge(right[[c for c in right.columns if c == "contig_id" or c not in left.columns]],
                      on="contig_id", how="left")


def validate_controls(ev, quality, novelty, tab_dir):
    def status(oid):
        r = ev[ev["original_id"] == oid]
        return r["candidate_status"].iloc[0] if len(r) else "absent"

    q = quality.merge(ev[["contig_id", "original_id"]], on="contig_id") if not quality.empty else quality
    pos_q = q[q["original_id"] == "CTRL_POS_TAILED_PHAGE_1"] if not q.empty else q
    pos_n = novelty.merge(ev[["contig_id", "original_id"]], on="contig_id") if not novelty.empty else novelty
    pos_n = pos_n[pos_n["original_id"] == "CTRL_POS_TAILED_PHAGE_1"] if not pos_n.empty else pos_n
    checks = [
        ("positive phage is a priority candidate", status("CTRL_POS_TAILED_PHAGE_1") == "heuristic_priority_candidate"),
        ("cellular control is excluded", status("CTRL_NEG_CELLULAR_HOUSEKEEPING_1") == "excluded_cellular_or_host_island"),
        ("short contig is archived (<5 kb), not screened", status("CTRL_SHORT_CONTIG_LT5KB") == "absent"),
        ("positive phage DTR detected", (not pos_q.empty) and pos_q["tr_type"].iloc[0] == "DTR"),
        ("positive phage matched itself as rediscovery",
         (not pos_n.empty) and pos_n["novelty_interpretation"].iloc[0] == "known_species_level_match_or_rediscovery"),
    ]
    df = pd.DataFrame([{"check": c, "result": "PASS" if ok else "FAIL"} for c, ok in checks])
    df.to_csv(tab_dir / "control_validation.tsv", sep="\t", index=False)
    for _, r in df.iterrows():
        logging.info(f"  [{r['result']}] {r['check']}")
    return bool((df["result"] == "PASS").all())


def main():
    cfg = load_cfg()
    setup_logging("07_compile_results", cfg)
    tab_dir = rdir(cfg, "tables")
    rep_dir = Path(cfg["paths"]["report_dir"])
    rep_dir.mkdir(parents=True, exist_ok=True)
    mode = cfg["run"]["mode"]

    prep = read_or_empty(tab_dir / "01_preprocessing_counts.tsv")
    ev = read_or_empty(tab_dir / "04_contig_candidate_evidence.tsv")
    short = read_or_empty(tab_dir / "04_shortlist_candidates.tsv")
    qual = read_or_empty(tab_dir / "05_candidate_quality.tsv")
    nov = read_or_empty(tab_dir / "06_novelty_comparisons.tsv")

    merged = short.copy()
    if not merged.empty:
        merged = merge_unique(merge_unique(merged, qual), nov).astype(object).fillna("NA")
    merged.to_csv(tab_dir / "main_candidate_results.tsv", sep="\t", index=False)

    ok = True
    if mode == "sanity":
        logging.info("Validating synthetic control panel:")
        ok = validate_controls(ev, qual, nov, tab_dir)

    g = lambda c: prep[c].iloc[0] if (not prep.empty and c in prep.columns) else "NA"  # noqa: E731
    L = []
    if mode == "sanity":
        L += ["> **SANITY RUN: synthetic control data only. No biological or scientific claims are made here.**", ""]
    L += ["# Computational Screen for Candidate Tailed dsDNA Bacteriophages in an Assembled Metagenome", "",
          f"**Generated:** {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}  |  **Mode:** `{mode}`", "",
          "## 1. Research question", "",
          "Can a reproducible computational workflow recover well-supported tailed dsDNA phage sequences "
          "(*Caudoviricetes*) from an assembled freshwater metagenome that lack close matches in specified "
          "public reference databases? Finding zero such candidates is a valid result.", "",
          "## 2. Dataset and preprocessing", "",
          f"- Raw records evaluated: {g('total_raw_records')}",
          f"- Retained (>=5 kb): {g('retained_ge5kb')} (priority tier >=10 kb: {g('priority_ge10kb')})",
          f"- Archived (<5 kb, not declared non-viral): {g('archived_short_lt5kb')}",
          f"- Excluded (corrupt/duplicate), see `01_exclusion_log.tsv`: {g('excluded_corrupt_or_duplicate')}", "",
          "## 3. Candidates", ""]
    if merged.empty:
        L += ["No contig met the conservative hallmark criteria. This is a valid null result; "
              "see `04_contig_candidate_evidence.tsv` for per-contig evidence.", ""]
    else:
        L += [f"{len(merged)} contig(s) met the criteria (>=2 independent hallmark categories, or 1 hallmark "
              "plus >=3 supporting genes, without cellular housekeeping contamination).", "",
              "| Candidate | Original ID | Length (bp) | GC (%) | Coding density | Hallmarks | Closure | Closest ref | Approx WGS (%) | Interpretation |",
              "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |"]
        for _, r in merged.iterrows():
            L.append("| " + " | ".join(str(r.get(c, "NA")) for c in [
                "contig_id", "original_id", "length_bp", "gc_percent", "coding_density",
                "hallmark_categories_detected", "closure_status", "closest_reference_accession",
                "approx_wgs_pct", "novelty_interpretation"]) + " |")
        L.append("")
    L += ["## 4. Limitations and required follow-up", "",
          "1. Without geNomad/CheckV, candidates are heuristic and completeness is `NA`, not an estimate.",
          "2. Direct terminal repeats support a terminally redundant assembly, not proof of a circular virion genome.",
          "3. Novelty is judged only against the frozen reference snapshots in the manifest. Sketch and exact-anchor "
          "methods detect only closely related genomes; divergent relatives need protein-level comparison "
          "(e.g. vConTACT/geNomad taxonomy) before any novelty claim.",
          "4. Candidates are assembly-only; read mapping across any repeat junction and coverage uniformity are required "
          "before formal taxonomic proposals.",
          "5. A candidate is 'provisionally novel relative to the snapshot', never 'confirmed new phage'."]
    (rep_dir / "scientific_report.md").write_text("\n".join(L) + "\n")
    logging.info(f"Wrote {tab_dir / 'main_candidate_results.tsv'} and {rep_dir / 'scientific_report.md'}")
    if not ok:
        logging.error("Control validation FAILED.")
        sys.exit(1)


if __name__ == "__main__":
    main()
