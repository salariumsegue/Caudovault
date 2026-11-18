#!/usr/bin/env python3
"""
src/01b_prepare_pharokka_phrogs.py
Alternative to 01_build_phrog_hmms.py when you have the Pharokka database bundle, which ships
PHROGs already built as a PyHMMER-readable .h3m (all 38,880 profiles). That skips building HMMs
from MSAs entirely.

The shipped profiles carry a PHROG id but usually no functional description, and this pipeline
classifies profiles BY description, so this script injects 'category | annot' from
phrog_annot_v4.tsv into each profile, then writes:
  data/references/phrogs_hmm/all_phrogs.hmm   (all profiles, binary h3m)
  data/references/phrogs_hmm/hallmarks.hmm    (hallmark + supporting subset, binary h3m)

Usage:
  python3 src/01b_prepare_pharokka_phrogs.py --pharokka-db /path/to/pharokka_db
  python3 src/01b_prepare_pharokka_phrogs.py --hmm <all_phrogs.h3m> --annot <phrog_annot_v4.tsv>
"""

import argparse
import logging
import re
import sys
from pathlib import Path

import pandas as pd
from pyhmmer.plan7 import HMMFile

from common import kw_match, load_cfg, setup_logging, text

NUM = re.compile(r"(\d+)")


def find_inputs(db_dir):
    """Locate the all-PHROGs profile file and annotation table inside a Pharokka db directory."""
    hmm = sorted([p for p in db_dir.rglob("*.h3m") if "phrog" in p.name.lower()] or
                 list(db_dir.rglob("*.h3m")), key=lambda p: -p.stat().st_size)
    annot = (list(db_dir.rglob("phrog_annot*.tsv")) or list(db_dir.rglob("*annot*.tsv")))
    return (hmm[0] if hmm else None), (annot[0] if annot else None)


def load_annotations(path):
    """phrog_annot_v4.tsv -> {phrog number as str: 'category | annot'}; tolerant of column naming."""
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("unknown")
    cols = {c.lower(): c for c in df.columns}
    id_col = cols.get("phrog") or list(df.columns)[0]
    annot_col = cols.get("annot") or cols.get("annotation") or cols.get("product")
    cat_col = cols.get("category") or cols.get("cat")
    out = {}
    for _, r in df.iterrows():
        m = NUM.search(str(r[id_col]))
        if not m:
            continue
        a = str(r[annot_col]) if annot_col else "hypothetical protein"
        c = str(r[cat_col]) if cat_col else "unknown function"
        out[m.group(1)] = f"{c} | {a}"
    return out


def is_first_pass(desc, cfg):
    for kws in cfg["hmmer_search"]["hallmark_categories"].values():
        if any(kw_match(k, desc) for k in kws):
            return True
    return any(kw_match(k, desc) for k in cfg["hmmer_search"]["supporting_categories"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pharokka-db", type=Path, help="Directory of the unpacked Pharokka database")
    ap.add_argument("--hmm", type=Path, help="all-PHROGs .h3m (overrides --pharokka-db)")
    ap.add_argument("--annot", type=Path, help="phrog_annot_v4.tsv (overrides --pharokka-db)")
    args = ap.parse_args()

    cfg = load_cfg()
    setup_logging("01b_prepare_pharokka_phrogs", cfg)

    hmm_path, annot_path = args.hmm, args.annot
    if (not hmm_path or not annot_path) and args.pharokka_db:
        f_hmm, f_annot = find_inputs(args.pharokka_db)
        hmm_path = hmm_path or f_hmm
        annot_path = annot_path or f_annot
    if not hmm_path or not hmm_path.exists():
        sys.exit("Could not find the PHROGs .h3m. Pass --hmm explicitly.")
    if not annot_path or not annot_path.exists():
        sys.exit("Could not find phrog_annot_v4.tsv. Pass --annot explicitly.")
    logging.info(f"Profiles : {hmm_path}")
    logging.info(f"Annotation: {annot_path}")

    annot = load_annotations(annot_path)
    logging.info(f"Loaded {len(annot)} PHROG annotations.")

    out_all = Path(cfg["paths"]["phrogs_all_hmm"])
    out_hall = Path(cfg["paths"]["phrogs_hallmark_hmm"])
    out_all.parent.mkdir(parents=True, exist_ok=True)

    n = n_desc = n_hall = 0
    with HMMFile(str(hmm_path)) as hf, open(out_all, "wb") as fa, open(out_hall, "wb") as fh:
        for hmm in hf:
            name = text(hmm.name)
            m = NUM.search(name)
            desc = annot.get(m.group(1)) if m else None
            if desc:
                hmm.description = desc
                n_desc += 1
            else:
                desc = text(hmm.description) or "unknown function | hypothetical protein"
            if not name.startswith("phrog_") and m:
                hmm.name = f"phrog_{m.group(1)}"
            hmm.write(fa, binary=True)
            n += 1
            if is_first_pass(desc, cfg):
                hmm.write(fh, binary=True)
                n_hall += 1
            if n % 5000 == 0:
                logging.info(f"  ...{n} profiles processed")

    logging.info(f"Wrote {n} profiles -> {out_all} ({n_desc} given descriptions from the annotation table).")
    logging.info(f"Wrote {n_hall} hallmark/supporting profiles -> {out_hall}")
    if n_desc == 0:
        logging.error("No profile matched the annotation table; classification would find nothing. "
                      "Check that the .h3m profile names contain PHROG numbers.")
        sys.exit(1)
    if n_hall == 0:
        logging.error("Hallmark subset is empty; check the annotation table's category/annot columns.")
        sys.exit(1)


if __name__ == "__main__":
    main()
