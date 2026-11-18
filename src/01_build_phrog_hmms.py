#!/usr/bin/env python3
"""
src/01_build_phrog_hmms.py
Converts downloaded PHROG alignments into HMMER3 profiles (pyhmmer Builder.build_msa).
Writes all_phrogs.hmm plus hallmarks.hmm (hallmark + supporting profiles for the first-pass screen).
"""

import logging
import sys
from pathlib import Path

import pandas as pd
from pyhmmer.easel import Alphabet, MSAFile
from pyhmmer.plan7 import Background, Builder

from common import kw_match, load_cfg, setup_logging


def load_annotations(tsv_path):
    if not tsv_path.exists():
        logging.warning(f"PHROG annotation TSV not found at {tsv_path}; descriptions will be 'unknown'.")
        return {}
    df = pd.read_csv(tsv_path, sep="\t", dtype=str).fillna("unknown")
    meta = {}
    for _, row in df.iterrows():
        raw = str(row.get("phrog", "")).strip()
        pid = raw if raw.startswith("phrog_") else f"phrog_{raw}"
        meta[pid] = {"annot": row.get("annot", "hypothetical protein"),
                     "category": row.get("category", "unknown function")}
    return meta


def is_first_pass(annot, category, cfg):
    desc = f"{category} {annot}"
    for kws in cfg["hmmer_search"]["hallmark_categories"].values():
        if any(kw_match(k, desc) for k in kws):
            return True
    return any(kw_match(k, desc) for k in cfg["hmmer_search"]["supporting_categories"])


def main():
    cfg = load_cfg()
    setup_logging("01_build_phrog_hmms", cfg)

    msa_dir = Path(cfg["paths"]["phrogs_msa_dir"])
    out_hall = Path(cfg["paths"]["phrogs_hallmark_hmm"])
    out_all = Path(cfg["paths"]["phrogs_all_hmm"])
    files = sorted(p for ext in ("*.faa", "*.fasta", "*.afa", "*.msa", "*.fma") for p in msa_dir.glob(ext))
    if not files:
        logging.warning(f"No PHROG MSA files in {msa_dir}. Keeping existing HMMs if present "
                        f"(hallmark HMM exists: {out_hall.exists()}).")
        return

    meta = load_annotations(Path(cfg["paths"]["phrogs_annotation_tsv"]))
    alpha = Alphabet.amino()
    builder, bg = Builder(alpha), Background(alpha)
    out_hall.parent.mkdir(parents=True, exist_ok=True)

    n_all = n_hall = n_err = 0
    with open(out_all, "wb") as fa, open(out_hall, "wb") as fh:
        for path in files:
            pid = path.stem.split(".")[0]
            if not pid.startswith("phrog_"):
                pid = f"phrog_{pid}"
            try:
                with MSAFile(str(path), format="afa", digital=True, alphabet=alpha) as mf:
                    msa = mf.read()
                info = meta.get(pid, {"annot": "hypothetical protein", "category": "unknown function"})
                msa.name = pid
                msa.accession = pid
                msa.description = f"{info['category']} | {info['annot']}"
                hmm, _, _ = builder.build_msa(msa, bg)
                hmm.write(fa)
                n_all += 1
                if is_first_pass(info["annot"], info["category"], cfg):
                    hmm.write(fh)
                    n_hall += 1
            except Exception as e:  # noqa: BLE001
                logging.warning(f"Skipping {path.name}: {e}")
                n_err += 1
    logging.info(f"Built {n_all} PHROG HMMs ({n_hall} hallmark/supporting, {n_err} skipped).")
    if n_all == 0:
        sys.exit("No HMMs were built; check MSA format.")


if __name__ == "__main__":
    main()
