#!/usr/bin/env python3
"""
src/04_assess_quality.py
Terminal-repeat (DTR/ITR) analysis with low-complexity filtering, host-flank flag, and optional
CheckV integration. If CheckV is unavailable, completeness is recorded as NA (never 0%).
"""

import logging
import math
import shutil
import subprocess
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO

from common import load_cfg, rdir, setup_logging

_RC = str.maketrans("ACGTacgtN", "TGCAtgcaN")


def revcomp(s):
    return s.translate(_RC)[::-1]


def entropy(seq):
    if not seq:
        return 0.0
    n = len(seq)
    return -sum((c / n) * math.log2(c / n) for c in Counter(seq).values())


def top_base_fraction(seq):
    return max(Counter(seq).values()) / len(seq) if seq else 1.0


def detect_terminal_repeats(seq, min_len, max_len, max_homo, min_ent):
    n = len(seq)
    limit = min(max_len, n // 2)
    dtr = next((k for k in range(limit, min_len - 1, -1) if seq[:k] == seq[n - k:]), 0)
    rc_end = revcomp(seq[n - limit:])
    itr = next((k for k in range(limit, min_len - 1, -1) if seq[:k] == rc_end[:k]), 0)

    for kind, k in (("DTR", dtr), ("ITR", itr)):
        if k == 0:
            continue
        rep = seq[:k]
        ent, homo = entropy(rep), top_base_fraction(rep)
        if ent >= min_ent and homo <= max_homo:
            return {"tr_type": kind, "tr_length_bp": k, "tr_shannon_entropy": round(ent, 3),
                    "closure_status": f"terminal-repeat-supported assembly ({kind})"}
        return {"tr_type": f"low_complexity_{kind}_artifact", "tr_length_bp": k,
                "tr_shannon_entropy": round(ent, 3),
                "closure_status": "linear_or_unclosed (low-complexity end repeat excluded)"}
    return {"tr_type": "none", "tr_length_bp": 0, "tr_shannon_entropy": np.nan,
            "closure_status": "no_terminal_repeat_detected"}


def run_checkv(fna, cfg):
    exe = shutil.which("checkv")
    db = Path(cfg["paths"]["checkv_db"])
    out = Path(cfg["paths"]["results_dir"]) / "checkv"
    if not (exe and db.exists() and any(db.iterdir()) and fna.stat().st_size > 0):
        return {}
    try:
        subprocess.run([exe, "end_to_end", str(fna), str(out), "-d", str(db),
                        "-t", str(cfg["project"]["threads"]), "--restart"],
                       check=True, capture_output=True, text=True)
    except Exception as e:  # noqa: BLE001
        logging.warning(f"CheckV failed ({e}); using pip-only fallback (completeness = NA).")
        return {}
    summ = out / "quality_summary.tsv"
    if not summ.exists():
        return {}
    df = pd.read_csv(summ, sep="\t")
    return {str(r["contig_id"]): {"checkv_quality": r.get("checkv_quality", "NA"),
                                  "completeness_pct": r.get("completeness", np.nan),
                                  "completeness_method": r.get("completeness_method", "NA"),
                                  "contamination_pct": r.get("contamination", np.nan),
                                  "checkv_provirus": r.get("provirus", "NA")}
            for _, r in df.iterrows()}


def main():
    cfg = load_cfg()
    setup_logging("04_assess_quality", cfg)
    tab_dir, seq_dir = rdir(cfg, "tables"), rdir(cfg, "sequences")
    short_tsv, short_fna = tab_dir / "04_shortlist_candidates.tsv", seq_dir / "shortlist_candidates.fna"
    out_tsv = tab_dir / "05_candidate_quality.tsv"

    short = pd.read_csv(short_tsv, sep="\t")
    if short.empty:
        logging.info("No shortlisted candidates; nothing to assess.")
        pd.DataFrame(columns=["contig_id"]).to_csv(out_tsv, sep="\t", index=False)
        return

    cv = run_checkv(short_fna, cfg)
    q = cfg["quality"]
    seqs = {r.id: str(r.seq).upper() for r in SeqIO.parse(short_fna, "fasta")}
    rows = []
    for _, row in short.iterrows():
        cid = row["contig_id"]
        tr = detect_terminal_repeats(seqs[cid], int(q["tr_min_length_bp"]), int(q["tr_max_length_bp"]),
                                     float(q["tr_max_homopolymer_frac"]), float(q["tr_min_shannon_entropy"]))
        c = cv.get(cid, {})
        limits = []
        if not c:
            limits.append("CheckV_not_run_completeness_unknown")
        if tr["tr_type"] == "none":
            limits.append("unclosed_contig_fragment")
        limits.append("read_mapping_junction_validation_not_performed")
        rows.append({"contig_id": cid,
                     "checkv_quality_tier": c.get("checkv_quality", "NA"),
                     "completeness_estimate_pct": c.get("completeness_pct", np.nan),
                     "completeness_method": c.get("completeness_method", "not_measured_pip_fallback"),
                     "contamination_estimate_pct": c.get("contamination_pct", np.nan),
                     "checkv_provirus": c.get("checkv_provirus", "NA"),
                     "host_flank_flag": ("potential_host_flank_detected" if int(row["cellular_housekeeping_genes"]) > 0
                                         else "no_cellular_housekeeping_detected"),
                     **tr, "read_support_status": "NA (assembly_only_not_read_validated)",
                     "quality_limitations": ";".join(limits)})
    pd.DataFrame(rows).to_csv(out_tsv, sep="\t", index=False)
    logging.info(f"Quality assessment complete for {len(rows)} candidates (CheckV used: {bool(cv)}).")


if __name__ == "__main__":
    main()
