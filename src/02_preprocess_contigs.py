#!/usr/bin/env python3
"""
src/02_preprocess_contigs.py
Validates the input FASTA (sample assembly, or the synthetic control panel in sanity mode),
excludes empty/duplicate/invalid/high-N records (each logged with a reason), de-duplicates
exact and reverse-complement copies, and splits >=5 kb (priority >=10 kb) from <5 kb.
"""

import gzip
import hashlib
import logging
import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from common import load_cfg, rdir, setup_logging

VALID = set("ACGTNRYSWKMBDHV")
_RC = str.maketrans("ACGTNRYSWKMBDHV", "TGCANYRSWMKVHDB")


def open_fasta(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def canonical_hash(seq, include_rc=True):
    if include_rc:
        seq = min(seq, seq.translate(_RC)[::-1])
    return hashlib.sha256(seq.encode("ascii")).hexdigest()


def main():
    cfg = load_cfg()
    setup_logging("02_preprocess_contigs", cfg)

    mode = cfg["run"]["mode"]
    if mode == "sanity":
        in_path, label = Path(cfg["paths"]["controls_fasta"]), "sanity_control"
    else:
        in_path, label = Path(cfg["paths"]["raw_contigs"]), "sample"
    if not in_path.exists():
        sys.exit(f"Input FASTA not found: {in_path} (mode={mode}). "
                 "Stage the assembly in data/raw/ or run `bash run_all.sh sanity`.")

    pp = cfg["preprocessing"]
    min_len, prio_len = int(pp["min_contig_length_bp"]), int(pp["priority_contig_length_bp"])
    max_n, dedup_rc = float(pp["max_ambiguous_fraction"]), bool(pp["deduplicate_rc"])

    seq_dir, tab_dir = rdir(cfg, "sequences"), rdir(cfg, "tables")
    counts = {"mode": mode, "total_raw_records": 0, "excluded_corrupt_or_duplicate": 0,
              "archived_short_lt5kb": 0, "retained_ge5kb": 0, "priority_ge10kb": 0}
    seen_ids, seen_hash = set(), {}
    inventory, exclusions = [], []

    def exclude(oid, length, reason, matched="NA"):
        exclusions.append({"original_id": oid, "source": label, "length_bp": length,
                           "reason": reason, "matched_id": matched})
        counts["excluded_corrupt_or_duplicate"] += 1

    n_ctg = 0
    with open_fasta(in_path) as fin, \
            open(seq_dir / "contigs_retained_ge5kb.fna", "w") as fkeep, \
            open(seq_dir / "contigs_short_lt5kb.fna", "w") as fshort:
        for rec in SeqIO.parse(fin, "fasta"):
            counts["total_raw_records"] += 1
            oid, seq = rec.id, str(rec.seq).upper()
            length = len(seq)

            if oid in seen_ids:
                exclude(oid, length, "duplicate_id", oid)
                continue
            seen_ids.add(oid)
            if length == 0:
                exclude(oid, 0, "empty_sequence")
                continue
            if not set(seq) <= VALID:
                exclude(oid, length, "invalid_characters")
                continue

            a, c, g, t = (seq.count(x) for x in "ACGT")
            acgt = a + c + g + t
            n_frac = (length - acgt) / length
            gc = (g + c) / acgt * 100 if acgt else 0.0
            if n_frac > max_n:
                exclude(oid, length, f"ambiguous_fraction_gt_{max_n}")
                continue

            h = canonical_hash(seq, dedup_rc)
            if h in seen_hash:
                exclude(oid, length, "exact_or_rc_duplicate", seen_hash[h])
                continue
            seen_hash[h] = oid

            if length < min_len:
                fshort.write(f">{oid} source={label}\n{seq}\n")
                counts["archived_short_lt5kb"] += 1
                internal, tier = "NA", "archived_lt5kb"
            else:
                n_ctg += 1
                internal = f"CTG_{n_ctg:06d}"
                tier = "priority_ge10kb" if length >= prio_len else "screened_5kb_to_10kb"
                counts["retained_ge5kb"] += 1
                counts["priority_ge10kb"] += int(length >= prio_len)
                fkeep.write(f">{internal} original_id={oid} source={label} length={length}\n{seq}\n")

            inventory.append({"internal_id": internal, "original_id": oid, "source_type": label,
                              "length_bp": length, "gc_percent": round(gc, 3),
                              "ambiguous_base_fraction": round(n_frac, 6),
                              "canonical_sha256": h, "length_tier": tier})

    pd.DataFrame(inventory).to_csv(tab_dir / "01_sequence_inventory.tsv", sep="\t", index=False)
    pd.DataFrame(exclusions, columns=["original_id", "source", "length_bp", "reason", "matched_id"]
                 ).to_csv(tab_dir / "01_exclusion_log.tsv", sep="\t", index=False)
    pd.DataFrame([counts]).to_csv(tab_dir / "01_preprocessing_counts.tsv", sep="\t", index=False)
    logging.info(f"Preprocessing summary: {counts}")


if __name__ == "__main__":
    main()
