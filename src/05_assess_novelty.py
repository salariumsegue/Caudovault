#!/usr/bin/env python3
"""
src/05_assess_novelty.py
Novelty screen of shortlisted candidates against FROZEN reference snapshots.
  A. Reference de-duplication with alias preservation.
  B. Sourmash k=21/31 sketch prefilter (shared-hash ranking).
  C. Alignment of each candidate to its top-N sketch hits: BLAST+ if installed, otherwise an exact
     31-mer anchor aligner. Reports non-redundant aligned bases, query/reference aligned fractions,
     and an APPROXIMATE whole-genome similarity (VIRIDIC-like; not VIRIDIC).

Honesty rules: with no references staged, novelty is "not_assessed" (never "novel"). Sketch/exact-anchor
methods only see ~>85-90% identity; divergent relatives need protein-level comparison (follow-up).
"""

import gzip
import hashlib
import logging
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import sourmash
from Bio import SeqIO

from common import load_cfg, rdir, setup_logging

_RC = str.maketrans("ACGTN", "TGCAN")


def revcomp(s):
    return s.translate(_RC)[::-1]


def open_maybe_gz(p):
    return gzip.open(p, "rt") if str(p).endswith(".gz") else open(p, "rt")


def load_references(cfg):
    sources = [("RefSeq_Phage", cfg["paths"]["refseq_fasta"]), ("INPHARED", cfg["paths"]["inphared_fasta"]),
               ("Env_Catalog", cfg["paths"]["env_viral_fasta"])]
    refs = {}
    for label, p in sources:
        p = Path(p)
        if not p.exists():
            continue
        with open_maybe_gz(p) as fh:
            for rec in SeqIO.parse(fh, "fasta"):
                s = str(rec.seq).upper()
                h = hashlib.sha256(s.encode()).hexdigest()
                if h in refs:
                    refs[h]["aliases"].append(f"{label}:{rec.id}")
                else:
                    refs[h] = {"ref_id": rec.id, "source_db": label, "aliases": [f"{label}:{rec.id}"],
                               "seq": s, "length_bp": len(s)}
        logging.info(f"Loaded {label}: running total {len(refs)} unique references.")
    if not refs and cfg["novelty"].get("use_control_as_reference"):
        for rec in SeqIO.parse(cfg["paths"]["controls_fasta"], "fasta"):
            if "POS_TAILED_PHAGE" in rec.id:
                s = str(rec.seq).upper()
                refs[hashlib.sha256(s.encode()).hexdigest()] = {
                    "ref_id": f"REF_{rec.id}", "source_db": "Sanity_Control_Ref",
                    "aliases": [f"Sanity_Control_Ref:{rec.id}"], "seq": s, "length_bp": len(s)}
        logging.warning("SANITY MODE: using the synthetic positive control as the only reference.")
    return refs


def merge_len(intervals):
    if not intervals:
        return 0
    iv = sorted((min(a, b), max(a, b)) for a, b in intervals)
    merged = [iv[0]]
    for s, e in iv[1:]:
        if s <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
        else:
            merged.append((s, e))
    return sum(e - s + 1 for s, e in merged)


def _summarize(method, q_len, r_len, nr_q, nr_r, ident):
    wgs = ((nr_q + nr_r) * ident / 100.0) / (q_len + r_len) * 100.0
    return {"method": method, "aligned_identity_pct": round(ident, 2),
            "nonredundant_aligned_bases_query": nr_q,
            "query_aligned_fraction": round(nr_q / q_len, 4) if q_len else 0.0,
            "ref_aligned_fraction": round(nr_r / r_len, 4) if r_len else 0.0,
            "genome_length_ratio": round(min(q_len, r_len) / max(q_len, r_len), 4),
            "approx_wgs_pct": round(wgs, 2)}


def align_blast(q, r, exe):
    with tempfile.TemporaryDirectory() as d:
        qf, rf = Path(d) / "q.fna", Path(d) / "r.fna"
        qf.write_text(f">q\n{q}\n")
        rf.write_text(f">r\n{r}\n")
        res = subprocess.run([exe, "-query", str(qf), "-subject", str(rf), "-task", "dc-megablast",
                              "-outfmt", "6 qstart qend sstart send pident length"],
                             capture_output=True, text=True)
    q_iv, r_iv, wsum, tot = [], [], 0.0, 0
    for line in res.stdout.strip().splitlines():
        p = line.split("\t")
        q_iv.append((int(p[0]), int(p[1])))
        r_iv.append((int(p[2]), int(p[3])))
        wsum += float(p[4]) * int(p[5])
        tot += int(p[5])
    ident = wsum / tot if tot else 0.0
    return _summarize("blastn_dc_megablast", len(q), len(r), merge_len(q_iv), merge_len(r_iv), ident)


def align_kmer(q, r, k=31):
    """Exact k-mer anchors on both strands; identity reported only for anchored bases (=100)."""
    r_pos = {}
    for i in range(len(r) - k + 1):
        r_pos.setdefault(r[i:i + k], i + 1)
    best = (0, 0)
    for qq in (q, revcomp(q)):
        q_iv, r_iv = [], []
        for i in range(len(qq) - k + 1):
            rp = r_pos.get(qq[i:i + k])
            if rp:
                q_iv.append((i + 1, i + k))
                r_iv.append((rp, rp + k - 1))
        nq, nr = merge_len(q_iv), merge_len(r_iv)
        if nq > best[0]:
            best = (nq, nr)
    return _summarize("exact_31mer_anchor_fallback", len(q), len(r), best[0], best[1], 100.0 if best[0] else 0.0)


def main():
    cfg = load_cfg()
    setup_logging("05_assess_novelty", cfg)
    tab_dir, seq_dir = rdir(cfg, "tables"), rdir(cfg, "sequences")
    short_fna = seq_dir / "shortlist_candidates.fna"
    out_n, out_m = tab_dir / "06_novelty_comparisons.tsv", tab_dir / "06_pairwise_similarity_matrix.tsv"
    if not short_fna.exists() or short_fna.stat().st_size == 0:
        logging.info("No shortlisted candidates; skipping novelty.")
        pd.DataFrame(columns=["contig_id"]).to_csv(out_n, sep="\t", index=False)
        pd.DataFrame(columns=["candidate_id"]).to_csv(out_m, sep="\t", index=False)
        return

    cands = {r.id: str(r.seq).upper() for r in SeqIO.parse(short_fna, "fasta")}
    refs = load_references(cfg)
    nv = cfg["novelty"]
    scaled, top_n = int(nv["sourmash_scaled"]), int(nv["align_top_n"])

    def sk(seq):
        out = []
        for k in nv["sourmash_k"]:
            mh = sourmash.MinHash(n=0, ksize=int(k), scaled=scaled)
            mh.add_sequence(seq, force=True)
            out.append(mh)
        return out

    logging.info(f"Sketching {len(refs)} references (k={nv['sourmash_k']}, scaled={scaled}).")
    ref_sk = {h: sk(r["seq"]) for h, r in refs.items()}
    exe = shutil.which("blastn")
    logging.info(f"Aligner: {'blastn dc-megablast' if exe else 'exact 31-mer anchor fallback (no BLAST+)'}")

    rows, pairs = [], []
    for cid, cseq in cands.items():
        base = {"contig_id": cid}
        if not refs:
            rows.append({**base, "closest_reference_accession": "none", "closest_reference_db": "none",
                         "closest_reference_aliases": "none", "novelty_interpretation": "not_assessed_no_references_staged"})
            continue
        csk = sk(cseq)
        ranked = []
        for h, rs in ref_sk.items():
            shared = sum(c.count_common(r) for c, r in zip(csk, rs))
            if shared > 0:
                ranked.append((shared, h))
        ranked.sort(reverse=True)

        best = None
        for shared, h in ranked[:top_n]:
            r = refs[h]
            aln = align_blast(cseq, r["seq"], exe) if exe else align_kmer(cseq, r["seq"])
            j = [round(c.jaccard(x), 4) for c, x in zip(csk, ref_sk[h])]
            pairs.append({"candidate_id": cid, "reference_id": r["ref_id"], "reference_db": r["source_db"],
                          "sourmash_jaccard_k21": j[0], "sourmash_jaccard_k31": j[-1],
                          "shared_hashes_total": shared, "approx_wgs_pct": aln["approx_wgs_pct"],
                          "query_aligned_fraction": aln["query_aligned_fraction"]})
            if best is None or aln["approx_wgs_pct"] > best[1]["approx_wgs_pct"]:
                best = (r, aln, j)

        if best is None:
            rows.append({**base, "closest_reference_accession": "none", "closest_reference_db": "none",
                         "closest_reference_aliases": "none", "sourmash_shared_hashes": 0,
                         "novelty_interpretation": "no_sketch_or_anchor_match_in_references_provisional_candidate"})
            continue
        r, aln, j = best
        w = aln["approx_wgs_pct"]
        if w >= float(nv["species_wgs_threshold"]):
            interp = "known_species_level_match_or_rediscovery"
        elif w >= float(nv["genus_wgs_threshold"]):
            interp = "provisional_novel_species_within_known_genus_threshold"
        else:
            interp = "no_close_match_in_frozen_references_provisional_candidate"
        rows.append({**base, "closest_reference_accession": r["ref_id"], "closest_reference_db": r["source_db"],
                     "closest_reference_aliases": ";".join(r["aliases"]),
                     "sourmash_jaccard_k21": j[0], "sourmash_jaccard_k31": j[-1], **aln,
                     "novelty_interpretation": interp})

    pd.DataFrame(rows).to_csv(out_n, sep="\t", index=False)
    pd.DataFrame(pairs, columns=["candidate_id", "reference_id", "reference_db", "sourmash_jaccard_k21",
                                 "sourmash_jaccard_k31", "shared_hashes_total", "approx_wgs_pct",
                                 "query_aligned_fraction"]).to_csv(out_m, sep="\t", index=False)
    logging.info(f"Novelty comparisons written for {len(rows)} candidates.")


if __name__ == "__main__":
    main()
