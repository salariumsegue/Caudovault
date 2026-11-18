#!/usr/bin/env python3
"""
src/03_identify_candidates.py
Calls genes (Pyrodigal, metagenomic mode), screens them against hallmark/supporting PHROG HMMs
(PyHMMER, fixed Z), computes gene-architecture metrics, checks for competing cellular/plasmid
explanations, optionally ingests geNomad scores, and writes the candidate shortlist.
"""

import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyhmmer
import pyrodigal
from Bio import SeqIO
from pyhmmer.easel import Alphabet, DigitalSequenceBlock, TextSequence
from pyhmmer.plan7 import HMMFile

from common import contig_of, kw_match, load_cfg, rdir, setup_logging, text

HOUSEKEEPING = ["ribosomal protein", "rrna", "trna synthetase", "dna gyrase", "ftsz", "secy",
                "rpob", "rpoc", "elongation factor tu", "chaperonin groel", "cellular_housekeeping"]
PLASMID = ["para", "parb", "repa", "plasmid partition", "conjugal transfer", "trag", "moba"]


def classify_profile(desc, cfg):
    for cat, kws in cfg["hmmer_search"]["hallmark_categories"].items():
        if any(kw_match(k, desc) for k in kws):
            return "hallmark", cat
    if any(h in desc.lower() for h in HOUSEKEEPING):
        return "cellular_housekeeping", "housekeeping"
    if any(kw_match(p, desc) for p in PLASMID):
        return "plasmid_associated", "plasmid"
    for supp in cfg["hmmer_search"]["supporting_categories"]:
        if kw_match(supp, desc):
            return "supporting", supp
    return "other", "other"


def architecture(contig_len, genes):
    if not genes:
        return {"gene_count": 0, "coding_density": 0.0, "median_gene_len_bp": 0.0,
                "plus_strand_fraction": 0.0, "max_colinear_run": 0, "median_intergenic_bp": np.nan}
    mask = np.zeros(contig_len, dtype=bool)
    lengths, strands, coords = [], [], []
    for g in genes:
        s0, e0 = max(0, g["start"] - 1), min(contig_len, g["end"])
        mask[s0:e0] = True
        lengths.append(e0 - s0)
        strands.append(g["strand"])
        coords.append((g["start"], g["end"]))
    run = best = 1
    for i in range(1, len(strands)):
        run = run + 1 if strands[i] == strands[i - 1] else 1
        best = max(best, run)
    coords.sort()
    ig = [coords[i][0] - coords[i - 1][1] - 1 for i in range(1, len(coords))]
    return {"gene_count": len(genes),
            "coding_density": round(float(mask.sum()) / contig_len, 4),
            "median_gene_len_bp": round(float(np.median(lengths)), 1),
            "plus_strand_fraction": round(sum(s == 1 for s in strands) / len(strands), 4),
            "max_colinear_run": int(best),
            "median_intergenic_bp": round(float(np.median(ig)), 1) if ig else 0.0}


def main():
    cfg = load_cfg()
    setup_logging("03_identify_candidates", cfg)
    seq_dir, tab_dir = rdir(cfg, "sequences"), rdir(cfg, "tables")

    hmm_path = Path(cfg["paths"]["phrogs_hallmark_hmm"])
    if not hmm_path.exists():
        sys.exit(f"Hallmark HMM file missing: {hmm_path}. Run src/01_build_phrog_hmms.py with PHROG MSAs "
                 "staged (real mode) or src/00_make_sanity_panel.py (sanity mode). "
                 "Refusing to screen without profiles.")

    contigs_fna = seq_dir / "contigs_retained_ge5kb.fna"
    inv = pd.read_csv(tab_dir / "01_sequence_inventory.tsv", sep="\t")
    inv_lookup = inv[inv["internal_id"].notna()].set_index("internal_id").to_dict("index")

    finder = pyrodigal.GeneFinder(meta=bool(cfg["gene_calling"]["meta_mode"]))
    alpha = Alphabet.amino()
    digital, gene_rows, contig_genes, contig_len = [], [], {}, {}

    with open(contigs_fna) as fin, open(seq_dir / "proteins.faa", "w") as ffaa, \
            open(seq_dir / "genes.fna", "w") as ffna, open(seq_dir / "genes.gff3", "w") as fgff:
        fgff.write("##gff-version 3\n")
        for rec in SeqIO.parse(fin, "fasta"):
            cid = rec.id
            contig_len[cid] = len(rec.seq)
            contig_genes[cid] = []
            preds = finder.find_genes(str(rec.seq).encode())
            preds.write_gff(fgff, sequence_id=cid, header=False)
            preds.write_translations(ffaa, sequence_id=cid)
            preds.write_genes(ffna, sequence_id=cid)
            for idx, g in enumerate(preds, 1):
                gid = f"{cid}_{idx}"
                aa = g.translate().replace("*", "")
                digital.append(TextSequence(name=gid, sequence=aa).digitize(alpha))
                info = {"gene_id": gid, "contig_id": cid, "gene_index": idx, "start": g.begin,
                        "end": g.end, "strand": g.strand, "partial_5prime": g.partial_begin,
                        "partial_3prime": g.partial_end, "aa_length": len(aa)}
                gene_rows.append(info)
                contig_genes[cid].append(info)
    logging.info(f"Called {len(gene_rows)} genes on {len(contig_len)} contigs.")
    pd.DataFrame(gene_rows).to_csv(tab_dir / "02_gene_to_contig_map.tsv", sep="\t", index=False)

    hit_rows = []
    if digital:
        block = DigitalSequenceBlock(alpha, digital)
        with HMMFile(str(hmm_path)) as hf:
            hmms = list(hf)
        logging.info(f"Searching {len(hmms)} hallmark/supporting profiles (fixed Z={cfg['hmmer_search']['fixed_Z']}).")
        for top in pyhmmer.hmmsearch(hmms, block, cpus=int(cfg["project"]["threads"]),
                                     Z=int(cfg["hmmer_search"]["fixed_Z"])):
            pid, plen, pdesc = text(top.query.name), top.query.M, text(top.query.description)
            rtype, rcat = classify_profile(pdesc, cfg)
            for hit in top:
                gid = text(hit.name)
                for dom in hit.domains:
                    if dom.i_evalue > float(cfg["hmmer_search"]["max_domain_evalue"]):
                        continue
                    aln = dom.alignment
                    cov = (aln.hmm_to - aln.hmm_from + 1) / plen
                    if cov >= float(cfg["hmmer_search"]["min_profile_coverage"]):
                        hit_rows.append({"contig_id": contig_of(gid), "gene_id": gid, "profile_id": pid,
                                         "profile_description": pdesc, "role_type": rtype,
                                         "role_category": rcat, "bit_score": round(dom.score, 2),
                                         "domain_i_evalue": dom.i_evalue, "profile_coverage": round(cov, 4)})
    cols = ["contig_id", "gene_id", "profile_id", "profile_description", "role_type", "role_category",
            "bit_score", "domain_i_evalue", "profile_coverage"]
    hits = pd.DataFrame(hit_rows, columns=cols)
    hits.to_csv(tab_dir / "03_hmmer_hallmark_hits.tsv", sep="\t", index=False)

    gmap = {}
    gtsv = Path(cfg["paths"]["results_dir"]) / "genomad" / "assembly_virus_summary.tsv"
    if gtsv.exists():
        gdf = pd.read_csv(gtsv, sep="\t")
        gmap = {str(r["seq_name"]).split()[0]: float(r.get("virus_score", np.nan)) for _, r in gdf.iterrows()}
        logging.info(f"Loaded {len(gmap)} geNomad scores.")

    by_contig = {cid: grp for cid, grp in hits.groupby("contig_id")}
    rules = cfg["candidate_rules"]
    rows = []
    for cid, clen in contig_len.items():
        arch = architecture(clen, contig_genes[cid])
        sub = by_contig.get(cid)
        if sub is not None:
            hall = sub[sub["role_type"] == "hallmark"]
            cats = sorted(hall["role_category"].unique())
            hall_genes = set(hall["gene_id"])
            indep = min(len(cats), len(hall_genes))  # one gene cannot satisfy two categories
            supp = set(sub[sub["role_type"] == "supporting"]["gene_id"]) - hall_genes
            hk = set(sub[sub["role_type"] == "cellular_housekeeping"]["gene_id"])
            pl = set(sub[sub["role_type"] == "plasmid_associated"]["gene_id"])
            max_score = float(sub["bit_score"].max())
            summary = ";".join(f"{r['role_category']}:{r['profile_id']}({r['bit_score']})"
                               for _, r in hall.sort_values("bit_score", ascending=False)
                               .drop_duplicates("role_category").iterrows()) or "none"
        else:
            indep, cats, supp, hk, pl, max_score, summary = 0, [], set(), set(), set(), np.nan, "none"

        if len(hk) > int(rules["max_cellular_housekeeping_genes"]):
            status = "excluded_cellular_or_host_island"
        elif indep >= int(rules["priority_min_distinct_hallmarks"]):
            status = "heuristic_priority_candidate"
        elif indep >= int(rules["review_min_hallmarks"]) and len(supp) >= int(rules["review_min_supporting_genes"]):
            status = "heuristic_review_candidate"
        else:
            status = "unresolved"

        meta = inv_lookup.get(cid, {})
        orig = meta.get("original_id", cid)
        rows.append({"contig_id": cid, "original_id": orig, "source_type": meta.get("source_type", "NA"),
                     "length_bp": clen, "gc_percent": meta.get("gc_percent", np.nan),
                     "ambiguous_base_fraction": meta.get("ambiguous_base_fraction", 0.0), **arch,
                     "low_coding_density_flag": arch["coding_density"] < float(rules["min_coding_density_warning"]),
                     "independent_hallmark_categories": indep,
                     "hallmark_categories_detected": ";".join(cats) if cats else "none",
                     "hallmark_evidence_summary": summary, "supporting_phage_genes": len(supp),
                     "cellular_housekeeping_genes": len(hk), "plasmid_associated_genes": len(pl),
                     "max_hallmark_bit_score": max_score,
                     "genomad_virus_score": gmap.get(orig, gmap.get(cid, np.nan)),
                     "candidate_status": status})

    ev = pd.DataFrame(rows)
    ev.to_csv(tab_dir / "04_contig_candidate_evidence.tsv", sep="\t", index=False)
    short = ev[ev["candidate_status"].isin(["heuristic_priority_candidate", "heuristic_review_candidate"])]
    short.to_csv(tab_dir / "04_shortlist_candidates.tsv", sep="\t", index=False)

    keep = set(short["contig_id"])
    with open(contigs_fna) as fin, open(seq_dir / "shortlist_candidates.fna", "w") as fout:
        for rec in SeqIO.parse(fin, "fasta"):
            if rec.id in keep:
                fout.write(f">{rec.id}\n{rec.seq}\n")
    with open(seq_dir / "proteins.faa") as fin, open(seq_dir / "shortlist_proteins.faa", "w") as fout:
        for rec in SeqIO.parse(fin, "fasta"):
            if contig_of(rec.id) in keep:
                SeqIO.write(rec, fout, "fasta")

    logging.info(f"Status counts: {ev['candidate_status'].value_counts().to_dict()}")
    logging.info(f"Shortlisted {len(short)} of {len(ev)} contigs >=5 kb.")


if __name__ == "__main__":
    main()
