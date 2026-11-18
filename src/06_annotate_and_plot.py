#!/usr/bin/env python3
"""
src/06_annotate_and_plot.py
Annotates shortlisted proteins against the full PHROG profile set (unsupported genes stay
'hypothetical protein') and draws: (1) stage attrition, (2) linear genome maps with GC,
(3) candidate-reference similarity heatmap, (4) dual-threshold PHROG-sharing network.
No synthetic reference nodes are ever added to figures.
"""

import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as mpatches  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import networkx as nx  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pyhmmer  # noqa: E402
import seaborn as sns  # noqa: E402
from Bio import SeqIO  # noqa: E402
from pyhmmer.easel import Alphabet, DigitalSequenceBlock, TextSequence  # noqa: E402
from pyhmmer.plan7 import HMMFile  # noqa: E402

from common import load_cfg, rdir, setup_logging, text  # noqa: E402

COLORS = {"packaging": "#d95f02", "head/capsid": "#7570b3", "tail": "#1b9e77", "replication": "#e7298a",
          "lysis": "#66a61e", "regulation/integration": "#e6ab02", "unknown": "#bdbdbd"}


def functional_group(desc):
    """Descriptions are 'category | annotation'. Use the specific annotation first; fall back to the
    broad category only if the annotation gives no match (so 'head and packaging | major capsid
    protein' is head/capsid, not packaging)."""
    category, _, annot = desc.partition(" | ")
    if not annot:
        category, annot = "", category
    for part in (annot, category):
        grp = _group_of(part.lower())
        if grp != "unknown":
            return grp
    return "unknown"


def _group_of(d):
    groups = [("packaging", ["terminase", "portal", "packaging"]),
              ("head/capsid", ["capsid", "head", "coat", "scaffold"]),
              ("tail", ["tail", "baseplate", "tape measure", "fiber", "neck"]),
              ("replication", ["polymerase", "helicase", "primase", "recombinase", "exonuclease", "replication"]),
              ("lysis", ["holin", "endolysin", "lysin", "spanin", "lysis"]),
              ("regulation/integration", ["integrase", "repressor", "regulator", "transcription"])]
    for name, kws in groups:
        if any(k in d for k in kws):
            return name
    return "unknown"


def plot_retention(cfg, fig_dir, tab_dir):
    prep = pd.read_csv(tab_dir / "01_preprocessing_counts.tsv", sep="\t")
    ev = pd.read_csv(tab_dir / "04_contig_candidate_evidence.tsv", sep="\t")
    labels = ["Raw records", "Retained >=5 kb", "Priority >=10 kb", "Review candidates", "Priority candidates"]
    counts = [int(prep["total_raw_records"].iloc[0]), int(prep["retained_ge5kb"].iloc[0]),
              int(prep["priority_ge10kb"].iloc[0]),
              int((ev["candidate_status"] == "heuristic_review_candidate").sum()),
              int((ev["candidate_status"] == "heuristic_priority_candidate").sum())]
    fig, ax = plt.subplots(figsize=(8, 4.5), dpi=200)
    bars = ax.bar(labels, counts, color=["#6baed6", "#4292c6", "#2171b5", "#fd8d3c", "#d94801"])
    for b in bars:
        ax.text(b.get_x() + b.get_width() / 2, b.get_height(), str(int(b.get_height())),
                ha="center", va="bottom", fontsize=9, fontweight="bold")
    ax.set_ylabel("Contig count")
    ax.set_title("Contig retention across filtering and phage identification stages")
    plt.xticks(rotation=20, ha="right")
    plt.tight_layout()
    plt.savefig(fig_dir / "fig1_stage_retention.png")
    plt.savefig(fig_dir / "fig1_stage_retention.svg")
    plt.close()


def plot_genome_map(cid, seq, genes, fig_dir):
    fig, (axg, axc) = plt.subplots(2, 1, figsize=(11, 4.2), sharex=True, dpi=200,
                                   gridspec_kw={"height_ratios": [1.3, 1]})
    n = len(seq)
    for _, g in genes.iterrows():
        s, e, strand = int(g["start"]), int(g["end"]), int(g["strand"])
        width = max(50, e - s)
        y = 0.55 if strand == 1 else 0.25
        axg.add_patch(mpatches.FancyArrow(
            s if strand == 1 else e, y, width if strand == 1 else -width, 0, width=0.18, head_width=0.28,
            head_length=min(width * 0.35, 250), length_includes_head=True,
            facecolor=COLORS.get(g["functional_group"], "#bdbdbd"), edgecolor="black", linewidth=0.5))
    axg.set_ylim(0, 1)
    axg.set_xlim(0, n)
    axg.set_yticks([0.25, 0.55])
    axg.set_yticklabels(["- strand", "+ strand"])
    axg.set_title(f"Predicted linear genome architecture: {cid} ({n:,} bp) [computational prediction]")
    axg.legend(handles=[mpatches.Patch(color=c, label=k) for k, c in COLORS.items()],
               loc="upper right", ncol=4, fontsize=7.5)
    win, step = 500, 200
    xs = list(range(0, max(1, n - win), step))
    gcs = [(seq[i:i + win].count("G") + seq[i:i + win].count("C")) / len(seq[i:i + win]) * 100 for i in xs]
    if xs:
        axc.plot([x + win // 2 for x in xs], gcs, color="#2b5c8f", linewidth=1.2)
        axc.axhline(np.mean(gcs), color="gray", linestyle="--", linewidth=0.8, label=f"mean GC {np.mean(gcs):.1f}%")
        axc.legend(loc="upper right", fontsize=8)
    axc.set_ylabel("GC (%)")
    axc.set_xlabel("Contig coordinate (bp)")
    plt.tight_layout()
    plt.savefig(fig_dir / f"fig2_genome_map_{cid}.png")
    plt.savefig(fig_dir / f"fig2_genome_map_{cid}.svg")
    plt.close()


def plot_heatmap(tab_dir, fig_dir):
    p = tab_dir / "06_pairwise_similarity_matrix.tsv"
    if not p.exists():
        return
    df = pd.read_csv(p, sep="\t")
    if df.empty:
        logging.info("No aligned candidate-reference pairs; heatmap skipped.")
        return
    piv = df.pivot_table(index="candidate_id", columns="reference_id", values="approx_wgs_pct")
    fig, ax = plt.subplots(figsize=(max(6, 0.9 * piv.shape[1] + 3), max(3.5, 0.6 * piv.shape[0] + 2)), dpi=200)
    sns.heatmap(piv, annot=True, fmt=".1f", cmap="YlGnBu", vmin=0, vmax=100,
                cbar_kws={"label": "Approx. whole-genome similarity (%)"}, ax=ax)
    ax.set_title("Candidate vs. closest references (blank = not aligned)")
    plt.tight_layout()
    plt.savefig(fig_dir / "fig3_similarity_heatmap.png")
    plt.close()


def plot_network(annot, cfg, fig_dir):
    fams = {cid: set(g.loc[g["best_phrog_id"] != "NA", "best_phrog_id"]) for cid, g in annot.groupby("contig_id")}
    if len(fams) < 2:
        logging.info("Fewer than 2 candidates; PHROG network skipped.")
        return
    thresholds = [float(cfg["novelty"]["network_jaccard_primary"]), float(cfg["novelty"]["network_jaccard_secondary"])]
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), dpi=200)
    nodes = list(fams)
    for ax, t in zip(axes, thresholds):
        G = nx.Graph()
        G.add_nodes_from(nodes)
        for i in range(len(nodes)):
            for j in range(i + 1, len(nodes)):
                a, b = fams[nodes[i]], fams[nodes[j]]
                jac = len(a & b) / len(a | b) if (a | b) else 0.0
                if jac >= t:
                    G.add_edge(nodes[i], nodes[j], weight=round(jac, 3))
        nx.draw_networkx(G, pos=nx.spring_layout(G, seed=int(cfg["project"]["random_seed"])), ax=ax,
                         node_color="#d95f02", node_size=550, font_size=8, edge_color="#636363")
        ax.set_title(f"PHROG family sharing (Jaccard >= {t})")
        ax.axis("off")
    plt.tight_layout()
    plt.savefig(fig_dir / "fig4_phrog_network.png")
    plt.close()


def main():
    cfg = load_cfg()
    setup_logging("06_annotate_and_plot", cfg)
    tab_dir, seq_dir, fig_dir = rdir(cfg, "tables"), rdir(cfg, "sequences"), rdir(cfg, "figures")
    plot_retention(cfg, fig_dir, tab_dir)

    short_fna, short_faa = seq_dir / "shortlist_candidates.fna", seq_dir / "shortlist_proteins.faa"
    if not short_fna.exists() or short_fna.stat().st_size == 0:
        logging.info("No shortlisted candidates to annotate.")
        pd.DataFrame().to_csv(tab_dir / "07_shortlist_gene_annotations.tsv", sep="\t", index=False)
        return

    alpha = Alphabet.amino()
    prots = [TextSequence(name=r.id, sequence=str(r.seq).replace("*", "")).digitize(alpha)
             for r in SeqIO.parse(short_faa, "fasta")]
    best, conflicts = {}, {}
    all_hmm = Path(cfg["paths"]["phrogs_all_hmm"])
    if all_hmm.exists() and prots:
        block = DigitalSequenceBlock(alpha, prots)
        with HMMFile(str(all_hmm)) as hf:
            hmms = list(hf)
        for top in pyhmmer.hmmsearch(hmms, block, cpus=int(cfg["project"]["threads"]),
                                     Z=int(cfg["hmmer_search"]["fixed_Z"])):
            pid, plen, pdesc = text(top.query.name), top.query.M, text(top.query.description)
            for hit in top:
                gid = text(hit.name)
                for dom in hit.domains:
                    if dom.i_evalue > float(cfg["hmmer_search"]["max_domain_evalue"]):
                        continue
                    cov = (dom.alignment.hmm_to - dom.alignment.hmm_from + 1) / plen
                    if cov < float(cfg["hmmer_search"]["min_profile_coverage"]):
                        continue
                    cand = {"phrog_id": pid, "desc": pdesc, "score": round(dom.score, 2),
                            "evalue": dom.i_evalue, "cov": round(cov, 4)}
                    if gid not in best or cand["score"] > best[gid]["score"]:
                        if gid in best:
                            conflicts.setdefault(gid, []).append(best[gid]["phrog_id"])
                        best[gid] = cand
                    else:
                        conflicts.setdefault(gid, []).append(pid)

    seqs = {r.id: str(r.seq).upper() for r in SeqIO.parse(short_fna, "fasta")}
    gmap = pd.read_csv(tab_dir / "02_gene_to_contig_map.tsv", sep="\t")
    genes = gmap[gmap["contig_id"].isin(seqs)].copy()
    rows = []
    for _, g in genes.iterrows():
        bh = best.get(g["gene_id"])
        rows.append({**g.to_dict(),
                     "best_phrog_id": bh["phrog_id"] if bh else "NA",
                     "predicted_function": f"putative {bh['desc']}" if bh else "hypothetical protein",
                     "functional_group": functional_group(bh["desc"]) if bh else "unknown",
                     "bit_score": bh["score"] if bh else np.nan,
                     "domain_evalue": bh["evalue"] if bh else np.nan,
                     "profile_coverage": bh["cov"] if bh else np.nan,
                     "secondary_conflicting_profiles": ";".join(conflicts.get(g["gene_id"], [])) or "none"})
    annot = pd.DataFrame(rows)
    annot.to_csv(tab_dir / "07_shortlist_gene_annotations.tsv", sep="\t", index=False)

    short = pd.read_csv(tab_dir / "04_shortlist_candidates.tsv", sep="\t")
    short["_prio"] = (short["candidate_status"] != "heuristic_priority_candidate").astype(int)
    for cid in short.sort_values(["_prio", "length_bp"], ascending=[True, False])["contig_id"].head(3):
        plot_genome_map(cid, seqs[cid], annot[annot["contig_id"] == cid], fig_dir)
    plot_heatmap(tab_dir, fig_dir)
    plot_network(annot, cfg, fig_dir)
    logging.info("Annotation and figures complete.")


if __name__ == "__main__":
    main()
