"""Merge geNomad output onto the candidate table. Adds host/prophage screening.
Usage: python src/08_merge_genomad.py <genomad_out_dir> [results_dir]
Verdicts are screening labels, not confirmation."""
import sys, glob, os
import pandas as pd

gdir = sys.argv[1]; res = sys.argv[2] if len(sys.argv) > 2 else "results"
def one(pat):
    f = glob.glob(os.path.join(gdir, pat), recursive=True)
    if not f: sys.exit(f"missing geNomad file: {pat}")
    return f[0]
rd = lambda p: pd.read_csv(p, sep="\t", keep_default_na=False)
agg = rd(one("**/*_aggregated_classification.tsv"))
genes = rd(one("**/*_genes.tsv"))
try: vs = rd(one("**/*_virus_summary.tsv"))
except SystemExit: vs = pd.DataFrame(columns=["seq_name","virus_score","n_hallmarks","topology"])

genes["contig_id"] = genes["gene"].str.rsplit("_", n=1).str[0]
genes["uscg"] = (~genes["uscg"].astype(str).str.strip().isin(["", "NA", "0", "0.0", "nan"])).astype(int)
uscg = genes.groupby("contig_id")["uscg"].sum().rename("genomad_uscg_genes")
vs["contig_id"] = vs["seq_name"].str.split("|").str[0]
vs["is_provirus"] = vs["seq_name"].str.contains("|provirus", regex=False)
vsum = vs.groupby("contig_id").agg(genomad_virus_score=("virus_score", "max"),
        genomad_provirus_calls=("is_provirus", "sum"))
agg = agg.rename(columns={"seq_name": "contig_id", "chromosome_score": "genomad_chromosome_score",
                          "plasmid_score": "genomad_plasmid_score", "virus_score": "genomad_agg_virus_score"})
keep = ["contig_id", "genomad_chromosome_score", "genomad_plasmid_score", "genomad_agg_virus_score"]
cand = pd.read_csv(os.path.join(res, "tables/main_candidate_results.tsv"), sep="\t", keep_default_na=False)
out = (cand.merge(agg[keep], on="contig_id", how="left")
           .merge(uscg, on="contig_id", how="left").merge(vsum, on="contig_id", how="left"))
for c in ["genomad_uscg_genes", "genomad_provirus_calls"]: out[c] = pd.to_numeric(out[c], errors="coerce").fillna(0)
for c in ["genomad_chromosome_score","genomad_agg_virus_score","genomad_virus_score"]:
    out[c] = pd.to_numeric(out[c], errors="coerce")

def verdict(r):
    if pd.isna(r.genomad_agg_virus_score): return "not_scored_by_genomad"
    if r.genomad_uscg_genes >= 2: return "likely_cellular_or_host_island"
    if r.genomad_provirus_calls > 0: return "contains_provirus_region_check_boundaries"
    if r.genomad_agg_virus_score >= 0.8 and r.genomad_chromosome_score < 0.5: return "virus_supported_by_genomad"
    return "weak_or_ambiguous_virus_signal"
out["host_screen_verdict"] = out.apply(verdict, axis=1)
p = os.path.join(res, "tables/08_candidates_with_host_screen.tsv")
out.to_csv(p, sep="\t", index=False)
print(out["host_screen_verdict"].value_counts().to_string())
print("wrote", p)
