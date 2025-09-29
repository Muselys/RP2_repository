#!/usr/bin/env python3
"""
Marker-centric validity analysis.
- marker_inclusivity_exclusivity_pre.tsv
- marker_inclusivity_exclusivity_post.tsv
- genus_leakage_summary_pre.tsv
- genus_leakage_summary_post.tsv
"""

import argparse, os, sys, re
import pandas as pd
import numpy as np
from collections import defaultdict

REQ_COLS = [
    "qseqid","gene_name","true_species",
    "sseqid","bitscore","qcovs","pident","evalue",
    "ref_species"
]

def norm_sp(x):
    if x is None or (isinstance(x,float) and pd.isna(x)): return None
    s = str(x).strip()
    if not s: return None
    s = re.sub(r"\s+"," ", s)
    return s.replace(" ","_")

def genus_of(sp):
    return sp.split("_",1)[0] if sp else None

def read_required_cols(path, chunksize):
    dtype = {"qseqid":str, "sseqid":str, "gene_name":str}
    for chunk in pd.read_csv(path, sep="\t", chunksize=chunksize, dtype=dtype, low_memory=False):
        missing = [c for c in REQ_COLS if c not in chunk.columns]
        if missing:
            raise ValueError(f"{path} is missing columns: {missing}")
        chunk["ref_species"]  = chunk["ref_species"].map(norm_sp)
        chunk["true_species"] = chunk["true_species"].map(norm_sp)
        yield chunk

def collect_hits(path, chunksize=500_000):
    species_refs = defaultdict(set)
    hits_by_q = defaultdict(lambda: defaultdict(set))
    truth_map = {}
    for chunk in read_required_cols(path, chunksize):
        for r in chunk.itertuples(index=False):
            truth_map[r.qseqid] = (r.true_species, r.gene_name)
            if r.ref_species:
                species_refs[r.ref_species].add(r.sseqid)
                hits_by_q[r.qseqid][r.ref_species].add(r.sseqid)
    return species_refs, hits_by_q, truth_map

def marker_table(path, chunksize=500_000):
    species_refs, hits_by_q, truth_map = collect_hits(path, chunksize)
    rows = []
    for q, (true_sp, gname) in truth_map.items():
        target_refs = species_refs.get(true_sp, set())
        non_target_refs = set().union(*[refs for sp, refs in species_refs.items() if sp != true_sp])
        hits = hits_by_q.get(q, {})
        target_hits = hits.get(true_sp, set())
        non_target_hits = set().union(*[refs for sp, refs in hits.items() if sp != true_sp])

        incl = len(target_hits)/len(target_refs) if len(target_refs)>0 else np.nan
        excl = len(non_target_hits)/len(non_target_refs) if len(non_target_refs)>0 else np.nan

        tgt_genus = genus_of(true_sp)
        same_genus_refs = set().union(*[refs for sp, refs in species_refs.items() if genus_of(sp)==tgt_genus and sp!=true_sp])
        cross_genus_refs = set().union(*[refs for sp, refs in species_refs.items() if genus_of(sp)!=tgt_genus])

        same_genus_hits = set().union(*[hits.get(sp,set()) for sp in species_refs if genus_of(sp)==tgt_genus and sp!=true_sp])
        cross_genus_hits = set().union(*[hits.get(sp,set()) for sp in species_refs if genus_of(sp)!=tgt_genus])

        same_leak = len(same_genus_hits)/len(same_genus_refs) if len(same_genus_refs)>0 else np.nan
        cross_leak = len(cross_genus_hits)/len(cross_genus_refs) if len(cross_genus_refs)>0 else np.nan

        rows.append({
            "qseqid": q,
            "gene_name": gname,
            "true_species": true_sp.replace("_"," ") if true_sp else None,
            "inclusivity": incl,
            "exclusivity_leak": excl,
            "same_genus_leak": same_leak,
            "cross_genus_leak": cross_leak,
            "target_refs": len(target_refs),
            "target_refs_hit": len(target_hits),
            "non_target_refs": len(non_target_refs),
            "non_target_refs_hit": len(non_target_hits)
        })
    return pd.DataFrame(rows)

def genus_summary(df):
    return pd.DataFrame({
        "median_inclusivity":[df["inclusivity"].median()],
        "mean_inclusivity":[df["inclusivity"].mean()],
        "median_excl_leak":[df["exclusivity_leak"].median()],
        "mean_excl_leak":[df["exclusivity_leak"].mean()],
        "median_same_genus_leak":[df["same_genus_leak"].median()],
        "mean_same_genus_leak":[df["same_genus_leak"].mean()],
        "median_cross_genus_leak":[df["cross_genus_leak"].median()],
        "mean_cross_genus_leak":[df["cross_genus_leak"].mean()],
        "n_markers":[len(df)]
    })

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pre", required=True)
    ap.add_argument("--post", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--chunksize", type=int, default=500_000)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    print("[INFO] Calculating marker metrics (PRE)...")
    pre_df = marker_table(args.pre, args.chunksize)
    pre_df.to_csv(os.path.join(args.outdir, "marker_inclusivity_exclusivity_pre.tsv"), sep="\t", index=False)
    genus_summary(pre_df).to_csv(os.path.join(args.outdir, "genus_leakage_summary_pre.tsv"), sep="\t", index=False)

    print("[INFO] Calculating marker metrics (POST)...")
    post_df = marker_table(args.post, args.chunksize)
    post_df.to_csv(os.path.join(args.outdir, "marker_inclusivity_exclusivity_post.tsv"), sep="\t", index=False)
    genus_summary(post_df).to_csv(os.path.join(args.outdir, "genus_leakage_summary_post.tsv"), sep="\t", index=False)

    print("[OK] Done.")

if __name__ == "__main__":
    main()

"""
python /data/pam/team230/sm71/scratch/rp2/blast_run/ref_db/results/marker_validity.py \
  --pre /data/pam/team230/sm71/scratch/rp2/blast_run/ref_db/results/blast_results_pre.tsv \
  --post /data/pam/team230/sm71/scratch/rp2/blast_run/ref_db/results/blast_results_post.tsv \
  --outdir /data/pam/team230/sm71/scratch/rp2/blast_run/ref_db/results/marker_analysis \
  --chunksize 200000
"""


