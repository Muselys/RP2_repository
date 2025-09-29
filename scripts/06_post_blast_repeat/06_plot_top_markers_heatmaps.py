#!/usr/bin/env python3
import argparse, os, re
import pandas as pd
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

METRICS = ["inclusivity","exclusivity_leak","same_genus_leak","cross_genus_leak"]

def norm_sp(x):
    if x is None or (isinstance(x,float) and pd.isna(x)): return None
    s = str(x).strip()
    s = re.sub(r"\s+"," ", s)
    return s.replace(" ","_") if s else None

def main():
    ap = argparse.ArgumentParser(description="Top-50 markers per species (POST) heatmaps.")
    ap.add_argument("--post_marker_tsv", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--topn", type=int, default=50)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    df = pd.read_csv(args.post_marker_tsv, sep="\t")
    # Expect columns: qseqid, gene_name, true_species, inclusivity, exclusivity_leak, same_genus_leak, cross_genus_leak, ...
    # Normalize species to underscores (the file likely has spaces)
    df["true_species_norm"] = df["true_species"].map(norm_sp)

    # ensure metric numeric
    for m in METRICS:
        df[m] = pd.to_numeric(df[m], errors="coerce")

    species_list = sorted(df["true_species_norm"].dropna().unique())

    for sp in species_list:
        sub = df[df["true_species_norm"] == sp].copy()

        # rank: high inclusivity, then low cross-genus leak, then low exclusivity leak, then low same-genus leak
        sub = sub.sort_values(
            by=["inclusivity","cross_genus_leak","exclusivity_leak","same_genus_leak"],
            ascending=[False, True, True, True]
        ).head(args.topn)

        if sub.empty:
            continue

        mat = sub[METRICS].copy()
        # clip to [0,1] range for color stability
        mat = mat.clip(lower=0.0, upper=1.0)
        arr = mat.values

        # plot
        h, w = arr.shape
        fig_h = max(4, 0.3 * h)
        fig_w = 6
        fig, ax = plt.subplots(figsize=(fig_w, fig_h))
        im = ax.imshow(arr, origin="lower", aspect="auto", vmin=0.0, vmax=1.0, cmap="viridis")
        ax.set_yticks(range(h))
        # marker label: prefer gene_name then qseqid
        labels = []
        for _, r in sub.iterrows():
            g = str(r.get("gene_name") or "").strip()
            q = str(r.get("qseqid") or "").strip()
            labels.append(g if g else q)
        ax.set_yticklabels(labels, fontsize=7)
        ax.set_xticks(range(w))
        ax.set_xticklabels(["Inclusivity","Excl. leak","Same-genus leak","Cross-genus leak"], rotation=0, fontsize=8)
        ax.set_title(f"Top {len(sub)} markers — {sp.replace('_',' ')} (POST)", fontsize=11)
        fig.colorbar(im, ax=ax).ax.tick_params(labelsize=8)
        fig.tight_layout()

        out_png = os.path.join(args.outdir, f"heat_top_markers_{sp}.png")
        fig.savefig(out_png, dpi=150)
        plt.close(fig)

if __name__ == "__main__":
    main()

"""
python plot_top_markers_heatmaps.py \
  --post_marker_tsv /data/pam/team230/sm71/scratch/rp2/blast_run/ref_db/results/marker_analysis/marker_inclusivity_exclusivity_post.tsv \
  --outdir /data/pam/team230/sm71/scratch/rp2/blast_run/ref_db/results/marker_analysis \
  --topn 50
"""