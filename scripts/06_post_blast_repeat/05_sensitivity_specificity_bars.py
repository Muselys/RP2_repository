#!/usr/bin/env python3
import argparse, os
import pandas as pd
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def load_metrics(path, tag):
    df = pd.read_csv(path, sep="\t")
    # columns: species, TP, FP, FN, TN, Sensitivity, Specificity
    df["species"] = df["species"].astype(str)
    df = df[["species","Sensitivity","Specificity"]].copy()
    df = df.rename(columns={"Sensitivity":f"Sensitivity_{tag}",
                            "Specificity":f"Specificity_{tag}"})
    return df

def plot_bars(df, metric_cols, title, out_png):
    species = df["species"].tolist()
    x = np.arange(len(species))
    width = 0.38

    fig, ax = plt.subplots(figsize=(max(10, 0.5*len(species)), 6))
    ax.bar(x - width/2, df[metric_cols[0]], width, label=metric_cols[0].split("_")[-1])
    ax.bar(x + width/2, df[metric_cols[1]], width, label=metric_cols[1].split("_")[-1])

    ax.set_xticks(x)
    ax.set_xticklabels(species, rotation=60, ha="right", fontsize=8)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Score")
    ax.set_title(title)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)

def main():
    ap = argparse.ArgumentParser(description="Bar charts for Sensitivity/Specificity per species (pre vs post).")
    ap.add_argument("--pre", required=True)
    ap.add_argument("--post", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    pre = load_metrics(args.pre, "pre")
    post = load_metrics(args.post, "post")
    df = pd.merge(pre, post, on="species", how="outer").fillna(0.0)

    # Sort species for cleaner view (by post sensitivity desc)
    df = df.sort_values("Sensitivity_post", ascending=False)

    # Sensitivity bars
    plot_bars(df, ["Sensitivity_pre","Sensitivity_post"],
              "Sensitivity per species (pre vs post)",
              os.path.join(args.outdir, "bars_sensitivity.png"))

    # Specificity bars
    plot_bars(df, ["Specificity_pre","Specificity_post"],
              "Specificity per species (pre vs post)",
              os.path.join(args.outdir, "bars_specificity.png"))

if __name__ == "__main__":
    main()

"""
python plot_metrics_bars.py \
  --pre  /data/pam/team230/sm71/scratch/rp2/blast_run/ref_db/results/classification_analysis/per_species_metrics_pre.tsv \
  --post /data/pam/team230/sm71/scratch/rp2/blast_run/ref_db/results/classification_analysis/per_species_metrics_post.tsv \
  --outdir /data/pam/team230/sm71/scratch/rp2/blast_run/ref_db/results/classification_analysis

"""