#!/usr/bin/env python3
import argparse, os
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

def load_cm(path):
    df = pd.read_csv(path, sep="\t", index_col=0)
    df = df.apply(pd.to_numeric, errors="coerce").fillna(0).astype(int)
    return df

def save_heatmap(df, title, out_png, vmin=0.0, vmax=1.0, cmap="Blues"):
    arr = df.values
    h, w = arr.shape
    fig_w = min(max(8, w * 0.5), 40)
    fig_h = min(max(6, h * 0.5), 40)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    im = ax.imshow(arr, origin="lower", aspect="auto", vmin=vmin, vmax=vmax, cmap=cmap)
    ax.set_xticks(range(w)); ax.set_xticklabels([c.replace("_"," ") for c in df.columns], rotation=60, ha="right", fontsize=8)
    ax.set_yticks(range(h)); ax.set_yticklabels([i.replace("_"," ") for i in df.index], fontsize=8)
    ax.set_xlabel("Predicted species"); ax.set_ylabel("True species")
    ax.set_title(title, fontsize=12)
    # red boxes on diagonal
    n = min(h, w)
    for i in range(n):
        ax.add_patch(Rectangle((i-0.5, i-0.5), 1, 1, fill=False, ec="red", lw=1.5))
    cbar = fig.colorbar(im, ax=ax)
    cbar.ax.tick_params(labelsize=8)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)

def row_normalize(cm):
    rs = cm.sum(axis=1)
    return cm.div(rs.replace(0, np.nan), axis=0).fillna(0.0)

def main():
    ap = argparse.ArgumentParser(description="Confusion matrix heatmaps (row-normalized).")
    ap.add_argument("--pre", required=True)
    ap.add_argument("--post", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    cm_pre = load_cm(args.pre)
    cm_post = load_cm(args.post)

    # Row-normalize so each true species sums to 1
    rn_pre = row_normalize(cm_pre)
    rn_post = row_normalize(cm_post)

    save_heatmap(rn_pre,  "Confusion Matrix of Species Classification (Pre-filtering, Row-normalised)",  os.path.join(args.outdir, "heatmap_confusion_pre.png"))
    save_heatmap(rn_post, "Confusion Matrix of Species Classification (Post-filtering, Row-normalised)", os.path.join(args.outdir, "heatmap_confusion_post.png"))

if __name__ == "__main__":
    main()


"""



"""

python plot_marker_passfail.py \
  --pre  marker_inclusivity_exclusivity_pre.tsv \
  --post marker_inclusivity_exclusivity_post.tsv \
  --out  bar_marker_passcounts.png
