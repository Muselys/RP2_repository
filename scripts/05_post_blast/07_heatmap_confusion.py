#!/usr/bin/env python3
import argparse, os
import numpy as np
import pandas as pd

# headless matplotlib BEFORE pyplot import
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle  # used for diagonal boxes

# ---------- helpers ----------
def norm_sp(x):
    if pd.isna(x): return None
    s = str(x).strip().replace(" ", "_")
    return s or None

def load_predictions(path):
    usecols = ["qseqid","gene_name","true_species","pred_pre","pred_post"]
    df = pd.read_csv(path, sep="\t", usecols=usecols, low_memory=False)
    for c in ["true_species","pred_pre","pred_post"]:
        df[c] = df[c].map(norm_sp)
    return df

def read_species_list(path):
    with open(path) as fh:
        vals = [norm_sp(l.strip()) for l in fh if l.strip()]
    return [v for v in vals if v]

def save_heatmap(df, title, out_png, vmin=None, vmax=None,
                 rotate_xticks=60, xlabel=None, ylabel=None, cmap="Blues",
                 draw_diag_boxes=False, diag_band=0, diag_color="red", diag_lw=1.8,
                 force_square=False):
    arr = df.values; h, w = arr.shape
    if h == 0 or w == 0:
        fig, ax = plt.subplots(figsize=(6,3))
        ax.text(0.5,0.5,"No data", ha="center", va="center"); ax.axis("off")
        fig.suptitle(title, fontsize=12)
        fig.tight_layout()
        fig.savefig(out_png, dpi=150)
        plt.close(fig)
        return

    fig_w = min(max(8, w * 0.18), 60); fig_h = max(3.5, h * 0.35)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    im = ax.imshow(
        arr,
        aspect=("equal" if force_square else "auto"),
        vmin=vmin, vmax=vmax, cmap=cmap, origin="lower"
    )
    xlabels = [str(s).replace("_", " ") for s in df.columns]
    ylabels = [str(s).replace("_", " ") for s in df.index]
    ax.set_xticks(range(w)); ax.set_xticklabels(xlabels, rotation=rotate_xticks, ha="right", fontsize=8)
    ax.set_yticks(range(h)); ax.set_yticklabels(ylabels, fontsize=8)
    if xlabel: ax.set_xlabel(xlabel)
    if ylabel: ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=12)
    fig.colorbar(im, ax=ax).ax.tick_params(labelsize=8)

    if draw_diag_boxes:
        n = min(h, w)
        for i in range(n):
            for j in range(n):
                if abs(i - j) <= diag_band:
                    ax.add_patch(Rectangle((j-0.5, i-0.5), 1, 1, fill=False,
                                           ec=diag_color, lw=diag_lw))

    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)

# ---------- CLI ----------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--predictions", required=True)
    ap.add_argument("--species-file", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    pred = load_predictions(args.predictions)
    species = read_species_list(args.species_file)

    # fixed display order (drop any not present in species file)
    DESIRED_ORDER = [
        "Enterococcus_faecalis",
        "Enterococcus_faecium",
        "Staphylococcus_argenteus",
        "Staphylococcus_aureus",
        "Staphylococcus_capitis",
        "Staphylococcus_epidermidis",
        "Staphylococcus_haemolyticus",
        "Staphylococcus_pseudintermedius",
        "Staphylococcus_sciuri",
        "Streptococcus_agalactiae",
        "Streptococcus_dysgalactiae",
        "Streptococcus_equi",
        "Streptococcus_mitis",
        "Streptococcus_mutans",
        "Streptococcus_pneumoniae",
        "Streptococcus_pyogenes",
        "Streptococcus_suis",
        "Streptococcus_uberis",
    ]
    species = [s for s in DESIRED_ORDER if s in species]

    for mode, col in (("PRE","pred_pre"), ("POST","pred_post")):
        # build full confusion (including None/off-target) to get correct denominators per true species
        d_all = pred.copy()
        d_all[col] = d_all[col].fillna("None")
        cm_all = pd.crosstab(d_all["true_species"], d_all[col], dropna=False)

        # row totals = total queries per true species
        row_totals = cm_all.sum(axis=1)

        # restrict to 18×18 (targets only)
        cm_counts = cm_all.reindex(index=species, columns=species, fill_value=0)

        # row-normalize by FULL denominators
        denom = row_totals.reindex(species).replace(0, np.nan)
        cm_frac = cm_counts.div(denom, axis=0).fillna(0.0)

        # plot heatmap
        save_heatmap(
            cm_frac,
            f"Confusion ({mode}) — row-normalized",
            os.path.join(args.outdir, f"heat_confusion_{mode.lower()}.png"),
            vmin=0.0, vmax=1.0,
            xlabel="Predicted species", ylabel="True species",
            draw_diag_boxes=True, diag_band=0, force_square=True
        )

if __name__ == "__main__":
    main()

"""
EVAL=/lustre/scratch127/pam/teams/team230/sm71/rp2/blast_run/post_blast/eval

python heatmap_confusion.py \
  --predictions "$EVAL/predictions_pre_post.tsv" \
  --species-file /lustre/scratch127/pam/teams/team230/sm71/rp2/blast_run/post_blast/targets_18.txt \
  --outdir "$EVAL"
"""

"""
Confusion Matrix Info:
"""
#confusion_metrics.py
"""
Overall accuracy pre/post (+ delta)
Top misclassification pairs pre/post
Per-species recall with Δ recall (improved/no change/worse)
"""
#!/usr/bin/env python3
import pandas as pd
import numpy as np
import os, argparse

def load_cm(path):
    # expects TSV with row index = true species, columns = predicted species, values = counts
    df = pd.read_csv(path, sep="\t", index_col=0)
    # enforce numeric
    df = df.apply(pd.to_numeric, errors="coerce").fillna(0).astype(int)
    return df

def overall_accuracy(cm):
    correct = int(np.trace(cm.values))
    total = int(cm.values.sum())
    return (correct / total if total else 0.0), correct, total

def top_confusions(cm, top_n=5):
    df = cm.copy()
    np.fill_diagonal(df.values, 0)
    pairs = (
        df.stack()
          .sort_values(ascending=False)
          .reset_index()
          .rename(columns={"level_0":"true", "level_1":"pred", 0:"count"})
    )
    return pairs.head(top_n)

def per_species_recall(cm):
    rowsum = cm.sum(axis=1)
    diag = pd.Series(np.diag(cm.values), index=cm.index, name="tp")
    out = pd.DataFrame({"tp": diag, "total_true": rowsum})
    out["recall"] = out["tp"] / out["total_true"].replace(0, np.nan)
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pre", required=True)
    ap.add_argument("--post", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    cm_pre  = load_cm(args.pre)
    cm_post = load_cm(args.post)

    # overall accuracy
    acc_pre, tp_pre, tot_pre   = overall_accuracy(cm_pre)
    acc_post, tp_post, tot_post = overall_accuracy(cm_post)

    print(f"Overall accuracy PRE:  {acc_pre:.3f} ({tp_pre}/{tot_pre})")
    print(f"Overall accuracy POST: {acc_post:.3f} ({tp_post}/{tot_post})")
    print(f"Δ accuracy: {acc_post - acc_pre:+.3f}")

    # top off-diagonal pairs
    print("\nTop off-diagonal PRE:")
    print(top_confusions(cm_pre, 5).to_string(index=False))
    print("\nTop off-diagonal POST:")
    print(top_confusions(cm_post, 5).to_string(index=False))

    # species recall changes
    rec_pre  = per_species_recall(cm_pre).add_suffix("_pre")
    rec_post = per_species_recall(cm_post).add_suffix("_post")
    rec = rec_pre.join(rec_post, how="outer").fillna(0)
    rec["delta_recall"] = rec["recall_post"] - rec["recall_pre"]

    # save tables
    rec.reset_index(names="species").to_csv(os.path.join(args.outdir,"species_recall_deltas.tsv"),
                                            sep="\t", index=False)
    top_pre  = top_confusions(cm_pre, 20)
    top_post = top_confusions(cm_post, 20)
    top_pre.to_csv(os.path.join(args.outdir,"top_offdiag_pre.tsv"), sep="\t", index=False)
    top_post.to_csv(os.path.join(args.outdir,"top_offdiag_post.tsv"), sep="\t", index=False)

    # quick console highlights
    improved = rec.sort_values("delta_recall", ascending=False).head(5)
    worse    = rec[rec["delta_recall"] < 0].sort_values("delta_recall").head(5)
    flat     = rec[np.isclose(rec["delta_recall"], 0, atol=1e-9)]

    print("\nMost improved (top 5):")
    print(improved[["recall_pre","recall_post","delta_recall"]].to_string())
    if not flat.empty:
        print("\nNo improvement:")
        print(flat[["recall_pre","recall_post","delta_recall"]].to_string())
    if not worse.empty:
        print("\nGot worse (top 5):")
        print(worse[["recall_pre","recall_post","delta_recall"]].to_string())

if __name__ == "__main__":
    main()

#Get unique queries pre vs post from blast.tsv

python - <<'PY'
import pandas as pd

blast_path = "/data/pam/team230/sm71/scratch/rp2/blast_run/post_blast/blast.tsv"

usecols = ["qseqid","pident","qcovs"]  # matches your header
b = pd.read_csv(blast_path, sep="\t", usecols=usecols, low_memory=False)

pre_unique  = b["qseqid"].nunique()
post_unique = b.query("pident >= 80 and qcovs >= 90")["qseqid"].nunique()

print(f"Unique queries PRE={pre_unique}")
print(f"Unique queries POST={post_unique}")

# (optional) assignments/rows as well:
print(f"Total BLAST rows PRE={len(b)}")
print(f"Total BLAST rows POST={len(b.query('pident >= 80 and qcovs >= 90'))}")
PY
