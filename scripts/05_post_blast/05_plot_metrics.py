#!/usr/bin/env python3
"""
What this makes:
  1) Confusion heatmaps (1 PRE, 1 POST)
     - x = Predicted species, y = True species
  2) Sensitivity/Specificity
     - TSV tables (pre/post)
     - Bar charts (one PRE, one POST). Species on x, bars = Sens & Spec (0..1).
  3) Per-species candidate heatmaps (for each of the 18 species, PRE & POST)
     - counts: rows = species (pred), cols = that target’s top genes (up to 300)
     - fractions: same matrix column-normalized (per-gene)
     - red border on cells where the gene is ONLY present for the target species

Inputs:
  - predictions_pre_post.tsv (from metrics_top1.py)
  - --species-file (18 lines; spaces or underscores ok)

Outputs:
  outdir/
    heat_confusion_pre.png
    heat_confusion_post.png
    metrics_pre.tsv, metrics_post.tsv
    metrics_pre_bars.png, metrics_post_bars.png
    per_species/<species>_{pre,post}_{counts,frac}.{png,tsv}
"""

import argparse, os
import numpy as np
import pandas as pd

# headless matplotlib (cluster-safe) BEFORE pyplot import
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

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
        vmin=vmin,
        vmax=vmax,
        cmap=cmap,
        origin="lower"
    )
    xlabels = [str(s).replace("_", " ") for s in df.columns]
    ylabels = [str(s).replace("_", " ") for s in df.index]
    ax.set_xticks(range(w)); ax.set_xticklabels(xlabels, rotation=rotate_xticks, ha="right", fontsize=8)
    ax.set_yticks(range(h)); ax.set_yticklabels(ylabels, fontsize=8)
    if xlabel: ax.set_xlabel(xlabel)
    if ylabel: ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=12)
    fig.colorbar(im, ax=ax).ax.tick_params(labelsize=8)

    # NEW: draw diagonal boxes
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

def save_heatmap_with_boxes(df, title, out_png, box_mask,
                            vmin=None, vmax=None, rotate_xticks=60,
                            xlabel=None, ylabel=None, cmap="Blues",
                            force_square=False):
    arr = df.values; h, w = arr.shape
    if h == 0 or w == 0:
        fig, ax = plt.subplots(figsize=(6,3))
        ax.text(0.5, 0.5, "No data", ha="center", va="center"); ax.axis("off")
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

    # Pretty labels (underscores -> spaces)
    xlabels = [str(s).replace("_", " ") for s in df.columns]
    ylabels = [str(s).replace("_", " ") for s in df.index]
    ax.set_xticks(range(w)); ax.set_xticklabels(xlabels, rotation=rotate_xticks, ha="right", fontsize=8)
    ax.set_yticks(range(h)); ax.set_yticklabels(ylabels, fontsize=8)

    if xlabel: ax.set_xlabel(xlabel)
    if ylabel: ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=12)
    fig.colorbar(im, ax=ax).ax.tick_params(labelsize=8)

    # Draw red boxes where box_mask is True
    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            if box_mask[i, j]:
                ax.add_patch(Rectangle((j-0.5, i-0.5), 1, 1, fill=False, ec="red", lw=1.5))

    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


# ---------- confusion + metrics ----------
def confusion_species_only(pred, species, col_pred):
    # keep only rows where predicted species is one of the 18 targets
    d = pred[pred[col_pred].isin(species)].copy()
    # counts: rows=true species, cols=predicted species, both strictly the 18
    cm = pd.crosstab(d["true_species"], d[col_pred], dropna=False)
    cm = cm.reindex(index=species, columns=species, fill_value=0)
    return cm

def row_normalize(cm: pd.DataFrame) -> pd.DataFrame:
    row_sums = cm.sum(axis=1).replace(0, np.nan)
    return (cm.div(row_sums, axis=0)).fillna(0.0)

def accuracy_diagonal(cm_frac: pd.DataFrame) -> pd.DataFrame:
    # assumes columns == index == species in same order
    acc = np.diag(cm_frac.to_numpy())
    n = cm_frac.index.map(lambda s: cm_frac.loc[s].sum())  # should be 1.0 if normalized
    # also add counts to give context
    counts = cm_frac.index.map(lambda s: int(cm_frac.loc[s].to_numpy().sum() * 0 + cm_frac.loc[s].to_numpy().size))  # dummy to keep shape
    # better: pull counts from the unnormalized cm you used to normalize:
    return pd.DataFrame({"species": cm_frac.index, "accuracy": acc})


def metrics_table(pred, species, col_pred):
    ts = pred["true_species"].values
    ps = pred[col_pred].values
    rows=[]
    for S in species:
        t = (ts == S); p = (ps == S)
        TP = int(np.sum(t & p)); FN = int(np.sum(t & ~p))
        FP = int(np.sum(~t & p)); TN = int(np.sum(~t & ~p))
        sens = TP/(TP+FN) if (TP+FN)>0 else np.nan
        spec = TN/(TN+FP) if (TN+FP)>0 else np.nan
        rows.append(dict(species=S, TP=TP, FN=FN, FP=FP, TN=TN,
                         Sensitivity=sens, Specificity=spec))
    return pd.DataFrame(rows).set_index("species")

def metrics_barplot(mt, title, out_png):
    # mt: index=species, columns include Sensitivity, Specificity
    species = mt.index.tolist()
    x = np.arange(len(species))
    width = 0.42
    fig_w = max(10, len(species)*0.6); fig_h = 4.5
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.bar(x - width/2, mt["Sensitivity"].values, width)
    ax.bar(x + width/2, mt["Specificity"].values, width)
    ax.set_xticks(x); ax.set_xticklabels(species, rotation=60, ha="right", fontsize=8)
    ax.set_ylim(0, 1.0)
    ax.set_ylabel("Score (0–1)")
    ax.set_title(title)
    ax.legend(["Sensitivity", "Specificity"], fontsize=9)
    fig.tight_layout(); fig.savefig(out_png, dpi=150); plt.close(fig)

# ---------- per-species candidate heatmaps ----------
def build_gene_set_for_target(pred, target, col_pred, genes_top, max_genes, min_gene_count):
    d = pred.dropna(subset=["gene_name"]).copy()
    # rank genes by how often they are predicted AS the target (in this mode)
    top_as_target = (d[d[col_pred] == target]
                       .groupby("gene_name")["qseqid"].size()
                       .sort_values(ascending=False))
    genes = top_as_target.index.tolist()
    if genes_top: genes = genes[:genes_top]
    genes = genes[:max_genes]
    if not genes: return []
    # optional: keep genes with enough total hits (across all species in this mode)
    counts_all = (d[d["gene_name"].isin(genes)]
                    .groupby(["gene_name", col_pred])["qseqid"].size()
                    .reset_index(name="n"))
    tot = counts_all.groupby("gene_name")["n"].sum()
    genes = [g for g in genes if tot.get(g, 0) >= min_gene_count]
    return genes

def species_by_genes_matrix(pred, species_list, genes, col_pred):
    if not genes:
        return pd.DataFrame(index=species_list, columns=[], dtype=int)
    d = pred.dropna(subset=["gene_name"]).copy()
    d = d[d["gene_name"].isin(genes)]
    grp = d.groupby([col_pred, "gene_name"])["qseqid"].size().reset_index(name="n")
    M = grp.pivot(index=col_pred, columns="gene_name", values="n").fillna(0).astype(int)
    M = M.reindex(index=species_list, columns=genes, fill_value=0)
    return M  # rows = predicted species, cols = genes

def unique_mask_for_target(counts_mat, target):
    present = (counts_mat > 0).astype(int)
    col_nz = present.sum(axis=0)  # per gene
    mask = np.zeros_like(present.values, dtype=bool)
    if target in present.index:
        trow = present.index.get_loc(target)
        for j, g in enumerate(present.columns):
            if col_nz.iloc[j] == 1 and present.iloc[trow, j] == 1:
                mask[trow, j] = True
    return mask

# ---------- CLI ----------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--predictions", required=True)
    ap.add_argument("--species-file", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--genes-top", type=int, default=300, help="max top genes per target (x-axis)")
    ap.add_argument("--max-genes", type=int, default=300, help="hard cap")
    ap.add_argument("--min-gene-count", type=int, default=1, help="min total hits per gene to keep")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    per_dir = os.path.join(args.outdir, "per_species"); os.makedirs(per_dir, exist_ok=True)

    pred = load_predictions(args.predictions)
    species = read_species_list(args.species_file)

    # Lock the axis order exactly as requested (x & y will match this order)
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
    # Use the desired order; (optionally) drop anything not in predictions
    species = [s for s in DESIRED_ORDER if s in species]


    for mode, col in (("PRE","pred_pre"), ("POST","pred_post")):
        # 1.1) Full confusion (incl. None + off-target) for correct denominators
        d_all = pred.copy()
        d_all[col] = d_all[col].fillna("None")
        cm_all = pd.crosstab(d_all["true_species"], d_all[col], dropna=False)

        # row totals = total queries per true species (what you want as denominator)
        row_totals = cm_all.sum(axis=1)

        # 1.2) Species-only 18×18 counts (drop None + off-target columns, keep only target rows)
        cm_counts = cm_all.reindex(index=species, columns=species, fill_value=0)
        cm_counts.to_csv(os.path.join(args.outdir, f"confusion_{mode.lower()}_counts.tsv"), sep="\t")

        # 1.3) Row-normalized fractions in [0,1], using FULL row totals
        denom = row_totals.reindex(species).replace(0, np.nan)
        cm_frac = cm_counts.div(denom, axis=0).fillna(0.0)
        cm_frac.to_csv(os.path.join(args.outdir, f"confusion_{mode.lower()}_fraction.tsv"),
                    sep="\t", float_format="%.6f")

        # 1.4) Per-species accuracy (diagonal) + n_queries for context
        acc = pd.DataFrame({
            "species": species,
            "accuracy": np.diag(cm_frac.to_numpy()),
            "n_queries": denom.fillna(0).astype(int).values,
        })
        acc.to_csv(os.path.join(args.outdir, f"per_species_accuracy_{mode.lower()}.tsv"),
                sep="\t", index=False, float_format="%.6f")

        # 1.5) Heatmap from fractions (0..1)
        save_heatmap(
            cm_frac,
            f"Confusion ({mode}) — row-normalized",
            os.path.join(args.outdir, f"heat_confusion_{mode.lower()}.png"),
            vmin=0.0, vmax=1.0,
            xlabel="Predicted species", ylabel="True species",
            draw_diag_boxes=True, diag_band=0, force_square=True
        )



    # 2) Sensitivity/Specificity (tables + bar charts)
    mt_pre  = metrics_table(pred, species, "pred_pre")
    mt_post = metrics_table(pred, species, "pred_post")
    mt_pre.to_csv( os.path.join(args.outdir, "metrics_pre.tsv"),  sep="\t")
    mt_post.to_csv(os.path.join(args.outdir, "metrics_post.tsv"), sep="\t")
    metrics_barplot(mt_pre,  "Sensitivity & Specificity (PRE)",  os.path.join(args.outdir, "metrics_pre_bars.png"))
    metrics_barplot(mt_post, "Sensitivity & Specificity (POST)", os.path.join(args.outdir, "metrics_post_bars.png"))

    # 3) Per-species candidate heatmaps (PRE + POST)
    for mode, col in (("pre","pred_pre"), ("post","pred_post")):
        for target in species:
            genes = build_gene_set_for_target(pred, target, col,
                                              genes_top=args.genes_top,
                                              max_genes=args.max_genes,
                                              min_gene_count=args.min_gene_count)
            counts = species_by_genes_matrix(pred, species, genes, col)
            # fractions: per-gene (column) normalize
            if counts.shape[1] > 0:
                colsum = counts.sum(axis=0).replace(0, np.nan)
                frac = (counts / colsum).fillna(0.0)
            else:
                frac = counts.astype(float)

            # red border on genes unique to the target species
            mask = unique_mask_for_target(counts, target)

            # save TSVs
            base = f"{target}_{mode}"
            counts.to_csv(os.path.join(per_dir, f"{base}_counts.tsv"), sep="\t")
            frac.to_csv(  os.path.join(per_dir, f"{base}_frac.tsv"),   sep="\t")

            # plots
            title_c = f"{target} candidates ({mode.upper()} counts)"
            title_f = f"{target} candidates ({mode.upper()} fractions)"
            png_c = os.path.join(per_dir, f"{base}_counts.png")
            png_f = os.path.join(per_dir, f"{base}_frac.png")

            save_heatmap_with_boxes(counts, title_c, png_c, mask,
                        vmin=0, xlabel="Genes (top for target)", ylabel="Predicted species",
                        force_square=True)
            save_heatmap_with_boxes(frac, title_f, png_f, mask,
                        vmin=0.0, vmax=1.0, xlabel="Genes (top for target)", ylabel="Predicted species",
                        force_square=True)


if __name__ == "__main__":
    main()




