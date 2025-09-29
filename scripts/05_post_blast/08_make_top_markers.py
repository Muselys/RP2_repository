#!/usr/bin/env python3
import argparse, os
import numpy as np
import pandas as pd

def norm_sp(x):
    if pd.isna(x): return None
    s = str(x).strip().replace(" ", "_")
    return s or None

def read_species_list(path):
    with open(path) as fh:
        vals = [norm_sp(l.strip()) for l in fh if l.strip()]
    return [v for v in vals if v]

def load_predictions(path):
    usecols = ["qseqid","gene_name","true_species","pred_pre","pred_post"]
    df = pd.read_csv(path, sep="\t", usecols=usecols, low_memory=False)
    for c in ["true_species","pred_pre","pred_post"]:
        df[c] = df[c].map(norm_sp)
    # keep only rows with a gene_name and a POST prediction (species label present)
    df = df.dropna(subset=["gene_name", "pred_post"]).copy()
    return df

def species_by_genes_matrix(pred, species_list, genes, col_pred="pred_post"):
    if not genes:
        return pd.DataFrame(index=species_list, columns=[], dtype=int)
    d = pred[pred["gene_name"].isin(genes)].copy()
    grp = d.groupby([col_pred, "gene_name"])["qseqid"].size().reset_index(name="n")
    M = grp.pivot(index=col_pred, columns="gene_name", values="n").fillna(0).astype(int)
    M = M.reindex(index=species_list, columns=genes, fill_value=0)
    return M

def main():
    ap = argparse.ArgumentParser(
        description="Build a flat table of top-K candidate markers per target species (POST only)."
    )
    ap.add_argument("--predictions", required=True, help="predictions_pre_post.tsv")
    ap.add_argument("--species-file", required=True, help="One species per line (18 lines)")
    ap.add_argument("--outdir", required=True, help="Where to write topK_markers_post.tsv")
    ap.add_argument("--topk", type=int, default=20, help="Top K markers per species (default 20)")
    ap.add_argument("--min-gene-count", type=int, default=1,
                    help="Drop genes with fewer than this many total hits across species (default 1)")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    pred = load_predictions(args.predictions)
    species = read_species_list(args.species_file)

    # precompute total hits per gene across all species (POST)
    totals = (pred.groupby("gene_name")["qseqid"]
                .size()
                .rename("total_gene_hits")
                .astype(int))

    # optional filter by min total hits
    keep_genes = set(totals[totals >= args.min_gene_count].index.tolist())
    pred = pred[pred["gene_name"].isin(keep_genes)].copy()
    totals = totals.loc[sorted(keep_genes)] if keep_genes else totals.iloc[0:0]

    rows = []
    for target in species:
        # rank genes by how often they’re predicted AS this target (POST)
        counts_target = (pred[pred["pred_post"] == target]
                         .groupby("gene_name")["qseqid"].size()
                         .sort_values(ascending=False))
        if counts_target.empty:
            continue

        ranked_genes = counts_target.index.tolist()[:args.topk]

        # build counts across species for these ranked genes to assess uniqueness + fractions
        counts_mat = species_by_genes_matrix(pred, species, ranked_genes, col_pred="pred_post")
        total_gene_hits = counts_mat.sum(axis=0)

        # for each ranked gene, assemble output row
        for rank, g in enumerate(ranked_genes, start=1):
            count_in_target = int(counts_mat.loc[target, g])
            total_hits = int(total_gene_hits[g])
            frac_in_target = (count_in_target / total_hits) if total_hits > 0 else 0.0
            unique_to_target = (counts_mat[g] > 0).sum() == 1 and count_in_target > 0

            rows.append({
                "species": target,
                "rank": rank,
                "gene_name": g,
                "count_in_target": count_in_target,
                "total_gene_hits": total_hits,
                "frac_in_target": round(frac_in_target, 6),
                "unique_to_target": unique_to_target
            })

    out_path = os.path.join(args.outdir, f"top{args.topk}_markers_post.tsv")
    df_out = pd.DataFrame(rows).sort_values(["species", "rank"])
    df_out.to_csv(out_path, sep="\t", index=False)
    print(f"[OK] wrote {out_path} ({len(df_out)} rows)")

if __name__ == "__main__":
    main()

"""
EVAL=/lustre/scratch127/pam/teams/team230/sm71/rp2/blast_run/post_blast/eval
python make_top_markers_post.py \
  --predictions $EVAL/predictions_pre_post.tsv \
  --species-file /lustre/scratch127/pam/teams/team230/sm71/rp2/blast_run/post_blast/targets_18.txt \
  --outdir $EVAL/top20_candidates \
  --topk 20 \
  --min-gene-count 1
"""

