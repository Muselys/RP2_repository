#!/usr/bin/env python3
"""
Query-level classification analysis (TOP-1 per query):
- predictions_pre_post.tsv
- per_species_metrics_pre.tsv, per_species_metrics_post.tsv
- confusion_pre.tsv, confusion_post.tsv

Ranking (desc unless noted):
  bitscore, qcovs, pident, -evalue (smaller is better), sseqid (lex asc)

Notes:
- Builds aligned confusion matrices: rows = the 18 target species (true),
  columns = union of all predicted labels across PRE+POST + 'None'.
"""

import argparse, os, re
import numpy as np
import pandas as pd

TARGET_SPECIES_DEFAULT = [
    "Staphylococcus aureus",
    "Staphylococcus epidermidis",
    "Staphylococcus pseudintermedius",
    "Staphylococcus haemolyticus",
    "Staphylococcus capitis",
    "Staphylococcus sciuri",
    "Staphylococcus argenteus",
    "Streptococcus pneumoniae",
    "Streptococcus pyogenes",
    "Streptococcus agalactiae",
    "Streptococcus suis",
    "Streptococcus equi",
    "Streptococcus dysgalactiae",
    "Streptococcus uberis",
    "Streptococcus mutans",
    "Streptococcus mitis",
    "Enterococcus faecium",
    "Enterococcus faecalis",
]

REQ_COLS = [
    "qseqid","gene_name","true_species",
    "sseqid","bitscore","qcovs","pident","evalue",
    "ref_species"
]

def norm_sp(x):
    if x is None or (isinstance(x, float) and pd.isna(x)): return None
    s = str(x).strip()
    if not s: return None
    s = re.sub(r"\s+"," ", s)
    return s.replace(" ", "_")

def top1_predictions(path, chunksize=500_000):
    """
    Stream TSV and keep the single best hit (global TOP-1) per qseqid.
    Returns:
      - truth_map: qseqid -> (true_species, gene_name)
      - pred_map:  qseqid -> ref_species (best hit) or None
    """
    truth_map = {}
    best = {}  # q -> (score_tuple, ref_species)

    usecols = ["qseqid","sseqid","bitscore","qcovs","pident","evalue",
               "ref_species","true_species","gene_name"]
    dtype = {"qseqid":str, "sseqid":str}

    for chunk in pd.read_csv(path, sep="\t", chunksize=chunksize, usecols=usecols,
                             dtype=dtype, low_memory=False):
        # numeric coercion
        for c in ("bitscore","qcovs","pident","evalue"):
            chunk[c] = pd.to_numeric(chunk[c], errors="coerce")
        # normalize labels
        chunk["ref_species"]  = chunk["ref_species"].map(norm_sp)
        chunk["true_species"] = chunk["true_species"].map(norm_sp)

        for r in chunk.itertuples(index=False):
            q = r.qseqid
            if q not in truth_map:
                truth_map[q] = (r.true_species, getattr(r, "gene_name", None))

            sp = r.ref_species
            if sp is None or pd.isna(r.bitscore):
                continue

            # build ranking tuple (higher is better for all except evalue)
            bs   = r.bitscore if pd.notna(r.bitscore) else -np.inf
            qc   = r.qcovs    if pd.notna(r.qcovs)    else -np.inf
            pid  = r.pident   if pd.notna(r.pident)   else -np.inf
            ev   = r.evalue   if pd.notna(r.evalue)   else np.inf
            sseq = r.sseqid or ""
            score = (bs, qc, pid, -ev, sseq)

            cur = best.get(q)
            if (cur is None) or (score > cur[0]):
                best[q] = (score, sp)

    pred_map = {q: sp for q, (_, sp) in best.items()}
    return truth_map, pred_map

def predictions_df(truth_map, pre_pred, post_pred):
    rows = []
    for q, (tsp, gname) in truth_map.items():
        rows.append((q, gname, tsp, pre_pred.get(q), post_pred.get(q)))
    return pd.DataFrame(rows, columns=["qseqid","gene_name","true_species","pred_pre","pred_post"])

def metrics_ovr(pred, true, species_list):
    ps = np.asarray(pred)
    ts = np.asarray(true)
    rows = []
    for S in species_list:
        t = (ts == S)
        p = (ps == S)
        TP = int(np.sum(t & p))
        FN = int(np.sum(t & ~p))
        FP = int(np.sum(~t & p))
        TN = int(np.sum(~t & ~p))
        sens = TP/(TP+FN) if (TP+FN)>0 else np.nan
        spec = TN/(TN+FP) if (TN+FP)>0 else np.nan
        rows.append(dict(species=S.replace("_"," "),
                         TP=TP, FP=FP, FN=FN, TN=TN,
                         Sensitivity=sens, Specificity=spec))
    return pd.DataFrame(rows)

def confusion_table_full(true, pred, row_order, col_order, none_label="None"):
    pred = pd.Series(pred, copy=True).fillna(none_label)
    true = pd.Series(true, copy=True)
    cm = pd.crosstab(true, pred, dropna=False)
    cm = cm.reindex(index=row_order, columns=col_order, fill_value=0)
    return cm

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pre", required=True)
    ap.add_argument("--post", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--chunksize", type=int, default=500_000)
    ap.add_argument("--none-label", default="None")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # Build TOP-1 predictions independently for PRE and POST
    truth_pre,  pred_pre  = top1_predictions(args.pre,  chunksize=args.chunksize)
    truth_post, pred_post = top1_predictions(args.post, chunksize=args.chunksize)

    # unify truth (prefer POST, then PRE)
    truth_map = dict(truth_pre)
    truth_map.update(truth_post)

    # predictions table
    pred_df = predictions_df(truth_map, pred_pre, pred_post)
    pred_df.to_csv(os.path.join(args.outdir, "predictions_pre_post.tsv"), sep="\t", index=False)

    # per-species metrics (OVR) over your 18 targets
    species_list = [norm_sp(s) for s in TARGET_SPECIES_DEFAULT]
    m_pre  = metrics_ovr(pred_df["pred_pre"],  pred_df["true_species"], species_list)
    m_post = metrics_ovr(pred_df["pred_post"], pred_df["true_species"], species_list)
    m_pre.to_csv( os.path.join(args.outdir, "per_species_metrics_pre.tsv"),  sep="\t", index=False)
    m_post.to_csv(os.path.join(args.outdir, "per_species_metrics_post.tsv"), sep="\t", index=False)

    # aligned confusion matrices: rows = 18 targets; cols = union of predictions + None
    rows = species_list
    preds_union = set(pred_df["pred_pre"].dropna()) | set(pred_df["pred_post"].dropna())
    cols = sorted(preds_union | {args.none_label})

    cm_pre  = confusion_table_full(pred_df["true_species"], pred_df["pred_pre"],  rows, cols, args.none_label)
    cm_post = confusion_table_full(pred_df["true_species"], pred_df["pred_post"], rows, cols, args.none_label)

    # pretty headers for human TSVs
    cm_pre_out  = cm_pre.copy();  cm_pre_out.index  = [s.replace("_"," ") for s in cm_pre.index]
    cm_pre_out.columns = [s.replace("_"," ") for s in cm_pre.columns]
    cm_post_out = cm_post.copy(); cm_post_out.index = [s.replace("_"," ") for s in cm_post.index]
    cm_post_out.columns = [s.replace("_"," ") for s in cm_post.columns]

    cm_pre_out.to_csv( os.path.join(args.outdir, "confusion_pre.tsv"),  sep="\t")
    cm_post_out.to_csv(os.path.join(args.outdir, "confusion_post.tsv"), sep="\t")

if __name__ == "__main__":
    main()


"""

python /data/pam/team230/sm71/scratch/rp2/blast_run/ref_db/results/classification_metrics.py \
  --pre  /data/pam/team230/sm71/scratch/rp2/blast_run/ref_db/results/blast_results_pre.tsv \
  --post /data/pam/team230/sm71/scratch/rp2/blast_run/ref_db/results/blast_results_post.tsv \
  --outdir /data/pam/team230/sm71/scratch/rp2/blast_run/ref_db/results/classification_analysis/ \
  --chunksize 1000000

"""