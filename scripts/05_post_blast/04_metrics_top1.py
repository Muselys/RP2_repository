#!/usr/bin/env python3
"""
Compute per-species TP/FP/FN/TN + Sensitivity/Specificity using TOP-1 hit per query.

- PRE  = from the annotated file (no thresholds)
- POST = from the filtered file   (after thresholds)

Ranking tie-breaks (desc unless noted):
  bitscore, qcovs, pident, -evalue (smaller is better), sseqid (lex asc)

Outputs:
  - predictions_pre_post.tsv     (qseqid, gene_name, true_species, pred_pre, pred_post)
  - per_species_metrics_pre.tsv  (TP,FP,FN,TN,Sensitivity,Specificity)
  - per_species_metrics_post.tsv
  - confusion_pre.tsv            (S x S matrix of true vs predicted PRE, incl. "None")
  - confusion_post.tsv           (S x S matrix of true vs predicted POST, incl. "None")
"""

import argparse, os, sys
import numpy as np
import pandas as pd

def norm_sp(x):
    if pd.isna(x): return None
    s = str(x).strip().replace(" ", "_")
    return s if s else None

def better(score_a, score_b):
    """Return True if score_a is better than score_b (lexicographic)."""
    if score_b is None:
        return True
    return score_a > score_b

def top1_from_chunks(path, chunksize=500_000):
    """
    Stream a (possibly huge) TSV and keep the best hit per qseqid.
    Returns:
      - top:   dict qseqid -> (score_tuple, ref_species)
      - truth: dict qseqid -> (true_species, gene_name)
    """
    top = {}
    truth = {}
    usecols = [
        "qseqid","sseqid","bitscore","qcovs","pident","evalue",
        "ref_species","true_species","gene_name"
    ]
    dtype = {"qseqid":str, "sseqid":str}

    for chunk in pd.read_csv(path, sep="\t", chunksize=chunksize,
                             dtype=dtype, usecols=usecols, low_memory=False):
        # numeric coercion
        for c in ("bitscore","qcovs","pident","evalue"):
            chunk[c] = pd.to_numeric(chunk[c], errors="coerce")

        # normalize species labels
        chunk["ref_species"]  = chunk["ref_species"].map(norm_sp)
        chunk["true_species"] = chunk["true_species"].map(norm_sp)

        # iterate rows (dict keeps RAM small and is fast enough)
        for r in chunk.itertuples(index=False):
            q = r.qseqid
            if q not in truth:
                truth[q] = (r.true_species, getattr(r, "gene_name", None))

            bs   = r.bitscore if pd.notna(r.bitscore) else -np.inf
            qc   = r.qcovs    if pd.notna(r.qcovs)    else -np.inf
            pid  = r.pident   if pd.notna(r.pident)   else -np.inf
            ev   = r.evalue   if pd.notna(r.evalue)   else np.inf
            sseq = r.sseqid or ""

            score = (bs, qc, pid, -ev, sseq)  # lexicographic tuple
            cur = top.get(q)                  # (score, ref_species) or None
            if better(score, cur[0] if cur else None):
                top[q] = (score, r.ref_species)

    return top, truth

def predictions_df(truth_map, pre_top, post_top):
    rows = []
    for q, (tsp, gname) in truth_map.items():
        pp = pre_top.get(q, (None, None))[1]   # predicted species (pre)
        po = post_top.get(q, (None, None))[1]  # predicted species (post)
        rows.append((q, gname, tsp, pp, po))
    return pd.DataFrame(rows, columns=["qseqid","gene_name","true_species","pred_pre","pred_post"])

def metrics_ovr(pred, true, species_list):
    """
    One-vs-rest metrics per species.
    Accepts pandas Series or NumPy arrays for pred/true.
    """
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
        rows.append(dict(species=S, TP=TP, FP=FP, FN=FN, TN=TN,
                         Sensitivity=sens, Specificity=spec))
    return pd.DataFrame(rows)

def read_species_list(path):
    with open(path) as fh:
        S = {norm_sp(line.strip()) for line in fh if line.strip()}
    S.discard(None)
    return sorted(S)

def main():
    ap = argparse.ArgumentParser(
        description="Top-1 per-query species metrics (pre vs post filtering)."
    )
    ap.add_argument("--annotated", required=True, help="blast.annotated.tsv(.gz)")
    ap.add_argument("--filtered",   required=True, help="blast.filtered.tsv(.gz)")
    ap.add_argument("--outdir",     required=True, help="Output directory")
    ap.add_argument("--species-file", default=None, help="Optional file with one target species per line")
    ap.add_argument("--chunksize", type=int, default=500_000, help="Rows per read_csv chunk (default 500k)")
    ap.add_argument("--none-label", default="None", help="Label to use for 'no prediction' in confusion tables")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    print("[INFO] scanning annotated (pre)...", file=sys.stderr)
    pre_top, truth = top1_from_chunks(args.annotated, chunksize=args.chunksize)

    print("[INFO] scanning filtered (post)...", file=sys.stderr)
    post_top, _ = top1_from_chunks(args.filtered, chunksize=args.chunksize)

    # build predictions
    pred = predictions_df(truth, pre_top, post_top)
    pred_path = os.path.join(args.outdir, "predictions_pre_post.tsv")
    pred.to_csv(pred_path, sep="\t", index=False)

    # species scope
    if args.species_file:
        species_list = read_species_list(args.species_file)
    else:
        species_list = sorted({s for s in pred["true_species"].dropna().unique()})

    # metrics per species (OVR)
    m_pre  = metrics_ovr(pred["pred_pre"],  pred["true_species"], species_list)
    m_post = metrics_ovr(pred["pred_post"], pred["true_species"], species_list)
    m_pre.to_csv( os.path.join(args.outdir, "per_species_metrics_pre.tsv"),  sep="\t", index=False)
    m_post.to_csv(os.path.join(args.outdir, "per_species_metrics_post.tsv"), sep="\t", index=False)

    # confusion matrices (include explicit None label for clarity)
    pred_ct = pred.copy()
    pred_ct["pred_pre"]  = pred_ct["pred_pre"].fillna(args.none_label)
    pred_ct["pred_post"] = pred_ct["pred_post"].fillna(args.none_label)

    cm_pre  = pd.crosstab(pred_ct["true_species"], pred_ct["pred_pre"],  dropna=False).fillna(0).astype(int)
    cm_post = pd.crosstab(pred_ct["true_species"], pred_ct["pred_post"], dropna=False).fillna(0).astype(int)
    cm_pre.to_csv( os.path.join(args.outdir, "confusion_pre.tsv"),  sep="\t")
    cm_post.to_csv(os.path.join(args.outdir, "confusion_post.tsv"), sep="\t")

    print("[OK] wrote outputs to", args.outdir, file=sys.stderr)

if __name__ == "__main__":
    main()
