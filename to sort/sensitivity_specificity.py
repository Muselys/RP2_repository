#STAGE 1:
#Merge HSPs on query coordinates
#figure out how much of the query is covered,
#without double-counting overlapping bases.

#Stage 2:
#Compute coverage and identity

#Stage 3:
#Call Presence/absence (Boolean)

#Stage 4:
#Confusion matrix per marker (per qseqid)

#Stage 5:
#Metrics (per qseqid)

#Step 6:
#Ranking score

#!/usr/bin/env python3
"""
Marker scoring pipeline from BLAST tabular output (outfmt 6).

Inputs
------
- BLAST TSV columns (exact order expected):
  qseqid sseqid pident length qlen qcovs evalue bitscore qstart qend sstart send
- Metadata CSV: columns required -> sseqid,is_ef,species

Outputs
-------
- <prefix>.presence.tsv : per (qseqid,sseqid): covered_len, qcov, pident_merged, present
- <prefix>.metrics.tsv  : per qseqid: TP,FN,FP,TN,Sensitivity,Specificity,PPV,NPV
- <prefix>.ranking.tsv  : per qseqid: ef_rate, max_non_ef_rate, reliability_score

Notes
-----
- "Merged identity" is computed as per-base identity over UNIQUE covered
  query bases: higher-identity HSPs claim overlapping bases first.
- Presence rule (default): pident >= 0.80 and qcov >= 0.90 (configurable).
"""

import argparse
import numpy as np
import pandas as pd
from typing import List, Tuple

# -------------------------
# Helpers
# -------------------------

def safe_div(a: float, b: float):
    return (a / b) if b else np.nan

def merge_and_score_group(group: pd.DataFrame) -> pd.Series:
    """
    Greedy merge of HSPs on query coordinates (qL,qR) for a (qseqid,sseqid) pair.
    - Sort HSPs by pident desc, then length desc.
    - Add only the non-overlapping portion of each HSP to avoid double-counting.
    Returns covered_len and matched_bases (sum of pident * contributed_len), and qlen.
    """
    qlen = int(group["qlen"].iloc[0])

    # Prepare HSP intervals with identity and length
    rows = group[["qL","qR","pident","hsplen"]].sort_values(
        ["pident","hsplen"], ascending=[False,False]
    ).itertuples(index=False, name=None)

    taken: List[Tuple[int,int]] = []  # list of merged non-overlapping intervals
    covered_len = 0
    matched_bases = 0.0

    for qL, qR, pid, hsplen in rows:
        L, R = int(qL), int(qR)
        if L > R:
            continue

        # New segments initially the full interval
        new_segments = [(L, R)]
        # Cut out overlaps with already taken intervals
        for a, b in taken:
            next_segments = []
            for x, y in new_segments:
                if y < a or x > b:
                    next_segments.append((x, y))  # no overlap
                else:
                    # left piece
                    if x < a:
                        next_segments.append((x, a - 1))
                    # right piece
                    if y > b:
                        next_segments.append((b + 1, y))
            new_segments = next_segments
            if not new_segments:
                break

        # Accumulate newly contributed bases
        for x, y in new_segments:
            add_len = y - x + 1
            if add_len > 0:
                covered_len += add_len
                matched_bases += (pid / 100.0) * add_len

        # Insert current interval, keep 'taken' merged
        taken.append((L, R))
        taken.sort()
        merged = []
        for s, e in taken:
            if not merged or s > merged[-1][1] + 1:
                merged.append([s, e])
            else:
                merged[-1][1] = max(merged[-1][1], e)
        taken = [(s, e) for s, e in merged]

    return pd.Series({
        "qlen": qlen,
        "covered_len": covered_len,
        "matched_bases": matched_bases
    })


def detect_ef_species_name(meta: pd.DataFrame, explicit_name: str = None) -> str:
    """
    Pick the E. faecalis species label.
    - If explicit_name is given, use it.
    - Else, take the most frequent 'species' among rows with is_ef==1.
    """
    if explicit_name:
        return explicit_name
    ef_species = meta.loc[meta["is_ef"] == 1, "species"]
    if ef_species.empty:
        raise ValueError("No rows with is_ef==1 found in metadata.")
    return ef_species.mode().iloc[0]


def main():
    ap = argparse.ArgumentParser(description="Score markers from BLAST outfmt 6.")
    ap.add_argument("--blast_tsv", required=True, help="BLAST results TSV (outfmt 6).")
    ap.add_argument("--meta_csv",  required=True, help="Metadata CSV with sseqid,is_ef,species.")
    ap.add_argument("--out_prefix", default="marker", help="Output filename prefix.")
    ap.add_argument("--pident_thr", type=float, default=0.80, help="Presence pident threshold (0-1).")
    ap.add_argument("--qcov_thr",   type=float, default=0.90, help="Presence qcov threshold (0-1).")
    ap.add_argument("--ef_species", default=None, help="Explicit species name for E. faecalis (optional).")
    args = ap.parse_args()

    # =========================
    # STAGE 1:
    # Merge HSPs on query coordinates
    # figure out how much of the query is covered,
    # without double-counting overlapping bases.
    # =========================

    cols = ["qseqid","sseqid","pident","length","qlen","qcovs","evalue","bitscore","qstart","qend","sstart","send"]
    dtypes = {
        "qseqid":"string",
        "sseqid":"string",
        "pident":"float32",
        "length":"int32",
        "qlen":"int32",
        "qcovs":"float32",
        "evalue":"float64",
        "bitscore":"float32",
        "qstart":"int32",
        "qend":"int32",
        "sstart":"int32",
        "send":"int32",
    }
    df = pd.read_csv(args.blast_tsv, sep="\t", header=None, names=cols, dtype=dtypes, low_memory=False)

    # Normalize query intervals
    df["qL"] = df[["qstart","qend"]].min(axis=1).astype("int32")
    df["qR"] = df[["qstart","qend"]].max(axis=1).astype("int32")
    df["hsplen"] = (df["qR"] - df["qL"] + 1).clip(lower=0).astype("int32")

    # Merge & score per (qseqid,sseqid)
    merged = (
        df.groupby(["qseqid","sseqid"], sort=False)
          .apply(merge_and_score_group)
          .reset_index()
    )

    # =========================
    # Stage 2:
    # Compute coverage and identity
    # =========================
    merged["qcov"] = merged.apply(lambda r: safe_div(r["covered_len"], r["qlen"]), axis=1)
    merged["pident_merged"] = merged.apply(
        lambda r: safe_div(r["matched_bases"], r["covered_len"]) if r["covered_len"] > 0 else 0.0, axis=1
    )

    # =========================
    # Stage 3:
    # Call Presence/absence (Boolean)
    # =========================
    merged["present"] = ((merged["pident_merged"] >= args.pident_thr) & (merged["qcov"] >= args.qcov_thr)).astype("int8")

    # Save presence table
    presence_cols = ["qseqid","sseqid","qlen","covered_len","qcov","pident_merged","present"]
    presence_path = f"{args.out_prefix}.presence.tsv"
    merged[presence_cols].to_csv(presence_path, sep="\t", index=False)

    # =========================
    # Stage 4:
    # Confusion matrix per marker (per qseqid)
    # =========================
    meta = pd.read_csv(args.meta_csv)
    required_meta = {"sseqid","is_ef","species"}
    if not required_meta.issubset(meta.columns):
        raise ValueError(f"Metadata missing required columns: {required_meta - set(meta.columns)}")

    meta["is_ef"] = meta["is_ef"].astype(int)
    pred = merged.merge(meta[["sseqid","is_ef","species"]], on="sseqid", how="left")

    def summarize_confusion(g: pd.DataFrame) -> pd.Series:
        TP = int(((g.is_ef==1) & (g.present==1)).sum())
        FN = int(((g.is_ef==1) & (g.present==0)).sum())
        FP = int(((g.is_ef==0) & (g.present==1)).sum())
        TN = int(((g.is_ef==0) & (g.present==0)).sum())
        return pd.Series({"TP":TP,"FN":FN,"FP":FP,"TN":TN})

    conf = pred.groupby("qseqid").apply(summarize_confusion).reset_index()

    # =========================
    # Stage 5:
    # Metrics (per qseqid)
    # =========================
    def add_metrics(row: pd.Series) -> pd.Series:
        TP, FN, FP, TN = row["TP"], row["FN"], row["FP"], row["TN"]
        Sens = safe_div(TP, TP + FN)
        Spec = safe_div(TN, TN + FP)
        PPV  = safe_div(TP, TP + FP)
        NPV  = safe_div(TN, TN + FN)
        return pd.Series({"Sensitivity":Sens, "Specificity":Spec, "PPV":PPV, "NPV":NPV})

    mets = conf.copy()
    mets = pd.concat([mets, conf.apply(add_metrics, axis=1)], axis=1)
    metrics_path = f"{args.out_prefix}.metrics.tsv"
    mets.to_csv(metrics_path, sep="\t", index=False)

    # =========================
    # Step 6:
    # Ranking score
    # =========================
    ef_name = detect_ef_species_name(meta, args.ef_species)

    # Detection rates per (qseqid, species)
    det = (pred
           .groupby(["qseqid","species"], as_index=False)
           .agg(n=("present","size"), hits=("present","sum")))
    det["detection_rate"] = det["hits"] / det["n"]

    ef_det = det[det["species"] == ef_name][["qseqid","detection_rate"]].rename(columns={"detection_rate":"ef_rate"})
    non_det = det[det["species"] != ef_name]

    max_non = (non_det
               .groupby("qseqid", as_index=False)["detection_rate"]
               .max()
               .rename(columns={"detection_rate":"max_non_ef_rate"}))

    rank = ef_det.merge(max_non, on="qseqid", how="left").fillna({"max_non_ef_rate": 0.0})
    rank["reliability_score"] = rank["ef_rate"] - rank["max_non_ef_rate"]
    rank = rank.sort_values(["reliability_score","ef_rate"], ascending=[False,False])

    ranking_path = f"{args.out_prefix}.ranking.tsv"
    rank.to_csv(ranking_path, sep="\t", index=False)

    print("Done.")
    print(f"Presence table : {presence_path}")
    print(f"Metrics table  : {metrics_path}")
    print(f"Ranking table  : {ranking_path}")


if __name__ == "__main__":
    main()
