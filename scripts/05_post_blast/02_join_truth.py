#!/usr/bin/env python3
import pandas as pd, argparse, os, gzip, sys

COLS17 = ["qseqid","sseqid","pident","length","mismatch","gapopen",
          "qstart","qend","sstart","send","evalue","bitscore",
          "qlen","qcovhsp","qcovs","qcovus","stitle"]
COLS16 = COLS17[:-1]  # no stitle

def open_maybe_gzip(path, mode="rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)

def norm_sp(x):
    return x if pd.isna(x) else str(x).strip().replace(" ", "_")

def detect_cols(path):
    with open_maybe_gzip(path, "rt") as fh:
        line = fh.readline().strip()
    parts = line.split("\t")
    has_header = (len(parts) >= 2 and parts[0] == "qseqid")
    if has_header:
        return None, True  # file has its own header
    n = len(parts)
    if n == 17: return COLS17, False
    if n == 16: return COLS16, False
    # fallback to 16 (common case here)
    return COLS16, False

def main():
    ap = argparse.ArgumentParser(description="Join truth to BLAST (streaming).")
    ap.add_argument("--blast", required=True, help="BLAST TSV (outfmt 6); .gz ok")
    ap.add_argument("--truth-queries", required=True, help="TSV with gene_name,species,cluster_id")
    ap.add_argument("--truth-refs", required=True, help="TSV with sseqid,sample_id,species")
    ap.add_argument("--out", required=True, help="Annotated output (.tsv or .tsv.gz)")
    ap.add_argument("--chunksize", type=int, default=1_000_000)
    args = ap.parse_args()

    outdir = os.path.dirname(args.out) or "."
    os.makedirs(outdir, exist_ok=True)

    # Load truth (small)
    tq = pd.read_csv(args.truth_queries, sep="\t", dtype={"cluster_id":str})
    tq = tq.rename(columns={"cluster_id":"qseqid","species":"true_species"})
    tq["qseqid"] = tq["qseqid"].astype(str)
    tq["true_species"] = tq["true_species"].map(norm_sp)
    tq = tq.drop_duplicates("qseqid", keep="first").set_index("qseqid")

    tr = pd.read_csv(args.truth_refs, sep="\t", dtype={"sseqid":str})
    tr = tr.rename(columns={"species":"ref_species"})
    tr["sseqid"] = tr["sseqid"].astype(str)
    tr["ref_species"] = tr["ref_species"].map(norm_sp)
    tr = tr.drop_duplicates("sseqid", keep="first").set_index("sseqid")

    # Detect BLAST header/columns
    names, has_header = detect_cols(args.blast)

    # Output writer
    outfh = gzip.open(args.out, "wt") if args.out.endswith(".gz") else open(args.out, "w")
    wrote_header = False
    first_chunk = True

    # Minimal dtype hints (only those we know are present with 16 cols)
    dtypes = {"qseqid":str, "sseqid":str,
              "pident":"float32", "qcovs":"float32",
              "bitscore":"float32", "evalue":"float64"}

    reader = pd.read_csv(
        args.blast, sep="\t",
        names=names, header=0 if has_header else None,
        dtype=dtypes, chunksize=args.chunksize,
        engine="c", low_memory=False
    )

    for chunk in reader:
        # validate first chunk has what we need
        if first_chunk:
            required = {"qseqid","sseqid","pident","qcovs","bitscore","evalue"}
            missing = required - set(chunk.columns)
            if missing:
                print(f"[FAIL] Missing BLAST columns: {sorted(missing)}", file=sys.stderr)
                sys.exit(2)
            first_chunk = False

        # ensure key types
        chunk["qseqid"] = chunk["qseqid"].astype(str)
        chunk["sseqid"] = chunk["sseqid"].astype(str)

        # map truth (vectorized)
        chunk["gene_name"] = tq["gene_name"].reindex(chunk["qseqid"]).values if "gene_name" in tq.columns else None
        chunk["true_species"] = tq["true_species"].reindex(chunk["qseqid"]).values
        chunk["ref_species"]  = tr["ref_species"].reindex(chunk["sseqid"]).values
        if "sample_id" in tr.columns:
            chunk["sample_id"] = tr["sample_id"].reindex(chunk["sseqid"]).values

        # write
        chunk.to_csv(outfh, sep="\t", index=False, header=not wrote_header)
        wrote_header = True

    outfh.close()
    print(f"[OK] wrote {args.out}", file=sys.stderr)

if __name__ == "__main__":
    main()
