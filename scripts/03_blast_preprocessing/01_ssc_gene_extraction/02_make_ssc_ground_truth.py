#!/usr/bin/env python3
"""
01_make_ssc_ground_truth.py  — build SSC ground-truth (gene_name, species, cluster_id)
Reproducible, CLI-driven, pure Polars (lazy), CSV inputs/outputs.

Inputs (required)
  --ssc-tab PATH        species_specific_core.tab (from make_ssc_full_tab.sh)
  --gene-data PATH      gene_data.csv  (Panaroo gene data; CSV)

Options
  --outdir PATH         output directory (default: ./blast_preprocessing/ssc_queries)
  --gene-col NAME       override gene column name in gene_data.csv
  --cluster-col NAME    override cluster column name in gene_data.csv
  --threads N           POLARS_MAX_THREADS (default: 8)
  -h/--help             usage

Output
  <outdir>/ssc_genes_species_clusterid.tab   (TSV; gene_name, species, cluster_id)
  <outdir>/ssc_ground_truth_manifest.txt     (provenance: md5s, columns, versions)
"""
from __future__ import annotations
import argparse, os, sys, hashlib
import polars as pl

def md5(path: str) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

def parse_args():
    p = argparse.ArgumentParser(add_help=True)
    p.add_argument("--ssc-tab", required=True)
    p.add_argument("--gene-data", required=True)          # CSV (not parquet)
    p.add_argument("--outdir", default="blast_preprocessing/ssc_queries")
    p.add_argument("--gene-col", default="", help="column in gene_data.csv for gene name")
    p.add_argument("--cluster-col", default="", help="column in gene_data.csv for cluster id")
    p.add_argument("--threads", type=int, default=8)
    return p.parse_args()

def resolve(col_hint: str, candidates: list[str], available: list[str]) -> str:
    if col_hint:
        for c in available:
            if c == col_hint or c.lower() == col_hint.lower():
                return c
    lower_av = {c.lower(): c for c in available}
    for cand in candidates:
        if cand in available: return cand
        if cand.lower() in lower_av: return lower_av[cand.lower()]
    return ""

def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    os.environ["POLARS_MAX_THREADS"] = str(args.threads)

    ssc_tab = args.ssc_tab
    gene_csv = args.gene_data
    out_tsv = os.path.join(args.outdir, "ssc_genes_species_clusterid.tab")
    manifest = os.path.join(args.outdir, "ssc_ground_truth_manifest.txt")

    # ---- validate inputs
    for p in (ssc_tab, gene_csv):
        if not os.path.isfile(p):
            sys.stderr.write(f"[ERROR] missing file: {p}\n")
            sys.exit(2)

    # ---- read SSC table (TSV, lazy)
    ssc_lf = pl.scan_csv(ssc_tab, separator="\t", infer_schema_length=0)
    first = ssc_lf.head(1).collect()
    ssc_cols = set(first.columns)
    for req in ("gene_name", "details"):
        if req not in ssc_cols:
            sys.stderr.write(f"[ERROR] {ssc_tab} missing column: {req}\n")
            sys.exit(3)

    # parse species from 'details' (Core: ...), drop Inter:/Rare:
    core = (
        pl.col("details")
        .str.replace_all(r"(?s).*Core:\s*", "")
        .str.replace_all(r"(?s)\s+Inter:.*$", "")
        .str.replace_all(r"(?s)\s+Rare:.*$", "")
        .str.strip_chars()
    )

    ssc_long = (
        ssc_lf
        .with_columns(core.alias("__core"))
        .with_columns(pl.when(pl.col("__core") == "").then(None).otherwise(pl.col("__core")).alias("__core"))
        .with_columns(pl.col("__core").str.split("+").alias("__species_list"))
        .with_columns(pl.col("__species_list").list.eval(pl.element().str.strip_chars()).alias("__species_list"))
        .explode("__species_list")
        .filter(pl.col("__species_list").is_not_null() & (pl.col("__species_list") != ""))
        .rename({"__species_list": "species"})
        .select(["gene_name", "species"])
    )

    # ---- read gene_data.csv (lazy)
    gene_lf = pl.scan_csv(gene_csv)
    gene_cols = gene_lf.columns

    gene_candidates = [c for c in [args.gene_col] if c] or ["gene", "gene_name", "annotation_id"]
    cluster_candidates = [c for c in [args.cluster_col] if c] or ["cluster_id", "clustering_id"]

    gene_col = resolve(args.gene_col, gene_candidates, gene_cols)
    cluster_col = resolve(args.cluster_col, cluster_candidates, gene_cols)
    if not gene_col or not cluster_col:
        sys.stderr.write(
            "[ERROR] Could not resolve required columns in gene_data.csv.\n"
            f" Available: {gene_cols}\n"
            f" Tried gene candidates: {gene_candidates}\n"
            f" Tried cluster candidates: {cluster_candidates}\n"
        )
        sys.exit(4)

    gdf = gene_lf.select([pl.col(gene_col), pl.col(cluster_col)]).unique()

    # ---- join & write
    merged_lf = (
        ssc_long.join(gdf, left_on="gene_name", right_on=gene_col, how="left")
        .rename({cluster_col: "cluster_id"})
        .select(["gene_name", "species", "cluster_id"])
        .unique()
    )

    df = merged_lf.collect(engine="streaming")
    df.write_csv(out_tsv, separator="\t")

    # ---- manifest (provenance)
    with open(manifest, "w") as fh:
        fh.write("# 01_make_ssc_ground_truth.py\n")
        fh.write(f"ssc_tab\t{ssc_tab}\n")
        fh.write(f"gene_data_csv\t{gene_csv}\n")
        fh.write(f"out_tsv\t{out_tsv}\n")
        fh.write(f"md5_ssc_tab\t{md5(ssc_tab)}\n")
        fh.write(f"md5_gene_data_csv\t{md5(gene_csv)}\n")
        fh.write(f"resolved_gene_col\t{gene_col}\n")
        fh.write(f"resolved_cluster_col\t{cluster_col}\n")
        try:
            fh.write(f"polars_version\t{pl.__version__}\n")
        except Exception:
            pass

    print(f"[DONE] Ground truth -> {out_tsv}")

if __name__ == "__main__":
    main()
