#!/usr/bin/env python3
# search_missing_groups.py (robust columns, pure Polars)

import os, sys
import polars as pl

INDIR   = "/data/pam/team230/sm71/scratch/rp2/blast_preprocessing/ssc_queries"
MISSING = f"{INDIR}/missing_group_genes.txt"
SEARCH  = "/data/pam/team230/sm71/scratch/rp2/panaroo_output/gene_data.parquet"
OUT_FOUND   = f"{INDIR}/found_matches.tsv"
OUT_MISSING = f"{INDIR}/still_missing.txt"

# optional hints via env
GENE_COL_HINT    = os.environ.get("GENE_COL", "").strip()
CLUSTER_COL_HINT = os.environ.get("CLUSTER_COL", "").strip()

# --- load queries ---
if not os.path.exists(MISSING):
    sys.stderr.write(f"[ERROR] missing file not found: {MISSING}\n")
    sys.exit(1)

queries = (
    pl.read_csv(MISSING, has_header=False, new_columns=["gene_query"])
      .with_columns(pl.col("gene_query").str.strip_chars())
      .filter(pl.col("gene_query") != "")
      .unique()
)
qset = set(queries.get_column("gene_query").to_list())
print(f"[INFO] loaded {len(qset)} query genes from {MISSING}")

# --- inspect parquet schema safely ---
pq = pl.scan_parquet(SEARCH)
available = pq.collect_schema().names()
print(f"[INFO] parquet columns: {available}")

def resolve_col(hint: str, candidates: list[str], available: list[str]) -> str:
    # exact or case-insensitive match on hint
    if hint:
        for c in available:
            if c == hint or c.lower() == hint.lower():
                return c
    # try candidates, case-insensitive
    lower_av = {c.lower(): c for c in available}
    for cand in candidates:
        if cand in available:
            return cand
        lc = cand.lower()
        if lc in lower_av:
            return lower_av[lc]
    return ""

gene_candidates   = [c for c in [GENE_COL_HINT, "gene name", "gene", "annotation_id"] if c]
cluster_candidates= [c for c in [CLUSTER_COL_HINT, "clustering_id", "cluster_id"] if c]

GENE_COL    = resolve_col(GENE_COL_HINT, gene_candidates, available)
CLUSTER_COL = resolve_col(CLUSTER_COL_HINT, cluster_candidates, available)

if not GENE_COL or not CLUSTER_COL:
    sys.stderr.write(
        "[ERROR] Could not resolve required Parquet columns.\n"
        f" Available: {available}\n"
        f" Tried gene candidates: {gene_candidates}\n"
        f" Tried cluster candidates: {cluster_candidates}\n"
    )
    sys.exit(2)

print(f"[INFO] using GENE_COL={GENE_COL}, CLUSTER_COL={CLUSTER_COL}")

# --- pass 1: exact join on the resolved gene column ---
pq_min = pq.select([pl.col(GENE_COL), pl.col(CLUSTER_COL)]).unique()

found1 = (
    queries.lazy()
    .join(pq_min, left_on="gene_query", right_on=GENE_COL, how="inner")
    .select([
        pl.col("gene_query"),
        pl.lit(GENE_COL).alias("matched_column"),
        pl.col(GENE_COL).alias("matched_value"),
        pl.col(CLUSTER_COL)
    ])
)

found1_df = found1.collect(engine="streaming")
found1_set = set(found1_df["gene_query"].to_list())
remaining = [g for g in qset if g not in found1_set]
print(f"[INFO] pass1 exact join matched {len(found1_set)}; remaining {len(remaining)}")

# --- pass 2: exact match across multiple columns (skip big sequences) ---
SEARCH_COLS = [c for c in ["gene name","annotation_id","description","scaffold_name","gff_file"] if c in available]

if remaining and SEARCH_COLS:
    rem_df = pl.DataFrame({"gene_query": remaining}).lazy()

    long = (
        pq.select(SEARCH_COLS + [CLUSTER_COL])
          .melt(id_vars=[CLUSTER_COL], variable_name="matched_column", value_name="matched_value")
    )

    found2 = (
        rem_df.join(long, left_on="gene_query", right_on="matched_value", how="inner")
              .select(["gene_query", "matched_column", "matched_value", CLUSTER_COL])
              .unique()
    )
    found2_df = found2.collect(engine="streaming")
else:
    found2_df = pl.DataFrame({"gene_query": [], "matched_column": [], "matched_value": [], CLUSTER_COL: []})

# --- combine + write ---
found_all = pl.concat([found1_df, found2_df], how="vertical").unique()
found_all.write_csv(OUT_FOUND, separator="\t")
print(f"[INFO] wrote matches -> {OUT_FOUND} ({found_all.shape[0]} rows)")

matched_set = set(found_all["gene_query"].to_list())
still = [g for g in qset if g not in matched_set]
with open(OUT_MISSING, "w") as fh:
    for g in sorted(still):
        fh.write(f"{g}\n")
print(f"[INFO] wrote still-missing -> {OUT_MISSING} ({len(still)} genes)")

print("[DONE]")
