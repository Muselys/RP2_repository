"""
run this on conda env:
export POLARS_MAX_THREADS=${LSB_DJOB_NUMPROC:-8}
python make_ssc_ground_truth.py
"""
#!/usr/bin/env python3
"""
make_ssc_ground_truth.py  (PURE POLARS, robust columns)

Env (optional):
  SSC_TAB, PARQUET, OUTDIR, GENE_COL, CLUSTER_COL
"""
import os, sys
import polars as pl

SSC_TAB = os.environ.get(
    "SSC_TAB",
    "/data/pam/team230/sm71/scratch/rp2/blast_preprocessing/ssc_queries/species_specific_core.tab",
)
PARQUET = os.environ.get(
    "PARQUET",
    "/data/pam/team230/sm71/scratch/rp2/panaroo_output/gene_data.parquet",
)
OUTDIR = os.environ.get(
    "OUTDIR",
    "/data/pam/team230/sm71/scratch/rp2/blast_preprocessing/ssc_queries",
)
# user hints (optional); we’ll still auto-resolve
GENE_COL_HINT = os.environ.get("GENE_COL", "").strip()
CLUSTER_COL_HINT = os.environ.get("CLUSTER_COL", "").strip()

os.makedirs(OUTDIR, exist_ok=True)
OUTFILE = os.path.join(OUTDIR, "ssc_genes_species_clusterid.tab")

# -------- read SSC TSV lazily --------
ssc_lf = pl.scan_csv(SSC_TAB, separator="\t", infer_schema_length=0)

# sanity check columns (non-deprecated)
ssc_cols = set(ssc_lf.head(1).collect().columns)
for req in ("gene_name", "details"):
    if req not in ssc_cols:
        sys.stderr.write(f"[ERROR] {SSC_TAB} missing column: {req}\n")
        sys.exit(2)

# -------- parse species from details --------
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

# -------- read Parquet lazily --------
pq = pl.scan_parquet(PARQUET)
pq_cols = [c for c in pq.columns]

# auto-resolve gene + cluster columns
def resolve(col_hint: str, candidates: list[str], available: list[str]) -> str:
    if col_hint:
        # exact or case-insensitive match
        for c in available:
            if c == col_hint or c.lower() == col_hint.lower():
                return c
    # try candidates
    lower_av = {c.lower(): c for c in available}
    for cand in candidates:
        if cand in available:
            return cand
        if cand.lower() in lower_av:
            return lower_av[cand.lower()]
    return ""  # not found

gene_candidates = [GENE_COL_HINT] if GENE_COL_HINT else ["gene", "gene_name", "annotation_id"]
clust_candidates = [CLUSTER_COL_HINT] if CLUSTER_COL_HINT else ["cluster_id", "clustering_id"]

GENE_COL = resolve(GENE_COL_HINT, gene_candidates, pq_cols)
CLUSTER_COL = resolve(CLUSTER_COL_HINT, clust_candidates, pq_cols)

if not GENE_COL or not CLUSTER_COL:
    sys.stderr.write(
        "[ERROR] Could not resolve required columns in Parquet.\n"
        f" Available: {pq_cols}\n"
        f" Tried gene candidates: {gene_candidates}\n"
        f" Tried cluster candidates: {clust_candidates}\n"
    )
    sys.exit(3)

gdf = pq.select([pl.col(GENE_COL), pl.col(CLUSTER_COL)]).unique()

# -------- join + write --------
merged_lf = (
    ssc_long.join(gdf, left_on="gene_name", right_on=GENE_COL, how="left")
    .rename({CLUSTER_COL: "cluster_id"})
    .select(["gene_name", "species", "cluster_id"])
    .unique()
)

# stream engine (no deprecation)
df = merged_lf.collect(engine="streaming")
df.write_csv(OUTFILE, separator="\t")
print(f"[DONE] Ground truth -> {OUTFILE}")
