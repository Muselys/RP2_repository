#!/usr/bin/env python3
"""
BSUB submission example:
bsub -q normal -n 1 -M 4000 -R "select[mem>4000] rusage[mem=4000]" \
   -o /data/pam/team230/sm71/scratch/rp2/logs/unmatched_genes.out \
   -e /data/pam/team230/sm71/scratch/rp2/logs/unmatched_genes.err \
   "module load ISG/conda && conda activate py && python /data/pam/team230/sm71/scratch/rp2/scripts/unmatched_genes.py"


"""
import polars as pl

# --- paths ---
tab_path = "/data/pam/team230/sm71/scratch/rp2/blast_preprocessing/unmatched_genes.txt "
pq_path  = "/data/pam/team230/sm71/scratch/rp2/panaroo_output/gene_data.parquet"
out_path = Path("/data/pam/team230/sm71/scratch/rp2/blast_preprocessing/matched_genes.csv")


# 1) read the input TXT (single column, no header)
targets = (
    pl.read_csv(
        txt_path.as_posix(),
        has_header=False,
        separator="\t",
        new_columns=["gene_name"],
        infer_schema_length=0,
    )
    .with_columns(pl.col("gene_name").cast(pl.Utf8).str.strip())
    .filter(pl.col("gene_name").str.len_bytes() > 0)
    .with_row_count("order")  # preserve original order
)

# 2) load parquet with the needed columns
scan = (
    pl.scan_parquet(pq_path)
    .select(["gene_name", "species", "gff_file", "scaffold_name", "clustering_id"])
    .with_columns(
        pl.col("gene_name").cast(pl.Utf8).str.strip(),
        pl.col("species").cast(pl.Utf8).str.strip(),
        pl.col("gff_file").cast(pl.Utf8).str.strip(),
        pl.col("scaffold_name").cast(pl.Utf8).str.strip(),
        pl.col("clustering_id").cast(pl.Utf8).str.strip(),
    )
)

# 3) LEFT JOIN on gene_name (only this key)
enriched = (
    targets.lazy()
    .join(scan, on="gene_name", how="left")
    .sort("order")
    .select(["gene_name", "species", "gff_file", "scaffold_name", "clustering_id"])
    .collect(streaming=True)
)

# 4) write output as CSV with header
enriched.write_csv(out_path.as_posix(), has_header=True)

print(f"Wrote matched results to: {out_path}")
print(enriched.head())