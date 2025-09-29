from pathlib import Path
import polars as pl
import sys

species_target = sys.argv[1]

base_run   = Path("/data/pam/team230/sm71/scratch/rp2/run_blast/results")
gene_fp    = Path("/data/pam/team230/sm71/scratch/rp2/panaroo_output/gene_data.parquet")
meta_fp    = Path("/data/pam/team230/sm71/scratch/rp2/metadata/File4_QC_characterisation_661K.tsv")
out_base   = base_run / "reports"

gene_join_col = "clustering_id"
gene_name_col = "gene_name"

# Load gene + metadata once
gene = pl.read_parquet(gene_fp)
meta = (
    pl.read_csv(meta_fp, separator="\t")
      .select(["sample_id","species"])
      .unique()
)

def run_species(species_target: str):
    species_slug = species_target.replace(" ", "_")
    blast_fp = base_run / f"{species_slug}_ssc_candidates.out"
    if not blast_fp.exists():
        print(f"WARNING: BLAST file not found for {species_target} -> {blast_fp}")
        return

    out_dir = out_base / species_slug
    out_dir.mkdir(parents=True, exist_ok=True)

    out_summary = out_dir / f"{species_slug}.per_qseqid.table.tsv"
    out_detail  = out_dir / f"{species_slug}.per_row.detail.tsv"

    cols = ["qseqid","sseqid","pident","length","mismatch","gapopen",
            "qstart","qend","sstart","send","evalue","bitscore","qcov"]

    df = (
        pl.read_csv(blast_fp, separator="\t", has_header=False, infer_schema_length=0)
          .rename({f"column_{i+1}": cols[i] for i in range(len(cols))})
    )

    uq_q = df.select(pl.col("qseqid").unique().alias("qseqid"))
    q2gene = (
        uq_q.join(gene, left_on="qseqid", right_on=gene_join_col, how="left")
            .group_by("qseqid")
            .agg(pl.col(gene_name_col).drop_nulls().unique().sort().alias("gene_names"))
            .with_columns(
                pl.when(pl.col("gene_names").is_null())
                  .then(pl.lit(""))
                  .otherwise(pl.col("gene_names").list.join(", "))
                  .alias("Genes it represents")
            )
            .select(["qseqid","Genes it represents"])
    )

    df_sp = (
        df.with_columns(pl.col("sseqid").str.splitn(".", 2).struct.field("field_0").alias("sample_id"))
          .join(meta, on="sample_id", how="left")
          .with_columns(pl.col("species").str.strip_chars().alias("species"))
    )

    per_q = (
        df_sp.group_by("qseqid")
             .agg([
                 pl.col("sseqid").n_unique().alias("No of sseq ids blasted against"),
                 pl.col("sseqid")
                   .filter(pl.col("species") == species_target)
                   .n_unique()
                   .alias("No of qseqids"),
             ])
             .with_columns(
                 (pl.when(pl.col("No of sseq ids blasted against") > 0)
                    .then(100.0 * pl.col("No of qseqids") / pl.col("No of sseq ids blasted against"))
                    .otherwise(0.0)
                  ).alias("% of hits")
             )
    )

    summary_tbl = (
        per_q.join(q2gene, on="qseqid", how="left")
             .select([
                 pl.col("qseqid").alias("qseqid"),
                 pl.col("Genes it represents").alias("Genes"),
                 pl.col("No of qseqids").alias("total matches"),
                 pl.col("No of sseq ids blasted against").alias("total hits"),
                 pl.col("% of hits"),
             ])
             .sort("% of hits", descending=True)
    )

    detail_tbl = (
        df.with_columns(pl.col("sseqid").str.splitn(".", 2).struct.field("field_0").alias("sample_id"))
          .join(meta, on="sample_id", how="left")
          .with_columns((pl.col("species") == species_target).alias("same_species"))
          .join(q2gene, on="qseqid", how="left")
    )

    summary_tbl.write_csv(out_summary, separator="\t", quote_style="never")
    detail_tbl.write_csv(out_detail,  separator="\t", quote_style="never")
    print(f"[OK] {species_target}: wrote {out_summary.name}, {out_detail.name}")

run_species(species_target)
