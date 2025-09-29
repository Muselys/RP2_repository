
./make_ssc_full_tab.sh \
  --class twilight_analysis/classification.tab \
  --outdir blast_preprocessing/ssc_queries

#############################################################################################
# optional threads control
export POLARS_MAX_THREADS=8

python scripts/03_blast_preprocessing/02_make_ssc_ground_truth.py \
  --ssc-tab blast_preprocessing/ssc_queries/species_specific_core.tab \
  --gene-data panaroo_output/gene_data.csv \
  --outdir blast_preprocessing/ssc_queries
# (if your gene_data.csv uses different column names)
#   --gene-col ORF --cluster-col gene_cluster


#############################################################################################

scripts/03_blast_preprocessing/02_make_queries_pool_all.sh \
  --ssc blast_preprocessing/ssc_with_annotations1.tab \
  --fasta panaroo_output/combined_DNA_CDS.fasta \
  --outdir run_blast/queries
# Submit to LSF:
#   bsub < scripts/03_blast_preprocessing/03_make_queries_pool_all.sh
