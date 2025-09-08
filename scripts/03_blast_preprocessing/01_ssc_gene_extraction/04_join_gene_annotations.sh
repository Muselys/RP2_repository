#!/usr/bin/env bash
#BSUB -q normal
#BSUB -J join_gene_annotations_all_hits
#BSUB -n 1
#BSUB -M 8000
#BSUB -R "select[mem>8000] rusage[mem=8000]"
#BSUB -o /data/pam/team230/sm71/scratch/rp2/logs/join_gene_annotations.%J.out
#BSUB -e /data/pam/team230/sm71/scratch/rp2/logs/join_gene_annotations.%J.err
set -euo pipefail

export LC_ALL=C  # faster parsing

BP="/data/pam/team230/sm71/scratch/rp2/blast_preprocessing"

SPECIES_FILE="$BP/species_specific_core.tab"   # cols: gene_name  species
GENE_FILE="$BP/gene_data_copy.tsv"             # cols: gff_file  scaffold_name  clustering_id  gene_name
OUT_FILE="$BP/ssc_with_annotations.tab"
UNMATCHED_FILE="$BP/unmatched_genes.tab"

# Clean old outputs if present
: > "$OUT_FILE"
: > "$UNMATCHED_FILE"

awk -F'\t' -v OFS='\t' -v OUT="$OUT_FILE" -v UNM="$UNMATCHED_FILE" '
BEGIN {
  print "gene_name","species","gff_file","scaffold_name","clustering_id" > OUT
  print "gene_name","species" > UNM
}

# PASS 1: read the SMALL file (species_specific_core.tab)
# Keep species per gene and preserve row order for unmatched output.
FNR==NR {
  if (FNR==1) next  # skip header
  gene = $1
  species = $2
  if (!(gene in species_of)) {
    order[++n] = gene
    species_of[gene] = species
    matched[gene] = 0
  }
  next
}

# PASS 2: stream the BIG file (gene_data_copy.tsv)
# Print a row for EVERY matching hit (no dedup).
FNR==1 { next }  # skip header
{
  gene = $4
  if (gene in species_of) {
    print gene, species_of[gene], $1, $2, $3 >> OUT
    matched[gene]++
  }
  next
}

# END: write genes that never matched, in original order
END {
  for (i = 1; i <= n; i++) {
    g = order[i]
    if (matched[g] == 0) {
      print g, species_of[g] >> UNM
    }
  }
}
' "$SPECIES_FILE" "$GENE_FILE"
