#!/bin/bash
#BSUB -J blast_run
#BSUB -o blast_run.out
#BSUB -q yesterday
#BSUB -e blast_run.err
#BSUB -n 4
#BSUB -M 32000
#BSUB -R "select[mem>32000] rusage[mem=32000]"

# === Load BLAST module ===
module load blast/2.14.1--pl5321h6f7f691_0

# === INPUT FILES AND PATHS ===
QUERY="/lustre/scratch127/pam/teams/team230/sm71/rp2/blast_run/queries/queries.fa"
DB="/lustre/scratch127/pam/teams/team230/sm71/rp2/blast_run/ref_db/combined_ref_db"
OUT="/lustre/scratch127/pam/teams/team230/sm71/rp2/blast_run/ref_db/results/blast_results.tsv"

# === Make sure output directory exists ===
mkdir -p "$(dirname "$OUT")"

echo "=== RUNNING BLAST ==="

# === Run BLAST ===
blastn \
  -task megablast \
  -query "$QUERY" \
  -db "$DB" \
  -evalue 1e-6 \
  -num_threads 8 \
  -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qcovhsp qcovs qcovus' \
  -out "$OUT"

echo "=== BLAST COMPLETE ==="
