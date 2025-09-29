#!/usr/bin/env bash
set -euo pipefail

# Simple post-BLAST merge + metadata join

# ---- Paths ----
Q_META="/data/pam/team230/sm71/scratch/rp2/blast_run/post_blast/truth_queries.tsv"
R_META="/data/pam/team230/sm71/scratch/rp2/blast_run/post_blast/truth_refs.tsv"

# ---- 1) Combine all BLAST results ----
echo "[INFO] Combining shard outputs into blast.tsv"

printf 'qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tqcovhsp\tqcovs\tqcovus\n' > blast.tsv

find . -maxdepth 1 -type f -name 'q_*.megablast.top10.tsv' -size +0c -print0 \
  | sort -zV \
  | xargs -0 cat >> blast.tsv

echo "[OK] Combined BLAST results → blast.tsv"

# ---- 2) Join truth tables ----
echo "[INFO] Adding query truth labels..."
awk -F'\t' 'NR==FNR{Q[$1]=$2; next} FNR==1{print $0"\ttrue_query"} {print $0"\t"Q[$1]}' "$Q_META" blast.tsv > blast.with_q.tsv

echo "[INFO] Adding reference truth labels..."
awk -F'\t' 'NR==FNR{R[$1]=$2; next} FNR==1{print $0"\ttrue_ref"} {print $0"\t"R[$2]}' "$R_META" blast.with_q.tsv > blast.truth.tsv

rm blast.with_q.tsv
echo "[OK] Final annotated results → blast.truth.tsv"
