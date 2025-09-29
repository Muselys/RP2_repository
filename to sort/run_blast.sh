split queries into chunky shards
blast shards in parallel

Write compact output (outfmt 6) with fields you need for post-filtering + metrics.
Merge results after.
Apply thresholds after BLAST to tune specificity/sensitivity properly.

# paths to edit
DB_PREFIX="/data/pam/team230/sm71/scratch/rp2/blast_run/ref18/ref18"   # makeblastdb output prefix (no extension)
QUERY_DIR="/data/pam/team230/sm71/scratch/rp2/blast_run/queries"                 # folder with your query FASTA(s)
WORK="/data/pam/team230/sm71/scratch/rp2/blast_jobs/ref18"
LOGS="/data/pam/team230/sm71/scratch/rp2/logs"
mkdir -p "$WORK"/{shards,outs,dbcache} "$LOGS" "$QUERY_DIR"

mv /data/pam/team230/sm71/scratch/rp2/run_blast/queries/queries.fa \
    /data/pam/team230/sm71/scratch/rp2/blast_run/queries/queries.fa 


#!/usr/bin/env bash
set -euo pipefail
BLOCK="300M"  # change to 200M/500M if you want

# cat all queries (top-level *.fa) and split by record (starts with '>')
find "$QUERY_DIR" -maxdepth 1 -type f -iname '*.fa' -print0 \
 | xargs -0 cat \
 | parallel --pipe --recstart '>' --block "$BLOCK" --jobs 1 '
     i=$(printf "%04d" {#})
     printf "%s%s\n" ">" "$(cat -)" > "'"$WORK"'/shards/q_'"${i}"'.fa"
   '

# Make a list file for the array
ls -1 "$WORK"/shards/q_*.fa | LC_ALL=C sort > "$WORK/shards.list"
wc -l "$WORK/shards.list"
