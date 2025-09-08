#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# Submits one BLAST job per SSC query FASTA to LSF (one file = one job).
#
# What it does:
# - Scans FA_DIR for *.fa files (your BLAST query FASTAs).
# - For each file, calls `bsub < blastn_one.sh` and passes QUERY=<file> via env,
#   so the worker script’s LSF headers and logging are respected.
# - Exits with an error if no query FASTAs are found.
#
# Prereqs:
# - Run the extractor first (e.g., `extract_ssc_query_fastas.sh`) to create *.fa queries.
# - `blastn_one.sh` must exist, be executable, include LSF headers, and read $QUERY.
#
# Usage:
#   chmod +x 02.1_submit_blast_jobs.sh
#   ./02.1_submit_blast_jobs.sh
#   (To include *.fasta too, update the find pattern.)
#
# Notes:
# - Uses `mapfile` with `-print0` (null-safe) — handles spaces in filenames.
# - Parallelism is handled by the LSF scheduler; this script just submits jobs.
# - Output locations/logs are defined inside `blastn_one.sh`.
# - Safe defaults: `set -euo pipefail` -> fail fast on errors/undefined vars.
# -----------------------------------------------------------------------------

set -euo pipefail

FA_DIR="/data/pam/team230/sm71/scratch/rp2/blast/ssc_gene_fastas"
SCRIPT="blastn_one.sh"   # worker script above (with LSF headers)

# find all .fa files (add -o -name '*.fasta' if you want those too)
mapfile -d '' files < <(find "$FA_DIR" -maxdepth 1 -type f -name '*.fa' -print0)

if (( ${#files[@]} == 0 )); then
  echo "No .fa files in $FA_DIR" >&2
  exit 1
fi

echo "Submitting ${#files[@]} BLAST jobs..."
for q in "${files[@]}"; do
  base="$(basename "${q%.*}")"
  echo "  -> $base"
  # Pass QUERY via env so the LSF headers are honored
  QUERY="$q" bsub < "$SCRIPT"
  # Alternatively (also works): bsub bash "$SCRIPT" "$q"
done

echo "All jobs submitted.