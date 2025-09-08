#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# Runs a single BLASTn job for one SSC query FASTA (one file = one job) on LSF.
#
# What it does:
# - Loads BLAST+ module and checks the nucleotide DB with `blastdbcmd -info`.
# - Runs `blastn` against DB for the provided QUERY (*.fa / *.fasta).
# - Writes tabular results (outfmt 6) to results/<query>.out and stderr to logs/.
# - Skips work if the expected .out already exists (idempotent).
#
# Inputs:
#   - QUERY: path to a query FASTA (passed as $1 or env VAR QUERY)
#   - DB: BLAST database prefix (nucl), built earlier with makeblastdb
#
# Outputs:
#   - references/results/<basename>.out    (TSV: qseqid sseqid pident length mismatch gapopen qstart qend sstart send qcovs evalue bitscore)
#   - references/logs/<basename>.err       (stderr from blastn)
#
# LSF resources:
#   - Queue: normal | Cores: 1 | Mem: 4 GB | Wall: 2h
#
# Usage:
#   # submitted by the batch submitter:
#   QUERY=/path/to/query.fa bsub < 02.2_blastn_one.sh
#   # or directly:
#   bsub "bash -lc 'QUERY=/path/to/query.fa ./02.2_blastn_one.sh'"
#
# Notes:
#   - Requires module system; falls back to sourcing common module init scripts.
#   - Outfmt includes qcovs and bitscore; identity threshold set with -perc_identity 80.
#   - Final echo lines show how many .out files exist and how many jobs named blastn_one are running.
#   - Companion script: 02.1_submit_blast_jobs.sh submits one of these per query FASTA.
# -----------------------------------------------------------------------------

#BSUB -J blastn_one
#BSUB -q normal
#BSUB -n 1
#BSUB -R "select[mem>4000] rusage[mem=4000]"
#BSUB -M 4000
#BSUB -W 02:00
#BSUB -o /data/pam/team230/sm71/scratch/rp2/blast/logs/blastn.%J.out
#BSUB -e /data/pam/team230/sm71/scratch/rp2/blast/logs/blastn.%J.err

# make 'module' available
if ! command -v module >/dev/null 2>&1; then
  [[ -f /etc/profile.d/modules.sh ]] && source /etc/profile.d/modules.sh
  [[ -f /usr/share/Modules/init/bash ]] && source /usr/share/Modules/init/bash
fi
module load blast/2.14.1--pl5321h6f7f691_0

set -euo pipefail

DB="/data/pam/team230/sm71/scratch/rp2/blast/references/ref_db/db"   # DB prefix (no extension)
OUT_BASE="/data/pam/team230/sm71/scratch/rp2/blast/references/results"
LOG_DIR="/data/pam/team230/sm71/scratch/rp2/blast/references/logs"
mkdir -p "$OUT_BASE" "$LOG_DIR"

#Take the query from CLI arg or env
QUERY="${1:-${QUERY:-}}"
: "${QUERY:?Must provide QUERY (.fa filepath) via arg or env}"
[[ -s "$QUERY" ]] || { echo "QUERY not found or empty: $QUERY" >&2; exit 1; }

# sanity check DB
blastdbcmd -info -db "$DB" | head -n 5 >/dev/null

# derive names
base="$(basename "$QUERY")"; base="${base%.fa}"; base="${base%.fasta}"
OUT="$OUT_BASE/${base}.out"
ERR="$LOG_DIR/${base}.err"

# skip if already done (idempotent)
if [[ -s "$OUT" ]]; then
  echo "Skip $base: $OUT exists."
  exit 0
fi

echo "Running BLAST for $base ..."
blastn \
  -db "$DB" \
  -query "$QUERY" \
  -perc_identity 80 \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send qcovs evalue bitscore" \
  -out - 2> >(tee "$ERR" >&2) | tee "$OUT"

echo "Done. Results: $OUT  |  Stderr: $ERR"



#count how many .out blast results are complete:
ls -1 /data/pam/team230/sm71/scratch/rp2/blast/references/results/*.out 2>/dev/null | wc -l

#View all blastn_one jobs running right not under my name:
bjobs -a -u sm71  | grep blastn_one | wc -l