#!/usr/bin/env bash
#BSUB -J blastn_one
#BSUB -q normal
#BSUB -n 2
#BSUB -R "span[hosts=1] select[mem>4000] rusage[mem=4000]"
#BSUB -M 4000
#BSUB -o /data/pam/team230/sm71/scratch/rp2/blast/logs/blastn.%J.out
#BSUB -e /data/pam/team230/sm71/scratch/rp2/blast/logs/blastn.%J.err

set -euo pipefail

# make `module` available, then load BLAST+
if ! command -v module >/dev/null 2>&1; then
  [[ -f /etc/profile.d/modules.sh ]] && source /etc/profile.d/modules.sh
  [[ -f /usr/share/Modules/init/bash ]] && source /usr/share/Modules/init/bash
fi
module load blast/2.14.1--pl5321h6f7f691_0

# config
DB="/data/pam/team230/sm71/scratch/rp2/blast/references/ref_db/db"    # no extension
OUT_BASE="/data/pam/team230/sm71/scratch/rp2/blast/references/results"
LOG_DIR="/data/pam/team230/sm71/scratch/rp2/blast/references/logs"
mkdir -p "$OUT_BASE" "$LOG_DIR"

# input
QUERY="${1:-${QUERY:-}}"
: "${QUERY:?Must provide QUERY (.fa/.fasta path) via arg or env}"
[[ -s "$QUERY" ]] || { echo "QUERY not found or empty: $QUERY" >&2; exit 1; }

# sanity: database exists
blastdbcmd -info -db "$DB" | head -n1 >/dev/null

base="$(basename "$QUERY")"; base="${base%.fa}"; base="${base%.fasta}"
OUT="$OUT_BASE/${base}.out"
ERR="$LOG_DIR/${base}.stderr"

if [[ -s "$OUT" ]]; then
  echo "Skip $base: $OUT exists."
  exit 0
fi

# threads from LSF
NT="${LSB_DJOB_NUMPROC:-2}"

echo "Running BLAST for $base with $NT threads…"
# std (12 cols) + qcovs (13th). Matches your downstream: qcov in col 13.

blastn \
  -db "$DB" \
  -query "$QUERY" \
  -task blastn-fast \
  -perc_identity 80 \
  -qcov_hsp_perc 90 \
  -num_threads "$NT" \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" \
  -out "$OUT" \
  2> "$ERR"



echo "Done. Results: $OUT  |  Stderr: $ERR"
