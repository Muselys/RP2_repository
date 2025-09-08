#!/bin/bash -l
# ------------------------------------------------------------------
# Script: concat_refs_parallel.sh
# Purpose:
#   Parallel concat of many FASTA files into one big FASTA.
#   - Supports *.fa/*.fasta/*.fna and .gz variants
#   - Runs N chunk workers in parallel on one node
#   - Writes per-chunk parts on /tmp (fast), then stitches them
# Output:
#   /data/pam/team230/sm71/scratch/rp2/blast/references/combined_references.fa
# ------------------------------------------------------------------

#BSUB -q normal
#BSUB -J concat_refs_parallel
#BSUB -n 8
#BSUB -M 32000
#BSUB -R "select[mem>32000] rusage[mem=32000] span[hosts=1]"
#BSUB -o /data/pam/team230/sm71/scratch/rp2/blast/logs/concat_refs_parallel.%J.out
#BSUB -e /data/pam/team230/sm71/scratch/rp2/blast/logs/concat_refs_parallel.%J.err

set -euo pipefail
export LC_ALL=C

# -------- config --------
REF_DIR="/data/pam/team230/sm71/scratch/rp2/blast/references/r_sequences/unzipped"
OUT_DIR="/data/pam/team230/sm71/scratch/rp2/blast/references"
LOG_DIR="/data/pam/team230/sm71/scratch/rp2/blast/logs"
COMBINED="${OUT_DIR}/combined_references.fa"
CORES="${LSB_DJOB_NUMPROC:-8}"      # how many parallel chunk workers
BATCH_SIZE=1500                     # files per chunk (tune)
# ------------------------

mkdir -p "$OUT_DIR" "$LOG_DIR"

WORK="/tmp/${USER}/concat_${LSB_JOBID}"
PARTS="${WORK}/parts"
mkdir -p "$PARTS"

echo "[INFO] WORK dir: $WORK"
echo "[INFO] Listing input FASTA files…"

LIST="${WORK}/files.lst"
find "$REF_DIR" -type f \
  \( -name "*.fa" -o -name "*.fasta" -o -name "*.fna" -o -name "*.fa.gz" -o -name "*.fasta.gz" -o -name "*.fna.gz" \) \
  -print0 > "$LIST"

NUM_FILES=$(tr -cd '\0' < "$LIST" | wc -c || true)
echo "[INFO] Found ${NUM_FILES} files"
if [[ "$NUM_FILES" -eq 0 ]]; then
  echo "ERROR: No FASTA files found under $REF_DIR" >&2
  exit 1
fi

echo "[INFO] Splitting into chunk lists (CORES=${CORES}, BATCH_SIZE=${BATCH_SIZE})…"
tr '\0' '\n' < "$LIST" > "${WORK}/files.newline.lst"
split -l "$BATCH_SIZE" "${WORK}/files.newline.lst" "${WORK}/chunk_"

echo "[INFO] Running chunk workers in parallel…"
PIDS=()
i=0
for chunk in "${WORK}"/chunk_*; do
  i=$((i+1))
  part="${PARTS}/part_${i}.fa"
  (
    : > "$part"
    while IFS= read -r fpath; do
      [[ -z "$fpath" ]] && continue
      if [[ "$fpath" == *.gz ]]; then
        gzip -cd -- "$fpath" >> "$part"
      else
        cat -- "$fpath" >> "$part"
      fi
    done < "$chunk"
  ) &
  PIDS+=($!)
  # throttle to $CORES concurrent jobs
  if [[ ${#PIDS[@]} -ge $CORES ]]; then
    wait "${PIDS[@]}"
    PIDS=()
  fi
done
# wait any stragglers
if [[ ${#PIDS[@]} -gt 0 ]]; then
  wait "${PIDS[@]}"
fi

echo "[INFO] Stitching parts into final combined FASTA…"
COMBINED_TMP="${WORK}/combined_references.fa"
: > "$COMBINED_TMP"
# cat parts in numeric order
for p in $(ls -1 "${PARTS}"/part_*.fa | sort -V); do
  cat "$p" >> "$COMBINED_TMP"
done

[[ -s "$COMBINED_TMP" ]] || { echo "ERROR: Combined FASTA is empty"; exit 1; }

echo "[INFO] Moving combined FASTA to shared storage…"
mv -f "$COMBINED_TMP" "$COMBINED"
echo "[OK] Combined FASTA: $COMBINED"
