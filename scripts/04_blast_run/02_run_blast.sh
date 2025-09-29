#!/bin/bash -l
#BSUB -J "blast_megablast[1-__NSHARDS__]"    # <-- replace __NSHARDS__ before submitting
#BSUB -o blast_megablast.%J.%I.out
#BSUB -e blast_megablast.%J.%I.err
#BSUB -n 8
#BSUB -M 16000
#BSUB -R "select[mem>16000] rusage[mem=16000] span[hosts=1]"
#BSUB -q normal

set -euo pipefail
echo "== [$(date)] host: $(hostname) job: ${LSB_JOBID:-NA} idx: ${LSB_JOBINDEX:-NA} =="

# Init modules (adjust if your module name differs)
if [ -f /etc/profile.d/modules.sh ]; then source /etc/profile.d/modules.sh; fi
module load blast/2.14.1 --silent || true

# ---- Config (use /data for inputs but resolve to /lustre canonical paths) ----
DB_DATA="/data/pam/team230/sm71/scratch/rp2/blast_run/ref19_norm/ref19_norm"
DB="$(readlink -f "$DB_DATA" || echo "$DB_DATA")"

SHARD_DIR_DATA="/data/pam/team230/sm71/scratch/rp2/blast_run/queries/shards"
SHARD_DIR="$(readlink -f "$SHARD_DIR_DATA" || echo "$SHARD_DIR_DATA")"

OUTDIR_DATA="/data/pam/team230/sm71/scratch/rp2/blast_run/results"
OUTDIR="$(readlink -f "$OUTDIR_DATA" || echo "$OUTDIR_DATA")"

# ---- Checklist helpers ----
pass(){ echo "[PASS] $1"; }
fail(){ echo "[FAIL] $1" >&2; exit 1; }
info(){ echo "[INFO] $1"; }

echo "=== CHECKLIST ==="

# 1) blastn present
BLASTN="$(command -v blastn || true)"
[ -n "$BLASTN" ] && pass "Found blastn: $BLASTN" || fail "blastn not in PATH"
"$BLASTN" -version >/dev/null || fail "blastn version check failed"

# 2) DB prefix exists (monolithic or sharded)
if ls -1 "${DB}.nin" "${DB}.nsq" "${DB}.nhr" >/dev/null 2>&1 || ls -1 ${DB}.0*.nin >/dev/null 2>&1; then
  pass "DB prefix found: $DB"
else
  fail "DB not found for prefix: $DB"
fi

# 3) Build a stable, sorted list of shard files (canonical /lustre paths)
mapfile -d '' SHARDS < <(find "$SHARD_DIR" -maxdepth 1 -type f -name 'q_*.fa' -print0 | sort -z)
NS=${#SHARDS[@]}
[ "$NS" -gt 0 ] && pass "Found $NS shard(s) under $SHARD_DIR" || fail "No q_*.fa shards found in $SHARD_DIR"

# 4) Pick the shard for this array index
IDX="${LSB_JOBINDEX:?Array index missing}"
(( IDX>=1 && IDX<=NS )) || fail "Array index $IDX out of range 1..$NS"
SHARD="${SHARDS[$((IDX-1))]}"
[ -s "$SHARD" ] && pass "Shard[$IDX]: $SHARD" || fail "Shard file missing/empty: $SHARD"

# 5) Output directory directly on /lustre (or canonical path)
mkdir -p "$OUTDIR" || true
[ -d "$OUTDIR" ] && [ -w "$OUTDIR" ] && pass "Output dir ok: $OUTDIR" || fail "Cannot write to $OUTDIR"

# 6) Optional: DB info (warn-only)
if command -v blastdbcmd >/dev/null 2>&1; then
  blastdbcmd -db "$DB" -dbtype nucl -info | head -n 5 || info "blastdbcmd info skipped"
fi

# ---- Derive output name per shard ----
BASE="$(basename "$SHARD")"
STEM="${BASE%.fa}"
OUT="$OUTDIR/${STEM}.megablast.top10.tsv"

echo "=== RUNNING BLAST ==="
time "$BLASTN" \
  -task megablast \
  -query "$SHARD" \
  -db "$DB" \
  -evalue 1e6 \
  -num_threads 8 \
  -max_target_seqs 10 \
  -max_hsps 1 \
  -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qcovhsp qcovs qcovus ' \
  -out "$OUT"

# Verify output
if [ -f "$OUT" ]; then
  LINES=$(wc -l < "$OUT" || echo 0); SIZE=$(du -h "$OUT" | cut -f1)
  [ "$LINES" -ge 0 ] && pass "Output ok: $OUT (lines=$LINES, size=$SIZE)" || fail "Output exists but empty: $OUT"
else
  fail "No output file produced: $OUT"
fi

echo "== [$(date)] shard done =="

