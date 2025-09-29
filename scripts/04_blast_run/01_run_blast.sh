#!/usr/bin/env bash
set -euo pipefail

# 01_run_blast.sh — simple, reproducible BLAST run (megablast by default)
#
# Usage:
#   ./01_run_blast.sh \
#     --query /path/to/queries.fa \
#     --db /path/to/ref_db_prefix \
#     --out /path/to/results/blast_results.tsv
#
# Optional:
#   --threads 8            # CPU threads (default 8)
#   --task megablast       # blastn task (megablast, blastn, dc-megablast)
#   --evalue 1e-6          # e-value threshold
#   --module 'blast/2.14.1--pl5321h6f7f691_0'  # module to load (if your cluster uses modules)
#
# Notes:
#   - --db is the BLAST DB prefix (without .nin/.nsq), e.g. /.../ref19_norm/ref19_norm
#   - Output is tabular (-outfmt 6) with useful coverage columns.

# -------- defaults --------
QUERY=""
DB=""
OUT=""
THREADS="${THREADS:-8}"
TASK="megablast"
EVALUE="1e-6"
MODULE_SPEC="${MODULE_SPEC:-}"   # e.g. 'blast/2.14.1--pl5321h6f7f691_0'

usage() {
  sed -n '1,80p' "$0" | sed -n '1,/^# -------- defaults --------/p' | sed 's/^# \{0,1\}//'
  cat <<EOF
Examples:
  ./01_run_blast.sh --query /lustre/.../queries.fa --db /lustre/.../ref19_norm/ref19_norm --out /lustre/.../results/blast.tsv
  ./01_run_blast.sh --query q.fa --db db/ref19_norm --out out.tsv --threads 16 --module 'blast/2.14.1--pl5321h6f7f691_0'
EOF
}

# -------- args --------
while (( "$#" )); do
  case "$1" in
    --query) QUERY="$2"; shift 2;;
    --db) DB="$2"; shift 2;;
    --out) OUT="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --task) TASK="$2"; shift 2;;
    --evalue) EVALUE="$2"; shift 2;;
    --module) MODULE_SPEC="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "[FATAL] unknown arg: $1" >&2; usage; exit 1;;
  esac
done

# -------- tiny helpers --------
need_cmd(){ command -v "$1" >/dev/null 2>&1 || { echo "[FATAL] missing: $1" >&2; exit 1; }; }
pass(){ echo "[PASS] $*"; }
fail(){ echo "[FAIL] $*" >&2; exit 1; }
info(){ echo "[INFO] $*"; }

# -------- checks --------
[[ -n "$QUERY" && -n "$DB" && -n "$OUT" ]] || { usage; exit 1; }

# modules (optional)
if [[ -n "$MODULE_SPEC" ]] && command -v module >/dev/null 2>&1; then
  module load $MODULE_SPEC || fail "module load $MODULE_SPEC failed"
fi

need_cmd blastn
blastn -version >/dev/null || fail "blastn not working"

[[ -s "$QUERY" ]] || fail "query FASTA missing/empty: $QUERY"

# db can be monolithic or sharded; check presence
if ls -1 "${DB}.nhr" "${DB}.nin" "${DB}.nsq" >/dev/null 2>&1 || ls -1 ${DB}.0*.nin >/dev/null 2>&1; then
  pass "DB prefix ok: $DB"
else
  fail "BLAST DB not found for prefix: $DB"
fi

# out dir
mkdir -p "$(dirname "$OUT")"

# -------- run --------
info "Running blastn (task=$TASK, threads=$THREADS, evalue=$EVALUE)"
START_TS=$(date -u +'%Y-%m-%dT%H:%M:%SZ')

/usr/bin/time -f "[TIME] %E elapsed, %M KB max RSS" \
blastn \
  -task "$TASK" \
  -query "$QUERY" \
  -db "$DB" \
  -evalue "$EVALUE" \
  -num_threads "$THREADS" \
  -max_target_seqs 10 \
  -max_hsps 1 \
  -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qcovhsp qcovs qcovus' \
  -out "$OUT"

END_TS=$(date -u +'%Y-%m-%dT%H:%M:%SZ')

# -------- verify --------
if [[ -f "$OUT" ]]; then
  LINES=$(wc -l < "$OUT" || echo 0)
  SIZE=$(du -h "$OUT" | awk '{print $1}')
  pass "Output: $OUT (lines=$LINES, size=$SIZE)"
else
  fail "No output written: $OUT"
fi

echo "=== BLAST COMPLETE ==="
echo "Started : $START_TS"
echo "Finished: $END_TS"
