#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# build_ref_db.sh — concat FASTAs, normalize headers, and make a BLAST DB
#
# Pipeline (LSF):
#   fa_cat[1-N] -> fa_join -> fa_norm -> fa_mkdb (optional)
#
# Deterministic ordering (LC_ALL=C), resumable/idempotent, provenance files.
#
# Usage (common):
#   ./build_ref_db.sh
#   ./build_ref_db.sh --include-gz --shard-size 1000 --max-parallel 12
#   ./build_ref_db.sh --normalize-prefix ref19 --makeblastdb --dbtype nucl
#   ./build_ref_db.sh --module 'blast/2.14.1--pl5321h6f7f691_0'
#
# Key outputs (defaults):
#   ${DESTDIR}/combined_fasta.fa
#   ${DESTDIR}/norm/${NORM_PREFIX}_norm.fa
#   ${DESTDIR}/norm/${NORM_PREFIX}_norm.map.tsv
#   ${DESTDIR}/norm/${NORM_PREFIX}_norm.*   (BLAST database files)
#   ${DESTDIR}/manifests/{filelist.sha256,run.log}
# ==============================================================================

# ---------- defaults ----------
SRC="${SRC:-/data/pam/team230/sm71/scratch/rp2/atb/ref/fasta}"
DESTDIR="${DESTDIR:-/data/pam/team230/sm71/scratch/rp2/blast}"
OUT_DEFAULT="$DESTDIR/combined_fasta.fa"
OUT="${OUT:-$OUT_DEFAULT}"
TMP="${TMP:-$DESTDIR/partials}"
LOGDIR="${LOG:-/data/pam/team230/sm71/scratch/rp2/logs}"

PATTERN="fa,fasta,fna"
INCLUDE_GZ=0
SHARD_SIZE=500
MAX_PARALLEL=8
QUEUE="transfer"
MEM_MB=8000
TIME="12:00" # hh:mm
FORCE=0

# normalization + mkdb
NORM_PREFIX="${NORM_PREFIX:-ref19}"    # header prefix; set empty to skip normalize (not recommended)
NORM_DIR_DEFAULT="$DESTDIR/norm"
NORM_DIR="${NORM_DIR:-$NORM_DIR_DEFAULT}"
MAKEBLASTDB=0
DBTYPE="nucl"
MKDB_TITLE=""                          # auto if empty -> "${NORM_PREFIX}_norm"
MKDB_MAXSZ="3GB"
PARSE_SEQIDS=1
MODULE_SPEC="${MODULE_SPEC:-}"         # e.g. 'blast/2.14.1--pl5321h6f7f691_0'

export LC_ALL=C

# ---------- args ----------
usage() {
  cat <<EOF
Usage:
  $0 [--src DIR] [--dest DIR] [--out FILE]
     [--include-gz] [--shard-size N] [--max-parallel N]
     [--queue Q] [--mem-mb MB] [--time HH:MM] [--force]
     [--normalize-prefix STR] [--norm-dir DIR]
     [--makeblastdb] [--dbtype nucl|prot] [--mkdb-title STR]
     [--mkdb-max-file-sz 3GB] [--no-parse-seqids]
     [--module 'blast/2.14.1--...']

Examples:
  $0
  $0 --include-gz --shard-size 1000 --max-parallel 12
  $0 --normalize-prefix ref19 --makeblastdb --dbtype nucl
  $0 --module 'blast/2.14.1--pl5321h6f7f691_0'
EOF
}

while (( "$#" )); do
  case "$1" in
    --src) SRC="$2"; shift 2;;
    --dest) DESTDIR="$2"; shift 2;;
    --out) OUT="$2"; shift 2;;
    --include-gz) INCLUDE_GZ=1; shift;;
    --shard-size) SHARD_SIZE="$2"; shift 2;;
    --max-parallel) MAX_PARALLEL="$2"; shift 2;;
    --queue) QUEUE="$2"; shift 2;;
    --mem-mb) MEM_MB="$2"; shift 2;;
    --time) TIME="$2"; shift 2;;
    --force) FORCE=1; shift;;
    --normalize-prefix) NORM_PREFIX="$2"; shift 2;;
    --norm-dir) NORM_DIR="$2"; shift 2;;
    --makeblastdb) MAKEBLASTDB=1; shift;;
    --dbtype) DBTYPE="$2"; shift 2;;
    --mkdb-title) MKDB_TITLE="$2"; shift 2;;
    --mkdb-max-file-sz) MKDB_MAXSZ="$2"; shift 2;;
    --no-parse-seqids) PARSE_SEQIDS=0; shift;;
    --module) MODULE_SPEC="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "[FATAL] unknown arg: $1" >&2; usage; exit 1;;
  esac
done

# ---------- helpers ----------
need_cmd(){ command -v "$1" >/dev/null 2>&1 || { echo "[FATAL] missing: $1" >&2; exit 1; }; }

for c in awk sed sort find xargs sha256sum bsub; do need_cmd "$c"; done
[[ -d "$SRC" ]] || { echo "[FATAL] src not found: $SRC" >&2; exit 1; }

mkdir -p "$DESTDIR" "$TMP" "$LOGDIR" "$DESTDIR/manifests" "$NORM_DIR"

FILELIST="$DESTDIR/filelist.txt"
FILELIST_SHA="$DESTDIR/manifests/filelist.sha256"
RUNLOG="$DESTDIR/manifests/run.log"

# Normalized outputs
NORM_FA="$NORM_DIR/${NORM_PREFIX}_norm.fa"
MAP_TSV="$NORM_DIR/${NORM_PREFIX}_norm.map.tsv"
DB_BASE="$NORM_DIR/${NORM_PREFIX}_norm"
[[ -z "$MKDB_TITLE" ]] && MKDB_TITLE="${NORM_PREFIX}_norm"

# If combined already exists and not forcing, skip concat/join; still allow normalize/mkdb
SKIP_CONCAT=0
if [[ -s "$OUT" && $FORCE -eq 0 ]]; then
  echo "[SKIP] $OUT exists. (use --force to rebuild concat)"
  SKIP_CONCAT=1
fi

# ---------- deterministic file list ----------
if [[ $SKIP_CONCAT -eq 0 ]]; then
  find "$SRC" -type f \( -name '*.fa' -o -name '*.fasta' -o -name '*.fna' -o -name '*.fa.gz' -o -name '*.fasta.gz' -o -name '*.fna.gz' \) \
    | { if [[ $INCLUDE_GZ -eq 0 ]]; then grep -E '\.(fa|fasta|fna)$'; else cat; fi; } \
    | sort > "$FILELIST"

  COUNT=$(wc -l < "$FILELIST" || echo 0)
  echo "[INFO] files found: $COUNT"
  [[ $COUNT -gt 0 ]] || { echo "[FATAL] no FASTA files found under $SRC" >&2; exit 1; }

  sha256sum "$FILELIST" > "$FILELIST_SHA"

  # shard plan
  rm -f "$DESTDIR"/shard_* "$TMP"/part_*.fa "$TMP"/part_*.done || true
  split -d -a 3 -l "$SHARD_SIZE" "$FILELIST" "$DESTDIR/shard_"
  N=$(ls -1 "$DESTDIR"/shard_* | wc -l | awk '{print $1}')
  echo "[INFO] shards: $N (size: $SHARD_SIZE)"

  JOBNAME="fa_cat"
  bsub -q "$QUEUE" \
       -J "${JOBNAME}[1-$N]%${MAX_PARALLEL}" \
       -n 1 \
       -M "$MEM_MB" -R "select[mem>$MEM_MB] rusage[mem=$MEM_MB]" \
       -W "$TIME" \
       -oo "$LOGDIR/part.%I.out" \
       -eo "$LOGDIR/part.%I.err" \
       "set -euo pipefail
        DESTDIR='$DESTDIR'
        TMP='$TMP'
        i=\$(printf '%03d' \$LSB_JOBINDEX)
        shard=\"\$DESTDIR/shard_\$i\"
        out=\"\$TMP/part_\$i.fa\"
        stamp=\"\$TMP/part_\$i.done\"
        rm -f \"\$out\" \"\$stamp\" || true
        while IFS= read -r f; do
          if [[ \"\$f\" =~ \.gz\$ ]]; then
            zcat -- \"\$f\" >> \"\$out\"
          else
            cat -- \"\$f\" >> \"\$out\"
          fi
        done < \"\$shard\"
        tail -c1 \"\$out\" | od -An -t u1 | grep -q '^ *10\$' || echo >> \"\$out\"
        touch \"\$stamp\"
       "

  JOINNAME="fa_join"
  bsub -q "$QUEUE" \
       -J "$JOINNAME" \
       -w "ended($JOBNAME)" \
       -n 1 \
       -M "$MEM_MB" -R "select[mem>$MEM_MB] rusage[mem=$MEM_MB]" \
       -W "02:00" \
       -oo "$LOGDIR/join.out" \
       -eo "$LOGDIR/join.err" \
       "set -euo pipefail
        DESTDIR='$DESTDIR'
        TMP='$TMP'
        OUT='$OUT'
        for s in \$(ls -1 \"\$DESTDIR\"/shard_*); do
          i=\$(basename \"\$s\" | sed 's/^shard_//')
          [[ -f \"\$TMP/part_\$i.done\" ]] || { echo 'missing part '\$i >&2; exit 2; }
        done
        ls -1 \"\$TMP\"/part_*.fa | sort -V | xargs cat -- > \"\$OUT\"
       "
else
  JOINNAME="fa_join_existing"  # fake name so downstream -w dependencies still work
  bsub -q "$QUEUE" -J "$JOINNAME" -n 1 -W "00:01" -oo "$LOGDIR/join_exists.out" -eo "$LOGDIR/join_exists.err" "true" >/dev/null
  COUNT="$(wc -l < "${FILELIST:-/dev/null}" 2>/dev/null || echo NA)"
fi

# ---------- normalize headers (always recommended) ----------
# Load module if requested (useful for makeblastdb later)
if [[ -n "$MODULE_SPEC" ]] && command -v module >/dev/null 2>&1; then
  module load "$MODULE_SPEC" || true
fi
need_cmd awk

NORMJOB="fa_norm"
bsub -q "$QUEUE" \
     -J "$NORMJOB" \
     -w "ended($JOINNAME)" \
     -n 1 \
     -M "$MEM_MB" -R "select[mem>$MEM_MB] rusage[mem=$MEM_MB]" \
     -W "02:00" \
     -oo "$LOGDIR/norm.out" \
     -eo "$LOGDIR/norm.err" \
     "set -euo pipefail
      OUT='$OUT'
      NORM_FA='$NORM_FA'
      MAP_TSV='$MAP_TSV'
      PREFIX='$NORM_PREFIX'
      mkdir -p \"$(dirname "$NORM_FA")\"
      : > \"\$MAP_TSV\"
      awk -v pref=\"\$PREFIX\" '\"'\"'
        BEGIN{FS=" "; OFS="\\t"}
        /^>/{ ++c; id=sprintf("%s|%09d", pref, c); hdr=substr(\$0,2); print id, hdr >> ENVIRON["MAP_TSV"]; print ">" id; next }
        {print}
      '\"'\"' \"\$OUT\" > \"\$NORM_FA\"
     "

# ---------- (optional) makeblastdb ----------
if [[ $MAKEBLASTDB -eq 1 ]]; then
  need_cmd makeblastdb
  MKDBJOB="fa_mkdb"
  bsub -q "$QUEUE" \
       -J "$MKDBJOB" \
       -w "ended($NORMJOB)" \
       -n 1 \
       -M "$MEM_MB" -R "select[mem>$MEM_MB] rusage[mem=$MEM_MB]" \
       -W "02:00" \
       -oo "$LOGDIR/makeblastdb.out" \
       -eo "$LOGDIR/makeblastdb.err" \
       "set -euo pipefail
        NORM_FA='$NORM_FA'
        DB_BASE='$DB_BASE'
        TITLE='$MKDB_TITLE'
        OPT_PARSE=''
        [[ $PARSE_SEQIDS -eq 1 ]] && OPT_PARSE='-parse_seqids'
        makeblastdb -dbtype '$DBTYPE' -input_type fasta -in \"\$NORM_FA\" -title \"\$TITLE\" -out \"\$DB_BASE\" \$OPT_PARSE -max_file_sz '$MKDB_MAXSZ' -logfile \"${NORM_DIR}/makeblastdb.log\"
       "
fi

# ---------- provenance ----------
{
  echo "===== build_ref_db run $(date -u +'%Y-%m-%dT%H:%M:%SZ') ====="
  echo "SRC=$SRC"
  echo "DESTDIR=$DESTDIR"
  echo "OUT=$OUT"
  echo "COUNT=${COUNT:-NA}"
  echo "SHARD_SIZE=$SHARD_SIZE  SHARDS=${N:-NA}"
  echo "QUEUE=$QUEUE  PARALLEL=$MAX_PARALLEL  MEM_MB=$MEM_MB  TIME=$TIME"
  echo "INCLUDE_GZ=$INCLUDE_GZ  FORCE=$FORCE"
  echo "NORM_PREFIX=$NORM_PREFIX  NORM_DIR=$NORM_DIR"
  echo "MAKEBLASTDB=$MAKEBLASTDB  DBTYPE=$DBTYPE  MKDB_TITLE=$MKDB_TITLE  MKDB_MAXSZ=$MKDB_MAXSZ  PARSE_SEQIDS=$PARSE_SEQIDS"
  [[ -f "$FILELIST_SHA" ]] && echo "FILELIST_SHA256=$(cut -d' ' -f1 "$FILELIST_SHA")"
} >> "$RUNLOG"

echo "[SUBMITTED] join=$JOINNAME, normalize=$NORMJOB, mkdb=${MKDBJOB:-skipped}"
echo "[NOTE] outputs:"
echo "  combined:   $OUT"
echo "  normalized: $NORM_FA"
echo "  map tsv:    $MAP_TSV"
echo "  blast db:   $DB_BASE.* (only if --makeblastdb)"
