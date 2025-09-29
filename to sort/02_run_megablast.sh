BASEDIR=/data/pam/team230/sm71/scratch/rp2/blast
QUERY_DIR="$BASEDIR/chunks"
ALIAS_DIR="$BASEDIR"                 # where ref_db.* lives
OUTDIR="$BASEDIR/results"
LOGS="$BASEDIR/logs"
mkdir -p "$OUTDIR" "$LOGS"

N=$(ls -1 "$QUERY_DIR"/*.fa 2>/dev/null | wc -l)
[ "$N" -gt 0 ] || { echo "No .fa chunks in $QUERY_DIR" >&2; exit 1; }

bsub -q normal -J "mb[1-$N]%64" -n 4 \
  -R 'span[hosts=1] select[mem>16000] rusage[mem=16000]' -M 16000 -W 60 \
  -o "$LOGS/%J_%I.out" -e "$LOGS/%J_%I.err" \
  bash -lc '
    set -euo pipefail
    module load blast/2.14.1--pl5321h6f7f691_0

    # -------- node-local scratch --------
    TMPDIR=${TMPDIR:-/scratch/$USER/$LSB_JOBID}
    [ -d "$TMPDIR" ] || TMPDIR=$(pwd)/tmp_$LSB_JOBID
    mkdir -p "$TMPDIR/db"

    # -------- bind host paths into the container --------
    export SINGULARITY_BINDPATH="$BASEDIR:$BASEDIR,$TMPDIR:$TMPDIR"

    # -------- stage DB locally (faster than NFS) --------
    if command -v rsync >/dev/null 2>&1; then
      rsync -a "$ALIAS_DIR"/ref_db.* "$TMPDIR/db/"
    else
      cp -p "$ALIAS_DIR"/ref_db.* "$TMPDIR/db/"
    fi
    export BLASTDB="$TMPDIR/db"
    DB="$TMPDIR/db/ref_db"

    # -------- pick my chunk --------
    INFILE=$(ls "$QUERY_DIR"/*.fa | sort | sed -n "${LSB_JOBINDEX}p")
    [ -n "$INFILE" ] || { echo "No input for index $LSB_JOBINDEX" >&2; exit 2; }
    BN=$(basename "$INFILE" .fa)
    OUT="$OUTDIR/${BN}.tsv.gz"

    # -------- blastn (balanced profile) --------
    blastn -task megablast \
      -query "$INFILE" -db "$DB" \
      -evalue 1e-10 -perc_identity 95 -qcov_hsp_perc 80 \
      -word_size 32 -culling_limit 1 -max_hsps 1 -max_target_seqs 5 \
      -num_threads 4 \
      -outfmt "6 qseqid sseqid pident length qlen qcovs evalue bitscore qstart qend sstart send" \
      | gzip -c > "$OUT"
  '
