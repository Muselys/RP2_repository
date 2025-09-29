mkdir -p /data/pam/team230/sm71/scratch/rp2/blast/chunks
# split into 128 parts (≈ ~300MB each for a 40GB file)
seqkit split -p 128 -O /data/pam/team230/sm71/scratch/rp2/blast/chunks \
  /data/pam/team230/sm71/scratch/rp2/run_blast/queries/queries.fa

BASEDIR=/data/pam/team230/sm71/scratch/rp2/blast
QUERY_DIR=$BASEDIR/chunks
OUTDIR=$BASEDIR/results_chunks
LOGS=$BASEDIR/logs_chunks
DB=/data/pam/team230/sm71/scratch/rp2/blast/ref_db
mkdir -p "$OUTDIR" "$LOGS"

N=$(ls -1 "$QUERY_DIR"/*.fa | wc -l)
bsub -q normal -J "mb[1-$N]%64" -n 4 \
  -R "span[hosts=1] select[mem>16000] rusage[mem=16000]" -M 16000 -W 60 \
  -o "$LOGS/%J_%I.out" -e "$LOGS/%J_%I.err" \
  bash -lc '
    set -eo pipefail
    module load blast/2.14.1--pl5321h6f7f691_0
    INFILE=$(ls '"$QUERY_DIR"'/*.fa | sort | sed -n "${LSB_JOBINDEX}p")
    BN=$(basename "$INFILE" .fa)
    OUT='"$OUTDIR"'/"${BN}.tsv.gz"

    blastn -task megablast \
      -query "$INFILE" -db '"$DB"' \
      -evalue 1e-10 -perc_identity 90 -qcov_hsp_perc 80 \
      -word_size 32 -culling_limit 1 -max_target_seqs 5 \
      -num_threads 4 \
      -outfmt "6 qseqid sseqid pident length qlen qcovs evalue bitscore qstart qend sstart send" \
      | gzip -c > "$OUT"
  '
#merge:
zcat $OUTDIR/*.tsv.gz > $BASEDIR/megablast_results.tsv

#run one chunk to estimate eta
/usr/bin/time -v blastn \
  -task megablast \
  -query /data/pam/team230/sm71/scratch/rp2/blast/chunks/queries.part_001.fa \
  -db /data/pam/team230/sm71/scratch/rp2/blast/ref_db \
  -perc_identity 90 \
  -qcov_hsp_perc 80 \
  -word_size 32 \
  -culling_limit 1 \
  -max_target_seqs 5 \
  -num_threads 4 \
  -outfmt "6 qseqid sseqid pident length qlen qcovs evalue bitscore qstart qend sstart send" \
  > /dev/null
