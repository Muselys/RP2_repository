# paths
BASEDIR=/data/pam/team230/sm71/scratch/rp2/blast
QUERY_DIR=$BASEDIR/chunks
RESULTS=$BASEDIR/megablast_results.tsv
DBBASE=$BASEDIR/ref_db         # has ref_db.nal etc.

mkdir -p "$BASEDIR/results" "$BASEDIR/logs"

# 64 chunks, max 24 running at once; 16 threads per job; ~8 GB/thread
bsub -q normal -J "mb[1-64]%24" -n 16 \
  -R 'span[hosts=1] select[mem>8000] rusage[mem=8000]' -M 8000 \
  -o "$BASEDIR/logs/%J_%I.out" -e "$BASEDIR/logs/%J_%I.err" \
  "bash -lc '
   set -euo pipefail
   module load blast/2.14.1--pl5321h6f7f691_0

   export BLASTDB=\"$BASEDIR\"   # use shared DB; no per-node copy
   DB=\"$DBBASE\"

   INFILE=\$(ls \"$QUERY_DIR\" | sort | sed -n \"\${LSB_JOBINDEX}p\")
   OUT=\"$BASEDIR/results/blast_\${LSB_JOBINDEX}.tsv\"

   blastn -task megablast \
     -query \"$QUERY_DIR/\$INFILE\" -db \"$DB\" \
     -evalue 1e-10 -perc_identity 90 -qcov_hsp_perc 80 \
     -word_size 32 -culling_limit 1 -max_target_seqs 5 \
     -outfmt \"6 qseqid sseqid pident length qlen qcovs evalue bitscore qstart qend sstart send\" \
     -num_threads 16 > \"\$OUT\"
  '"

bsub -w "ended(mb)" -q normal -J merge -n 1 \
  -o "$BASEDIR/logs/merge.%J.out" -e "$BASEDIR/logs/merge.%J.err" \
  "bash -lc '
   set -euo pipefail
   { printf \"qseqid\tsseqid\tpident\tlength\tqlen\tqcovs\tevalue\tbitscore\tqstart\tqend\tsstart\tsend\n\"
     cat \"$BASEDIR\"/results/*.tsv; } > \"$RESULTS\"
   echo Finished: $RESULTS
  '"

BASEDIR=/data/pam/team230/sm71/scratch/rp2/blast
QUERY_DIR=$BASEDIR/chunks
mkdir -p "$BASEDIR/results" "$BASEDIR/logs"

bsub -q normal -J mb_test -n 16 \
  -R 'span[hosts=1] select[mem>8000] rusage[mem=8000]' -M 8000 \
  -o "$BASEDIR/logs/mb_test.%J.out" -e "$BASEDIR/logs/mb_test.%J.err" \
  "bash -lc '
   set -euo pipefail
   module load blast/2.14.1--pl5321h6f7f691_0

   # make DB visible to the containerized BLAST
   export BLASTDB=$BASEDIR
   # (if your site needs strict binds, uncomment)
   # export SINGULARITY_BINDPATH=$BASEDIR:$BASEDIR

   INFILE=\$(ls \"$QUERY_DIR\" | sort | sed -n \"1p\")
   OUT=\"$BASEDIR/results/blast_1.tsv\"

   blastn -task megablast \
     -query \"$QUERY_DIR/\$INFILE\" -db ref_db \
     -evalue 1e-10 -perc_identity 90 -qcov_hsp_perc 80 \
     -word_size 32 -culling_limit 1 -max_target_seqs 5 \
     -outfmt \"6 qseqid sseqid pident length qlen qcovs evalue bitscore qstart qend sstart send\" \
     -num_threads 16 > \"\$OUT\"
  '"
