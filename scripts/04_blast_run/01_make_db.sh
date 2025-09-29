#REF Map ID to Sample_ID located here:
#/data/pam/team230/sm71/scratch/rp2/blast_run/ref19_norm/non_db_stuff.ref19_norm.map.tsv

#!/usr/bin/env bash
#BSUB -J makeblastdb_ref19_norm
#BSUB -q yesterday
#BSUB -n 4
#BSUB -R "span[hosts=1] select[mem>32000] rusage[mem=32000]"
#BSUB -M 32000
#BSUB -o /data/pam/team230/sm71/scratch/rp2/logs/makeblastdb_ref19_norm.%J.out
#BSUB -e /data/pam/team230/sm71/scratch/rp2/logs/makeblastdb_ref19_norm.%J.err

set -euo pipefail
module load blast/2.14.1--pl5321h6f7f691_0

# --- INPUTS ---
SRC="/data/pam/team230/sm71/scratch/rp2/blast_run/combined_ref.fa"
OUTDIR="/lustre/scratch127/pam/teams/team230/sm71/rp2/blast_run/ref19_norm"
NORM_FA="$OUTDIR/ref19_norm.fa"
MAP_TSV="$OUTDIR/ref19_norm.map.tsv"
DB_BASE="$OUTDIR/ref19_norm"

mkdir -p "$OUTDIR"

echo "[START] $(date) Normalizing headers"
awk 'BEGIN{FS=" "; OFS="\t"}
     /^>/{
        ++c;
        id=sprintf("ref19|%09d", c);
        hdr=substr($0,2);            # full original header (without ">")
        print id, hdr >> "'"$MAP_TSV"'";
        print ">" id; next
     }
     {print}
' "$SRC" > "$NORM_FA"

echo "[INFO] Normalized FASTA written: $NORM_FA"
echo "[INFO] Mapping table written: $MAP_TSV"

echo "[START] $(date) Building BLAST DB"
makeblastdb \
  -dbtype nucl \
  -input_type fasta \
  -in "$NORM_FA" \
  -title "combined_ref19_norm" \
  -out "$DB_BASE" \
  -parse_seqids \
  -max_file_sz 3GB \
  -logfile "$OUTDIR/makeblastdb.log"

echo "[DONE] $(date)"
ls -lh "$OUTDIR" | sed -n "1,200p"



