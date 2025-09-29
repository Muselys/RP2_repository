#!/usr/bin/env bash
set -euo pipefail

FA_DIR="/data/pam/team230/sm71/scratch/rp2/blast/references/queries"   # where your .fa live
LIST="/data/pam/team230/sm71/scratch/rp2/blast/references/queries.list"

# build list of queries (absolute paths), sorted and unique
find "$FA_DIR" -maxdepth 1 -type f \( -name "*.fa" -o -name "*.fasta" \) -print0 \
 | xargs -0 -I{} readlink -f "{}" \
 | sort -u > "$LIST"

N=$(wc -l < "$LIST")
[[ "$N" -gt 0 ]] || { echo "No FASTA found in $FA_DIR"; exit 1; }

echo "Submitting $N BLAST jobs…"

bsub -q normal \
  -J "blastn_one[1-$N]%8" \
  -n 2 \
  -M 4000 \
  -R "span[hosts=1] select[mem>4000] rusage[mem=4000]" \
  -o /data/pam/team230/sm71/scratch/rp2/logs/blastn.%I.out \
  -e /data/pam/team230/sm71/scratch/rp2/logs/blastn.%I.err \
  "FILE=\$(sed -n \${LSB_JOBINDEX}p $LIST); \
   echo \"[blastn_one] idx=\${LSB_JOBINDEX} -> \$FILE\"; \
   bash /data/pam/team230/sm71/scratch/rp2/blast/02.2_blastn_one.sh \"\$FILE\""
