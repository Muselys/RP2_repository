#Check Folder Storage Size:
# du -sh /data/pam/team230/sm71/scratch/rp2/atb/ref/fasta
# 933G

#Paths
SRC="/data/pam/team230/sm71/scratch/rp2/atb/ref/fasta"
DESTDIR="/data/pam/team230/sm71/scratch/rp2/blast"
OUT="$DESTDIR/combined_fasta.fa"
TMP="$DESTDIR/partials"
LOG="/data/pam/team230/sm71/scratch/rp2/logs"
mkdir -p "$TMP" "$LOG"

#Build a deterministic file list
find "$SRC" -type f -name '*.fa' | LC_ALL=C sort > "$DESTDIR/filelist.txt"
wc -l "$DESTDIR/filelist.txt"   # sanity: how many files?

#Split into shards (500 files per shard)
split -d -a 3 -l 500 "$DESTDIR/filelist.txt" "$DESTDIR/shard_"
# shards look like: shard_000, shard_001, ...
N=$(ls -1 "$DESTDIR"/shard_* | wc -l)
echo "Shards: $N"

#Submit array job on the transfer queue
bsub -q transfer -J "fa_cat[1-$N]%8" -n 1 -M 8000 -R "select[mem>8000] rusage[mem=8000]" -W 12:00 \
  -oo "$LOG/part.%I.out" -eo "$LOG/part.%I.err" '
    set -euo pipefail
    DESTDIR="'"$DESTDIR"'"
    TMP="'"$TMP"'"
    i=$(printf "%03d" $LSB_JOBINDEX)
    shard="$DESTDIR/shard_$i"
    out="$TMP/part_$i.fa"
    # Simple concat (fast, tiny RAM)
    xargs -a "$shard" -r cat -- > "$out"
  '

#join the shards only when the array ends
bsub -q transfer -J "fa_join" -w "ended(fa_cat)" -n 1 -M 8000 -R "select[mem>8000] rusage[mem=8000]" -W 02:00 \
  -oo "$LOG/join.out" -eo "$LOG/join.err" \
  'ls -1 "'"$TMP"'/part_"*.fa | sort -V | xargs cat -- > "'"$OUT"'"'

#Count headers in parts vs final (quick sanity)
grep -h -c "^>" "$TMP"/part_*.fa | awk '{s+=$1} END{print "headers_in_parts:", s}'
grep -c "^>" "$OUT"
