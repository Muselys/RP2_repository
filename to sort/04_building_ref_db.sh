#!/usr/bin/env bash
#BSUB -q normal
#BSUB -J ref_list
#BSUB -n 1
#BSUB -M 8000
#BSUB -R "select[mem>8000] rusage[mem=8000] span[hosts=1]"
#BSUB -o /data/pam/team230/sm71/scratch/rp2/logs/ref_list.%J.out
#BSUB -e /data/pam/team230/sm71/scratch/rp2/logs/ref_list.%J.err
set -euo pipefail
export LC_ALL=C

# Paths
RUNBLAST_BASE="/data/pam/team230/sm71/scratch/rp2/run_blast"
RTAB="/data/pam/team230/sm71/scratch/rp2/panaroo_output/gene_presence_absence.Rtab"

REFDB="${RUNBLAST_BASE}/ref_db"
OUT_URLS="${REFDB}/urls.txt"              # URL-only list
OUT_URLS_TSV="${REFDB}/urls.tsv"         # sample_id<TAB>url
OUT_MISSING="${REFDB}/missing_ids.txt"   # samples with no URL found
OUT_IDS="${REFDB}/accessions.txt"        # one sample/accession per line
DL_DIR="${REFDB}/r_sequences"            # downloads land here

# Make directories
mkdir -p "$REFDB" "$DL_DIR"

echo "[1/4] Extracting sample IDs from Rtab header → $OUT_IDS"
head -n1 "$RTAB" \
  | cut -f2- \
  | tr '\t' '\n' \
  | sed 's/^[[:space:]]*//; s/[[:space:]]*$//' \
  | awk 'NF' \
  | LC_ALL=C sort -u > "$OUT_IDS"
echo "   IDs: $(wc -l < "$OUT_IDS")"

# Clean outputs
: > "$OUT_URLS"
: > "$OUT_URLS_TSV"
: > "$OUT_MISSING"

API="https://www.ebi.ac.uk/ena/portal/api/search"

fetch_url_for_sample () {
  # Prints first matching HTTPS URL for the given sample ID (if any), else prints nothing.
  local S="$1"
  curl -fsS -G "$API" \
    --data-urlencode 'result=analysis' \
    --data-urlencode 'format=tsv' \
    --data-urlencode 'limit=0' \
    --data-urlencode 'fields=analysis_accession,sample_accession,analysis_type,submitted_ftp' \
    --data-urlencode "query=(sample_accession=\"$S\" OR secondary_sample_accession=\"$S\") AND analysis_type=\"SEQUENCE_ASSEMBLY\"" \
  | awk -F'\t' '
      NR>1 && $4!=""{
        split($4,a,";");
        for(i in a){
          if(a[i] ~ /(contigs|scaffolds).*\.(fa|fna|fasta)(\.gz)?$/){
            u=a[i]; gsub(/^ftp:\/\//,"https://",u);
            print u; exit
          }
        }
      }'
}

echo "[2/4] Querying ENA (analysis endpoint) per sample…"
n=0; found=0; miss=0
while IFS= read -r S; do
  n=$((n+1))
  url="$(fetch_url_for_sample "$S" || true)"
  if [[ -n "${url:-}" ]]; then
    echo "$url" >> "$OUT_URLS"
    printf "%s\t%s\n" "$S" "$url" >> "$OUT_URLS_TSV"
    found=$((found+1))
  else
    echo "$S" >> "$OUT_MISSING"
    miss=$((miss+1))
  fi
  # gentle pacing to avoid hammering the API
  sleep 0.15
  # lightweight progress every 1000
  if (( n % 1000 == 0 )); then
    echo "   …processed $n (urls: $found, missing: $miss)"
  fi
done < "$OUT_IDS"

# De-duplicate (just in case)
LC_ALL=C sort -u -o "$OUT_URLS" "$OUT_URLS"
LC_ALL=C sort -u -o "$OUT_URLS_TSV" "$OUT_URLS_TSV"
LC_ALL=C sort -u -o "$OUT_MISSING" "$OUT_MISSING"

echo "[3/4] URL list written:"
echo "   URL-only : $OUT_URLS  (count: $(wc -l < "$OUT_URLS"))"
echo "   TSV map  : $OUT_URLS_TSV"
echo "   Missing  : $OUT_MISSING  (count: $(wc -l < "$OUT_MISSING"))"

echo "[4/4] Download directory prepared:"
echo "   $DL_DIR"
echo "   (Use wget/curl or your parallel downloader to fetch into this folder.)"

echo "Done."
