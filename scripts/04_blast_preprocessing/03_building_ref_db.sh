#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# Builds HTTPS download URLs for contig FASTA files from ENA using sample IDs
# pulled from a Panaroo Rtab header, in safe chunks.
#
# What it does:
# - Inputs:
#     * gene_presence_absence.Rtab (Panaroo) — uses header to get sample IDs
# - Steps:
#     1) Grab all sample IDs from the Rtab header (everything after the first "gene" column).
#     2) Split the IDs into chunks (default 200) to avoid huge API requests.
#     3) Query ENA’s analysis endpoint (analysis_type=SEQUENCE_ASSEMBLY) per chunk.
#     4) Combine all ENA TSV responses into one file.
#     5) From the ENA “submitted_ftp” field, extract exactly one URL per sample
#        that ends with “.contigs.fa.gz” (strict filter), and convert ftp:// → https://.
#     6) Write both a SAMPLE↦URL TSV and a plain URL list for quick download tests.
#     7) Report which Panaroo header sample IDs did not get a contig URL.
#
# Outputs:
#   - all_sequence_assemblies.tsv : concatenated ENA responses (with header)
#   - contig_urls.tsv              : SAMPLE<TAB>https://…contigs.fa.gz (unique by sample)
#   - contig_urls.txt              : https://…contigs.fa.gz (URL-only list)
#   - missing_samples.txt          : sample IDs present in Panaroo header but no URL found
#
# Notes:
# - Uses ENA API fields: analysis_accession, sample_accession, analysis_type, submitted_ftp.
# - Only accepts entries ending with “contigs.fa.gz” (ignores other files in submitted_ftp).
# - Retries and chunking make it network- and memory-safe for large cohorts.
# - Prints summary counts (header samples, URLs found, missing) for quick QC.
# -----------------------------------------------------------------------------

#BSUB -q normal
#BSUB -J build_urls_all
#BSUB -n 1
#BSUB -M 2000
#BSUB -R "select[mem>2000] rusage[mem=2000]"
#BSUB -o /data/pam/team230/sm71/scratch/rp2/blast/logs/build_urls_all.%J.out
#BSUB -e /data/pam/team230/sm71/scratch/rp2/blast/logs/build_urls_all.%J.err

#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C


# --- config ---
PANAROO="/data/pam/team230/sm71/scratch/rp2/panaroo_output/gene_presence_absence.Rtab"

# Put all working dirs *inside* the current references directory
BASE="/data/pam/team230/sm71/scratch/rp2/blast/references/urlbuild"
API="https://www.ebi.ac.uk/ena/portal/api/search"
FIELDS="analysis_accession,sample_accession,analysis_type,submitted_ftp"
CHUNK_SIZE=200

# outputs
CSV="$BASE/reference_sampleids.csv"
WORK="$BASE/chunks"
RESP="$BASE/resp"
ALL_TSV="$RESP/all_sequence_assemblies.tsv"
URLS_TSV="$RESP/contig_urls.tsv"   # SAMPLE<TAB>URL (unique by sample)
URLS_TXT="$RESP/contig_urls.txt"   # URL only


mkdir -p "$BASE" "$WORK" "$RESP"

require() { command -v "$1" >/dev/null 2>&1 || { echo "[ERR] $1 not found"; exit 1; }; }
require curl; require awk; require tr; require sed; require split

# helper: join IDs as "A","B","C"
make_in_list() {
  awk '{printf (NR==1?"":"\",\"") $0} END{print ""}' "$1" \
  | sed '1s/^/"/; $s/$/"/'
}

echo "[STEP 1] Extract ALL sample IDs from Panaroo header -> $CSV"
head -n1 "$PANAROO" | cut -f2- | tr '\t' ',' > "$CSV"
SAMPLE_COUNT=$(awk -F',' '{print NF}' "$CSV")
echo "[INFO] Panaroo sample count (header): $SAMPLE_COUNT"

echo "[STEP 2] CSV -> one-per-line; split into chunks of $CHUNK_SIZE"
IDS="$BASE/sample_ids.txt"
tr ',' '\n' < "$CSV" \
 | sed 's/^[[:space:]]*//; s/[[:space:]]*$//' \
 | awk 'NF' > "$IDS"
echo "[INFO] IDs written: $(wc -l < "$IDS")"

# fresh workspace
find "$WORK" -maxdepth 1 -type f -name 'ids_chunk_*' -delete || true
find "$RESP" -maxdepth 1 -type f -name 'ids_chunk_*.tsv' -delete || true
: > "$ALL_TSV" || true

split -l "$CHUNK_SIZE" "$IDS" "$WORK/ids_chunk_"
echo "[INFO] chunks: $(ls -1 "$WORK"/ids_chunk_* | wc -l)"

echo "[STEP 3] Query ENA per chunk (analysis_type=SEQUENCE_ASSEMBLY)"
HEADER_WRITTEN=0
for CHUNK in "$WORK"/ids_chunk_*; do
  n=$(wc -l < "$CHUNK"); base=$(basename "$CHUNK")
  OUT_TSV="$RESP/${base}.tsv"
  echo "[INFO] querying $n IDs in $base"

  QLIST=$(make_in_list "$CHUNK")
  RAW_QUERY="sample_accession in (${QLIST}) AND analysis_type=\"SEQUENCE_ASSEMBLY\""

  curl -fsS --retry 5 --retry-delay 5 --retry-all-errors \
       --max-time 180 -G "$API" \
       --data-urlencode "result=analysis" \
       --data-urlencode "format=tsv" \
       --data-urlencode "limit=0" \
       --data-urlencode "fields=${FIELDS}" \
       --data-urlencode "query=${RAW_QUERY}" \
       > "$OUT_TSV"

  if [ $HEADER_WRITTEN -eq 0 ]; then
    cat "$OUT_TSV" > "$ALL_TSV"; HEADER_WRITTEN=1
  else
    awk 'NR>1' "$OUT_TSV" >> "$ALL_TSV"
  fi
done


echo "[STEP 3 DONE] Combined TSV -> $ALL_TSV"
echo "[STATS] rows excl header: $(awk 'END{print NR-1}' "$ALL_TSV")"

echo "[STEP 4] Build unique contig URL list (SAMPLE<TAB>https://…) -> $URLS_TSV"
# Only accept *.contigs.fa.gz files (strict requirement)
awk -F'\t' '
  function emit(sample,u){
    gsub(/^ftp:\/\//,"",u)
    printf "%s\thttps://%s\n", sample, u
  }
  NR>1 && $4!="" {
    if (seen[$2]++) next
    split($4,a,";")

    # strict: only contigs.fa.gz
    for(i in a) if (a[i] ~ /contigs\.fa\.gz$/) { emit($2,a[i]); next }
  }
' "$ALL_TSV" > "$URLS_TSV"

cut -f2 "$URLS_TSV" > "$URLS_TXT"

echo "[DONE] URL build complete."
echo "  - ENA rows:          $ALL_TSV"
echo "  - Contig URLs (TSV): $URLS_TSV"
echo "  - Contig URLs (TXT): $URLS_TXT"
echo "  - Counts:"
echo "      header samples: $SAMPLE_COUNT"
echo "      unique URLs:    $(wc -l < "$URLS_TSV")"

echo "[STEP 5] Identify missing sample IDs (in Panaroo but no contig URL)"

MISSING="$RESP/missing_samples.txt"

# get header IDs from Panaroo (skip the "gene" col, split on tabs)
head -n1 "$PANAROO" | cut -f2- | tr '\t' '\n' | sed 's/^[[:space:]]*//; s/[[:space:]]*$//' | sort -u > "$RESP/header_ids.txt"

# extract sample IDs that we actually got URLs for
cut -f1 "$URLS_TSV" | sort -u > "$RESP/url_ids.txt"

# compare: header IDs - URL IDs = missing
comm -23 "$RESP/header_ids.txt" "$RESP/url_ids.txt" > "$MISSING"

echo "[DONE] Missing sample IDs written to $MISSING"
echo "  header samples : $(wc -l < "$RESP/header_ids.txt")"
echo "  with URLs      : $(wc -l < "$RESP/url_ids.txt")"
echo "  missing        : $(wc -l < "$MISSING")"