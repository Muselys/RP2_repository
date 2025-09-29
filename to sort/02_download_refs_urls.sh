#!/usr/bin/env bash

# -----------------------------------------------------------------------------
# Parallel downloader for ENA contig URLs listed in a text file.
#
# What it does:
# - Reads URLs from urlbuild/resp/urls.txt (one per line).
# - Downloads each file into r_sequences/ using up to P workers (P = LSF slot count).
# - Skips files that already exist (resumable with wget -c when available).
# - Logs any failed URLs into logs/fails.txt (thread-safe via flock).
#
# Inputs:
#   - urlbuild/resp/urls.txt   (list of https://…contigs.fa.gz URLs)
#
# Outputs:
#   - r_sequences/<filename>   (downloaded contig archives)
#   - logs/fails.txt           (URLs that failed to download)
#
# Notes:
#   - Uses wget if present (quiet, continue, timeout/tries); falls back to curl.
#   - Converts Windows line endings in the URL list if dos2unix exists.
#   - Parallelism controlled by LSB_DJOB_NUMPROC or defaults to 20.
#   - All downloads happen on a single host (span[hosts=1]).
# -----------------------------------------------------------------------------


set -euo pipefail
export LC_ALL=C

BASE="/data/pam/team230/sm71/scratch/rp2/blast/references"
URLS="$BASE/urlbuild/resp/urls.txt"
OUT_DIR="$BASE/r_sequences"
LOG_DIR="$BASE/../logs"
FAILS="$LOG_DIR/fails.txt"
LOCK="$LOG_DIR/fails.lock"

mkdir -p "$OUT_DIR" "$LOG_DIR"
: > "$FAILS"

command -v dos2unix >/dev/null 2>&1 && dos2unix -q "$URLS" || true

# match -n to nslots (lightweight CPU usage)
P=${LSB_DJOB_NUMPROC:-20}

echo "[START] P=$P urls=$(wc -l < "$URLS")"

xargs -P "$P" -n 1 -I{} bash -c '
  url="{}"
  file=$(basename "$url")
  dest="'"$OUT_DIR"'/$file"

  if [ -e "$dest" ]; then
    echo "[SKIP] $file"
    exit 0
  fi

  echo "[DL] $file"
  if command -v wget >/dev/null 2>&1; then
    if ! wget -q -c --timeout=120 --tries=2 -O "$dest" "$url"; then
      (
        flock 9
        echo "$url" >> "'"$FAILS"'"
      ) 9>>"'"$LOCK"'"
      echo "[FAIL] $file"
      exit 1
    fi
  else
    if ! curl -s -L --retry 2 --retry-delay 3 -o "$dest" "$url"; then
      (
        flock 9
        echo "$url" >> "'"$FAILS"'"
      ) 9>>"'"$LOCK"'"
      echo "[FAIL] $file"
      exit 1
    fi
  fi
' _ < "$URLS"

echo "[DONE] Downloads attempted."
echo "[INFO] Fail list: $FAILS  (count: $(wc -l < "$FAILS"))"


S="SAMN10396905"  # <- put your sample (SAMEA*/SAMD*/SRS*/SRR*)
API="https://www.ebi.ac.uk/ena/portal/api/search"
curl -fsS -G "$API" \
  --data-urlencode 'result=analysis' \
  --data-urlencode 'format=tsv' \
  --data-urlencode 'limit=0' \
  --data-urlencode 'fields=analysis_accession,sample_accession,analysis_type,submitted_ftp' \
  --data-urlencode "query=(sample_accession=\"$S\" OR secondary_sample_accession=\"$S\") AND analysis_type=\"SEQUENCE_ASSEMBLY\"" \
| awk -F'\t' 'NR>1 && $4!=""{split($4,a,";"); for(i in a) if(a[i] ~ /(contigs|scaffolds).*\.(fa|fna|fasta)(\.gz)?$/){gsub(/^ftp:\/\//,"https://",a[i]); print a[i]; exit}}'
