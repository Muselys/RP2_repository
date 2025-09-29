#!/usr/bin/env bash
set -euo pipefail

IN="$1"      # annotated .tsv or .tsv.gz
OUT="$2"     # .tsv.gz
PIDENT="${3:-80}"   # default 80
QCOVS="${4:-90}"    # default 90
THREADS="${THREADS:-8}"

# choose a reader based on extension
reader() {
  case "$IN" in
    *.gz) pigz -dc -- "$IN" ;;
    *)    cat -- "$IN" ;;
  esac
}

reader | awk -v OFS='\t' -v PIDENT="$PIDENT" -v QCOVS="$QCOVS" -F'\t' '
  NR==1 {
    for(i=1;i<=NF;i++){h[$i]=i}
    # require columns exist
    if(!("pident" in h) || !("qcovs" in h)){
      print "Missing pident/qcovs columns" > "/dev/stderr"; exit 1
    }
    pid=h["pident"]; qcv=h["qcovs"];
    print $0; next
  }
  {
    # numeric compare; coerce empty/non-numeric to 0/NaN-safe
    p=$pid+0; q=$qcv+0;
    if (p >= PIDENT && q >= QCOVS) print $0
  }
' | pigz -p "$THREADS" -1 > "$OUT"

#FILTER_THREADS=16 ./filter_thresholds.sh blast.annotated.tsv.gz blast.filtered.tsv.gz 80 90


