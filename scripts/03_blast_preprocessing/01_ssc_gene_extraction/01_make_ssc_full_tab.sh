#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# make_ssc_full_tab.sh — extract species-specific core genes
#
# Produces (in OUTDIR):
#   species_specific_core.tab  (keeps all original columns; expands gene_name per token)
#   ssc_counts.txt             (totals + per-species before/after)
#
# Behaviour:
#   - Filters rows with specific_class == "Species specific core"
#   - Expands clustered gene_name on "~~~" into individual rows
#   - Parses species list from details column (Core: …), ignoring Inter:/Rare:
#
# Dependencies:
#   awk/mawk, bash ≥4
#
# Example:
#   ./make_ssc_full_tab.sh --class twilight_analysis/classification.tab \
#                          --outdir blast_preprocessing/ssc_queries
# -----------------------------------------------------------------------------
#BSUB -q normal
#BSUB -J make_ssc_full_tab
#BSUB -n 1
#BSUB -M 8000
#BSUB -R "select[mem>8000] rusage[mem=8000]"
#BSUB -o logs/make_ssc_full_tab.%J.out
#BSUB -e logs/make_ssc_full_tab.%J.err

set -euo pipefail
export LC_ALL=C

usage() {
  cat <<EOF
Usage: $0 --class classification.tab --outdir output_dir

Options:
  --class    Path to Twilight classification.tab (required)
  --outdir   Output directory (required)
  --help     Show this message
EOF
}

IN_CLASS=""
OUTDIR=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --class) IN_CLASS="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    --help|-h) usage; exit 0;;
    *) echo "Unknown option: $1"; usage; exit 2;;
  esac
done

[[ -n "$IN_CLASS" && -f "$IN_CLASS" ]] || { echo "ERROR: --class file missing"; exit 2; }
[[ -n "$OUTDIR" ]] || { echo "ERROR: --outdir required"; exit 2; }

mkdir -p "$OUTDIR"
TAB_OUT="$OUTDIR/species_specific_core.tab"
STATS_OUT="$OUTDIR/ssc_counts.txt"
MANIFEST="$OUTDIR/ssc_manifest.txt"

AWK_BIN="$(command -v mawk || command -v awk)"

# Provenance log
{
  echo "# make_ssc_full_tab.sh v1.0.0"
  echo "# date_utc: $(date -u +"%Y-%m-%dT%H:%M:%SZ")"
  command -v md5sum >/dev/null && md5sum "$IN_CLASS" | awk '{print "input_md5\t"$1}'
  echo -e "input_file\t$IN_CLASS"
  echo -e "awk_bin\t$AWK_BIN"
} > "$MANIFEST"

echo "[INFO] Input : $IN_CLASS"
echo "[INFO] Outdir: $OUTDIR"
echo "[INFO] Output: $TAB_OUT , $STATS_OUT"

"$AWK_BIN" -F'\t' -v OFS='\t' -v TAB="$TAB_OUT" '
NR==1{
  for(i=1;i<=NF;i++){
    h=$i; gsub(/^[ \t]+|[ \t]+$/, "", h)
    H[i]=h
    if(h=="specific_class") sc=i
    if(h=="gene_name")      gn=i
    if(h=="details")        det=i
  }
  if(!sc || !gn || !det){
    print "ERROR: required columns missing (specific_class, gene_name, details)" > "/dev/stderr"
    exit 3
  }
  for(i=1;i<=NF;i++) printf "%s%s", H[i], (i==NF?ORS:OFS) > TAB
  next
}
{
  sub(/\r$/, "")
  if($sc != "Species specific core") next
  ssc_rows_cls++
  d=$det
  sub(/.*Core:[[:space:]]*/, "", d)
  sub(/[[:space:]]+Inter:.*$/, "", d)
  sub(/[[:space:]]+Rare:.*$/,  "", d)
  gsub(/^[[:space:]]+|[[:space:]]+$/, "", d)
  if(d=="") next
  m=split(d, sp, /\+/)
  delete row_sp
  for (k=1; k<=m; k++) {
    s=sp[k]; gsub(/^[[:space:]]+|[[:space:]]+$/, "", s)
    if (s!=""){ row_sp[s]=1; seen_sp[s]=1 }
  }
  for(s in row_sp) before[s]++
  n=split($gn, a, /~~~/)
  for (i=1; i<=n; i++) {
    t=a[i]; gsub(/^[[:space:]]+|[[:space:]]+$/, "", t)
    if (t=="") continue
    for(s in row_sp){ after[s]++; total_after++ }
    for(j=1;j<=NF;j++){
      val = (j==gn ? t : $j)
      printf "%s%s", val, (j==NF?ORS:OFS) >> TAB
    }
    total_rows_tab++
  }
}
END{
  printf "ssc rows in classification.tab: %d\n", (ssc_rows_cls+0)
  printf "rows in species_specific_core.tab (after expansion): %d\n", (total_rows_tab+0)
  printf "expanded gene tokens total: %d\n", (total_after+0)
  print ""
  print "species\tbefore_rows\tafter_expanded_genes"
  for(s in seen_sp){
    printf "%s\t%d\t%d\n", s, (before[s]+0), (after[s]+0)
  }
}
' "$IN_CLASS" | tee "$STATS_OUT"

echo "[DONE] species_specific_core.tab and ssc_counts.txt written to $OUTDIR"
