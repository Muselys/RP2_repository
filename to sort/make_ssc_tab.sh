#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# make_species_specific_core.sh
#
# Builds species_specific_core.tab and ssc_counts.txt
# - "Before" = rows per species in species_specific_core.tab (not expanded)
# - "After"  = expanded gene count = sum of tokens after splitting gene_name on "~~~"
#
# Inputs (override via env):
#   IN_CLASS=/data/pam/team230/sm71/scratch/rp2/twilight_analysis/classification.tab
#   OUTDIR=/data/pam/team230/sm71/scratch/rp2/blast_preprocessing
#
# Outputs (in OUTDIR):
#   species_specific_core.tab         (gene_name, species)
#   ssc_counts.txt                    (totals + per-species before/after)
# -----------------------------------------------------------------------------
#BSUB -q normal
#BSUB -J make_ssc_tab
#BSUB -n 1
#BSUB -M 8000
#BSUB -R "select[mem>8000] rusage[mem=8000]"
#BSUB -o /data/pam/team230/sm71/scratch/rp2/logs/make_ssc_tab.%J.out
#BSUB -e /data/pam/team230/sm71/scratch/rp2/logs/make_ssc_tab.%J.err

set -euo pipefail
export LC_ALL=C

IN_CLASS="${IN_CLASS:-/data/pam/team230/sm71/scratch/rp2/twilight_analysis/classification.tab}"
OUTDIR="${OUTDIR:-/data/pam/team230/sm71/scratch/rp2/blast_preprocessing}"
mkdir -p "$OUTDIR"

RAW_CORE="$OUTDIR/species_specific_core.tab"
STATS_OUT="$OUTDIR/ssc_counts.txt"

AWK_BIN="$(command -v mawk || command -v awk)"

echo "[INFO] Input : $IN_CLASS"
echo "[INFO] Outdir: $OUTDIR"

"$AWK_BIN" -F'\t' -v OFS='\t' -v RAW="$RAW_CORE" '
NR==1{
  # detect columns
  for(i=1;i<=NF;i++){
    h=$i; gsub(/^[ \t]+|[ \t]+$/, "", h)
    if(h=="specific_class") sc=i
    if(h=="gene_name")      gn=i
    if(h=="details")        det=i
  }
  if(!sc || !gn || !det){
    print "ERROR: required columns missing (specific_class, gene_name, details)" > "/dev/stderr"
    exit 3
  }
  print "gene_name","species" > RAW
  next
}
{
  sub(/\r$/, "")  # handle CRLF
  if($sc != "Species specific core") next

  ssc_rows_cls++  # SSC rows in classification.tab

  # token count for this gene_name (used for "after expansion")
  n=split($gn, a, /~~~/)
  tok_n=0
  for(i=1;i<=n;i++){
    t=a[i]; gsub(/^[[:space:]]+|[[:space:]]+$/, "", t)
    if(t!="") tok_n++
  }

  # parse species list from details
  d=$det
  sub(/.*Core:[[:space:]]*/, "", d)
  sub(/[[:space:]]+Inter:.*$/, "", d)
  sub(/[[:space:]]+Rare:.*$/,  "", d)
  gsub(/^[[:space:]]+|[[:space:]]+$/, "", d)
  if(d=="") next

    m=split(d, sp, /\+/)
  for (k=1; k<=m; k++) {
    s=sp[k]; gsub(/^[[:space:]]+|[[:space:]]+$/, "", s)
    if (s=="") continue

    # keep the "before" (pre-expansion) row count per species
    before[s]++
    seen_sp[s]=1

    # write one line per expanded token and count them
    for (i=1; i<=n; i++) {
      t=a[i]; gsub(/^[[:space:]]+|[[:space:]]+$/, "", t)
      if (t=="") continue
      print t, s >> RAW
      after[s]++        # per-species expanded count
      total_rows_raw++  # lines written to species_specific_core.tab (expanded)
      total_after++     # total expanded tokens across species rows
    }
  }

END{
  # Totals
  printf "ssc rows in classification.tab: %d\n", (ssc_rows_cls+0)
  printf "ssc rows in species_specific_core.tab: %d\n", (total_rows_raw+0)
  printf "expanded gene count (rows in species_specific_core.tab after expanded): %d\n", (total_after+0)

  # Per-species
  print ""
  print "species\tbefore_rows\tafter_expanded_genes"
  for(s in seen_sp){
    printf "%s\t%d\t%d\n", s, (before[s]+0), (after[s]+0)
  }
}
' "$IN_CLASS" | tee "$STATS_OUT"

echo "[DONE] species_specific_core.tab -> $RAW_CORE"
echo "[DONE] ssc_counts.txt -> $STATS_OUT"
