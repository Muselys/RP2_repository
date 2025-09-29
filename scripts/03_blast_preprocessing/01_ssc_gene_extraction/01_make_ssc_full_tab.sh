#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# make_ssc_full_tab.sh
#
# Produces (in OUTDIR):
#   species_specific_core.tab  (KEEPS ALL ORIGINAL COLUMNS; gene_name expanded per token)
#   ssc_counts.txt             (totals + per-species before/after)
#
# Behavior:
#   - Filters rows to specific_class == "Species specific core" (case-sensitive)
#   - Expands clustered gene_name on "~~~" into individual rows
#   - Parses species from details (Core: ...), ignoring Inter:/Rare:
#
# Inputs (override via env):
#   IN_CLASS=/data/pam/team230/sm71/scratch/rp2/twilight_analysis/classification.tab
#   OUTDIR=/data/pam/team230/sm71/scratch/rp2/blast_preprocessing/ssc_queries
# -----------------------------------------------------------------------------
#BSUB -q normal
#BSUB -J make_ssc_full_tab
#BSUB -n 1
#BSUB -M 8000
#BSUB -R "select[mem>8000] rusage[mem=8000]"
#BSUB -o /data/pam/team230/sm71/scratch/rp2/logs/make_ssc_full_tab.%J.out
#BSUB -e /data/pam/team230/sm71/scratch/rp2/logs/make_ssc_full_tab.%J.err

set -euo pipefail
export LC_ALL=C

IN_CLASS="${IN_CLASS:-/data/pam/team230/sm71/scratch/rp2/twilight_analysis/classification.tab}"
OUTDIR="${OUTDIR:-/data/pam/team230/sm71/scratch/rp2/blast_preprocessing/ssc_queries}"
mkdir -p "$OUTDIR"

TAB_OUT="$OUTDIR/species_specific_core.tab"
STATS_OUT="$OUTDIR/ssc_counts.txt"

AWK_BIN="$(command -v mawk || command -v awk)"

echo "[STAGE 1] Starting make_ssc_full_tab.sh"
echo "[INFO] Input  : $IN_CLASS"
echo "[INFO] Outdir : $OUTDIR"
echo "[INFO] Output : $TAB_OUT , $STATS_OUT"
echo "[STAGE 2] Running AWK processing..."

"$AWK_BIN" -F'\t' -v OFS='\t' -v TAB="$TAB_OUT" '
NR==1{
  # capture header and column indices
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
  # write header unchanged
  for(i=1;i<=NF;i++) printf "%s%s", H[i], (i==NF?ORS:OFS) > TAB
  next
}
{
  sub(/\r$/, "")  # handle CRLF
  if($sc != "Species specific core") next

  ssc_rows_cls++  # rows retained from classification.tab

  # Parse species list from details -> Core: ... (strip Inter: and Rare:)
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
  # before = once per species for this input row
  for(s in row_sp) before[s]++

  # Expand gene_name tokens
  n=split($gn, a, /~~~/)

  for (i=1; i<=n; i++) {
    t=a[i]; gsub(/^[[:space:]]+|[[:space:]]+$/, "", t)
    if (t=="") continue

    # after counts incremented per species
    for(s in row_sp){ after[s]++; total_after++ }

    # Emit row with SAME columns, but gene_name replaced by token t
    for(j=1;j<=NF;j++){
      val = (j==gn ? t : $j)
      printf "%s%s", val, (j==NF?ORS:OFS) >> TAB
    }
    total_rows_tab++
  }
}
END{
  # Totals
  printf "ssc rows in classification.tab: %d\n", (ssc_rows_cls+0)
  printf "rows in species_specific_core.tab (after expansion): %d\n", (total_rows_tab+0)
  printf "expanded gene tokens total: %d\n", (total_after+0)

  # Per-species
  print ""
  print "species\tbefore_rows\tafter_expanded_genes"
  for(s in seen_sp){
    printf "%s\t%d\t%d\n", s, (before[s]+0), (after[s]+0)
  }
}
' "$IN_CLASS" | tee "$STATS_OUT"

echo "[STAGE 3] AWK processing complete"
echo "[STAGE 4] species_specific_core.tab written to: $TAB_OUT"
echo "[STAGE 5] ssc_counts.txt written to: $STATS_OUT"
echo "[DONE] make_ssc_full_tab.sh finished successfully"
