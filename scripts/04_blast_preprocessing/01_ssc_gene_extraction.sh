#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# This script extracts "Species specific core" (SSC) genes from a Twilight
# classification table and prepares them for downstream BLAST/QC analysis.
#
# Steps:
# 1. Filter classification.tab to keep only rows with "Species specific core".
# 2. Parse species information from the "details" column.
# 3. Expand gene_name fields split by "~~~" into individual gene entries.
# 4. Write a final SSC table (gene_name, species, specific_class).
# 5. Generate per-species counts of unique SSC genes for QC/summary.
#
# Outputs:
# - species_specific_core.tab   : expanded SSC gene table
# - genes_per_species_counts.tsv: per-species gene counts (sorted)
# -----------------------------------------------------------------------------
 
#BSUB -q normal
#BSUB -J twilight_ssc_extract
#BSUB -n 1
#BSUB -M 8000
#BSUB -R "select[mem>8000] rusage[mem=8000]"
#BSUB -W 00:10
#BSUB -o ssc_extract.%J.out
#BSUB -e ssc_extract.%J.err

set -euo pipefail
export LC_ALL=C

# ---- paths ----
in_class="/data/pam/team230/sm71/scratch/rp2/twilight_output/classification.tab"
outdir="/data/pam/team230/sm71/scratch/rp2/blast"
final_tab="$outdir/species_specific_core.tab"          # FINAL: gene_name  species  specific_class
counts_out="$outdir/genes_per_species_counts.tsv"

# ---- setup ----
[[ -s "$in_class" ]] || { echo "ERR: missing input $in_class" >&2; exit 1; }
mkdir -p "$outdir"

# temps
tmp_ssc="$(mktemp)"; tmp_final="$(mktemp)"
cleanup(){ rm -f "$tmp_ssc" "$tmp_final"; }
trap cleanup EXIT

########################################
# STEP 1: filter to "Species specific core" (TEMP)
########################################
awk -F'\t' -v OFS='\t' '
NR==1{
  for(i=1;i<=NF;i++){
    h=$i; gsub(/^[ \t]+|[ \t]+$/, "", h)
    if(h=="specific_class") sc=i
  }
  if(!sc){ print "ERROR: specific_class column not found" > "/dev/stderr"; exit 2 }
  print; next
}
{
  gsub(/\r$/,"")
  if($sc=="Species specific core") print
}' "$in_class" > "$tmp_ssc"

echo "Filtered SSC rows (temp): $(wc -l < "$tmp_ssc") lines (incl. header)"

########################################
# STEP 2+3: parse species from details AND expand gene_name -> FINAL
########################################
awk -F'\t' -v OFS='\t' '
NR==1{
  # find indices for needed columns in temp SSC table
  for(i=1;i<=NF;i++){
    h=$i
    if(h=="gene_name")      gn=i
    if(h=="details")        det=i
    if(h=="specific_class") sc=i
  }
  if(!gn || !det || !sc){ print "ERROR: required columns missing (gene_name/details/specific_class)" > "/dev/stderr"; exit 3 }
  print "gene_name","species","specific_class"; next
}
{
  # parse species from details
  d=$det
  sub(/\r$/,"", d)
  sub(/[[:space:]]+Inter:.*$/, "", d)
  sub(/[[:space:]]+Rare:.*$/,  "", d)
  sub(/^[[:space:]]*Core:[[:space:]]*/, "", d)
  sub(/[[:space:]]*\+.*$/, "", d)              # if "X+Y", keep "X"
  gsub(/^[[:space:]]+|[[:space:]]+$/, "", d)   # trim

  # expand gene_name on "~~~"
  n = split($gn, parts, /~~~/)
  for (i=1; i<=n; i++){
    g = parts[i]
    gsub(/^[[:space:]]+|[[:space:]]+$/, "", g)
    if (g != "") print g, d, $sc
  }
}
' "$tmp_ssc" > "$tmp_final"

# Overwrite the original SSC filename with the FINAL (expanded) table
mv -f "$tmp_final" "$final_tab"
echo "Wrote FINAL (expanded) -> $final_tab"
head -5 "$final_tab" || true

########################################
# STEP 4: count unique gene_name per species
########################################
awk -F'\t' -v OFS='\t' '
NR==1 { next }
{
  sp=$2; gene=$1
  key=sp SUBSEP gene
  if (!(key in seen)) { seen[key]=1; c[sp]++ }
}
END {
  print "species","gene_count"
  for (sp in c) print sp, c[sp]
}
' "$final_tab" | sort -k2,2nr > "$counts_out"
# For alphabetical instead:  | sort -k1,1

echo "Counts -> $counts_out"
head -10 "$counts_out" || true