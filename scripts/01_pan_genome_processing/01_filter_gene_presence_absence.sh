#!/usr/bin/env bash
#BSUB -J filter_rtab
#BSUB -o /data/pam/team230/sm71/scratch/rp2/logs/filter_rtab.%J.out
#BSUB -e /data/pam/team230/sm71/scratch/rp2/logs/filter_rtab.%J.err
#BSUB -n 1
#BSUB -M 8000
#BSUB -R "select[mem>8000] rusage[mem=8000]"

set -euo pipefail
export LC_ALL=C   # faster parsing

META_DIR="/data/pam/team230/sm71/scratch/rp2/metadata"
IN_DIR="/data/pam/team230/sm71/scratch/rp2/panaroo_output"
OUT_DIR="/data/pam/team230/sm71/scratch/rp2/panaroo_output"

FILE4="$META_DIR/File4_QC_characterisation_661K.tsv"
RTAB="$IN_DIR/gene_presence_absence.Rtab"
OUT="$OUT_DIR/gene_presence_absence_filtered.Rtab"
SUMMARY="$OUT_DIR/rtab_filter_summary.txt"

mkdir -p "$OUT_DIR"

# --- 18 target species (1 header + 18 rows = 19 lines) ---
TARGET="$META_DIR/target_species.tsv"
cat > "$TARGET" <<'EOF'
species
Staphylococcus aureus
Staphylococcus epidermidis
Staphylococcus pseudintermedius
Staphylococcus haemolyticus
Staphylococcus capitis
Staphylococcus sciuri
Staphylococcus argenteus
Streptococcus pneumoniae
Streptococcus pyogenes
Streptococcus agalactiae
Streptococcus suis
Streptococcus equi
Streptococcus dysgalactiae
Streptococcus uberis
Streptococcus mutans
Streptococcus mitis
Enterococcus faecium
Enterococcus faecalis
EOF

# --- Build keep-set from File4: (species ∧ high_quality==TRUE) ---
KEEP_IDS="$OUT_DIR/.keep_ids.txt"
TMP_SPEC_COUNT="$OUT_DIR/.tmp_species_candidates.count"
TMP_HQ_REMOVED="$OUT_DIR/.tmp_hq_removed.count"

awk -F'\t' -v KEEP="$KEEP_IDS" -v SPEC="$TMP_SPEC_COUNT" -v HQREM="$TMP_HQ_REMOVED" '
FNR==NR {
  # TARGET: single col "species" (skip header)
  if (NR>1 && $1!="") tgt[$1]=1
  next
}
FNR==1 { next }  # skip FILE4 header
{
  # Assume: FILE4 col1=sample_id, col2=species, last column=high_quality
  sid = $1
  sp  = $2
  hq  = $NF
  if (sid=="" || sp=="") next

  if (sp in tgt) {
    spec_candidates++
    if (hq=="TRUE") {          # strict TRUE; change to tolower(hq)=="true" for case-insensitive
      keep[sid]=1
    } else {
      hq_removed++
    }
  }
}
END {
  for (k in keep) print k > KEEP
  print (spec_candidates+0) > SPEC
  print (hq_removed+0) > HQREM
}
' "$TARGET" "$FILE4"

# Dedup keep IDs
sort -u -o "$KEEP_IDS" "$KEEP_IDS"

# --- Counts from mapping step ---
species_candidates=$(cat "$TMP_SPEC_COUNT")
hq_removed=$(cat "$TMP_HQ_REMOVED")
post_hq_candidates=$(wc -l < "$KEEP_IDS")

# --- Original RTAB counts (before any filtering) ---
orig_samples=$(head -n1 "$RTAB" | awk -F'\t' '{print NF-1}')
orig_genes=$(( $(wc -l < "$RTAB") - 1 ))

# --- How many of the post-HQ candidates actually exist in RTAB header? ---
matched_in_header=$(
  awk -F'\t' -v KEEP="$KEEP_IDS" '
  BEGIN { while ((getline k < KEEP) > 0) keep[k]=1; close(KEEP) }
  NR==1 { m=0; for (i=2;i<=NF;i++) if ($i in keep) m++; print m; exit }
  ' "$RTAB"
)
header_nontarget=$(( orig_samples - matched_in_header ))

# --- Filter RTAB (no pv) + count all-zero drops ---
GENE_DROP_FILE="$OUT_DIR/.tmp_genes_dropped.count"
AWK_BIN=$(command -v mawk || command -v awk)

"$AWK_BIN" -F'\t' -v OFS='\t' -v KEEP="$KEEP_IDS" -v DROPFILE="$GENE_DROP_FILE" '
BEGIN{
  while ((getline k < KEEP) > 0) keep[k]=1
  close(KEEP)
  dropped=0
}
NR==1{
  m=1; cols[m]=1           # always keep "Gene"
  for (i=2;i<=NF;i++) if ($i in keep) cols[++m]=i

  if (m==1) {
    print "ERROR: 0 columns matched keep-set. Check species/HQ filtering vs RTAB header." > "/dev/stderr"
    exit 2
  }

  # header
  for (j=1;j<=m;j++) printf "%s%s", $(cols[j]), (j<m?OFS:ORS)
  next
}
{
  present=0
  for (j=2;j<=m;j++){ v=$(cols[j]); if (v!="" && v!="0"){ present=1; break } }
  if (!present){ dropped++; next }

  for (j=1;j<=m;j++) printf "%s%s", $(cols[j]), (j<m?OFS:ORS)
}
END{
  print dropped > DROPFILE
}
' "$RTAB" > "$OUT"

# --- Filtered counts (after species+HQ, and zero-row drop) ---
filt_samples=$(head -n1 "$OUT" | awk -F'\t' '{print NF-1}')
filt_genes=$(( $(wc -l < "$OUT") - 1 ))
genes_dropped_all_zero=$(cat "$GENE_DROP_FILE")

# --- Write summary file ---
{
  echo -e "Metric\tCount"
  echo -e "Original_Rtab_Samples\t$orig_samples"
  echo -e "Original_Rtab_Genes\t$orig_genes"
  echo -e "Target_Species_Listed\t18"
  echo -e "Species_Candidates_in_File4\t$species_candidates"
  echo -e "Removed_by_high_quality_FALSE\t$hq_removed"
  echo -e "Candidates_after_HQ_TRUE\t$post_hq_candidates"
  echo -e "Matched_in_Rtab_Header\t$matched_in_header"
  echo -e "Header_NonTarget_Samples\t$header_nontarget"
  echo -e "Filtered_Rtab_Samples\t$filt_samples"
  echo -e "Genes_Dropped_AllZero\t$genes_dropped_all_zero"
  echo -e "Filtered_Rtab_Genes\t$filt_genes"
  echo -e "Filtered_Rtab_Path\t$OUT"
} > "$SUMMARY"

# --- Console output (human friendly) ---
printf "Kept sample IDs (post-HQ TRUE): %s\n" "$post_hq_candidates"
printf "Removed due to high_quality!=TRUE: %s\n" "$hq_removed"
printf "Header non-target (not in keep-set): %s\n" "$header_nontarget"
printf "Genes dropped (all-zero across kept samples): %s\n" "$genes_dropped_all_zero"
printf "Original Rtab: %s samples, %s genes\n" "$orig_samples" "$orig_genes"
printf "Filtered Rtab: %s samples, %s genes\n" "$filt_samples" "$filt_genes"
printf "Summary: %s\n" "$SUMMARY"
