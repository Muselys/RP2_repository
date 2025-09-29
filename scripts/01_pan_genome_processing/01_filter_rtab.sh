#!/usr/bin/env bash

# -----------------------------------------------------------------------------
# filter_rtab.sh — filter Panaroo/Roary Rtab to target species & HQ genomes
#
# What this script does
# 1) Creates metadata/target_species.tsv (hard-coded list of 18 species).
# 2) Reads File4_QC_characterisation_661K.tsv (TSV; col1=sample_id, col2=species,
#    last column=high_quality) and builds a keep-set:
#       keep_ids = { sample_id | species ∈ target_species AND high_quality == TRUE }.
#    (IDs are de-duplicated.)
# 3) Reads gene_presence_absence.Rtab and writes a filtered matrix that keeps:
#       • the first column ("Gene"), and
#       • only those sample columns whose header matches keep_ids.
# 4) Drops any gene row that is all zeros/empty across the kept sample columns.
# 5) Outputs:
#       • gene_presence_absence_filtered.Rtab  (filtered matrix)
#       • rtab_filter_summary.txt              (counts + output path)
# -----------------------------------------------------------------------------

#BSUB -J filter_rtab
#BSUB -o /data/pam/team230/sm71/scratch/rp2/logs/filter_rtab.%J.out
#BSUB -e /data/pam/team230/sm71/scratch/rp2/logs/filter_rtab.%J.err
#BSUB -n 1
#BSUB -M 4000
#BSUB -R "select[mem>4000] rusage[mem=4000]"

set -euo pipefail
export LC_ALL=C   # faster parsing

META_DIR="/data/pam/team230/sm71/scratch/rp2/metadata"
IN_DIR="/data/pam/team230/sm71/scratch/rp2/panaroo_output"
OUT_DIR="/data/pam/team230/sm71/scratch/rp2/pan_genome_processing"

FILE4="$META_DIR/File4_QC_characterisation_661K.tsv"
RTAB="$IN_DIR/gene_presence_absence.Rtab"
OUT="$OUT_DIR/gene_presence_absence_filtered.Rtab"
SUMMARY="$OUT_DIR/rtab_filter_summary.txt"

echo "Stage 0: Setup — creating output directory and confirming paths"
mkdir -p "$OUT_DIR"
echo "Stage 0: META_DIR=$META_DIR"
echo "Stage 0: IN_DIR=$IN_DIR"
echo "Stage 0: OUT_DIR=$OUT_DIR"
echo "Stage 0: FILE4=$FILE4"
echo "Stage 0: RTAB=$RTAB"

# --- 18 target species (1 header + 18 rows = 19 lines) ---
TARGET="$META_DIR/target_species.tsv"

echo "Stage 1: Writing target species list → $TARGET"
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
echo "Stage 1: Done (18 target species listed)."

# --- Build set of sample IDs present in the RTAB header ---
TMP_HEADER="$OUT_DIR/.tmp_rtab_header_ids.txt"

echo "Stage 1.5: Extracting sample IDs from RTAB header -> $TMP_HEADER"
head -n1 "$RTAB" | cut -f2- | tr '\t' '\n' | sort -u > "$TMP_HEADER"
echo "Stage 1.5: Header IDs: $(wc -l < "$TMP_HEADER")"

# --- Build keep-set from File4, but only for IDs also present in RTAB header ---
KEEP_IDS="$OUT_DIR/.keep_ids.txt"
TMP_SPEC_COUNT="$OUT_DIR/.tmp_species_candidates.count"
TMP_HQ_REMOVED="$OUT_DIR/.tmp_hq_removed.count"


echo "Stage 2: Building keep-set (target species ∧ high_quality==TRUE) ∩ Rtab header"
awk -F'\t' -v KEEP="$KEEP_IDS" -v SPEC="$TMP_SPEC_COUNT" -v HQREM="$TMP_HQ_REMOVED" '
# 1) load allowlist species
FNR==NR {
  if (NR>1 && $1!="") tgt[$1]=1
  next
}
# 2) load RTAB header sample IDs (one per line)
FILENAME==ARGV[2] {
  hdr[$1]=1
  next
}
# 3) scan File4 and count only those rows whose sample_id is in hdr
FILENAME==ARGV[3] {
  if (FNR==1) next  # skip header
  sid=$1; sp=$2; hq=$NF
  if (sid=="" || sp=="") next
  if ((sp in tgt) && (sid in hdr)) {
    spec_candidates++            # candidates *that actually exist in the RTAB header*
    if (tolower(hq)=="true") {
      keep[sid]=1
    } else {
      hq_removed++
    }
  }
  next
}
END {
  for (k in keep) print k > KEEP
  print (spec_candidates+0) > SPEC
  print (hq_removed+0) > HQREM
}
' "$TARGET" "$TMP_HEADER" "$FILE4"

echo "Stage 3: De-duplicating keep IDs"
sort -u -o "$KEEP_IDS" "$KEEP_IDS"

# --- Counts from mapping step ---
species_candidates=$(cat "$TMP_SPEC_COUNT")
hq_removed=$(cat "$TMP_HQ_REMOVED")
post_hq_candidates=$(wc -l < "$KEEP_IDS")
echo "Stage 3: Candidates present in RTAB header (target species): $species_candidates"
echo "Stage 3: Removed by high_quality!=TRUE: $hq_removed"
echo "Stage 3: Kept sample IDs after HQ filter: $post_hq_candidates"

# --- Original RTAB counts (before any filtering) ---
echo "Stage 4: Inspecting original RTAB dimensions"
orig_samples=$(head -n1 "$RTAB" | awk -F'\t' '{print NF-1}')
orig_genes=$(( $(wc -l < "$RTAB") - 1 ))
echo "Stage 4: Original RTAB — samples: $orig_samples, genes: $orig_genes"

# --- How many of the post-HQ candidates actually exist in RTAB header? ---
echo "Stage 5: Matching keep-set to RTAB header"
matched_in_header=$(
  awk -F'\t' -v KEEP="$KEEP_IDS" '
  BEGIN { while ((getline k < KEEP) > 0) keep[k]=1; close(KEEP) }
  NR==1 { m=0; for (i=2;i<=NF;i++) if ($i in keep) m++; print m; exit }
  ' "$RTAB"
)
header_nontarget=$(( orig_samples - matched_in_header ))
echo "Stage 5: Matched in RTAB header: $matched_in_header"
echo "Stage 5: Header non-target samples: $header_nontarget"

# --- Filter RTAB + count all-zero drops ---
GENE_DROP_FILE="$OUT_DIR/.tmp_genes_dropped.count"
AWK_BIN=$(command -v mawk || command -v awk)
echo "Stage 6: Filtering RTAB (columns ∈ keep-set + drop all-zero gene rows)"
echo "Stage 6: Using AWK binary: $AWK_BIN"
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
echo "Stage 6: Filtering complete → $OUT"

# --- Filtered counts (after species+HQ, and zero-row drop) ---
echo "Stage 7: Computing filtered RTAB dimensions"
filt_samples=$(head -n1 "$OUT" | awk -F'\t' '{print NF-1}')
filt_genes=$(( $(wc -l < "$OUT") - 1 ))
genes_dropped_all_zero=$(cat "$GENE_DROP_FILE")
echo "Stage 7: Filtered RTAB — samples: $filt_samples, genes: $filt_genes"
echo "Stage 7: Genes dropped (all-zero across kept samples): $genes_dropped_all_zero"

# --- Write summary file ---
echo "Stage 8: Writing summary → $SUMMARY"
{
  echo -e "Metric\tCount"
  echo -e "Original_Rtab_Samples\t$orig_samples"
  echo -e "Original_Rtab_Genes\t$orig_genes"
  echo -e "Target_Species_Listed\t18"
  echo -e "Species_Candidates_in_RtabHeader\t$species_candidates"
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
echo "Stage 9: Done — final recap"
printf "Kept sample IDs (post-HQ TRUE): %s\n" "$post_hq_candidates"
printf "Removed due to high_quality!=TRUE: %s\n" "$hq_removed"
printf "Header non-target (not in keep-set): %s\n" "$header_nontarget"
printf "Genes dropped (all-zero across kept samples): %s\n" "$genes_dropped_all_zero"
printf "Original Rtab: %s samples, %s genes\n" "$orig_samples" "$orig_genes"
printf "Filtered Rtab: %s samples, %s genes\n" "$filt_samples" "$filt_genes"
printf "Summary: %s\n" "$SUMMARY"
