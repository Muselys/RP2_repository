#!/usr/bin/env bash
# filter_rtab.sh — filter Panaroo/Roary Rtab to target species & HQ genomes (v1.0.0)

set -euo pipefail
export LC_ALL=C

usage() {
  cat <<EOF
Usage: $(basename "$0") \
  --rtab gene_presence_absence.Rtab \
  --metadata File4_QC_characterisation_661K.tsv \
  --outdir pan_genome_processing \
  [--species-file target_species.tsv]

Required:
  --rtab           Path to Panaroo/Roary gene_presence_absence.Rtab
  --metadata       File4_QC_characterisation_661K.tsv (Blackwell et al., 2021b)
  --outdir         Output directory

Optional:
  --species-file   TSV with header 'species' and 18 rows (default: built-in list)
  --help           Show this help
EOF
}

# --- parse args ---
RTAB="" ; META="" ; OUTDIR="" ; SPECIES_FILE=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --rtab) RTAB="$2"; shift 2;;
    --metadata) META="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    --species-file) SPECIES_FILE="$2"; shift 2;;
    --help|-h) usage; exit 0;;
    *) echo "ERROR: Unknown arg $1"; usage; exit 2;;
  esac
done

[[ -z "${RTAB}" || -z "${META}" || -z "${OUTDIR}" ]] && { usage; exit 2; }
[[ -f "$RTAB" ]] || { echo "ERROR: RTAB not found: $RTAB" >&2; exit 2; }
[[ -f "$META" ]] || { echo "ERROR: metadata not found: $META" >&2; exit 2; }
mkdir -p "$OUTDIR"

# Paths
OUT="$OUTDIR/gene_presence_absence_filtered.Rtab"
SUMMARY="$OUTDIR/rtab_filter_summary.txt"
TMP_HEADER="$OUTDIR/.tmp_rtab_header_ids.txt"
KEEP_IDS="$OUTDIR/.keep_ids.txt"
TMP_SPEC_COUNT="$OUTDIR/.tmp_species_candidates.count"
TMP_HQ_REMOVED="$OUTDIR/.tmp_hq_removed.count"
GENE_DROP_FILE="$OUTDIR/.tmp_genes_dropped.count"

cleanup() { rm -f "$TMP_HEADER" "$KEEP_IDS" "$TMP_SPEC_COUNT" "$TMP_HQ_REMOVED" "$GENE_DROP_FILE" 2>/dev/null || true; }
trap cleanup EXIT

# Species list
TARGET="${SPECIES_FILE:-$OUTDIR/target_species.tsv}"
if [[ -z "${SPECIES_FILE}" ]]; then
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
fi

echo "# filter_rtab.sh v1.0.0" | tee "$SUMMARY"
echo "# date: $(date -u +"%Y-%m-%dT%H:%M:%SZ")" | tee -a "$SUMMARY"
echo -e "# rtab_md5\t$(md5sum "$RTAB" | awk '{print $1}')" | tee -a "$SUMMARY"
echo -e "# meta_md5\t$(md5sum "$META" | awk '{print $1}')" | tee -a "$SUMMARY"
echo -e "# awk_bin\t$(command -v mawk || command -v awk)" | tee -a "$SUMMARY"

# Header IDs
head -n1 "$RTAB" | cut -f2- | tr '\t' '\n' | sort -u > "$TMP_HEADER"
header_n=$(wc -l < "$TMP_HEADER")
(( header_n > 0 )) || { echo "ERROR: RTAB has no sample columns" >&2; exit 2; }

# Build keep set
awk -F'\t' -v KEEP="$KEEP_IDS" -v SPEC="$TMP_SPEC_COUNT" -v HQREM="$TMP_HQ_REMOVED" '
FNR==NR { if (NR>1 && $1!="") tgt[$1]=1; next }                # target species
FILENAME==ARGV[2] { hdr[$1]=1; next }                          # RTAB header IDs
FILENAME==ARGV[3] {
  if (FNR==1) next
  sid=$1; sp=$2; hq=$NF
  if (sid=="" || sp=="") next
  if ((sp in tgt) && (sid in hdr)) {
    spec_candidates++
    if (tolower(hq)=="true") keep[sid]=1; else hq_removed++
  }
  next
}
END {
  for (k in keep) print k > KEEP
  print (spec_candidates+0) > SPEC
  print (hq_removed+0) > HQREM
}
' "$TARGET" "$TMP_HEADER" "$META"

sort -u -o "$KEEP_IDS" "$KEEP_IDS"
post_hq_candidates=$(wc -l < "$KEEP_IDS")
(( post_hq_candidates > 0 )) || { echo "ERROR: 0 samples passed species+HQ filters." >&2; exit 2; }

# Original dimensions
orig_samples=$(head -n1 "$RTAB" | awk -F'\t' '{print NF-1}')
orig_genes=$(( $(wc -l < "$RTAB") - 1 ))

# Match count in header
matched_in_header=$(
  awk -F'\t' -v KEEP="$KEEP_IDS" '
  BEGIN { while ((getline k < KEEP) > 0) keep[k]=1; close(KEEP) }
  NR==1 { m=0; for (i=2;i<=NF;i++) if ($i in keep) m++; print m; exit }
  ' "$RTAB"
)
header_nontarget=$(( orig_samples - matched_in_header ))

# Filtering
AWK_BIN=$(command -v mawk || command -v awk)
"$AWK_BIN" -F'\t' -v OFS='\t' -v KEEP="$KEEP_IDS" -v DROPFILE="$GENE_DROP_FILE" '
BEGIN{ while ((getline k < KEEP) > 0) keep[k]=1; close(KEEP); dropped=0 }
NR==1{
  m=1; cols[m]=1
  for (i=2;i<=NF;i++) if ($i in keep) cols[++m]=i
  if (m==1) { print "ERROR: no matching keep columns" > "/dev/stderr"; exit 2 }
  for (j=1;j<=m;j++) printf "%s%s", $(cols[j]), (j<m?OFS:ORS)
  next
}
{
  present=0
  for (j=2;j<=m;j++){ v=$(cols[j]); if (v!="" && v!="0"){ present=1; break } }
  if (!present){ dropped++; next }
  for (j=1;j<=m;j++) printf "%s%s", $(cols[j]), (j<m?OFS:ORS)
}
END{ print dropped > DROPFILE }
' "$RTAB" > "$OUT"

# Filtered dimensions
filt_samples=$(head -n1 "$OUT" | awk -F'\t' '{print NF-1}')
filt_genes=$(( $(wc -l < "$OUT") - 1 ))
genes_dropped_all_zero=$(cat "$GENE_DROP_FILE")
species_candidates=$(cat "$TMP_SPEC_COUNT")
hq_removed=$(cat "$TMP_HQ_REMOVED")

# Summary (tabular)
{
  echo -e "Metric\tCount"
  echo -e "Original_Rtab_Samples\t$orig_samples"
  echo -e "Original_Rtab_Genes\t$orig_genes"
  echo -e "Target_Species_Listed\t$(($(wc -l < "$TARGET")-1))"
  echo -e "Species_Candidates_in_RtabHeader\t$species_candidates"
  echo -e "Removed_by_high_quality_FALSE\t$hq_removed"
  echo -e "Candidates_after_HQ_TRUE\t$post_hq_candidates"
  echo -e "Matched_in_Rtab_Header\t$matched_in_header"
  echo -e "Header_NonTarget_Samples\t$header_nontarget"
  echo -e "Filtered_Rtab_Samples\t$filt_samples"
  echo -e "Genes_Dropped_AllZero\t$genes_dropped_all_zero"
  echo -e "Filtered_Rtab_Genes\t$filt_genes"
  echo -e "Filtered_Rtab_Path\t$OUT"
} >> "$SUMMARY"

echo "Done. Summary: $SUMMARY"
