#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# Script Name: make_groups_tab.sh
#
# Purpose:
#   Generate a Twilight-compatible groups.tab file.
#   - Extracts sample IDs from the header of gene_presence_absence_filtered.Rtab
#   - Maps each sample ID to its species using File4_QC_characterisation_661K.tsv
#   - Outputs a two-column tab-delimited file: sample_id <TAB> species
#   - Skips any samples that are not found in File4 (with a warning to stderr)
#
# Input:
#   1. gene_presence_absence_filtered.Rtab  (Panaroo output, filtered)
#   2. File4_QC_characterisation_661K.tsv   (metadata; col1=sample_id, col2=species)
#
# Output:
#   groups.tab (in twilight_input directory)
#   - No header row
#   - One line per sample kept
#BSUB -J make_groups_tab
#BSUB -o /data/pam/team230/sm71/scratch/rp2/logs/make_groups_tab.%J.out
#BSUB -e /data/pam/team230/sm71/scratch/rp2/logs/make_groups_tab.%J.err
#BSUB -n 1
#BSUB -M 2000
#BSUB -R "select[mem>2000] rusage[mem=2000]"
# -----------------------------------------------------------------------------

set -euo pipefail
export LC_ALL=C

META_DIR="/data/pam/team230/sm71/scratch/rp2/metadata"
IN_RTAB="/data/pam/team230/sm71/scratch/rp2/pan_genome_processing/gene_presence_absence_filtered.Rtab"
FILE4="$META_DIR/File4_QC_characterisation_661K.tsv"
OUT_DIR="/data/pam/team230/sm71/scratch/rp2/pan_genome_processing"
OUT="$OUT_DIR/groups.tab"

mkdir -p "$OUT_DIR"

AWK_BIN="$(command -v mawk || command -v awk)"

"$AWK_BIN" -F'\t' -v OFS='\t' -v F4="$FILE4" '
BEGIN{
  while ((getline line < F4) > 0) {
    sub(/\r$/, "", line)
    if (line=="") continue
    split(line, a, "\t")
    if (++row_f4==1) continue
    sid=a[1]; sp=a[2]
    if (sid!="") species[sid]=sp
  }
  close(F4)
}
NR==1{
  for (i=2; i<=NF; i++) {
    sid=$i
    if (sid in species) {
      print sid, species[sid]
      kept++
    } else {
      printf("WARN: sample_id %s not found in File4; skipping\n", sid) > "/dev/stderr"
    }
  }
  exit
}
' "$IN_RTAB" | sed 's/\r$//' > "$OUT"

rows=$(wc -l < "$OUT" || true)
echo "groups.tab written: $OUT (${rows} rows)."
