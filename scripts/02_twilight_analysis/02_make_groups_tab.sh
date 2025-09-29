#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# make_groups_tab.sh  — generate Twilight-compatible groups.tab     (v1.0.0)
#
# Purpose:
#   - Extract sample IDs from the header of a filtered RTAB
#   - Map each sample ID to its species using File4_QC_characterisation_661K.tsv
#   - Emit a two-column, tab-delimited file:  sample_id<TAB>species  (no header)
#
# Notes:
#   - Metadata file (File4_QC_characterisation_661K.tsv) is from Blackwell et al., 2021 (Figshare).
#   - This script is original work (not adapted from other code).
# -----------------------------------------------------------------------------

set -euo pipefail
export LC_ALL=C

usage() {
  cat <<'EOF'
Usage:
  make_groups_tab.sh \
    --rtab path/to/gene_presence_absence_filtered.Rtab \
    --metadata path/to/File4_QC_characterisation_661K.tsv \
    --outdir path/to/output \
    [--outfile groups.tab] [--sort-by id|species]

Required:
  --rtab        Filtered RTAB (first column 'Gene', remaining columns are sample IDs)
  --metadata    File4_QC_characterisation_661K.tsv (col1=sample_id, col2=species)
  --outdir      Output directory

Optional:
  --outfile     Output filename (default: groups.tab)
  --sort-by     Sort rows by 'id' or 'species' (default: id)
  --help        Show this help
EOF
}

# -------- parse args --------
RTAB="" ; META="" ; OUTDIR="" ; OUTFILE="groups.tab" ; SORT_BY="id"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --rtab) RTAB="$2"; shift 2 ;;
    --metadata) META="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --outfile) OUTFILE="$2"; shift 2 ;;
    --sort-by) SORT_BY="$2"; shift 2 ;;
    --help|-h) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage; exit 2 ;;
  endesac || true
  esac
done

# -------- validate --------
[[ -n "$RTAB" && -f "$RTAB" ]] || { echo "ERROR: --rtab not found: $RTAB" >&2; exit 2; }
[[ -n "$META" && -f "$META" ]] || { echo "ERROR: --metadata not found: $META" >&2; exit 2; }
[[ -n "$OUTDIR" ]] || { echo "ERROR: --outdir is required" >&2; exit 2; }
mkdir -p "$OUTDIR"

OUT="$OUTDIR/$OUTFILE"
TMP_IDS="$OUTDIR/.tmp_header_ids.$$"
TMP_RAW="$OUTDIR/.tmp_groups_raw.$$"
trap 'rm -f "$TMP_IDS" "$TMP_RAW"' EXIT

AWK_BIN="$(command -v mawk || command -v awk)"

# -------- provenance log --------
{
  echo "# make_groups_tab.sh v1.0.0"
  echo "# date_utc: $(date -u +"%Y-%m-%dT%H:%M:%SZ")"
  command -v md5sum >/dev/null 2>&1 && {
    echo -e "# rtab_md5\t$(md5sum "$RTAB" | awk "{print \$1}")"
    echo -e "# meta_md5\t$(md5sum "$META" | awk "{print \$1}")"
  } || true
  echo -e "# awk_bin\t$AWK_BIN"
} > "$OUTDIR/groups_tab_summary.txt"

# -------- extract header sample IDs from RTAB --------
head -n1 "$RTAB" | cut -f2- | tr '\t' '\n' | sed 's/\r$//' | sort -u > "$TMP_IDS"
header_n=$(wc -l < "$TMP_IDS")
(( header_n > 0 )) || { echo "ERROR: RTAB header contains no sample IDs" >&2; exit 2; }

# -------- build groups table (sample_id<TAB>species) --------
"$AWK_BIN
" -F'\t' -v OFS='\t' -v F4="$META" '
BEGIN{
  # load metadata: col1=sample_id, col2=species (skip header)
  while ((getline line < F4) > 0) {
    sub(/\r$/, "", line)
    if (line=="") continue
    n = split(line, a, "\t")
    if (++row_f4==1) continue
    sid=a[1]; sp=a[2]
    if (sid!="") species[sid]=sp
  }
  close(F4)
}
# read one sample_id per line from header list and emit mapping if present
{
  sid=$0
  if (sid in species) {
    print sid, species[sid]
    kept++
    counts[species[sid]]++
  } else {
    # warn to stderr for transparency
    printf("WARN: sample_id %s not found in metadata; skipping\n", sid) > "/dev/stderr"
    skipped++
  }
}
END{
  # print a small summary to stdout (caller can ignore)
  printf("kept=%d\tskipped=%d\n", (kept+0), (skipped+0)) > "/dev/stdout"
  for (sp in counts) printf("species_count\t%s\t%d\n", sp, counts[sp]) > "/dev/stdout"
}
' "$TMP_IDS" > "$TMP_RAW"

# -------- optional sorting --------
case "$SORT_BY" in
  id)      sort -t$'\t' -k1,1V "$TMP_RAW" > "$OUT" ;;
  species) sort -t$'\t' -k2,2V -k1,1V "$TMP_RAW" > "$OUT" ;;
  *) echo "WARN: unknown --sort-by '$SORT_BY' (using 'id')" >&2; sort -t$'\t' -k1,1V "$TMP_RAW" > "$OUT" ;;
esac

rows=$(wc -l < "$OUT" || true)
echo "groups.tab written: $OUT (${rows} rows)."

# -------- append counts per species to summary --------
{
  echo -e "Output\t$OUT"
  echo -e "Rows\t$rows"
  echo -e "SortBy\t$SORT_BY"
  echo -e "Species\tCount"
  cut -f2 "$OUT" | sort | uniq -c | awk '{print $2 "\t" $1}'
} >> "$OUTDIR/groups_tab_summary.txt"
