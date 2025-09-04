#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# This script filters a Panaroo gene_presence_absence.Rtab file
# It keeps only the samples listed in an "allowed samples" file
# It removes genes that are absent (0/blank) across all kept samples
# It writes out a new filtered Rtab file
# The script is set up to run on an LSF cluster with memory and time limits
# It also performs sanity checks (input files exist, sample IDs match)
# If no sample IDs match, it fails early instead of producing an empty/useless file
# -----------------------------------------------------------------------------

#BSUB -J filter_rtab
#BSUB -o logs/filter_rtab.%J.out
#BSUB -e logs/filter_rtab.%J.err
#BSUB -n 1
#BSUB -M 4000
#BSUB -R "select[mem>4000] rusage[mem=4000]"
#BSUB -W 01:00
set -euo pipefail

RTAB="/data/pam/team230/sm71/scratch/rp2/panaroo_output/gene_presence_absence.Rtab"
ALLOW="/data/pam/team230/sm71/scratch/rp2/twilight_input/samples_step2_allowed.tsv"
OUT="/data/pam/team230/sm71/scratch/rp2/panaroo_output/gene_presence_absence.filtered.Rtab"

[[ -s "$RTAB" ]]  || { echo "ERR: RTAB not found: $RTAB" >&2; exit 1; }
[[ -s "$ALLOW" ]] || { echo "ERR: allow-list not found: $ALLOW" >&2; exit 1; }

awk -v ALLOW="$ALLOW" -F'\t' -v OFS='\t' '
function trim(x){ gsub(/^[[:space:]]+|[[:space:]]+$/, "", x); return x }

BEGIN{
  # Read allowed sample IDs (first column only). Skip header if present.
  while ((getline line < ALLOW) > 0) {
    gsub(/\r$/, "", line)                 # drop CR
    if (line == "") continue
    n = split(line, a, /\t/)              # handle multi-column allow files
    id = trim(a[1])
    gsub(/^"|"$/, "", id)                 # strip quotes if any
    if (id == "" || tolower(id) == "sample" || tolower(id) == "sample_id") continue
    allow[id] = 1
    allow_count++
  }
  close(ALLOW)
}

NR==1{
  # Build keep list: always col 1 (gene), then only allowed sample IDs
  keepN = 0
  keep[++keepN] = 1

  total_samples = 0
  matched = 0

  for (i=2; i<=NF; i++) {
    total_samples++
    h = $i
    if ((h) in allow) {
      keep[++keepN] = i
      matched++
    } else {
      # uncomment to see missing:
      # print "MISS\t" h > "/dev/stderr"
    }
  }

  # Summary to stderr
  print "allow_ids:", allow_count, " header_samples:", total_samples, " matched:", matched > "/dev/stderr"

  # If nothing matched, hard fail so you don’t get a useless file
  if (matched == 0) {
    print "ERROR: 0 header columns matched any allowed IDs. Check naming/first column of allow file." > "/dev/stderr"
    exit 2
  }

  # Print filtered header
  for (k=1; k<=keepN; k++) {
    printf "%s%s", $(keep[k]), (k<keepN?OFS:ORS)
  }
  next
}

# Data rows: drop genes that are all 0/blank across kept sample columns
{
  present = 0
  for (k=2; k<=keepN; k++) {
    v = $(keep[k])
    if (v != "" && v != "0") { present = 1; break }
  }
  if (present) {
    for (k=1; k<=keepN; k++) {
      printf "%s%s", $(keep[k]), (k<keepN?OFS:ORS)
    }
  }
}
' "$RTAB" > "$OUT"

echo "done ✅ wrote: $OUT"
