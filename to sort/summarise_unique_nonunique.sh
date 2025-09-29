#!/usr/bin/env bash
# summarize_unique_nonunique_after_removal.sh  (LINE-BASED)
#
# Input:
#   /data/pam/team230/sm71/scratch/rp2/blast_preprocessing/species_specific_core.tab
#
# Output:
#   /data/pam/team230/sm71/scratch/rp2/blast_preprocessing/species_gene_summary.txt
#
# Behavior:
#   - Classify each gene as unique (appears in exactly 1 species) or non-unique (>=2 species).
#   - Per-species tallies are LINE-BASED so they match expanded counts:
#       * For each row (gene,species), if that gene is unique → unique_lines[species]++,
#         else nonunique_lines[species]++.
#   - "total_genes" per species = unique_lines + nonunique_lines.

set -euo pipefail
export LC_ALL=C

INFILE="/data/pam/team230/sm71/scratch/rp2/blast_preprocessing/species_specific_core.tab"
OUTFILE="/data/pam/team230/sm71/scratch/rp2/blast_preprocessing/species_gene_summary.txt"

awk -F'\t' -v OUTFILE="$OUTFILE" '
NR==1{ next }  # skip header

{
  g=$1; s=$2

  # Track distinct species per gene (for uniqueness classification)
  gs_set[g, s]=1

  # Track line counts per (gene,species) so we can tally line-based later
  pair_count[g, s]++

  # Remember species list
  species_set[s]=1
  overall_total_lines++
}

END{
  # 1) Decide uniqueness per gene by counting distinct species it occurs in
  for (k in gs_set) {
    split(k, a, SUBSEP); g=a[1]; s=a[2]
    species_count[g]++
  }
  for (g in species_count) {
    is_unique[g] = (species_count[g] == 1) ? 1 : 0
  }

  # 2) Line-based tallies per species
  for (k in pair_count) {
    split(k, a, SUBSEP); g=a[1]; s=a[2]
    c = pair_count[k]
    if (is_unique[g]) {
      unique_lines[s] += c
      overall_unique_lines += c
    } else {
      nonunique_lines[s] += c
      overall_nonunique_lines += c
    }
  }

  # 3) Write output
  out = OUTFILE
  print "Overall totals (line-based):"                 >  out
  print "  Total genes (lines): " overall_total_lines  >> out
  print "  Unique genes (kept): " overall_unique_lines >> out
  print "  Non-unique genes (removed): " overall_nonunique_lines >> out
  print "" >> out

  print "Per-species tallies (line-based):"           >> out
  print "species\ttotal_genes\tunique_lines\tnon_unique_lines" >> out
  for (s in species_set) {
    u  = (unique_lines[s]+0)
    nu = (nonunique_lines[s]+0)
    tot = u + nu
    print s "\t" tot "\t" u "\t" nu                    >> out
  }
  close(out)

  # sanity check to stderr
  if ((overall_unique_lines + overall_nonunique_lines) != overall_total_lines) {
    print "[WARN] Totals mismatch: unique+non-unique != total lines" > "/dev/stderr"
  }
  print "[DONE] Wrote summary to " out > "/dev/stderr"
}' "$INFILE"
