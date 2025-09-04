#!/bin/bash
# -----------------------------------------------------------------------------
# This script compares gene counts per species after Twilight filtering vs after
# QC filtering. It:
# - Reads the species_specific_core.tab file from BLAST output.
# - Counts all genes assigned to each species (Twilight results).
# - Identifies genes unique to exactly one species (QC definition).
# - Builds a comparison table per species with both counts, plus overall totals.
# - Saves the summary table to summary/gene_count_twilight_vs_qc.txt.
# -----------------------------------------------------------------------------
#BSUB -q normal
#BSUB -J twilight_qc_compare
#BSUB -n 1
#BSUB -M 8000
#BSUB -R "select[mem>8000] rusage[mem=8000]"
#BSUB -W 00:10
#BSUB -o %J.out
#BSUB -e %J.err
set -euo pipefail
export LC_ALL=C

SSC_TAB="/data/pam/team230/sm71/scratch/rp2/blast/species_specific_core.tab"
OUTDIR="/data/pam/team230/sm71/scratch/rp2/summary"
mkdir -p "$OUTDIR"
OUT="$OUTDIR/gene_count_twilight_vs_qc.txt"

# Build final table in one pass (no saved intermediates)
{
  echo -e "### Twilight vs QC gene counts per species"
  echo -e "Species\tAfter_Twilight\tAfter_QC"

  # One awk pass:
  # - raw[s] counts all SSC rows per species (After Twilight)
  # - gspec[g] tracks which species a gene_name belongs to; if it shows up in >1 species, mark MULTI
  # - uniq[s] counts genes that occur in exactly ONE species (After QC definition)
  awk -F'\t' 'NR>1{
      g=$1; s=$2
      raw[s]++
      if (!(g in gspec)) {
          gspec[g]=s
      } else if (gspec[g]!=s && gspec[g]!="MULTI") {
          gspec[g]="MULTI"
      }
  }
  END{
      # accumulate unique-only genes per species
      for (g in gspec) if (gspec[g]!="MULTI") uniq[gspec[g]]++

      # print per-species rows, sorted by species
      # collect species names first to sort deterministically
      for (s in raw) sp[++n]=s
      PROCINFO["sorted_in"]="@ind_str_asc"
      for (i=1;i<=n;i++){
          s=sp[i]
          at = raw[s]+0
          aq = (s in uniq ? uniq[s] : 0)
          print s "\t" at "\t" aq
          tot_at += at
          tot_aq += aq
      }

      # blank line + totals
      print ""
      print "TOTAL (all species):"
      print "After Twilight\t" tot_at
      print "After QC\t"       tot_aq
  }' "$SSC_TAB"
} > "$OUT"

echo "[DONE] Comparison written to $OUT"
