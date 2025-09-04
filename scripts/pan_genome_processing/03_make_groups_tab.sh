# -----------------------------------------------------------------------------
# Creates a Twilight-compatible groups.tab from samples_twilight.tsv by keeping
# only the first two columns (sample_id, species) and converting TSV → TAB.
# Removes any Windows carriage returns and drops the header row.
# -----------------------------------------------------------------------------


IN="/data/pam/team230/sm71/scratch/rp2/twilight_input/samples_twilight.tsv"
OUT="/data/pam/team230/sm71/scratch/rp2/twilight_input/groups.tab"

# Rebuild the 2 columns cleanly and save as .tab
awk -F'\t' 'NR==1{print "sample_id\tspecies"; next} {print $1 "\t" $2}' "$IN" \
  | sed 's/\r$//' > "$OUT"
#remove headers
sed -i '1d' /data/pam/team230/sm71/scratch/rp2/twilight_input/groups.tab




