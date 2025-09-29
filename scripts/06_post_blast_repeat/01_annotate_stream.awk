#!/usr/bin/env bash
# process_blast.sh
# Simple end-to-end script for splitting, sorting, and joining BLAST results

cd /data/pam/team230/sm71/scratch/rp2/blast_run/ref_db/results || exit 1

# STEP 1: Make a copy of the blast results
cp blast_results.tsv blast_results2.tsv

# STEP 2: Create parts folder
mkdir -p parts

# STEP 3: Split into 50 parts
split -d -n l/50 -a 4 --additional-suffix=.tsv \
  blast_results2.tsv parts/part_

# STEP 4: Sort and join each part with query truth (gene name + species)
for file in part_[0-9][0-9][0-9][0-9].tsv; do
  sorted="${file%.tsv}.sorted.tsv"
  joined="${file%.tsv}_with_species.tsv"

  # Skip if final joined file already exists
  if [[ -f "$joined" ]]; then
    echo "Skipping $file (already processed)"
    continue
  fi

  echo "Processing $file..."
  # Sort by qseqid (column 1)
  sort -t $'\t' -k1,1 "$file" > "$sorted"

  # Join with truth queries file
  join -t $'\t' -1 1 -2 1 -a 2 -e '' -o auto truth_queries.sorted.tsv  "$sorted" > "$joined"
done

# STEP 5: Join with reference truth to add ref species
for file in part_*_with_species.tsv; do
  final="${file%.tsv}_with_ref_species.tsv"

  # Skip if final ref-species file already exists
  if [[ -f "$final" ]]; then
    echo "⏩ Skipping $file (ref species already added)"
    continue
  fi

  echo "⚡ tagging ref species -> $file"
  awk -F'\t' -v OFS='\t' '
    FNR==NR {ref[$1]=$2; next}   # build lookup from truth_refs: ref[SAMD]=species
    {print $0, ref[$4]}          # append ref species by sample id in col 4
  ' truth_refs.sorted.tsv "$file" > "$final"
done

echo "✅ All done!"

# STEP 6: Merge all part_*_with_species_with_ref_species.tsv files into one:
cat part_*_with_species_with_ref_species.tsv > blast_results_pre.tsv

#check:
cat part_*_with_species_with_ref_species.tsv | wc -l
wc -l blast_results_pre.tsv

#STEP 7: add headers:
#Create a temp file with the header
echo -e "qseqid\tgene_name\ttrue_species\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tqcovhsp\tqcovs\tqcovus\tref_species" > blast_results_pre.tmp

#Append the original content
cat blast_results_pre.tsv >> blast_results_pre.tmp

#Replace the old file with the new one
mv blast_results_pre.tmp blast_results_pre.tsv

# Step 8: Produce blast results output with threshold filters:
awk -F'\t' -v OFS='\t' 'NR==1 || ($5 >= 80 && $17 >= 90)' blast_results_pre.tsv > blast_results_pre.filtered && mv blast_results_pre.filtered blast_results_pre.tsv
