#To count all RefSeq complete/chromosome assemblies per genus:
for G in Staphylococcus Enterococcus Streptococcus; do
  echo -n "$G: "
  datasets summary genome taxon "$G" \
    --assembly-source RefSeq \
    --assembly-level complete,chromosome \
    --as-json-lines \
  | wc -l
done


# Get all accessions into a file
for G in Staphylococcus Enterococcus Streptococcus; do
  datasets summary genome taxon "$G" \
    --assembly-source RefSeq \
    --assembly-level complete,chromosome \
    --report ids_only \
    --as-json-lines
done | jq -r '.accession' > all_refseq_accessions.txt

split -l 500 all_refseq_accessions.txt acc_batch_

for batch in acc_batch_*; do
  datasets download genome accession --inputfile "$batch" \
    --include genome,seq-report \
    --dehydrated \
    --filename "01_zips/${batch}.zip"
done

cd /data/pam/team230/sm71/scratch/rp2/run_blast/ref_db

# sanity: where are the batch files?
ls 01_zips/acc_batch_* | head

# download dehydrated zips for each batch
for batch in 01_zips/acc_batch_*; do
  out="01_zips/$(basename "$batch").zip"
  datasets download genome accession \
    --inputfile "$batch" \
    --include genome,seq-report \
    --dehydrated \
    --filename "$out"
done

# check you actually have .zip files now
ls -lh 01_zips | head

mkdir -p 02_raw
for z in 01_zips/*.zip; do
  d="02_raw/$(basename "${z%.zip}")"
  mkdir -p "$d"
  unzip -q "$z" -d "$d"
  datasets rehydrate --directory "$d"
done

mkdir -p 03_fa
cd 02_raw
find . -type f -name "*_genomic.fna" | sort | \
while read -r fa; do
  acc=$(basename "$(dirname "$fa")")
  awk -v p="${acc}|" 'BEGIN{OFS=""} /^>/{sub(/^>/,">" p); print; next} {print}' "$fa"
done > ../03_fa/refdb.fa
cd ..

cd 03_fa
makeblastdb -in refdb.fa -dbtype nucl -parse_seqids -title refdb_ncbi
minimap2 -d refdb.mmi refdb.fa
samtools faidx refdb.fa