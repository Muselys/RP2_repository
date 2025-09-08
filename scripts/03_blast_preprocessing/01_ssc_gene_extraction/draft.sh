#/data/pam/team230/sm71/scratch/rp2/
#    ├── blast_preprocessing/
#   │   ├── setup/
#  │   │   ├── fasta/                # FASTA + .fai (index)
#    │   │   ├── tables/               # TSVs: lists + meta
#   │   │   ├── ids/                  # 18 per-species cluster-id lists (*.ids)
#  │   │   └── logs/                 # logs for the prep steps
# │   └── README.md                 # what lives here / how to regenerate
#  └── run_blast/
#       └── queries/                  # per-species *_ssc_candidates*.fa (outputs)
 
#!/usr/bin/env bash
#BSUB -q normal
#BSUB -J make_lists
#BSUB -n 1
#BSUB -M 16000
#BSUB -R "rusage[mem=16000]"
#BSUB -o /data/pam/team230/sm71/scratch/rp2/logs/make_lists.%J.out
#BSUB -e /data/pam/team230/sm71/scratch/rp2/logs/make_lists.%J.err
set -euo pipefail
export LC_ALL=C

BASE=/data/pam/team230/sm71/scratch/rp2
BP="$BASE/blast_preprocessing/setup"
PAN="$BASE/panaroo_output"
TABDIR="$BP/tables"
IDDIR="$BP/ids"
SSC="$TABDIR/ssc_with_annotations1.tab"

mkdir -p "$BP"/{fasta,tables,ids,logs}
mkdir -p "$TABDIR" "$IDDIR"

# Put FASTA here so the .fai index sits alongside (simplest + avoids permission quirks)
ln -s "$PAN/combined_DNA_CDS.fasta" "$BP/fasta/combined_DNA_CDS.fasta"

module load samtools-1.19.2
module load python-3.11.6/perl-5.38.0
samtools faidx "$BP/fasta/combined_DNA_CDS.fasta"   # creates combined_DNA_CDS.fasta.fai

# Bring in the joined annotations (symlink or copy—your call)
ln -s "$BASE/blast_preprocessing/ssc_with_annotations1.tab" \
      "$BP/tables/ssc_with_annotations1.tab"




# Prefer gawk (better for large inputs)
AWK="$(command -v gawk || command -v awk)"

[ -s "$SSC" ] || { echo "ERROR: $SSC not found or empty." >&2; exit 1; }

echo "[1/4] Find clusters present in EXACTLY one species"
"$AWK" -F'\t' '
NR==1{
  for(i=1;i<=NF;i++){
    n=tolower($i)
    if(n=="clustering_id") cid=i
    if(n=="species")      sp=i
  }
  if(!cid||!sp){ print "ERROR: need clustering_id and species columns" > "/dev/stderr"; exit 2 }
  next
}
{
  k = $cid SUBSEP $sp
  if(!(k in seen)){ seen[k]=1; spcount[$cid]++; lastsp[$cid]=$sp }
}
END{
  for(c in spcount){
    if (spcount[c]==1) print c "\t" lastsp[c];
    else               print c "\t" spcount[c] > "'"$TABDIR/clusters_dropped_multi_species.tsv"'"
  }
}
' "$SSC" | sort -t $'\t' -k1,1 > "$TABDIR/clusters_species_specific.tsv"

kept=$(wc -l < "$TABDIR/clusters_species_specific.tsv" || true)
drop=$(wc -l < "$TABDIR/clusters_dropped_multi_species.tsv" 2>/dev/null || echo 0)
echo "   Kept species-specific clusters : $kept"
echo "   Dropped multi-species clusters : $drop  → $TABDIR/clusters_dropped_multi_species.tsv"

echo "[2/4] Build META (cluster → species, gff, scaffold, ALL gene names joined with ~~~)"
"$AWK" -F'\t' -v OFS='\t' '
FNR==NR { keep[$1]=$2; next }    # cluster_id -> species (kept only)
NR==1{
  for(i=1;i<=NF;i++){
    n=tolower($i)
    if(n=="gene_name")      gn=i
    if(n=="species")        sp=i
    if(n=="gff_file")       gf=i
    if(n=="scaffold_name")  sc=i
    if(n=="clustering_id")  cid=i
  }
  if(!gn||!sp||!gf||!sc||!cid){ print "ERROR: missing required columns" > "/dev/stderr"; exit 3 }
  next
}
($cid in keep){
  cl=$cid
  # collect unique gene names per cluster (for header; join with ~~~)
  g=$gn
  key=cl SUBSEP g
  if(!(key in seen_g)){ seen_g[key]=1; genes[cl]=(genes[cl]==""?g:genes[cl]"~~~"g) }
  # record first seen gff/scaffold as example
  if(!(cl in have)){ gff[cl]=$gf; scf[cl]=$sc; have[cl]=1 }
}
END{
  print "clustering_id","species","gff_file","scaffold_name","all_gene_names"
  for (cl in keep) print cl, keep[cl], gff[cl], scf[cl], genes[cl]
}
' "$TABDIR/clusters_species_specific.tsv" "$SSC" \
| sort -t $'\t' -k2,2 -k1,1 > "$TABDIR/species_specific_clusters.meta.tsv"

echo "[3/4] Species list & per-species ID files"
cut -f2 "$TABDIR/clusters_species_specific.tsv" | sort -u > "$TABDIR/species_list.txt"
while read -r S; do
  [[ -z "$S" ]] && continue
  SF=${S// /_}
  awk -F'\t' -v s="$S" '$2==s{print $1}' "$TABDIR/clusters_species_specific.tsv" > "$IDDIR/${SF}.ids"
done < "$TABDIR/species_list.txt"

echo "[4/4] Summary"
ls -lh "$TABDIR/clusters_species_specific.tsv" \
       "$TABDIR/clusters_dropped_multi_species.tsv" \
       "$TABDIR/species_specific_clusters.meta.tsv" \
       "$TABDIR/species_list.txt" || true
echo "Per-species ID lists in: $IDDIR"
