#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C

# -----------------------------------------------------------------------------
# Build per-species query FASTAs of species-specific clusters.
# Inputs:
#   - ssc_with_annotations1.tab (TSV; cols: gene_name, species, gff_file, scaffold_name, clustering_id)
#   - combined_DNA_CDS.fasta    (FASTA; first token after '>' is gene_name)
# Outputs (under OUTDIR):
#   - <Species>_ssc_candidates.fa           (human headers with spaces; all gene names joined by ~~~)
#   - <Species>_ssc_candidates.blastsafe.fa (single-token headers; no spaces, for BLAST parsing)
#   - clusters_species_specific.tsv, clusters_dropped_multi_species.tsv, species_specific_clusters.meta.tsv
# -----------------------------------------------------------------------------

SSC="${SSC:-/data/pam/team230/sm71/scratch/rp2/blast_preprocessing/ssc_with_annotations1.tab}"
FASTA="${FASTA:-/data/pam/team230/sm71/scratch/rp2/panaroo_output/combined_DNA_CDS.fasta}"
OUTDIR="${OUTDIR:-/data/pam/team230/sm71/scratch/rp2/run_blast/queries}"

mkdir -p "$OUTDIR"

LEN="$OUTDIR/gene_len.tsv"
SPECIFIC="$OUTDIR/clusters_species_specific.tsv"
DROPPED="$OUTDIR/clusters_dropped_multi_species.tsv"
META="$OUTDIR/species_specific_clusters.meta.tsv"

[[ -s "$SSC" && -s "$FASTA" ]] || { echo "ERROR: missing $SSC or $FASTA" >&2; exit 1; }

echo "[1/4] Computing per-gene CDS lengths → $LEN"
awk '
  /^>/ { id=$1; sub(/^>/,"",id); curr=id; next }
  { len[curr]+=length($0) }
  END { for (g in len) print g "\t" len[g] }
' "$FASTA" | LC_ALL=C sort -u > "$LEN"

echo "[2/4] Selecting clusters seen in EXACTLY one species → $SPECIFIC"
awk -F'\t' '
NR==1{
  for(i=1;i<=NF;i++){
    n=tolower($i);
    if(n=="clustering_id") cid=i;
    if(n=="species")      sp=i;
  }
  if(!cid||!sp){ print "ERROR: header must include clustering_id and species" > "/dev/stderr"; exit 2 }
  next
}
{
  key=$cid SUBSEP $sp
  if(!(key in seen)){ seen[key]=1; spcount[$cid]++; lastsp[$cid]=$sp }
}
END{
  for(c in spcount){
    if (spcount[c]==1) print c "\t" lastsp[c];
    else               print c "\t" spcount[c] > "'"$DROPPED"'"
  }
}
' "$SSC" | LC_ALL=C sort -t $'\t' -k1,1 > "$SPECIFIC"

kept=$(wc -l < "$SPECIFIC" || true)
dropped=$(wc -l < "$DROPPED" 2>/dev/null || echo 0)
echo "   Kept (species-specific clusters): $kept"
echo "   Dropped (>1 species):             $dropped  → $DROPPED"

echo "[3/4] For kept clusters: choose representative (longest CDS) and aggregate ALL gene names → $META"
awk -F'\t' -v OFS='\t' '
FNR==NR { L[$1]=$2; next }                      # lengths: gene_name -> len
FNR==NR+NR { keep[$1]=1; sp_of[$1]=$2; next }   # kept clusters: cluster_id -> species
NR==1{
  for(i=1;i<=NF;i++){
    n=tolower($i)
    if(n=="gene_name")      gn=i
    if(n=="species")        sp=i
    if(n=="gff_file")       gf=i
    if(n=="scaffold_name")  sc=i
    if(n=="clustering_id")  cid=i
  }
  if(!gn||!sp||!gf||!sc||!cid){
    print "ERROR: SSC must include gene_name,species,gff_file,scaffold_name,clustering_id" > "/dev/stderr"; exit 3
  }
  next
}
($cid in keep){
  cl=$cid; s=$sp; g=$gn; gff=$gf; scf=$sc

  # collect unique gene names per cluster (for header, joined by ~~~)
  k=cl SUBSEP g
  if(!(k in seen_g)){ seen_g[k]=1; genes[cl]=(genes[cl]==""?g:genes[cl]"~~~"g) }

  # pick representative gene by longest CDS (fallback: first seen with len=-1)
  len=(g in L ? L[g]+0 : -1)
  if(!(cl in best_len) || len > best_len[cl]){
    best_len[cl]=len; rep_gene[cl]=g; rep_gff[cl]=gff; rep_scaf[cl]=scf
  }
}
END{
  print "clustering_id","species","rep_gene","rep_gff","rep_scaffold","all_gene_names"
  for (cl in rep_gene){
    print cl, sp_of[cl], rep_gene[cl], rep_gff[cl], rep_scaf[cl], genes[cl]
  }
}
' "$LEN" "$SPECIFIC" "$SSC" | LC_ALL=C sort -t $'\t' -k2,2 -k1,1 > "$META"

echo "[4/4] Extract sequences and write per-species FASTAs under: $OUTDIR"
awk -v META="$META" -v OUTDIR="$OUTDIR" '
BEGIN{
  FS=OFS="\t"
  # read meta rows (skip header)
  while((getline < META)>0){
    if(++r==1) continue
    cl=$1; sp=$2; rg=$3; gf=$4; sc=$5; allg=$6
    info[rg]=cl OFS sp OFS gf OFS sc OFS allg
    reps[rg]=1
    spfile=sp; gsub(/ /,"_",spfile)
    out_human[rg]=OUTDIR "/" spfile "_ssc_candidates.fa"
    out_safe[rg] =OUTDIR "/" spfile "_ssc_candidates.blastsafe.fa"
  }
  close(META)
}
# Stream FASTA once and route each rep gene to its species files
/^>/{
  id=$1; sub(/^>/,"",id)             # first token after '>'
  if(id in reps){
    split(info[id], a, OFS)
    cl=a[1]; sp=a[2]; gf=a[3]; sc=a[4]; genes=a[5]

    # human-readable header with spaces; ALL gene names joined by ~~~
    human = ">" cl " " id " " sp " " gf " " sc " " genes

    # BLAST-safe header (single token; no spaces)
    sp_s=sp; gsub(/ /,"_",sp_s)
    gf_s=gf; gsub(/ /,"_",gf_s)
    sc_s=sc; gsub(/ /,"_",sc_s)
    genes_s=genes; gsub(/ /,"_",genes_s)
    safe = ">cl:" cl "|gene:" id "|sp:" sp_s "|gff:" gf_s "|scaf:" sc_s "|genes:" genes_s

    OH = out_human[id]; OS = out_safe[id]
    print human >> OH
    print safe  >> OS
    emit_h=OH; emit_s=OS; write=1
  } else {
    write=0
  }
  next
}
write{
  print >> emit_h
  print >> emit_s
}
' "$FASTA"

echo "Done."
echo "  Meta table                 : $META"
echo "  Kept clusters (species-only): $SPECIFIC"
echo "  Dropped clusters (>1 sp)   : $DROPPED"
echo "  Per-species FASTAs in      : $OUTDIR"
echo "  NOTE: Use *_blastsafe.fa for BLAST (no spaces in qseqid)."
