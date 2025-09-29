#!/usr/bin/env bash
#BSUB -q hugemem
#BSUB -J make_queries
#BSUB -n 1
#BSUB -M 400000
#BSUB -R "select[mem>400000] rusage[mem=400000]"
#BSUB -o /data/pam/team230/sm71/scratch/rp2/blast_preprocessing/setup/logs/make_queries_pool.%J.out
#BSUB -e /data/pam/team230/sm71/scratch/rp2/blast_preprocessing/setup/logs/make_queries_pool.%J.err
set -euo pipefail
export LC_ALL=C

# -----------------------------------------------------------------------------
# Build per-species query FASTAs from species-specific clusters (cluster-ID keyed FASTA).
#
# Inputs:
#   - blast_preprocessing/ssc_with_annotations1.tab
#       (TSV; cols: gene_name, species, gff_file, scaffold_name, clustering_id)
#   - panaroo_output/combined_DNA_CDS.fasta
#       (FASTA; headers are CLUSTER IDs like >0_0_0)
#
# Outputs:
#   - run_blast/queries/clusters_species_specific.tsv
#   - run_blast/queries/clusters_dropped_multi_species.tsv
#   - run_blast/queries/species_specific_clusters.meta.tsv
#   - run_blast/queries/<Species>_ssc_candidates.fa
#       (human headers: >CLUSTER_ID GENE1~~~GENE2~~~ SPECIES GFF SCAFFOLD)
#   - run_blast/queries/<Species>_ssc_candidates.blastsafe.fa
#       (BLAST-safe headers: >cl:ID|genes:...|sp:...|gff:...|scaf:...)
#
# Optional:
#   export SPECIES_FILTER="Enterococcus faecalis"   # restrict outputs to one species
# -----------------------------------------------------------------------------

BASE="/data/pam/team230/sm71/scratch/rp2"
PAN="$BASE/panaroo_output"
SETUP="$BASE/blast_preprocessing/setup"          # where we place links (optional)
TABDIR="$SETUP/tables"
FASTA_LINK="$SETUP/fasta/combined_DNA_CDS.fasta" # symlink to Panaroo FASTA
SSC_SRC="$BASE/blast_preprocessing/ssc_with_annotations1.tab"
OUTDIR="$BASE/run_blast/queries"
SPECIES_FILTER="${SPECIES_FILTER:-}"             # empty = all species

SPECIFIC="$OUTDIR/clusters_species_specific.tsv"
DROPPED="$OUTDIR/clusters_dropped_multi_species.tsv"
META="$OUTDIR/species_specific_clusters.meta.tsv"

mkdir -p "$SETUP"/{fasta,tables,logs} "$OUTDIR"

# --- Modules (robust in non-interactive shells)
type module &>/dev/null || source /etc/profile.d/modules.sh 2>/dev/null || true
module load samtools-1.19.2 || module load samtools/1.19.2 || true
module load python-3.11.6/perl-5.38.0 || true
command -v samtools >/dev/null || { echo "ERROR: samtools not on PATH"; exit 127; }

# --- Inputs: link FASTA & TSV (idempotent) -----------------------------------
ln -sfn "$PAN/combined_DNA_CDS.fasta" "$FASTA_LINK"
ln -sfn "$SSC_SRC" "$TABDIR/ssc_with_annotations1.tab"

# Index FASTA if needed
[ -s "$FASTA_LINK.fai" ] || samtools faidx "$FASTA_LINK"

# Sanitize TSV (strip CRs just in case)
SSC="$TABDIR/ssc_with_annotations1.tab"
SSC_SAN="$OUTDIR/.ssc_sanitized.$$"
sed 's/\r$//' "$SSC" > "$SSC_SAN"
trap 'rm -f "$SSC_SAN"' EXIT

# Prefer gawk if present
AWK="$(command -v gawk || command -v awk)"

# --- [1/2] Build META directly from SSC (no pre-filtering) --------------------
echo "[1/2] Building meta (all clusters, per species) → $META"
"$AWK" -F'\t' -v OFS='\t' '
NR==1{
  for(i=1;i<=NF;i++){
    h=$i; gsub(/\r/,"",h); gsub(/^[[:space:]]+|[[:space:]]+$/,"",h); n=tolower(h)
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
{
  cl=$cid; s=$sp; gname=$gn; gff=$gf; scf=$sc
  # Trim
  gsub(/^[[:space:]]+|[[:space:]]+$/, "", cl)
  gsub(/^[[:space:]]+|[[:space:]]+$/, "", s)
  gsub(/^[[:space:]]+|[[:space:]]+$/, "", gname)
  gsub(/^[[:space:]]+|[[:space:]]+$/, "", gff)
  gsub(/^[[:space:]]+|[[:space:]]+$/, "", scf)

  # Optional species restriction
  if ("'"$SPECIES_FILTER"'" != "" && s != "'"$SPECIES_FILTER"'") next

  key = cl SUBSEP s
  keys[key]=1

  # accumulate unique gene names per (cl,species)
  if (gname != ""){
    u=key SUBSEP gname
    if(!(u in seen_g)){
      seen_g[u]=1
      genes[key]=(genes[key]==""?gname:genes[key]"~~~"gname)
    }
  }

  # prefer most recent non-empty gff/scaffold
  if (gff != "") gff_by[key]=gff
  if (scf != "") scf_by[key]=scf
}
END{
  print "clustering_id","species","gff_file","scaffold_name","all_gene_names"
  for (k in keys){
    split(k,a,SUBSEP); cl=a[1]; s=a[2]
    g=(k in genes   && genes[k]   != "" ? genes[k]   : "NA")
    f=(k in gff_by  && gff_by[k]  != "" ? gff_by[k]  : "NA")
    c=(k in scf_by  && scf_by[k]  != "" ? scf_by[k]  : "NA")
    print cl, s, f, c, g
  }
}
' "$SSC_SAN" \
| sort -t $'\t' -k2,2 -k1,1 > "$META"


# --- [2/2] Extract sequences by cluster-id; write per-species FASTAs ----------
echo "[2/2] Extract sequences by cluster-id from FASTA → $OUTDIR"

# Clean previous outputs for targeted species/all species to avoid appends
if [ -n "${SPECIES_FILTER}" ]; then
  sp_und="${SPECIES_FILTER// /_}"
  rm -f "$OUTDIR/${sp_und}_ssc_candidates.fa" "$OUTDIR/${sp_und}_ssc_candidates.blastsafe.fa" 2>/dev/null || true
else
  rm -f "$OUTDIR"/*_ssc_candidates.fa "$OUTDIR"/*_ssc_candidates.blastsafe.fa 2>/dev/null || true
fi

"$AWK" -v META="$META" -v OUTDIR="$OUTDIR" '
BEGIN{
  FS=OFS="\t"
  SUBSEP_char = SUBSEP

  # Load META: allow multiple species per cluster
  while((getline < META)>0){
    if(++r==1) continue
    cl=$1; sp=$2; gf=$3; sc=$4; genes=$5

    # species filter handled upstream; still guard here
    if ("'"$SPECIES_FILTER"'" != "" && sp != "'"$SPECIES_FILTER"'") continue

    key = cl SUBSEP sp
    gf_of[key]=gf; sc_of[key]=sc; genes_of[key]=genes

    # per-species output paths
    spkey=sp; gsub(/ /,"_",spkey)
    if(!(sp in human_path)){ human_path[sp]=OUTDIR "/" spkey "_ssc_candidates.fa" }
    if(!(sp in safe_path )){ safe_path [sp]=OUTDIR "/" spkey "_ssc_candidates.blastsafe.fa" }

    # register mapping cluster -> list of species (SUBSEP-joined)
    cl_to_sp[cl] = (cl in cl_to_sp ? cl_to_sp[cl] SUBSEP sp : sp)
  }
  close(META)
}

# FASTA stream; headers are cluster IDs
/^>/{
  id=$1; sub(/^>/,"",id)
  if (id in cl_to_sp){
    # Build list of species for this cluster
    n = split(cl_to_sp[id], splist, SUBSEP)

    # Collect file lists for this cluster so we can append sequence lines
    filelist_h=""; filelist_s=""

    for (i=1; i<=n; i++){
      sp = splist[i]
      key = id SUBSEP sp
      gf = gf_of[key]; sc = sc_of[key]; genes = genes_of[key]
      OH = human_path[sp]; OS = safe_path[sp]

      # human header: >CLUSTER_ID GENE1~~~GENE2~~~ SPECIES GFF SCAFFOLD
      print ">" id " " genes " " sp " " gf " " sc >> OH

      # BLAST-safe header
      sps=sp; gsub(/ /,"_",sps); gfs=gf; gsub(/ /,"_",gfs)
      scs=sc; gsub(/ /,"_",scs); gens=genes; gsub(/ /,"_",gens)
      print ">cl:" id "|genes:" gens "|sp:" sps "|gff:" gfs "|scaf:" scs >> OS

      filelist_h = (filelist_h=="" ? OH : filelist_h SUBSEP OSUB filelist_h) # placeholder
      filelist_h = (filelist_h=="" ? OH : filelist_h SUBSEP OH)
      filelist_s = (filelist_s=="" ? OS : filelist_s SUBSEP OS)
    }

    H_paths[id]=filelist_h
    S_paths[id]=filelist_s
    cur_id=id; emit=1
  } else {
    emit=0
  }
  next
}

emit{
  nh = split(H_paths[cur_id], hp, SUBSEP)
  ns = split(S_paths[cur_id], sp, SUBSEP)
  for (i=1; i<=nh; i++) print >> hp[i]
  for (i=1; i<=ns; i++) print >> sp[i]
}
' "$FASTA_LINK"


echo "Done."
echo "  Kept clusters         : $SPECIFIC"
echo "  Dropped clusters      : $DROPPED"
echo "  Meta                  : $META"
echo "  Per-species FASTAs in : $OUTDIR"
echo "  NOTE: Use *_blastsafe.fa for BLAST (no spaces in qseqid)."
