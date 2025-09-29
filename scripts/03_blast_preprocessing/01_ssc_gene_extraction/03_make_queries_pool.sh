#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# 02_make_queries_pool_all.sh — build per-species SSC query FASTAs (all species)  v1.0.0
# Author: Sahra Musse | License: MIT
#
# Inputs (required)
#   --ssc   PATH   TSV with columns: gene_name, species, gff_file, scaffold_name, clustering_id
#   --fasta PATH   Panaroo combined_DNA_CDS.fasta (headers are cluster IDs like >0_0_0)
#   --outdir DIR   Output directory (created if absent)
#
# Outputs (in <outdir>)
#   clusters_species_specific.tsv
#   clusters_dropped_multi_species.tsv          # empty placeholder (reserved)
#   species_specific_clusters.meta.tsv          # cluster_id, species, gff, scaffold, all_gene_names
#   <Species>_ssc_candidates.fa                 # human headers
#   <Species>_ssc_candidates.blastsafe.fa       # BLAST-safe headers (no spaces)
#
# Requires: bash, awk/gawk, samtools
# Submit on LSF:  bsub < scripts/03_blast_preprocessing/02_make_queries_pool_all.sh
# -----------------------------------------------------------------------------

#BSUB -q hugemem
#BSUB -J make_queries_all
#BSUB -n 1
#BSUB -M 400000
#BSUB -R "select[mem>400000] rusage[mem=400000]"
#BSUB -o logs/make_queries_pool_all.%J.out
#BSUB -e logs/make_queries_pool_all.%J.err

set -euo pipefail
export LC_ALL=C

usage() {
  cat <<'EOF'
Usage:
  02_make_queries_pool_all.sh \
    --ssc blast_preprocessing/ssc_with_annotations1.tab \
    --fasta panaroo_output/combined_DNA_CDS.fasta \
    --outdir run_blast/queries
EOF
}

# ---- parse CLI ---------------------------------------------------------------
SSC="" ; FASTA="" ; OUTDIR=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --ssc)    SSC="$2"; shift 2;;
    --fasta)  FASTA="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    --help|-h) usage; exit 0;;
    *) echo "ERROR: unknown arg $1" >&2; usage; exit 2;;
  esac
done

[[ -n "$SSC"   && -f "$SSC"   ]] || { echo "ERROR: --ssc file not found: $SSC" >&2; exit 2; }
[[ -n "$FASTA" && -f "$FASTA" ]] || { echo "ERROR: --fasta file not found: $FASTA" >&2; exit 2; }
[[ -n "$OUTDIR" ]] || { echo "ERROR: --outdir is required" >&2; exit 2; }

mkdir -p "$OUTDIR" "$OUTDIR/logs"

AWK="$(command -v gawk || command -v awk)"
command -v samtools >/dev/null || { echo "ERROR: samtools not found on PATH" >&2; exit 127; }

# ---- provenance --------------------------------------------------------------
MANIFEST="$OUTDIR/make_queries_manifest.txt"
{
  echo "# 02_make_queries_pool_all.sh v1.0.0"
  echo "# date_utc: $(date -u +"%Y-%m-%dT%H:%M:%SZ")"
  command -v md5sum >/dev/null 2>&1 && {
    echo -e "md5_ssc\t$(md5sum "$SSC" | awk '{print $1}')"
    echo -e "md5_fasta\t$(md5sum "$FASTA" | awk '{print $1}')"
  }
  echo -e "ssc\t$SSC"
  echo -e "fasta\t$FASTA"
  echo -e "outdir\t$OUTDIR"
  echo -e "awk_bin\t$AWK"
} > "$MANIFEST"

SPECIFIC="$OUTDIR/clusters_species_specific.tsv"
DROPPED="$OUTDIR/clusters_dropped_multi_species.tsv"
META="$OUTDIR/species_specific_clusters.meta.tsv"

# ---- sanitize SSC (strip CR) -------------------------------------------------
SSC_SAN="$OUTDIR/.ssc_sanitized.$$"
trap 'rm -f "$SSC_SAN"' EXIT
sed 's/\r$//' "$SSC" > "$SSC_SAN"

# ---- ensure FASTA indexed ----------------------------------------------------
[ -s "$FASTA.fai" ] || samtools faidx "$FASTA"

# ---- [1/2] Build META from SSC ----------------------------------------------
echo "[1/2] Building meta → $META"
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
    print "ERROR: SSC requires gene_name,species,gff_file,scaffold_name,clustering_id" > "/dev/stderr"; exit 3
  }
  next
}
{
  cl=$cid; s=$sp; gname=$gn; gff=$gf; scf=$sc
  # trims
  gsub(/^[[:space:]]+|[[:space:]]+$/, "", cl); gsub(/^[[:space:]]+|[[:space:]]+$/, "", s)
  gsub(/^[[:space:]]+|[[:space:]]+$/, "", gname); gsub(/^[[:space:]]+|[[:space:]]+$/, "", gff)
  gsub(/^[[:space:]]+|[[:space:]]+$/, "", scf)

  key = cl SUBSEP s
  keys[key]=1

  # accumulate unique gene names per (cluster,species)
  if (gname != ""){
    u=key SUBSEP gname
    if(!(u in seen_g)){
      seen_g[u]=1
      genes[key]=(genes[key]==""?gname:genes[key]"~~~"gname)
    }
  }

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
' "$SSC_SAN" | sort -t $'\t' -k2,2 -k1,1 > "$META"

# empty placeholder (reserved for future filtering)
: > "$DROPPED"

# convenience pairs (cluster_id<TAB>species)
cut -f1,2 "$META" | tail -n +2 | sort -u > "$SPECIFIC"

# ---- [2/2] Stream FASTA and write per-species files --------------------------
echo "[2/2] Extracting sequences by cluster-id → per-species FASTAs in $OUTDIR"

# clean any previous outputs to avoid accidental append
rm -f "$OUTDIR"/*_ssc_candidates.fa "$OUTDIR"/*_ssc_candidates.blastsafe.fa" 2>/dev/null || true

"$AWK" -v FS="\t" -v OFS="\t" -v META="$META" -v OUTDIR="$OUTDIR" '
BEGIN{
  # load meta into maps
  while((getline < META)>0){
    if(++r==1) continue
    cl=$1; sp=$2; gf=$3; sc=$4; genes=$5

    key = cl SUBSEP sp
    gf_of[key]=gf; sc_of[key]=sc; genes_of[key]=genes

    spkey=sp; gsub(/ /,"_",spkey)
    if(!(sp in human_path)) human_path[sp]=OUTDIR "/" spkey "_ssc_candidates.fa"
    if(!(sp in safe_path )) safe_path [sp]=OUTDIR "/" spkey "_ssc_candidates.blastsafe.fa"

    cl_to_sp[cl] = (cl in cl_to_sp ? cl_to_sp[cl] SUBSEP sp : sp)
  }
  close(META)
}

# FASTA headers: >CLUSTER_ID
/^>/{
  id=$1; sub(/^>/,"",id)
  if (!(id in cl_to_sp)) { emit=0; next }
  emit=1; cur_id=id

  # write headers for each species using this cluster
  nh = split(cl_to_sp[id], splist, SUBSEP)
  Hn=0; Sn=0
  for (i=1;i<=nh;i++){
    sp = splist[i]
    key = id SUBSEP sp
    OH = human_path[sp]
    OS = safe_path[sp]

    # human header (spaces allowed)
    print ">" id " " genes_of[key] " " sp " " gf_of[key] " " sc_of[key] >> OH

    # BLAST-safe header (no spaces)
    sps=sp; gsub(/ /,"_",sps)
    gfs=gf_of[key]; gsub(/ /,"_",gfs)
    scs=sc_of[key]; gsub(/ /,"_",scs)
    gens=genes_of[key]; gsub(/ /,"_",gens)
    print ">cl:" id "|genes:" gens "|sp:" sps "|gff:" gfs "|scaf:" scs >> OS

    H_paths[++Hn]=OH
    S_paths[++Sn]=OS
  }
  next
}

emit{
  for (i=1;i<=Hn;i++) print >> H_paths[i]
  for (i=1;i<=Sn;i++) print >> S_paths[i]
}
' "$FASTA"

echo "[DONE]"
echo "  Meta                   : $META"
echo "  Cluster-species pairs  : $(wc -l < "$SPECIFIC" || echo 0)"
echo "  Dropped (placeholder)  : $DROPPED"
echo "  Per-species FASTAs     : $OUTDIR"
echo "  Tip: Use *_blastsafe.fa for BLAST (no spaces in qseqid)."
