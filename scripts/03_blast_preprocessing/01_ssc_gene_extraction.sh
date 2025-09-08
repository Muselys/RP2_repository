#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# Twilight SSC post-processing (final)
#
# What this does
# 1) Filter classification.tab to rows where specific_class == "Species specific core".
#    Parse species from 'details' (take part after "Core:", split on '+').
#    Output: species_specific_core.tab  [gene_names, species]
#
# 2) expand gene names, Tokenize gene_names by "~~~" to get individual tokens (unique list).
#
# 3) Parse gene_data.csv once (CSV-safe via gawk FPAT) and annotate tokens:
#      find gene_name in gene_data.csv first, if it matches with gene_name in species_specific_core.tab,
#      extract: gff_file, scaffold_name, clustering_id; keep ALL clustering_ids if multi-mapped
#    Output: species_specific_core_annotations.tab
#            [gene_name (token), species, clustering_id, gff_file, scaffold_name]
#
# 4) Cluster-level non-unique report (cluster_id appears in >1 species):
#    Output: nonunique_genes_report.tsv
#            [gene_names (union), cluster_id, species_list, species_count]
#    THEN remove those cluster_ids from species_specific_core_annotations.tab
#
# 5) Per-species counts of distinct clusters, before vs unique-only:
#    Output: genes_per_species_counts.MERGED.tsv
#            [species, count_before, count_after, delta, retained_fraction]
#
# Logs
# - unmapped_gene_names.txt      : list of token genes that don't match the gene_names in gene_data.csv
# - ambiguous_gene_names.txt     : gene-name tokens mapping to >1 clustering_id (we KEEP ALL)
# -----------------------------------------------------------------------------

#BSUB -q normal
#BSUB -J twilight_ssc_extract
#BSUB -n 1
#BSUB -M 32000
#BSUB -R "select[mem>32000] rusage[mem=32000]"
#BSUB -o /data/pam/team230/sm71/scratch/rp2/logs/ssc_extract.%J.out
#BSUB -e /data/pam/team230/sm71/scratch/rp2/logs/ssc_extract.%J.err

set -euo pipefail
export LC_ALL=C

# inputs
IN_CLASS="/data/pam/team230/sm71/scratch/rp2/twilight_analysis/classification.tab"    # TSV
GD_CSV="/data/pam/team230/sm71/scratch/rp2/panaroo_output/gene_data.csv"              # CSV

# output dir (absolute)
OUTDIR="/data/pam/team230/sm71/scratch/rp2/blast_preprocessing"
mkdir -p "$OUTDIR"

# outputs (now absolute paths under OUTDIR)
RAW_CORE="$OUTDIR/species_specific_core.tab"                          # [gene_names, species]
ANNOT="$OUTDIR/species_specific_core_annotations.tab"                 # [gene_name, species, clustering_id, gff_file, scaffold_name]
NONUNIQ="$OUTDIR/nonunique_genes_report.tsv"                          # [gene_names (union), cluster_id, species_list, species_count]
COUNTS="$OUTDIR/genes_per_species_counts.MERGED.tsv"                  # summary
UNMAPPED="$OUTDIR/unmapped_gene_names.txt"
AMBIG="$OUTDIR/ambiguous_gene_names.txt"

# tooling
GAWK="$(command -v gawk || true)"           # gawk for CSV FPAT
MAWK="$(command -v mawk || command -v awk)" # mawk/awk for TSV
if [[ -z "$GAWK" ]]; then
  echo "ERROR: need gawk (CSV FPAT). Try: module load gawk" >&2
  exit 127
fi

# temps
TOKENS="$(mktemp)"         # unique gene-name tokens
MAP="$(mktemp)"            # name -> (cid,gff,scaf) rows
C2N="$(mktemp)"            # cid -> union gene names
ANNOT_ALL="$(mktemp)"      # pre-filter annotations (for counts-before & non-unique)
SLIM="$(mktemp)"           # [gene_names, species]
trap 'rm -f "$TOKENS" "$MAP" "$C2N" "$ANNOT_ALL" "$SLIM"' EXIT

# -------------------------------
# STEP 1 — SSC filter -> RAW_CORE
# -------------------------------
echo "[STEP 1] start: filter SSC + parse species → $SLIM"
"$MAWK" -F'\t' -v OFS='\t' '
NR==1{
  for(i=1;i<=NF;i++){
    h=$i; gsub(/^[ \t]+|[ \t]+$/,"",h)
    if(h=="specific_class") sc=i
    if(h=="gene_name")      gn=i
    if(h=="details")        det=i
  }
  if(!sc||!gn||!det){ print "ERROR: missing specific_class/gene_name/details" > "/dev/stderr"; exit 3 }
  print "gene_names","species"; next
}
{
  gsub(/\r$/,"")
  if($sc!="Species specific core") next

  d=$det
  sub(/.*Core:[[:space:]]*/, "", d)
  sub(/[[:space:]]+Inter:.*$/, "", d)
  sub(/[[:space:]]+Rare:.*$/,  "", d)
  gsub(/^[[:space:]]+|[[:space:]]+$/,"", d)
  if(d=="") next

  n=split(d, sp, /\+/)
  for(k=1;k<=n;k++){
    s=sp[k]; gsub(/^[[:space:]]+|[[:space:]]+$/, "", s)
    if(s!="") print $gn, s
  }
}
' "$IN_CLASS" > "$SLIM"
cp -f "$SLIM" "$RAW_CORE"
echo "[STEP 1] done : rows(incl header)=$(wc -l < "$RAW_CORE") → $RAW_CORE"

# -------------------------------------------
# STEP 2 — expand gene_names -> unique TOKENS
# -------------------------------------------
echo "[STEP 2] start: expand gene_names on ~~~ → unique tokens"
"$MAWK" -F'\t' '
NR==1{next}
{
  n=split($1, a, /~~~/)
  for(i=1;i<=n;i++){
    t=a[i]; gsub(/^[[:space:]]+|[[:space:]]+$/, "", t)
    if(t!="") print t
  }
}
' "$SLIM" | sort -u > "$TOKENS"
echo "[STEP 2] done : tokens=$(wc -l < "$TOKENS")"

# -------------------------------------------------------------------
# STEP 3 — CSV-safe parse gene_data.csv; build name -> (cid,gff,scaf)
#          keep ALL cids; log ambiguous names (>1 distinct cid)
# -------------------------------------------------------------------
: > "$UNMAPPED"
: > "$AMBIG"
echo "[STEP 3a] start: parse CSV (gawk FPAT) & build maps"
"$GAWK" -v OFS='\t' -v TOK="$TOKENS" -v MAP_OUT="$MAP" -v C2N_OUT="$C2N" -v AMBIG="$AMBIG" '
BEGIN{
  FPAT = "([^,]+)|(\"([^\"]|\"\")*\")"  # CSV-safe fields
  while((getline t < TOK)>0) need[t]=1
  close(TOK)
}
NR==1{
  for(i=1;i<=NF;i++){
    h=$i; gsub(/^"|"$/,"",h); gsub(/^[ \t]+|[ \t]+$/,"",h)
    if(h=="gff_file")       gi=i
    if(h=="scaffold_name")  si=i
    if(h=="clustering_id")  ci=i
    if(h=="gene_name")      ni=i
  }
  if(!gi||!si||!ci||!ni){ print "ERROR: gene_data.csv headers missing" > "/dev/stderr"; exit 5 }
  next
}
{
  gff=$gi; scf=$si; cid=$ci; name=$ni
  gsub(/^"|"$/,"",gff); gsub(/^[ \t]+|[ \t]+$/,"",gff)
  gsub(/^"|"$/,"",scf); gsub(/^[ \t]+|[ \t]+$/,"",scf)
  gsub(/^"|"$/,"",cid); gsub(/^[ \t]+|[ \t]+$/, "",cid)
  gsub(/^"|"$/,"",name); gsub(/^[ \t]+|[ \t]+$/, "",name)

  if(!(name in need)) next

  key=name SUBSEP cid SUBSEP gff SUBSEP scf
  if(!(seen[key]++)){
    print name, cid, gff, scf >> MAP_OUT
    if(!(seen_c2n[cid,name]++)) c2n[cid] = (c2n[cid]!=""? c2n[cid] "~~~" name : name)
  }
  if(!(hit[name,cid]++)){ cnt[name]++; list[name]=(list[name]!=""?list[name] "," cid:cid) }
}
END{
  for(c in c2n) print c, c2n[c] > C2N_OUT
  for(n in cnt) if(cnt[n]>1) print n "\t" cnt[n] "\t" list[n] > AMBIG
}
' "$GD_CSV"
echo "[STEP 3a] done : name→rows=$(wc -l < "$MAP") ; cid→names=$(wc -l < "$C2N") ; ambiguous_names=$(wc -l < "$AMBIG")"

# ----------------------------------------------------------
# annotate: species_specific_core_annotations.tab (ALL rows)
# ----------------------------------------------------------
echo "[STEP 3b] start: annotate tokens per species (ALL mappings)"
"$MAWK" -F'\t' -v OFS='\t' -v MAP="$MAP" -v UNM="$UNMAPPED" '
BEGIN{
  while((getline < MAP) > 0){
    n=$1; cid=$2; g=$3; s=$4
    i=++idx[n]; C[n,i]=cid; G[n,i]=g; S[n,i]=s
  }
  print "gene_name","species","clustering_id","gff_file","scaffold_name"
}
NR==1{ next }
{
  gn=$1; sp=$2
  m=split(gn, a, /~~~/)
  for(i=1;i<=m;i++){
    t=a[i]; gsub(/^[[:space:]]+|[[:space:]]+$/, "", t)
    if(t=="") continue
    if(!(t in idx)){ if(!logged[t]++){ print t > UNM } ; continue }
    for(k=1;k<=idx[t];k++){
      print t, sp, C[t,k], G[t,k], S[t,k]
    }
  }
}
' "$RAW_CORE" > "$ANNOT_ALL"
echo "[STEP 3b] done : annotated_rows(incl header)=$(wc -l < "$ANNOT_ALL") ; unmapped_tokens=$(wc -l < "$UNMAPPED")"

# ----------------------------------------------------------------------
# STEP 4 — cluster-level non-unique report + remove from annotations
# ----------------------------------------------------------------------
echo "[STEP 4] start: build non-unique report & filter annotations"
NONUNIQ_CID_SET="$(mktemp)"
"$MAWK" -F'\t' -v OFS='\t' -v C2N="$C2N" -v CIDSET="$NONUNIQ_CID_SET" '
function addset(k,v,   kk){ kk=k SUBSEP v; if(!(seen[kk]++)) list[k]=(k in list && list[k]!=""? list[k] "," v : v) }
BEGIN{
  while((getline < C2N) > 0){ cid=$1; g=$2; c2names[cid]=g }
  print "gene_names","cluster_id","species_list","species_count"
}
NR==1{next}
{
  cid=$3; sp=$2
  addset(cid, sp)
}
END{
  for(c in list){
    tmp=list[c]; gsub(/[^,]/,"",tmp); cnt=(tmp==""?1:length(tmp)+1)
    if(cnt>1){
      print (c2names[c]?c2names[c]:""), c, list[c], cnt
      print c > CIDSET
    }
  }
}
' "$ANNOT_ALL" > "$NONUNIQ"

# Filter annotations to remove non-unique CIDs
"$MAWK" -F'\t' -v OFS='\t' -v CIDSET="$NONUNIQ_CID_SET" '
BEGIN{ while((getline c < CIDSET)>0) bad[c]=1; close(CIDSET) }
NR==1{ print; next }
{ if(!($3 in bad)) print }
' "$ANNOT_ALL" > "$ANNOT"
rm -f "$NONUNIQ_CID_SET"
echo "[STEP 4] done : nonunique_rows(incl header)=$(wc -l < "$NONUNIQ") ; kept_annotations(incl header)=$(wc -l < "$ANNOT")"

# ---------------------------------------------------------
# STEP 5 — counts (distinct CIDs per species), before vs after
# ---------------------------------------------------------
echo "[STEP 5] start: compute before/after counts"
CB="$(mktemp)"
"$MAWK" -F'\t' -v OFS='\t' '
NR==1{next}
{ k=$2 SUBSEP $3; if(!(seen[k]++)) c[$2]++ }
END{ print "species","count_before"; for(s in c) print s, c[s] }
' "$ANNOT_ALL" | sort -t$'\t' -k1,1 > "$CB"

CA="$(mktemp)"
"$MAWK" -F'\t' -v OFS='\t' '
NR==1{next}
{ k=$2 SUBSEP $3; if(!(seen[k]++)) c[$2]++ }
END{ print "species","count_after"; for(s in c) print s, c[s] }
' "$ANNOT" | sort -t$'\t' -k1,1 > "$CA"

"$MAWK" -F'\t' -v OFS='\t' '
FNR==NR{ if($1!="species") before[$1]=$2+0; next }
$1!="species"{
  sp=$1; ca=$2+0; seenA[sp]=1
  cb=(sp in before? before[sp] : 0)
  delta=ca-cb
  rf=(cb>0? ca/cb : (ca>0?1:0))
  printf "%s\t%d\t%d\t%d\t%.6f\n", sp, cb, ca, delta, rf
}
END{
  for(sp in before) if(!(sp in seenA)){
    cb=before[sp]
    printf "%s\t%d\t%d\t%d\t%.6f\n", sp, cb, 0, -cb, 0
  }
}
' "$CB" "$CA" | sort -t$'\t' -k1,1 > "$COUNTS"
rm -f "$CB" "$CA"
echo "[STEP 5] done : wrote $COUNTS"

# ----------------
# final summary
# ----------------
echo "Wrote to $OUTDIR:"
echo "  - $RAW_CORE"
echo "  - $ANNOT   (non-unique cluster_ids removed)"
echo "  - $NONUNIQ"
echo "  - $COUNTS"
[[ -s "$UNMAPPED" ]] && echo "  - $UNMAPPED (unmapped tokens present)" || echo "  - no unmapped tokens"
[[ -s "$AMBIG"   ]] && echo "  - $AMBIG (gene-name tokens with multiple CIDs)" || echo "  - no ambiguous tokens"
