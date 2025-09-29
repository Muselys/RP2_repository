#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# trim_gene_data.csv.sh  (TSV-only)
#
# Makes a slim TSV from Panaroo gene_data.csv containing only SSC tokens and
# only the columns: gff_file, scaffold_name, clustering_id, gene_name.
#
# Inputs:
#   classification.tab (TSV)  : /data/pam/team230/sm71/scratch/rp2/twilight_analysis/classification.tab
#   gene_data.csv      (CSV)  : /data/pam/team230/sm71/scratch/rp2/panaroo_output/gene_data.csv
#
# Outputs (TSV only, plus a tiny log):
#   /data/pam/team230/sm71/scratch/rp2/blast_preprocessing/gene_data_copy.tsv
#   /data/pam/team230/sm71/scratch/rp2/blast_preprocessing/gene_data_copy.unmapped_tokens.txt
#
# Notes:
#   - mawk for TSV steps; gawk for CSV FPAT.
#   - Very low RAM; I/O bound.
# -----------------------------------------------------------------------------

#BSUB -q normal
#BSUB -J trim_gene_data_csv
#BSUB -n 1
#BSUB -M 8000
#BSUB -R "select[mem>8000] rusage[mem=8000]"
#BSUB -o /data/pam/team230/sm71/scratch/rp2/logs/trim_gene_data.%J.out
#BSUB -e /data/pam/team230/sm71/scratch/rp2/logs/trim_gene_data.%J.err

set -euo pipefail
export LC_ALL=C

log(){ echo "[$(date +'%F %T')] $*"; }

# ---- Fixed paths ----
IN_CLASS="/data/pam/team230/sm71/scratch/rp2/twilight_analysis/classification.tab"
GD_CSV="/data/pam/team230/sm71/scratch/rp2/panaroo_output/gene_data.csv"
OUTDIR="/data/pam/team230/sm71/scratch/rp2/blast_preprocessing"
mkdir -p "$OUTDIR"

# ---- Tools ----
MAWK="$(command -v mawk || command -v awk)"
GAWK="$(command -v gawk || true)"
if [[ -z "$GAWK" ]]; then
  echo "ERROR: need gawk (for FPAT CSV parsing). Try: module load gawk" >&2
  exit 127
fi

# ---- Outputs ----
GD_COPY_TSV="$OUTDIR/gene_data_copy.tsv"
UNMAPPED_TOKENS="$OUTDIR/gene_data_copy.unmapped_tokens.txt"

# ---- Temps ----
TOKENS="$(mktemp)"
trap 'rm -f "$TOKENS"' EXIT

# ------------------------------------------------------------------
# Step 1 — Build SSC token list from classification.tab (TSV)
# ------------------------------------------------------------------
log "[STEP 1] Extract SSC tokens from classification.tab"
[[ -s "$IN_CLASS" ]] || { echo "ERR: missing $IN_CLASS" >&2; exit 1; }

"$MAWK" -F'\t' '
NR==1{
  for(i=1;i<=NF;i++){
    h=$i; gsub(/^[ \t]+|[ \t]+$/,"",h)
    if(h=="specific_class") sc=i
    if(h=="gene_name")      gn=i
  }
  if(!sc||!gn){ print "ERROR: need columns specific_class and gene_name" > "/dev/stderr"; exit 3 }
  next
}
{
  gsub(/\r$/,"")
  if($sc!="Species specific core") next
  n=split($gn, a, /~~~/)
  for(i=1;i<=n;i++){
    t=a[i]; gsub(/^[[:space:]]+|[[:space:]]+$/, "", t)
    if(t!="") print t
  }
}
' "$IN_CLASS" | sort -u > "$TOKENS"

log "  SSC tokens: $(wc -l < "$TOKENS")"
if [[ ! -s "$TOKENS" ]]; then
  echo "No SSC tokens found; nothing to do." >&2
  exit 0
fi

# ------------------------------------------------------------------
# Step 2 — Stream gene_data.csv → slim TSV + unmapped tokens log
# ------------------------------------------------------------------
log "[STEP 2] Stream gene_data.csv and write slim TSV to $OUTDIR"
: > "$GD_COPY_TSV"
: > "$UNMAPPED_TOKENS"

# TSV header
echo -e "gff_file\tscaffold_name\tclustering_id\tgene_name" > "$GD_COPY_TSV"

"$GAWK" -v TOK="$TOKENS" -v OUT_TSV="$GD_COPY_TSV" -v UNM="$UNMAPPED_TOKENS" '
BEGIN{
  FPAT = "([^,]+)|(\"([^\"]|\"\")*\")"  # CSV-safe fields
  while((getline t < TOK)>0){ need[t]=1 }
  close(TOK)
}
function dequote(x){ gsub(/^"|"$/,"",x); return x }
function trim(x){ gsub(/^[ \t]+|[ \t]+$/,"",x); return x }
NR==1{
  for(i=1;i<=NF;i++){
    h=trim(dequote($i))
    if(h=="gff_file")       gi=i
    else if(h=="scaffold_name")  si=i
    else if(h=="clustering_id")  ci=i
    else if(h=="gene_name")      ni=i
  }
  if(!gi||!si||!ci||!ni){
    print "ERROR: gene_data.csv missing required headers" > "/dev/stderr"; exit 5
  }
  next
}
{
  gff=trim(dequote($gi))
  scf=trim(dequote($si))
  cid=trim(dequote($ci))
  nam=trim(dequote($ni))
  if(!(nam in need)) next

  print gff "\t" scf "\t" cid "\t" nam >> OUT_TSV
  seen[nam]=1
}
END{
  for(n in need) if(!(n in seen)) print n > UNM
}
' "$GD_CSV"

log "  slim TSV rows (incl header): $(wc -l < "$GD_COPY_TSV")"
if [[ -s "$UNMAPPED_TOKENS" ]]; then
  log "  unmapped tokens: $(wc -l < "$UNMAPPED_TOKENS") (see $UNMAPPED_TOKENS)"
else
  log "  no unmapped tokens"
fi

log "[DONE] Wrote:"
echo "  $GD_COPY_TSV"
[[ -s "$UNMAPPED_TOKENS" ]] && echo "  $UNMAPPED_TOKENS"
