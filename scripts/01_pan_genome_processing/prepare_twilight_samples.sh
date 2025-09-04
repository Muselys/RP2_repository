#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# This script prepares input sample lists for twilight analysis.
# It filters Panaroo gene presence/absence data using metadata,
# keeps only samples from a fixed 18-species allowlist,
# and applies a high-quality (HQ==TRUE) filter.
# It outputs a final sample list, species counts at each step,
# and a summary log for QC.
# -----------------------------------------------------------------------------
#BSUB -J twilight_prep
#BSUB -o logs/twilight_prep.%J.out
#BSUB -e logs/twilight_prep.%J.err
#BSUB -R 'select[mem>8000] rusage[mem=8000]'
#BSUB -n 1
#BSUB -M 8000
#BSUB -W 01:00
set -euo pipefail
export LC_ALL=C

# ---- paths ----
Rtab="/data/pam/team230/sm71/scratch/rp2/panaroo_output/gene_presence_absence.Rtab"
META="/data/pam/team230/sm71/scratch/rp2/File4_QC_characterisation_661K.tsv"
OUT_DIR="/data/pam/team230/sm71/scratch/rp2/twilight_input"

STEP1="$OUT_DIR/samples_step1_ids.tsv"
STEP2="$OUT_DIR/samples_step2_allowed.tsv"
STEP3="$OUT_DIR/samples_step3_HQ.tsv"
FINAL="$OUT_DIR/samples_twilight.tsv"

SUMMARY="$OUT_DIR/samples_twilight.summary.txt"
COUNTS_TXT="$OUT_DIR/species_counts_after.txt"
COUNTS_STEP2="$OUT_DIR/species_counts_step2.tsv"
COUNTS_STEP3="$OUT_DIR/species_counts_step3.tsv"

mkdir -p "$OUT_DIR" logs
: > "$SUMMARY"

[[ -s "$Rtab" ]] || { echo "ERR: missing $Rtab" >&2; exit 1; }
[[ -s "$META" ]] || { echo "ERR: missing $META" >&2; exit 1; }

# ---- 18 species allowlist ----
read -r -d '' allowed <<'EOF' || true
Enterococcus faecalis
Enterococcus faecium
Staphylococcus argenteus
Staphylococcus aureus
Staphylococcus capitis
Staphylococcus epidermidis
Staphylococcus haemolyticus
Staphylococcus pseudintermedius
Staphylococcus sciuri
Streptococcus agalactiae
Streptococcus dysgalactiae
Streptococcus equi
Streptococcus mitis
Streptococcus mutans
Streptococcus pneumoniae
Streptococcus pyogenes
Streptococcus suis
Streptococcus uberis
EOF

############################################
# STEP 1: extract sample IDs from Rtab header
############################################
awk -F '\t' 'NR==1{for(i=2;i<=NF;i++) print $i}' "$Rtab" \
 | awk 'BEGIN{print "sample_id"}1' > "$STEP1"

ids_step1=$(( $(wc -l < "$STEP1") - 1 ))
echo "[STEP1] IDs found in Rtab header: $ids_step1 → $STEP1" | tee -a "$SUMMARY"

############################################
# STEP 2: join species from META; keep only 18 species
# (auto-detect sample_id & species by header name)
############################################
TMP_LOG="$(mktemp)"
awk -F '\t' -v OFS='\t' -v out="$STEP2" -v logf="$TMP_LOG" -v allow_block="$allowed" '
function trim(x){gsub(/^[ \t]+|[ \t]+$/,"",x); return x}

BEGIN{
  print "sample_id\tspecies" > out

  # Build allow-set after normalizing spaces -> underscores
  n = split(allow_block, L, /\n/)
  for(i=1;i<=n;i++){
    s = trim(L[i]); if(s=="") continue
    ns = s; gsub(/[[:space:]]+/, "_", ns)
    allow[ns]=1
  }
}

# First input (META): header → locate columns
FNR==1 && NR==FNR {
  sidcol=0; spcol=0
  for(i=1;i<=NF;i++){
    h=tolower($i)
    if(h=="sample_id") sidcol=i
    if(h=="species")   spcol=i
  }
  if(!sidcol) sidcol=1
  if(!spcol)  spcol=2
  next
}

# First input (META): map ID -> species
NR==FNR {
  id=trim($(sidcol)); sp=trim($(spcol))
  if(id!=""){
    nsp=sp; gsub(/[[:space:]]+/, "_", nsp)
    SP[id]=sp; NSP[id]=nsp
  }
  next
}

# Second input (STEP1): keep only allowed species
FNR==1 && NR!=FNR { next } # skip STEP1 header
{
  id=$1; sp=SP[id]; nsp=NSP[id]
  if(sp==""){ missing[id]++; next }
  if(nsp in allow) { print id, sp >> out; kept++ }
  else { removed_by_species[sp]++; }
}
END{
  print "[STEP2] removed not-allowed per species:" > logf
  total=0
  for(s in removed_by_species){ print s "\t" removed_by_species[s] >> logf; total+=removed_by_species[s] }
  if(total==0) print "(none)\t0" >> logf
  m=0; for(x in missing) m++
  print "[STEP2] samples with missing species (treated as removed):\t" m >> logf
  print "[STEP2] kept after allowlist:\t" kept >> logf
}
' "$META" "$STEP1"

remain_step2=$(( $(wc -l < "$STEP2") - 1 ))
removed_step2=$(( ids_step1 - remain_step2 ))
{
  echo -e "species\tremoved_not_allowed"
  awk 'NR>1' "$TMP_LOG" | sed -n '1!p' | grep -v '^\[STEP2\]'
  echo "—"
  echo "[STEP2] removed (total): $removed_step2"
  echo "[STEP2] remaining: $remain_step2"
} | tee -a "$SUMMARY" >/dev/null
rm -f "$TMP_LOG"

# Per-species counts after STEP2
awk -F '\t' 'NR>1 { cnt[$2]++ } END{ for(s in cnt) print s "\t" cnt[s] }' "$STEP2" \
 | sort > "$COUNTS_STEP2"
echo "[STEP2] species counts → $COUNTS_STEP2" | tee -a "$SUMMARY"
cat "$COUNTS_STEP2" | tee -a "$SUMMARY" >/dev/null

############################################
# STEP 3: add high_quality from META; keep HQ==TRUE
# (detect high_quality by header name; no hardcoded index)
############################################
TMP_ANN="$(mktemp)"
awk -F '\t' -v OFS='\t' -v out="$TMP_ANN" '
function trim(x){gsub(/^[ \t]+|[ \t]+$/,"",x); return x}

# META header: find sample_id + high_quality
FNR==1 && NR==FNR {
  sidcol=hqcol=0
  for(i=1;i<=NF;i++){
    h=tolower($i)
    if(h=="sample_id")   sidcol=i
    if(h=="high_quality") hqcol=i
  }
  if(!sidcol){ sidcol=1 }
  if(!hqcol){  print "ERROR: high_quality column not found" > "/dev/stderr"; exit 2 }
  next
}

# META rows: map ID -> HQ
NR==FNR {
  id=trim($(sidcol)); v=trim($(hqcol))
  if(id!="") HQ[id]=v
  next
}

# Second file (STEP2): append HQ
FNR==1 && NR!=FNR { print "sample_id\tspecies\thigh_quality" > out; next }
{
  id=$1; sp=$2; v=HQ[id]; print id, sp, v >> out
}
' "$META" "$STEP2"

# Filter HQ==TRUE (case-insensitive) + log removals per species
TMP_FIL="$(mktemp)"; TMP_REMOVED="$(mktemp)"
awk -F '\t' -v OFS='\t' -v out="$TMP_FIL" -v logf="$TMP_REMOVED" '
FNR==1 { print $0 > out; next }
{
  v=tolower($3)
  if(v=="true"){ print $1, $2, $3 >> out; kept++ }
  else { removed[$2]++ }
}
END{
  print "[STEP3] removed HQ!=TRUE per species:" > logf
  total=0; for(s in removed){ print s "\t" removed[s] >> logf; total+=removed[s] }
  if(total==0) print "(none)\t0" >> logf
  print "[STEP3] kept after HQ filter:\t" kept >> logf
}
' "$TMP_ANN"

mv "$TMP_FIL" "$STEP3"; rm -f "$TMP_ANN"

# Per-species counts after STEP3
awk -F '\t' 'NR>1 { cnt[$2]++ } END{ for(s in cnt) print s "\t" cnt[s] }' "$STEP3" \
 | sort > "$COUNTS_STEP3"

{
  echo -e "species\tremoved_high_quality_not_true"
  awk 'NR>1' "$TMP_REMOVED"
  echo "—"
  echo -e "species\tremaining_after_HQ_filter"
  cat "$COUNTS_STEP3"
} | tee -a "$SUMMARY" >/dev/null
rm -f "$TMP_REMOVED"

############################################
# FINAL: sample_id \t species + full 18-species counts
############################################
awk -F'\t' 'NR==1{print "sample_id\tspecies"; next} {print $1 "\t" $2}' "$STEP3" > "$FINAL"
kept=$(( $(wc -l < "$FINAL") - 1 ))
echo "[FINAL] kept $kept samples → $FINAL" | tee -a "$SUMMARY"

# Expand counts to all 18 species (fill zeros)
ALLOWED_ZERO="$OUT_DIR/_allowed_zero.tsv"
COUNTS_ACTUAL="$OUT_DIR/_counts_actual.tsv"
echo "$allowed" | awk -v OFS='\t' '{print $0, 0}' | sort > "$ALLOWED_ZERO"
sort -o "$COUNTS_ACTUAL" "$COUNTS_STEP3" 2>/dev/null || : > "$COUNTS_ACTUAL"
join -t $'\t' -a1 -e 0 -o 1.1,2.2 "$ALLOWED_ZERO" "$COUNTS_ACTUAL" > "$COUNTS_TXT"
rm -f "$ALLOWED_ZERO" "$COUNTS_ACTUAL"

echo "[COUNTS] final species counts (18 species, zeros included) → $COUNTS_TXT" | tee -a "$SUMMARY"
echo "Summary written to $SUMMARY"
