#!/usr/bin/env bash
set -euo pipefail

# -------- usage --------
usage() {
  cat <<EOF
Usage:
  $(basename "$0") -t REF_FASTA_DIR -s SUBSET_TSV -o OUT_DIR  BLAST1.out[.gz] [BLAST2 ...]
  or with globs:
  $(basename "$0") -t REF_FASTA_DIR -s SUBSET_TSV -o OUT_DIR  /path/to/*/*_ssc_candidates.out.gz

Required:
  -t  Directory with reference FASTAs (filenames are sample IDs; .fa/.fasta)
  -s  subset.tsv (header: sample<tab>species at least)
  -o  Output directory to write per-file results + summary

Notes:
  * BLAST files must be outfmt 6 with sseqid in column 2
  * Handles .out and .out.gz
EOF
  exit 1
}

# -------- args --------
REF_FASTA_DIR=""
SUBSET=""
OUT_DIR=""
while getopts ":t:s:o:" opt; do
  case "$opt" in
    t) REF_FASTA_DIR="$OPTARG" ;;
    s) SUBSET="$OPTARG" ;;
    o) OUT_DIR="$OPTARG" ;;
    *) usage ;;
  }
done

shift $((OPTIND-1))
[[ -d "$REF_FASTA_DIR" && -f "$SUBSET" && -n "${OUT_DIR:-}" ]] || usage
[[ $# -ge 1 ]] || usage

mkdir -p "$OUT_DIR"

# -------- build truth once --------
SUMMARY_DIR="$OUT_DIR/summary"
mkdir -p "$SUMMARY_DIR"

SAMPLES_IN_REF="$SUMMARY_DIR/samples_in_ref.txt"
S2S_ALL="$SUMMARY_DIR/sample2species_all.tsv"
TRUTH_ALL="$SUMMARY_DIR/sample_truth.tsv"
SUMMARY_TSV="$SUMMARY_DIR/confusion_summary.tsv"

# samples present in ref (strip .fa/.fasta)
ls "$REF_FASTA_DIR" | sed -E 's/\.(fa|fasta)$//' | sort -u > "$SAMPLES_IN_REF"

# sample -> species from subset (skip header)
awk -F'\t' 'NR>1 {print $1"\t"$2}' "$SUBSET" | sort -u > "$S2S_ALL"

# keep only samples present in ref
awk 'NR==FNR{keep[$1]=1; next} ($1 in keep)' "$SAMPLES_IN_REF" "$S2S_ALL" > "$TRUTH_ALL"

# header for summary
echo -e "blast_file\tspecies_inferred\tTP\tFP\tTN\tFN\tSensitivity\tSpecificity" > "$SUMMARY_TSV"

# helper: infer species name from parent dir or filename; canonicalize like "genus species"
canon_species() {
  local s="$1"
  echo "$s" | tr '_' ' ' | awk '{print tolower($0)}' | sed -E 's/[[:space:]]+/ /g; s/^ | $//g'
}

infer_species_from_path() {
  local p="$1" bn dir
  bn="$(basename "$p")"
  dir="$(basename "$(dirname "$p")")"
  # prefer directory name if it looks like Genus_species
  if [[ "$dir" =~ ^[A-Za-z]+_[A-Za-z0-9]+$ ]]; then
    echo "$dir"
  else
    # fallback: grab leading Genus_species-like token from filename
    echo "$bn" | sed -E 's/^([A-Za-z]+_[A-Za-z0-9]+).*/\1/'
  fi
}

# -------- process each BLAST file --------
for BLAST in "$@"; do
  if [[ ! -s "$BLAST" ]]; then
    echo "[warn] Skipping missing/empty: $BLAST" >&2
    continue
  fi

  species_name="$(infer_species_from_path "$BLAST")"
  target_canon="$(canon_species "$species_name")"
  out_base="$OUT_DIR/$(basename "$species_name")"
  mkdir -p "$out_base"

  echo "==> $BLAST  (species: $species_name → $target_canon)"

  # choose reader
  if [[ "$BLAST" == *.gz ]]; then READER="zcat"; else READER="cat"; fi

  # 1) predicted positive samples from BLAST (sseqid col2 → take prefix before first '.')
  PRED_SAMP="$out_base/predicted_pos_samples.txt"
  $READER "$BLAST" \
    | awk -F'\t' 'NF>=2 {split($2,a,"."); print a[1]}' \
    | sort -u > "$PRED_SAMP"

  # 2) join truth + predictions to label POS/NEG
  LABELLED="$out_base/sample_truth_with_pred.tsv"
  awk 'BEGIN{FS=OFS="\t"} NR==FNR{pos[$1]=1; next} {print $1, $2, ( ($1 in pos) ? "POS" : "NEG")}' \
    "$PRED_SAMP" "$TRUTH_ALL" > "$LABELLED"

  # 3) confusion matrix using your canon logic
  CM_OUT="$out_base/confusion_matrix.txt"
  awk -v FS="\t" -v target="$target_canon" '
  function canon(s,t){ t=tolower(s); gsub(/[_]+/," ",t); gsub(/[[:space:]]+/, " ", t); sub(/^ | $/,"",t); return t }
  function genus_strip(g,h){ h=g; sub(/_[a-z]+$/,"",h); return h }
  function is_target(label,part,i){
    n=split(label,part,/;/);
    for(i=1;i<=n;i++){
      tok=canon(part[i]); split(tok,a," ");
      if(length(a)<2) continue;
      gen=genus_strip(a[1]); spp=a[2];
      if(gen" "spp==target) return 1
    }
    return 0
  }
  BEGIN{tp=fp=tn=fn=0}
  {
    st=$2; pred=($3=="POS"); tgt=is_target(st);
    if (tgt && pred) tp++;
    else if (tgt && !pred) fn++;
    else if (!tgt && pred) fp++;
    else tn++;
  }
  END{
    printf "TP\t%d\nFP\t%d\nTN\t%d\nFN\t%d\n", tp, fp, tn, fn;
    sens = (tp+fn)? tp/(tp+fn) : 0;
    spec = (tn+fp)? tn/(tn+fp) : 0;
    printf "Sensitivity\t%.6f\nSpecificity\t%.6f\n", sens, spec;
  }' "$LABELLED" | tee "$CM_OUT" >/dev/null

  # append to summary
  TP=$(awk '/^TP\t/{print $2}' "$CM_OUT");   FP=$(awk '/^FP\t/{print $2}' "$CM_OUT")
  TN=$(awk '/^TN\t/{print $2}' "$CM_OUT");   FN=$(awk '/^FN\t/{print $2}' "$CM_OUT")
  SENS=$(awk '/^Sensitivity\t/{print $2}' "$CM_OUT"); SPEC=$(awk '/^Specificity\t/{print $2}' "$CM_OUT")
  echo -e "$(basename "$BLAST")\t$species_name\t$TP\t$FP\t$TN\t$FN\t$SENS\t$SPEC" >> "$SUMMARY_TSV"
done

echo "[ok] Summary written: $SUMMARY_TSV"
