REF_FASTA_DIR=/data/pam/team230/sm71/scratch/rp2/atb/ref/fasta

ls "$REF_FASTA_DIR" | sed -E 's/\.(fa|fasta)$//' | sort -u > samples_in_ref.txt
wc -l samples_in_ref.txt
head samples_in_ref.txt



#make a truth set:
SUBSET=/data/pam/team230/sm71/scratch/rp2/atb/subset.tsv

# sample \t species (skip header), then keep only samples present in ref
awk -F'\t' 'NR>1 {print $1"\t"$2}' "$SUBSET" | sort -u > sample2species_all.tsv
awk 'NR==FNR{keep[$1]=1; next} ($1 in keep)' samples_in_ref.txt sample2species_all.tsv \
  > sample_truth.tsv

wc -l sample_truth.tsv
head sample_truth.tsv

awk -F'\t' '{
  sp=tolower($2); gsub(/[_]+/," ",sp); gsub(/[[:space:]]+/," ",sp); sub(/^ | $/,"",sp)
  print $1, sp
}' sample_truth.tsv > sample_truth.norm.tsv
mv sample_truth.norm.tsv sample_truth.tsv



#3) Build the prediction set from BLAST:
zcat /data/pam/team230/sm71/scratch/rp2/run_blast/results/Enterococcus_faecium/Enterococcus_faecium_ssc_candidates.out.gz \
| cut -f2 \
| sed 's/\.[Cc][Oo][Nn][Tt][Ii][Gg].*$//' \
| sort -u > predicted_pos_samples.txt

wc -l predicted_pos_samples.txt
head -50 predicted_pos_samples.txt



#4) Join truth and predictions, label POS/NEG:
awk 'BEGIN{FS=OFS="\t"}
     NR==FNR {pos[$1]=1; next}
     {print $1, $2, ( ($1 in pos) ? "POS" : "NEG")}
' predicted_pos_samples.txt sample_truth.tsv > sample_truth_with_pred.tsv

head -5 sample_truth_with_pred.tsv
# sample_id \t species_truth \t POS/NEG


#Confusion matrix
TARGET_CANON="enterococcus faecium"

awk -v FS="\t" -v target="$TARGET_CANON" '
function canon(s,t){ t=tolower(s); gsub(/[_]+/," ",t); gsub(/[[:space:]]+/, " ", t); sub(/^ | $/,"",t); return t }
function genus_strip(g,h){ h=g; sub(/_[a-z]+$/,"",h); return h }
function is_target(label,part,i,tok,gen,spp){
  n=split(label,part,/;/)
  for(i=1;i<=n;i++){
    tok=canon(part[i]); split(tok,a," ")
    if(length(a)<2) continue
    gen=genus_strip(a[1]); spp=a[2]
    if(gen" "spp==target) return 1
  }
  return 0
}
BEGIN{tp=fp=tn=fn=0}
{
  st=$2; pred=($3=="POS"); tgt=is_target(st)
  if (tgt && pred) tp++
  else if (tgt && !pred) fn++
  else if (!tgt && pred) fp++
  else tn++
}
END{
  printf "TP\t%d\nFP\t%d\nTN\t%d\nFN\t%d\n", tp, fp, tn, fn
  sens = (tp+fn)? tp/(tp+fn) : 0
  spec = (tn+fp)? tn/(tn+fp) : 0
  printf "Sensitivity\t%.6f\nSpecificity\t%.6f\n", sens, spec
}
' sample_truth_with_pred.tsv




./data/pam/team230/sm71/scratch/rp2/run_blast/results/eval_one_blast.sh \
  -t /data/pam/team230/sm71/scratch/rp2/atb/ref/fasta \
  -s /data/pam/team230/sm71/scratch/rp2/atb/subset.tsv \
  -b /data/pam/team230/sm71/scratch/rp2/run_blast/results/Enterococcus_faecalis/Enterococcus_faecalis_ssc_candidates.out.gz \
  -o /data/pam/team230/sm71/scratch/rp2/run_blast/results/
# (optional) -x "enterococcus faecium"  # if you want to override autodetected species
