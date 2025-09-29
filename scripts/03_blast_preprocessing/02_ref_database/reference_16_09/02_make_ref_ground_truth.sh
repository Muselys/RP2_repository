cd /data/pam/team230/sm71/scratch/rp2/atb/

awk -F'\t' '{
  split($3,a,"/"); 
  id=a[length(a)]; 
  sub(/\.fa$/,"",id); 
  print id "\t" $1
}' subset_18.tsv > ref_ground_truth.tsv
