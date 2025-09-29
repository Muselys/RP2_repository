#concatenate blast outputs into one file
printf 'qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tqcovhsp\tqcovs\tqcovus\n' > blast.tsv
find . -maxdepth 1 -type f -name 'q_*.megablast.top10.tsv' -size +0c -print0 \
| sort -zV \
| xargs -0 cat >> blast.tsv

#need truth tables for qseqid and sseqid
Q_META="/data/pam/team230/sm71/scratch/rp2/blast_run/post_blast/truth_queries.tsv"
R_META="/data/pam/team230/sm71/scratch/rp2/blast_run/post_blast/truth_refs.tsv"
#join truth table columns to blast output
#apply thresholds and save a blast copy
#Classify 
#calculate specificity and sensitivity per species in 18 target species


