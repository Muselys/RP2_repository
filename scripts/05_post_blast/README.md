python3 02_join_truth.py \
  --blast /data/pam/team230/sm71/scratch/rp2/blast_run/results/blast.tsv \
  --truth-queries /data/pam/team230/sm71/scratch/rp2/blast_run/post_blast/truth_queries.tsv \
  --truth-refs /data/pam/team230/sm71/scratch/rp2/blast_run/post_blast/truth_refs.tsv \
  --out /data/pam/team230/sm71/scratch/rp2/blast_run/post_blast/blast_annotated.tsv


./03_filter_thresholds.sh IN.tsv[.gz] OUT.tsv.gz [PIDENT] [QCOVS]

python3 04_compute_metrics.py \
  --annotated /path/to/blast_annotated.tsv[.gz] \
  --filtered  /path/to/blast_filtered.tsv[.gz] \
  --outdir    /path/to/eval_dir \
  --species-file /path/to/targets_18.txt


python3 05_plot_metrics.py

python3 plot_metrics.py \
  --predictions /lustre/scratch127/pam/teams/team230/sm71/rp2/blast_run/post_blast/eval/predictions_pre_post.tsv \
  --species-file /lustre/scratch127/pam/teams/team230/sm71/rp2/blast_run/post_blast/targets_18.txt \
  --outdir /lustre/scratch127/pam/teams/team230/sm71/rp2/blast_run/post_blast/eval \
  --genes-top 300 \
  --max-genes 300 \
  --min-gene-count 1
