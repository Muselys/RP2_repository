./01_run_blast.sh \
  --query /path/to/queries.fa \
  --db /path/to/ref19_norm/ref19_norm \
  --out /path/to/results/blast_results.tsv \
  [--threads 8] \
  [--task megablast] \
  [--evalue 1e-6] \
  [--module blast/2.14.1--pl5321h6f7f691_0]
