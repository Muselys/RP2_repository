#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# 03_run_twilight.sh — submit Twilight analysis as an LSF job
#
# Purpose:
#   Run classify_genes.R (Twilight species-level analysis) on the cluster.
#   Uses the filtered Panaroo RTAB and the groups.tab file.
#   Requests ~300 GB RAM on a single host, activates R conda env, and writes
#   outputs to twilight_output with logs in logs/.
# -----------------------------------------------------------------------------

#bsub -q normal \
#     -o logs/twilight_%J.out \
#     -e logs/twilight_%J.err \
#     -R "select[mem>300000] rusage[mem=300000] span[hosts=1]" \
#     -M 300000 \
#     bash -lc 'module load ISG/conda && conda activate rconda && \
#     Rscript scripts/02_twilight_analysis/01_classify_genes.R \
#         -p pan_genome_processing/gene_presence_absence_filtered.Rtab \
#         -g pan_genome_processing/groups.tab \
#         -o twilight_output \
#         -s 10 -c 0.95 -r 0.15'
