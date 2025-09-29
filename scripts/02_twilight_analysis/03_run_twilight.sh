#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# Submits an LSF job to run classify_genes_1.R (species-level twilight analysis)
# using the filtered Panaroo Rtab and the new groups.tab.
# Requests ~300 GB RAM on a single host, sets up the R conda env, and writes
# outputs to twilight_output with stdout/stderr logs captured.
# -----------------------------------------------------------------------------

#Twilight analysis: species-level
bsub -q normal \
    -o /data/pam/team230/sm71/scratch/rp2/logs/twilight_%J.out \
    -e /data/pam/team230/sm71/scratch/rp2/logs/twilight_%J.err \
    -R "select[mem>300000] rusage[mem=300000] span[hosts=1]" \
    -M 300000 \
    bash -lc 'module load ISG/conda && conda activate rconda && \
    Rscript /data/pam/team230/sm71/scratch/rp2/twilight_analysis/classify_genes.R \
        -p /data/pam/team230/sm71/scratch/rp2/pan_genome_processing/gene_presence_absence_filtered.Rtab \
        -g /data/pam/team230/sm71/scratch/rp2/pan_genome_processing/groups.tab \
        -o s'

#122085 genes found in Twilight analysis