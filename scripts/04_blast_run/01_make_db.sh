#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# Builds a BLAST+ nucleotide database from a combined reference FASTA on LSF.
#
# What it does:
# - Submits an LSF job that loads BLAST 2.14.1 and runs `makeblastdb`.
# - Input FASTA: /data/.../blast/references/combined_references.fa
# - Output DB prefix: /data/.../blast/references/db  (BLAST+ will create .nhr/.nin/.nsq,
#   and multi-volume files like db.00.nsq, db.01.nsq, etc. if large)
# - Logs stdout/stderr to: /data/.../blast/logs/makeblastdb.<JOBID>.out|.err
#
# LSF resources:
# - Queue: normal
# - Cores: 1
# - Memory: 4 GB (select[mem>4000] rusage[mem=4000], -M 4000)
#
# Notes:
# - `-dbtype nucl` builds a nucleotide DB (use this for nucleotide BLASTs).
# - Use the DB by passing `-db /data/.../blast/references/db` to blastn, etc.
# - Very large inputs will be split into multiple volumes automatically.
# - Module path/version may differ on other clusters; adjust `module load` as needed.
#
# Run metadata (example from a previous run):
# - Start→Finish: ~56 min walltime
# - Host: node-14-17 (queue: normal)
# - Sequences: ~10M
# - Output DB size: ~65 GB (23 volumes)
# - CPU time: ~56 min; Max RSS: ~107 MB
# -----------------------------------------------------------------------------

bsub -q normal\
  -J makeblastdb \
  -n 1 -R "select[mem>4000] rusage[mem=4000]" \
  -M 4000 \
  -o /data/pam/team230/sm71/scratch/rp2/blast/logs/makeblastdb.%J.out \
  -e /data/pam/team230/sm71/scratch/rp2/blast/logs/makeblastdb.%J.err \
  "bash -lc 'module load blast/2.14.1--pl5321h6f7f691_0 && \
    makeblastdb -in /data/pam/team230/sm71/scratch/rp2/blast/references/combined_references.fa \
      -dbtype nucl \
      -out /data/pam/team230/sm71/scratch/rp2/blast/references/db'"

