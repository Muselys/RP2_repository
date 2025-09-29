""" Extract qseqid and sseqid from the blast output
"""
import sys

if len(sys.argv) != 3:
    print(f"Usage: python {sys.argv[0]} <blast_output.txt> <output_hits.txt>")
    sys.exit(1)

blast_file = sys.argv[1]
output_file = sys.argv[2]

with open(blast_file) as f, open(output_file, 'w') as out:
    for line in f:
        if line.strip():
            parts = line.strip().split()
            if len(parts) >= 2:
                out.write(f"{parts[0]}\t{parts[1]}\n")


""" Create an empty sensitivity/specificity matrix for the blast output, 
only label the columns (sseqid) and the first row (sseqid)"""

import sys
import csv

if len(sys.argv) != 3:
    print(f"Usage: python {sys.argv[0]} blast_hits.txt empty_matrix.tsv")
    sys.exit(1)

hits_file = sys.argv[1]
output_file = sys.argv[2]

genes = set()
genomes = set()

with open(hits_file) as f:
    for line in f:
        if line.strip():
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                genes.add(parts[0])
                genomes.add(parts[1])

genes = sorted(genes)
genomes = sorted(genomes)

with open(output_file, "w", newline='') as out:
    writer = csv.writer(out, delimiter='\t')
    # header
    writer.writerow(["Gene"] + genomes)
    # empty cells
    for gene in genes:
        row = [gene] + [""] * len(genomes)
        writer.writerow(row)

print(f"Empty matrix saved to {output_file}")