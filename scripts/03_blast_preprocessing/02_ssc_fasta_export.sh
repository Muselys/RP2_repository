#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# 02_ssc_fasta_export.sh  —  export per-gene FASTAs per species
#
# INPUTS
#   species_specific_core.tab    (gene_name, species)
#   PAN/combined_DNA_CDS.fasta                  (headers: >clustering_id)
# gene_data_copy.tsv 
# OUTPUTS
#   /data/pam/team230/sm71/scratch/rp2/run_blast/queries/<species>/<gene>.fa
#       FASTA header: >gff_file, scaffold_name, cluster_id, species, gene_name
#   $BP/ssc_export_summary.tsv
#   $BP/ssc_missing_clusters.tsv   (gene_name, species, clustering_id not found in FASTA)
#
# NOTES
#   - Uses mawk where sensible; Python only to stream FASTA once & write files.
#   - Handles multiple genes mapping to the same cluster_id (writes same seq to
#     multiple gene files, as requested).
# -----------------------------------------------------------------------------

#BSUB -q normal
#BSUB -J ssc_fasta_export
#BSUB -n 1
#BSUB -M 4000
#BSUB -R "select[mem>4000] rusage[mem=4000]"
#BSUB -W 01:00
#BSUB -o /data/pam/team230/sm71/scratch/rp2/logs/ssc_fasta_export.%J.out
#BSUB -e /data/pam/team230/sm71/scratch/rp2/logs/ssc_fasta_export.%J.err

set -euo pipefail
export LC_ALL=C

# ---- paths ----
BP="/data/pam/team230/sm71/scratch/rp2/blast_preprocessing"             # reports saved here
PAN="/data/pam/team230/sm71/scratch/rp2/panaroo_output"
RUNROOT="/data/pam/team230/sm71/scratch/rp2/run_blast/queries"

ANNOT_TSV="$BP/species_specific_core_annotations.tab"
FASTA="$PAN/combined_DNA_CDS.fasta"

SUMMARY="$BP/ssc_export_summary.tsv"
MISSING="$BP/ssc_missing_clusters.tsv"

AWK_BIN="$(command -v mawk || command -v awk)"

echo "[INFO] Annotations : $ANNOT_TSV"
echo "[INFO] FASTA       : $FASTA"
echo "[INFO] Output root : $RUNROOT"
echo "[INFO] Reports to  : $BP"

# ---- sanity checks ----
[[ -s "$ANNOT_TSV" ]] || { echo "ERR: missing $ANNOT_TSV" >&2; exit 1; }
[[ -s "$FASTA"     ]] || { echo "ERR: missing $FASTA" >&2; exit 1; }
mkdir -p "$RUNROOT" "$BP"

# ---- build minimal map: cluster_id -> rows (gene_name, species, gff_file, scaffold_name)
# also emit species list
TMPDIR_LOCAL="$(mktemp -d)"
MAP_TSV="$TMPDIR_LOCAL/cluster_map.tsv"   # cluster_id \t gene_name \t species \t gff_file \t scaffold_name
SP_LIST="$TMPDIR_LOCAL/species.list"

echo "[STEP 1] Parse annotations → cluster map + species list"
$AWK_BIN -F'\t' -v OFS='\t' '
NR==1{
  for(i=1;i<=NF;i++){
    h=$i
    if(h=="gene_name")      gn=i
    else if(h=="species")   sp=i
    else if(h=="clustering_id") cid=i
    else if(h=="gff_file")  gf=i
    else if(h=="scaffold_name") sc=i
  }
  if(!gn||!sp||!cid||!gf||!sc){ print "ERR: missing required columns in annotations" > "/dev/stderr"; exit 2 }
  next
}
{
  print $cid, $gn, $sp, $gf, $sc
}
' "$ANNOT_TSV" > "$MAP_TSV"

cut -f3 "$MAP_TSV" | sort -u > "$SP_LIST"
echo "[STEP 1] done : species=$(wc -l < "$SP_LIST")  map_rows=$(wc -l < "$MAP_TSV")"

# ---- FASTA export (stream once; write per-gene files under species dirs)
echo "[STEP 2] Stream FASTA and write per-gene FASTAs per species"

python3 - "$MAP_TSV" "$FASTA" "$RUNROOT" "$MISSING" "$SUMMARY" << 'PY'
import sys, os, re
from collections import defaultdict

map_tsv, fasta, out_root, missing_path, summary_path = sys.argv[1:6]

# cluster_id -> list of (gene_name, species, gff_file, scaffold_name)
cid2rows = defaultdict(list)
species_set = set()

def sanitize_filename(s: str) -> str:
    # keep filename safe; collapse weird chars to underscore
    s = re.sub(r'[^\w.\-]+', '_', s.strip())
    s = re.sub(r'_{2,}', '_', s).strip('_')
    return s or "NA"

with open(map_tsv, 'r') as m:
    for line in m:
        cid, gn, sp, gff, scf = line.rstrip('\n').split('\t')
        cid2rows[cid].append((gn, sp, gff, scf))
        species_set.add(sp)

# Counters & tracking
written_per_species = defaultdict(int)
seen_cids = set()
missing_records = []  # (gene_name, species, cluster_id)

def emit_record(cid, seq):
    # write one record per (gene, species) mapping
    rows = cid2rows[cid]
    for gn, sp, gff, scf in rows:
        sp_dir = os.path.join(out_root, sp)
        os.makedirs(sp_dir, exist_ok=True)
        gene_file = os.path.join(sp_dir, f"{sanitize_filename(gn)}.fa")
        header = f">{gff}, {scf}, {cid}, {sp}, {gn}\n"
        with open(gene_file, 'a') as out:
            out.write(header)
            # wrap to 60 chars
            for i in range(0, len(seq), 60):
                out.write(seq[i:i+60] + "\n")
        written_per_species[sp] += 1

# FASTA streaming by cluster_id
cur_id = None
cur_seq = []

def flush():
    global cur_id, cur_seq
    if cur_id is None:
        return
    if cur_id in cid2rows:
        emit_record(cur_id, ''.join(cur_seq).replace('\n',''))
        seen_cids.add(cur_id)
    cur_id, cur_seq = None, []

def parse_cid(h):
    # header like: >clustering_id [optional rest]
    h = h[1:].strip()
    return h.split()[0] if h else ""

with open(fasta, 'r') as fin:
    for line in fin:
        if line.startswith('>'):
            flush()
            cur_id = parse_cid(line)
            cur_seq = []
        else:
            if cur_id is not None:
                cur_seq.append(line.strip())
    flush()

# Missing = any cid present in the map that we never saw in FASTA
for cid, rows in cid2rows.items():
    if cid not in seen_cids:
        for gn, sp, gff, scf in rows:
            missing_records.append((gn, sp, cid))

# Write missing report
with open(missing_path, 'w') as miss:
    miss.write("gene_name\tspecies\tclustering_id\n")
    for gn, sp, cid in missing_records:
        miss.write(f"{gn}\t{sp}\t{cid}\n")

# Summary
with open(summary_path, 'w') as summ:
    summ.write("metric\tvalue\n")
    total = sum(written_per_species.values())
    summ.write(f"seqs_written_total\t{total}\n")
    for sp in sorted(written_per_species):
        summ.write(f"seqs_written[{sp}]\t{written_per_species[sp]}\n")
    summ.write(f"missing_cluster_rows\t{len(missing_records)}\n")
PY

echo "[STEP 2] done."

echo "[DONE] Summary: $SUMMARY"
echo "[DONE] Missing : $MISSING"
echo "[DONE] FASTAs under: $RUNROOT"
rm -rf "$TMPDIR_LOCAL"
