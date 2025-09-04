#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# QC-filters partial-match (pm) cluster IDs and merges them with exact IDs to
# produce a clean per-SSC gene clustering_ids column.
#
# What it does:
# - Inputs:
#     1) species_specific_core_with_clustering_id.tab
#        (5 cols: gene_name, species, specific_class, clustering_ids, pm_clustering_ids)
#     2) species_specific_core.tab
#        (defines the SSC gene set: gene_name, species, specific_class)
#     3) gene_data.partial_matches.txt
#        (Stage 1 output: matched_target, gene_name(csv), clustering_id, dna_sequence)
#
# - QC rules on pm_clustering_ids (per gene *per species*):
#     (a) Drop pm CIDs whose mapped SSC gene occurs in >1 species
#         (not species-specific).
#     (b) Drop pm CIDs where the CSV gene name is itself an SSC gene
#         (indicating it’s not a true partial).
#
# - Merge step:
#     After pruning, merge exact 'clustering_ids' with pruned 'pm_clustering_ids'
#     → a single, de-duplicated, comma-joined clustering_ids column.
#     (clustering_ids = union(exact_cids, pruned_pm_cids))
#
# - Outputs:
#     * species_specific_core_with_clustering_id_qc.tab
#       (FINAL, 4 cols: gene_name, species, specific_class, clustering_ids)
#     * qc_prune_audit.txt
#       (per-gene audit: pm_before, pm_after, pruned_rule_a, pruned_rule_b)
#
# Usage:
#   bsub < qc_prune_pm_ids.sh
#
# Monitor (LSF):
#   tail -f logs/qc_prune_pm_ids.<JOBID>.out
#   tail -f logs/qc_prune_pm_ids.<JOBID>.err
# -----------------------------------------------------------------------------

#BSUB -q normal
#BSUB -J qc_prune_pm_ids
#BSUB -n 1
#BSUB -M 8000
#BSUB -R "select[mem>8000] rusage[mem=8000]"
#BSUB -o logs/qc_prune_pm_ids.%J.out
#BSUB -e logs/qc_prune_pm_ids.%J.err

set -euo pipefail
export LC_ALL=C
export PYTHONUNBUFFERED=1

# ---- paths ----
BLAST_DIR="/data/pam/team230/sm71/scratch/rp2/blast"
TMPDIR="/data/pam/team230/sm71/scratch/rp2/tmp"
FINAL_IN="$BLAST_DIR/species_specific_core_with_clustering_id.tab"      # 5 cols
SSC_TAB="$BLAST_DIR/species_specific_core.tab"                           # 3 cols
PARTIAL_LOG="$TMPDIR/gene_data.partial_matches.txt"                   # Stage 1 output
FINAL_OUT="$BLAST_DIR/species_specific_core_with_clustering_id_qc.tab"   # 4 cols after merge
AUDIT_OUT="$BLAST_DIR/qc_prune_audit.txt"

mkdir -p "$BLAST_DIR" logs
[[ -s "$FINAL_IN"    ]] || { echo "[ERR] missing $FINAL_IN"    >&2; exit 1; }
[[ -s "$SSC_TAB"     ]] || { echo "[ERR] missing $SSC_TAB"     >&2; exit 1; }
[[ -s "$PARTIAL_LOG" ]] || { echo "[ERR] missing $PARTIAL_LOG" >&2; exit 1; }

python3 - "$FINAL_IN" "$SSC_TAB" "$PARTIAL_LOG" "$FINAL_OUT" "$AUDIT_OUT" << 'PY'
import sys, csv

final_in, ssc_tab, partial_log, final_out, audit_out = sys.argv[1:]

def split_csv_field(s):
    s = (s or "").strip()
    if not s:
        return []
    return [x for x in (p.strip() for p in s.split(",")) if x]

# --- 1) SSC gene list + species multiplicity per gene ------------------------
ssc_genes = set()
gene_to_species = {}
with open(ssc_tab, newline="") as fin:
    fin.readline()  # header
    for line in fin:
        if not line.strip():
            continue
        parts = line.rstrip("\n").split("\t")
        parts += [""] * (3 - len(parts))
        gene, species = parts[0].strip(), parts[1].strip()
        if gene:
            ssc_genes.add(gene)
            gene_to_species.setdefault(gene, set()).add(species)

# --- 2) Exact maps: gene -> exact CIDs; CID -> genes (from FINAL_IN) ---------
exact_map = {}         # gene -> set(CIDs)
cluster_to_genes = {}  # CID  -> set(genes)
with open(final_in, newline="") as fin:
    fin.readline()  # header
    for line in fin:
        if not line.strip():
            continue
        parts = line.rstrip("\n").split("\t")
        parts += [""] * (5 - len(parts))
        gene, _, _, exact_cids, _ = parts[:5]
        gene = gene.strip()
        ex_set = set(split_csv_field(exact_cids))
        exact_map[gene] = ex_set
        for cid in ex_set:
            cluster_to_genes.setdefault(cid, set()).add(gene)

# --- 3) Rule (b) support via partial log -------------------------------------
# Drop pm CID for target T if the CSV gene name is an SSC gene.
# Partial log header is:
#   matched_target  gene_name(csv)  clustering_id  dna_sequence
prune_by_partial_name = {}  # matched_target -> set(CIDs to drop)
with open(partial_log, newline="") as fin:
    r = csv.reader(fin, delimiter="\t")
    header = next(r, None) or []
    col = {name: i for i, name in enumerate(header)}
    i_mt  = col.get("matched_target", 0)
    i_csv = col.get("gene_name(csv)", 1)
    i_cid = col.get("clustering_id", 2)

    for row in r:
        if not row:
            continue
        need = max(i_mt, i_csv, i_cid)
        if need >= len(row):
            row = row + [""] * (need + 1 - len(row))
        mt  = row[i_mt].strip()
        csv_name = row[i_csv].strip()
        cid = row[i_cid].strip()
        if cid and csv_name and (csv_name in ssc_genes):
            prune_by_partial_name.setdefault(mt, set()).add(cid)

# --- 4) QC filter + MERGE + audit (per gene *per species*) -------------------
audit_rows = []  # (gene, species, pm_before, pm_after, pruned_rule_a, pruned_rule_b)

with open(final_in, newline="") as fin, open(final_out, "w", newline="") as fout:
    # Write 4-col header for merged deliverable
    fout.write("\t".join(["gene_name","species","specific_class","clustering_ids"]) + "\n")

    # Read source rows
    fin_header = fin.readline()

    for line in fin:
        if not line.strip():
            continue
        parts = line.rstrip("\n").split("\t")
        parts += [""] * (5 - len(parts))
        gene, species, spec_class, exact_cids, pm_cids = [p.strip() for p in parts[:5]]

        exact_set = set(split_csv_field(exact_cids))
        pm_set = set(split_csv_field(pm_cids))
        pm_before = len(pm_set)
        pruned_a = set()
        pruned_b = set()

        if pm_set:
            # ---- Rule (a): prune pm CIDs whose mapped gene appears in >1 species
            for cid in list(pm_set):
                genes_for_cid = cluster_to_genes.get(cid, set())
                if any(len(gene_to_species.get(g, set())) > 1 for g in genes_for_cid):
                    pruned_a.add(cid)
            pm_set -= pruned_a

            # ---- Rule (b): prune pm CIDs where the CSV gene name is an SSC gene
            pruned_b = prune_by_partial_name.get(gene, set())
            pm_set -= pruned_b

        pm_after = len(pm_set)

        # ---- Merge exact + pruned partials into one clustering_ids column
        merged = sorted(exact_set | pm_set)
        merged_out = ",".join(merged)

        # 4-column output row
        fout.write("\t".join([gene, species, spec_class, merged_out]) + "\n")

        audit_rows.append((gene, species, pm_before, pm_after, len(pruned_a), len(pruned_b)))

# --- 5) Audit report ----------------------------------------------------------
with open(audit_out, "w", newline="") as aout:
    aout.write("gene_name\tspecies\tpm_before\tpm_after\tpruned_rule_a\tpruned_rule_b\n")
    for row in audit_rows:
        aout.write("\t".join(map(str, row)) + "\n")
PY

echo "[QC] Done."
echo "[QC] Output written -> $FINAL_OUT"
echo "[QC] Audit report   -> $AUDIT_OUT"
head -5 "$AUDIT_OUT" || true
