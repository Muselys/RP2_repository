#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# Extracts SSC gene sequences from Panaroo’s combined_DNA_CDS.fasta using
# clustering_ids from the QC table, writing one FASTA per (species, gene).
#
# IMPORTANT ORDER:
#
# What it does:
# 1) Read species_specific_core_with_clustering_id_qc.tab (4 cols) and map:
#      cluster_id → [(out_path=OUT_DIR/species_gene.fa, header_prefix=species_gene_)]
# 2) Stream combined_DNA_CDS.fasta once (RAM-safe):
#      - Exact match: first FASTA header token == cluster_id → write record.
#      - Fallback: if no exact match, search header for any cluster_id substring.
#      - Write FASTA with header “>species_gene_<cluster_id>”, wrap sequence at 60 chars.
# 3) Sanitize filenames (keep letters/numbers/.-_/; replace others with “_”).
# 4) Audit any cluster_ids not found in the FASTA.
#
# Inputs:
#   - species_specific_core_with_clustering_id_qc.tab
#   - combined_DNA_CDS.fasta  (Panaroo)
#
# Outputs:
#   - OUT_DIR/species_gene.fa  (per gene; may contain multiple cluster IDs)
#   - extract_sequences_missing_ids.txt  (cluster_ids not found)
#
# Notes:
#   - Appends to existing per-gene FASTAs (safe for re-runs).
#   - Case-sensitive matching on cluster_ids.
#   - Uses a union regex for up to ~50k IDs; falls back to substring scan if needed.
#   - Prints summary stats (rows parsed, files written, top files by record count).
# -----------------------------------------------------------------------------

#BSUB -q normal
#BSUB -J extract_gene_fastas
#BSUB -n 1
#BSUB -M 64000
#BSUB -R "select[mem>64000] rusage[mem=64000]"
#BSUB -o /data/pam/team230/sm71/scratch/rp2/blast/ssc_gene_fastas/logs/extract_gene_fastas.%J.out
#BSUB -e /data/pam/team230/sm71/scratch/rp2/blast/ssc_gene_fastas/logs/extract_gene_fastas.%J.err

set -euo pipefail
export LC_ALL=C
export PYTHONUNBUFFERED=1

# ---- paths ----
BLAST_DIR="/data/pam/team230/sm71/scratch/rp2/blast"
QC_TAB="$BLAST_DIR/species_specific_core_with_clustering_id_qc.tab"  # 4-col input
PANAROO_FASTA="/data/pam/team230/sm71/scratch/rp2/panaroo_output/combined_DNA_CDS.fasta"

OUT_DIR="$BLAST_DIR/ssc_gene_fastas"
MISSING_AUDIT="$BLAST_DIR/extract_sequences_missing_ids.txt"

mkdir -p "$BLAST_DIR" "$OUT_DIR" logs
[[ -s "$QC_TAB"        ]] || { echo "[ERR] missing $QC_TAB"        >&2; exit 1; }
[[ -s "$PANAROO_FASTA" ]] || { echo "[ERR] missing $PANAROO_FASTA" >&2; exit 1; }

python3 - "$QC_TAB" "$PANAROO_FASTA" "$OUT_DIR" "$MISSING_AUDIT" << 'PY'
import sys, os, re

qc_tab, fasta_path, out_dir, missing_audit = sys.argv[1:]

def split_ids(s):
    s = (s or "").strip()
    if not s:
        return []
    return [x for x in (p.strip() for p in s.split(",")) if x]

def sanitize_filename(s):
    # Keep letters, numbers, .-_ ; replace others with _
    return re.sub(r'[^A-Za-z0-9.\-_/]+', '_', s)

def main():
    # --- 1) Build mapping: cluster_id -> list of output files (species_gene.fa) ---
    cid_to_targets = {}     # cid -> list of (outfile_path, species_gene_prefix)
    file_to_cids   = {}     # outfile_path -> set(cids)
    rows = 0
    assigned_ids = 0

    with open(qc_tab) as fin:
        header_cols = fin.readline().rstrip("\n").split("\t")
        if len(header_cols) < 4:
            raise SystemExit(f"[ERR] {qc_tab} must have 4 columns: gene_name, species, specific_class, clustering_ids")
        col = {name:i for i,name in enumerate(header_cols)}
        i_gene = col.get("gene_name", 0)
        i_species = col.get("species", 1)
        i_spec_class = col.get("specific_class", 2)
        i_cids = col.get("clustering_ids", 3)

        for line in fin:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            parts += [""] * (max(i_gene, i_species, i_spec_class, i_cids) + 1 - len(parts))
            gene    = parts[i_gene].strip()
            species = parts[i_species].strip()
            cids    = split_ids(parts[i_cids])
            rows += 1
            if not gene or not species or not cids:
                continue

            base = f"{species}_{gene}.fa"
            out_path = os.path.join(out_dir, sanitize_filename(base))
            file_to_cids.setdefault(out_path, set()).update(cids)

            prefix = f"{species}_{gene}_"
            for cid in set(cids):
                cid_to_targets.setdefault(cid, []).append((out_path, prefix))
                assigned_ids += 1

    print(f"[INFO] parsed {rows} rows; {len(cid_to_targets)} unique cluster IDs to fetch; {assigned_ids} total (with duplicates per gene).", file=sys.stderr)

    # --- 2) Stream the FASTA and write sequences as we find them ------------------
    found = {}  # cid -> True if written at least once
    written_counts = {}  # outfile_path -> number of records written

    def write_record(cid, header_line, seq, targets):
        for out_path, prefix in targets:
            os.makedirs(os.path.dirname(out_path), exist_ok=True)
            with open(out_path, "a") as fh:
                fh.write(f">{prefix}{cid}\n")
                for i in range(0, len(seq), 60):
                    fh.write(seq[i:i+60] + "\n")
            written_counts[out_path] = written_counts.get(out_path, 0) + 1

    # locals captured by flush_record via nonlocal
    header = None
    seq_chunks = []
    total_records = 0

    # Optional regex union for substring fallback
    cid_list = list(cid_to_targets.keys())
    union_rx = None
    if 0 < len(cid_list) <= 50000:
        try:
            union_rx = re.compile("|".join(re.escape(x) for x in cid_list))
        except Exception:
            union_rx = None

    def flush_record():
        nonlocal header, seq_chunks, total_records
        if header is None:
            return
        total_records += 1
        seq = "".join(seq_chunks).replace(" ", "").replace("\r", "").replace("\n", "")
        if not seq:
            header = None
            seq_chunks = []
            return

        first_token = header[1:].split(None, 1)[0] if header.startswith(">") else header.split(None,1)[0]
        targets = cid_to_targets.get(first_token)
        if targets is not None:
            write_record(first_token, header, seq, targets)
            found[first_token] = True
        else:
            # substring fallback
            m_cids = []
            if union_rx:
                for m in union_rx.finditer(header):
                    m_cids.append(m.group(0))
            else:
                for cid in cid_to_targets.keys():
                    if cid in header:
                        m_cids.append(cid)
            if m_cids:
                for cid in set(m_cids):
                    write_record(cid, header, seq, cid_to_targets[cid])
                    found[cid] = True

        header = None
        seq_chunks = []

    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                flush_record()
                header = line.rstrip("\n")
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
        flush_record()

    print(f"[INFO] streamed FASTA; wrote to {len(written_counts)} files; matched {len(found)} unique cluster IDs.", file=sys.stderr)

    # --- 3) Report missing cluster IDs -------------------------------------------
    missing = sorted([cid for cid in cid_to_targets.keys() if cid not in found])
    with open(missing_audit, "w") as aout:
        aout.write("# cluster_ids_not_found_in_combined_DNA_CDS.fasta\n")
        for cid in missing:
            aout.write(cid + "\n")

    print(f"[INFO] missing IDs: {len(missing)} -> {missing_audit}", file=sys.stderr)
    top = sorted(written_counts.items(), key=lambda x: x[1], reverse=True)[:5]
    for p,c in top:
        print(f"[INFO] top file: {os.path.basename(p)} -> {c} seqs", file=sys.stderr)

if __name__ == "__main__":
    main()

PY

echo "[EXTRACT] Done."
echo "[EXTRACT] Output FASTAs dir -> $OUT_DIR"
echo "[EXTRACT] Missing IDs audit -> $MISSING_AUDIT"
ls -1 "$OUT_DIR" | head -10 || true


