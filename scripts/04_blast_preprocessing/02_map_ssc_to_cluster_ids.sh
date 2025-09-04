#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# Extracts SSC targets and maps them to clustering IDs using a single-pass,
# RAM-safe stream over gene_data.csv. Outputs a per-target summary and a final
# SSC table joined with exact + partial clustering IDs.
#
# Inputs:
#   BLAST_DIR/species_specific_core.tab
#   BLAST_DIR/gene_data.csv
#
# Outputs:
#   TMPDIR/_targets.list
#   TMPDIR/_ssc_exact_map.tsv
#   TMPDIR/gene_data.partial_matches.txt
#   TMPDIR/_exact.counts, TMPDIR/_partial.counts, TMPDIR/_partial.names (Stage 2)
#   BLAST_DIR/ssc_match_summary.tsv
#   TMPDIR/_clusters.agg
#   TMPDIR/_pm_agg.tsv
#   BLAST_DIR/species_specific_core_with_clustering_id.tab  (FINAL)
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# WHAT THIS JOB DOES (step-by-step)
# 1) Collect SSC targets (gene_name col 1) from species_specific_core.tab.
#    -> writes: _targets.list
# 2) Stream gene_data.csv ONCE for ALL targets (RAM-safe):
#    - EXACT hits (CSV gene == SSC target)  -> _ssc_exact_map.tsv
#      Columns: matched_target  gene_name(csv)  clustering_id  dna_sequence
#    - PARTIAL hits (target substring appears anywhere in CSV gene, but not exact)
#      -> gene_data.partial_matches.txt (same 4 columns)
#    - Progress printed every PROG_EVERY rows; final DONE summary.
# 3) Build per-gene summary (Stage 2):
#    - exact_count, partial_count, comma-joined partial_names per target
#    -> ssc_match_summary.tsv
#    - QA audits: exact/partial hits with empty clustering_id
# 4) Aggregate:
#    - exact clustering_ids per gene (comma-joined, empty CIDs skipped)
#      -> _clusters.agg
#    - partial clustering_ids per gene (comma-joined, empty CIDs skipped)
#      -> _pm_agg.tsv
# 5) Left-join back onto SSC to produce final deliverable:
#    species_specific_core_with_clustering_id.tab
#    Columns (with header):
#      gene_name  species  specific_class  clustering_ids  pm_clustering_ids
# -----------------------------------------------------------------------------

#BSUB -q normal
#BSUB -J ssc_map_clusters
#BSUB -n 1
#BSUB -M 32000
#BSUB -R "select[mem>32000] rusage[mem=32000]"
#BSUB -o /data/pam/team230/sm71/scratch/rp2/blast/logs/ssc_map_clusters.%J.out
#BSUB -e /data/pam/team230/sm71/scratch/rp2/blast/logs/ssc_map_clusters.%J.err

set -euo pipefail
export LC_ALL=C
export PYTHONUNBUFFERED=1

# Optional modules (don’t fail if absent)
if command -v module &>/dev/null; then
  module purge || true
  module load python/3.12.0 2>/dev/null || module load python/3.10 2>/dev/null || true
fi
PYTHON="$(command -v python3)"
echo "[ENV] Using $($PYTHON --version) @ $PYTHON"
which "$PYTHON" || true

# ---- config ----
export PROG_EVERY=${PROG_EVERY:-200000}   # print progress every N CSV rows

# ---- paths ----
BLAST_DIR="/data/pam/team230/sm71/scratch/rp2/blast"
GD_CSV="$BLAST_DIR/gene_data.csv"                       # huge CSV
SSC_TAB="$BLAST_DIR/species_specific_core.tab"          # 3-col: gene_name  species  specific_class
OUT_JOIN="$BLAST_DIR/species_specific_core_with_clustering_id.tab"

# ---- temp dirs ----
TMPDIR="/data/pam/team230/sm71/scratch/rp2/tmp"
mkdir -p "$TMPDIR"/{stage0,stage1,stage2,stage3,stage4} "$BLAST_DIR" logs
export TMPDIR

# ---- temps + side outputs (kept, organized) ----

# Stage 0: target list
TARGETS_LIST="$TMPDIR/stage0/_targets.list"

# Stage 1: streaming results
MAP_TMP="$TMPDIR/stage1/_ssc_exact_map.tsv"                 # matched_target  gene_name(csv)  clustering_id  dna_sequence
PARTIAL_LOG="$TMPDIR/stage1/gene_data.partial_matches.txt"  # same 4 cols

# Stage 2: summary + audits
EXACT_COUNT="$TMPDIR/stage2/_exact.counts"
PARTIAL_COUNT="$TMPDIR/stage2/_partial.counts"
PARTIAL_PAIRS="$TMPDIR/stage2/_partial.pairs"      # matched_target \t gene_name(csv)
PARTIAL_AGG="$TMPDIR/stage2/_partial.names"        # matched_target \t name1,name2,...
EXACT_EMPTY="$TMPDIR/stage2/_exact_empty_cid.tsv"
PM_EMPTY="$TMPDIR/stage2/_partial_empty_cid.tsv"
SUMMARY="$BLAST_DIR/ssc_match_summary.tsv"         # final summary stays in BLAST_DIR

# Stage 3: aggregation
CLUSTERS_PAIRS="$TMPDIR/stage3/_clusters.pairs"
CLUSTERS_AGG="$TMPDIR/stage3/_clusters.agg"
PM_PAIRS="$TMPDIR/stage3/_pm_pairs.tsv"
PM_AGG="$TMPDIR/stage3/_pm_agg.tsv"

# Stage 4 doesn’t really need temps, except maybe the intermediate join
TMP_JOIN="$TMPDIR/stage4/_ssc_with_cids.tmp"

# input checks
[[ -s "$GD_CSV" ]] || { echo "[ERR] missing $GD_CSV" >&2; exit 1; }
[[ -s "$SSC_TAB" ]] || { echo "[ERR] missing $SSC_TAB" >&2; exit 1; }

# -----------------------------
# Stage 0: collect targets
# -----------------------------
echo "[STAGE 0] $(date +%T) Collecting target gene names from SSC…"
cut -f1 "$SSC_TAB" | sed '1d' | sort -u > "$TARGETS_LIST"
echo "[STAGE 0] $(date +%T) Targets: $(wc -l < "$TARGETS_LIST")"

# -----------------------------
# Stage 1: single-pass streamer (ALL targets) — compiled regex
# -----------------------------
echo "[STAGE 1] $(date +%T) Streaming gene_data.csv once for ALL targets (exact + partial)…"
$PYTHON - "$GD_CSV" "$TARGETS_LIST" "$MAP_TMP" "$PARTIAL_LOG" "$PROG_EVERY" << 'PY'
import sys, csv, time, re

gd_csv, targets_list, map_tmp, partial_log, prog_every = sys.argv[1:]
prog_every = int(prog_every)

# Load targets
targets = []
with open(targets_list, 'r', newline='') as f:
    for line in f:
        g = line.strip()
        if g:
            targets.append(g)

# Fast exact check
targets_set = set(targets)

# Build compiled regexes for partials:
# - Escape every target (literal match)
# - Sort by length desc (helps avoid pathological alternation)
# - Chunk to keep patterns manageable (e.g., 400–800 per chunk)
def build_regex_chunks(items, chunk_size=600):
    esc = [re.escape(x) for x in items]
    esc.sort(key=len, reverse=True)
    chunks = []
    for i in range(0, len(esc), chunk_size):
        pat = "|".join(esc[i:i+chunk_size])
        # Named group 't' so we can pull the matched target easily
        chunks.append(re.compile(f"(?P<t>{pat})"))
    return chunks

regex_chunks = build_regex_chunks(targets)

print(f"[INFO] partial regex built in {len(regex_chunks)} chunk(s) covering {len(targets)} targets.",
      file=sys.stderr, flush=True)

# CSV prep
try:
    csv.field_size_limit(10**9)
except Exception:
    pass

rows = exact = partial = 0
t0 = time.time()

with open(map_tmp, "w", newline="") as m, \
     open(partial_log, "w", newline="") as p, \
     open(gd_csv, newline="") as fin:

    exact_w = csv.writer(m, delimiter="\t")
    part_w  = csv.writer(p, delimiter="\t")

    # unified 4-col headers
    exact_w.writerow(["matched_target","gene_name(csv)","clustering_id","dna_sequence"])
    part_w.writerow (["matched_target","gene_name(csv)","clustering_id","dna_sequence"])

    r = csv.reader(fin)
    header = next(r)
    col = {name:i for i,name in enumerate(header)}
    for need in ("gene_name","clustering_id","dna_sequence"):
        if need not in col:
            print(f"[ERROR] Missing column '{need}' in gene_data.csv", file=sys.stderr, flush=True)
            sys.exit(2)
    max_idx = max(col.values())

    for row in r:
        rows += 1
        if len(row) <= max_idx:
            continue

        gn  = row[col["gene_name"]]     or ""
        cid = row[col["clustering_id"]] or ""
        dna = row[col["dna_sequence"]]  or ""

        # EXACT: perfect name match (fast set membership)
        if gn in targets_set:
            # matched_target == gn for exacts
            exact_w.writerow([gn, gn, cid, dna])
            exact += 1
        else:
            # PARTIALS: run compiled regex chunks; collect all matching targets in this gn
            hits = set()
            for rx in regex_chunks:
                for m in rx.finditer(gn):
                    t = m.group('t')
                    # keep only true partials: 'gn' is NOT exactly the target
                    if t != gn:
                        hits.add(t)
                # small optimization: if we found lots already, keep going but avoid duplicates
            if hits:
                for t in sorted(hits):
                    part_w.writerow([t, gn, cid, dna])
                    partial += 1

        if rows % prog_every == 0:
            elapsed = time.time()-t0
            rate = rows/elapsed if elapsed>0 else 0
            print(f"[PROG] {rows:,} rows | exact:{exact:,} partial:{partial:,} | {rate:,.0f} rows/s | {elapsed/60:.1f} min",
                  file=sys.stderr, flush=True)

elapsed = time.time()-t0
rate = rows/elapsed if elapsed>0 else 0
print(f"[DONE] Streamed {rows:,} rows | exact:{exact:,} partial:{partial:,} | {rate:,.0f} rows/s | {elapsed/60:.1f} min",
      file=sys.stderr, flush=True)
PY
echo "[STAGE 1] $(date +%T) Done. Exact map: $MAP_TMP ; Partial audit: $PARTIAL_LOG"
head -3 "$MAP_TMP" || true
head -5 "$PARTIAL_LOG" || true


# -----------------------------
# Stage 2: build per-gene summary
# -----------------------------
echo "[STAGE 2] $(date +%T) Building summary (exact/partial counts + partial names)…"

# exact counts per matched_target (from _ssc_exact_map.tsv col1)
awk -F'\t' 'NR>1{c[$1]++} END{for(k in c) print k"\t"c[k]}' "$MAP_TMP" \
  | sort -t$'\t' -k1,1 > "$EXACT_COUNT"

# partial counts per matched_target (from partial log col1)
awk -F'\t' 'NR>1{c[$1]++} END{for(k in c) print k"\t"c[k]}' "$PARTIAL_LOG" \
  | sort -t$'\t' -k1,1 > "$PARTIAL_COUNT"

# unique partial names per target (matched_target \t gene_name(csv))
cut -f1,2 "$PARTIAL_LOG" | sed '1d' \
  | sort -t$'\t' -k1,1 -k2,2 -u > "$PARTIAL_PAIRS"

# aggregate partial names into comma list per target
awk -F'\t' '{a[$1]=(a[$1]?a[$1]","$2:$2)} END{for(k in a) print k"\t"a[k]}' \
  "$PARTIAL_PAIRS" | sort -t$'\t' -k1,1 > "$PARTIAL_AGG"

# compose summary over ALL targets; fill missing with 0 / empty
printf "matched_target\texact_count\tpartial_count\tpartial_names\n" > "$SUMMARY"
awk -F'\t' '
  FNR==NR { ec[$1]=$2; next }                           # EXACT_COUNT
  FILENAME==ARGV[2] { pc[$1]=$2; next }                 # PARTIAL_COUNT
  FILENAME==ARGV[3] { pn[$1]=$2; next }                 # PARTIAL_AGG
  {
    g=$0;
    print g "\t" (g in ec?ec[g]:0) "\t" (g in pc?pc[g]:0) "\t" (g in pn?pn[g]:"")
  }' \
  "$EXACT_COUNT" "$PARTIAL_COUNT" "$PARTIAL_AGG" \
  "$TARGETS_LIST" >> "$SUMMARY"

echo "[STAGE 2] $(date +%T) Summary written -> $SUMMARY"
head -5 "$SUMMARY" || true

# QA audits of empty clustering_id (exact/partial)
awk -F'\t' 'NR>1 && $3=="" {print $1"\t"$2}' "$MAP_TMP" > "$EXACT_EMPTY" || true
awk -F'\t' 'NR>1 && $3=="" {print $1"\t"$2}' "$PARTIAL_LOG" > "$PM_EMPTY" || true
echo "[AUDIT] exact hits with empty clustering_id:   $(wc -l < "$EXACT_EMPTY") -> $EXACT_EMPTY"
echo "[AUDIT] partial hits with empty clustering_id: $(wc -l < "$PM_EMPTY")   -> $PM_EMPTY"


# -----------------------------
# Stage 3: aggregate exact IDs
# -----------------------------
echo "[STAGE 3] $(date +%T) Aggregating exact cluster IDs per gene…"
# unique gene->cluster pairs (skip empty clustering_id values), using cols 1,3
cut -f1,3 "$MAP_TMP" | sed '1d' | awk -F'\t' 'NF==2 && $2!=""' \
  | sort -t$'\t' -k1,1 -k2,2 -u > "$CLUSTERS_PAIRS"

# aggregate to comma list per matched_target
awk -F'\t' '{a[$1]=(a[$1]?a[$1]","$2:$2)} END{for(k in a) print k"\t"a[k]}' \
  "$CLUSTERS_PAIRS" > "$CLUSTERS_AGG"

echo "[STAGE 3] $(date +%T) Exact aggregation done: $(wc -l < "$CLUSTERS_AGG") genes"
head -5 "$CLUSTERS_AGG" || true

# -----------------------------
# Stage 3.5: aggregate partial IDs
# -----------------------------
echo "[STAGE 3.5] $(date +%T) Aggregating partial-match cluster IDs per gene…"
# matched_target \t clustering_id from partial log (cols 1,3); drop empty CIDs
awk -F'\t' 'NR>1 && $3!="" {print $1"\t"$3}' "$PARTIAL_LOG" \
  | sort -t$'\t' -k1,1 -k2,2 -u > "$PM_PAIRS"

awk -F'\t' '{a[$1]=(a[$1]?a[$1]","$2:$2)} END{for(k in a) print k"\t"a[k]}' \
  "$PM_PAIRS" > "$PM_AGG"

echo "[STAGE 3.5] $(date +%T) Partial aggregation done: $(wc -l < "$PM_AGG") genes"
head -5 "$PM_AGG" || true

# -----------------------------
# Stage 4: join onto SSC
# -----------------------------
echo "[STAGE 4] $(date +%T) Joining clusters + partial clusters onto SSC…"


# header
printf "gene_name\tspecies\tspecific_class\tclustering_ids\tpm_clustering_ids\n" > "$OUT_JOIN"

# left-join exact clustering_ids
join -t $'\t' -a1 -e "" -o 1.1,1.2,1.3,2.2 \
  -1 1 -2 1 \
  <(sed '1d' "$SSC_TAB" | sort -t$'\t' -k1,1) \
  <(sort -t$'\t' -k1,1 "$CLUSTERS_AGG") \
  > "$TMP_JOIN"

# left-join pm_clustering_ids onto previous result
join -t $'\t' -a1 -e "" -o 1.1,1.2,1.3,1.4,2.2 \
  -1 1 -2 1 \
  <(sort -t$'\t' -k1,1 "$TMP_JOIN") \
  <(sort -t$'\t' -k1,1 "$PM_AGG") \
  >> "$OUT_JOIN"

echo "[STAGE 4] $(date +%T) Wrote -> $OUT_JOIN"
head -5 "$OUT_JOIN" || true
echo -n "[CHECK] final rows (incl header): "; wc -l < "$OUT_JOIN"
echo -n "[CHECK] SSC rows (incl header): "; wc -l < "$SSC_TAB"

echo "[DONE] $(date +%T) All stages complete."

# Optional cleanup (keep commented if you want to reuse intermediates)
# echo "[CLEANUP] Removing temp files in $TMPDIR"
# rm -f "$TMPDIR"/_*

