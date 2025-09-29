#!/usr/bin/env bash
set -euo pipefail

# Build a stable ground-truth mapping: <reference_id>\t<species>
# from subset_18.tsv produced by get_refs.sh
#
# Columns expected in subset_18.tsv:
#   1: species (normalized)
#   2: raw taxonomy
#   3: file path (inside tar)
#   4: tarball name
#   5: URL
#
# Usage:
#   ./make_ref_ground_truth.sh [--in subset_18.tsv] [--out ref_ground_truth.tsv] [--force]
#
# Notes:
# - Output is sorted (LC_ALL=C) by reference_id and deduplicated for reproducibility.
# - We strip known FASTA extensions: .fa, .fasta, .fna and optional .gz
# - We DO NOT strip text after underscores in the ID (keeps full basename).
#   If you need pre-underscore IDs, see the optional flag below.

export LC_ALL=C

# ---------- CLI ----------
IN="subset_18.tsv"
OUT="ref_ground_truth.tsv"
FORCE=0
STRIP_UNDERSCORE=0   # optional: turn on to mimic comm/_wanted ids behavior

usage() {
  cat <<EOF
Usage: $0 [--in FILE] [--out FILE] [--force] [--strip-after-underscore]

Options:
  --in FILE                    Input TSV (default: subset_18.tsv)
  --out FILE                   Output TSV (default: ref_ground_truth.tsv)
  --force                      Overwrite OUT if it already exists
  --strip-after-underscore     Trim ID at first underscore (e.g. GCF123_foo -> GCF123)
EOF
}

while (( $# )); do
  case "$1" in
    --in) IN="$2"; shift 2;;
    --out) OUT="$2"; shift 2;;
    --force) FORCE=1; shift;;
    --strip-after-underscore) STRIP_UNDERSCORE=1; shift;;
    -h|--help) usage; exit 0;;
    *) echo "[FATAL] unknown arg: $1" >&2; usage; exit 1;;
  endesac
done

# ---------- checks ----------
[[ -f "$IN" ]] || { echo "[FATAL] input not found: $IN" >&2; exit 1; }
if [[ -s "$OUT" && $FORCE -eq 0 ]]; then
  echo "[SKIP] $OUT exists (use --force to overwrite)"; exit 0
fi

# ---------- build ----------
# Normalize CRLF just in case, extract final path component, strip extensions,
# optional underscore-trim, then emit: <id>\t<species>, sort unique.
tmp="$(mktemp)"
trap 'rm -f "$tmp"' EXIT

# Use awk to parse and clean
# shellcheck disable=SC2016
awk -F'\t' -v strip_us="$STRIP_UNDERSCORE" '
  {
    gsub(/\r$/, "", $0)                    # normalize CRLF
    species = $1
    p = $3
    # take last path component
    gsub(/^.*\//, "", p)
    # strip common fasta extensions (optionally with .gz)
    sub(/\.(fa|fasta|fna)(\.gz)?$/, "", p)
    if (strip_us == 1) { sub(/_.*/, "", p) }  # optional
    if (p != "" && species != "") {
      print p "\t" species
    }
  }
' "$IN" > "$tmp"

# stable, unique output
sort -u -t$'\t' -k1,1 "$tmp" > "$OUT"

echo "[OK] wrote $OUT (rows: $(wc -l < "$OUT"))"
