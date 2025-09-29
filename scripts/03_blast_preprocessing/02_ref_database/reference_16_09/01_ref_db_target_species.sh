#!/usr/bin/env bash
set -euo pipefail

# =========================================
# Usage: ./get_refs_from_zero.sh [--force]
# - Idempotent: skips steps if outputs exist
# - Use --force to rebuild artifacts
# Env overrides:
#   ATB=/path/to/workdir
#   REFDIR=/path/to/workdir/ref/fasta
#   TARDIR=/path/to/workdir/tarballs
#   JOBS=6
#
#bsub -q normal -J get_refs -n 8 \
#     -R "select[mem>64000] rusage[mem=64000]" -M 64000 \
#     -o /data/pam/team230/sm71/scratch/rp2/logs/get_refs.%J.out \
#     -e /data/pam/team230/sm71/scratch/rp2/logs/get_refs.%J.err \
#     bash get_refs_from_atb.sh
#
# =========================================

FORCE=0
if [[ "${1:-}" == "--force" ]]; then FORCE=1; fi

# ========== 0) setup ==========
ATB="${ATB:-/data/pam/team230/sm71/scratch/rp2/atb}"
REFDIR="${REFDIR:-$ATB/ref/fasta}"
TARDIR="${TARDIR:-$ATB/tarballs}"
JOBS="${JOBS:-6}"

mkdir -p "$REFDIR" "$TARDIR"
cd "$ATB"
export LC_ALL=C
alias zcat='zcat'  # keep explicit

# tools
command -v wget >/dev/null 2>&1 || { echo "[FATAL] wget not found"; exit 1; }
command -v tar  >/dev/null 2>&1 || { echo "[FATAL] tar not found";  exit 1; }
if ! command -v parallel >/dev/null 2>&1; then
  echo "[WARN] GNU parallel not found; will extract serially"
  USE_PARALLEL=0
else
  USE_PARALLEL=1
fi

# helper: run step if file missing or --force
need() { [[ $FORCE -eq 1 || ! -s "$1" ]]; }

# ========== 1) targets ==========
if need "target_species.tsv"; then
  cat > target_species.tsv <<'EOF'
Staphylococcus aureus
Staphylococcus epidermidis
Staphylococcus pseudintermedius
Staphylococcus haemolyticus
Staphylococcus capitis
Staphylococcus sciuri
Staphylococcus argenteus
Streptococcus pneumoniae
Streptococcus pyogenes
Streptococcus agalactiae
Streptococcus suis
Streptococcus equi
Streptococcus dysgalactiae
Streptococcus uberis
Streptococcus mutans
Streptococcus mitis
Enterococcus faecium
Enterococcus faecalis
EOF
  echo "[OK] wrote target_species.tsv"
else
  echo "[SKIP] target_species.tsv"
fi

if need "target_species.norm.tsv"; then
  awk '{print tolower($0)}' target_species.tsv > target_species.norm.tsv
  echo "[OK] wrote target_species.norm.tsv"
else
  echo "[SKIP] target_species.norm.tsv"
fi

# ========== 2) master list ==========
if need "file_list.all.latest.tsv.gz"; then
  wget -q https://osf.io/download/4yv85/ -O file_list.all.latest.tsv.gz
  echo "[OK] downloaded file_list.all.latest.tsv.gz"
else
  echo "[SKIP] file_list.all.latest.tsv.gz"
fi
echo "[INFO] AllTheBacteria rows: $(zcat file_list.all.latest.tsv.gz | tail -n +2 | wc -l)"

# ========== 3) filter to the 18 ==========
if need "subset_18.tsv"; then
  zcat file_list.all.latest.tsv.gz \
  | awk -F'\t' '
    BEGIN{
      while ((getline l < "target_species.norm.tsv") > 0) if(l!="") T[l]=1
    }
    NR==1{next}
    {
      raw=$2; fn=$4; tar=$5; url=$6
      n=split(raw, parts, / *; */)
      keep=0; tag=""
      for(i=1;i<=n && !keep;i++){
        sp=tolower(parts[i])
        split(sp,a,/ +/)
        if(!(a[1] && a[2])) continue
        g=a[1]; s=a[2]
        sub(/_.*/,"",g)   # Enterococcus_B -> Enterococcus
        sub(/_.*/, "",s)  # faecium_A -> faecium
        if(g=="mammaliicoccus" && s=="sciuri") norm="staphylococcus sciuri"
        else norm=g" "s
        if(norm in T){ keep=1; tag=norm }
      }
      if(keep) print tag"\t"$2"\t"fn"\t"tar"\t"url
    }
  ' > subset_18.tsv
  echo "[OK] wrote subset_18.tsv"
else
  echo "[SKIP] subset_18.tsv"
fi

# summary (cheap, always refresh)
{
  echo "Targets present: $(cut -f1 subset_18.tsv | sort -u | wc -l) of 18"
  cut -f1 subset_18.tsv | sort | uniq -c | sort -nr | sed "s/^ *//"
} > target_counts.txt
echo "[OK] wrote target_counts.txt"

# ========== 4) plan files ==========
if need "tar2url.tsv"; then
  awk -F'\t' '{print $4"\t"$5}' subset_18.tsv | sort -u > tar2url.tsv
  echo "[OK] wrote tar2url.tsv"
else
  echo "[SKIP] tar2url.tsv"
fi

if need "files_by_tar.tsv"; then
  awk -F'\t' '{print $4"\t"$3}' subset_18.tsv | sort > files_by_tar.tsv
  echo "[OK] wrote files_by_tar.tsv"
else
  echo "[SKIP] files_by_tar.tsv"
fi

# ========== 5) desired IDs (SampleID only) ==========
if need "wanted.ids.txt"; then
  cut -f2 files_by_tar.tsv \
  | awk -F'/' '{print $NF}' \
  | sed 's/\.[^.]*$//' | sed 's/_.*$//' \
  | sort -u > wanted.ids.txt
  echo "[OK] wrote wanted.ids.txt"
else
  echo "[SKIP] wanted.ids.txt"
fi
echo "[INFO] total target IDs: $(wc -l < wanted.ids.txt)"

# current have (recompute every run)
find "$REFDIR" -type f \( -name '*.fa' -o -name '*.fa.gz' \) -printf '%f\n' \
| sed 's/\.[^.]*$//' | sed 's/_.*$//' | sort -u > have.ids.txt

comm -23 wanted.ids.txt have.ids.txt > missing.ids.txt
echo "[INFO] missing IDs: $(wc -l < missing.ids.txt)"

# rows to extract (only missing)
awk -v OFS='\t' '
  NR==FNR {miss[$1]=1; next}
  {
    bn=$2
    sub(/^.*\//,"",bn)
    sub(/\.[^.]*$/,"",bn)
    sub(/_.*/, "", bn)
    if (bn in miss) print $1, $2
  }
' missing.ids.txt files_by_tar.tsv > files_by_tar.todo.tsv

cut -f1 files_by_tar.todo.tsv | sort -u > tar.todo.list
echo "[INFO] TODO members: $(wc -l < files_by_tar.todo.tsv)"
echo "[INFO] TODO tarballs: $(wc -l < tar.todo.list)"

# bail early if nothing to do
if [[ ! -s files_by_tar.todo.tsv ]]; then
  echo "[DONE] nothing to do; all target IDs already present in $REFDIR"
  exit 0
fi

# ========== 6) map URLs + download into TARDIR ==========
sed -i 's/\r$//' tar.todo.list tar2url.tsv  # normalize CRLF

# refresh tar2url from master if missing/forced (extra safety)
if need "tar2url.tsv"; then
  zcat file_list.all.latest.tsv.gz | awk -F'\t' 'NR>1{print $5"\t"$6}' | sort -u > tar2url.tsv
  echo "[OK] refreshed tar2url.tsv from master"
fi

# join needed tars to urls
sort -u tar.todo.list > tar.todo.sorted
sort -u -k1,1 tar2url.tsv > tar2url.sorted.tsv
join -t $'\t' -1 1 -2 1 tar.todo.sorted tar2url.sorted.tsv > tar.todo.urls.tsv

HAVE_URLS=$(wc -l < tar.todo.urls.tsv)
NEED_TARS=$(wc -l < tar.todo.sorted)
echo "[INFO] URL matches: $HAVE_URLS / $NEED_TARS"
if [[ "$HAVE_URLS" -ne "$NEED_TARS" ]]; then
  echo "[WARN] tar names with no URL mapping:"
  comm -23 tar.todo.sorted <(cut -f1 tar.todo.urls.tsv | sort -u) | head -n 10 | sed 's/^/  - /' || true
fi

echo "[INFO] downloading tarballs to $TARDIR (resume/retry)"
while IFS=$'\t' read -r tar url; do
  if [[ ! -s "$TARDIR/$tar" ]]; then
    echo "[DL] $tar"
    wget -c --retry-connrefused --waitretry=5 --tries=10 \
         --read-timeout=60 --timeout=60 \
         --no-verbose --show-progress \
         -O "$TARDIR/$tar" "$url" || echo -e "$tar\t$url" >> download_failures.log
  fi
done < tar.todo.urls.tsv

# present tarballs (non-empty) — keep list as **names** (not paths)
: > tar.todo.present.list
while read -r t; do
  [[ -s "$TARDIR/$t" ]] && echo "$t" >> tar.todo.present.list
done < tar.todo.list
echo "[INFO] present tarballs: $(wc -l < tar.todo.present.list)"

# ========== 7) extract into REFDIR (flatten to basenames) ==========
echo "[INFO] extracting into $REFDIR (flatten paths)"
if [[ "$USE_PARALLEL" -eq 1 ]]; then
  parallel -j "$JOBS" '
    awk -F"\t" -v t="{1}" '"'"'$1==t{print $2}'"'"' files_by_tar.todo.tsv > list.{#}.txt
    echo "[EXTRACT] {1} members: $(wc -l < list.{#}.txt)"
    tar -xJf "'"$TARDIR"'/{1}" \
        -T list.{#}.txt \
        -C "'"$REFDIR"'" \
        --transform="s|.*/||" \
        --skip-old-files
  ' :::: tar.todo.present.list
else
  while read -r TARF; do
    awk -F'\t' -v t="$TARF" '$1==t{print $2}' files_by_tar.todo.tsv > list.txt
    echo "[EXTRACT] $TARF members: $(wc -l < list.txt)"
    tar -xJf "$TARDIR/$TARF" \
        -T list.txt \
        -C "$REFDIR" \
        --transform="s|.*/||" \
        --skip-old-files
  done < tar.todo.present.list
fi








bsub -q transfer \
  -J move_refs2 \
  -n 1 \
  -R "span[hosts=1] select[mem>2000] rusage[mem=2000]" \
  -M 2000 \
  -o /data/pam/team230/sm71/scratch/rp2/logs/references_2.%J.out \
  -e /data/pam/team230/sm71/scratch/rp2/logs/references_2.%J.err \
  bash -lc 'cd /data/pam/team230/sm71/scratch/rp2/atb; ./references_2.sh'

