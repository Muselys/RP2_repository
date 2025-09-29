#Step 7: This script is the bridge between your messy collection of FASTA files and a clean, BLAST-ready reference database with traceable IDs.
#!/bin/bash
#BSUB -q hugemem
#BSUB -J mkdb_ref_ids
#BSUB -n 4
#BSUB -R "rusage[mem=128GB]"
#BSUB -o mkdb_ids.out
#BSUB -e mkdb_ids.err

set -euo pipefail
module load blast/2.14.1--pl5321h6f7f691_0 || true

BASE="/data/pam/team230/sm71/scratch/rp2/atb"
REF_DIR="$BASE/ref"                          # where the .fa live after extraction
SUBSET="$BASE/subset.tsv"                    # your filtered table
OUT_DIR="$BASE/blastdb"
mkdir -p "$OUT_DIR"

WORK="${TMPDIR:-/tmp}/mkdb_ref"
mkdir -p "$WORK"
MEGA="$WORK/ref_db.fasta"
MAP="$WORK/ref_db.idmap.tsv"

# Build an associative array from subset.tsv for species and tar info
# subset.tsv header: sample  species_sylph  species_miniphy  filename_in_tar_xz  tar_xz  tar_xz_url  tar_xz_md5  tar_xz_size_MB
awk -F'\t' 'NR>1{gsub(/ /,"_",$2); sp[$1]=$2; fn[$1]=$4; tar[$1]=$5} END{
  for (k in sp) { print k"\t"sp[k]"\t"fn[k]"\t"tar[k] }
}' "$SUBSET" | sort > "$WORK/meta.tsv"

# Create MEGA FASTA with unique IDs and a mapping file
# seqid = "sample|species_sylph|filename_in_tar_xz" (all underscores, no spaces)
echo -e "seqid\tsample\tspecies_sylph\tfilename_in_tar_xz\ttar_xz" > "$MAP"

# We’ll read each .fa, detect its sample from the filename (basename without .fa), and rewrite headers.
# Then join to meta.tsv to get species and tar info.
join -t $'\t' -1 1 -2 1 \
  <(awk -F'\t' '{print $1"\t"$0}' "$WORK/meta.tsv" | sort -k1,1) \
  <(
     find "$REF_DIR" -type f -name '*.fa' -print0 \
     | while IFS= read -r -d '' f; do
         sid="$(basename "$f" .fa)"
         echo -e "${sid}\t${f}"
       done | sort -k1,1
   ) \
| awk -F'\t' -v MEGA="$MEGA" -v MAP="$MAP" '
  {
    # Fields after join:
    # 1:sid  2:sid  3:species  4:filename_in_tar_xz  5:tar_xz  6:full_path_to_fa
    sid=$1; species=$3; fin_tar=$4; tarxz=$5; file=$6;
    # Build a safe, unique seqid
    # (species already underscored; sid has no spaces; fin_tar has slashes—keep them for uniqueness if you like)
    seqid = sid "|" species "|" fin_tar;
    gsub(/[ \t]/,"_",seqid);  # just in case
    # Emit mapping
    print seqid "\t" sid "\t" species "\t" fin_tar "\t" tarxz >> MAP;
    # Rewrite FASTA on the fly
    while ((getline line < file) > 0) {
      if (line ~ /^>/) {
        # replace defline with our unique seqid; keep any extra info after first space as description if desired
        print ">" seqid;
      } else {
        print line;
      }
    }
    close(file);
  }
' > "$MEGA"

#STEP 8: make reference database:
