#STEP 1
#sample, species_sylph, species_miniphy, filename_in_tar_xz, tar_xz, tar_xz_url, tar_xz_md5, tar_xz_size_MB.
#We’ll use that to filter by genus and figure out which tarballs to fetch.
#Make a directory called atb (AllTheBacteria)
mkdir -p atb && cd atb

#STEP 2
# get the latest file list (tab-delimited; includes species_sylph + tar info)
# note: URL may change; this is the one the docs show.
wget https://osf.io/download/4yv85/ -O file_list.all.latest.tsv.gz

zcat /data/pam/team230/sm71/scratch/rp2/atb/file_list.all.latest.tsv.gz |  tail -n +2 | wc -l

#2440377 actual assemblies in the AllTheBacteria Genome database

#STEP 3.1
tail -n +2 /data/pam/team230/sm71/scratch/rp2/atb/subset.tsv | wc -l
#410304 subset.tsv actual assemblies in the AllTheBacterua Genome database that are from the Genera Staphylococcus, Streptococcus and Enterococcus

#Work out which archives to download (unique tar name + URL):
awk -F'\t' 'NR>1 {print $5"\t"$6}' subset.tsv | sort -u > tar2url.tsv



#STEP 3.2
#check the presence of the target species in the AllTheBacteria Genome database using:
#/data/pam/team230/sm71/scratch/rp2/metadata/target_species.tsv
awk -F'\t' '
  NR==FNR {targets[$1]; next}         # read target species into array
  FNR==1 {next}                       # skip header in subset.tsv
  {
    # clean the genus part of the observed species
    split($2, a, " ");                # split "Genus species..." on space
    g = a[1];                         # first token = Genus or Genus_suffix
    sub(/_.*/, "", g);                # drop any suffix after underscore
    observed[$2];                     # keep full species string
    observed_genus[g];                # also store genus-only form
  }
  END {
    for (t in targets) {
      if (t in observed) {
        print t "\tTRUE";             # exact match in subset.tsv
      } else {
        # Check genus-only match (handles _A, _B cases)
        split(t, b, " ");             # genus of target
        tg = b[1];
        sub(/_.*/, "", tg);
        if (tg in observed_genus) print t "\tTRUE";
        else print t "\tFALSE";
      }
    }
  }
' /data/pam/team230/sm71/scratch/rp2/metadata/target_species.tsv \
  /data/pam/team230/sm71/scratch/rp2/atb/subset.tsv \
> target_species_presence.txt

#All target species are present!

#STEP 5: Download those archives (rename correctly for OSF):

# no parallel
while IFS=$'\t' read -r tar url; do
  echo "Downloading $tar"
  wget -O "$tar" "$url"
done < tar2url.tsv

# optional: speed up if you have GNU parallel
# parallel -j 4 --colsep '\t' 'wget -O {1} {2}' :::: tar2url.tsv

#STEP 6: Extract only your target FASTAs from each tar.xz (no wasted I/O):
# Make a per-archive file list of the members we want to extract
#$5 = tar_xz → the archive file name (e.g. atb.assembly.r0.2.batch.110.tar.xz)
#$4 = filename_in_tar_xz → the path inside the archive (e.g. atb.assembly.r0.2.batch.110/SAMEA...fa)
#after the list of requires fasta files are made in lists.txt, it is then downloaded in the the ref directory
awk -F'\t' 'NR>1 {print $5"\t"$4}' subset.tsv | sort > files_by_tar.tsv



#EXTRACT_REFS.SH
#!/bin/bash
#BSUB -q hugemem              # use hugemem queue
#BSUB -J extract_refs         # job name
#BSUB -n 1                    # number of cores
#BSUB -R "rusage[mem=128GB]"  # adjust memory request if needed
#BSUB -o extract_refs.out     # stdout log
#BSUB -e extract_refs.err     # stderr log

# Load any modules you need for tar (if your system uses modules)
# module load tar

# Move into your working directory
cd /data/pam/team230/sm71/scratch/rp2/atb

# Loop over tarballs and extract only target FASTAs into ./ref/
cut -f1 files_by_tar.tsv | sort -u | while read -r tarfile; do
  echo "Extracting from $tarfile"
  grep -P "^${tarfile}\t" files_by_tar.tsv | cut -f2 > list.txt
  # Extract into ./ref
  tar -xJf "$tarfile" -T list.txt -C ref
done



#when it finishes i'll have a bunch of paths like this: atb.assembly.incr_release.202408.batch.X/SAMxxxxx.fa


#EXTRACT_REFS.SH
#!/bin/bash
#BSUB -q hugemem              # use hugemem queue
#BSUB -J extract_refs         # job name
#BSUB -n 1                    # number of cores
#BSUB -R "rusage[mem=128GB]"  # adjust memory request if needed
#BSUB -o extract_refs.out     # stdout log
#BSUB -e extract_refs.err     # stderr log

# Load any modules you need for tar (if your system uses modules)
# module load tar

# Move into your working directory
cd /data/pam/team230/sm71/scratch/rp2/atb

# Loop over tarballs and extract only target FASTAs into ./ref/
cut -f1 files_by_tar.tsv | sort -u > tar.list
export FILES_BY_TAR=files_by_tar.tsv
parallel -j 6 '
  grep -P "^{1}\t" "$FILES_BY_TAR" | cut -f2 > list.{#}.txt &&
  tar --use-compress-program="xz -T0 -dc" -xf {1} -T list.{#}.txt -C ref --skip-old-files
' :::: tar.list

#sanity checks once finished:
# Expect these to match (or be very close)
echo "Expected: $(tail -n +2 subset.tsv | wc -l)"
echo "Found:    $(find ref -type f -name '*.fa' | wc -l)"

# Spot check a file looks like FASTA
head -n 3 "$(find ref -type f -name '*.fa' | head -n1)"
