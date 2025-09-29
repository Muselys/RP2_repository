

# remove the blocked defaults
conda config --remove-key default_channels

# create env with only conda-forge + bioconda + nodefaults
conda create -y -n phylign --override-channels \
  -c conda-forge -c bioconda -c nodefaults \
  snakemake python=3.11

  # activate + sanity check
conda activate phylign
conda install -y -c bioconda pigz
module load minimap2/2.28--h577a1d6_4
module load samtools-1.19.2/ 
which snakemake && snakemake --version

ASMS=/data/pam/team230/sm71/scratch/rp2/Phylign/asms
COBS=/data/pam/team230/sm71/scratch/rp2/Phylign/cobs
INPUT="/data/pam/team230/sm71/scratch/rp2/run_blast/queries"
BATCHLIST=/data/pam/team230/sm71/scratch/rp2/Phylign/data/batches_2m.txt

ASM_URL=https://ftp.ebi.ac.uk/pub/databases/AllTheBacteria/Releases/0.2/assembly/
COBS_URL=https://ftp.ebi.ac.uk/pub/databases/AllTheBacteria/Releases/0.2/indexes/phylign/

ASM_SRC=/data/pam/software/AllTheBacteria/Releases/0.2/assembly
COB_SRC=/data/pam/software/AllTheBacteria/Releases/0.2/indexes/phylign

#Dry run - preview, no changes
# Assemblies 
curl -s https://ftp.ebi.ac.uk/pub/databases/AllTheBacteria/Releases/0.2/assembly/ \
| grep -Eio 'href="[^"]+\.asm\.tar\.xz"' \
| sed 's/^href="//; s/"$//' \
| grep -Ei '^(streptococcus|staphylococcus|enterococcus)'

# COBS 
curl -s https://ftp.ebi.ac.uk/pub/databases/AllTheBacteria/Releases/0.2/indexes/phylign/ \
| grep -Eio 'href="[^"]+\.cobs_classic\.xz"' \
| sed 's/^href="//; s/"$//' \
| grep -Ei '^(streptococcus|staphylococcus|enterococcus)'



#Symlink assemblies → asms/
cd /data/pam/team230/sm71/scratch/rp2/Phylign/asms

for f in $(curl -s https://ftp.ebi.ac.uk/pub/databases/AllTheBacteria/Releases/0.2/assembly/ \
  | grep -Eio 'href="[^"]+\.asm\.tar\.xz"' \
  | sed 's/^href="//; s/"$//' \
  | grep -Ei '^(streptococcus|staphylococcus|enterococcus)'); do
  ln -s /data/pam/software/AllTheBacteria/Releases/0.2/assembly/"$f" .
done

#Symlink indexes → cobs/
cd /data/pam/team230/sm71/scratch/rp2/Phylign/cobs

for f in $(curl -s https://ftp.ebi.ac.uk/pub/databases/AllTheBacteria/Releases/0.2/indexes/phylign/ \
  | grep -Eio 'href="[^"]+\.cobs_classic\.xz"' \
  | sed 's/^href="//; s/"$//' \
  | grep -Ei '^(streptococcus|staphylococcus|enterococcus)'); do
  ln -s /data/pam/software/AllTheBacteria/Releases/0.2/indexes/phylign/"$f" .
done

# Symlink all *_ssc_candidates.fa files into input/
cd /data/pam/team230/sm71/scratch/rp2/Phylign/input
ln -s /data/pam/team230/sm71/scratch/rp2/run_blast/queries/*_ssc_candidates.fa .

#Querying a subset of the AllTheBacteria dataset
ls -1 /data/pam/team230/sm71/scratch/rp2/Phylign/asms/*.asm.tar.xz \
| xargs -n1 basename \
| sed 's/\.asm\.tar\.xz$//' \
| grep -E '^(streptococcus|staphylococcus|enterococcus)' \
| sort \
> /data/pam/team230/sm71/scratch/rp2/Phylign/data/batches_2m.txt


#Sanity checks (counts + mismatches)
ls -1 /data/pam/team230/sm71/scratch/rp2/Phylign/asms | wc -l
ls -1 /data/pam/team230/sm71/scratch/rp2/Phylign/cobs | wc -l
ls -1 /data/pam/team230/sm71/scratch/rp2/Phylign/input | wc -l

#Clean u intermediates from previous runs:
make clean

#Run the pipeline
make              #execute Snakemake with the corresponding parameters
make match        
make map



