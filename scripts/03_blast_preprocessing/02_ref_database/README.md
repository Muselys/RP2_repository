#After a successful run your atb directory will look like this:
atb/
├── file_list.all.latest.tsv.gz      # Master metadata from OSF
├── target_species.tsv               # Original species list (editable)
├── target_species.norm.tsv          # Normalized species names (lowercased)
├── subset_18.tsv                    # Filtered metadata for target species
├── tar2url.tsv                      # Tarball-to-URL mapping
├── files_by_tar.tsv                 # File members per tarball
├── wanted.ids.txt                  # All target sample IDs
├── have.ids.txt                    # Already-downloaded sample IDs
├── missing.ids.txt                 # Missing IDs (will be downloaded)
├── files_by_tar.todo.tsv           # Tar members to extract
├── tar.todo.list                   # Tarballs to download
├── tarballs/                       # Cached tarballs (resumable)
├── ref/
│   └── fasta/                      # Extracted FASTA files (flattened)
└── manifests/
    ├── versions.txt               # Tool versions used in the run
    ├── checksums.txt              # SHA256 of downloaded files
    └── run.log                    # Run summary (counts, timestamps, etc.)

# default paths
./01_ref_db_target_species.sh

# force rebuild everything
./01_ref_db_target_species.sh --force

# custom work dir + more threads
./01_ref_db_target_species.sh \
  --atb /data/pam/team230/sm71/scratch/rp2/atb \
  --refdir /data/pam/team230/sm71/scratch/rp2/atb/ref/fasta \
  --tardir /data/pam/team230/sm71/scratch/rp2/atb/tarballs \
  --jobs 12
  
# downlaod fastas then normalize + BLAST DB
./build_ref_db.sh --include-gz --shard-size 1000 --max-parallel 12 \
  --normalize-prefix ref19 --makeblastdb --dbtype nucl \
  --module 'blast/2.14.1--pl5321h6f7f691_0'