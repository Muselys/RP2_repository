#BSUB -J csv2parquet
#BSUB -q normal
#BSUB -n 8                   # CPU threads (Polars will use them)
#BSUB -R "span[hosts=1] select[mem>8000] rusage[mem=8000]"  # per-core mem request
#BSUB -M 8000
#BSUB -o /data/pam/team230/sm71/scratch/rp2/logs/csv2parquet.%J.out
#BSUB -e /data/pam/team230/sm71/scratch/rp2/logs/csv2parquet.%J.err

set -euo pipefail

# 1) activate your env
source ~/.bashrc
conda activate py  # or whatever your env is called
conda install polars -c conda-forge

import polars as pl
# 2) let Polars use all LSF cores
export POLARS_MAX_THREADS=${LSB_DJOB_NUMPROC:-8}

# 3) paths (edit if needed)
CDIR="/data/pam/team230/sm71/scratch/rp2/panaroo_output"
IN="${CDIR}/gene_data.csv"
OUT="${CDIR}/gene_data.parquet"

# 4) convert CSV -> Parquet (ZSTD compression)
python - <<PY
import polars as pl, os
inp = "${IN}"
out = "${OUT}"
print(f"Converting {inp} -> {out}")
lf = pl.scan_csv(inp)
# write once; columnar + compressed
lf.sink_parquet(out, compression="zstd", statistics=True)
print("Done.")
PY