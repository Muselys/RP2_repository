python3 - <<'PY'
from pathlib import Path
import polars as pl
import pandas as pd
import numpy as np
import plotly.express as px

BASE = Path("/data/pam/team230/sm71/scratch/rp2/run_blast/results/reports")

SPECIES = [
    "Enterococcus_faecalis",
    "Enterococcus_faecium",
    "Staphylococcus_aureus",
    "Staphylococcus_epidermidis",
    "Staphylococcus_capitis",
    "Staphylococcus_haemolyticus",
    "Staphylococcus_pseudintermedius",
    "Staphylococcus_argenteus",
    "Staphylococcus_sciuri",
    "Streptococcus_agalactiae",
    "Streptococcus_dysgalactiae",
    "Streptococcus_equi",
    "Streptococcus_mitis",
    "Streptococcus_mutans",
    "Streptococcus_pneumoniae",
    "Streptococcus_pyogenes",
    "Streptococcus_suis",
    "Streptococcus_uberis",
]

# collect per (gene, species): presence, max %hits, n_qseqid
rows = []
for sp in SPECIES:
    f = BASE / sp / f"{sp}.per_qseqid.table.tsv"
    if not f.exists():
        print(f"[SKIP] missing: {f}")
        continue

    df = pl.read_csv(f, separator="\t")

    # keep only rows with % of hits > 0
    if "% of hits" not in df.columns:
        raise SystemExit(f"Missing '% of hits' column in {f}")
    df = df.filter(pl.col("% of hits") > 0)

    if df.is_empty():
        continue

    # explode comma-separated "Genes" into one-per-row
    exploded = (
        df.select([
            pl.col("qseqid"),
            pl.col("Genes"),
            pl.col("% of hits").alias("pct_hits"),
        ])
        .with_columns(
            pl.when(pl.col("Genes").is_null() | (pl.col("Genes").str.strip_chars() == ""))
              .then(pl.lit([]))
              .otherwise(pl.col("Genes").str.split(","))
              .alias("genes_list")
        )
        .explode("genes_list")
        .with_columns(pl.col("genes_list").str.strip_chars().alias("gene"))
        .filter(pl.col("gene") != "")
    )

    if exploded.is_empty():
        continue

    # aggregate per gene within species
    agg = (
        exploded.group_by("gene")
        .agg([
            pl.col("qseqid").n_unique().alias("n_qseqid"),
            pl.col("pct_hits").max().alias("pct_hits_max"),
        ])
    )

    # append to list
    for r in agg.iter_rows(named=True):
        rows.append({
            "gene": r["gene"],
            "species": sp,
            "present": 1,
            "n_qseqid": int(r["n_qseqid"]),
            "pct_hits_max": float(r["pct_hits_max"]),
        })

if not rows:
    raise SystemExit("No data after filtering (% of hits > 0).")

df_all = pd.DataFrame(rows)

# pivot to matrices
presence = df_all.pivot_table(index="gene", columns="species", values="present",
                              aggfunc="max", fill_value=0).astype(int)

pct_max = df_all.pivot_table(index="gene", columns="species", values="pct_hits_max",
                             aggfunc="max")  # keep NaN where absent

n_qseq = df_all.pivot_table(index="gene", columns="species", values="n_qseqid",
                            aggfunc="sum")  # NaN where absent

# build hover text using pct_max and n_qseq
# fill NaNs for display (but keep presence for color)
pct_disp = pct_max.copy()
nq_disp = n_qseq.copy()
pct_disp = pct_disp.where(~pct_disp.isna(), other=0.0)
nq_disp = nq_disp.where(~nq_disp.isna(), other=0)

hover = pct_disp.copy()
for c in hover.columns:
    hover[c] = (
        "Gene: " + hover.index.astype(str) +
        "<br>Species: " + c +
        "<br>Present: " + presence[c].map({0:"No", 1:"Yes"}).astype(str) +
        "<br>Max % of hits: " + pct_disp[c].round(2).astype(str) +
        "<br># qseqid (>% hits 0): " + nq_disp[c].astype(int).astype(str)
    )

# plotly heatmap: color by presence (0/1), hover shows % and counts
fig = px.imshow(
    presence.values,
    x=presence.columns.tolist(),
    y=presence.index.tolist(),
    color_continuous_scale="Viridis",
    aspect="auto",
    origin="upper",
)

# attach custom hover text
fig.update_traces(
    hovertemplate="%{customdata}",
    customdata=hover.values,
    hoverongaps=False,
    zmin=0, zmax=1,
    showscale=True,
    colorbar_title="Presence",
)

fig.update_layout(
    title="Gene presence across species (hover shows max % of hits)",
    xaxis_title="Species",
    yaxis_title="Gene",
)

out_html = BASE / "gene_species_heatmap_presence_hover.html"
fig.write_html(str(out_html), include_plotlyjs="cdn")
print(f"[OK] Interactive heatmap -> {out_html}")
PY
