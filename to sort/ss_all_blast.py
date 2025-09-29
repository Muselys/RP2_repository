#!/usr/bin/env python3
import argparse, os, sys, gzip, glob, gc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from matplotlib.colors import PowerNorm


def log(m): print(m, file=sys.stderr)


def list_blast_files(results_dir):
    root = Path(results_dir).resolve()
    files = []
    for pat in ("*.out", "*.out.gz"):
        for f in root.rglob(pat):
            if not f.is_file():
                continue
            # Skip top-level .out files (usually scheduler logs)
            if f.parent == root:
                continue
            # Skip anything under reports/
            if "reports" in f.parts:
                continue
            # Skip known job-log patterns
            if f.name.startswith(("build_reports.", "lsf.", "bsub.")):
                continue
            files.append(str(f))
    return sorted(files)


def list_samples_from_ref_dir(ref_dir):
    return sorted({os.path.splitext(f)[0] for f in os.listdir(ref_dir)
                   if f.endswith((".fa", ".fasta"))})

def load_subset_sample2species(subset_path):
    sub = pd.read_csv(subset_path, sep="\t", dtype=str)
    if sub.shape[1] < 2:
        raise ValueError("subset.tsv must have at least two columns: sample, species.")
    sub = sub.iloc[:, :2]
    sub.columns = ["sample", "species"]
    return sub.drop_duplicates("sample")

def canon_spaces(s):
    s = (s or "").replace("_", " ").strip()
    return " ".join(s.split())

def primary_binomial(s):
    s = canon_spaces(s)
    first = s.split(";", 1)[0].strip()
    toks = first.split()
    if len(toks) < 2: return first
    genus = toks[0].split("_")[0]  # strip _A/_B
    species = toks[1]
    return f"{genus} {species}"

def _extract_sample_fast(s: pd.Series) -> pd.Series:
    """
    Cheap, vectorised: take prefix before first '.'.
    Example: 'SAMPLE01.contig00042_len...' -> 'SAMPLE01'
    """
    # str.partition avoids regex; returns (head, sep, tail)
    return s.str.partition('.', expand=True)[0]

def read_blast_qs(paths, samples_in_ref, chunksize=1_000_000, compact_every=25):
    """
    Stream BLAST .out/.out.gz files and return unique (qseqid, sample) pairs
    limited to samples_in_ref. Avoid regex; avoid large in-memory spikes.
    """
    samples_in_ref = set(samples_in_ref)  # O(1) membership
    partial = []
    total_rows = 0
    total_kept = 0

    for i, p in enumerate(paths, 1):
        try:
            # pandas handles gz automatically
            for chunk in pd.read_csv(
                p, sep="\t", header=None, usecols=[0,1],
                names=["qseqid", "sseqid"], dtype=str,
                engine="c", on_bad_lines="skip",
                chunksize=chunksize
            ):
                total_rows += len(chunk)
                if chunk.empty:
                    continue

                # FAST sample extraction (no regex)
                chunk["sample"] = _extract_sample_fast(chunk["sseqid"])

                # Keep only columns we actually need
                chunk = chunk[["qseqid", "sample"]]

                # Filter early to shrink memory
                chunk = chunk[chunk["sample"].isin(samples_in_ref)]
                if chunk.empty:
                    continue

                # Drop dupes within chunk
                chunk = chunk.drop_duplicates(["qseqid", "sample"])

                # Use smaller dtypes to save RAM
                chunk["qseqid"] = chunk["qseqid"].astype("string")
                chunk["sample"] = chunk["sample"].astype("category")

                total_kept += len(chunk)
                partial.append(chunk)

        except Exception as e:
            log(f"[warn] failed {p}: {e}")
            continue

        # Periodic compaction to cap memory
        if i % compact_every == 0 and partial:
            df = pd.concat(partial, ignore_index=True)
            df = df.drop_duplicates(["qseqid", "sample"])
            partial = [df]  # keep one compacted frame
            gc.collect()
            log(f"  read {i}/{len(paths)} files; rows seen: {total_rows:,}; kept so far: {total_kept:,}")

    if not partial:
        return pd.DataFrame(columns=["qseqid","sample"])

    # Final compact
    qs = pd.concat(partial, ignore_index=True).drop_duplicates(["qseqid","sample"])
    # Ensure plain dtypes for downstream merge/pivot
    qs["qseqid"] = qs["qseqid"].astype(str)
    qs["sample"] = qs["sample"].astype(str)
    return qs

def plot_heatmap(df, out_png, title):
    if df.empty:
        log(f"[warn] empty matrix for plot {out_png}"); return
    # Guard against absurdly wide/tall plots
    if df.shape[0] > 5000 or df.shape[1] > 500:
        log(f"[warn] matrix too large for plotting ({df.shape[0]}x{df.shape[1]}). Skipping {out_png}.")
        return
    fig, ax = plt.subplots(figsize=(max(8, df.shape[1]*0.4), max(6, df.shape[0]*0.2)))
    norm = PowerNorm(gamma=0.6, vmin=0.0, vmax=1.0)
    im = ax.imshow(df.values, aspect='auto', interpolation='nearest', cmap='Blues', norm=norm)
    ax.set_xticks(np.arange(df.shape[1])); ax.set_xticklabels(df.columns, rotation=90)
    ax.set_yticks(np.arange(df.shape[0])); ax.set_yticklabels(df.index)
    ax.set_xlabel("Reference species"); ax.set_ylabel("qseqid (genes)")
    ax.set_title(title)
    cbar = plt.colorbar(im, ax=ax); cbar.set_label("Detection rate")
    fig.tight_layout(); fig.savefig(out_png, dpi=200); plt.close(fig)

def main():
    ap = argparse.ArgumentParser(description="Build gene×species detection matrix + heatmaps (no focus species).")
    ap.add_argument("--results-dir", required=True, help="Root dir holding species subdirs with BLAST .out/.out.gz files.")
    ap.add_argument("--subset", required=True, help="ATB subset.tsv (header). col1=sample, col2=species.")
    ap.add_argument("--ref-fasta-dir", required=True, help="Directory with reference FASTAs (one file per sample).")
    ap.add_argument("--reports-dir", default="reports", help="Output directory.")
    ap.add_argument("--chunksize", type=int, default=1_000_000, help="Rows per chunk when streaming BLAST outputs.")
    ap.add_argument("--compact-every", type=int, default=25, help="Compact partial frames every N files.")
    args = ap.parse_args()

    reports = Path(args.reports_dir); reports.mkdir(parents=True, exist_ok=True)

    log(">>> scanning BLAST files...")
    blast_files = list_blast_files(args.results_dir)
    if not blast_files: sys.exit("No BLAST .out files found.")
    (reports / "blast_files_used.txt").write_text("\n".join(blast_files) + "\n")
    log(f"  found {len(blast_files)} files")

    log(">>> listing samples in reference...")
    samples_in_ref = set(list_samples_from_ref_dir(args.ref_fasta_dir))
    log(f"  samples_in_ref: {len(samples_in_ref)}")

    log(">>> loading subset (sample->species)...")
    sub = load_subset_sample2species(args.subset)
    sub = sub[sub["sample"].isin(samples_in_ref)].copy()
    if sub.empty: sys.exit("No overlap between subset.tsv and ref FASTA filenames.")
    sub["species_norm"] = sub["species"].map(primary_binomial)

    # denominators (genome counts per species)
    denom = sub.groupby("species_norm", as_index=False)["sample"].nunique()
    denom.columns = ["species_norm","n_genomes"]

    log(">>> reading BLAST qseqid/sseqid across all species...")
    qs = read_blast_qs(blast_files, samples_in_ref, chunksize=args.chunksize, compact_every=args.compact_every)
    if qs.empty: sys.exit("No usable BLAST hits after filtering to ref samples.")

    # presence per gene×sample
    qs["present"] = 1
    qs = qs.drop_duplicates(["qseqid","sample"])

    # map to species
    qs = qs.merge(sub[["sample","species_norm"]], on="sample", how="left").dropna(subset=["species_norm"])

    # aggregate to gene×species counts
    agg = qs.groupby(["qseqid","species_norm"], as_index=False)["present"].sum()
    mat = agg.pivot(index="qseqid", columns="species_norm", values="present").fillna(0)

    # ensure all species columns present
    all_species = denom["species_norm"].tolist()
    missing = [sp for sp in all_species if sp not in mat.columns]
    if missing:
        mat = pd.concat([mat, pd.DataFrame(0, index=mat.index, columns=missing)], axis=1)

    # consistent column order
    mat = mat[sorted(mat.columns)]

    # convert to detection rates (float32 to save RAM/disk)
    denom_s = denom.set_index("species_norm")["n_genomes"]
    mat = mat.divide(denom_s, axis=1).fillna(0.0).astype(np.float32)

    # write full matrix
    matrix_csv = reports / "detection_rate_matrix.csv"
    mat.to_csv(matrix_csv)
    log(f"[ok] wrote {matrix_csv}")

    # full heatmap (order rows by overall mean detection rate descending)
    row_order = mat.mean(axis=1).sort_values(ascending=False).index
    mat_full_plot = mat.loc[row_order]
    full_png = reports / "heatmap_detection_rates_full.png"
    plot_heatmap(mat_full_plot, full_png, "Detection rate per gene × species (full)")
    log(f"[ok] wrote {full_png}")

    # per-species reports for ALL species columns
    for sp in mat.columns:
        others = mat.drop(columns=[sp], errors="ignore")
        nonmax = others.max(axis=1) if others.shape[1] else 0.0
        rel = (mat[sp] - nonmax).rename("reliability").sort_values(ascending=False)

        mat_rows = mat.loc[rel.index]
        other_cols = [c for c in mat.columns if c != sp]
        mat_plot = mat_rows[[sp] + other_cols]

        base = reports / sp.replace(" ", "_")
        mat_plot.to_csv(f"{base}_matrix.csv")
        rel.to_csv(f"{base}_gene_reliability.tsv", sep="\t", header=True)
        plot_heatmap(mat_plot, f"{base}_heatmap.png", f"Detection rates — focus: {sp}")
        log(f"[ok] {sp}: wrote {base}_matrix.csv, {base}_gene_reliability.tsv, {base}_heatmap.png")

    log(">>> done.")

if __name__ == "__main__":
    main()
