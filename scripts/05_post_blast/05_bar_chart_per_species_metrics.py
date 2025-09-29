import pandas as pd
import matplotlib.pyplot as plt

# load metrics
pre = pd.read_csv("per_species_metrics_pre.tsv", sep="\t")
post = pd.read_csv("per_species_metrics_post.tsv", sep="\t")

def plot_metrics(df, title, outpng):
    df = df.sort_values("species")
    x = range(len(df))
    width = 0.35
    fig, ax = plt.subplots(figsize=(14,6))
    ax.bar([i - width/2 for i in x], df["Sensitivity"], width, label="Sensitivity")
    ax.bar([i + width/2 for i in x], df["Specificity"], width, label="Specificity")
    ax.set_xticks(x)
    ax.set_xticklabels(df["species"], rotation=90)
    ax.set_ylim(0,1.05)
    ax.set_ylabel("Value")
    ax.set_title(title)
    ax.legend()
    plt.tight_layout()
    plt.savefig(outpng, dpi=300)

plot_metrics(pre, "Per-species metrics (pre-filter)", "metrics_pre.png")
plot_metrics(post, "Per-species metrics (post-filter)", "metrics_post.png")
