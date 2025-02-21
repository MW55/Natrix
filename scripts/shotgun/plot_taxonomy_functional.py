import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load inputs
taxonomy_df = pd.read_csv(snakemake.input.taxonomy_matrix, sep='\t', index_col=0)
functional_df = pd.read_csv(snakemake.input.functional_matrix, sep='\t', index_col=0)

# Plot Taxonomy Bar Plot
plt.figure(figsize=(12, 6))
taxonomy_df.sum(axis=1).sort_values(ascending=False).plot(kind='bar', color='skyblue')
plt.xlabel("Samples")
plt.ylabel("Abundance")
plt.title("Taxonomic Composition")
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig(snakemake.output.taxonomy_plot)
plt.close()

# Plot Functional Heatmap
plt.figure(figsize=(12, 8))
sns.heatmap(functional_df, cmap="viridis", xticklabels=True, yticklabels=False)
plt.xlabel("Samples")
plt.ylabel("Functions")
plt.title("Functional Abundance Heatmap")
plt.tight_layout()
plt.savefig(snakemake.output.functional_heatmap)
plt.close()

print("Taxonomy and functional visualization completed successfully.")