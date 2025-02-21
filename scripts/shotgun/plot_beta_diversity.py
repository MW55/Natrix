import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load beta diversity results
beta_df = pd.read_csv(snakemake.input.beta_diversity, sep='\t', index_col=0)

# Extract first two principal coordinates for visualization
pc1 = beta_df.iloc[:, 0]
pc2 = beta_df.iloc[:, 1]

# Create PCoA scatter plot
plt.figure(figsize=(8, 6))
sns.scatterplot(x=pc1, y=pc2)
plt.xlabel("PCoA 1")
plt.ylabel("PCoA 2")
plt.title("PCoA Plot of Beta Diversity")
plt.grid(True)
plt.tight_layout()
plt.savefig(snakemake.output.pcoa_plot)
plt.close()

print("Beta diversity visualization completed successfully.")