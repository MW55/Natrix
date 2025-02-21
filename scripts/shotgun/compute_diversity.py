import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist, squareform
from skbio.diversity import alpha_diversity
from skbio.stats.ordination import pcoa
from skbio.stats.distance import DistanceMatrix

df = pd.read_csv(snakemake.input.taxonomy_matrix, sep='\t', index_col=0)

df = df.apply(pd.to_numeric, errors='coerce')

if df.empty or df.isnull().all().all():
    raise ValueError("Error: The input taxonomy matrix is empty or contains only non-numeric values.")

# Normalize data to relative abundances (important for Bray-Curtis)
df = df.div(df.sum(axis=1), axis=0)

# Compute Alpha Diversity Metrics Individually
alpha_metrics = ['shannon', 'simpson']
alpha_results = {}

for metric in alpha_metrics:
    alpha_results[metric] = alpha_diversity(metric=metric, counts=df.to_numpy(), ids=df.index.to_list())

# Convert to DataFrame
alpha_df = pd.DataFrame(alpha_results, index=df.index)
alpha_df.to_csv(snakemake.output.alpha_diversity, sep='\t')

# Compute Beta Diversity (PCoA using Bray-Curtis)
beta_matrix = squareform(pdist(df.to_numpy(), metric='braycurtis'))

"""
pcoa_results = pcoa(beta_matrix)
beta_df = pd.DataFrame(pcoa_results.samples, index=df.index)
beta_df.to_csv(snakemake.output.beta_diversity, sep='\t')
"""

distmat = DistanceMatrix(beta_matrix, list(df.index))
pcoa_results = pcoa(distmat)
beta_df = pcoa_results.samples

beta_df.to_csv(snakemake.output.beta_diversity, sep='\t')

print("Diversity analysis completed successfully.")
