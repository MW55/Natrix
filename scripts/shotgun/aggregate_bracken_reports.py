import os
import pandas as pd

# Read inputs and outputs from Snakemake
input_files = snakemake.input
output_file = snakemake.output.taxonomy_matrix

print("Input files:", list(input_files))


all_dfs = []

for filepath in input_files:
    # Derive the sample name from the filename
    # e.g. "A071SE_A_bracken_report.txt" -> "A071SE_A"
    sample_name = os.path.basename(filepath).replace("_bracken_report.txt", "")

    # Read the Bracken report
    df = pd.read_csv(filepath, sep="\t")

    # Typical Bracken columns include:
    #   name, taxonomy_id, taxonomy_lvl, kraken_assigned_reads,
    #   added_reads, new_est_reads, fraction_total_reads
    # We'll just keep the core columns needed (taxonomy_id and new_est_reads).
    df = df[["taxonomy_id", "new_est_reads"]]

    # Set 'taxonomy_id' as the index, so we can transpose and get columns = taxonomy IDs
    df.set_index("taxonomy_id", inplace=True)

    # Rename 'new_est_reads' column to the sample name
    # After transpose, row = sample, col = taxonomy_id
    df.rename(columns={"new_est_reads": sample_name}, inplace=True)

    # Transpose so we have exactly one row (the sample) and many columns (the taxonomy IDs)
    df = df.T

    # Collect for merging
    all_dfs.append(df)

# Concatenate all single-row dataframes into one multi-row dataframe
merged_df = pd.concat(all_dfs, axis=0)

# Replace missing values with 0
merged_df.fillna(0, inplace=True)

# Save to file
merged_df.to_csv(output_file, sep="\t")

print(f"Aggregated {len(input_files)} Bracken reports.")
print(f"Resulting matrix has {merged_df.shape[0]} samples (rows) and {merged_df.shape[1]} taxonomy IDs (columns).")
print(f"Output written to: {output_file}")
