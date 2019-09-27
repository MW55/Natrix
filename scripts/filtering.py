import pandas as pd
import tables
import numpy as np

df = pd.read_hdf(str(snakemake.input), 'df')

if snakemake.params.filter_method == "split_sample":
    samples =  set(map(lambda x: x.split("_")[0], list(df.columns)))
    for sample in samples:
        pair = [ sample+'_A', sample+'_B']
        # replace value by zero where XOR is True for A and B pair of samples,
        # since where() replace False positions, replace the invert of xor
        df[pair] = df[pair].where(np.invert(np.logical_xor(df[pair[0]]>=1, df[pair[1]]>=1)),0)
    # finally remove all empty rows
    selected_rows = df.sum(axis=1)>0
    df_filtered = df[selected_rows]
    df_filtered_out = df[np.invert(selected_rows)]
elif snakemake.params.filter_method == "not_split":
    df_filtered = df[(df > snakemake.params.cutoff).all(1)]
    df_filtered_out = df[(df <= snakemake.params.cutoff).all(1)]
else:
    raise ValueError("Valid filter methods are 'split_sample' and 'not_split'")

df_filtered.to_csv(snakemake.output[0])
df_filtered_out.to_csv(snakemake.output[1])
