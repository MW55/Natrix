import pandas as pd
import tables

df = pd.read_hdf(str(snakemake.input), 'df')

if snakemake.params.filter_method == "split_sample":
    # If more samples contain a value for a sequence than the unique list of
    # samples without unit identifier, than the sequence is present in at least
    # one _A, _B combination
    uniq_samples =  len(set(sample[:-2] for sample in list(df.columns)))
    df_filtered = df[(df >= 1).sum(axis=1) > uniq_samples]
    df_filtered_out = df[(df >= 1).sum(axis=1) <= uniq_samples]
elif snakemake.params.filter_method == "not_split":
    df_filtered = df[(df > snakemake.params.cutoff).all(1)]
    df_filtered_out = df[(df <= snakemake.params.cutoff).all(1)]
else:
    raise ValueError("Valid filter methods are 'split_sample' and 'not_split'")

df_filtered.to_csv(snakemake.output[0])
df_filtered_out.to_csv(snakemake.output[1])
