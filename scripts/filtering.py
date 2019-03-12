import json
import pandas as pd
import collections as col

with open(str(snakemake.input), 'r') as f_:
    seq_dict = json.load(f_)

filtered_out = col.OrderedDict()
filtered = col.OrderedDict()

if snakemake.params.filter_method == "split_sample":
    # If the list of samples without unit identifier contains only
    # unique strings, the sequence does not apper in both samples for
    # any sample pair.
    for sequence in seq_dict.keys():
        sample_lcomp = [sample[:-2] for sample in seq_dict[sequence].keys()]
        if len(sample_lcomp) == len(set(sample_lcomp)):
            filtered_out[sequence] = seq_dict[sequence]
        else:
            filtered[sequence] = seq_dict[sequence]
elif snakemake.params.filter_method == 'not_split':
    for sequence in seq_dict.keys():
        if all([int(value) <= snakemake.params.cutoff for value
                in seq_dict[sequence].values()]):
            filtered_out[sequence] = seq_dict[sequence]
        else:
            filtered[sequence] = seq_dict[sequence]
else:
    raise ValueError('Valid filter methods are "split_sample" and "not_split"')

df_filtered = pd.DataFrame.from_dict(filtered, orient='index').fillna(0)
df_filtered.index.name = 'sequences'
df_filtered.to_csv(snakemake.output[0])
df_filtered_out = pd.DataFrame.from_dict(filtered_out, orient='index').fillna(0)
df_filtered_out.index.name = 'sequences'
df_filtered_out.to_csv(snakemake.output[1])

with open(str(snakemake.output[2]), 'w') as f_:
    json.dump(filtered, f_)
with open(str(snakemake.output[3]), 'w') as g_:
    json.dump(filtered_out, g_)
