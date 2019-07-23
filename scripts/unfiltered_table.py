import dinopy
import pandas as pd
import collections as col
import tables
import numpy as np
# Script to create a nested dictionary out of all fasta files,
# containing each sequence as a key with a dict as value which,
# in turn, cointains all samples as key and the abundance of
# the particular sequence as value.

f_names = []
seq_dict = col.OrderedDict()
for i in range(len(snakemake.input)):
  seqs = dinopy.FastaReader(snakemake.input[i])
  f_name = snakemake.input[i].split("/")[-1].split(".")[0]
  f_names.append(f_name)
  for entry in seqs.entries():
      if entry.sequence.decode() in seq_dict:
          seq_dict[entry.sequence.decode()][f_name] =(
                  int(str(entry.name).split("size=")[1].split(";")[0]))
      else:
          seq_dict[entry.sequence.decode()] =(
                  {f_name:int(str(entry.name).split("size=")[1].split(";")[0])})

# replace all non-occuring sequences with zero, otherwise will be replaced with nan in dataframe (float64 when replaced)
for key in seq_dict.keys():
        for i in f_names:
            if not i in seq_dict[key]:
                seq_dict[key][i] = 0

# Export the dict in hdf5 and csv format for further processing.
df = pd.DataFrame.from_dict(seq_dict, orient="index", dtype=np.int32) # use int32 to get smaller memmory footprint
df.index.name = "sequences"
df.to_hdf(snakemake.output[1], key='df', mode='w')
df.to_csv(snakemake.output[0])
