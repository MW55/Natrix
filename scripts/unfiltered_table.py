import dinopy
import pandas as pd
import tables
import numpy as np

f_names = []
uniq_seqs = []
for i in range(len(snakemake.input)):
    seqs = dinopy.FastaReader(snakemake.input[i])
    uniq_seqs = uniq_seqs + list(set([entry.sequence.decode() for entry in seqs.entries()]))
uniq_seqs = set(uniq_seqs)

# create empty matrix and fill, all other solutions cost too much memory
sample_names = [i.split("/")[-1].split(".")[0] for i in snakemake.input]
df = pd.DataFrame(0, index=uniq_seqs, columns=sample_names, dtype=np.uint16)

del uniq_seqs

# fill matrix
for i in range(len(snakemake.input)):
    sample_name = sample_names[i]
    seqs = dinopy.FastaReader(snakemake.input[i])
    for entry in seqs.entries():
        seq = entry.sequence.decode()
        value = np.uint16(entry.name.decode().split("size=")[1].split(";")[0])
        df.at[seq, sample_name] = value

# save to file
df.index.name = "sequences"
df.to_hdf(snakemake.output[1], key='df', mode='w')
df.to_csv(snakemake.output[0])
