import dinopy
import pandas as pd

# Function to create a nested dictionary out of all fasta files,
# containing each sequence as a key with a dict as value which,
# in turn, cointains all samples as key and the abundance of
# the particular sequence as value. It's quicker than it looks.

def seq_dict_creator(data_files):
    seq_dict = {}
    for i in range(len(data_files)):
        seqs = dinopy.FastaReader(data_files[i])
        f_name = data_files[i].split('/')[-1].split('.')[0]
        for entry in seqs.entries():
            if entry.sequence in seq_dict:
                (seq_dict[entry.sequence.decode()][f_name]
                    = str(entry.name).split('size=')[1].split(';')[0])
            else:
                (seq_dict[entry.sequence.decode()]
                    = {f_name:str(entry.name).split('size=')[1].split(';')[0]})

# Writes the dict to a csv file, index is the sequence with each
# column being the abundance of each sequence for a particular
# sample. Same structure than the old pipeline.

df = pd.DataFrame.from_dict(seq_dict_creator(snakemake.input),
        orient='index').fillna(0)
df.index.name = 'sequences'
df.to_csv(snakmake.output)
