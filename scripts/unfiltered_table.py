import dinopy
import json
import pandas as pd
import collections as col
# Function to create a nested dictionary out of all fasta files,
# containing each sequence as a key with a dict as value which,
# in turn, cointains all samples as key and the abundance of
# the particular sequence as value. It's quicker than it looks.

def seq_dict_creator(data_files):
    seq_dict = col.OrderedDict()
    for i in range(len(data_files)):
        seqs = dinopy.FastaReader(data_files[i])
        f_name = data_files[i].split('/')[-1].split('.')[0]
        for entry in seqs.entries():
            if entry.sequence.decode() in seq_dict:
                seq_dict[entry.sequence.decode()][f_name] =(
                        str(entry.name).split('size=')[1].split(';')[0])
            else:
                seq_dict[entry.sequence.decode()] =(
                        {f_name:str(entry.name).split('size=')[1].split(';')[0]})
    return seq_dict

seq_dict = seq_dict_creator(snakemake.input)

# Export the dict im yaml format for further processing.
# This increases the modularity of the pipeline but cost time to
# read/write the file. I could put the other functions working
# with the dict in this script to work on the dict in memory,
# should the reading/writing take up too much space

#with open(str(snakemake.output[1]), 'w') as f_:
#    yaml.dump(seq_dict, stream = f_)

with open(str(snakemake.output[1]), 'w') as f_:
    json.dump(seq_dict, f_)


# Write the dict to a csv file, index is the sequence with each
# column being the abundance of each sequence for a particular
# sample. Same structure than the old pipeline.

df = pd.DataFrame.from_dict(seq_dict, orient='index').fillna(0)
df.index.name = 'sequences'
df.to_csv(snakemake.output[0])


# Export the dict im yaml format for further processing.
# This increases the modularity of the pipeline but cost time to
# read/write the file. I could put the other functions working
# with the dict in this script to work on the dict in memory,
# should the reading/writing take up too much space
