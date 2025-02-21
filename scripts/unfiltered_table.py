import pandas as pd
import numpy as np
import os
import re
import yaml
import glob
import sys
import shutil
import subprocess
import pathlib

# FASTA parser
def read_fasta(file):
    sequences = {}
    with open(file, 'r') as f:
        header = None
        seq = []
        for line in f:
            if line.startswith('>'):
                if header:
                    sequences[header] = ''.join(seq)
                header = line.strip()[1:]
                seq = []
            else:
                seq.append(line.strip())
        if header:
            sequences[header] = ''.join(seq)
    return sequences

p_table = pd.read_csv(snakemake.params.primertable, index_col='Probe')
primertable = p_table.to_dict('index')
data_folder = str(snakemake.params.filename)
file_path_list = sorted(glob.glob(data_folder + "/*.fasta*"))

def process_fasta_files(file_list):
    f_names = []
    uniq_seqs = set()
    for file in file_list:
        sequences = read_fasta(file)
        uniq_seqs.update(sequences.values())
    return list(uniq_seqs)

uniq_seqs = process_fasta_files(file_path_list)

sample_names = [i.split("/")[-1].split(".")[0] for i in snakemake.input]
df = pd.DataFrame(0, index=uniq_seqs, columns=sample_names, dtype=np.uint16)

for i, file in enumerate(snakemake.input):
    sample_name = sample_names[i]
    sequences = read_fasta(file)
    for header, seq in sequences.items():
        value = np.uint16(header.split("size=")[1].split(";")[0]) if "size=" in header else 1
        df.at[seq, sample_name] = value

df.index.name = "sequences"
df.to_hdf(snakemake.output[1], key='df', mode='w')
df.to_csv(snakemake.output[0])