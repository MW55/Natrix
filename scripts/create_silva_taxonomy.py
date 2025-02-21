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
import logging

# FASTA parser
def read_fasta(file):
    sequences = {}
    headers = []
    with open(file, 'r') as f:
        header = None
        seq = []
        for line in f:
            if line.startswith('>'):
                if header:
                    sequences[header] = ''.join(seq)
                    headers.append(header.split(" ", 1))
                header = line.strip()[1:]
                seq = []
            else:
                seq.append(line.strip())
        if header:
            sequences[header] = ''.join(seq)
            headers.append(header.split(" ", 1))
    return sequences, headers

p_table = pd.read_csv(snakemake.params.primertable, index_col='Probe')
primertable = p_table.to_dict('index')
data_folder = str(snakemake.params.filename)
file_path_list = sorted(glob.glob(data_folder + "/*.fasta*"))

uniq_seqs = set()
headers = []
for file in file_path_list:
    sequences, file_headers = read_fasta(file)
    uniq_seqs.update(sequences.values())
    headers.extend(file_headers)

df = pd.DataFrame(headers, columns=["id", "taxonomy"])
df = df.set_index(keys="id", drop=True)

df.to_hdf(snakemake.output[0], key='df', mode='w')
