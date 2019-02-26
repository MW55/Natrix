import yaml
import pandas as pd
import numpy as np
from glob import glob
import sys
# Create the datatable containing the samples, units and paths of all
# fastq files formatted correctly. This is vital for the snakemake
# pipeline, without it, the wildcards can't be created.

# might be better to glob abspaths
with open(sys.argv[1]) as f_:
    config = yaml.load(f_)

fpl = sorted(glob('demultiplexed/' + "*.gz"))
fl = sorted([file_.split('/')[-1] for file_ 
                    in glob('demultiplexed/' + "*.gz")])
if config['merge']['paired_End']:
    df = pd.DataFrame(columns=['sample', 'unit', 'fq1', 'fq2'],
        index =range(int(len(fl)/2)), dtype=str)
    i, j = (0, 0)
    while i < len(fl)/2:
        df.loc[i]['sample'] = fl[j].split('_')[0]
        df.loc[i]['unit'] = fl[j].split('_')[1]
        df.loc[i]['fq1'] = fpl[j][:-3]
        df.loc[i]['fq2'] = fpl[j+1][:-3]
        j += 2
        i += 1
else:
    df = pd.DataFrame(columns=['sample', 'unit', 'fq1', 'fq2'], 
        index = range(int(len(fl))), dtype=str)
    i = 0
    while i < len(fl):
        df.loc[i]['sample'] = fl[i].split('_')[0]
        df.loc[i]['unit'] = fl[i].split('_')[1]
        df.loc[i]['fq1'] = fpl[i][:-3]
        df.loc[i]['fq2'] = np.nan
        i += 1

df.to_csv('units.tsv', sep='\t')
