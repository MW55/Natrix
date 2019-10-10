import yaml
import pandas as pd
import numpy as np
from glob import glob
import sys

# Create the datatable containing the samples, units and paths of all
# fastq files formatted correctly. This is vital for the snakemake
# pipeline, without it, the wildcards can't be created.

with open(sys.argv[1]) as f_:
    config = yaml.load(f_, Loader=yaml.FullLoader)

def create_dataframe(fl, fpl, config, slice):
    if config['merge']['paired_End'] and not config['general']['already_assembled']:
        df = pd.DataFrame(columns=['sample', 'unit', 'fq1', 'fq2'],
            index =range(int(len(fl)/2)), dtype=str)
        i, j = (0, 0)
        while i < len(fl)/2:
            df.loc[i]['sample'] = fl[j].split('_')[0]
            df.loc[i]['unit'] = fl[j].split('_')[1]
            df.loc[i]['fq1'] = fpl[j][:slice]
            df.loc[i]['fq2'] = fpl[j+1][:slice]
            j += 2
            i += 1
    else:
        df = pd.DataFrame(columns=['sample', 'unit', 'fq1', 'fq2'], 
            index = range(int(len(fl))), dtype=str)
        i = 0
        while i < len(fl):
            df.loc[i]['sample'] = fl[i].split('_')[0]
            df.loc[i]['unit'] = fl[i].split('_')[1]
            df.loc[i]['fq1'] = fpl[i][:slice]
            df.loc[i]['fq2'] = np.nan
            i += 1
    return df

if __name__ == '__main__':
    if not config['general']['already_assembled']:
        file_path_list = sorted(glob('demultiplexed/*.gz'))
        file_list = sorted([file_.split('/')[-1] for file_ 
                    in file_path_list])
        slice = -3 # Remove the .gz extension from the file paths.
    else:
        file_path_list = sorted(glob('results/assembly/*/*.fastq'))
        file_list = sorted([file_.split('/')[-1] for file_ 
                    in file_path_list])
        slice = None
    df = create_dataframe(file_list, file_path_list, config, slice)
    df.to_csv('units.tsv', sep='\t')
