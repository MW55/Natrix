import yaml
import os
import pandas as pd
import numpy as np

with open('test_data.yaml') as f_:
    config = yaml.load(f_)

file_list = os.listdir(config['general']['filename'])
file_list.sort()
df = pd.DataFrame(columns=['sample', 'unit', 'fq1', 'fq2'], 
        index =range(int(len(file_list)/2)))

if config['merge']['paired_End']:
    df = pd.DataFrame(columns=['sample', 'unit', 'fq1', 'fq2'], 
            index =range(int(len(file_list)/2)))
    i, j = (0, 0)
    while i < len(file_list)/2:
        df.loc[i] = file_list[j].split('_')[0], file_list[j].split('_')[1], \
            file_list[j], file_list[j+1]
        j += 2
        i += 1
else:
    df = pd.DataFrame(columns=['sample', 'unit', 'fq1', 'fq2'], 
            index = range(int(len(file_list))))
    i = 0
    while i < len(file_list):
        df.loc[i] = file_list[i].split('_')[0], file_list[i].split('_')[1], \
            file_list[i], np.nan
        i += 1

df.to_csv('units.tsv', sep='\t', index=False)
pd.DataFrame(df['sample']).to_csv('samples.tsv', sep='\t', index=False)
