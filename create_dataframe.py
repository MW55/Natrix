import pathlib

import yaml
import pandas as pd
import numpy as np
from glob import glob
import os
import sys

# Create the datatable containing the samples, units and paths of all
# fastq files formatted correctly. This is vital for the snakemake
# pipeline, without it, the wildcards can't be created.
# Additionally, options will be checked.

with open(sys.argv[1]) as f_:
    config = yaml.load(f_, Loader=yaml.FullLoader)


def create_dataframe(fl, fpl, config, slice):
    if config['merge']['paired_End'] and not config['general']['already_assembled']:
        df = pd.DataFrame(columns=['sample', 'unit', 'fq1', 'fq2'],
            index =range(int(len(fl)/2)), dtype=str)
        i, j = (0, 0)
        while i < len(fl)/2:
            # last split needs to be fwd or rev read
            # second last can be unit
            unit = fl[j].split('_')[-2]
            if unit in ['A', 'B']:
            	df.loc[i]['unit'] = unit
            	df.loc[i]['sample'] = '_'.join(fl[j].split('_')[:-2])
            else:
            	df.loc[i]['unit'] = ''
            	df.loc[i]['sample'] = '_'.join(fl[j].split('_')[:-1])
            df.loc[i]['fq1'] = fpl[j][:slice]
            df.loc[i]['fq2'] = fpl[j+1][:slice]
            j += 2
            i += 1
    else:
        df = pd.DataFrame(columns=['sample', 'unit', 'fq1', 'fq2'],
            index = range(int(len(fl))), dtype=str)
        i = 0
        while i < len(fl):
            # last split needs to be fwd or rev read
            # second last can be unit
            unit = fl[i].split('_')[-2]
            if unit in ['A', 'B']:
            	df.loc[i]['unit'] = unit
            	df.loc[i]['sample'] = '_'.join(fl[i].split('_')[:-2])
            else:
            	df.loc[i]['unit'] = ''
            	df.loc[i]['sample'] = '_'.join(fl[i].split('_')[:-1])
            df.loc[i]['fq1'] = fpl[i][:slice]
            df.loc[i]['fq2'] = np.nan
            i += 1
    return df


if __name__ == '__main__':
    # check config options
    if "-" in config["general"]["output_dir"]:
    	sys.exit("Please rename output folder, do not use a dash in the folder name")
    if config["classify"]["mothur"] and config["blast"]["blast"]:
        sys.exit("Please decide whether to use blast or classification with mothur. Both config options cannot be set to TRUE")
    if config["general"]["seq_rep"] == "ASV" and config["postcluster"]["mumu"]:
        print("Postclustering with mumu is not supported for ASVs.")
        changeopt = input("To proceed and set the mumu config option to FALSE, type 'yes': ")
        if changeopt == "yes":
            print("Proceeding")
            config["postcluster"]["mumu"] = False
        else:
            sys.exit("Workflow will be aborted")
    if not config['general']['already_assembled']:
        file_path_list = [os.path.join(config["general"]["output_dir"],'demultiplexed/' + name.split('/')[-1]) for name in
                          sorted(glob(config['general']['filename'].rstrip("/") + '/*.gz'))]
        file_list = sorted([file_.split('/')[-1] for file_
                    in file_path_list])
        slice = -3 # Remove the .gz extension from the file paths.
    else:
        file_path_list = sorted(glob(os.path.join(config["general"]["output_dir"],'assembly/*/*.fastq')))
        file_list = sorted([file_.split('/')[-1] for file_
                    in file_path_list])
        slice = None
    # create dataframe
    df = create_dataframe(file_list, file_path_list, config, slice)

    pathlib.Path(config["general"]["output_dir"]).mkdir(parents=True, exist_ok=True)
    df.to_csv(os.path.join(config["general"]["output_dir"],config["general"]['units']), sep='\t')
