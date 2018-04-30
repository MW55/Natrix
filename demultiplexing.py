import pandas as pd
import numpy as np
import dinopy
import os
import re
import yaml
import glob

#IMPORTANT: The name formatting differs between the new and the old style: now an additional integer itentifier is right
# next to the unit identifier, have to ask is that is the planned name formatting and if so, adjust the string formatting 
# during the sample table generation so that the integer identifier is in the sample column, not the unit column

# config and p_table has to be adjusted for the corresponding path to the primertable with config config['general']['filename'] + '.csv'
with open('../../workflow_test/amplicon_snakemake_pipeline/demultiplexer/test_data.yaml') as f_:
    config = yaml.load(f_)
p_table = pd.read_csv('../../workflow_test/amplicon_snakemake_pipeline/demultiplexer/test_data.csv', index_col='Probe')
primertable = p_table.to_dict('index')
data_folder = f_.name.rsplit('/', 1)[0] + '/' + config['general']['filename']
file_path_list = sorted(glob.glob(data_folder + "/*.gz"))

### Demultiplexing part of the script, this is only needed for the new samples with the triplet barcode and the "manual" pooling in-lab ###
### might be better to put it in a seperate script ###

# regex substitution dict
iupac_dict_regex = {'M':'[AC]', 'R':'[AG]', 'W':'[AT]', 'S':'[CG]', 'Y':'[CT]',
'K':'[GT]', 'V':'[ACG]', 'H':'[ACT]', 'D':'[AGT]', 'B':'[CGT]', 'X':'[ACGT]', 'N':'[ACGT]'}

# substitution helper function
def iupac_replace(sequence, iupac_dict):
    for i, j in iupac_dict_regex.items():
        sequence = sequence.replace(i, j)
    return sequence

# Higher order function (define_direction) which returns a closure (a function together with an environment, here check_for_match),
# which has access to the environment of the higher order function, this allows creation of similar functions
# while limiting code redundancy
def define_direction(polyN, prim, barcode):
    def check_for_match(sequence, sample):
        poly_prim_bar = [primertable[sample][key] for key in primertable[sample].keys() if key in [polyN, prim, barcode]]
        prim_bar = re.compile(poly_prim_bar[1] + iupac_replace(poly_prim_bar[2], iupac_dict_regex))
        for i in [0, 1, -1, 2, -2]:
            start = np.clip(len(primertable[sample][polyN]) + i, a_min=0, a_max=None)
            end = np.clip(len(''.join(poly_prim_bar)) + i, a_min=0, a_max=None)
            if prim_bar.match(sequence[start : end]):
                return True
        else:
            return False
    return check_for_match

# These are the variations of the check_for_match closure for the forward and reverse primer, the arguments are the column indices
# of the corresponding polyN col, the col after the primer (as the slice beginning is inclusive, end is exclusive) and the 
# col index of the barcode
check_for_match_fwd = define_direction('poly_N', 'specific_forward_primer', 'Barcode_forward')
check_for_match_rev = define_direction('poly_N_rev', 'specific_reverse_primer', 'Barcode_reverse')

# Dinopy
def demultiplexer(file_path_list):
    samples = []
    output_filepaths = []
    os.mkdir('demultiplexed')
    for sample in primertable.keys():
        samples.append(sample + '_R1')
        samples.append(sample + '_R2')
        output_filepaths.append('demultiplexed/' + sample + '_R1.fastq.gz')
        output_filepaths.append('demultiplexed/' + sample + '_R2.fastq.gz')

    # create a dict of writers
    writers = {name: dinopy.FastqWriter(path) for name, path in 
               zip(samples, output_filepaths)}
    #open all writers
    for writer in writers.values():
        writer.open()
    
    #start writing
    for sample in file_path_list:
        sequence = dinopy.FastqReader(sample)
        for read in sequence.reads(quality_values=True):
            for sample in primertable.keys():
                if check_for_match_fwd(read.sequence.decode(), sample):
                    writers[sample + '_R1'].write(read.sequence, read.name, read.quality)
                elif check_for_match_rev(read.sequence.decode(), sample):
                    writers[sample + '_R2'].write(read.sequence, read.name, read.quality)
                else:
                    pass
            
    # close all writers
    for writer in writers.values():
        writer.close()

# Run the demultiplexing script
demultiplexer(file_path_list)

### Create the datatable containing the samples, units and paths of all fastq files formatted correctly ###
### this is vital for the snakemake pipeline, without it, the wildcards can't be created ###

# might be better to glob abspaths

file_path_list = sorted(glob.glob('demultiplexed/' + "/*.gz"))
file_list = sorted([file_.split('/')[-1] for file_ in glob.glob('demultiplexed/' + "/*.gz")])
df = pd.DataFrame(columns=['sample', 'unit', 'fq1', 'fq2'],
    index =range(int(len(file_list)/2)))
if config['merge']['paired_End']:
    df = pd.DataFrame(columns=['sample', 'unit', 'fq1', 'fq2'],
        index =range(int(len(file_list)/2)))
    i, j = (0, 0)
    while i < len(file_list)/2:
        df.loc[i] = file_list[j].split('_')[0], file_list[j].split('_')[1], \
            file_path_list[j][:-3], file_path_list[j+1][:-3]
        j += 2
        i += 1
    reads = [1,2]
else:
    df = pd.DataFrame(columns=['sample', 'unit', 'fq1', 'fq2'], 
        index = range(int(len(file_list))))
    i = 0
    while i < len(file_list):
        df.loc[i] = file_list[i].split('_')[0], file_list[i].split('_')[1], \
            file_list[i][:-3], np.nan
        i += 1
    reads = 1
df.to_csv('units.tsv', sep='\t')
pd.DataFrame(df['sample']).to_csv('samples.tsv', sep='\t')
samples = pd.read_table(config["general"]["samples"], index_col="sample")
units = pd.read_table(config["general"]["units"], index_col=["sample", "unit"],
    dtype=str)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])
