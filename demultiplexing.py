import pandas as pd
import numpy as np
import dinopy
import os
import re
import yaml
import glob
import sys
import shutil
# IMPORTANT: The name formatting differs between the new and the old
# style: now an additional integer itentifier is right next to the
# unit identifier, have to ask is that is the planned name formatting
# and if so, adjust the string formatting during the sample table
# generation so that the integer identifier is in the sample column,
# not the unit column

# config and p_table has to be adjusted for the corresponding path to the 
# primertable with config config['general']['filename'] + '.csv'
with open(sys.argv[1] + '.yaml') as f_:
    config = yaml.load(f_)
p_table = pd.read_csv(sys.argv[1] + '.csv', index_col='Probe') #test_data.csv
primertable = p_table.to_dict('index')
data_folder = config['general']['filename'] #f_.name.rsplit('/', 1)[0] + '/' + 
file_path_list = sorted(glob.glob(data_folder + "/*.fastq*"))

# Demultiplexing part of the script, this is only needed for the new
# samples with the triplet barcode and the "manual" pooling in-lab

# regex substitution dict
iupac_dict_regex = {'M':'[AC]', 'R':'[AG]', 'W':'[AT]', 'S':'[CG]', 'Y':'[CT]',
                    'K':'[GT]', 'V':'[ACG]', 'H':'[ACT]', 'D':'[AGT]',
                    'B':'[CGT]', 'X':'[ACGT]', 'N':'[ACGT]'}

# substitution helper function
def iupac_replace(sequence, iupac_dict):
    for i, j in iupac_dict_regex.items():
        sequence = sequence.replace(i, j)
    return sequence

# Higher order function (define_direction) which returns a closure
# (a function together with an environment, here check_for_match),
# which has access to the environment of the higher order function,
# this allows creation of similar functions while limiting code
# redundancy
def define_direction_demulti(polyN, prim, barcode):
    def check_for_match_demulti(sequence, sample):
        poly_prim_bar = [primertable[sample][key] for key
                        in primertable[sample].keys() if key
                        in [polyN, prim, barcode]]
        prim_bar = re.compile(poly_prim_bar[1] + iupac_replace(poly_prim_bar[2],
                                iupac_dict_regex))
        for i in [0, 1, -1, 2, -2]:
            start = np.clip(len(primertable[sample][polyN]) + i, a_min=0,
                                a_max=None)
            end = np.clip(len(''.join(poly_prim_bar)) + i, a_min=0,
                                a_max=None)
            if prim_bar.match(sequence[start : end]):
                return True
        else:
            return False
    return check_for_match_demulti

# These are the variations of the check_for_match closure for the
# forward and reverse primer, the arguments are the column indices
# of the corresponding polyN col, the col after the primer
# (as the slice beginning is inclusive, end is exclusive) and the
# col index of the barcode.
check_for_match_fwd_demulti = define_direction_demulti('poly_N', 'specific_forward_primer',
                                       'Barcode_forward')
check_for_match_rev_demulti = define_direction_demulti('poly_N_rev', 'specific_reverse_primer',
                                       'Barcode_reverse')

def define_direction_sort(prim):
    def check_for_match_sort(sequence, sample):
        prim_regex = re.compile(iupac_replace(primertable[sample][prim],
                                              iupac_dict_regex))
        if prim_regex.match(sequence[:len(primertable[sample][prim])]):
            return True
        else:
            return False
    return check_for_match_sort

check_for_match_sort_fwd = define_direction_sort('specific_forward_primer')
check_for_match_sort_rev = define_direction_sort('specific_reverse_primer')

# Create a dict of Dinopy writer instances and write the sequences
# according to their barcode and primer sequence in the corresponding
# files defined in the primertable.
def demultiplexer(file_path_list):
    samples = []
    output_filepaths = []
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
                if check_for_match_fwd_demulti(read.sequence.decode(), sample):
                    writers[sample + '_R1'].write(read.sequence, read.name,
                            read.quality)
                elif check_for_match_rev_demulti(read.sequence.decode(), sample):
                    writers[sample + '_R2'].write(read.sequence, read.name,
                            read.quality)
                else:
                    pass

    # close all writers
    for writer in writers.values():
        writer.close()

if not os.path.exists('demultiplexed'):
    os.mkdir('demultiplexed')

if not os.path.exists('demultiplexed/not_sorted'):
    os.mkdir('demultiplexed/not_sorted')

def read_sorter(primertable):
    samples = []
    output_filepaths = []
    for sample in primertable.keys():
        samples.append(sample + config['merge']['name_ext'][:-1] + '1')
        samples.append(sample + config['merge']['name_ext'][:-1] + '2')
        samples.append(sample + '_not_sorted')
        output_filepaths.append('demultiplexed/' + sample + '_R1.fastq.gz')
        output_filepaths.append('demultiplexed/' + sample + '_R2.fastq.gz')
        output_filepaths.append('demultiplexed/not_sorted/' + sample + '_not_sorted.fastq.gz')

    # create a dict of writers
    writers = {name: dinopy.FastqWriter(path) for name, path in 
               zip(samples, output_filepaths)}

    #open all writers
    for writer in writers.values():
        writer.open()

    #start writing
    for sample in primertable.keys():
        fwd = dinopy.FastqReader(data_folder + '/' + sample + config['merge']['name_ext'][:-1] + '1.fastq.gz')
        rev = dinopy.FastqReader(data_folder + '/' + sample + config['merge']['name_ext'][:-1] + '2.fastq.gz')
        for read_f, read_r in zip(fwd.reads(quality_values=True),rev.reads(quality_values=True)):
            if check_for_match_sort_fwd(read_f.sequence.decode(), sample.split('/')[-1]) and check_for_match_sort_rev(read_r.sequence.decode(), sample.split('/')[-1]):
                writers[sample + '_R1'].write(read_f.sequence, read_f.name,
                                read_f.quality)
                writers[sample + '_R2'].write(read_r.sequence, read_r.name,
                                read_r.quality)
            elif check_for_match_sort_rev(read_f.sequence.decode(), sample.split('/')[-1]) and check_for_match_sort_fwd(read_r.sequence.decode(), sample.split('/')[-1]):
                writers[sample + '_R2'].write(read_f.sequence, read_f.name,
                                read_f.quality)
                writers[sample + '_R1'].write(read_r.sequence, read_r.name,
                                read_r.quality)
            else:
                writers[sample + '_not_sorted'].write(read_f.sequence, read_f.name,
                                read_f.quality)
                writers[sample + '_not_sorted'].write(read_r.sequence, read_r.name,
                                read_r.quality)

    # close all writers
    for writer in writers.values():
        writer.close()

# If the files do not need demultiplexing, just move them to
# the demultiplexed folder.
if config['general']['demultiplexing']:
    demultiplexer(file_path_list)
elif config['general']['read_sorting']:
    read_sorter(primertable)
else:
    # Run the demultiplexing script
    for file in file_path_list:
        shutil.move(file, 'demultiplexed/')
