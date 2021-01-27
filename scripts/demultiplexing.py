import pandas as pd
import numpy as np
import dinopy
import os
import re
import yaml
import glob
import sys
import shutil
import subprocess
import pathlib
# This script is used for moving already assembled fastq files to the correct folder for further processing,
# demultiplexing samples that were pooled in-lab (instead of at the sequencing company) and
# to manually sort reads depending on their primersequence.

p_table = pd.read_csv(snakemake.params.primertable, index_col='Probe')
primertable = p_table.to_dict('index')
data_folder = str(snakemake.params.filename)
file_path_list = sorted(glob.glob(data_folder + "/*.fast*"))

# Regex substitution dict.
iupac_dict_regex = {'M':'[AC]', 'R':'[AG]', 'W':'[AT]', 'S':'[CG]', 'Y':'[CT]',
                    'K':'[GT]', 'V':'[ACG]', 'H':'[ACT]', 'D':'[AGT]',
                    'B':'[CGT]', 'X':'[ACGT]', 'N':'[ACGT]'}


# Substitution helper function.
def iupac_replace(sequence, iupac_dict):
    for i, j in iupac_dict_regex.items():
        sequence = sequence.replace(i, j)
    return sequence


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

    # Create a dict of writers.
    writers = {name: dinopy.FastqWriter(path) for name, path in
               zip(samples, output_filepaths)}

    # Open all writers.
    for writer in writers.values():
        writer.open()

    # Start writing.
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

    # Close all writers.
    for writer in writers.values():
        writer.close()

def read_sorter(primertable):
    if not os.path.exists('demultiplexed/not_sorted'):
        os.mkdir('demultiplexed/not_sorted')
    samples = []
    output_filepaths = []
    for sample in primertable.keys():
        samples.append(sample + snakemake.params.name_ext[:-1] + '1')
        samples.append(sample + snakemake.params.name_ext[:-1] + '2')
        samples.append(sample + '_not_sorted')
        output_filepaths.append('demultiplexed/' + sample + '_R1.fastq.gz')
        output_filepaths.append('demultiplexed/' + sample + '_R2.fastq.gz')
        output_filepaths.append('demultiplexed/not_sorted/' + sample + '_not_sorted.fastq.gz')

    # Create a dict of writers.
    writers = {name: dinopy.FastqWriter(path) for name, path in
               zip(samples, output_filepaths)}

    # Open all writers.
    for writer in writers.values():
        writer.open()

    # Start writing.
    for sample in primertable.keys():
        fwd = dinopy.FastqReader('../' + data_folder + '/' + sample + str(snakemake.params.name_ext)[:-1] + '1.fastq.gz')
        rev = dinopy.FastqReader('../' + data_folder + '/' + sample + str(snakemake.params.name_ext)[:-1] + '2.fastq.gz')
        for read_f, read_r in zip(fwd.reads(quality_values=True), rev.reads(quality_values=True)):
            if check_for_match_sort_fwd(read_f.sequence.decode(),
                sample.split('/')[-1]) and check_for_match_sort_rev(read_r.sequence.decode(),
                    sample.split('/')[-1]):
                writers[sample + '_R1'].write(read_f.sequence, read_f.name,
                                read_f.quality)
                writers[sample + '_R2'].write(read_r.sequence, read_r.name,
                                read_r.quality)
            elif check_for_match_sort_rev(read_f.sequence.decode(),
                sample.split('/')[-1]) and check_for_match_sort_fwd(read_r.sequence.decode(),
                    sample.split('/')[-1]):
                writers[sample + '_R2'].write(read_f.sequence, read_f.name,
                                read_f.quality)
                writers[sample + '_R1'].write(read_r.sequence, read_r.name,
                                read_r.quality)
            else:
                writers[sample + '_not_sorted'].write(read_f.sequence, read_f.name,
                                read_f.quality)
                writers[sample + '_not_sorted'].write(read_r.sequence, read_r.name,
                                read_r.quality)

    # Close all writers.
    for writer in writers.values():
        writer.close()


def already_assembled(primertable, file_path_list):
    for f_ in file_path_list:
        if '.gz' in f_:
            subprocess.run(['gunzip', f_])
            f_ = f_.split('.gz')[0]
        for sample in primertable.keys():
            pathlib.Path('../results/assembly/' + sample).mkdir(parents=True, exist_ok=True)
            if sample in f_:
                shutil.copy(f_, '../results/assembly/' + sample + '/' + sample + '_assembled.fastq')
                shutil.copy(f_, '../demultiplexed/')


# Run the demultiplexing / read sorting script.
if snakemake.params.demultiplexing:
    print('1')
    demultiplexer(file_path_list)
elif snakemake.params.read_sorting:
    read_sorter(primertable)
    print('2')
elif snakemake.params.assembled:
    already_assembled(primertable, file_path_list)
    print('3')
else:
    print('4')
    # If the files do not need demultiplexing / sorting, just copy them to
    # the demultiplexed folder. Leave original files in input folder
    for file in file_path_list:
        print(file)
        shutil.copy(file, 'demultiplexed/')
