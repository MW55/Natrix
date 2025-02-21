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

# FASTQ parser
def read_fastq(file):
    with open(file, 'r') as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            sequence = f.readline().strip()
            plus = f.readline().strip()
            quality = f.readline().strip()
            yield header, sequence, quality

# FASTQ writer
def write_fastq(file, records):
    with open(file, 'w') as f:
        for header, sequence, quality in records:
            f.write(f"{header}\n{sequence}\n+\n{quality}\n")

p_table = pd.read_csv(snakemake.params.primertable, index_col='Probe')
primertable = p_table.to_dict('index')
data_folder = str(snakemake.params.filename)
file_path_list = sorted(glob.glob(data_folder + "/*.fast*"))

iupac_dict_regex = {'M':'[AC]', 'R':'[AG]', 'W':'[AT]', 'S':'[CG]', 'Y':'[CT]',
                    'K':'[GT]', 'V':'[ACG]', 'H':'[ACT]', 'D':'[AGT]',
                    'B':'[CGT]', 'X':'[ACGT]', 'N':'[ACGT]'}

def iupac_replace(sequence, iupac_dict):
    for i, j in iupac_dict.items():
        sequence = sequence.replace(i, j)
    return sequence

def define_direction_demulti(polyN, prim, barcode):
    def check_for_match_demulti(sequence, sample):
        poly_prim_bar = [primertable[sample][key] for key in primertable[sample].keys() if key in [polyN, prim, barcode]]
        prim_bar = re.compile(poly_prim_bar[1] + iupac_replace(poly_prim_bar[2], iupac_dict_regex))
        for i in [0, 1, -1, 2, -2]:
            start = np.clip(len(primertable[sample][polyN]) + i, a_min=0, a_max=None)
            end = np.clip(len(''.join(poly_prim_bar)) + i, a_min=0, a_max=None)
            if prim_bar.match(sequence[start : end]):
                return True
        return False
    return check_for_match_demulti

check_for_match_fwd_demulti = define_direction_demulti('poly_N', 'specific_forward_primer', 'Barcode_forward')
check_for_match_rev_demulti = define_direction_demulti('poly_N_rev', 'specific_reverse_primer', 'Barcode_reverse')

def define_direction_sort(prim):
    def check_for_match_sort(sequence, sample):
        prim_regex = re.compile(iupac_replace(primertable[sample][prim], iupac_dict_regex))
        return bool(prim_regex.match(sequence[:len(primertable[sample][prim])]))
    return check_for_match_sort

check_for_match_sort_fwd = define_direction_sort('specific_forward_primer')
check_for_match_sort_rev = define_direction_sort('specific_reverse_primer')

def demultiplexer(file_path_list):
    for sample in file_path_list:
        output_records = {name: [] for name in primertable.keys()}
        for header, sequence, quality in read_fastq(sample):
            for sample_name in primertable.keys():
                if check_for_match_fwd_demulti(sequence, sample_name):
                    output_records[sample_name + '_R1'].append((header, sequence, quality))
                elif check_for_match_rev_demulti(sequence, sample_name):
                    output_records[sample_name + '_R2'].append((header, sequence, quality))
        for name, records in output_records.items():
            write_fastq(f'demultiplexed/{name}.fastq.gz', records)

def already_assembled(primertable, file_path_list):
    for f_ in file_path_list:
        rm_unzipped = False
        if '.gz' in f_:
            subprocess.run(['gunzip', '-k', f_])
            f_ = f_.split('.gz')[0]
            rm_unzipped = True
        for sample in primertable.keys():
            pathlib.Path(f'results/assembly/{sample}').mkdir(parents=True, exist_ok=True)
            if sample in f_:
                shutil.copy(f_, f'results/assembly/{sample}/{sample}_assembled.fastq')
        if rm_unzipped:
            os.remove(f_)

if snakemake.params.demultiplexing:
    demultiplexer(file_path_list)
elif snakemake.params.assembled:
    already_assembled(primertable, file_path_list)
else:
    for file in file_path_list:
        shutil.copy(file, 'demultiplexed/')
