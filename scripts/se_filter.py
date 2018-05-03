import dinopy
import pandas as pd
import yaml
import re
import numpy as np
from glob import glob

primertable = pd.read_csv('../../workflow_test/amplicon_snakemake_pipeline/fwd_read_fake/fwd_read_fake.csv', index_col='Probe').to_dict('index')
iupac_dict_regex = {'M':'[AC]', 'R':'[AG]', 'W':'[AT]', 'S':'[CG]', 'Y':'[CT]',
        'K':'[GT]', 'V':'[ACG]', 'H':'[ACT]', 'D':'[AGT]', 'B':'[CGT]', 'X':'[ACGT]', 'N':'[ACGT]'}

def primer_len_filter(path, sample):
        sequence = dinopy.FastqReader(path)
        assembled = dinopy.FastqWriter(path + '_assembled.fastq')
        filt_out = dinopy.FastqWriter(path + '_filtered_out.fastq')
        assembled_counter = 0
        filt_out_counter = 0
        #open all writers
        assembled.open()
        filt_out.open()
        #check for match, remove primer,barcode & polyN and compare the sequence length to the cutoffs defined in the config
        for read in sequence.reads(quality_values=True):
            seq = check_for_match(read.sequence.decode(), sample)
            if seq[0] and maxlen >= len(seq[1]) >= minlen:
                assembled.write(seq[1].encode(), read.name, read.quality)
                assembled_counter += 1
            else:
                filt_out.write(read.sequence, read.name, read.quality)
                filt_out_counter += 1
        #this goes into the logfile
        print('{}: {} sequences were kept, {} sequences were filtered out'.format(sample, assembled_counter, filt_out_counter))
        # close all writers
        assembled.close()
        filt_out.close()

# helper function
def iupac_replace(sequence, iupac_dict):
    for i, j in iupac_dict_regex.items():
        sequence = sequence.replace(i, j)
    return sequence

def check_for_match(sequence, sample):
    if prim_rm:
        return (True, sequence)
    else:
        poly_prim_bar = [primertable[sample][key] for key in primertable[sample].keys() if key in ['poly_N', 'specific_forward_primer', 'Barcode_forward']]
        prim_bar = re.compile(poly_prim_bar[1] + iupac_replace(poly_prim_bar[2], iupac_dict_regex))
        for i in [0, 1, -1, 2, -2]:
            start = np.clip(len(primertable[sample]['poly_N']) + i, a_min=0, a_max=None)
            end = np.clip(len(''.join(poly_prim_bar)) + i, a_min=0, a_max=None)
            if prim_bar.match(sequence[start : end]):
                return (True, sequence.replace(sequence[ : end], ''))
            else:
                return (False, sequence)


primer_len_filter(sample, sample.split('/')[-1].rsplit('_', 1)[0]):
