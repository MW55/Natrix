import pandas as pd
import subprocess
import dinopy
import yaml
import re
import numpy as np
from glob import glob
import logging

primer_table = pd.read_csv(snakemake.input.primer_t, index_col='Probe').to_dict('index')

def primer_len_filter(path, sample):
    sequence = dinopy.FastqReader(path)
    assembled = dinopy.FastqWriter(path.rsplit('_', 1)[0] + '_assembled.fastq')
    filt_out = dinopy.FastqWriter(path.rsplit('_', 1)[0] + '_filtered_out.fastq')
    assembled_counter = 0
    filt_out_counter = 0
    #open all writers
    assembled.open()
    filt_out.open()
    #check for match, remove primer,barcode & polyN and compare the sequence length to the cutoffs defined in the config
    for read in sequence.reads(quality_values=True):
        seq = check_for_match(read.sequence.decode(), sample)
        if seq[0] and snakemake.params.maxlen >= len(seq[1]) >= snakemake.params.minlen:
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
    if snakemake.params.prim_rm:
        return (True, sequence)
    else:
        poly_prim_bar = [primer_table[sample][key] for key in primer_table[sample].keys() if key in ['poly_N', 'specific_forward_primer', 'Barcode_forward']]
        prim_bar = re.compile(poly_prim_bar[1] + iupac_replace(poly_prim_bar[2], iupac_dict_regex))
        for i in [0, 1, -1, 2, -2]:
            start = np.clip(len(primer_table[sample]['poly_N']) + i, a_min=0, a_max=None)
            end = np.clip(len(''.join(poly_prim_bar)) + i, a_min=0, a_max=None)
            if prim_bar.match(sequence[start : end]):
                return (True, sequence.replace(sequence[ : end], ''))
            else:
                return (False, sequence)

if snakemake.params.paired_end:
    r1_primer = primer_table[snakemake.wildcards.sample + '_' + snakemake.wildcards.unit]['specific_forward_primer']
    r2_primer = primer_table[snakemake.wildcards.sample + '_' + snakemake.wildcards.unit]['specific_reverse_primer']
    if snakemake.params.prim_rm:
        subprocess.call('pandaseq -f {r1} -r {r2} -B -a -F -g {log} -w {output} -N -t  {threshold}, -o {minoverlap} -l {minlen} -L {maxlen} -C min_phred:{minqual}'.format(r1=snakemake.input[0],r2=snakemake.input[1],log=snakemake.log, output=snakemake.output, r1_primer=r1_primer,r2_primer=r2_primer,threshold=snakemake.params.threshold, minoverlap=snakemake.params.minoverlap, minlen=snakemake.params.minlen, maxlen=snakemake.params.maxlen, minqual=snakemake.params.minqual), shell = True)
    else:
        subprocess.call('pandaseq -f {r1} -r {r2} -B -a -F -g {log} -w {output} -N -p "{r1_primer}" -q "{r2_primer}" -t  {threshold}, -o {minoverlap} -l {minlen} -L {maxlen} -C min_phred:{minqual}'.format(r1=snakemake.input[0],r2=snakemake.input[1],log=snakemake.log, output=snakemake.output, r1_primer=r1_primer,r2_primer=r2_primer,threshold=snakemake.params.threshold, minoverlap=snakemake.params.minoverlap, minlen=snakemake.params.minlen, maxlen=snakemake.params.maxlen, minqual=snakemake.params.minqual), shell = True)
else:
    logging.basicConfig(filename=str(snakemake.log), level = logging.DEBUG)
    iupac_dict_regex = {'M':'[AC]', 'R':'[AG]', 'W':'[AT]', 'S':'[CG]', 'Y':'[CT]',
        'K':'[GT]', 'V':'[ACG]', 'H':'[ACT]', 'D':'[AGT]', 'B':'[CGT]', 'X':'[ACGT]', 'N':'[ACGT]'}

    def primer_len_filter(path, sample):
        sequence = dinopy.FastqReader(path)
        assembled = dinopy.FastqWriter(path.rsplit('_', 1)[0] + '_assembled.fastq')
        filt_out = dinopy.FastqWriter(path.rsplit('_', 1)[0] + '_filtered_out.fastq')
        assembled_counter = 0
        filt_out_counter = 0
        #open all writers
        assembled.open()
        filt_out.open()
        #check for match, remove primer,barcode & polyN and compare the sequence length to the cutoffs defined in the config
        for read in sequence.reads(quality_values=True):
            seq = check_for_match(read.sequence.decode(), sample)
            if seq[0] and snakemake.params.maxlen >= len(seq[1]) >= snakemake.params.minlen:
                assembled.write(seq[1].encode(), read.name, read.quality)
                assembled_counter += 1
            else:
                filt_out.write(read.sequence, read.name, read.quality)
                filt_out_counter += 1
        #this goes into the logfile
        logging.info('{}: {} sequences were kept, {} sequences were filtered out'.format(sample, assembled_counter, filt_out_counter))
        # close all writers
        assembled.close()
        filt_out.close()

    # helper function
    def iupac_replace(sequence, iupac_dict):
        for i, j in iupac_dict_regex.items():
            sequence = sequence.replace(i, j)
        return sequence

    def check_for_match(sequence, sample):
        if snakemake.params.prim_rm:
            return (True, sequence)
        else:
            poly_prim_bar = [primer_table[sample][key] for key in primer_table[sample].keys() if key in ['poly_N', 'specific_forward_primer', 'Barcode_forward']]
            prim_bar = re.compile(poly_prim_bar[1] + iupac_replace(poly_prim_bar[2], iupac_dict_regex))
            for i in [0, 1, -1, 2, -2]:
                start = np.clip(len(primer_table[sample]['poly_N']) + i, a_min=0, a_max=None)
                end = np.clip(len(''.join(poly_prim_bar)) + i, a_min=0, a_max=None)
                if prim_bar.match(sequence[start : end]):
                    return (True, sequence.replace(sequence[ : end], ''))
                else:
                    return (False, sequence)

    primer_len_filter(snakemake.input[0], snakemake.input[0].split('/')[-1].rsplit('_', 1)[0])


