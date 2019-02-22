import re
import yaml
import dinopy
import logging
import subprocess
import numpy as np
import pandas as pd
from glob import glob

# Script to assemble paired end reads using PandaSeq.
# If the primer were previously not removed, PandaSeq will remove them.
# If the sequences are single end, the primer will be cut off
# and the remaining length of the sequences will be compared to the
# cutoff values defined in the config.

primer_table = pd.read_csv(snakemake.input.primer_t, index_col='Probe',
        na_filter=False).to_dict('index')

if snakemake.params.paired_end:
    #old way using the old primertable
#    r1_primer = primer_table[snakemake.wildcards.sample + '_'
#            + snakemake.wildcards.unit]['specific_forward_primer']
#    r2_primer = primer_table[snakemake.wildcards.sample + '_'
#            + snakemake.wildcards.unit]['specific_reverse_primer']
# new way, now further down so that primertables without
# primer work (when everything is already removed)
#    r1_primer = primer_table[snakemake.wildcards.sample + '_'
#            + snakemake.wildcards.unit]['f_primer']
#    r2_primer = primer_table[snakemake.wildcards.sample + '_'
#            + snakemake.wildcards.unit]['r_primer']

    if snakemake.params.prim_rm:
        subprocess.call(['pandaseq',
            '-f',snakemake.input[0], '-r', snakemake.input[1], '-B', '-a', '-F',
            '-g', str(snakemake.log),
            '-w', str(snakemake.output), '-N',
            '-T', str(snakemake.threads),
            '-t', str(snakemake.params.threshold),
            '-o', str(snakemake.params.minoverlap),
            '-l', str(snakemake.params.minlen),
            '-L', str(snakemake.params.maxlen),
            '-C' 'min_phred:' + str(snakemake.params.minqual)])
    else:
        r1_primer = primer_table[snakemake.wildcards.sample + '_'
            + snakemake.wildcards.unit]['f_primer']
        r2_primer = primer_table[snakemake.wildcards.sample + '_'
            + snakemake.wildcards.unit]['r_primer']

        subprocess.call(['pandaseq',
            '-f', snakemake.input[0], '-r', snakemake.input[1], '-B', '-a', '-F',
            '-g', str(snakemake.log),
            '-w', str(snakemake.output),'-N',
            '-p', r1_primer, '-q', r2_primer,
            '-T', str(snakemake.threads),
            '-t', str(snakemake.params.threshold),
            '-o', str(snakemake.params.minoverlap),
            '-l', str(snakemake.params.minlen),
            '-L', str(snakemake.params.maxlen),
            '-C' 'min_phred:' + str(snakemake.params.minqual)])
else:
    logging.basicConfig(filename=str(snakemake.log),
            level = logging.DEBUG)
    iupac_dict_regex = {'M':'[AC]', 'R':'[AG]', 'W':'[AT]', 'S':'[CG]',
            'Y':'[CT]','K':'[GT]', 'V':'[ACG]', 'H':'[ACT]',
            'D':'[AGT]', 'B':'[CGT]', 'X':'[ACGT]', 'N':'[ACGT]'}

    # Function to remove the primer, barcode & polyN from the sequences
    # and compare the sequence length to the cutoffs defined in the config.
    def primer_len_filter(path, sample):
        sequence = dinopy.FastqReader(path)
        assembled = dinopy.FastqWriter(path.rsplit('_', 1)[0]
                + '_assembled.fastq')
        filt_out = dinopy.FastqWriter(path.rsplit('_', 1)[0]
                + '_filtered_out.fastq')
        assembled_counter = 0
        filt_out_counter = 0
        assembled.open()
        filt_out.open()
        for read in sequence.reads(quality_values=True):
            seq = check_for_match(read.sequence.decode(), sample)
            if seq[0] and snakemake.params.maxlen >= len(seq[1]) >= snakemake.params.minlen:
                assembled.write(seq[1].encode(), read.name, read.quality)
                assembled_counter += 1
            else:
                filt_out.write(read.sequence, read.name, read.quality)
                filt_out_counter += 1
        logging.info('{}: {} sequences were kept, \
                {} sequences were filtered out'.format(sample,
                    assembled_counter, filt_out_counter))
        assembled.close()
        filt_out.close()

    # Helper function for the IUPAC extended nucleotide base code
    def iupac_replace(sequence, iupac_dict):
        for i, j in iupac_dict_regex.items():
            sequence = sequence.replace(i, j)
        return sequence

    # Function to search and remove primer and barcode sequences.
    # It will also search for matching sequences two bases shifted
    # to the left and right from the position defined by the primertable,
    # to account for potential errors in the polyN column of the primertable.
    def check_for_match(sequence, sample):
        if snakemake.params.prim_rm:
            return (True, sequence)
        else:
            poly_prim_bar = [primer_table[sample][key] for key
                    in primer_table[sample].keys() if key in ['poly_N',
                        'specific_forward_primer', 'Barcode_forward']]
            prim_bar = re.compile(poly_prim_bar[1]
                    + iupac_replace(poly_prim_bar[2], iupac_dict_regex))
            for i in [0, 1, -1, 2, -2]:
                start = np.clip(len(primer_table[sample]['poly_N']) + i,
                        a_min=0, a_max=None)
                end = np.clip(len(''.join(poly_prim_bar)) + i, a_min=0,
                        a_max=None)
                if prim_bar.match(sequence[start : end]):
                    return (True, sequence.replace(sequence[ : end], ''))
                else:
                    return (False, sequence)

    primer_len_filter(snakemake.input[0],
            snakemake.input[0].split('/')[-1].rsplit('_', 1)[0])
