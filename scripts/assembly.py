import re
import yaml
import logging
import subprocess
import numpy as np
import pandas as pd
from glob import glob


# FASTQ file parser
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


# FASTQ file writer
def write_fastq(file, records):
    with open(file, 'w') as f:
        for header, sequence, quality in records:
            f.write(f"{header}\n{sequence}\n+\n{quality}\n")


# Read primer table
primer_table = pd.read_csv(snakemake.input.primer_t, index_col="Probe", na_filter=False).to_dict("index")

if snakemake.params.paired_end:
    if snakemake.params.prim_rm:
        subprocess.call([
            "pandaseq", "-f", snakemake.input[0], "-r", snakemake.input[1], "-B", "-a", "-F",
            "-g", str(snakemake.log), "-w", str(snakemake.output), "-N",
            "-T", str(snakemake.threads), "-t", str(snakemake.params.threshold),
            "-o", str(snakemake.params.minoverlap), "-l", str(snakemake.params.minlen),
            "-L", str(snakemake.params.maxlen), "-C", "min_phred:" + str(snakemake.params.minqual)
        ])
    else:
        r1_primer = primer_table[snakemake.wildcards.sample + "_" + snakemake.wildcards.unit]["specific_forward_primer"]
        r2_primer = primer_table[snakemake.wildcards.sample + "_" + snakemake.wildcards.unit]["specific_reverse_primer"]

        subprocess.call([
            "pandaseq", "-f", snakemake.input[0], "-r", snakemake.input[1], "-B", "-a", "-F",
            "-g", str(snakemake.log), "-w", str(snakemake.output), "-N",
            "-p", r1_primer, "-q", r2_primer,
            "-T", str(snakemake.threads), "-t", str(snakemake.params.threshold),
            "-o", str(snakemake.params.minoverlap), "-l", str(snakemake.params.minlen),
            "-L", str(snakemake.params.maxlen), "-C", "min_phred:" + str(snakemake.params.minqual)
        ])
else:
    logging.basicConfig(filename=str(snakemake.log), level=logging.DEBUG)
    iupac_dict_regex = {
        "M": "[AC]", "R": "[AG]", "W": "[AT]", "S": "[CG]", "Y": "[CT]", "K": "[GT]",
        "V": "[ACG]", "H": "[ACT]", "D": "[AGT]", "B": "[CGT]", "X": "[ACGT]", "N": "[ACGT]"
    }


    def iupac_replace(sequence, iupac_dict):
        for i, j in iupac_dict.items():
            sequence = sequence.replace(i, j)
        return sequence


    def check_for_match(sequence, sample):
        if snakemake.params.prim_rm:
            return True, sequence
        else:
            poly_prim_bar = [
                primer_table[sample][key] for key in primer_table[sample]
                if key in ["poly_N", "specific_forward_primer", "Barcode_forward"]
            ]
            prim_bar = re.compile(poly_prim_bar[1] + iupac_replace(poly_prim_bar[2], iupac_dict_regex))
            for i in [0, 1, -1, 2, -2]:
                start = np.clip(len(primer_table[sample]["poly_N"]) + i, a_min=0, a_max=None)
                end = np.clip(len("".join(poly_prim_bar)) + i, a_min=0, a_max=None)
                if prim_bar.match(sequence[start:end]):
                    return True, sequence.replace(sequence[:end], "")
            return False, sequence


    def primer_len_filter(path, sample):
        assembled_records = []
        filt_out_records = []
        assembled_counter = 0
        filt_out_counter = 0

        for header, sequence, quality in read_fastq(path):
            seq = check_for_match(sequence, sample)
            if seq[0] and snakemake.params.maxlen >= len(seq[1]) >= snakemake.params.minlen:
                assembled_records.append((header, seq[1], quality))
                assembled_counter += 1
            else:
                filt_out_records.append((header, sequence, quality))
                filt_out_counter += 1

        write_fastq(path.rsplit("_", 1)[0] + "_assembled.fastq", assembled_records)
        write_fastq(path.rsplit("_", 1)[0] + "_filtered_out.fastq", filt_out_records)

        logging.info(
            f"{sample}: {assembled_counter} sequences were kept, {filt_out_counter} sequences were filtered out")


    primer_len_filter(snakemake.input[0], snakemake.input[0].split("/")[-1].rsplit("_", 1)[0])
