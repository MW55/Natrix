import re
import yaml
import dinopy
import logging
import subprocess
import numpy as np
import pandas as pd
from glob import glob

primer_table = pd.read_csv(snakemake.input.primer_t, index_col="Probe",
        na_filter=False).to_dict("index")

if snakemake.params.paired_end:
    if snakemake.params.prim_rm:
        subprocess.call(["mv", snakemake.input[0], snakemake.output[0]])
        subprocess.call(["mv", snakemake.input[1], snakemake.output[1]])
    else:
        r1_primer = primer_table[snakemake.wildcards.sample + "_"
            + snakemake.wildcards.unit]["f_primer"]
        r2_primer = primer_table[snakemake.wildcards.sample + "_"
            + snakemake.wildcards.unit]["r_primer"]

        subprocess.call(["cutadapt", "-g", r1_primer, "-G", r2_primer, "-o",
                         snakemake.output[0], "-p", snakemake.output[1], snakemake.input[0], snakemake.output[1]])
else:
    r1_primer = primer_table[snakemake.wildcards.sample + "_"
                             + snakemake.wildcards.unit]["f_primer"]
    subprocess.call(["cutadapt", "-g", r1_primer, "-o", snakemake.output[0], snakemake.input[0]])