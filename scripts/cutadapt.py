import subprocess
import pandas as pd

primer_table = pd.read_csv(snakemake.input.primer_t, index_col="Probe",
        na_filter=False).to_dict("index")

r1_barcode = primer_table[snakemake.wildcards.sample + "_"
+ snakemake.wildcards.unit]["Barcode_forward"]
r2_barcode = primer_table[snakemake.wildcards.sample + "_"
             + snakemake.wildcards.unit]["Barcode_reverse"]
r1_primer = primer_table[snakemake.wildcards.sample + "_"
                         + snakemake.wildcards.unit]["specific_forward_primer"]
r2_primer = primer_table[snakemake.wildcards.sample + "_"
                         + snakemake.wildcards.unit]["specific_reverse_primer"]


if snakemake.params.paired_end:
    if snakemake.params.prim_rm:
        subprocess.call(["mv", snakemake.input[0], snakemake.output[0]])
        subprocess.call(["mv", snakemake.input[1], snakemake.output[1]])
    elif not snakemake.params.bar_removed:
        subprocess.call(["cutadapt","-j 0", "-m", str(snakemake.params.minlen), "-M", str(snakemake.params.maxlen),
                         "-g", r1_primer + r1_barcode, "-G", r2_primer + r2_barcode, "-o",
                         snakemake.output[0], "-p", snakemake.output[1], snakemake.input[0], snakemake.input[1]])
    else:
        subprocess.call(["cutadapt","-j 0", "-m", str(snakemake.params.minlen), "-M", str(snakemake.params.maxlen),
                         "-g", r1_primer, "-G", r2_primer, "-o", snakemake.output[0], "-p",
                         snakemake.output[1], snakemake.input[0], snakemake.input[1]])
else:
    if not snakemake.params.bar_removed:
        subprocess.call(["cutadapt","-j 0", "-m", str(snakemake.params.minlen), "-M", str(snakemake.params.maxlen),
                         "-g", r1_primer + r1_barcode, "-o", snakemake.output[0], snakemake.input[0]])
    else:
        subprocess.call(["cutadapt","-j 0", "-m", str(snakemake.params.minlen), "-M", str(snakemake.params.maxlen),
                         "-g", r1_primer, "-o", snakemake.output[0], snakemake.input[0]])