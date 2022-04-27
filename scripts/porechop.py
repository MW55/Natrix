import pandas as pd
import logging
import subprocess
from Bio.Seq import Seq

primer_table = pd.read_csv(snakemake.input.primer_t, index_col="Probe",na_filter=False).to_dict("index")

f_primers = primer_table[snakemake.wildcards.sample + "_" + snakemake.wildcards.unit]["specific_forward_primer"].split("-")
while len(f_primers) < 2:
	f_primers.append("")
r_primers = primer_table[snakemake.wildcards.sample + "_" + snakemake.wildcards.unit]["specific_reverse_primer"].split("-")
while len(r_primers) < 2:
	r_primers.append("")
f_barcode = primer_table[snakemake.wildcards.sample + "_" + snakemake.wildcards.unit]["Barcode_forward"]
r_barcode = primer_table[snakemake.wildcards.sample + "_" + snakemake.wildcards.unit]["Barcode_reverse"]

f_adapter = "".join([f_primers[0],f_barcode,f_primers[1]])
r_adapter = "".join([r_primers[0],r_barcode,r_primers[1]])

f_adapter_rev_complement = str(Seq(f_adapter).reverse_complement())
r_adapter_rev_complement = str(Seq(r_adapter).reverse_complement())


process = subprocess.run(["porechop",
    "-i", snakemake.input[0], "-o", snakemake.output[0],
    "--adapter_forward", "-".join([f_adapter, r_adapter]),
    "--adapter_reverse", "-".join([f_adapter_rev_complement, r_adapter_rev_complement]),
    "-t", str(snakemake.threads),
    "--trimmed_only", "true",
    "--tail_crop", str(snakemake.params["tail_crop"]),
    "--head_crop", str(snakemake.params["head_crop"]),
    "--min_length", str(snakemake.params["min_length"]),
    "--max_length", str(snakemake.params["max_length"]),
    "--custom_adapter_names", "-".join(["forward","reverse"]),
    "--discard_middle",
    "--correct_read_direction"
    ])


