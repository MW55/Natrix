import subprocess
import pandas as pd

# Read the primer table for completeness (used here only for logging or potential future use)
primer_table = pd.read_csv(snakemake.input.primer_t, index_col="Probe", na_filter=False).to_dict("index")
r1_primer = primer_table[snakemake.wildcards.sample + "_" + snakemake.wildcards.unit]["specific_forward_primer"]
r2_primer = primer_table[snakemake.wildcards.sample + "_" + snakemake.wildcards.unit]["specific_reverse_primer"]

logfile = open(str(snakemake.log), "w")

minlen = snakemake.params.minlen
maxlen = snakemake.params.maxlen  # Note: fastp does not directly support maximum length filtering.

# fastp performs automatic adapter detection even in single-end mode.
if snakemake.params.paired_end:
    if snakemake.params.prim_rm:
        # If primer removal is set, mimic a no-op (or simply rename the input files)
        subprocess.call(["mv", snakemake.input[0], snakemake.output[0]])
        subprocess.call(["mv", snakemake.input[1], snakemake.output[1]])
    else:
        # Build fastp command for paired-end reads.
        cmd = [
            "fastp",
            "-i", snakemake.input[0],
            "-I", snakemake.input[1],
            "-o", snakemake.output[0],
            "-O", snakemake.output[1],
            "--length_required", str(minlen),
            "--detect_adapter_for_pe"
        ]
        # Optionally log the primer sequences (if needed)
        logfile.write("Paired-end fastp command: " + " ".join(cmd) + "\n")
        subprocess.call(cmd, stdout=logfile)
else:
    if snakemake.params.prim_rm:
        subprocess.call(["mv", snakemake.input[0], snakemake.output[0]])
    else:
        # Build fastp command for single-end reads.
        cmd = [
            "fastp",
            "-i", snakemake.input[0],
            "-o", snakemake.output[0],
            "--length_required", str(minlen),
            "--detect_adapter_for_pe"
        ]
        logfile.write("Single-end fastp command: " + " ".join(cmd) + "\n")
        subprocess.call(cmd, stdout=logfile)

logfile.close()
